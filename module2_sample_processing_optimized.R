# module2_sample_processing_optimized.R
# 模块2: 单样品处理 (降维、聚类、标记、注释、可视化)
# 版本日期：2025.05.14 (已修正合并逻辑)

# ============================================================================
# 模块2主要功能：
# 1. (已在模块1完成) 合并样品复制
# 2. 样品标准化、特征选择、降维、聚类
# 3. 标记基因查找
# 4. 细胞类型注释
# 5. 可视化
# ============================================================================

# 加载必要的库（只保留模块2特有的库，其他通过主脚本和utils.R加载）
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  # 其他库不需要重复加载
})

# 确保 utils.R 和 clustering_evaluation.R 已由主脚本加载
# source("utils.R") # 由主脚本加载
# source("clustering_evaluation.R")  # 由主脚本加载

#' @title 处理单个样品 (已由模块1合并技术复制)
#' @description 对单个样品进行标准化、变量特征查找、数据缩放和PCA。
#'              模块1的输出已经是合并了技术复制的Seurat对象。
#' @param seu Seurat对象 (单个样品，技术复制已合并)
#' @param npcs 主成分个数
#' @param nfeatures 高变基因数量
#' @return 处理后的Seurat对象
process_single_sample <- function(seu, npcs, nfeatures = 3000) {
  message("  Step 2.1: Processing sample (Normalize, FindVariableFeatures, Scale, PCA)...")
  # 确保输入是有效的Seurat对象
  if (is.null(seu) || !inherits(seu, "Seurat")) {
    stop("process_single_sample: 输入的seu不是有效的Seurat对象。")
  }
  
  seu <- NormalizeData(seu, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) %>%
    ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE) %>% # 考虑从config读取vars.to.regress
    RunPCA(npcs = npcs, verbose = FALSE)
  message("    Processing complete.")
  return(seu)
}

#' @title 确定最佳PCA维度
#' @description 基于累积方差贡献率确定最佳PCA维度
#' @param seu Seurat对象，需要包含PCA结果
#' @param variance_threshold 方差贡献阈值 (%)
#' @param min_dim 最小维度
#' @param max_dim 最大维度
#' @return 最佳维度数
determine_optimal_dim <- function(seu, variance_threshold = 75, min_dim = 10, max_dim = 35) {
  message("  Step 2.2: Determining optimal PCA dimensions...")
  if (!("pca" %in% names(seu@reductions))) {
    stop("PCA results not found in Seurat object.")
  }
  pca_stdev <- seu@reductions$pca@stdev
  if (length(pca_stdev) == 0) {
    stop("PCA standard deviations are empty.")
  }
  pca_var <- pca_stdev^2
  pca_var_pct <- pca_var / sum(pca_var) * 100
  pca_cumvar <- cumsum(pca_var_pct)
  
  optimal_dim <- which(pca_cumvar >= variance_threshold)[1]
  if (is.na(optimal_dim)) {
    optimal_dim <- max_dim # 如果阈值未达到，使用max_dim
    message("    Variance threshold (", variance_threshold, "%) not met with available PCs. Using max_dim: ", max_dim)
  } else {
    optimal_dim <- min(max(optimal_dim, min_dim), max_dim) # 将维度限制在范围内
    message(paste0("    Optimal PCA dimension determined: ", optimal_dim, 
                   " (explains ", round(pca_cumvar[optimal_dim], 1), "% variance)"))
  }
  return(optimal_dim)
}

#' @title 运行降维（UMAP、邻居图）
#' @description 在确定最佳PCA维度后运行UMAP降维和FindNeighbors
#' @param seu Seurat对象
#' @param optimal_dim 最佳PC维度
#' @param viz_params 可视化参数 (来自config)
#' @return 执行降维后的Seurat对象
run_dim_reduction <- function(seu, optimal_dim, viz_params) {
  message("  Step 2.3: Running dimensionality reduction (UMAP, FindNeighbors)...")
  
  reduction_to_use <- "pca" # UMAP和Neighbors基于PCA结果
  
  cell_count <- ncol(seu)
  message(paste0("    Cell count for UMAP/Neighbors: ", cell_count))
  
  # UMAP参数设置
  umap_params <- list(
    reduction = reduction_to_use,
    dims = 1:optimal_dim,
    reduction.name = "umap", # 标准UMAP名称
    verbose = FALSE
  )
  # FindNeighbors参数设置
  neighbor_params <- list(
    reduction = reduction_to_use,
    dims = 1:optimal_dim,
    verbose = FALSE
  )
  
  # 根据细胞数量调整UMAP和FindNeighbors参数 (来自config$common_viz_params或$module2$viz_params)
  if (cell_count > (viz_params$large_dataset_threshold %||% 15000) && !is.null(viz_params$large_umap_neighbors)) {
    message("    Using large dataset UMAP/Neighbors parameters...")
    umap_params$n.neighbors <- viz_params$large_umap_neighbors
    umap_params$min.dist <- viz_params$large_umap_min_dist
    umap_params$repulsion.strength <- viz_params$large_repulsion %||% 1 # uwot的默认值
    umap_params$spread <- viz_params$large_umap_spread %||% 1 # uwot的默认值
    
    neighbor_params$k.param <- viz_params$large_neighbor_k %||% 20 # FindNeighbors的k参数
  }
  
  # 运行UMAP
  message("    Running UMAP with dims: 1:", optimal_dim, 
          if(!is.null(umap_params$n.neighbors)) paste0(", n.neighbors: ", umap_params$n.neighbors) else "",
          if(!is.null(umap_params$min.dist)) paste0(", min.dist: ", umap_params$min.dist) else "")
  seu <- do.call(RunUMAP, c(object = seu, umap_params))
  
  # 运行FindNeighbors
  message("    Running FindNeighbors with dims: 1:", optimal_dim,
          if(!is.null(neighbor_params$k.param)) paste0(", k.param: ", neighbor_params$k.param) else "")
  seu <- do.call(FindNeighbors, c(object = seu, neighbor_params))
  
  message("    Dimensionality reduction complete.")
  return(seu)
}

#' @title 确定最佳聚类分辨率
#' @description 使用clustering_evaluation.R中的函数确定最佳聚类分辨率
#' @param seu Seurat对象
#' @param optimal_dim 最佳PC维度
#' @param resolution_range 要评估的分辨率范围 (来自config)
#' @param clustering_metric_weights 聚类评估指标权重 (来自config)
#' @param sample_id 当前样品ID (用于绘图和保存)
#' @param sample_dir 当前样品的结果目录 (用于保存图表)
#' @return 最佳分辨率值
determine_optimal_res_with_eval <- function(seu, optimal_dim, resolution_range, clustering_metric_weights, sample_id, sample_dir) {
  message("  Step 2.4: Determining optimal clustering resolution using evaluation metrics...")
  
  # 确保FindNeighbors已运行，并获取图名称
  # FindNeighbors默认创建 AssayName_snn 和 AssayName_knn
  # 我们需要SNN图进行聚类
  assay_name <- DefaultAssay(seu) %||% "RNA" # 获取当前默认Assay
  graph_name_to_use <- paste0(assay_name, "_snn") 
  
  if (!graph_name_to_use %in% names(seu@graphs)) {
    warning(paste0("SNN graph '", graph_name_to_use, "' not found. Will attempt to use the first available graph or re-run FindNeighbors."))
    if (length(names(seu@graphs)) > 0) {
      graph_name_to_use <- names(seu@graphs)[1] # 使用第一个可用的图
      message(paste0("    Using first available graph: ", graph_name_to_use))
    } else {
      message("    No graphs found. Re-running FindNeighbors with default PCA reduction.")
      seu <- FindNeighbors(seu, reduction = "pca", dims = 1:optimal_dim, verbose = FALSE)
      graph_name_to_use <- paste0(assay_name, "_snn") # 再次尝试默认
      if (!graph_name_to_use %in% names(seu@graphs)) {
        stop("Failed to find or create an SNN graph for clustering.")
      }
    }
  }
  message(paste0("    Using graph '", graph_name_to_use, "' for clustering evaluation."))
  
  # 评估不同分辨率
  # clustering_evaluation.R 中的 evaluate_clustering 函数
  clustering_metrics <- evaluate_clustering(
    seu_obj = seu, 
    resolutions = resolution_range,
    dims = 1:optimal_dim, # 确保传递正确的dims给内部的FindClusters
    reduction_method_for_graph = "pca", # 假设SNN图是基于PCA构建的
    graph_name = graph_name_to_use, 
    compute_silhouette = TRUE,
    compute_modularity = TRUE,
    compute_db_index = TRUE
  )
  
  # 确定最佳分辨率
  # clustering_evaluation.R 中的 determine_optimal_resolution 函数
  optimal_res <- determine_optimal_resolution(
    eval_metrics = clustering_metrics,
    min_clusters = 5, # 可从config读取
    prefer_higher_clusters = TRUE, # 可从config读取
    weights = clustering_metric_weights 
  )
  
  # 可视化评估指标
  # clustering_evaluation.R 中的 plot_clustering_metrics 函数
  p_metrics <- plot_clustering_metrics(
    eval_metrics = clustering_metrics, 
    optimal_res = optimal_res,
    title = paste0(sample_id, " - Clustering Evaluation Metrics")
  )
  
  # 保存评估图表 (使用utils.R中的save_plot_safe)
  metrics_plot_path <- file.path(sample_dir, paste0(sample_id, "_clustering_metrics.pdf"))
  save_plot_safe(metrics_plot_path, p_metrics, width = 12, height = 10)
  
  message("    Optimal resolution determined: ", optimal_res)
  return(optimal_res)
}


#' @title 执行最终聚类
#' @description 使用确定的最佳分辨率执行最终聚类
#' @param seu Seurat对象
#' @param optimal_res 最佳分辨率
#' @return 聚类后的Seurat对象
run_final_clustering <- function(seu, optimal_res) {
  message("  Step 2.5: Running final clustering with resolution ", optimal_res)
  
  assay_name <- DefaultAssay(seu) %||% "RNA"
  graph_name_to_use <- paste0(assay_name, "_snn")
  
  if (!graph_name_to_use %in% names(seu@graphs)) {
    stop(paste0("SNN graph '", graph_name_to_use, "' not found for final clustering."))
  }
  message(paste0("    Using graph '", graph_name_to_use, "' for final clustering."))
  
  # 执行聚类
  seu <- FindClusters(seu, resolution = optimal_res, graph.name = graph_name_to_use, verbose = FALSE)
  
  # Seurat v5 FindClusters 会将结果存储在 seu$seurat_clusters (如果resolution是单个值)
  # 或 seu[[paste0(graph_name, "_res.", resolution)]]
  # 我们需要确保 seu$seurat_clusters 和 Idents(seu) 被正确设置
  
  cluster_col_name <- paste0(graph_name_to_use, "_res.", optimal_res) # Seurat v3/v4 风格
  
  if (cluster_col_name %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- seu@meta.data[[cluster_col_name]]
    message(paste0("    'seurat_clusters' column updated from '", cluster_col_name, "'."))
  } else if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    # 如果 FindClusters (v5) 直接更新了 seurat_clusters, 确保它对应于 optimal_res
    # 这部分可能需要更复杂的检查，但通常如果 optimal_res 是单个值，FindClusters会处理好
    message(paste0("    'seurat_clusters' column already exists. Assuming it corresponds to resolution ", optimal_res, "."))
  } else {
    warning(paste0("Could not find cluster column for the chosen optimal resolution '", cluster_col_name, 
                   "' and 'seurat_clusters' is also missing. Clustering might not have run as expected."))
    # 作为回退，尝试使用FindClusters可能创建的任何默认聚类列
    default_cluster_cols <- grep(paste0(graph_name_to_use, "_res\\."), colnames(seu@meta.data), value = TRUE)
    if (length(default_cluster_cols) > 0) {
      seu$seurat_clusters <- seu[[default_cluster_cols[1]]] # 使用第一个找到的
      warning(paste0("    Fallback: Using first available resolution column '", default_cluster_cols[1], "' as 'seurat_clusters'."))
    } else {
      # 如果实在没有，这可能是个问题
      seu$seurat_clusters <- factor(rep(1, ncol(seu))) # 创建一个虚拟聚类
      warning("    Critical: No clustering columns found. Created a dummy 'seurat_clusters' with all cells in cluster 1.")
    }
  }
  
  # 确保seurat_clusters为因子，按数值顺序排序水平
  if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    current_clusters <- seu$seurat_clusters
    if (!is.factor(current_clusters)) {
      unique_vals <- sort(as.numeric(as.character(unique(current_clusters))))
      seu$seurat_clusters <- factor(current_clusters, levels = as.character(unique_vals))
      message("    Converted 'seurat_clusters' to factor with numerically sorted levels.")
    } else {
      current_levels <- levels(current_clusters)
      numeric_levels <- suppressWarnings(as.numeric(current_levels))
      if (!anyNA(numeric_levels) && !identical(current_levels, as.character(sort(numeric_levels)))) {
        seu$seurat_clusters <- factor(current_clusters, levels = as.character(sort(numeric_levels)))
        message("    Re-sorted factor levels for 'seurat_clusters' numerically.")
      }
    }
    Idents(seu) <- "seurat_clusters" # 设置主身份
  }
  
  n_clusters <- length(unique(Idents(seu)))
  message(paste0("    Final clustering complete. Found ", n_clusters, " clusters."))
  
  # 尝试应用JoinLayers (Seurat v5)
  tryCatch({
    if(utils::packageVersion("Seurat") >= "5.0.0" && 
       (exists("JoinLayers", envir = asNamespace("Seurat")) || exists("JoinLayers", envir = asNamespace("SeuratObject")))) {
      message("    Applying JoinLayers (Seurat v5)...")
      if (exists("JoinLayers", envir = asNamespace("SeuratObject"))) {
        seu <- SeuratObject::JoinLayers(seu)
      } else {
        seu <- Seurat::JoinLayers(seu)
      }
    } else {
      message("    Skipping JoinLayers (not Seurat v5 or function not found).")
    }
  }, error = function(e) {
    message("    Skipping JoinLayers due to error: ", e$message)
  })
  
  return(seu)
}


#' @title 查找样品标记基因
#' @description 为样品中的每个聚类查找标记基因 (使用utils.R中的find_markers)
#' @param seu Seurat对象
#' @param marker_params 标记基因参数 (来自config)
#' @return 包含标记基因的数据框
find_sample_markers_util <- function(seu, marker_params) {
  message("  Step 2.6: Finding sample markers using utility function...")
  # 调用 utils.R 中的 find_markers 函数
  # 注意：utils.R中的find_markers函数名可能与此不同，需确认
  # 假设 utils.R 中有一个名为 find_all_markers_robust 的函数
  if (exists("find_markers", where = "package:Seurat")) { # 避免与Seurat::FindMarkers冲突
    # 如果 utils.R 定义了更健壮的 find_markers, 确保它被调用
    # 假设 utils.R 中的函数是 find_sample_markers_core 或类似名称
    if (exists("find_markers_core_from_utils", envir = .GlobalEnv)) { # 假设 utils.R 中函数叫这个
      all_markers <- find_markers_core_from_utils(seu, marker_params)
    } else {
      # 如果 utils.R 中没有特定的 find_markers, 则使用标准 FindAllMarkers
      message("    Using Seurat::FindAllMarkers as no custom utility function was found.")
      all_markers <- Seurat::FindAllMarkers(
        seu,
        only.pos = marker_params$only.pos %||% TRUE,
        min.pct = marker_params$min.pct %||% 0.25,
        logfc.threshold = marker_params$logfc.threshold %||% 0.25,
        test.use = marker_params$test.use %||% "wilcox",
        verbose = FALSE,
        max.cells.per.ident = marker_params$max.cells.per.ident %||% Inf # 默认不限制
      )
    }
  } else {
    stop("Seurat package or FindAllMarkers function not available.")
  }
  
  if(!is.null(all_markers) && nrow(all_markers) > 0) {
    if(!"gene" %in% colnames(all_markers) && !is.null(rownames(all_markers))) {
      all_markers$gene <- rownames(all_markers) # 确保gene列存在
    }
    # 进一步过滤和处理可以在这里或find_markers_core_from_utils内部完成
    message(paste0("    Found ", nrow(all_markers), " potential markers."))
  } else {
    message("    No markers found.")
    return(data.frame()) # 返回空数据框
  }
  return(all_markers)
}

#' @title 为标记基因添加注释信息 (使用utils.R中的函数)
#' @description 将MSU注释信息添加到标记基因数据框
#' @param markers_df 标记基因数据框
#' @param msu_annot MSU注释数据框 (已加载)
#' @return 添加注释信息后的标记基因数据框
add_annotation_to_markers_util <- function(markers_df, msu_annot) {
  message("    Adding annotation to markers using utility function...")
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    message("    No markers to annotate.")
    return(markers_df)
  }
  if (is.null(msu_annot) || nrow(msu_annot) == 0) {
    message("    MSU annotation data is missing. Skipping annotation.")
    return(markers_df)
  }
  # 调用 utils.R 中的 add_annotation_to_markers 函数
  # 假设 utils.R 中的函数名就是 add_annotation_to_markers
  if (exists("add_annotation_to_markers", envir = .GlobalEnv)) {
    return(add_annotation_to_markers(markers_df, msu_annot))
  } else {
    warning("Utility function 'add_annotation_to_markers' not found. Skipping annotation.")
    return(markers_df) # 原样返回
  }
}


#' @title 查找、注释和可视化标记基因
#' @description 为样品找到标记基因，注释聚类，生成可视化图表
#' @param seu Seurat对象
#' @param sample_id 样品ID
#' @param marker_params 标记基因参数 (来自config)
#' @param annotation_params 注释参数 (来自config)
#' @param msu_annot MSU注释数据 (已加载)
#' @param viz_params 可视化参数 (来自config)
#' @param sample_dir 样品结果目录
#' @return 包含注释的Seurat对象
find_annotate_visualize_markers <- function(seu, sample_id, marker_params, annotation_params, msu_annot, viz_params, sample_dir) {
  message("  Step 2.7: Finding, annotating, and visualizing markers for ", sample_id)
  
  # ---- 生成并保存带有聚类的UMAP图 (使用utils.R中的函数) ----
  # utils.R 中的 generate_visualization 函数
  if (exists("generate_visualization", envir = .GlobalEnv)) {
    p_umap <- generate_visualization(
      seu, 
      type = "umap", 
      params = list(annotation_col="seurat_clusters", title=paste0(sample_id, " - Clusters (res: ", seu@meta.data$seurat_clusters_resolution[1] %||% "NA", ")")), # 假设分辨率存储在元数据中
      viz_params = viz_params
    )
    umap_plot_path <- file.path(sample_dir, paste0(sample_id, "_cluster_umap.pdf"))
    save_plot_safe(umap_plot_path, p_umap, width = 10, height = 8)
  } else {
    warning("generate_visualization function not found. Skipping UMAP plot.")
  }
  
  # ---- 查找标记基因 ----
  message("    Finding markers for sample ", sample_id)
  markers <- find_sample_markers_util(seu, marker_params) # 使用封装的工具函数
  
  # 如果找到标记基因，进行处理、注释和可视化
  if(nrow(markers) > 0) {
    markers$gene <- fix_gene_labels(markers$gene) # utils.R
    
    markers_to_save <- markers
    if("cluster" %in% colnames(markers_to_save)) {
      markers_to_save$cluster_num <- suppressWarnings(as.numeric(as.character(markers_to_save$cluster)))
      if(!anyNA(markers_to_save$cluster_num)) {
        markers_to_save <- markers_to_save %>% arrange(cluster_num, p_val_adj, desc(avg_log2FC %||% avg_logFC))
      } else {
        markers_to_save <- markers_to_save %>% arrange(cluster, p_val_adj, desc(avg_log2FC %||% avg_logFC))
      }
      markers_to_save$cluster_num <- NULL 
    }
    
    # 添加基因注释 (使用封装的工具函数)
    markers_anno <- add_annotation_to_markers_util(markers, msu_annot)
    markers_to_save_anno <- add_annotation_to_markers_util(markers_to_save, msu_annot)
    
    marker_file <- file.path(sample_dir, paste0(sample_id, "_markers_annotated.csv"))
    write.csv(markers_to_save_anno, marker_file, row.names = FALSE)
    message("    Annotated markers saved to: ", marker_file)
    
    # ---- 聚类注释 (使用utils.R中的函数) ----
    message("    Annotating clusters...")
    # utils.R 中的 annotate_clusters 函数
    if (exists("annotate_clusters", envir = .GlobalEnv)) {
      cluster_annotations <- annotate_clusters(markers_anno, msu_annot, annotation_params)
      
      # 将注释添加到Seurat对象元数据
      if (nrow(cluster_annotations) > 0 && "cluster" %in% colnames(cluster_annotations) && "label" %in% colnames(cluster_annotations)) {
        # 创建映射：聚类ID -> 注释标签
        annot_map <- setNames(cluster_annotations$label, as.character(cluster_annotations$cluster))
        
        # 应用映射到细胞
        # 确保 seu$seurat_clusters 是字符或因子，并且与 annot_map 的键类型一致
        current_cell_clusters <- as.character(seu$seurat_clusters)
        seu$custom_annotation <- annot_map[current_cell_clusters]
        seu$custom_annotation[is.na(seu$custom_annotation)] <- "Unannotated" # 处理未匹配到的聚类
        
        message("    'custom_annotation' column added to Seurat object metadata.")
        
        # 保存聚类注释表
        annotation_table_path <- file.path(sample_dir, paste0(sample_id, "_cluster_annotations.csv"))
        write.csv(cluster_annotations, annotation_table_path, row.names = FALSE)
        message("    Cluster annotation table saved to: ", annotation_table_path)
        
      } else {
        message("    Annotation results are empty or malformed. Skipping adding to Seurat metadata.")
        seu$custom_annotation <- "Unannotated" # 默认值
      }
    } else {
      warning("annotate_clusters function not found. Skipping cluster annotation.")
      seu$custom_annotation <- "Unannotated" # 默认值
    }
    
    # ---- 标记可视化 (使用utils.R中的函数) ----
    if(nrow(markers) > 0 && exists("generate_visualization", envir = .GlobalEnv)) {
      tryCatch({
        # 热图
        top_n_heatmap <- annotation_params$markers_per_cluster_heatmap %||% 3 # 从config获取
        message(paste0("    Generating markers heatmap (top ", top_n_heatmap, " per cluster)..."))
        p_heatmap <- generate_visualization(
          seu, "heatmap", 
          list(
            markers_df = markers, 
            n_markers_per_cluster = top_n_heatmap,
            annotation_col = "seurat_clusters", # 或 "custom_annotation" 如果想用注释名
            title = paste0(sample_id, " - Top ", top_n_heatmap, " Markers Per Cluster")
          ), 
          viz_params
        )
        heatmap_plot_path <- file.path(sample_dir, paste0(sample_id, "_markers_heatmap.pdf"))
        save_plot_safe(heatmap_plot_path, p_heatmap, 
                       width = 12, 
                       height = min(15, 4 + 0.25 * length(unique(markers$cluster)) * top_n_heatmap)) # 动态高度
        
        # 气泡图
        top_n_bubble <- annotation_params$markers_per_cluster_bubble %||% 3 # 从config获取
        message(paste0("    Generating markers bubble plot (top ", top_n_bubble, " per cluster)..."))
        top_markers_bubble <- markers %>%
          filter(!is.na(gene) & gene %in% rownames(seu)) %>% # 确保基因在对象中
          group_by(cluster) %>% 
          slice_max(order_by = avg_log2FC %||% avg_logFC, n = top_n_bubble) %>% 
          ungroup() 
        
        if(nrow(top_markers_bubble) > 0) {
          p_bubble <- generate_visualization(
            seu, "bubble", 
            list(
              features = unique(top_markers_bubble$gene), 
              annotation_col = "seurat_clusters", # 或 "custom_annotation"
              title = paste0(sample_id, " - Top Markers Bubble Plot")
            ), 
            viz_params
          )
          bubble_plot_path <- file.path(sample_dir, paste0(sample_id, "_markers_bubbleplot.pdf"))
          save_plot_safe(bubble_plot_path, p_bubble, width = max(10, 0.5 * length(unique(top_markers_bubble$gene))), height = 8) # 动态宽度
        } else {
          message("    No valid markers found for bubble plot after filtering.")
        }
        
        # 特征图 (所有找到的顶尖标记基因，可能会很多，考虑限制数量)
        top_n_feature <- annotation_params$markers_per_cluster_feature %||% 1 # 从config获取
        message(paste0("    Generating top marker feature plots (top ", top_n_feature, " per cluster)..."))
        top_markers_feature <- markers %>%
          filter(!is.na(gene) & gene %in% rownames(seu)) %>%
          group_by(cluster) %>%
          slice_max(order_by = avg_log2FC %||% avg_logFC, n = top_n_feature) %>%
          ungroup() %>%
          distinct(gene) %>% # 取不重复的基因
          slice_head(n = viz_params$max_features_for_featureplot %||% 12) # 限制总基因数
        
        if(nrow(top_markers_feature) > 0) {
          p_feature <- generate_visualization(
            seu, "feature", 
            list(
              features = top_markers_feature$gene, 
              title = paste0(sample_id, " - Top Markers Expression")
            ), 
            viz_params
          )
          feature_plot_path <- file.path(sample_dir, paste0(sample_id, "_markers_featureplot.pdf"))
          save_plot_safe(feature_plot_path, p_feature, width=12, height=max(8, 2*ceiling(nrow(top_markers_feature)/(viz_params$feature_ncol %||% 3)))) # 动态高度
        } else {
          message("    No valid markers found for feature plot after filtering.")
        }
      }, error = function(e) {
        message("    Error during marker visualization: ", e$message)
      })
    } else {
      message("    No markers found or generate_visualization function missing. Skipping marker visualizations.")
    }
  } else {
    message("    No markers found to annotate or visualize.")
    # 创建一个空的markers文件以避免下游处理错误
    empty_markers <- data.frame(p_val=numeric(0), avg_log2FC=numeric(0), pct.1=numeric(0), pct.2=numeric(0), p_val_adj=numeric(0), cluster=character(0), gene=character(0))
    marker_file <- file.path(sample_dir, paste0(sample_id, "_markers_annotated.csv"))
    write.csv(empty_markers, marker_file, row.names = FALSE)
    message("    Empty markers file saved to: ", marker_file)
    seu$custom_annotation <- "Unannotated" # 默认值
  }
  
  return(seu)
}

#' @title 模块2: 单样品处理主函数
#' @description 对通过质控和双细胞过滤的样本进行标准分析流程。
#'              模块1的输出 seu_clean_list 已经是合并了技术复制的样本列表。
#' @param seu_clean_list 通过质控和双细胞过滤的样本列表 (sample_id -> SeuratObject)
#' @param results_dir 结果保存主目录
#' @param processing_params 处理参数 (来自config)
#' @param marker_params 标记基因参数 (来自config)
#' @param annotation_params 注释参数 (来自config)
#' @param msu_annot MSU注释数据 (已加载)
#' @param viz_params 可视化参数 (来自config)
#' @return 处理后的Seurat对象列表 (瘦身版)
run_sample_processing <- function(seu_clean_list, 
                                  results_dir, 
                                  processing_params, 
                                  marker_params, 
                                  annotation_params, 
                                  msu_annot, 
                                  viz_params) {
  
  message("=== Running Module 2: Single Sample Processing ===")
  
  processed_samples_list <- list() # 用于存储最终处理好的（瘦身版）Seurat对象
  
  # 遍历每个样品 (seu_clean_list 的每个元素是一个已合并技术复制的Seurat对象)
  for (sample_id in names(seu_clean_list)) {
    
    current_seu <- seu_clean_list[[sample_id]] # 这是单个Seurat对象
    message(paste0("\nProcessing Sample: ", sample_id, " - Start Time: ", Sys.time()))
    
    # 为这个样品创建结果目录 (在主results_dir下按sample_id分子目录)
    sample_result_dir <- file.path(results_dir, "Module2_SingleSampleOutput", sample_id) # 修改了路径
    dir.create(sample_result_dir, recursive = TRUE, showWarnings = FALSE)
    
    processed_seu_for_sample <- NULL # 初始化当前样品的处理结果
    
    tryCatch({
      # --- 步骤 2.0: 检查输入的Seurat对象 ---
      if (is.null(current_seu) || !inherits(current_seu, "Seurat")) {
        message("  Skipping sample ", sample_id, ": Invalid Seurat object received from Module 1.")
        next # 跳到下一个sample_id
      }
      message("  Input Seurat object for sample ", sample_id, ": ", ncol(current_seu), " cells, ", nrow(current_seu), " features.")
      
      # --- 基因名标准化 (以防万一) ---
      current_seu <- standardize_gene_names(current_seu) # utils.R
      
      # --- 步骤 2.1: 核心处理 (Normalize, FindVariableFeatures, Scale, PCA) ---
      message("  Step 2.1: process_single_sample - Start Time: ", Sys.time())
      current_seu <- process_single_sample(current_seu, 
                                           npcs = processing_params$npcs, 
                                           nfeatures = processing_params$nfeatures)
      message("  Step 2.1: process_single_sample - End Time: ", Sys.time())
      
      # --- 步骤 2.2: 确定最佳PCA维度 --- 
      message("  Step 2.2: Determine/Set optimal_dim - Start Time: ", Sys.time())
      # 允许在config中为特定样品预设optimal_dim，或全局设置
      sample_specific_optimal_dim <- processing_params$sample_optimal_dims[[sample_id]]
      global_optimal_dim_setting <- processing_params$optimal_dim 
      
      if (!is.null(sample_specific_optimal_dim) && is.numeric(sample_specific_optimal_dim)) {
        optimal_dim <- sample_specific_optimal_dim
        message("    Using pre-defined optimal PCA dimension for ", sample_id, ": ", optimal_dim)
      } else if (!is.null(global_optimal_dim_setting) && is.numeric(global_optimal_dim_setting)) {
        optimal_dim <- global_optimal_dim_setting
        message("    Using globally pre-defined optimal PCA dimension: ", optimal_dim)
      } else {
        optimal_dim <- determine_optimal_dim(current_seu, 
                                             variance_threshold = processing_params$variance_threshold %||% 75,
                                             min_dim = processing_params$min_dim %||% 10, 
                                             max_dim = processing_params$max_dim %||% 35)
      }
      message("  Step 2.2: Determine/Set optimal_dim - End Time: ", Sys.time())
      current_seu@misc$optimal_dim <- optimal_dim # 存储选择的维度
      
      # --- 步骤 2.3: 运行降维 (UMAP, FindNeighbors) ---
      message("  Step 2.3: run_dim_reduction - Start Time: ", Sys.time())
      current_seu <- run_dim_reduction(current_seu, optimal_dim, viz_params)
      message("  Step 2.3: run_dim_reduction - End Time: ", Sys.time())
      
      # --- 步骤 2.4 & 2.5: 确定最佳分辨率并聚类 ---
      message("  Step 2.4/2.5: Determine optimal resolution and run clustering - Start Time: ", Sys.time())
      optimal_res <- NULL
      # 允许在config中为特定样品预设resolution，或全局设置
      sample_specific_resolution <- processing_params$sample_resolutions[[sample_id]]
      global_resolution_setting <- processing_params$fixed_resolution 
      
      if (!is.null(sample_specific_resolution) && is.numeric(sample_specific_resolution)) {
        optimal_res <- sample_specific_resolution
        message("    Using pre-defined resolution for ", sample_id, ": ", optimal_res)
        current_seu <- run_final_clustering(current_seu, optimal_res)
      } else if (!is.null(global_resolution_setting) && is.numeric(global_resolution_setting)) {
        optimal_res <- global_resolution_setting
        message("    Using globally pre-defined resolution: ", optimal_res)
        current_seu <- run_final_clustering(current_seu, optimal_res)
      } else if (!is.null(processing_params$enable_resolution_optimization) && 
                 !processing_params$enable_resolution_optimization) {
        optimal_res <- processing_params$default_resolution %||% 0.6 # 使用config中的默认值
        message("    Resolution optimization disabled. Using default resolution: ", optimal_res)
        current_seu <- run_final_clustering(current_seu, optimal_res)
      } else {
        optimal_res <- determine_optimal_res_with_eval(
          current_seu, 
          optimal_dim, 
          resolution_range = processing_params$resolution_range,
          clustering_metric_weights = processing_params$clustering_metric_weights, # 来自config
          sample_id = sample_id,
          sample_dir = sample_result_dir # 传递目录用于保存图表
        )
        current_seu <- run_final_clustering(current_seu, optimal_res)
      }
      message("  Step 2.4/2.5: Clustering - End Time: ", Sys.time())
      current_seu@misc$optimal_resolution <- optimal_res # 存储选择的分辨率
      # 将实际使用的分辨率也存储在元数据中，方便后续查看
      if("seurat_clusters" %in% colnames(current_seu@meta.data)){
        current_seu$seurat_clusters_resolution <- optimal_res
      }
      
      
      # --- 步骤 2.6 & 2.7: 查找、注释和可视化标记 --- 
      message("  Step 2.6/2.7: find_annotate_visualize_markers - Start Time: ", Sys.time())
      current_seu <- find_annotate_visualize_markers(
        current_seu, sample_id, 
        marker_params, annotation_params, msu_annot, 
        viz_params, sample_result_dir
      )
      message("  Step 2.6/2.7: find_annotate_visualize_markers - End Time: ", Sys.time())
      
      processed_seu_for_sample <- current_seu # 保存处理好的对象
      
    }, error = function(e) {
      message(paste0("  ERROR processing sample ", sample_id, ": ", e$message))
      # 可以在这里记录更详细的错误日志，或者将空的/部分处理的对象存起来
      # e.g., saveRDS(current_seu, file.path(sample_result_dir, paste0(sample_id, "_error_partial.rds")))
      processed_seu_for_sample <- NULL 
    })
    
    # --- 对象瘦身与保存 ---
    if (!is.null(processed_seu_for_sample) && inherits(processed_seu_for_sample, "Seurat")) {
      message("  Slimming down processed Seurat object for sample: ", sample_id)
      # 使用 utils.R 中的 slim_seurat 函数
      if (exists("slim_seurat", envir = .GlobalEnv)) {
        slimmed_seu <- slim_seurat(
          processed_seu_for_sample,
          assays_to_keep = processing_params$slim_assays_to_keep %||% "RNA",
          reductions_to_keep = processing_params$slim_reductions_to_keep %||% c("pca", "umap"),
          clean_scale_data = processing_params$slim_clean_scale_data %||% TRUE,
          keep_graphs = processing_params$slim_keep_graphs %||% FALSE
        )
      } else {
        warning("slim_seurat function not found in utils.R. Using internal slimming logic (less configurable).")
        # Fallback to a basic slimming if slim_seurat is not available
        # (可以复制/粘贴模块2原始脚本中的瘦身逻辑到这里作为备用)
        # 这里为了简洁，假设slim_seurat总是可用的
        slimmed_seu <- processed_seu_for_sample # 如果没有slim_seurat，则不瘦身或使用简单瘦身
        if ("scale.data" %in% slotNames(slimmed_seu[['RNA']])) {
          slimmed_seu[['RNA']]@scale.data <- methods::new("dgCMatrix")
        }
        slimmed_seu@graphs <- list()
      }
      
      # 保存瘦身后的RDS文件到 Module2_RDS_Files_Slimmed 目录
      slim_rds_dir <- file.path(results_dir, "Module2_RDS_Files_Slimmed")
      dir.create(slim_rds_dir, showWarnings = FALSE, recursive = TRUE)
      slim_rds_path <- file.path(slim_rds_dir, paste0(sample_id, "_module2_slim.rds"))
      
      message("    Saving slimmed individual sample RDS to: ", slim_rds_path)
      tryCatch({
        saveRDS(slimmed_seu, file = slim_rds_path)
        message("      -> Successfully saved slimmed RDS for ", sample_id, 
                " (Cells: ", ncol(slimmed_seu), ", Features: ", nrow(slimmed_seu), ")")
        processed_samples_list[[sample_id]] <- slimmed_seu # 将瘦身对象存入返回列表
      }, error = function(e){
        message("      -> ERROR saving slimmed RDS for sample ", sample_id, ": ", e$message)
        processed_samples_list[[sample_id]] <- NULL # 保存失败则不加入列表
      })
    } else {
      message(paste0("  Sample ", sample_id, " processing failed or did not produce a valid Seurat object. Skipping slimming and saving."))
      processed_samples_list[[sample_id]] <- NULL
    }
    
    gc() # 在每个样品处理后进行垃圾回收
  } # 结束样品循环
  
  message("\n=== Module 2 Finished ===")
  
  # 返回处理后的瘦身Seurat对象列表
  return(processed_samples_list)
}
