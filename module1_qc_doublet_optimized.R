# module1_qc_doublet_optimized.R
# 模块1: 数据质控与双细胞检测 (优化版)
# 适配 config.R 和 utils.R

# --- 1. 加载必要的库 ---
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("模块1错误: Seurat包未安装。")
  library(Seurat)
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) stop("模块1错误: DoubletFinder包未安装。")
  library(DoubletFinder)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("模块1错误: dplyr包未安装。")
  library(dplyr)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("模块1错误: ggplot2包未安装。")
  library(ggplot2)
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("模块1错误: patchwork包未安装。")
  library(patchwork)
})

# --- 2. 辅助函数 ---

#' @title 计算并应用QC指标 (优化版)
#' @description 计算线粒体/核糖体比例，生成QC图，并过滤细胞。
#' @param seu_rep 单个Seurat对象 (代表一个复制)。
#' @param qc_params QC参数列表 (来自config$module1$qc_params)。
#' @param rep_name 复制的名称 (例如 "SampleA_rep1")。
#' @param rep_plot_dir 此复制的图表输出目录。
#' @param utils_path utils.R的路径 (如果不在默认搜索路径中)。
#' @return 经过QC过滤的Seurat对象。
calculate_and_apply_qc_optimized <- function(seu_rep, qc_params, rep_name, rep_plot_dir, utils_path = "utils.R") {
  # 确保 utils.R 中的 save_plot_safe 可用
  if (file.exists(utils_path)) source(utils_path, local = TRUE)
  if (!exists("save_plot_safe")) stop("calculate_and_apply_qc_optimized: save_plot_safe function not found from utils.R.")
  
  message(paste0("      Step 2.1: 计算QC指标 for ", rep_name))
  
  # 计算线粒体基因比例
  if (!is.null(qc_params$mt_pattern) && nzchar(qc_params$mt_pattern)) {
    tryCatch({
      seu_rep[['percent.mt']] <- PercentageFeatureSet(seu_rep, pattern = qc_params$mt_pattern)
    }, error = function(e) {
      message(paste0("        Warning: 计算线粒体基因比例失败 for ", rep_name, " using pattern '", qc_params$mt_pattern, "'. Error: ", e$message))
      seu_rep[['percent.mt']] <- 0 # 设置为0以避免后续错误
    })
  } else {
    message(paste0("        Skipping mitochondrial calculation for ", rep_name, " (mt_pattern is null or empty)."))
    seu_rep[['percent.mt']] <- 0
  }
  
  # 计算核糖体基因比例
  if (!is.null(qc_params$ribo_pattern) && nzchar(qc_params$ribo_pattern)) {
    tryCatch({
      seu_rep[['percent.ribo']] <- PercentageFeatureSet(seu_rep, pattern = qc_params$ribo_pattern)
    }, error = function(e) {
      message(paste0("        Warning: 计算核糖体基因比例失败 for ", rep_name, " using pattern '", qc_params$ribo_pattern, "'. Error: ", e$message))
      seu_rep[['percent.ribo']] <- 0 # 设置为0以避免后续错误
    })
  } else {
    message(paste0("        Skipping ribosomal calculation for ", rep_name, " (ribo_pattern is null or empty)."))
    # 如果不存在，则添加一个空的列以避免后续的subset错误
    if (!"percent.ribo" %in% colnames(seu_rep@meta.data)) {
      seu_rep[['percent.ribo']] <- 0
    }
  }
  
  # 生成QC小提琴图
  qc_features_to_plot <- c("nFeature_RNA", "nCount_RNA")
  if ("percent.mt" %in% colnames(seu_rep@meta.data) && sum(seu_rep$percent.mt > 0) > 0) qc_features_to_plot <- c(qc_features_to_plot, "percent.mt")
  if ("percent.ribo" %in% colnames(seu_rep@meta.data) && sum(seu_rep$percent.ribo > 0) > 0) qc_features_to_plot <- c(qc_features_to_plot, "percent.ribo")
  
  if (length(qc_features_to_plot) > 0) {
    vln_plot <- VlnPlot(seu_rep, features = qc_features_to_plot, ncol = length(qc_features_to_plot), pt.size = 0.1) +
      plot_annotation(title = paste0("QC Metrics (Before Filtering) - ", rep_name))
    save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_qc_vlnplot_before.jpeg")), vln_plot, width = max(8, 3*length(qc_features_to_plot)), height = 6)
  }
  
  # 生成QC散点图
  scatter_plot1 <- FeatureScatter(seu_rep, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + labs(title = rep_name)
  if ("percent.mt" %in% colnames(seu_rep@meta.data)) {
    scatter_plot2 <- FeatureScatter(seu_rep, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(title = rep_name)
    combined_scatter <- (scatter_plot1 + scatter_plot2) + plot_annotation(title = paste0("QC Scatter (Before Filtering) - ", rep_name))
    save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_qc_scatter_before.jpeg")), combined_scatter, width = 10, height = 5)
  } else {
    save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_qc_scatter_ncount_nfeature_before.jpeg")), scatter_plot1, width = 6, height = 5)
  }
  
  message(paste0("      Step 2.2: 应用QC过滤 for ", rep_name))
  cells_before_qc <- ncol(seu_rep)
  
  # 构建一个逻辑向量，初始全为TRUE
  pass_qc <- rep(TRUE, ncol(seu_rep))
  conditions_applied_count <- 0

  if (!is.null(qc_params$min_features) && "nFeature_RNA" %in% colnames(seu_rep@meta.data)) {
    pass_qc <- pass_qc & (seu_rep$nFeature_RNA > qc_params$min_features)
    conditions_applied_count <- conditions_applied_count + 1
    message(paste0("        Applying min_features > ", qc_params$min_features))
  }
  if (!is.null(qc_params$max_features) && "nFeature_RNA" %in% colnames(seu_rep@meta.data)) {
    pass_qc <- pass_qc & (seu_rep$nFeature_RNA < qc_params$max_features)
    conditions_applied_count <- conditions_applied_count + 1
    message(paste0("        Applying max_features < ", qc_params$max_features))
  }
  if (!is.null(qc_params$max_mt) && "percent.mt" %in% colnames(seu_rep@meta.data)) {
    pass_qc <- pass_qc & (seu_rep$percent.mt < qc_params$max_mt)
    conditions_applied_count <- conditions_applied_count + 1
    message(paste0("        Applying percent.mt < ", qc_params$max_mt))
  }
  if (!is.null(qc_params$max_ribo) && "percent.ribo" %in% colnames(seu_rep@meta.data) && sum(seu_rep$percent.ribo) > 0) {
    pass_qc <- pass_qc & (seu_rep$percent.ribo < qc_params$max_ribo)
    conditions_applied_count <- conditions_applied_count + 1
    message(paste0("        Applying percent.ribo < ", qc_params$max_ribo))
  }
  
  if (conditions_applied_count > 0) {
    seu_filtered <- seu_rep[, pass_qc]
  } else {
    message("        No QC filters applied as parameters are not set appropriately or all thresholds are NULL.")
    seu_filtered <- seu_rep 
  }
  
  cells_after_qc <- ncol(seu_filtered)
  message(paste0("      Step 2.3: QC过滤完成. Before=", cells_before_qc, ", After=", cells_after_qc, " for ", rep_name))
  
  if (cells_after_qc == 0 && cells_before_qc > 0) {
    warning(paste0("QC Alert for ", rep_name, ": 0 cells remaining after QC filtering. Check QC parameters and data quality."))
    # 返回原始对象以避免下游全空错误，但标记它
    seu_rep$qc_passed <- FALSE
    return(seu_rep) 
  }
  seu_filtered$qc_passed <- TRUE
  
  # 过滤后QC图
  if (length(qc_features_to_plot) > 0 && cells_after_qc > 0) {
    vln_plot_after <- VlnPlot(seu_filtered, features = qc_features_to_plot, ncol = length(qc_features_to_plot), pt.size = 0.1) +
      plot_annotation(title = paste0("QC Metrics (After Filtering) - ", rep_name))
    save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_qc_vlnplot_after.jpeg")), vln_plot_after, width = max(8, 3*length(qc_features_to_plot)), height = 6)
  }
  
  return(seu_filtered)
}

#' @title 运行DoubletFinder进行双细胞检测 (对齐双胞2.R逻辑)
run_doublet_detection_optimized <- function(seu_rep_qc, doublet_params, rep_name, rep_plot_dir, sweep_dir, utils_path = "utils.R") {
  if (file.exists(utils_path)) source(utils_path, local = TRUE)
  if (!exists("save_plot_safe")) stop("run_doublet_detection_optimized: save_plot_safe function not found from utils.R.")
  
  message(paste0("      Step 3.1: 为 ", rep_name, " 运行DoubletFinder (对齐双胞2.R逻辑)..." ) )
  ncells <- ncol(seu_rep_qc)
  
  if (ncells <= (doublet_params$min_cells_for_doubletfinder %||% 20)) {
    message(paste0("        Skipping DoubletFinder for ", rep_name, " due to low cell count (", ncells, "). Assigning all as Singlets."))
    seu_rep_qc$doublet_score <- 0
    seu_rep_qc$doublet_class <- "Singlet_low_cell_count"
    return(seu_rep_qc)
  }
  
  # --- 预处理: 恢复为标准流程 (LogNormalize, etc.) ---
  message(paste0("        DEBUG: Preprocessing for DoubletFinder using LogNormalize workflow for ", rep_name))
  DefaultAssay(seu_rep_qc) <- "RNA"
  
  vars_to_regress_param <- doublet_params$vars_to_regress_align 
  if (!is.null(vars_to_regress_param) && !vars_to_regress_param %in% colnames(seu_rep_qc@meta.data)){
      message(paste0("        WARNING: vars.to.regress for ScaleData '", vars_to_regress_param, "' not found in meta.data. Proceeding without regression."))
      vars_to_regress_param <- NULL
  }
  
  message("        DEBUG: Normalizing Data...")
  seu_rep_qc <- NormalizeData(seu_rep_qc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  
  nfeatures_to_find <- doublet_params$nfeatures_align %||% 2500
  message(paste0("        DEBUG: Finding Variable Features (nfeatures=", nfeatures_to_find, ")..."))
  seu_rep_qc <- FindVariableFeatures(seu_rep_qc, nfeatures = nfeatures_to_find, verbose = FALSE)
  
  message(paste0("        DEBUG: Scaling Data (vars.to.regress=", ifelse(is.null(vars_to_regress_param), "NULL", vars_to_regress_param), ")..."))
  seu_rep_qc <- ScaleData(seu_rep_qc, vars.to.regress = vars_to_regress_param, verbose = FALSE, features = VariableFeatures(seu_rep_qc))
  
  npcs_to_run <- doublet_params$npcs_align %||% 30
  message(paste0("        DEBUG: Running PCA (npcs=", npcs_to_run, ", on assay: RNA)..."))
  tryCatch({
    seu_rep_qc <- RunPCA(seu_rep_qc, npcs = npcs_to_run, verbose = FALSE, assay = "RNA", features = VariableFeatures(seu_rep_qc))
  }, error = function(e_pca){
    message(paste0("        ERROR during RunPCA for ", rep_name, " on assay RNA: ", e_pca$message))
    stop(paste0("RunPCA failed for ", rep_name, ", cannot run DoubletFinder."))
  })
  # --- 预处理结束 ---

  # --- Sanitize Seurat Object (Post-PCA, Pre-paramSweep/DoubletFinder) ---
  message("        DEBUG: Sanitizing Seurat object globally before paramSweep & DoubletFinder...")

  # 1. Sanitize PCA embeddings (the one we'll keep)
  if ("pca" %in% Reductions(seu_rep_qc)) {
    message("        DEBUG: Sanitizing PCA cell embeddings...")
    current_pca_reduc <- seu_rep_qc@reductions$pca
    embeds <- current_pca_reduc@cell.embeddings
    
    if (!is.matrix(embeds) || !is.numeric(embeds)) {
      message("        WARNING: PCA embeddings not a numeric matrix. Attempting conversion.")
      embeds <- as.matrix(embeds)
      if(!is.numeric(embeds)) stop("PCA embeddings cannot be coerced to numeric matrix for ", rep_name)
    }
    # Ensure standard colnames like PC_1, PC_2 using the key from the reduction
    pca_key <- Key(current_pca_reduc)
    colnames(embeds) <- paste0(pca_key, 1:ncol(embeds))
    seu_rep_qc@reductions$pca@cell.embeddings <- embeds
    message("        DEBUG: PCA cell embeddings sanitized.")
  } else {
    stop(paste0("PCA reduction not found after RunPCA for ", rep_name, ", cannot proceed."))
  }

  # 2. Set DefaultAssay and remove others
  assay_to_keep <- DefaultAssay(seu_rep_qc) # Should be SCT if SCTransform succeeded
  message(paste0("        DEBUG: Current DefaultAssay to keep is: ", assay_to_keep))
  assays_present <- Assays(seu_rep_qc)
  for (assay_name in assays_present) {
    if (assay_name != assay_to_keep) {
      message(paste0("        DEBUG: Removing assay: ", assay_name))
      seu_rep_qc[[assay_name]] <- NULL
    }
  }

  # 3. Remove other reductions
  reductions_present <- Reductions(seu_rep_qc)
  for (reduc_name in reductions_present) {
    if (reduc_name != "pca") {
      message(paste0("        DEBUG: Removing reduction: ", reduc_name))
      seu_rep_qc@reductions[[reduc_name]] <- NULL
    }
  }

  # 4. Deep sanitize meta.data
  message("        DEBUG: Performing deep sanitization of meta.data...")
  if (!is.null(seu_rep_qc) && "meta.data" %in% slotNames(seu_rep_qc) && ncol(seu_rep_qc@meta.data) > 0) {
    md <- seu_rep_qc@meta.data
    original_md_rownames <- rownames(md)
    md_df <- as.data.frame(md)
    for (col_name in colnames(md_df)) {
      if (is.list(md_df[[col_name]]) || !is.atomic(md_df[[col_name]])) {
        message(paste0("        DEBUG: Column '", col_name, "' in meta.data is complex. Converting to character by pasting elements."))
        md_df[[col_name]] <- sapply(md_df[[col_name]], function(x) paste(as.character(x), collapse=";"))
      } else if (is.factor(md_df[[col_name]])) {
        message(paste0("        DEBUG: Column '", col_name, "' in meta.data is factor. Converting to character."))
        md_df[[col_name]] <- as.character(md_df[[col_name]])
      }
    }
    seu_rep_qc@meta.data <- md_df
    if (!is.null(original_md_rownames) && length(original_md_rownames) == nrow(seu_rep_qc@meta.data)) {
        rownames(seu_rep_qc@meta.data) <- original_md_rownames
    } else {
        message("        WARNING: Could not restore rownames to meta.data after deep sanitization for ", rep_name)
    }
    message("        DEBUG: Deep meta.data sanitization complete.")
  } else {
    message("        DEBUG: meta.data is NULL or has no columns. Skipping deep sanitization for ", rep_name)
  }
  message("        DEBUG: Global Seurat object sanitization complete for ", rep_name, ".")
  # --- End Global Sanitization ---

  # --- DIAGNOSTIC: Direct DoubletFinder Call ---
  message("        DEBUG: --- Starting Diagnostic Direct DoubletFinder Call ---")
  pK_test_diagnostic <- doublet_params$pK_fallback_align %||% 0.09 # Use fallback pK or a common value
  nExp_poi_diagnostic <- round(ncol(seu_rep_qc) * (doublet_params$expected_doublet_rate_align %||% 0.009))
  pN_diagnostic <- doublet_params$pN_align %||% 0.25
  pcs_for_diagnostic <- 1:(doublet_params$pcs_align %||% 30)
  sct_diagnostic <- FALSE # Explicitly FALSE for this diagnostic path

  message(paste0("        DEBUG_DIAGNOSTIC: Calling doubletFinder directly with pN=", pN_diagnostic, 
                 ", pK=", pK_test_diagnostic, ", PCs=1:", max(pcs_for_diagnostic), 
                 ", nExp=", nExp_poi_diagnostic, ", sct=", sct_diagnostic))
  
  seu_obj_copy_for_diagnostic <- seu_rep_qc # Work on a copy to not affect the main flow yet
  diagnostic_call_successful <- FALSE
  tryCatch({
    seu_obj_copy_for_diagnostic <- doubletFinder(
      seu = seu_obj_copy_for_diagnostic,
      PCs = pcs_for_diagnostic, 
      pN = pN_diagnostic,
      pK = pK_test_diagnostic,   
      nExp = nExp_poi_diagnostic,
      reuse.pANN = FALSE, 
      sct = sct_diagnostic
    )
    diagnostic_call_successful <- TRUE
    message("        DEBUG_DIAGNOSTIC: Direct doubletFinder call SUCCEEDED.")
    # Optional: check for DF classification columns
    df_cols_diag <- colnames(seu_obj_copy_for_diagnostic@meta.data)[grep(paste0("^(DF.classifications|pANN)_", pN_diagnostic, "_", pK_test_diagnostic), colnames(seu_obj_copy_for_diagnostic@meta.data))]
    message(paste0("        DEBUG_DIAGNOSTIC: DoubletFinder output columns found: ", paste(df_cols_diag, collapse=", ")))

  }, error = function(e_diag) {
    message(paste0("        ERROR_DIAGNOSTIC: Direct doubletFinder call FAILED: ", e_diag$message))
  })
  message("        DEBUG: --- Finished Diagnostic Direct DoubletFinder Call (Success: ", diagnostic_call_successful, ") ---")
  # --- END DIAGNOSTIC --- 

  # --- DoubletFinder参数优化 (paramSweep) --- 
  pcs_for_sweep_and_df <- 1:(doublet_params$pcs_align %||% 30)
  message(paste0("        DEBUG: PCs for Sweep & DF: 1 to ", max(pcs_for_sweep_and_df)))
  
  sweep_file <- file.path(sweep_dir, paste0(rep_name, "_sweep_results.rds"))
  optimal_pK <- NULL
  
  if (!(doublet_params$recalculate_sweep %||% FALSE) && file.exists(sweep_file)) {
    message(paste0("        Loading existing sweep results for ", rep_name, " from: ", sweep_file))
    sweep.res_list <- readRDS(sweep_file)
    # Ensure sweep.res_list is not NULL or empty if loaded
    if (!is.null(sweep.res_list) && length(sweep.res_list) > 0) {
        sweep.stats <- summarizeSweep(sweep.res_list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
        if (nrow(bcmvn) > 0 && "BCmetric" %in% colnames(bcmvn)) {
    optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message(paste0("        Optimal pK from loaded sweep: ", optimal_pK))
  } else {
             message("        Warning: BCmvn from loaded sweep did not yield optimal pK. Will run paramSweep or use default.")
        }
    } else {
        message("        Warning: Loaded sweep file was empty or NULL. Will run paramSweep or use default pK.")
    }
  }
  
  # If optimal_pK still NULL (not foundPrince or recalculate_sweep is TRUE)
  if (is.null(optimal_pK) || (doublet_params$recalculate_sweep %||% FALSE) ) {
    message(paste0("        Running paramSweep for ", rep_name, " (PCs: 1:", max(pcs_for_sweep_and_df), ")"))
    sweep.res_list <- paramSweep(seu_rep_qc, PCs = pcs_for_sweep_and_df, sct = doublet_params$sct %||% FALSE, num.cores = max(1, parallel::detectCores() - 2))
    if (!is.null(sweep.res_list) && length(sweep.res_list) > 0) {
      tryCatch({
        saveRDS(sweep.res_list, file = sweep_file)
        message(paste0("        Sweep results saved to: ", sweep_file))
      }, error = function(e_save) {
        message(paste0("        Warning: Failed to save sweep results for ", rep_name, ". Error: ", e_save$message))
      })
      sweep.stats <- summarizeSweep(sweep.res_list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      if (nrow(bcmvn) > 0 && "BCmetric" %in% colnames(bcmvn)) {
        optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        message(paste0("        Optimal pK from new sweep: ", optimal_pK))
        
        pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_point() + geom_line() +
          geom_vline(xintercept = optimal_pK, linetype="dashed", color="red") +
                     labs(title = paste0("pK Optimization - ", rep_name), x = "pK", y = "BCmetric") + theme_bw()
        save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_pk_optimization.jpeg")), pk_plot, width = 6, height = 4)
      } else {
        message("        Warning: BCmvn from new sweep did not yield optimal pK. Using default pK.")
      }
    } else {
      message("        Warning: paramSweep returned NULL or empty. Using default pK.")
    }
  }
  
  # Fallback pK if still NULL
  if (is.null(optimal_pK) || is.na(optimal_pK) || length(optimal_pK) == 0) {
    optimal_pK <- doublet_params$pK_fallback_align %||% 0.09 # Configurable fallback pK
    message(paste0("        Using fallback pK: ", optimal_pK, " for ", rep_name))
  }
  
  # --- 计算预期双细胞数 (nExp) ---
  # modelHomotypic part from example - requires clusters. We don't have them here yet.
  # For now, keep the original nExp calculation. If SCT works, we can refine this later.
  doublet_rate_to_use <- doublet_params$expected_doublet_rate_align %||% 0.009 
  nExp_poi <- round(ncells * doublet_rate_to_use)
  message(paste0("        Expected doublets (nExp_poi, rate=", doublet_rate_to_use*100, "%): ", nExp_poi, " for ", rep_name))
  
  # --- 执行DoubletFinder --- 
  pN_to_use <- doublet_params$pN_align %||% 0.25 # Configurable, default 0.25
  # Ensure sct_to_use is FALSE for the main path now
  sct_to_use <- FALSE # DefaultAssay will be RNA, and config$sct is FALSE
  message(paste0("        Running doubletFinder with pN=", pN_to_use, ", pK=", optimal_pK, ", PCs=1:", max(pcs_for_sweep_and_df), ", nExp=", nExp_poi, ", sct=", sct_to_use))
  
  # DEBUG 信息
  message(paste0("        DEBUG: Class of seu_rep_qc: ", class(seu_rep_qc)))
  if ("pca" %in% names(seu_rep_qc@reductions)) {
    message(paste0("        DEBUG: Dimensions of PCA embedding: ", paste(dim(seu_rep_qc@reductions$pca@cell.embeddings), collapse=" x ")))
  }
  message(paste0("        DEBUG: PCs to use for DF: 1:", max(pcs_for_sweep_and_df)))
  
  tryCatch({
    seu_rep_qc <- doubletFinder(
      seu = seu_rep_qc,
      PCs = pcs_for_sweep_and_df, 
      pN = pN_to_use,
      pK = optimal_pK,
      nExp = nExp_poi,
      reuse.pANN = FALSE, 
      sct = sct_to_use # This will now be correctly FALSE
    )
    
    df_score_col <- colnames(seu_rep_qc@meta.data)[grep(paste0("pANN_",pN_to_use,"_", optimal_pK, "_"), colnames(seu_rep_qc@meta.data))]
    df_class_col <- colnames(seu_rep_qc@meta.data)[grep(paste0("DF.classifications_",pN_to_use,"_", optimal_pK, "_"), colnames(seu_rep_qc@meta.data))]
    
    if (length(df_score_col) == 1 && length(df_class_col) == 1) {
      seu_rep_qc$doublet_score <- seu_rep_qc@meta.data[[df_score_col]]
      seu_rep_qc$doublet_class <- seu_rep_qc@meta.data[[df_class_col]]
      
      doublet_count <- sum(seu_rep_qc$doublet_class == "Doublet", na.rm = TRUE)
      message(paste0("        Step 3.2: Detected ", doublet_count, " doublets (", 
                     round(doublet_count/ncells*100, 2), "%) for ", rep_name))
      
      if ("umap" %in% names(seu_rep_qc@reductions)) {
        plot_df_class <- DimPlot(seu_rep_qc, reduction = "umap", group.by = "doublet_class", cols = c("Singlet" = "grey", "Doublet" = "red", "Singlet_low_cell_count"="blue", "Singlet_df_exception"="orange", "Singlet_df_error"="purple")) +
          labs(title = paste0("Doublet Classification - ", rep_name))
        save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_doublet_class_umap.jpeg")), plot_df_class, width = 7, height = 6)
        
        plot_df_score <- FeaturePlot(seu_rep_qc, features = "doublet_score", reduction = "umap", order = TRUE, pt.size = 0.5) +
          scale_color_viridis_c() +
          labs(title = paste0("Doublet Score - ", rep_name))
        save_plot_safe(file.path(rep_plot_dir, paste0(rep_name, "_doublet_score_umap.jpeg")), plot_df_score, width = 7, height = 6)
      }
      
    } else {
      message(paste0("        Warning: DoubletFinder output columns not found as expected for ", rep_name, ". Assigning all as Singlets. Scores found: ", length(df_score_col), ", Classes found: ", length(df_class_col)))
      seu_rep_qc$doublet_score <- 0
      seu_rep_qc$doublet_class <- "Singlet_df_error"
    }
  }, error = function(e) {
    message(paste0("        ERROR during DoubletFinder for ", rep_name, ": ", e$message, ". Assigning all as Singlets."))
    seu_rep_qc$doublet_score <- 0
    seu_rep_qc$doublet_class <- "Singlet_df_exception"
  })
  
  return(seu_rep_qc)
}

#' @title 过滤双细胞 (优化版)
filter_doublets_optimized <- function(seu_rep_doublet, rep_name) {
  message(paste0("      Step 4.1: 过滤双细胞 for ", rep_name))
  cells_before_removal <- ncol(seu_rep_doublet)
  
  if (!"doublet_class" %in% colnames(seu_rep_doublet@meta.data)) {
    message(paste0("        Warning: 'doublet_class' column not found in metadata for ", rep_name, ". Skipping doublet filtering."))
    return(seu_rep_doublet)
  }
  
  if (!is.factor(seu_rep_doublet$doublet_class)) {
    seu_rep_doublet$doublet_class <- factor(seu_rep_doublet$doublet_class)
  }
  expected_levels <- c("Singlet", "Doublet", "Singlet_low_cell_count", "Singlet_df_exception", "Singlet_df_error")
  current_values <- unique(as.character(seu_rep_doublet$doublet_class))
  all_expected_levels <- unique(c(levels(seu_rep_doublet$doublet_class), expected_levels, current_values))
  seu_rep_doublet$doublet_class <- factor(seu_rep_doublet$doublet_class, levels = all_expected_levels)

  singlet_categories <- c("Singlet", "Singlet_low_cell_count", "Singlet_df_exception", "Singlet_df_error")
  seu_clean <- subset(seu_rep_doublet, subset = doublet_class %in% singlet_categories)
  cells_after_removal <- ncol(seu_clean)
  
  if (cells_after_removal == 0 && cells_before_removal > 0) {
    message(paste0("        Warning: All cells were filtered out as doublets for ", rep_name, ". Returning the object before doublet filtering to avoid errors."))
    return(seu_rep_doublet) 
  }
  
  message(paste0("      Step 4.2: 双细胞过滤完成. Before=", cells_before_removal, ", After=", cells_after_removal, " for ", rep_name))
  return(seu_clean)
}


# --- 3. 主模块函数 ---
#' @title 执行QC和双细胞去除 (优化版)
#' @description 对原始Seurat对象列表执行QC和双细胞检测/去除。
#' @param seu_raw_list 原始Seurat对象列表 (结构: sample_id -> rep_id -> seu_object)。
#' @param results_dir 主结果输出目录。
#' @param qc_params QC参数列表 (来自config$module1$qc_params)。
#' @param doublet_params DoubletFinder参数列表 (来自config$module1$doublet_params)。
#' @param utils_path utils.R的路径。
#' @return 处理后的Seurat对象列表。
run_qc_and_doublet_removal <- function(seu_raw_list, results_dir, qc_params, doublet_params, utils_path = "utils.R") {
  
  module_start_time <- Sys.time()
  message(paste0("\n=== 开始执行 Module 1: QC 与双细胞过滤 (优化版) (", module_start_time, ") ==="))
  
  # 创建此模块的输出子目录
  module1_output_dir <- file.path(results_dir, "Module1_QC_Doublet")
  if (!dir.exists(module1_output_dir)) dir.create(module1_output_dir, recursive = TRUE)
  
  plots_output_dir <- file.path(module1_output_dir, "Plots")
  if (!dir.exists(plots_output_dir)) dir.create(plots_output_dir, recursive = TRUE)
  
  sweep_output_dir <- file.path(module1_output_dir, "DoubletFinder_SweepData")
  if (!dir.exists(sweep_output_dir)) dir.create(sweep_output_dir, recursive = TRUE)
  
  # 加载utils.R中的辅助函数 (如果需要)
  if (file.exists(utils_path)) {
    source(utils_path, local = TRUE) # local=TRUE 避免污染全局环境
    message(paste("    辅助函数脚本 '", utils_path, "' 已加载到模块环境中。"))
  } else {
    warning(paste("    警告: 辅助函数脚本 '", utils_path, "' 未找到。某些功能可能受限。"))
  }
  
  # 初始化存储处理后Seurat对象的列表
  seu_processed_list <- list()
  
  # 遍历每个样品组
  for (sample_id in names(seu_raw_list)) {
    message(paste0("\n  处理样品组: ", sample_id))
    seu_processed_list[[sample_id]] <- list()
    
    sample_replicates_raw <- seu_raw_list[[sample_id]]
    if (is.null(sample_replicates_raw) || length(sample_replicates_raw) == 0) {
      message(paste0("    样品组 ", sample_id, " 不包含任何原始复制数据，跳过。"))
      next
    }
    
    # 遍历样品组中的每个复制
    for (rep_id in names(sample_replicates_raw)) {
      seu_rep_raw <- sample_replicates_raw[[rep_id]]
      full_rep_name <- paste0(sample_id, "_", rep_id) # 例如 "BC-0h_rep1"
      message(paste0("    处理复制: ", full_rep_name))
      
      # 为当前复制创建特定的绘图子目录
      rep_plot_dir_specific <- file.path(plots_output_dir, sample_id, rep_id)
      if (!dir.exists(rep_plot_dir_specific)) dir.create(rep_plot_dir_specific, recursive = TRUE)
      
      if (is.null(seu_rep_raw) || !inherits(seu_rep_raw, "Seurat")) {
        message(paste0("      复制 ", full_rep_name, " 的原始数据为NULL或不是Seurat对象，跳过。"))
        seu_processed_list[[sample_id]][[rep_id]] <- NULL
        next
      }
      if (ncol(seu_rep_raw) == 0) {
        message(paste0("      复制 ", full_rep_name, " 原始数据中细胞数为0，跳过。"))
        seu_processed_list[[sample_id]][[rep_id]] <- NULL
        next
      }
      
      seu_rep_current <- seu_rep_raw # 当前处理的Seurat对象
      
      tryCatch({
        # --- Step 1: (数据加载已在主脚本完成) ---
        # 添加原始样本和复制信息到元数据 (如果主脚本未添加)
        if (!"original_sample_id" %in% colnames(seu_rep_current@meta.data)) seu_rep_current$original_sample_id <- sample_id
        if (!"replicate_id" %in% colnames(seu_rep_current@meta.data)) seu_rep_current$replicate_id <- rep_id
        
        # --- Step 2: 质量控制 ---
        message(paste0("      Starting QC for ", full_rep_name))
        seu_rep_qc <- calculate_and_apply_qc_optimized(
          seu_rep = seu_rep_current, 
          qc_params = qc_params, 
          rep_name = full_rep_name,
          rep_plot_dir = rep_plot_dir_specific,
          utils_path = utils_path
        )
        
        # 检查QC后细胞数
        if (is.null(seu_rep_qc) || ncol(seu_rep_qc) == 0 || (!is.null(seu_rep_qc$qc_passed) && !all(seu_rep_qc$qc_passed))) {
          message(paste0("      QC失败或导致0细胞 for ", full_rep_name, ". 跳过此复制的后续步骤。"))
          seu_processed_list[[sample_id]][[rep_id]] <- NULL # 标记为处理失败
          next # 跳到下一个复制
        }
        seu_rep_current <- seu_rep_qc
        
        # --- Step 3: 双细胞检测 ---
        if (!is.null(doublet_params$expected_doublet_rate) && doublet_params$expected_doublet_rate > 0) {
          message(paste0("      Starting Doublet Detection for ", full_rep_name))
          seu_rep_doublet_info <- run_doublet_detection_optimized(
            seu_rep_qc = seu_rep_current,
            doublet_params = doublet_params,
            rep_name = full_rep_name,
            rep_plot_dir = rep_plot_dir_specific,
            sweep_dir = sweep_output_dir, # 传递sweep结果保存目录
            utils_path = utils_path
          )
          seu_rep_current <- seu_rep_doublet_info
          
          # --- Step 4: 双细胞过滤 ---
          message(paste0("      Starting Doublet Filtering for ", full_rep_name))
          seu_rep_final <- filter_doublets_optimized(
            seu_rep_doublet = seu_rep_current,
            rep_name = full_rep_name
          )
        } else {
          message(paste0("      Skipping Doublet Detection and Filtering for ", full_rep_name, " as expected_doublet_rate is <= 0 or not set."))
          # 如果跳过双细胞检测，确保 doublet_class 列存在且为 "Singlet"
          if (!"doublet_class" %in% colnames(seu_rep_current@meta.data)) {
            seu_rep_current$doublet_class <- "Singlet_not_run"
          }
          if (!"doublet_score" %in% colnames(seu_rep_current@meta.data)) {
            seu_rep_current$doublet_score <- 0
          }
          seu_rep_final <- seu_rep_current
        }
        
        # 存储最终处理好的对象 (如果还有细胞)
        if (ncol(seu_rep_final) > 0) {
          seu_processed_list[[sample_id]][[rep_id]] <- seu_rep_final
          message(paste0("    复制 ", full_rep_name, " 处理成功，最终细胞数: ", ncol(seu_rep_final)))
        } else {
          message(paste0("    警告: 复制 ", full_rep_name, " 在处理后细胞数为0。此复制将被跳过。"))
          seu_processed_list[[sample_id]][[rep_id]] <- NULL
        }
        
      }, error = function(e) {
        message(paste0("    错误: 处理复制 ", full_rep_name, " 时发生严重错误: ", e$message))
        # 记录错误，并将此复制的结果设为NULL
        seu_processed_list[[sample_id]][[rep_id]] <- NULL
      })
      
      # 内存管理
      if (exists("seu_rep_raw")) rm(seu_rep_raw)
      if (exists("seu_rep_qc")) rm(seu_rep_qc)
      if (exists("seu_rep_doublet_info")) rm(seu_rep_doublet_info)
      if (exists("seu_rep_final")) rm(seu_rep_final)
      if (exists("seu_rep_current") && !is.null(seu_processed_list[[sample_id]][[rep_id]])) {
        # 如果成功，seu_rep_current 可能是 seu_rep_final 的引用，不需要删除
      } else if (exists("seu_rep_current")) {
        rm(seu_rep_current)
      }
      gc(verbose = FALSE) # 在每个复制后进行垃圾回收
      
    } # 结束复制循环
    
    # 清理当前样品组中失败的复制 (NULL元素)
    if (length(seu_processed_list[[sample_id]]) > 0) {
      seu_processed_list[[sample_id]] <- Filter(Negate(is.null), seu_processed_list[[sample_id]])
      if (length(seu_processed_list[[sample_id]]) == 0) { # 如果所有复制都失败了
        message(paste0("  样品组 ", sample_id, " 的所有复制处理失败或产生0细胞。"))
        seu_processed_list[[sample_id]] <- NULL # 将整个样品组设为NULL
      } else {
        message(paste0("  样品组 ", sample_id, " 处理完成，包含 ", length(seu_processed_list[[sample_id]]), " 个有效复制。"))
      }
    } else { # 如果最初就没有复制数据或都跳过了
      seu_processed_list[[sample_id]] <- NULL
    }
    
  } # 结束样品组循环
  
  # 清理掉完全失败的样品组 (那些值为NULL的)
  final_seu_clean_list <- Filter(Negate(is.null), seu_processed_list)
  
  if (length(final_seu_clean_list) == 0) {
    warning("Module 1 警告: 所有样品组处理失败或未产生任何有效数据。返回空列表。")
  } else {
    message(paste0("\nModule 1 有效处理了 ", length(final_seu_clean_list), " 个样品组。"))
  }
  
  module_end_time <- Sys.time()
  message(paste0("=== Module 1: QC 与双细胞过滤 (优化版) 完成 (总耗时: ", 
                 format(difftime(module_end_time, module_start_time, units = "auto")), ") ==="))
  
  # 尝试调用 utils.R 中的 memory_usage (如果已加载)
  if (exists("memory_usage") && is.function(memory_usage)) {
    try(memory_usage(), silent = TRUE)
  }
  
  return(final_seu_clean_list) 
}

# 示例用法 (在主脚本中调用):
# source("module1_qc_doublet_optimized.R") # 确保此行在主脚本中
# seu_clean_list_m1 <- run_qc_and_doublet_removal(
#   seu_raw_list = seu_raw_list_from_main_script, # 这是主脚本加载的原始数据列表
#   results_dir = config$paths$results_dir,      # 来自config对象
#   qc_params = config$module1$qc_params,        # 来自config对象
#   doublet_params = config$module1$doublet_params, # 来自config对象
#   utils_path = "utils.R" # 确保utils.R的路径正确
# )
# # 之后可以保存 seu_clean_list_m1
# saveRDS(seu_clean_list_m1, file.path(config$paths$results_dir, "Module1_QC_Doublet", "seu_clean_list_m1_output.rds"))

