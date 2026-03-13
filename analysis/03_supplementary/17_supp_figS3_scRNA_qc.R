# =========================================================================
# Script: 17_supp_figS3_scRNA_qc.R
# Purpose: Generate Supplementary Figure S3 (QC + batch correction diagnostics)
# Author: Zuo
# Date: 2026-03-05
#
# Inputs:
#   - data/processed_full/single_cell/gse274583/gse274583_cd4_integrated_reviewer_diet.rds
#   - R/00_theme_setup.R   # Theme & palette (optional; will fall back to theme_bw)
#   - R/03_index_logger.R  # Logging helper (optional)
#
# Outputs (to output/supplement):
#   - FigS3A_scRNA_qc_by_batch1234.pdf
#   - FigS3B_harmony_umap_by_batch1234.pdf
#   - output/logs/17_supp_figS3_scRNA_qc.log
#   - output/logs/17_supp_figS3_scRNA_qc_sessionInfo.txt
#   - output/logs/17_supp_figS3_scRNA_qc_audit.tsv
# =========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# Project root + path setup
# ----------------------------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop(
    "Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).",
    call. = FALSE
  )
}

DATA_PROCESSED_FULL_DIR <- file.path(base_dir, "data", "processed_full")
SUPP_DIR               <- file.path(base_dir, "output", "supplement")
LOG_DIR                <- file.path(base_dir, "output", "logs")
THEME_FILE             <- file.path(base_dir, "R", "00_theme_setup.R")
LOGGER_FILE            <- file.path(base_dir, "R", "03_index_logger.R")

SC_REVIEWER_RDS <- file.path(
  DATA_PROCESSED_FULL_DIR,
  "single_cell", "gse274583",
  "gse274583_cd4_integrated_reviewer_diet.rds"
)

dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE         <- file.path(LOG_DIR, "17_supp_figS3_scRNA_qc.log")
SESSIONINFO_FILE <- file.path(LOG_DIR, "17_supp_figS3_scRNA_qc_sessionInfo.txt")

# Start logging
sink(LOG_FILE, append = FALSE, split = TRUE)
cat("=== 17_supp_figS3_scRNA_qc.R ===\n")
cat("Base dir: ", base_dir, "\n", sep = "")
cat("Input RDS: ", SC_REVIEWER_RDS, "\n", sep = "")
cat("Outputs:\n")
cat("  - ", file.path(SUPP_DIR, "FigS3A_scRNA_qc_by_batch1234.pdf"), "\n", sep = "")
cat("  - ", file.path(SUPP_DIR, "FigS3B_harmony_umap_by_batch1234.pdf"), "\n", sep = "")
cat("Log: ", LOG_FILE, "\n\n", sep = "")

# ----------------------------
# Theme setup (optional)
# ----------------------------
if (file.exists(THEME_FILE)) {
  source(THEME_FILE)
}
if (exists("set_plot_theme", mode = "function")) {
  set_plot_theme()
}

# Logger helper (optional)
if (file.exists(LOGGER_FILE)) {
  source(LOGGER_FILE)
}

# ----------------------------
# Load object
# ----------------------------
if (!file.exists(SC_REVIEWER_RDS)) {
  stop("Missing reviewer RDS: ", SC_REVIEWER_RDS, call. = FALSE)
}

obj <- readRDS(SC_REVIEWER_RDS)

# ----------------------------
# Ensure batch1234 exists (DO NOT use library_id for plotting)
# ----------------------------
valid_batches <- paste0("Batch", 1:4)

# Priority: use existing batch1234 if present
if ("batch1234" %in% colnames(obj@meta.data)) {
  cat("Using existing batch1234 field from object\n")
  
  # Validate levels
  bad_levels <- setdiff(unique(obj$batch1234), valid_batches)
  if (length(bad_levels) > 0) {
    stop(
      "Unexpected batch1234 levels found: ",
      paste(bad_levels, collapse = ", "),
      "\nValid levels are: ", paste(valid_batches, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Ensure no NA values
  if (any(is.na(obj$batch1234))) {
    stop("batch1234 contains NA values", call. = FALSE)
  }
  
  # Set factor levels in fixed order
  obj$batch1234 <- factor(obj$batch1234, levels = valid_batches)
  
} else {
  # Fallback: derive from library_id
  cat("Deriving batch1234 from library_id\n")
  
  if (!"library_id" %in% colnames(obj@meta.data)) {
    stop("Missing meta column: library_id (needed to derive batch1234)", call. = FALSE)
  }
  
  # Map library_id -> Batch1/2/3/4
  obj$batch1234 <- dplyr::case_when(
    obj$library_id %in% c("Library1") ~ "Batch1",
    obj$library_id %in% c("Library2") ~ "Batch2",
    obj$library_id %in% c("Library3") ~ "Batch3",
    grepl("^Library4", obj$library_id) ~ "Batch4",
    TRUE ~ NA_character_
  )
  
  # Hard guardrail: only allow Batch1-4
  bad_levels <- setdiff(unique(obj$batch1234), valid_batches)
  if (length(bad_levels) > 0) {
    stop(
      "Unexpected batch1234 levels found: ",
      paste(bad_levels, collapse = ", "),
      "\nCheck library_id values in meta.data.",
      call. = FALSE
    )
  }
  
  # Ensure no NA values
  if (any(is.na(obj$batch1234))) {
    stop("batch1234 contains NA values after mapping from library_id", call. = FALSE)
  }
  
  # Set factor levels in fixed order
  obj$batch1234 <- factor(obj$batch1234, levels = valid_batches)
}

# ----------------------------
# Audit output
# ----------------------------
AUDIT_FILE <- file.path(LOG_DIR, "17_supp_figS3_scRNA_qc_audit.tsv")
cat("\n=== Audit Tables ===\n")
cat("Batch1234 distribution:\n")
batch_dist <- table(obj$batch1234)
print(batch_dist)

cat("\nBatch1234 vs library_id:\n")
batch_lib_table <- table(obj$batch1234, obj$library_id)
print(batch_lib_table)

cat("\nBatch1234 vs seurat_clusters:\n")
if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  batch_cluster_table <- table(obj$batch1234, obj$seurat_clusters)
  print(batch_cluster_table)
} else {
  cat("seurat_clusters not found in meta.data\n")
}

# Write audit file
sink(AUDIT_FILE)
cat("=== 17_supp_figS3_scRNA_qc Audit ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Batch1234 distribution:\n")
print(batch_dist)

cat("\nBatch1234 vs library_id:\n")
print(batch_lib_table)

if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  cat("\nBatch1234 vs seurat_clusters:\n")
  print(batch_cluster_table)
}

sink()

# percent.mt: compute if missing
if (!"percent.mt" %in% colnames(obj@meta.data)) {
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  cat("Computed percent.mt\n")
}

# ----------------------------
# Helper: sampling for plotting (50k) with stratified sampling by batch1234
# ----------------------------
set.seed(2026)
sample_cells <- function(seu, max_cells = 50000, strata = "batch1234") {
  n <- ncol(seu)
  if (n <= max_cells) return(colnames(seu))
  
  # Stratified sampling by batch1234
  if (!is.null(strata) && strata %in% colnames(seu@meta.data)) {
    # Calculate batch proportions
    batch_props <- prop.table(table(seu@meta.data[[strata]]))
    
    # Calculate number of cells per batch
    batch_counts <- round(batch_props * max_cells)
    
    # Ensure total is exactly max_cells
    batch_counts[1] <- batch_counts[1] + (max_cells - sum(batch_counts))
    
    # Sample cells from each batch
    sampled_cells <- c()
    for (batch in names(batch_counts)) {
      batch_cells <- colnames(seu)[seu@meta.data[[strata]] == batch]
      n_batch <- min(batch_counts[batch], length(batch_cells))
      sampled_cells <- c(sampled_cells, sample(batch_cells, size = n_batch, replace = FALSE))
    }
    
    return(sampled_cells)
  } else {
    # Fallback to simple random sampling
    return(sample(colnames(seu), size = max_cells, replace = FALSE))
  }
}

cells_50k <- sample_cells(obj, 50000, strata = "batch1234")
obj_50k <- subset(obj, cells = cells_50k)

# Log batch proportions in sample
cat("\nBatch proportions in full object:\n")
print(prop.table(table(obj$batch1234)))
cat("Batch proportions in 50k sample:\n")
print(prop.table(table(obj_50k$batch1234)))

# ----------------------------
# Fig S3A: QC metrics by batch1234
# ----------------------------
cat("\n[FigS3A] QC plots by batch1234...\n")

# Reduce extreme compression for nCount by trimming high outliers for plotting only
# (keeps distributions readable; does not change object saved elsewhere)
cat("QC plot subset: nCount_RNA < 50000 (plot-only)\n")
obj_qc_plot <- subset(obj, subset = nCount_RNA < 50000)

# 分别创建三个QC指标的图，统一y轴范围
p_nFeature <- VlnPlot(
  obj_qc_plot,
  features = "nFeature_RNA",
  group.by = "batch1234",
  pt.size = 0,
  cols = trae_colors_celltype  # 使用主题配色
) +
  theme_allergy() +  # 使用主题配色
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),  # 去掉x轴标题
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  ylim(0, 4000)  # 统一nFeature_RNA的y轴上限

p_nCount <- VlnPlot(
  obj_qc_plot,
  features = "nCount_RNA",
  group.by = "batch1234",
  pt.size = 0,
  cols = trae_colors_celltype  # 使用主题配色
) +
  theme_allergy() +  # 使用主题配色
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),  # 去掉x轴标题
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  ylim(0, 50000)  # 统一nCount_RNA的y轴上限

p_percent_mt <- VlnPlot(
  obj_qc_plot,
  features = "percent.mt",
  group.by = "batch1234",
  pt.size = 0,
  cols = trae_colors_celltype  # 使用主题配色
) +
  theme_allergy() +  # 使用主题配色
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),  # 去掉x轴标题
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  ylim(0, 15)  # 统一percent.mt的y轴上限

# 组合三个图，添加共用图例
p_qc <- p_nFeature + p_nCount + p_percent_mt +
  plot_layout(ncol = 3, guides = "collect") +  # 收集图例
  theme(legend.position = "bottom")  # 图例放在底部

figS3A_out <- file.path(SUPP_DIR, "FigS3A_scRNA_qc_by_batch1234.pdf")
ggsave(figS3A_out, p_qc, width = 11, height = 4.2, units = "in")
cat("Saved: ", figS3A_out, "\n", sep = "")

if (exists("log_output", mode = "function")) {
  log_output(
    file_path = figS3A_out,
    file_type = "figure",
    figure_id = "Fig S3A",
    description = "scRNA QC metrics (nFeature, nCount, percent.mt) by batch1234",
    script_name = "17_supp_figS3_scRNA_qc.R",
    resolution = "vector pdf"
  )
}

# ----------------------------
# Fig S3B: Batch correction diagnostics (before/after) by batch1234
# ----------------------------
cat("\n[FigS3B] Harmony diagnostics UMAP by batch1234...\n")

if (!"pca" %in% names(obj@reductions)) {
  warning("PCA reduction not found. Will plot Harmony UMAP only.")
}
if (!"umap" %in% names(obj@reductions)) {
  stop("UMAP reduction not found in object. Ensure integration script saved umap.", call. = FALSE)
}

# Panel 1: PCA colored by batch (if available)
p_pca <- NULL
if ("pca" %in% names(obj_50k@reductions)) {
  p_pca <- DimPlot(
    obj_50k,
    reduction = "pca",
    group.by = "batch1234",
    raster = FALSE,  # 与主文Fig3保持一致
    pt.size = 0.2,  # 与主文Fig3保持一致
    cols = trae_colors_celltype  # 使用主题配色
  ) +
    theme_allergy() +  # 使用主题配色
    ggtitle("Before Harmony (PCA) by batch1234") +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold")
    )
}

# Panel 2: UMAP (Harmony space) colored by batch
p_umap_batch <- DimPlot(
  obj_50k,
  reduction = "umap",
  group.by = "batch1234",
  raster = FALSE,  # 与主文Fig3保持一致
  pt.size = 0.2,  # 与主文Fig3保持一致
  cols = trae_colors_celltype  # 使用主题配色
) +
  theme_allergy() +  # 使用主题配色
  ggtitle("After Harmony (UMAP based on Harmony embeddings) by batch1234") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  )

# Combine
p_harmony_diag <- if (!is.null(p_pca)) {
  p_pca + p_umap_batch + plot_layout(ncol = 2)
} else {
  p_umap_batch
}

figS3B_out <- file.path(SUPP_DIR, "FigS3B_harmony_umap_by_batch1234.pdf")
ggsave(figS3B_out, p_harmony_diag, width = 14, height = 6, units = "in")
cat("Saved: ", figS3B_out, "\n", sep = "")

if (exists("log_output", mode = "function")) {
  log_output(
    file_path = figS3B_out,
    file_type = "figure",
    figure_id = "Fig S3B",
    description = "Batch correction diagnostics (PCA before + UMAP after) by batch1234",
    script_name = "17_supp_figS3_scRNA_qc.R",
    resolution = "vector pdf"
  )
}

# ----------------------------
# Optional: library_id对照图（仅用于自查，不用于投稿）
# ----------------------------
cat("\nGenerating library_id对照图（仅用于自查）...\n")
lib_umap_out <- file.path(LOG_DIR, "17_tmp_umap_by_library_id.pdf")
p_lib_umap <- DimPlot(
  obj_50k,
  reduction = "umap",
  group.by = "library_id",
  raster = TRUE
) +
  theme_bw() +
  ggtitle("UMAP by library_id (for internal review only)") +
  theme(legend.position = "bottom")
ggsave(lib_umap_out, p_lib_umap, width = 10, height = 8, units = "in")
cat("Saved: ", lib_umap_out, "\n", sep = "")

# ----------------------------
# Session info
# ----------------------------
capture.output(sessionInfo(), file = SESSIONINFO_FILE)
cat("\nSessionInfo saved: ", SESSIONINFO_FILE, "\n", sep = "")
cat("\n=== Done ===\n")

# Close logging
sink()