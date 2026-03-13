# =========================================================================
# Script: 09_fig3_scRNA.R
# Purpose: Render Fig3A–E from reviewer DietSeurat object (GSE274583 CD4+ T cells).
# Author: Zuo
# Date: 2026-03-04
#
# Inputs:
#   - data/processed_full/single_cell/gse274583/gse274583_cd4_integrated_reviewer_diet.rds
#   - R/00_theme_setup.R  # Theme & color palette (single source for figure styling)
#
# Outputs (to output/figures_main):
#   - Fig3A_umap_clusters.pdf
#   - Fig3B_umap_sig11_score.pdf
#   - Fig3C_umap_th2_score.pdf
#   - Fig3D_vln_scores_by_cluster.pdf
#   - Fig3E_dotplot_locked_genes.pdf
#   - output/logs/09_fig3_scRNA_sessionInfo.txt
# =========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ----------------------------
# Project root + path setup
# ----------------------------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).",
       call. = FALSE)
}

DATA_PROCESSED_FULL_DIR <- file.path(base_dir, "data", "processed_full")
FIG_DIR                <- file.path(base_dir, "output", "figures_main")
LOG_DIR                <- file.path(base_dir, "output", "logs")
THEME_FILE             <- file.path(base_dir, "R", "00_theme_setup.R")

SC_REVIEWER_RDS <- file.path(
  DATA_PROCESSED_FULL_DIR,
  "single_cell", "gse274583",
  "gse274583_cd4_integrated_reviewer_diet.rds"
)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

SESSIONINFO_FILE <- file.path(LOG_DIR, "09_fig3_scRNA_sessionInfo.txt")

if (!file.exists(SC_REVIEWER_RDS)) stop(paste0("Missing reviewer RDS: ", SC_REVIEWER_RDS), call. = FALSE)
if (!file.exists(THEME_FILE)) stop(paste0("Missing theme file: ", THEME_FILE), call. = FALSE)

# ----------------------------
# Load theme (must define theme + palettes used in successful build)
# ----------------------------
source(THEME_FILE)
# Set plot theme and color palette
set_plot_theme()

# ----------------------------
# Read reviewer diet object
# ----------------------------
obj <- readRDS(SC_REVIEWER_RDS)

# Required fields sanity checks
req_meta <- c("library_id", "batch1234", "percent.mt", "seurat_clusters", "Sig11_Score1", "Th2_Score1")
missing_meta <- setdiff(req_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required meta columns in reviewer RDS: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

# Validate batch1234 values
if (any(is.na(obj$batch1234))) {
  stop("batch1234 contains NA values. Please ensure 04c writes batch1234 to the reviewer RDS.", call. = FALSE)
}

valid_batches <- paste0("Batch", 1:4)
invalid_batches <- setdiff(unique(obj$batch1234), valid_batches)
if (length(invalid_batches) > 0) {
  stop("Invalid batch1234 values found: ", paste(invalid_batches, collapse = ", "), ". Valid values are: ", paste(valid_batches, collapse = ", "), call. = FALSE)
}

# Create audit log
AUDIT_FILE <- file.path(LOG_DIR, "09_fig3_scRNA_meta_audit.tsv")
audit_data <- list(
  "Full object batch distribution" = table(obj$batch1234),
  "Full object batch vs library_id" = table(obj$batch1234, obj$library_id),
  "Full object batch vs clusters" = table(obj$batch1234, obj$seurat_clusters)
)

# Write audit log
sink(AUDIT_FILE)
cat("=== 09_fig3_scRNA Meta Audit ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

for (name in names(audit_data)) {
  cat(name, "\n")
  print(audit_data[[name]])
  cat("\n")
}

cat("=== End Audit ===\n")
sink()

if (!("umap" %in% names(obj@reductions))) {
  stop("UMAP reduction not found in reviewer RDS. Ensure Script 04c saved umap in DietSeurat.", call. = FALSE)
}

locked_genes <- c("CLCA1","SERPINB2","CPA3","CCR3","HDC","MUC5AC","IL4","MUC5B","IL13","CCL24","POSTN")
# For CD4+ T cell scRNA, some epithelial mucins may be absent; keep list but handle missing gracefully.
locked_genes_present <- intersect(locked_genes, rownames(obj))

# ----------------------------
# Sampling for rendering only (50,000) with stratified sampling by batch1234
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

cells_50k <- sample_cells(obj, max_cells = 50000, strata = "batch1234")
obj_50k <- subset(obj, cells = cells_50k)

# Add batch1234 audit for sampled object
sink(AUDIT_FILE, append = TRUE)
cat("\n=== Sampled (50k) Object Audit ===\n")
cat("Batch distribution in sampled object:\n")
print(table(obj_50k$batch1234))
cat("Batch proportions in sampled object:\n")
print(prop.table(table(obj_50k$batch1234)))
cat("\n=== End Sampled Audit ===\n")
sink()

# ----------------------------
# Figure output helpers
# ----------------------------
save_pdf <- function(p, filename, width = 7, height = 6) {
  ggsave(filename = file.path(FIG_DIR, filename), plot = p, width = width, height = height, units = "in")
}

# ----------------------------
# Fig 3A: UMAP clusters
# ----------------------------
# Set seed for consistent label placement
set.seed(2026)

# 1) 画底图（不带任何 label）
p3a <- DimPlot(
  obj_50k,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.2,  # 增大点的大小，使细胞图更清晰
  raster = FALSE,  # 关闭光栅化，提高清晰度
  cols = trae_colors_celltype
) +
  theme_allergy() +
  NoLegend() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  coord_cartesian(clip = "off") +  # Prevent labels from being clipped
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  # Add margin to accommodate labels

# 2) 叠加 cluster 编号（黑字，不要框）
# 从 UMAP embedding + meta 里算每个 cluster 的中心点（使用 median，更稳）
emb <- as.data.frame(Embeddings(obj_50k, reduction = "umap"))
# 确保列名正确，使用实际的列名
colnames(emb) <- c("UMAP_1", "UMAP_2")  # 强制设置列名
emb$seurat_clusters <- obj_50k$seurat_clusters

centers <- emb %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# 使用 ggrepel::geom_text_repel() 添加纯文字标签
p3a <- p3a +
  ggrepel::geom_text_repel(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters),
    color = "black",
    size = 4,
    seed = 2026,  # 固定种子，确保标签位置一致
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2
  )

# 保存PDF和PNG
save_pdf(p3a, "Fig3A_umap_clusters.pdf", width = 7.5, height = 6.2)
# 保存PNG
png(file.path(FIG_DIR, "Fig3A_umap_clusters.png"), width = 7.5, height = 6.2, units = "in", res = 300)
print(p3a)
dev.off()

# ----------------------------
# Calculate score quantiles for consistent color scales
# ----------------------------
# 使用2%-98%的分位数，增强热点可见性
sig11_quantiles <- quantile(obj$Sig11_Score1, probs = c(0.02, 0.98), na.rm = TRUE)
th2_quantiles <- quantile(obj$Th2_Score1, probs = c(0.02, 0.98), na.rm = TRUE)

# ----------------------------
# 通用函数：绘制UMAP灰底+top fraction热点
# ----------------------------
plot_umap_top_frac_rank <- function(seu, score_col, 
                                     top_frac = 0.05, 
                                     base_color = "grey85", 
                                     top_color  = "#1f4e99", 
                                     pt_base = 0.20, 
                                     pt_top  = 0.40, 
                                     reduction = "umap") {

  emb <- as.data.frame(Seurat::Embeddings(seu, reduction = reduction))
  stopifnot(ncol(emb) >= 2)
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")

  score <- seu@meta.data[[score_col]]
  if (is.null(score)) stop(paste0("Missing meta column: ", score_col))

  df <- cbind(emb[, c("UMAP_1", "UMAP_2")], score = score)

  # strict top fraction by rank (robust to ties)
  n <- nrow(df)
  cutoff_rank <- floor((1 - top_frac) * n)
  df$top <- rank(df$score, ties.method = "first") > cutoff_rank

  ggplot(df, aes(UMAP_1, UMAP_2)) + 
    geom_point(data = df[df$top == FALSE, ], color = base_color, size = pt_base) + 
    geom_point(data = df[df$top == TRUE,  ], color = top_color,  size = pt_top) + 
    theme_allergy() + 
    labs(x = "UMAP_1", y = "UMAP_2")
}

# ----------------------------
# Fig 3B: UMAP Sig11_Score1 (base gray + strict top 5% by rank)
# ----------------------------
p3b <- plot_umap_top_frac_rank(obj_50k, score_col = "Sig11_Score1", top_frac = 0.05)

save_pdf(p3b, "Fig3B_umap_sig11_score.pdf", width = 7.5, height = 6.2)
png(file.path(FIG_DIR, "Fig3B_umap_sig11_score.png"),
    width = 7.5, height = 6.2, units = "in", res = 300)
print(p3b)
dev.off()

# ----------------------------
# Fig 3C: UMAP Th2_Score1 (base gray + strict top 5% by rank)
# ----------------------------
p3c <- plot_umap_top_frac_rank(obj_50k, score_col = "Th2_Score1", top_frac = 0.05)

save_pdf(p3c, "Fig3C_umap_th2_score.pdf", width = 7.5, height = 6.2)
png(file.path(FIG_DIR, "Fig3C_umap_th2_score.png"),
    width = 7.5, height = 6.2, units = "in", res = 300)
print(p3c)
dev.off()

# 快速自检：确认top5真的是5%
message("Fig3B top fraction (Sig11): ", mean(rank(obj_50k$Sig11_Score1, ties.method="first") > floor(0.95*ncol(obj_50k))))
message("Fig3C top fraction (Th2): ",   mean(rank(obj_50k$Th2_Score1,   ties.method="first") > floor(0.95*ncol(obj_50k))))

# ----------------------------
# Fig 3D: Violin (scores by cluster) - only Sig11_Score1
# ----------------------------
# Keep panel style simple and stable but more violin-like.
# Calculate 2nd and 98th percentiles to avoid extreme values compressing the main distribution
sig11_02nd <- quantile(obj_50k$Sig11_Score1, probs = 0.02, na.rm = TRUE)
sig11_98th <- quantile(obj_50k$Sig11_Score1, probs = 0.98, na.rm = TRUE)

# 创建violin plot
p3d <- VlnPlot(
  obj_50k,
  features = "Sig11_Score1",
  group.by = "seurat_clusters",
  pt.size = 0,  # 不画散点
  adjust = 2,  # 增加平滑度
  cols = trae_colors_celltype  # 使用与Fig3A相同的配色
) +
  theme_allergy() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +  # 标签正向显示，字体大小调整为8
  coord_cartesian(ylim = c(sig11_02nd, 0.01)) +  # 使用coord_cartesian截断显示，上限设置为0.01
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +  # 叠加薄boxplot，不显示outlier
  NoLegend()  # 移除图例

# 保存PDF和PNG
save_pdf(p3d, "Fig3D_vln_scores_by_cluster.pdf", width = 8.5, height = 5.0)
# 保存PNG
png(file.path(FIG_DIR, "Fig3D_vln_scores_by_cluster.png"), width = 8.5, height = 5.0, units = "in", res = 300)
print(p3d)
dev.off()

# ----------------------------
# Fig 3E: DotPlot for locked genes (use FULL object, not sampled)
# ----------------------------
message(paste("Locked genes present:", length(locked_genes_present)))
message(paste("Locked genes:", paste(locked_genes_present, collapse=", ")))
message(paste("Number of clusters:", length(unique(obj$seurat_clusters))))

p3e_genes_candidate <- locked_genes_present
if (length(p3e_genes_candidate) == 0) {
  message("Fig3E: no locked genes present in object. Writing placeholder.")
  p3e <- ggplot() + theme_void() + ggtitle("DotPlot skipped: no locked genes present")
  save_pdf(p3e, "Fig3E_dotplot_locked_genes.pdf", width = 10.5, height = 6)
} else {
  # 1) 先生成 DotPlot 的数据（让 Seurat 自己处理 layer/assay）
  dp <- tryCatch({
    Seurat::DotPlot(obj, features = p3e_genes_candidate, group.by = "seurat_clusters")$data
  }, error = function(e) {
    message("Fig3E: DotPlot data generation failed: ", e$message)
    NULL
  })
  
  if (is.null(dp) || nrow(dp) == 0) {
    message("Fig3E: DotPlot$data is empty. Writing placeholder.")
    p3e <- ggplot() + theme_void() + ggtitle("DotPlot skipped: empty DotPlot data")
    save_pdf(p3e, "Fig3E_dotplot_locked_genes.pdf", width = 10.5, height = 6)
  } else {
    # 2) 从 DotPlot$data 中筛基因：至少一个 cluster pct.exp > 0
    genes_to_plot <- dp %>%
      dplyr::group_by(features.plot) %>%
      dplyr::summarise(mx_pct = max(pct.exp, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(mx_pct > 0) %>%
      dplyr::pull(features.plot)
    
    message("Fig3E genes_to_plot: ", paste(genes_to_plot, collapse = ", "))
    
    if (length(genes_to_plot) == 0) {
      message("Fig3E: all candidate genes have 0 expression across clusters. Writing placeholder.")
      p3e <- ggplot() + theme_void() + ggtitle("DotPlot skipped: no expressed genes detected")
      save_pdf(p3e, "Fig3E_dotplot_locked_genes.pdf", width = 10.5, height = 6)
    } else {
      # 防御性检查：防止seurat_clusters变为NA或只有一个cluster
      stopifnot(!all(is.na(obj$seurat_clusters)))
      stopifnot(length(unique(obj$seurat_clusters[!is.na(obj$seurat_clusters)])) > 1)
      
      # 固定基因顺序：免疫相关在左，上皮/黏蛋白在右
      immune_genes_ordered <- c("IL4", "IL13", "CCR3", "HDC", "CCL24", "SERPINB2", "CPA3")
      epithelial_genes_ordered <- c("CLCA1", "MUC5B")
      
      # 按固定顺序排序基因
      genes_sorted <- c(
        intersect(immune_genes_ordered, genes_to_plot),
        intersect(epithelial_genes_ordered, genes_to_plot)
      )
      
      # 显式固定assay/layer（Seurat v5）
      DefaultAssay(obj) <- "RNA"  # 假设使用RNA assay
      
      # 不尝试排序cluster，直接创建DotPlot
      message("Creating DotPlot without cluster sorting")
      
      # 创建DotPlot
      p3e <- Seurat::DotPlot(
        obj, 
        features = genes_sorted, 
        group.by = "seurat_clusters",
        cols = c("lightgrey", "blue"),  # 统一颜色方案
        dot.scale = 8
      ) +
        RotatedAxis() +
        theme_allergy() +
        ggtitle("Locked genes expression across clusters") +
        labs(caption = "Epithelial mucin genes are expected to be low/absent in CD4+ T cells") +
        theme(
          plot.caption = element_text(size = 6, hjust = 0.5, face = "italic")  # 字号更小
        )
      
      # 修改size legend刻度为0 / 5 / 10 / 20
      p3e <- p3e +
        guides(size = guide_legend(breaks = c(0, 5, 10, 20), title = "Percent expressed"))
      
      # 保存PDF和PNG
      save_pdf(p3e, "Fig3E_dotplot_locked_genes.pdf", width = 10.5, height = 6)
      # 保存PNG
      png(file.path(FIG_DIR, "Fig3E_dotplot_locked_genes.png"), width = 10.5, height = 6, units = "in", res = 300)
      print(p3e)
      dev.off()
    }
  }
}

message("DotPlot creation completed")

# ----------------------------
# Session info
# ----------------------------
capture.output(sessionInfo(), file = SESSIONINFO_FILE)

message("Done. Figures saved to: ", FIG_DIR)
message("SessionInfo saved to: ", SESSIONINFO_FILE)