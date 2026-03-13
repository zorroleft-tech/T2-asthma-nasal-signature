# =========================================================================
# Script: 15_supp_figS1_qc_fixed.R
# Purpose: Generate Figure S1 for supplementary materials
#   - Figure S1A: Volcano plot (Discovery cohort DE results; highlight locked 11 genes)
#   - Figure S1B: PCA plot (Discovery cohort post-QC expression; T2-high vs T2-low + 95% ellipses)
# Author: Zuo
# Date: 2026-03-01
#
# Inputs:
#   - 02_data_processed/phase2_outputs/GSE152004_full_transcriptome_deg.csv  # DE results
#   - 02_data_processed/phase2_outputs/gse152004_qc_passed.rds  # Post-QC expression matrix
#   - 02_data_processed/phase3_outputs/GSE152004_discovery_labels_final.csv  # T2 status metadata
#
# Outputs:
#   - outputs/figures/supp/FigS1A_GSE152004_volcano.pdf  # Volcano plot
#   - outputs/figures/supp/FigS1A_GSE152004_volcano.png  # Volcano plot
#   - outputs/figures/supp/FigS1B_GSE152004_PCA.pdf  # PCA plot
#   - outputs/figures/supp/FigS1B_GSE152004_PCA.png  # PCA plot
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
DATA_DIR <- file.path(base_dir, "data")
OUTPUT_DIR <- file.path(base_dir, "output")
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures_main")
SUPP_DIR <- file.path(FIGURES_DIR, "supp")

# Create directory if not exist
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

# Log base directory
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Supp figures directory: ", SUPP_DIR, "\n\n"))

# INPUT PATHS (adjust if your actual files differ)
deg_file <- file.path(DATA_DIR, "derived/limma_full_stats.csv")
expr_file <- file.path(DATA_DIR, "processed/expr_gse152004&GPL11154__FULL.rds")
meta_file <- file.path(DATA_DIR, "derived/GSE152004_T2_labels.csv")

# OUTPUT PATHS
out_volcano_base <- file.path(SUPP_DIR, "FigS1A_GSE152004_volcano")
out_pca_base     <- file.path(SUPP_DIR, "FigS1B_GSE152004_PCA")

# Locked 11 genes (do not change)
locked_11 <- c(
  "CLCA1","SERPINB2","CPA3","CCR3","HDC",
  "MUC5AC","IL4","MUC5B","IL13","CCL24","POSTN"
)

# Significance thresholds (match your original plot logic)
fdr_thr   <- 0.05
logfc_thr <- 1

# -----------------------------
# 1) Theme + color system
# -----------------------------
# Try to use your existing theme setup if present. Otherwise fall back to a clean minimal theme.
theme_setup_file <- file.path(base_dir, "R/00_theme_setup.R")

if (file.exists(theme_setup_file)) {
  source(theme_setup_file)
} else {
  message("Theme setup not found at: ", theme_setup_file)
  message("Falling back to minimal theme_trae_publication() and default colors.")
  theme_trae_publication <- function() theme_minimal(base_size = 14)
  trae_colors_main <- c("#D64550", "#2E86AB", "#333333")
}

# Utility: save both png and pdf using ggsave
save_png_pdf <- function(p, file_base, width_in = 12, height_in = 10, dpi = 300) {
  ggsave(paste0(file_base, ".pdf"), p, width = width_in, height = height_in, units = "in")
  ggsave(paste0(file_base, ".png"), p, width = width_in, height = height_in, units = "in", dpi = dpi)
}

# -----------------------------
# 2) Figure S1A: Volcano plot
# -----------------------------
if (!file.exists(deg_file)) stop("DEG file not found: ", deg_file)
deg_tbl <- readr::read_csv(deg_file, show_col_types = FALSE)

# Standardize column names (supports a few common variants)
# Required columns after mapping: symbol, logFC, adj.P.Val
if ("FDR" %in% colnames(deg_tbl)) deg_tbl <- dplyr::rename(deg_tbl, adj.P.Val = FDR)
if ("padj" %in% colnames(deg_tbl)) deg_tbl <- dplyr::rename(deg_tbl, adj.P.Val = padj)
if ("log2FoldChange" %in% colnames(deg_tbl)) deg_tbl <- dplyr::rename(deg_tbl, logFC = log2FoldChange)
if ("Symbol" %in% colnames(deg_tbl)) deg_tbl <- dplyr::rename(deg_tbl, symbol = Symbol)
if ("gene_symbol" %in% colnames(deg_tbl)) deg_tbl <- dplyr::rename(deg_tbl, symbol = gene_symbol)
if ("Gene" %in% colnames(deg_tbl)) deg_tbl <- dplyr::rename(deg_tbl, symbol = Gene)

required_cols <- c("symbol", "logFC", "adj.P.Val")
missing_cols <- setdiff(required_cols, colnames(deg_tbl))
if (length(missing_cols) > 0) {
  stop("DEG table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

deg_tbl <- deg_tbl %>%
  mutate(
    sig_status = case_when(
      adj.P.Val < fdr_thr & logFC >  logfc_thr ~ "Up",
      adj.P.Val < fdr_thr & logFC < -logfc_thr ~ "Down",
      TRUE ~ "NS"
    ),
    is_locked11 = ifelse(symbol %in% locked_11, "Yes", "No"),
    neglog10_fdr = -log10(adj.P.Val)
  )

bg_tbl   <- deg_tbl %>% filter(is_locked11 == "No")
core_tbl <- deg_tbl %>% filter(is_locked11 == "Yes")

# Volcano: keep original style
p_volcano <- ggplot() +
  geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(
    data = bg_tbl,
    aes(x = logFC, y = neglog10_fdr, color = sig_status),
    alpha = 0.30, size = 1.2, stroke = 0
  ) +
  scale_color_manual(values = c("Up" = "#F3D1D9", "Down" = "#D9D2F0", "NS" = "grey85")) +
  geom_point(
    data = core_tbl,
    aes(x = logFC, y = neglog10_fdr),
    fill = trae_colors_main[1], color = "black", shape = 21, size = 3.5, stroke = 0.8
  ) +
  geom_label_repel(
    data = core_tbl,
    aes(x = logFC, y = neglog10_fdr, label = symbol),
    size = 4, fontface = "bold",
    box.padding = 0.8, point.padding = 0.5,
    segment.color = "black", segment.size = 0.5,
    fill = scales::alpha("white", 0.8), color = "black",
    max.overlaps = Inf
  ) +
  labs(
    x = bquote(~log[2]~"(Fold Change)"),
    y = bquote(~-log[10]~"(Adjusted P-value)")
  ) +
  theme_trae_publication() +
  theme(legend.position = "none")

save_png_pdf(p_volcano, out_volcano_base, width_in = 12, height_in = 10, dpi = 300)

cat("✅ Fig S1A saved:\n")
cat("   ", paste0(out_volcano_base, ".pdf"), "\n")
cat("   ", paste0(out_volcano_base, ".png"), "\n\n")

# -----------------------------
# 3) Figure S1B: PCA plot
# -----------------------------
if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file)
if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)

expr_mat <- readRDS(expr_file)
meta_data <- readr::read_csv(meta_file, show_col_types = FALSE)

# Expect expr_mat: samples x genes with rownames = sample IDs
if (is.null(rownames(expr_mat))) stop("Expression matrix has no rownames (sample IDs).")

# Expect meta_data: has Sample_ID and T2_Status
if (! ("Sample_ID" %in% colnames(meta_data))) stop("Metadata missing 'Sample_ID' column: ", meta_file)
if (! ("T2_Status" %in% colnames(meta_data))) stop("Metadata missing 'T2_Status' column: ", meta_file)

# Align meta to expr order
meta_data <- meta_data[match(rownames(expr_mat), meta_data$Sample_ID), ]

if (any(is.na(meta_data$Sample_ID))) {
  stop("Some samples in expression matrix were not found in metadata Sample_ID.")
}
if (any(rownames(expr_mat) != meta_data$Sample_ID)) {
  stop("Sample ID mismatch between expression matrix and metadata (after matching).")
}

# Normalize label levels to T2-high / T2-low
meta_data$T2_Status <- as.character(meta_data$T2_Status)
meta_data$T2_Status[grepl("high", meta_data$T2_Status, ignore.case = TRUE)] <- "T2-high"
meta_data$T2_Status[grepl("low",  meta_data$T2_Status, ignore.case = TRUE)] <- "T2-low"
meta_data$T2_Status <- factor(meta_data$T2_Status, levels = c("T2-low", "T2-high"))

# PCA (scale+center as in your original)
pca_res <- prcomp(expr_mat, scale. = TRUE, center = TRUE)

pc_df <- tibble(
  sample_id = rownames(expr_mat),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  T2_Status = meta_data$T2_Status
)

var_explained <- summary(pca_res)$importance[2, ] * 100
pc1_var <- var_explained[1]
pc2_var <- var_explained[2]

p_pca <- ggplot(pc_df, aes(x = PC1, y = PC2, color = T2_Status, fill = T2_Status)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, linetype = "dashed", linewidth = 0.8) +
  geom_point(shape = 21, size = 1.8, color = "black", stroke = 0.15, alpha = 0.7) +
  scale_fill_manual(values = trae_colors_main[c(1, 2)]) +
  scale_color_manual(values = trae_colors_main[c(1, 2)]) +
  labs(
    x = sprintf("Principal Component 1 (%.1f%% Variance)", pc1_var),
    y = sprintf("Principal Component 2 (%.1f%% Variance)", pc2_var)
  ) +
  theme_trae_publication() +
  theme(
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "black", linewidth = 0.5),
    legend.margin = margin(t = 3, r = 6, b = 3, l = 6)
  )

save_png_pdf(p_pca, out_pca_base, width_in = 12, height_in = 10, dpi = 300)

cat("✅ Fig S1B saved:\n")
cat("   ", paste0(out_pca_base, ".pdf"), "\n")
cat("   ", paste0(out_pca_base, ".png"), "\n\n")

# -----------------------------
# 4) Quick audit outputs (for caption alignment)
# -----------------------------
cat("=== Figure S1 QC summary (for caption alignment) ===\n")
cat("Post-QC discovery sample size used for PCA: n = ", nrow(expr_mat), "\n", sep = "")
cat("Volcano highlighted genes (locked 11): ", paste(locked_11, collapse = ", "), "\n", sep = "")
cat("PC1 variance explained: ", round(pc1_var, 1), "%\n", sep = "")
cat("PC2 variance explained: ", round(pc2_var, 1), "%\n", sep = "")
cat("All done.\n")
