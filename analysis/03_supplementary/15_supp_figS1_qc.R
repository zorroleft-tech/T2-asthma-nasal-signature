# =========================================================================
# Script: 15_supp_figS1_qc.R
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
  library(limma)
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
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
DIET_DIR <- file.path(DATA_DIR, "processed_diet")
OUTPUT_DIR <- file.path(base_dir, "output")
SUPP_DIR <- file.path(OUTPUT_DIR, "supplement")

# Create directory if not exist
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

# Log base directory
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Supp figures directory: ", SUPP_DIR, "\n\n"))

# INPUT PATHS (using existing files)
expr_file <- file.path(PROCESSED_DIR, "expr_gse152004&GPL11154__FULL.rds")
diet_file <- file.path(DIET_DIR, "GSE152004_diet_expr.rds")
meta_file <- file.path(DATA_DIR, "derived", "GSE152004_T2_labels.csv")

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

# -----------------------------# 1) Theme + color system# -----------------------------
# Try to use your existing theme setup if present. Otherwise fall back to a clean minimal theme.
theme_setup_file <- file.path(base_dir, "R/00_theme_setup.R")

if (file.exists(theme_setup_file)) {
  source(theme_setup_file)
} else {
  message("Theme setup not found at: ", theme_setup_file)
  message("Falling back to minimal theme_trae_publication() and default colors.")
  theme_trae_publication <- function() theme_minimal(base_size = 24) +
    theme(
      axis.title = element_text(size = 28, face = "bold"),
      axis.text = element_text(size = 22, face = "bold"),
      plot.title = element_text(size = 32, face = "bold")
    )
  trae_colors_main <- c("#D64550", "#2E86AB", "#333333")
}

# Utility: save both png and pdf using ggsave
save_png_pdf <- function(p, file_base, width_in = 12, height_in = 10, dpi = 300) {
  ggsave(paste0(file_base, ".pdf"), p, width = width_in, height = height_in, units = "in")
  ggsave(paste0(file_base, ".png"), p, width = width_in, height = height_in, units = "in", dpi = dpi)
}

# -----------------------------# 2) Read data for both plots# -----------------------------if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file)
if (!file.exists(diet_file)) stop("Diet expression file not found: ", diet_file)
if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)

# Read expression matrix and metadata
expr_full <- readRDS(expr_file)
expr_diet <- readRDS(diet_file)

# Check that both matrices have the same number of columns
stopifnot(ncol(expr_full) == ncol(expr_diet))

# Replace FULL column names with DIET column names (HRxxxx format)
colnames(expr_full) <- colnames(expr_diet)

# Add safety check to prevent sample order mismatch
common_genes <- intersect(rownames(expr_diet), rownames(expr_full))
common_genes <- common_genes[!is.na(common_genes)]

# Use first 10 common genes as sentinel
sentinel <- head(common_genes, 10)

cors <- sapply(sentinel, function(g) {
  suppressWarnings(cor(as.numeric(expr_full[g, ]), as.numeric(expr_diet[g, ]), use = "pairwise.complete.obs"))
})

cat("Sentinel gene correlations (FULL vs DIET):\n")
print(round(cors, 3))

if (median(cors, na.rm = TRUE) < 0.95) {
  stop("Column order mismatch suspected: FULL columns may not align with DIET columns. Do not proceed to DEG.")
}

# Use the FULL matrix for analysis
expr_mat <- expr_full

# Add log2 transformation to fix large logFC values
if(max(expr_mat, na.rm=TRUE) > 50) { 
  expr_mat <- log2(expr_mat + 1) 
}

meta_data <- read.csv(meta_file, header = TRUE)

# Print sample size for figure caption
cat("Post-QC discovery sample size: n = ", ncol(expr_mat), "\n\n", sep = "")

# Handle sample ID column compatibility
sample_id_col <- NULL
possible_sample_cols <- c("Sample_ID", "sample_id", "GSM", "sample")
for (col in possible_sample_cols) {
  if (col %in% colnames(meta_data)) {
    sample_id_col <- col
    break
  }
}
if (is.null(sample_id_col)) {
  stop("No sample ID column found in metadata. Expected one of: ", paste(possible_sample_cols, collapse = ", "))
}

# Handle T2 status column compatibility
t2_col <- NULL
possible_t2_cols <- c("T2_Status", "t2_status", "discovery_label", "T2_label", "t2_label")
for (col in possible_t2_cols) {
  if (col %in% colnames(meta_data)) {
    t2_col <- col
    break
  }
}
if (is.null(t2_col)) {
  stop("No T2 status column found in metadata. Expected one of: ", paste(possible_t2_cols, collapse = ", "))
}

# Align meta to expr order
expr_samples <- colnames(expr_mat)
meta_data <- meta_data[match(expr_samples, meta_data[[sample_id_col]]), ]

if (any(is.na(meta_data[[sample_id_col]]))) {
  stop("Some samples in expression matrix were not found in metadata ", sample_id_col, ".")
}
if (any(expr_samples != meta_data[[sample_id_col]])) {
  stop("Sample ID mismatch between expression matrix and metadata (after matching).")
}

# Normalize label levels to T2_high / T2_low (valid R names)
meta_data$T2_Status <- as.character(meta_data[[t2_col]])
meta_data$T2_Status[grepl("high", meta_data$T2_Status, ignore.case = TRUE)] <- "T2_high"
meta_data$T2_Status[grepl("low",  meta_data$T2_Status, ignore.case = TRUE)] <- "T2_low"
meta_data$T2_Status <- factor(meta_data$T2_Status, levels = c("T2_low", "T2_high"))

# -----------------------------# 3) Figure S1A: Volcano plot (calculate DEG using limma)# -----------------------------
# Prepare design matrix for limma
design <- model.matrix(~ 0 + meta_data$T2_Status)
colnames(design) <- levels(meta_data$T2_Status)

# Fit linear model using limma
fit <- lmFit(expr_mat, design)

# Define contrast (T2_high vs T2_low)
contrast <- makeContrasts("T2_high - T2_low", levels = design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Extract results
deg_results <- topTable(fit, adjust = "BH", sort.by = "P", n = Inf)

# Convert to tibble and rename columns
deg_tbl <- as_tibble(deg_results, rownames = "symbol") %>%
  rename(logFC = logFC, adj.P.Val = adj.P.Val)

# Add columns for plotting
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
    alpha = 0.60, size = 1.2, stroke = 0
  ) +
  scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2C7FB8", "NS" = "grey75")) +
  geom_point(
    data = core_tbl,
    aes(x = logFC, y = neglog10_fdr),
    fill = "#E63946", color = "black", shape = 21, size = 3.5, stroke = 0.8
  ) +
  geom_label_repel(
    data = core_tbl,
    aes(x = logFC, y = neglog10_fdr, label = symbol),
    size = 8, fontface = "bold",
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
  theme(
    legend.position = "none",
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26, face = "bold")
  )

save_png_pdf(p_volcano, out_volcano_base, width_in = 12, height_in = 10, dpi = 300)

cat("✅ Fig S1A saved:\n")
cat("   ", paste0(out_volcano_base, ".pdf"), "\n")
cat("   ", paste0(out_volcano_base, ".png"), "\n\n")

# -----------------------------# 4) Figure S1B: PCA plot# -----------------------------

# Filter out constant or zero-variance genes for PCA
variances <- apply(expr_mat, 1, var)
non_constant_genes <- which(variances > 0)
expr_mat_pca <- expr_mat[non_constant_genes, ]

# PCA (scale+center as in your original) - transpose matrix to analyze samples
pca_res <- prcomp(t(expr_mat_pca), scale. = TRUE, center = TRUE)

pc_df <- tibble(
  sample_id = expr_samples,
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
    legend.text = element_text(size = 24, face = "bold"),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "black", linewidth = 0.5),
    legend.margin = margin(t = 3, r = 6, b = 3, l = 6),
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26, face = "bold")
  )

save_png_pdf(p_pca, out_pca_base, width_in = 12, height_in = 10, dpi = 300)

cat("✅ Fig S1B saved:\n")
cat("   ", paste0(out_pca_base, ".pdf"), "\n")
cat("   ", paste0(out_pca_base, ".png"), "\n\n")

# -----------------------------
# 4) Quick audit outputs (for caption alignment)
# -----------------------------
cat("=== Figure S1 QC summary (for caption alignment) ===\n")
cat("Post-QC discovery sample size used for PCA: n = ", ncol(expr_mat), "\n", sep = "")
cat("Volcano highlighted genes (locked 11): ", paste(locked_11, collapse = ", "), "\n", sep = "")
cat("PC1 variance explained: ", round(pc1_var, 1), "%\n", sep = "")
cat("PC2 variance explained: ", round(pc2_var, 1), "%\n", sep = "")
cat("All done.\n")
