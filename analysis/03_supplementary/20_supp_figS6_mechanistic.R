# =========================================================================
# Script: 20_supp_figS6_mechanistic.R
# Purpose: Generate Supplementary Figure S6 with mechanistic figures
# Author: Zuo 
# Date: 2026-03-01
#
# Inputs:
#   - data/processed/expr_gse65204&GPL14550__FULL.rds  # GSE65204 expression data
#   - data/processed/phase3_outputs/GSE65204_pheno_raw.rds  # GSE65204 phenotype data
#   - R/00_theme_setup.R                    # Theme setup script
#   - R/03_index_logger.R                   # Logging script
#
# Outputs:
#   - output/supplement/FigS6_IgE_stratified_boxplots.pdf  # 5 representative genes IgE-stratified boxplots
#   - output/supplement/FigS6_IgE_stratified_boxplots.png  # (PNG version)
#   - output/logs/20_supp_figS6_mechanistic.log            # Log file
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# 1. Directory Setup
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "analysis"))) {
  stop("Please run this script from the project root directory.", call. = FALSE)
}

processed_dir <- file.path(base_dir, "data", "processed")
supplement_dir <- file.path(base_dir, "output", "supplement")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "20_supp_figS6_mechanistic.log")
sink(log_file, append = FALSE, split = TRUE)

cat("=======================================================\n")
cat("=== Generating Supplementary Figure S6 (Mechanistic) ===\n")
cat("=======================================================\n")

# Source theme setup and utility functions
source(file.path(base_dir, "R", "03_index_logger.R"))
theme_setup_file <- file.path(base_dir, "R", "00_theme_setup.R")
if (file.exists(theme_setup_file)) {
  source(theme_setup_file)
  if (exists("set_plot_theme")) set_plot_theme()
}

# Save function
save_png_pdf <- function(p, file_base, width_in = 12, height_in = 4.5, dpi = 300) {
  ggsave(paste0(file_base, ".pdf"), p, width = width_in, height = height_in, device = "pdf")
  ggsave(paste0(file_base, ".png"), p, width = width_in, height = height_in, dpi = dpi, bg = "white")
}

# 2. Load Data (GSE65204 clinical anchor)
cat("\n[1/3] Loading GSE65204 expression and phenotype data...\n")
# Use the actual file paths
phase3_dir <- file.path(processed_dir, "phase3_outputs")
expr_file <- file.path(processed_dir, "expr_gse65204&GPL14550__FULL.rds")
pheno_file <- file.path(phase3_dir, "GSE65204_pheno_raw.rds")

if (!file.exists(expr_file) || !file.exists(pheno_file)) {
  stop("Missing required GSE65204 processed data. Please check data/processed/ directory.")
}

expr_mat <- readRDS(expr_file)
pheno_data <- readRDS(pheno_file)

# 3. Extract IgE and define Tertiles
cat("\n[2/3] Extracting IgE values and calculating tertiles...\n")

# Sample ID column with explicit priority
if ("std_sample_id" %in% names(pheno_data)) {
  sample_id_col <- "std_sample_id"
} else if ("Sample_ID" %in% names(pheno_data)) {
  sample_id_col <- "Sample_ID"
} else if ("sample_id" %in% names(pheno_data)) {
  sample_id_col <- "sample_id"
} else if ("geo_accession" %in% names(pheno_data)) {
  sample_id_col <- "geo_accession"
} else if (any(grepl("^GSM", names(pheno_data)))) {
  sample_id_col <- grep("^GSM", names(pheno_data), value = TRUE)[1]
} else {
  stop("No valid sample ID column found. Please manually specify the sample ID column.")
}

# Extract lnige from raw_characteristics_all
if ("raw_characteristics_all" %in% names(pheno_data)) {
  ige_col <- "raw_characteristics_all"
  # Extract lnige values from the string
  pheno_data$lnige <- as.numeric(gsub(".*lnige: ([0-9.]+).*", "\\1", pheno_data[[ige_col]]))
} else {
  stop("No IgE data found in GSE65204 phenotype data.")
}

# Align samples
common_samples <- intersect(colnames(expr_mat), pheno_data[[sample_id_col]])
if(length(common_samples) < 10) stop("Not enough overlapping samples between expression and phenotype.")

# Log验收信息
cat(sprintf("  - Expression matrix dimensions (before filtering): %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  - IgE column: %s (extracted lnige)\n", ige_col))
cat(sprintf("  - Sample ID column: %s\n", sample_id_col))
cat(sprintf("  - Number of common samples: %d\n", length(common_samples)))

expr_mat <- expr_mat[, common_samples]
pheno_data <- pheno_data[match(common_samples, pheno_data[[sample_id_col]]), ]

# Process IgE
ige_num <- pheno_data$lnige
valid_idx <- !is.na(ige_num)
cat(sprintf("  - Number of samples with valid IgE values: %d\n", sum(valid_idx)))

expr_mat <- expr_mat[, valid_idx]
ige_num <- ige_num[valid_idx]

# Use the extracted lnige values directly
ln_ige <- ige_num

# Calculate IgE tertiles with fallback for non-unique breaks
tryCatch({
  ige_tertiles <- quantile(ln_ige, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  ige_group <- cut(ln_ige, breaks = ige_tertiles, labels = c("Low", "Medium", "High"), include.lowest = TRUE)
}, error = function(e) {
  cat("  - Warning: Using ntile() instead of cut() due to non-unique quantile breaks\n")
  ige_group <- factor(dplyr::ntile(ln_ige, 3), levels = c(1, 2, 3), labels = c("Low", "Medium", "High"))
})

cat("  - Samples grouped by IgE tertiles:\n")
print(table(ige_group))

# 4. Define genes and module mapping
cat("\n[3/3] Generating Boxplots for Representative Genes...\n")
# Select 5 highly representative genes from the locked 11-gene signature
plot_genes <- c("SERPINB2", "MUC5AC", "MUC5B", "CCL24", "CPA3")
available_genes <- intersect(plot_genes, rownames(expr_mat))

cat(sprintf("  - Requested genes: %s\n", paste(plot_genes, collapse = ", ")))
cat(sprintf("  - Available genes in expression matrix: %s\n", paste(available_genes, collapse = ", ")))

# Check if all requested genes are available
if (length(available_genes) != length(plot_genes)) {
  missing_genes <- setdiff(plot_genes, available_genes)
  stop(sprintf("Missing required genes: %s", paste(missing_genes, collapse = ", ")))
}

# Assign generic modules for aesthetic coloring
module_def <- data.frame(
  Gene = c("SERPINB2", "MUC5AC", "MUC5B", "CCL24", "CPA3"),
  Module = c("Epithelial/Mucus", "Epithelial/Mucus", "Epithelial/Mucus", "Immune/Chemotaxis", "Immune/Chemotaxis"),
  stringsAsFactors = FALSE
)

plot_df <- map_dfr(available_genes, function(g) {
  tibble(
    Gene = g,
    IgE_Group = ige_group,
    Expression = as.numeric(expr_mat[g, ])
  )
}) %>%
  left_join(module_def, by = "Gene") %>%
  mutate(
    Gene = factor(Gene, levels = available_genes),
    IgE_Group = factor(IgE_Group, levels = c("Low", "Medium", "High"))
  )

# Color palette
mod_colors <- c("Epithelial/Mucus" = "#E64B35", "Immune/Chemotaxis" = "#4DBBD5")

# 5. Build and Save Plot
p_s6 <- ggplot(plot_df, aes(x = IgE_Group, y = Expression, fill = Module)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.6, width = 0.6) +
  facet_wrap(~Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = mod_colors) +
  labs(
    x = "IgE Tertile",
    y = "Expression (Z-score)",
    fill = "Functional Module"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_text(face = "bold")
  )

out_base <- file.path(supplement_dir, "FigS6_IgE_stratified_boxplots")
save_png_pdf(p_s6, out_base)

cat(sprintf("\n✅ Figure S6 successfully saved to:\n  - %s.pdf\n", out_base))

# 6. Final logging
if(exists("log_output")) {
  log_output(
    file_path = paste0(out_base, ".pdf"),
    file_type = "figure",
    figure_id = "Fig S6",
    description = "5 representative genes IgE-stratified expression boxplots",
    script_name = "20_supp_figS6_mechanistic.R",
    resolution = "300 dpi"
  )
}

cat("\nNOTE: Fig S6 genes are restricted to locked 11; no substitution from legacy gene sets is allowed.\n")
cat("\n=== Supplementary Figure S6 generation completed ===\n")

# Memory cleanup
rm(list = ls(pattern = "^(plot|expr|pheno|module|p_s6)", all.names = TRUE))
gc()

sink()