#!/usr/bin/env Rscript

# =========================================================================
# Script: 06_table2_model_coef.R
# Purpose: Generate Table 2 (locked 11-gene weight table) for manuscript
# Author: Zuo
# Date: 2026-02-27
#
# Inputs:
#   - data/derived/locked_weights.csv
#
# Outputs:
#   - output/tables_main/Table2_model_coefficients.csv  # Model coefficients table
#   - output/tables_main/Table2_model_coefficients.xlsx  # Excel version (optional)
#   - output/logs/06_table2_model_coef.log  # Log file
# =========================================================================

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
required_dirs <- c("data", "analysis", "data_preparation")
required_paths <- file.path(base_dir, required_dirs)
dirs_exist <- all(dir.exists(required_paths))
if (!dirs_exist) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define paths using file.path()
locked_weights_file <- file.path(base_dir, "data", "derived", "locked_weights.csv")
tables_dir <- file.path(base_dir, "output", "tables_main")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "06_table2_model_coef.log")

# Only set up sink if not already set
if (sink.number(type="output") == 0) {
  sink(log_file, append = TRUE)
}
if (sink.number(type="message") == 0) {
  sink(log_file, append = TRUE, type = "message")
}

# Ensure sink is always closed even if script stops early
on.exit({
  if (sink.number(type="message") > 0) sink(type="message")
  if (sink.number(type="output") > 0) sink(type="output")
}, add = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Input file: ", locked_weights_file, "\n")
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(tables_dir, "Table2_model_coefficients.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Table2_model_coefficients.xlsx")))
cat(sprintf("  - %s\n", log_file))

# Load required packages
library(tidyverse)

# Try to load openxlsx for Excel output
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  cat("Warning: openxlsx package not installed, Excel output will be skipped\n")
  create_xlsx <- FALSE
} else {
  library(openxlsx)
  create_xlsx <- TRUE
}

# Read locked_weights.csv
cat("\n=== Reading locked_weights.csv ===\n")
if (!file.exists(locked_weights_file)) {
  stop(sprintf("Locked weights file not found: %s", locked_weights_file))
}

weights_df <- read_csv(locked_weights_file, show_col_types = FALSE, progress = FALSE)
cat(sprintf("Read %d rows and %d columns\n", nrow(weights_df), ncol(weights_df)))
cat("Columns found:", names(weights_df), "\n")

# Identify gene column
gene_col_candidates <- c("Gene", "gene", "Gene_Symbol", "symbol")
gene_col <- intersect(names(weights_df), gene_col_candidates)
if (length(gene_col) != 1) {
  stop(sprintf("Could not identify unique gene column. Found: %s", paste(gene_col, collapse = ", ")))
}

# Identify weight column
weight_col_candidates <- c("Weight", "weight", "coef", "Coefficient", "beta")
weight_col <- intersect(names(weights_df), weight_col_candidates)
if (length(weight_col) != 1) {
  stop(sprintf("Could not identify unique weight column. Found: %s", paste(weight_col, collapse = ", ")))
}

# Rename columns to standard names
weights_df <- weights_df %>%
  rename(Gene = all_of(gene_col),
         Weight = all_of(weight_col))

# Convert Weight to numeric
weights_df$Weight <- as.numeric(weights_df$Weight)
if (any(is.na(weights_df$Weight))) {
  stop("Weight column contains non-numeric values")
}

# Check for duplicate genes
if (n_distinct(weights_df$Gene) != nrow(weights_df)) {
  stop("Duplicate genes found in locked_weights.csv")
}

# Define expected 11 genes
expected_genes <- c("CLCA1", "SERPINB2", "CPA3", "CCR3", "HDC", "MUC5AC", "IL4", "MUC5B", "IL13", "CCL24", "POSTN")

# Check if gene set is complete
actual_genes <- sort(weights_df$Gene)
expected_genes_sorted <- sort(expected_genes)

if (!identical(actual_genes, expected_genes_sorted)) {
  missing_genes <- setdiff(expected_genes_sorted, actual_genes)
  extra_genes <- setdiff(actual_genes, expected_genes_sorted)
  cat("Actual genes found:", actual_genes, "\n")
  if (length(missing_genes) > 0) {
    cat("Missing genes:", missing_genes, "\n")
  }
  if (length(extra_genes) > 0) {
    cat("Extra genes:", extra_genes, "\n")
  }
  stop("Gene set does not match expected locked 11 genes")
}

cat("\n=== Gene set validation passed ===\n")
cat("Found all 11 locked genes:\n")
print(actual_genes)

# Create functional notes dictionary
functional_notes <- c(
  "SERPINB2" = "IL-13-responsive marker; T2-high airway inflammation context",
  "CLCA1" = "Epithelial IL-13-responsive marker; mucus-related program",
  "POSTN" = "Extracellular matrix remodeling marker; T2-high context",
  "MUC5AC" = "Airway goblet cell mucin; mucus hypersecretion context",
  "MUC5B" = "Airway mucin; mucus/remodeling context",
  "IL4" = "Type 2 cytokine; Th2 signaling axis",
  "IL13" = "Type 2 cytokine; mucus program driver",
  "CCL24" = "Eotaxin-2; eosinophil recruitment axis",
  "CCR3" = "Eosinophil chemokine receptor; trafficking axis",
  "CPA3" = "Mast cell protease marker",
  "HDC" = "Histamine synthesis; mast cell/basophil axis"
)

# Generate Table 2
# Note: Weight sign only indicates coefficient direction in multivariate model, not equivalent to differential expression up/downregulation
table2 <- weights_df %>%
  mutate(
    Abs_Weight = abs(Weight),
    Direction = case_when(
      Weight > 0 ~ "Up",
      Weight < 0 ~ "Down",
      TRUE ~ "Zero"
    ),
    Functional_note = functional_notes[Gene]
  ) %>%
  arrange(desc(Abs_Weight))

# Verify row count
if (nrow(table2) != 11) {
  stop(sprintf("Expected 11 rows in Table 2, but got %d", nrow(table2)))
}

# Save Table 2 as CSV
csv_output <- file.path(tables_dir, "Table2_model_coefficients.csv")
write_csv(table2, csv_output)
cat(sprintf("\nSaved Table 2 CSV: %s\n", basename(csv_output)))

# Save Table 2 as Excel if openxlsx is available
if (create_xlsx) {
  xlsx_output <- file.path(tables_dir, "Table2_model_coefficients.xlsx")
  write.xlsx(table2, xlsx_output, rowNames = FALSE)
  cat(sprintf("Saved Table 2 Excel: %s\n", basename(xlsx_output)))
}

# Print Table 2 summary
cat("\n=== Table 2 Summary ===\n")
cat(sprintf("Total rows: %d\n", nrow(table2)))
cat("Sorted by Abs_Weight descending\n")
cat("Note: Weight sign only indicates coefficient direction in multivariate model, not equivalent to differential expression up/downregulation\n")
print(table2[, c("Gene", "Weight", "Abs_Weight", "Direction")])

# Print session info for reproducibility
cat("\n=== Session Info ===\n")
sessionInfo()

# Check output file sizes
cat("\n=== Output File Sizes ===\n")
if (file.exists(csv_output)) {
  csv_size <- file.size(csv_output)
  cat(sprintf("Table2_model_coefficients.csv: %.2f KB\n", csv_size / 1024))
}
if (create_xlsx && file.exists(xlsx_output)) {
  xlsx_size <- file.size(xlsx_output)
  cat(sprintf("Table2_model_coefficients.xlsx: %.2f KB\n", xlsx_size / 1024))
}

# Memory cleanup
rm(list = ls(pattern = "^(table|weight|gene)", all.names = TRUE))
gc()

cat("\n=== Table 2 generation completed ===\n")