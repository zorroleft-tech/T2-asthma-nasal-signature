# =========================================================================
# Script: 14_supp_tableS9_platform.R
# Purpose: Generate Supplementary Table S9 with gene coverage across platforms
# Author: Zuo
# Date: 2026-02-27
#
# Inputs:
#   - data/derived/locked_weights.csv  # Locked 11 genes
#   - data/processed_diet/*_diet_expr.rds  # Diet expression matrices
#
# Outputs:
#   - output/supplement/TableS9_platform_long.csv  # Long format table
#   - output/supplement/TableS9_platform_wide.csv  # Wide format table
#   - output/logs/14_supp_tableS9_platform.log  # Log file
# =========================================================================

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "analysis")) || !dir.exists(file.path(base_dir, "data_preparation"))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define paths using file.path()
data_derived_dir <- file.path(base_dir, "data", "derived")
data_processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
supplement_dir <- file.path(base_dir, "output", "supplement")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "14_supp_tableS9_platform.log")
sink(log_file, append = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Input files:\n")
cat(sprintf("  - %s\n", file.path(data_derived_dir, "locked_weights.csv")))
cat(sprintf("  - %s/*_diet_expr.rds\n", data_processed_diet_dir))
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS9_platform_long.csv")))
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS9_platform_wide.csv")))
cat(sprintf("  - %s\n", log_file))

# Load required packages
library(tidyverse)

# =========================================================================
# Step 1: Read locked genes
# =========================================================================

cat("\n[Step 1] Reading locked genes...\n")

# Read locked weights file
locked_weights_file <- file.path(data_derived_dir, "locked_weights.csv")
if (!file.exists(locked_weights_file)) {
  stop("Locked weights file not found: ", locked_weights_file, call. = FALSE)
}

locked_genes_df <- read.csv(locked_weights_file)

# Find gene column (case-insensitive)
gene_col <- grep("gene", colnames(locked_genes_df), ignore.case = TRUE)
if (length(gene_col) > 0) {
  locked_genes <- as.character(locked_genes_df[, gene_col[1]])
  cat(paste0("Using gene column: ", colnames(locked_genes_df)[gene_col[1]], "\n"))
} else {
  # Fallback to first column
  locked_genes <- as.character(locked_genes_df[, 1])
  cat("No gene column found, using first column\n")
}

# Validate locked genes
if (length(locked_genes) != 11) {
  stop("Locked genes count is not 11: ", length(locked_genes), call. = FALSE)
}

if (!all(c("CLCA1", "POSTN") %in% locked_genes)) {
  stop("Locked genes must include CLCA1 and POSTN", call. = FALSE)
}

cat(paste0("Locked genes (11): ", paste(locked_genes, collapse = ", "), "\n"))

# =========================================================================
# Step 2: Define cohorts to process
# =========================================================================

cat("\n[Step 2] Defining cohorts to process...\n")

# Define cohorts with their roles
cohorts <- list(
  # Anchors
  list(gse = "GSE65204", role = "Anchor"),
  list(gse = "GSE201955", role = "Anchor"),
  list(gse = "GSE45111", role = "Anchor"),
  # Replication-only
  list(gse = "GSE103166", role = "Replication-only"),
  list(gse = "GSE40888", role = "Replication-only"),
  list(gse = "GSE43696", role = "Replication-only"),
  list(gse = "GSE115770", role = "Replication-only"),
  list(gse = "GSE118761", role = "Replication-only"),
  list(gse = "GSE123750", role = "Replication-only"),
  list(gse = "GSE230048", role = "Replication-only")
)

# =========================================================================
# Step 3: Process each cohort
# =========================================================================

cat("\n[Step 3] Processing each cohort...\n")

# Initialize results lists
long_data <- list()
wide_data <- list()

for (cohort in cohorts) {
  gse <- cohort$gse
  role <- cohort$role
  
  cat(paste0("\nProcessing cohort: ", gse, " (", role, ")\n"))
  
  # Read diet expression matrix
  diet_expr_file <- file.path(data_processed_diet_dir, paste0(gse, "_diet_expr.rds"))
  if (!file.exists(diet_expr_file)) {
    stop("Diet expression file not found: ", diet_expr_file, call. = FALSE)
  }
  
  diet_expr <- readRDS(diet_expr_file)
  present_genes <- rownames(diet_expr)
  
  cat(paste0("  Diet expression matrix dimensions: ", nrow(diet_expr), " x ", ncol(diet_expr), "\n"))
  
  # Extract platform from processed_full filename if possible
  platform <- NA
  processed_full_files <- list.files(file.path(base_dir, "data", "processed_full"), pattern = paste0("expr_.*", tolower(gse), ".*__FULL.rds"), full.names = TRUE)
  if (length(processed_full_files) > 0) {
    filename <- basename(processed_full_files[1])
    gpl_match <- str_extract(filename, "GPL\\d+")
    if (!is.na(gpl_match)) {
      platform <- gpl_match
    }
  }
  
  # Calculate present status for each locked gene
  present_status <- locked_genes %in% present_genes
  coverage <- sum(present_status)
  missing_genes <- locked_genes[!present_status]
  missing_genes_str <- if (length(missing_genes) > 0) paste(missing_genes, collapse = ", ") else "None"
  
  cat(paste0("  Coverage: ", coverage, "/11\n"))
  cat(paste0("  Missing genes: ", missing_genes_str, "\n"))
  
  # Build long format data
  long_cohort_data <- data.frame(
    Dataset = gse,
    Gene = locked_genes,
    Present = as.integer(present_status),
    Source = "diet_expr",
    Notes = NA
  )
  long_data[[gse]] <- long_cohort_data
  
  # Build wide format data
  wide_cohort_data <- data.frame(
    Dataset = gse,
    Platform = platform,
    Role = role
  )
  
  # Add columns for each locked gene
  for (gene in locked_genes) {
    wide_cohort_data[[gene]] <- as.integer(gene %in% present_genes)
  }
  
  # Add coverage and missing genes columns
  wide_cohort_data$Coverage_X_of_11 <- coverage
  wide_cohort_data$Missing_Genes <- missing_genes_str
  
  wide_data[[gse]] <- wide_cohort_data
}

# =========================================================================
# Step 4: Generate and save tables
# =========================================================================

cat("\n[Step 4] Generating and saving tables...\n")

# Combine long data
long_table <- bind_rows(long_data)
long_output <- file.path(supplement_dir, "TableS9_platform_long.csv")
write_csv(long_table, long_output)
cat(paste0("Saved long format table: ", basename(long_output), " (", nrow(long_table), " rows)\n"))

# Combine wide data in the order of cohorts
wide_table <- bind_rows(wide_data[names(wide_data)])
# Reorder columns to match locked genes order
wide_cols <- c("Dataset", "Platform", "Role", locked_genes, "Coverage_X_of_11", "Missing_Genes")
wide_table <- wide_table[, wide_cols]

wide_output <- file.path(supplement_dir, "TableS9_platform_wide.csv")
write_csv(wide_table, wide_output)
cat(paste0("Saved wide format table: ", basename(wide_output), " (", nrow(wide_table), " rows)\n"))

# =========================================================================
# Step 5: Final validation
# =========================================================================

cat("\n[Step 5] Final validation...\n")

# Validate wide table rows
if (nrow(wide_table) != 10) {
  stop("Wide table should have 10 rows (3 anchors + 7 replication cohorts)", call. = FALSE)
}

# Validate coverage values
coverage_values <- wide_table$Coverage_X_of_11
if (any(coverage_values < 0) || any(coverage_values > 11)) {
  stop("Coverage values must be between 0 and 11", call. = FALSE)
}

# Validate specific cohorts
if (wide_table[wide_table$Dataset == "GSE65204", "Coverage_X_of_11"] != 6) {
  stop("GSE65204 should have coverage 6/11", call. = FALSE)
}

if (wide_table[wide_table$Dataset == "GSE40888", "Coverage_X_of_11"] != 10) {
  stop("GSE40888 should have coverage 10/11", call. = FALSE)
}

cat("✓ All validations passed\n")

# =========================================================================
# Final report
# =========================================================================

cat("\n=== Supplementary Table S9 generation completed ===\n")
cat(paste0("Long format: ", nrow(long_table), " rows\n"))
cat(paste0("Wide format: ", nrow(wide_table), " rows\n"))
cat("\nKey findings:\n")
for (i in 1:nrow(wide_table)) {
  row <- wide_table[i, ]
  cat(paste0("  ", row$Dataset, ": ", row$Coverage_X_of_11, "/11, Missing: ", row$Missing_Genes, "\n"))
}

# Memory cleanup
rm(list = ls())
gc()

# Close logging
sink()
