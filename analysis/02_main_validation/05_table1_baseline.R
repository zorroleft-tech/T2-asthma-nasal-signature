#!/usr/bin/env Rscript

# =========================================================================
# Script: 05_table1_baseline.R
# Purpose: Generate Table 1 with baseline characteristics of all 11 bulk cohorts
# Author: Zuo
# Date: 2026-02-27
#
# Inputs:
#   - data/processed_diet/*_pheno.rds  # Phenotype data for each cohort
#
# Outputs:
#   - output/tables_main/Table1_baseline.csv  # Baseline characteristics table
#   - output/tables_main/Table1_baseline.xlsx  # Excel version (optional)
#   - output/logs/05_table1_baseline.log  # Log file
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
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
tables_dir <- file.path(base_dir, "output", "tables_main")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "05_table1_baseline.log")

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
cat("Input directory: ", processed_diet_dir, "\n")
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(tables_dir, "Table1_baseline.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Table1_baseline.xlsx")))
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

# Define cohorts to include (11 bulk cohorts)
cohorts <- c(
  "GSE152004",  # Discovery cohort
  "GSE65204",    # Anchor cohort (IgE)
  "GSE201955",   # Anchor cohort (FeNO)
  "GSE45111",    # Anchor cohort (Sputum eos)
  "GSE103166",   # Replication cohort
  "GSE115770",   # Replication cohort
  "GSE118761",   # Replication cohort
  "GSE123750",   # Replication cohort
  "GSE230048",   # Replication cohort
  "GSE40888",    # Replication cohort
  "GSE43696"     # Replication cohort
)

# Cohort metadata map (authoritative text source)
cohort_metadata_map <- tibble::tribble(
  ~`GSE ID`,     ~Role,              ~Tissue,           ~`Platform (GPL)`,                    ~`N (Total)`,            ~`No. of Cases`, ~`No. of Controls`, ~Age,                                                                 ~Sex,
  "GSE152004",   "Discovery",        "Nasal",           "GPL11154 (RNA-seq)",                  "695",                   "347",          "348",              "8–21 years (Puerto Rican pediatric cohort; GALA II)",                "~54% male (from Sajuthi et al. 2022)",
  "GSE65204",    "Clinical anchor",  "Nasal",           "GPL14550",                            "46",                    "23",           "23",               "9-12 years (mean 11.0 ± 0.9)",                                       "51% female (35/69)",
  
  # Anchor cohorts missing from current Table 1 in manuscript (must be added)
  "GSE201955",   "Clinical anchor",  "Bronchial",       "GPL20301 (RNA-seq)",                  "74",                    NA,              NA,                  "21-58 years",                                                       "~51.4% female (38/74)",
  "GSE45111",    "Clinical anchor",  "Induced sputum",  "GPL6104 (microarray)",                "47",                    NA,              NA,                  "21-78 years",                                                       "~63.8% female (30/47)",
  
  "GSE103166",   "Replication-only", "Nasal",           "GPL23961 (Illumina HumanHT-12 v4.0)", "106",                   "41",           "65",               "0–16 years (mean ~4.9 years; WE: 4.35 ± 3.31 y, Conv: 4.48 ± 1.62 y, Ctrl: 6.40 ± 4.41 y)", "53.8% male (57/106; from Khoo et al. 2019 Table I)",
  "GSE115770",   "Replication-only", "Nasal",           "RNA-seq",                             "222 (samples)^b",       "132",          "90",               "Median 11 years (IQR: 8–13)",                                       "49.1% female (52/106 participants)^b",
  "GSE118761",   "Replication-only", "Paired airway",   "GPL11154 (RNA-seq)",                  "104",                   "17",           "87",               "2.0-16.4 years (children subset; median 3.78; mean 4.99)",           "Approximately balanced (~1:1)",
  "GSE123750",   "Replication-only", "Blood",           "GPL13158",                            "216",                   "91",           "125",              "1-17 years",                                                         "40% female (96/216)",
  "GSE230048",   "Replication-only", "Blood + neutrophils", "GPL16791",                         "88",                    "44",           "44",               "Mean ~7.8 years (NR: 7.6 ± 2.5 y; R: 7.9 ± 2.2 y; from Tsai et al. 2023 Table 1, TCCAS cohort)", "67.3% male (37/55; NR: 17/26, R: 20/29; from Tsai et al. 2023 Table 1)",
  "GSE40888",    "Replication-only", "PBMC",            "GPL6244",                             "105",                   "39",           "66",               "Median 9.5 years (IQR: 6.5–11.4) for AA; median 7.7 years (IQR: 6.6–9.3) for NA; range 4–15 years (from Raedler et al. 2015 Table I)", "63.1% male (AA+NA combined: 58/92; AA: 50/74, 67.6%; NA: 8/18, 44.4%; from Raedler et al. 2015 Table I)",
  "GSE43696",    "Replication-only", "Bronchial",       "GPL6480 (Agilent-014850 4×44K)",     "108",                   "12",           "96",               "18.2-62.1 years (mean 37.1 ± 12.4)",                                 "68.5% female (74/108)"
)

# Function to find pheno file for a cohort
find_pheno_file <- function(gse_id, base_dir) {
  # Try different naming patterns
  patterns <- c(
    paste0(gse_id, "_pheno.rds"),
    paste0(gse_id, "_diet_pheno.rds"),
    paste0(gse_id, "_pheno.csv")
  )
  
  for (pattern in patterns) {
    file_path <- file.path(base_dir, pattern)
    if (file.exists(file_path)) {
      return(file_path)
    }
  }
  return(NULL)
}

# Process each cohort to check pheno files
cat("\n=== Checking pheno files ===\n")
pheno_files <- list()
pheno_n <- list()
for (cohort in cohorts) {
  pheno_file <- find_pheno_file(cohort, processed_diet_dir)
  if (!is.null(pheno_file)) {
    pheno_files[[cohort]] <- pheno_file
    # Read pheno to get nrow
    pheno <- readRDS(pheno_file)
    n_rows <- nrow(pheno)
    pheno_n[[cohort]] <- n_rows
    cat(sprintf("%s: %s (nrow = %d)\n", cohort, basename(pheno_file), n_rows))
  } else {
    cat(sprintf("%s: No pheno file found\n", cohort))
  }
}

# Validate GSE230048 N=88
cat("\n=== Validating GSE230048 N=88 ===\n")
gse230048_blood_file <- file.path(processed_diet_dir, "GSE230048_pheno_blood.rds")
gse230048_neutro_file <- file.path(processed_diet_dir, "GSE230048_pheno_neutrophil.rds")
gse230048_combined_file <- file.path(processed_diet_dir, "GSE230048_pheno.rds")

if (file.exists(gse230048_blood_file) && file.exists(gse230048_neutro_file)) {
  # Method 1: Dual pheno files
  blood_pheno <- readRDS(gse230048_blood_file)
  neutro_pheno <- readRDS(gse230048_neutro_file)
  blood_n <- nrow(blood_pheno)
  neutro_n <- nrow(neutro_pheno)
  total_n <- blood_n + neutro_n
  
  stopifnot(blood_n == 46, neutro_n == 42, total_n == 88)
  cat(sprintf("GSE230048: blood_n=%d + neutrophil_n=%d => table_n=88\n", blood_n, neutro_n))
  pheno_n[['GSE230048']] <- total_n
} else if (file.exists(gse230048_combined_file)) {
  # Method 2: Single pheno file with compartment
  combined_pheno <- readRDS(gse230048_combined_file)
  if ('compartment' %in% names(combined_pheno)) {
    blood_n <- sum(combined_pheno$compartment == 'blood')
    neutro_n <- sum(combined_pheno$compartment == 'neutrophil')
    total_n <- blood_n + neutro_n
    
    stopifnot(blood_n == 46, neutro_n == 42, total_n == 88)
    cat(sprintf("GSE230048: blood_n=%d + neutrophil_n=%d => table_n=88\n", blood_n, neutro_n))
    pheno_n[['GSE230048']] <- total_n
  } else {
    stop("GSE230048 requires blood + neutrophil pheno to justify N=88")
  }
} else {
  stop("GSE230048 requires blood + neutrophil pheno to justify N=88")
}

# Generate Table 1 from metadata map
cat("\n=== Generating Table 1 ===\n")
table1 <- cohort_metadata_map

# Validate row count
if (nrow(table1) != 11) {
  stop(sprintf("Expected 11 rows in Table 1, but got %d", nrow(table1)), call. = FALSE)
}

# Validate anchor cohorts
anchor_cohorts <- c("GSE65204", "GSE201955", "GSE45111")
if (!all(anchor_cohorts %in% table1$`GSE ID`)) {
  missing_anchors <- setdiff(anchor_cohorts, table1$`GSE ID`)
  stop(sprintf("Missing anchor cohorts: %s", paste(missing_anchors, collapse = ", ")), call. = FALSE)
}

# Consistency validation block
cat("\n=== Consistency Validation ===\n")

# Helper function to extract numeric N from Table 1
get_table_n <- function(gse_id) {
  if (gse_id == "GSE230048") {
    return(88)  # Fixed value
  } else if (gse_id == "GSE115770") {
    return(222)  # Fixed value
  } else if (gse_id == "GSE65204") {
    return(46)  # Fixed value
  } else if (gse_id == "GSE152004") {
    return(695)  # Fixed value
  } else if (gse_id == "GSE201955") {
    return(74)  # Fixed value
  } else if (gse_id == "GSE45111") {
    return(47)  # Fixed value
  } else if (gse_id == "GSE103166") {
    return(106)  # Fixed value
  } else if (gse_id == "GSE118761") {
    return(104)  # Fixed value
  } else if (gse_id == "GSE123750") {
    return(216)  # Fixed value
  } else if (gse_id == "GSE40888") {
    return(105)  # Fixed value
  } else if (gse_id == "GSE43696") {
    return(108)  # Fixed value
  } else {
    stop("Unknown cohort: ", gse_id)
  }
}

# Strong validation (must match)
strong_validation_cohorts <- c(
  "GSE201955", "GSE45111", "GSE103166", 
  "GSE118761", "GSE123750", "GSE40888", "GSE43696"
)

strong_validation_passed <- c()
for (cohort in strong_validation_cohorts) {
  table_n <- get_table_n(cohort)
  pheno_n_val <- pheno_n[[cohort]]
  
  if (table_n != pheno_n_val) {
    stop(sprintf("Pheno nrow does not match Table 1 N for %s: pheno_n=%d, table_n=%d", cohort, pheno_n_val, table_n))
  }
  strong_validation_passed <- c(strong_validation_passed, cohort)
  cat(sprintf("%s: pheno_n=%d, table_n=%d (PASS)\n", cohort, pheno_n_val, table_n))
}

# Exception cohorts (explain in log)
exception_cohorts <- c("GSE65204", "GSE115770")
for (cohort in exception_cohorts) {
  table_n <- get_table_n(cohort)
  pheno_n_val <- pheno_n[[cohort]]
  
  if (cohort == "GSE65204") {
    cat(sprintf("%s: pheno_n=%d, table_n=%d (IgE extreme tertiles; middle tertile excluded)\n", cohort, pheno_n_val, table_n))
  } else if (cohort == "GSE115770") {
    cat(sprintf("%s: pheno_n=%d, table_n=%d (sample-level transcriptomes; pheno contains longitudinal records)\n", cohort, pheno_n_val, table_n))
  }
}

# Save Table 1 as CSV
csv_output <- file.path(tables_dir, "Table1_baseline.csv")
write_csv(table1, csv_output)
cat(sprintf("Saved Table 1 CSV: %s\n", basename(csv_output)))

# Save Table 1 as Excel if openxlsx is available
if (create_xlsx) {
  xlsx_output <- file.path(tables_dir, "Table1_baseline.xlsx")
  write.xlsx(table1, xlsx_output, rowNames = FALSE)
  cat(sprintf("Saved Table 1 Excel: %s\n", basename(xlsx_output)))
}

# Print Table 1 summary
cat("\n=== Table 1 Summary ===\n")
cat(sprintf("Total rows: %d\n", nrow(table1)))
cat("Cohorts included:\n")
print(table1$`GSE ID`)

# Log validation summary
cat("\n=== Validation Summary ===\n")
cat("Strong validation passed for:", paste(strong_validation_passed, collapse = ", "), "\n")
cat("Exceptions explained: GSE65204, GSE115770, GSE230048\n")

# Memory cleanup
rm(list = ls(pattern = "^(table|cohort|pheno|output)", all.names = TRUE))
gc()

cat("\n=== Table 1 generation completed ===\n")