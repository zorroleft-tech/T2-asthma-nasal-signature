# =========================================================================
# Script: 08_fig2_replication.R
# Purpose: Generate Figure 2 with replication cohorts analysis
# Author: Zuo
# Date: 2024-01-01
#
# Inputs:
#   - data/derived/locked_weights.csv  # Locked model weights
#   - data/processed_full/expr_*__FULL.rds  # Replication cohort expression data (FULL version)
#   - data/processed_full/expr_*__RAW.rds  # Replication cohort expression data (fallback)
#   - data/processed_diet/*_diet_expr.rds  # Replication cohort expression data (diet version)
#   - data/processed_diet/*_pheno.rds  # Replication cohort phenotype data
#   - R/00_theme_setup.R  # Theme setup script
#   - R/01_scoring_logic.R  # Scoring logic script
#   - R/02_stats_funcs.R  # Statistics functions script
#   - R/03_index_logger.R  # Logging script
#
# Outputs:
#   - output/figures_main/Fig2A_log2FC.pdf  # log2FC plot
#   - output/figures_main/Fig2B_heatmap.pdf  # Heatmap
#   - output/figures_main/Fig2C_attenuation.pdf  # Attenuation effect plot
#   - output/figures_main/Fig2D_paired_correlation.pdf  # Paired correlation plot
#   - output/figures_main/Fig2E_PBMC_boxplot.pdf  # PBMC boxplot
#   - output/logs/08_fig2_replication.log  # Log file
# =========================================================================

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "analysis")) || !dir.exists(file.path(base_dir, "data_preparation"))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define paths using file.path()
processed_dir <- file.path(base_dir, "data", "processed")
derived_dir <- file.path(base_dir, "data", "derived")
figures_dir <- file.path(base_dir, "output", "figures_main")
tables_dir <- file.path(base_dir, "output", "tables_main")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "08_fig2_replication.log")
sink(log_file, append = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Input files:\n")
cat(sprintf("  - %s\n", file.path(derived_dir, "locked_weights.csv")))
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig2A_log2FC.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig2B_heatmap.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig2C_attenuation.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig2D_paired_correlation.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig2E_PBMC_boxplot.pdf")))
cat(sprintf("  - %s\n", log_file))

# Load required packages
library(tidyverse)
library(ggplot2)

# Source theme setup and utility functions
source(file.path(base_dir, "R", "00_theme_setup.R"))
source(file.path(base_dir, "R", "01_scoring_logic.R"))
source(file.path(base_dir, "R", "02_stats_funcs.R"))
source(file.path(base_dir, "R", "03_index_logger.R"))

# Source Figure 2 functions
source(file.path(base_dir, "R", "figures", "fig2A_logFC_concordance.R"))
source(file.path(base_dir, "R", "figures", "fig2B_circular_heatmap.R"))
source(file.path(base_dir, "R", "figures", "fig2C_attenuation.R"))

set_plot_theme()

# Load locked weights
locked_weights_file <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(locked_weights_file)) {
  stop(paste0("Locked weights file not found: ", locked_weights_file))
}

cat("Loading locked weights...\n")
locked_weights <- read_csv(locked_weights_file)

# Read replication contrasts manifest
manifest_path <- file.path(base_dir, "analysis", "02_main_validation", "00_replication_contrasts.yaml")
if (!file.exists(manifest_path)) {
  stop(paste0("Replication contrasts manifest not found: ", manifest_path))
}

library(yaml)
manifest <- read_yaml(manifest_path)

# Expected replication cohorts for Fig2A and Fig2B
expected_cohorts <- c("GSE103166", "GSE115770", "GSE118761", "GSE43696", "GSE123750", "GSE40888", "GSE230048", "GSE152004")

# Validate cohort_id collection in yaml
# Exclude global keys like minimal_stats_report_path
global_keys <- c("minimal_stats_report_path")
manifest_cohorts <- names(manifest)[!names(manifest) %in% global_keys]

# Check if all expected cohorts are present
missing_cohorts <- setdiff(expected_cohorts, manifest_cohorts)
if (length(missing_cohorts) > 0) {
  warning(paste0("Missing expected cohorts in manifest: ", paste(missing_cohorts, collapse = ", ")))
  cat("Continuing with available cohorts...\n")
}

# Check if there are extra cohorts
extra_cohorts <- setdiff(manifest_cohorts, expected_cohorts)
if (length(extra_cohorts) > 0) {
  warning(paste0("Extra cohorts found in manifest: ", paste(extra_cohorts, collapse = ", ")))
  cat("Continuing with expected cohorts...\n")
}

cat("Manifest validation completed\n")

# Define 11-gene signature
signature_genes <- c("CLCA1", "SERPINB2", "CPA3", "CCR3", "HDC", 
                    "MUC5AC", "IL4", "MUC5B", "IL13", "CCL24", "POSTN")

# Extract paths from manifest for Fig2A
  # Use FULL expression matrices instead of diet matrices
  gse152004_expr_path <- file.path(base_dir, "data", "processed_full", "expr_gse152004&GPL11154__FULL.rds")
  # Fallback to RAW if FULL is not available
  if (!file.exists(gse152004_expr_path)) {
    gse152004_expr_path <- file.path(base_dir, "data", "processed_full", "expr_gse152004&GPL11154__RAW.rds")
  }
  gse152004_pheno_path <- file.path(base_dir, "data", "processed_diet", "GSE152004_pheno.rds")
  gse40888_expr_path <- file.path(base_dir, "data", "processed_full", "expr_gse40888&GPL6244__FULL.rds")
  # Fallback to RAW if FULL is not available
  if (!file.exists(gse40888_expr_path)) {
    gse40888_expr_path <- file.path(base_dir, "data", "processed_full", "expr_gse40888&GPL6244__RAW.rds")
  }
  gse40888_pheno_path <- file.path(base_dir, "data", "processed_diet", "GSE40888_pheno.rds")

# Extract minimal stats report path for Fig2C
minimal_stats_report_path <- file.path(base_dir, "data", "processed", "phase3_outputs", "minimal_stats_report.csv")

# Create cohort manifest data frame for Fig2B
cohort_manifest_df <- data.frame(
  cohort_id = names(manifest)[names(manifest) != "minimal_stats_report_path"],
  tissue_group = sapply(manifest[names(manifest) != "minimal_stats_report_path"], function(x) x$tissue_group),
  is_anchor = sapply(manifest[names(manifest) != "minimal_stats_report_path"], function(x) x$is_anchor),
  expr_z_path = sapply(manifest[names(manifest) != "minimal_stats_report_path"], function(x) x$expr_z_path)
)

# Generate Fig 2A: log2FC concordance plot
cat("Generating Fig 2A: log2FC concordance plot...\n")
fig2a_result <- make_fig2A_logfc_concordance(
  gse152004_expr_path = gse152004_expr_path,
  gse152004_pheno_path = gse152004_pheno_path,
  gse40888_expr_path = gse40888_expr_path,
  gse40888_pheno_path = gse40888_pheno_path,
  output_dir = figures_dir,
  tables_dir = tables_dir,
  base_dir = base_dir
)

# Generate Fig 2B: Circular heatmap
cat("Generating Fig 2B: Circular heatmap...\n")
fig2b_result <- make_fig2B_circular_heatmap(
  cohort_manifest = cohort_manifest_df,
  signature_genes = signature_genes,
  output_dir = figures_dir,
  tables_dir = tables_dir,
  base_dir = base_dir
)

# ======================================================================
# Fig2C Contract Builder: Generate Fig2C_cohort_effect_sizes.csv
# ======================================================================

# 强制输出到控制台，确保我们能看到执行日志
cat("\n=== Starting Fig2C Contract Builder ===\n", file = stdout())

# 立即执行Fig2C Contract Builder部分，确保它不会被跳过
cat("Replication-only cohorts: GSE103166, GSE115770, GSE118761, GSE43696, GSE40888, GSE123750, GSE230048\n", file = stdout())
cat("Anchor cohorts: GSE65204, GSE45111, GSE201955\n", file = stdout())

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(yaml)
})

# Define fixed vectors for replication-only and anchor cohorts
replication_only_cohorts <- c("GSE103166", "GSE115770", "GSE118761", "GSE43696", "GSE40888", "GSE123750", "GSE230048")
anchor_cohorts <- c("GSE65204", "GSE45111", "GSE201955")

cat(sprintf("Replication-only cohorts: %s\n", paste(replication_only_cohorts, collapse = ", ")))
cat(sprintf("Anchor cohorts: %s\n", paste(anchor_cohorts, collapse = ", ")))

# Load locked weights
locked_weights_path <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(locked_weights_path)) {
  stop(paste0("Locked weights file not found: ", locked_weights_path))
}
locked_weights <- read_csv(locked_weights_path)
cat("Loaded locked weights\n")

# Define 11-gene signature
signature_genes <- c("CLCA1", "SERPINB2", "CPA3", "CCR3", "HDC", 
                    "MUC5AC", "IL4", "MUC5B", "IL13", "CCL24", "POSTN")
cat(sprintf("Signature genes: %s\n", paste(signature_genes, collapse = ", ")))

# Define tissue group mapping function
map_tissue_group <- function(tissue_raw) {
  if (is.na(tissue_raw)) {
    return("Other")
  }
  tissue_lower <- tolower(tissue_raw)
  if (tissue_lower == "nasal") {
    return("Nasal")
  } else if (tissue_lower %in% c("bronchial", "sputum", "airway", "paired airway")) {
    return("Airway")
  } else if (tissue_lower %in% c("whole blood", "pbmc", "blood", "blood + neutrophils")) {
    return("Blood")
  } else {
    return("Other")
  }
}

# Read Table1_baseline.csv for tissue information
table1_path <- file.path(tables_dir, "Table1_baseline.csv")
if (!file.exists(table1_path)) {
  stop("Table1_baseline.csv not found. Please run the script that generates this file first.")
}

table1_data <- read_csv(table1_path, show_col_types = FALSE)
cat("Loaded Table1_baseline.csv\n")

# Read replication contrasts manifest
manifest_path <- file.path(base_dir, "analysis", "02_main_validation", "00_replication_contrasts.yaml")
if (!file.exists(manifest_path)) {
  stop(paste0("Replication contrasts manifest not found: ", manifest_path))
}
manifest <- read_yaml(manifest_path)
cat("Loaded replication contrasts manifest\n")

# Initialize results lists
effect_sizes_list <- list()
drop_report_list <- list()

# Process each replication-only cohort
for (cohort_id in replication_only_cohorts) {
  cat(sprintf("Processing cohort: %s\n", cohort_id))
  
  # Check if cohort is in anchor list (should not happen, but just in case)
  if (cohort_id %in% anchor_cohorts) {
    stop("Anchor cohort must not be included in Fig2C replication-only contract")
  }
  
  # Initialize status and fail reason
  status <- "PASS"
  fail_reason <- NA
  expected_inputs <- c("expression file", "phenotype file", "contrast definition", "sample groups")
  found_inputs <- c()
  notes <- NA
  
  # Get cohort info from manifest
  if (!cohort_id %in% names(manifest)) {
    status <- "FAIL"
    fail_reason <- "Cohort not found in manifest"
    found_inputs <- c()
  } else {
    cohort_info <- manifest[[cohort_id]]
    
    # Get expression and phenotype paths
    # Directly use diet_expr.rds instead of expr_z_path
    expr_path <- cohort_info$expr_path  # Use diet expression data
    pheno_path <- cohort_info$meta_path
    
    cat(sprintf("  Using diet expression file: %s\n", expr_path))
    cat(sprintf("  File exists: %s\n", file.exists(expr_path)))
    
    if (!file.exists(expr_path)) {
      status <- "FAIL"
      fail_reason <- paste0("Expression file not found: ", expr_path)
      found_inputs <- c()
      next
    } else {
      found_inputs <- c(found_inputs, "expression file")
    }
    
    if (!file.exists(pheno_path)) {
      status <- "FAIL"
      fail_reason <- paste0("Phenotype file not found: ", pheno_path)
    } else {
      found_inputs <- c(found_inputs, "phenotype file")
      
      # Load data
      expr_data <- readRDS(expr_path)
      pheno_data <- readRDS(pheno_path)
      
      # Ensure expr_data is gene x sample
      if (nrow(expr_data) < ncol(expr_data) && any(grepl("^GSM", rownames(expr_data)))) {
        expr_data <- t(expr_data)
      }
      
      # Get sample IDs
      sample_ids <- colnames(expr_data)
      
      # Get case/control labels
      if (!all(c("group_col", "case_levels", "control_levels") %in% names(cohort_info))) {
        status <- "FAIL"
        fail_reason <- "Contrast definition incomplete in manifest"
      } else {
        found_inputs <- c(found_inputs, "contrast definition")
        
        group_col <- cohort_info$group_col
        case_label <- cohort_info$case_levels
        control_label <- cohort_info$control_levels
        
        if (!group_col %in% colnames(pheno_data)) {
          status <- "FAIL"
          fail_reason <- paste0("Group column '", group_col, "' not found in pheno data")
        } else {
          # Determine the correct sample ID column name
          sample_id_col <- if("std_sample_id" %in% colnames(pheno_data)) {
            "std_sample_id"
          } else if("SampleID" %in% colnames(pheno_data)) {
            "SampleID"
          } else {
            # Try to find a column that looks like sample ID
            sample_id_col <- colnames(pheno_data)[grepl("sample|id|GSM", colnames(pheno_data), ignore.case = TRUE)][1]
            if(is.na(sample_id_col)) {
              status <- "FAIL"
              fail_reason <- "No sample ID column found in pheno data"
              next
            }
            sample_id_col
          }
          
          # Filter pheno data to sample IDs in expr data
          pheno_data <- pheno_data %>% filter(!!sym(sample_id_col) %in% sample_ids)
          
          # Ensure order matches
          pheno_data <- pheno_data[match(sample_ids, pheno_data[[sample_id_col]]), ]
          
          # Create case/control vector
          groups <- pheno_data[[group_col]]
          case_indices <- groups == case_label
          control_indices <- groups == control_label
          
          N_case <- sum(case_indices)
          N_control <- sum(control_indices)
          N <- N_case + N_control
          
          if (N_case == 0 || N_control == 0) {
            status <- "FAIL"
            fail_reason <- "No cases or controls found"
          } else {
            found_inputs <- c(found_inputs, "sample groups")
            
            # Calculate locked score
            avail_genes <- intersect(signature_genes, rownames(expr_data))
            coverage_count <- length(avail_genes)
            coverage_missing <- length(signature_genes) - coverage_count
            
            if (coverage_count == 0) {
              status <- "FAIL"
              fail_reason <- "No signature genes found"
            } else {
              # Filter weights to available genes
              cohort_weights <- locked_weights %>% filter(Gene %in% avail_genes)
              
              # Calculate score
              score_matrix <- expr_data[avail_genes, ]
              weights_vector <- cohort_weights$Weight[match(avail_genes, cohort_weights$Gene)]
              scores <- as.numeric(crossprod(weights_vector, score_matrix))
              
              # Z-score normalize within cohort
              score_z <- as.numeric(scale(scores))
              
              # Calculate Cohen's d using raw scores (d_raw)
              case_scores_raw <- scores[case_indices]
              control_scores_raw <- scores[control_indices]
              
              mean_case_raw <- mean(case_scores_raw)
              mean_control_raw <- mean(control_scores_raw)
              mean_diff_raw <- mean_case_raw - mean_control_raw
              sd_pooled_raw <- sqrt(((N_case - 1) * var(case_scores_raw) + (N_control - 1) * var(control_scores_raw)) / (N - 2))
              d_raw <- mean_diff_raw / sd_pooled_raw
              
              # Calculate Cohen's d using Z-scored scores (d_z)
              case_scores_z <- score_z[case_indices]
              control_scores_z <- score_z[control_indices]
              
              mean_case_z <- mean(case_scores_z)
              mean_control_z <- mean(control_scores_z)
              mean_diff_z <- mean_case_z - mean_control_z
              sd_pooled_z <- sqrt(((N_case - 1) * var(case_scores_z) + (N_control - 1) * var(control_scores_z)) / (N - 2))
              d_z <- mean_diff_z / sd_pooled_z
              
              # Calculate Wilcoxon p-value (using Z-scored scores)
              wilcoxon_p <- wilcox.test(case_scores_z, control_scores_z)$p.value
              
              # Get tissue information from Table1
              table1_row <- table1_data[table1_data$`GSE ID` == cohort_id, ]
              if (nrow(table1_row) > 0) {
                tissue_raw <- table1_row$Tissue[1]
                tissue_group <- map_tissue_group(tissue_raw)
              } else {
                # Fallback to manifest if not found in Table1
                tissue_raw <- cohort_info$tissue_group
                tissue_group <- map_tissue_group(tissue_raw)
              }
              cat(sprintf("  Tissue from Table1: %s -> %s\n", tissue_raw, tissue_group))
            }
          }
        }
      }
    }
  }
  
  # Create drop report entry
  drop_report_entry <- tibble(
    cohort_id = cohort_id,
    status = status,
    fail_reason = fail_reason,
    expected_inputs = paste(expected_inputs, collapse = "; "),
    found_inputs = paste(found_inputs, collapse = "; "),
    notes = notes
  )
  drop_report_list[[cohort_id]] <- drop_report_entry
  
  # Create effect size entry
  if (status == "PASS") {
    effect_size_entry <- tibble(
      cohort_id = cohort_id,
      n_total = N,
      n_case = N_case,
      n_control = N_control,
      mean_diff_raw = mean_diff_raw,
      mean_diff_z = mean_diff_z,
      d_raw = d_raw,
      d_z = d_z,
      cohens_d = d_z,  # Use d_z as cohens_d for consistency
      p_value = wilcoxon_p,
      coverage_count = coverage_count,
      coverage_missing = coverage_missing,
      tissue_group = tissue_group,
      status = status,
      fail_reason = fail_reason,
      notes = notes
    )
  } else {
    # Create entry with NA values for failed cohorts
    effect_size_entry <- tibble(
      cohort_id = cohort_id,
      n_total = NA,
      n_case = NA,
      n_control = NA,
      mean_diff_raw = NA,
      mean_diff_z = NA,
      d_raw = NA,
      d_z = NA,
      cohens_d = NA,
      p_value = NA,
      coverage_count = NA,
      coverage_missing = NA,
      tissue_group = NA,
      status = status,
      fail_reason = fail_reason,
      notes = notes
    )
  }
  
  # Print the structure to debug
  cat(sprintf("  Effect size entry structure: %s\n", paste(colnames(effect_size_entry), collapse = ", ")))
  
  effect_sizes_list[[cohort_id]] <- effect_size_entry
  cat(sprintf("  Status: %s\n", status))
  if (status == "FAIL") {
    cat(sprintf("  Fail reason: %s\n", fail_reason))
  }
}

# Combine all effect sizes
effect_sizes_df <- bind_rows(effect_sizes_list)
cat(sprintf("Combined effect sizes for %d cohorts\n", nrow(effect_sizes_df)))

# Combine all drop reports
drop_report_df <- bind_rows(drop_report_list)
cat(sprintf("Combined drop reports for %d cohorts\n", nrow(drop_report_df)))

# Validate that all replication-only cohorts are present
if (!setequal(unique(effect_sizes_df$cohort_id), replication_only_cohorts)) {
  cat("ERROR: Missing or extra cohorts in output\n")
  cat(sprintf("Expected: %s\n", paste(replication_only_cohorts, collapse = ", ")))
  cat(sprintf("Actual: %s\n", paste(unique(effect_sizes_df$cohort_id), collapse = ", ")))
  stop("All replication-only cohorts must be present in the output")
} else {
  cat("Validation passed: All replication-only cohorts are present\n")
}

# Validate that no anchor cohorts are present
if (any(effect_sizes_df$cohort_id %in% anchor_cohorts)) {
  cat("ERROR: Anchor cohorts found in output\n")
  stop("Anchor cohort must not be included in Fig2C replication-only contract")
} else {
  cat("Validation passed: No anchor cohorts in output\n")
}

# Save the contract table
contract_output_path <- file.path(tables_dir, "Fig2C_cohort_effect_sizes.csv")
write_csv(effect_sizes_df, contract_output_path)
cat(paste0("Fig2C_cohort_effect_sizes.csv saved to: ", contract_output_path, "\n"))

# Save the drop report
drop_report_output_path <- file.path(tables_dir, "Fig2C_cohort_effect_sizes_drop_report.csv")
write_csv(drop_report_df, drop_report_output_path)
cat(paste0("Fig2C_cohort_effect_sizes_drop_report.csv saved to: ", drop_report_output_path, "\n"))

# Log the outputs
log_output(
  file_path = contract_output_path,
  file_type = "table",
  figure_id = "Fig2C_cohort_effect_sizes",
  description = "Effect sizes (Cohen's d) for Fig2C across 7 replication-only cohorts",
  script_name = "08_fig2_replication.R"
)

log_output(
  file_path = drop_report_output_path,
  file_type = "table",
  figure_id = "Fig2C_cohort_effect_sizes_drop_report",
  description = "Drop report for Fig2C replication-only cohorts",
  script_name = "08_fig2_replication.R"
)

cat("\n=== Fig2C Contract Builder completed ===\n")

# Memory cleanup
# Preserve output path variables
output_paths <- c("drop_report_output_path", "contract_output_path")
variables_to_remove <- setdiff(ls(pattern = "^(effect_|cohort_|expr_|score_|locked_|weight_|gse_|minimal_|signature_|table1_|manifest_)", all.names = TRUE), output_paths)
if (length(variables_to_remove) > 0) {
  rm(list = variables_to_remove)
}
gc()

# Verify the output files
cat("\n=== Verifying output files ===\n")
if (file.exists(contract_output_path)) {
  cat("✓ Fig2C_cohort_effect_sizes.csv exists\n")
  verify_df <- read_csv(contract_output_path, show_col_types = FALSE)
  cat(sprintf("  Rows: %d\n", nrow(verify_df)))
  cat(sprintf("  Columns: %s\n", paste(colnames(verify_df), collapse = ", ")))
  cat(sprintf("  Cohorts: %s\n", paste(verify_df$cohort_id, collapse = ", ")))
} else {
  cat("✗ Fig2C_cohort_effect_sizes.csv does not exist\n")
}

if (file.exists(drop_report_output_path)) {
  cat("✓ Fig2C_cohort_effect_sizes_drop_report.csv exists\n")
  verify_drop_df <- read_csv(drop_report_output_path, show_col_types = FALSE)
  cat(sprintf("  Rows: %d\n", nrow(verify_drop_df)))
  cat(sprintf("  Columns: %s\n", paste(colnames(verify_drop_df), collapse = ", ")))
} else {
  cat("✗ Fig2C_cohort_effect_sizes_drop_report.csv does not exist\n")
}

# Generate Fig 2C: Attenuation effect
cat("Generating Fig 2C: Attenuation effect...\n")
fig2c_result <- make_fig2C_attenuation(
  output_dir = figures_dir,
  tables_dir = tables_dir,
  base_dir = base_dir
)

# Log outputs
log_output(
  file_path = fig2a_result$output_pdf,
  file_type = "figure",
  figure_id = "Fig 2A",
  description = "log2FC concordance plot of 11-gene signature across tissues",
  script_name = "08_fig2_replication.R",
  resolution = "300 dpi"
)

log_output(
  file_path = fig2b_result$output_pdf,
  file_type = "figure",
  figure_id = "Fig 2B",
  description = "Circular heatmap of 11-gene signature expression across cohorts",
  script_name = "08_fig2_replication.R",
  resolution = "300 dpi"
)

log_output(
  file_path = fig2c_result$output_pdf,
  file_type = "figure",
  figure_id = "Fig 2C",
  description = "Tissue-specific attenuation of 11-gene signature effect size",
  script_name = "08_fig2_replication.R",
  resolution = "300 dpi"
)

cat("\n=== Figure 2 generation completed ===\n")

# Memory cleanup
rm(list = ls(pattern = "^(fig|cohort|expr|score|locked|weight|gse|minimal|signature)", all.names = TRUE))
gc()

# Close logging
sink()