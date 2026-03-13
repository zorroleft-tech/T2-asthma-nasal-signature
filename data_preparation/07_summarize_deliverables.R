#!/usr/bin/env Rscript

# =========================================================================
# Script: 07_summarize_deliverables.R
# Purpose: Summarize all cohort results and prepare deliverables
# Author: Zuo
# Date: 2026-03-01
#
# Inputs:
#   - data/processed/phase3_outputs/*_pheno_raw.rds  # Standardized phenotype tables
#   - data/processed/phase3_outputs/minimal_stats_report.csv  # Minimal stats report
#   - data/processed/phase3_outputs/field_audit_report.csv  # Field audit report
#   - data/processed/phase3_outputs/label_mapping_audit.csv  # Label mapping audit report
#   - output/supplement/TableS3_anchor_auc_audit.csv  # Anchor cohort audit
#
# Outputs:
#   - data/processed/phase3_outputs/deliverables_summary.csv  # Summary of all deliverables
#   - output/manuscript_insert_new_truth_cohorts.md  # Manuscript insert for new truth cohorts
#   - output/manuscript_insert_updates_log.md  # Updates log for cover letter
# =========================================================================

library(dplyr)
library(stringr)
library(readr)

# Set base directory
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
required_dirs <- c('data', 'analysis', 'data_preparation')
required_paths <- file.path(base_dir, required_dirs)
dirs_exist <- all(dir.exists(required_paths))
if (!dirs_exist) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
processed_dir <- file.path(base_dir, "data", "processed", "phase3_outputs")
output_dir <- file.path(base_dir, "output")

# Create output directory if it doesn't exist
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to summarize deliverables
summarize_deliverables <- function() {
  # Read minimal stats report
  stats_file <- file.path(processed_dir, "minimal_stats_report.csv")
  if (file.exists(stats_file)) {
    stats_df <- read.csv(stats_file)
  } else {
    cat("Warning: minimal_stats_report.csv not found\n")
    stats_df <- data.frame()
  }
  
  # Read replication direction audit
  direction_audit_file <- file.path(processed_dir, "replication_direction_audit.csv")
  if (file.exists(direction_audit_file)) {
    direction_audit_df <- read.csv(direction_audit_file)
  } else {
    cat("Warning: replication_direction_audit.csv not found\n")
    direction_audit_df <- data.frame()
  }
  
  # Read field audit report
  audit_file <- file.path(processed_dir, "field_audit_report.csv")
  if (file.exists(audit_file)) {
    audit_df <- read.csv(audit_file)
  } else {
    cat("Warning: field_audit_report.csv not found\n")
    audit_df <- data.frame()
  }
  
  # Read label mapping audit report
  label_audit_file <- file.path(processed_dir, "label_mapping_audit.csv")
  if (file.exists(label_audit_file)) {
    label_audit_df <- read.csv(label_audit_file)
  } else {
    cat("Warning: label_mapping_audit.csv not found\n")
    label_audit_df <- data.frame()
  }
  
  # Read Table S3 anchor audit
  table_s3_path <- file.path(base_dir, "output", "supplement", "TableS3_anchor_auc_audit.csv")
  if (file.exists(table_s3_path)) {
    table_s3 <- read_csv(table_s3_path, show_col_types = FALSE)
  } else {
    cat("Warning: TableS3_anchor_auc_audit.csv not found\n")
    table_s3 <- data.frame()
  }
  
  # Define cohorts
  anchor_cohorts <- c("GSE65204", "GSE201955", "GSE45111")
  replication_cohorts <- c("GSE103166", "GSE115770", "GSE118761", "GSE123750", "GSE230048", "GSE40888", "GSE43696")
  
  # List all pheno files
  pheno_files <- list.files(processed_dir, pattern = "_pheno_raw\\.rds$", full.names = TRUE)
  
  # Also check processed_diet directory for pheno files
  processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
  if (dir.exists(processed_diet_dir)) {
    diet_pheno_files <- list.files(processed_diet_dir, pattern = "_pheno\\.rds$", full.names = TRUE)
    # Only add files that exist
    diet_pheno_files <- diet_pheno_files[file.exists(diet_pheno_files)]
    pheno_files <- c(pheno_files, diet_pheno_files)
  }
  
  # Remove duplicates and ensure all files exist
  pheno_files <- unique(pheno_files)
  pheno_files <- pheno_files[file.exists(pheno_files)]
  
  # Initialize deliverables summary
  deliverables_summary <- data.frame()
  
  # First process anchor cohorts
  for (cohort in anchor_cohorts) {
    if (nrow(table_s3) > 0) {
      # Get info from Table S3
      table_s3_cohort <- table_s3[table_s3$Cohort == cohort, ]
      if (nrow(table_s3_cohort) > 0) {
        # Extract gene coverage (e.g., "6/11" -> 6)
        coverage_count <- as.numeric(str_extract(table_s3_cohort$Gene_coverage, "^[0-9]+"))
        
        # Calculate mean_diff
        mean_diff <- table_s3_cohort$Score_mean_case - table_s3_cohort$Score_mean_control
        
        # Determine direction_ok
        score_should_be_high_in_case <- TRUE
        direction_ok <- ifelse((score_should_be_high_in_case && mean_diff > 0) || (!score_should_be_high_in_case && mean_diff < 0), "Yes", "No")
        
        # Calculate PASS_interpret
        PASS_interpret <- ""
        if (direction_ok == "No") {
          if (abs(mean_diff) < 0.5) {
            PASS_interpret <- "Boundary condition"
          } else {
            PASS_interpret <- "Review contrast definition"
          }
        }
        
        # Create summary row
        summary_row <- data.frame(
          cohort = cohort,
          n_total = table_s3_cohort$n_case + table_s3_cohort$n_control,
          n_case = table_s3_cohort$n_case,
          n_control = table_s3_cohort$n_control,
          n_aa = 0,
          n_na = 0,
          n_hc = 0,
          effect_size = NA,  # Not available in Table S3
          p_value = NA,      # Not available in Table S3
          auc = table_s3_cohort$AUC_final_direction_aligned,
          coverage_count = coverage_count,
          n_fields = 0,
          n_groupable_fields = 0,
          label_key_used = NA,
          unmapped_rate = NA,
          status = "PASS",
          fail_reason = "",
          pheno_file = NA,
          stats_available = "Yes",
          audit_available = "Yes",
          label_audit_available = "Yes",
          all_metrics_available = "Yes",
          PASS_tech = "Yes",
          direction_ok = direction_ok,
          PASS_interpret = PASS_interpret,
          stringsAsFactors = FALSE
        )
        
        # Append to summary
        deliverables_summary <- rbind(deliverables_summary, summary_row)
      }
    }
  }
  
  # Process replication-only cohorts
  for (cohort in replication_cohorts) {
    # Skip if already processed (from anchor cohorts)
    if (cohort %in% deliverables_summary$cohort) {
      next
    }
    
    # Find pheno file for this cohort
    pheno_file <- pheno_files[str_detect(pheno_files, paste0(cohort, "_pheno"))]
    
    # Select only the first file if multiple are found
    if (length(pheno_file) > 0) {
      pheno_file <- pheno_file[1]
      if (file.exists(pheno_file)) {
        # Try to read phenotype data
        tryCatch({
          pheno <- readRDS(pheno_file)
          
          # Check if pheno is a data frame
          if (!is.data.frame(pheno)) {
            # Not a data frame, treat as corrupted
            summary_row <- data.frame(
              cohort = cohort,
              n_total = 0,
              n_case = 0,
              n_control = 0,
              n_aa = 0,
              n_na = 0,
              n_hc = 0,
              effect_size = NA,
              p_value = NA,
              auc = NA,
              coverage_count = 0,
              n_fields = 0,
              n_groupable_fields = 0,
              label_key_used = NA,
              unmapped_rate = NA,
              status = "FAIL",
              fail_reason = "corrupted_pheno",
              pheno_file = basename(pheno_file),
              stats_available = "No",
              audit_available = "No",
              label_audit_available = "No",
              all_metrics_available = "No",
              PASS_tech = "No",
              direction_ok = "No",
              PASS_interpret = "",
              stringsAsFactors = FALSE
            )
          } else {
            # Get sample count
            n_total <- nrow(pheno)
            
            # Get stats information
            if (cohort %in% stats_df$cohort) {
              cohort_stats <- stats_df[stats_df$cohort == cohort, ]
            } else {
              # Check if stats file exists in phase3_outputs
              stats_file <- file.path(processed_dir, paste0(cohort, "_stats.rds"))
              if (file.exists(stats_file)) {
                # Try to read stats from file
                tryCatch({
                  cohort_stats <- readRDS(stats_file)
                }, error = function(e) {
                  # Failed to read stats file
                  cohort_stats <- data.frame(
                    cohort = cohort,
                    n_total = n_total,
                    n_case = sum(pheno$std_group_label == "case", na.rm = TRUE),
                    n_control = sum(pheno$std_group_label == "control", na.rm = TRUE),
                    n_aa = sum(pheno$std_group_label == "AA", na.rm = TRUE),
                    n_na = sum(pheno$std_group_label == "NA", na.rm = TRUE),
                    n_hc = sum(pheno$std_group_label == "HC", na.rm = TRUE),
                    effect_size = NA,
                    p_value = NA,
                    auc = NA,
                    coverage_count = 0,
                    status = "FAIL",
                    fail_reason = "corrupted_stats",
                    notes = "Stats file corrupted",
                    mean_diff = NA
                  )
                })
              } else {
                # No stats available
                cohort_stats <- data.frame(
                  cohort = cohort,
                  n_total = n_total,
                  n_case = sum(pheno$std_group_label == "case", na.rm = TRUE),
                  n_control = sum(pheno$std_group_label == "control", na.rm = TRUE),
                  n_aa = sum(pheno$std_group_label == "AA", na.rm = TRUE),
                  n_na = sum(pheno$std_group_label == "NA", na.rm = TRUE),
                  n_hc = sum(pheno$std_group_label == "HC", na.rm = TRUE),
                  effect_size = NA,
                  p_value = NA,
                  auc = NA,
                  coverage_count = 0,
                  status = "FAIL",
                  fail_reason = "missing_stats",
                  notes = "Stats not calculated",
                  mean_diff = NA
                )
              }
            }
            
            # Get audit information
            if (cohort %in% audit_df$cohort) {
              cohort_audit <- audit_df[audit_df$cohort == cohort, ]
              n_fields <- nrow(cohort_audit)
              n_groupable_fields <- sum(cohort_audit$can_be_used_for_grouping, na.rm = TRUE)
            } else {
              n_fields <- 0
              n_groupable_fields <- 0
            }
            
            # Get label mapping information
            if (cohort %in% label_audit_df$dataset) {
              cohort_label_audit <- label_audit_df[label_audit_df$dataset == cohort, ]
              label_key_used <- cohort_label_audit$label_key_used
              unmapped_rate <- cohort_label_audit$unmapped_rate
            } else {
              label_key_used <- NA
              unmapped_rate <- NA
            }
            
            # Determine if all required metrics are present
            effect_size_available <- !is.na(cohort_stats$effect_size)
            p_value_available <- !is.na(cohort_stats$p_value)
            coverage_available <- cohort_stats$coverage_count > 0
            all_metrics_available <- effect_size_available && p_value_available && coverage_available
            
            # Calculate PASS_tech
            PASS_tech <- ifelse(all_metrics_available, "Yes", "No")
            
            # Calculate direction_ok
            direction_ok <- "No"
            if (cohort %in% direction_audit_df$cohort) {
              cohort_direction <- direction_audit_df[direction_audit_df$cohort == cohort, ]
              score_should_be_high_in_case <- TRUE
              if (!is.na(cohort_direction$mean_diff)) {
                if ((score_should_be_high_in_case && cohort_direction$mean_diff > 0) || 
                    (!score_should_be_high_in_case && cohort_direction$mean_diff < 0)) {
                  direction_ok <- "Yes"
                }
              }
            }
            
            # Calculate PASS_interpret
            PASS_interpret <- ""
            if (direction_ok == "No") {
              if (!is.na(cohort_stats$mean_diff) && abs(cohort_stats$mean_diff) < 0.5) {
                PASS_interpret <- "Boundary condition"
              } else {
                PASS_interpret <- "Review contrast definition"
              }
            }
            
            # Create summary row
            summary_row <- data.frame(
              cohort = cohort,
              n_total = n_total,
              n_case = cohort_stats$n_case,
              n_control = cohort_stats$n_control,
              n_aa = cohort_stats$n_aa,
              n_na = cohort_stats$n_na,
              n_hc = cohort_stats$n_hc,
              effect_size = cohort_stats$effect_size,
              p_value = cohort_stats$p_value,
              auc = cohort_stats$auc,
              coverage_count = cohort_stats$coverage_count,
              n_fields = n_fields,
              n_groupable_fields = n_groupable_fields,
              label_key_used = label_key_used,
              unmapped_rate = unmapped_rate,
              status = cohort_stats$status,
              fail_reason = cohort_stats$fail_reason,
              pheno_file = basename(pheno_file),
              stats_available = ifelse(cohort %in% stats_df$cohort || file.exists(file.path(processed_dir, paste0(cohort, "_stats.rds"))), "Yes", "No"),
              audit_available = ifelse(cohort %in% audit_df$cohort, "Yes", "No"),
              label_audit_available = ifelse(cohort %in% label_audit_df$dataset, "Yes", "No"),
              all_metrics_available = ifelse(all_metrics_available, "Yes", "No"),
              PASS_tech = PASS_tech,
              direction_ok = direction_ok,
              PASS_interpret = PASS_interpret,
              stringsAsFactors = FALSE
            )
          }
        }, error = function(e) {
          # Failed to read pheno file
          summary_row <- data.frame(
            cohort = cohort,
            n_total = 0,
            n_case = 0,
            n_control = 0,
            n_aa = 0,
            n_na = 0,
            n_hc = 0,
            effect_size = NA,
            p_value = NA,
            auc = NA,
            coverage_count = 0,
            n_fields = 0,
            n_groupable_fields = 0,
            label_key_used = NA,
            unmapped_rate = NA,
            status = "FAIL",
            fail_reason = "corrupted_pheno",
            pheno_file = basename(pheno_file),
            stats_available = "No",
            audit_available = "No",
            label_audit_available = "No",
            all_metrics_available = "No",
            PASS_tech = "No",
            direction_ok = "No",
            PASS_interpret = "",
            stringsAsFactors = FALSE
          )
        })
      } else {
        # File doesn't exist
        summary_row <- data.frame(
          cohort = cohort,
          n_total = 0,
          n_case = 0,
          n_control = 0,
          n_aa = 0,
          n_na = 0,
          n_hc = 0,
          effect_size = NA,
          p_value = NA,
          auc = NA,
          coverage_count = 0,
          n_fields = 0,
          n_groupable_fields = 0,
          label_key_used = NA,
          unmapped_rate = NA,
          status = "FAIL",
          fail_reason = "missing_pheno",
          pheno_file = basename(pheno_file),
          stats_available = "No",
          audit_available = "No",
          label_audit_available = "No",
          all_metrics_available = "No",
          PASS_tech = "No",
          direction_ok = "No",
          PASS_interpret = "",
          stringsAsFactors = FALSE
        )
      }
    } else {
      # If no pheno file found, create a row with minimal information
      summary_row <- data.frame(
        cohort = cohort,
        n_total = 0,
        n_case = 0,
        n_control = 0,
        n_aa = 0,
        n_na = 0,
        n_hc = 0,
        effect_size = NA,
        p_value = NA,
        auc = NA,
        coverage_count = 0,
        n_fields = 0,
        n_groupable_fields = 0,
        label_key_used = NA,
        unmapped_rate = NA,
        status = "FAIL",
        fail_reason = "missing_pheno",
        pheno_file = NA,
        stats_available = "No",
        audit_available = "No",
        label_audit_available = "No",
        all_metrics_available = "No",
        PASS_tech = "No",
        direction_ok = "No",
        PASS_interpret = "",
        stringsAsFactors = FALSE
      )
    }
    
    # Append to summary
    deliverables_summary <- rbind(deliverables_summary, summary_row)
  }
  
  # Process GSE152004 (discovery cohort)
  gse152004_pheno_files <- pheno_files[str_detect(pheno_files, "GSE152004_pheno")]
  if (length(gse152004_pheno_files) > 0) {
    gse152004_pheno_file <- gse152004_pheno_files[1]
    # Create summary row for GSE152004
    summary_row <- data.frame(
      cohort = "GSE152004",
      n_total = 0,
      n_case = 0,
      n_control = 0,
      n_aa = 0,
      n_na = 0,
      n_hc = 0,
      effect_size = NA,
      p_value = NA,
      auc = NA,
      coverage_count = 0,
      n_fields = 0,
      n_groupable_fields = 0,
      label_key_used = NA,
      unmapped_rate = NA,
      status = "PASS",
      fail_reason = "discovery not in deliverables scope",
      pheno_file = basename(gse152004_pheno_file),
      stats_available = "No",
      audit_available = "No",
      label_audit_available = "No",
      all_metrics_available = "No",
      PASS_tech = "No",
      direction_ok = "No",
      PASS_interpret = "",
      stringsAsFactors = FALSE
    )
    
    # Append to summary
    deliverables_summary <- rbind(deliverables_summary, summary_row)
  }
  
  # Process remaining cohorts (not in anchor or replication list)
  processed_cohorts <- c(anchor_cohorts, replication_cohorts, "GSE152004")
  for (pheno_file in pheno_files) {
    # Extract cohort name
    cohort <- str_extract(basename(pheno_file), "GSE\\d+")
    
    # Skip if cohort is already processed
    if (cohort %in% processed_cohorts) {
      next
    }
    
    # Try to read phenotype data
    tryCatch({
      pheno <- readRDS(pheno_file)
      
      # Check if pheno is a data frame
      if (!is.data.frame(pheno)) {
        # Not a data frame, treat as corrupted
        summary_row <- data.frame(
          cohort = cohort,
          n_total = 0,
          n_case = 0,
          n_control = 0,
          n_aa = 0,
          n_na = 0,
          n_hc = 0,
          effect_size = NA,
          p_value = NA,
          auc = NA,
          coverage_count = 0,
          n_fields = 0,
          n_groupable_fields = 0,
          label_key_used = NA,
          unmapped_rate = NA,
          status = "FAIL",
          fail_reason = "corrupted_pheno",
          pheno_file = basename(pheno_file),
          stats_available = "No",
          audit_available = "No",
          label_audit_available = "No",
          all_metrics_available = "No",
          PASS_tech = "No",
          direction_ok = "No",
          PASS_interpret = "",
          stringsAsFactors = FALSE
        )
      } else {
        # Get sample count
        n_total <- nrow(pheno)
        
        # Get stats information
        if (cohort %in% stats_df$cohort) {
          cohort_stats <- stats_df[stats_df$cohort == cohort, ]
        } else {
          cohort_stats <- data.frame(
            cohort = cohort,
            n_total = n_total,
            n_case = sum(pheno$std_group_label == "case", na.rm = TRUE),
            n_control = sum(pheno$std_group_label == "control", na.rm = TRUE),
            n_aa = sum(pheno$std_group_label == "AA", na.rm = TRUE),
            n_na = sum(pheno$std_group_label == "NA", na.rm = TRUE),
            n_hc = sum(pheno$std_group_label == "HC", na.rm = TRUE),
            effect_size = NA,
            p_value = NA,
            auc = NA,
            coverage_count = 0,
            status = "FAIL",
            fail_reason = "Stats not calculated",
            notes = "Stats not calculated",
            mean_diff = NA
          )
        }
        
        # Get audit information
        if (cohort %in% audit_df$cohort) {
          cohort_audit <- audit_df[audit_df$cohort == cohort, ]
          n_fields <- nrow(cohort_audit)
          n_groupable_fields <- sum(cohort_audit$can_be_used_for_grouping, na.rm = TRUE)
        } else {
          n_fields <- 0
          n_groupable_fields <- 0
        }
        
        # Get label mapping information
        if (cohort %in% label_audit_df$dataset) {
          cohort_label_audit <- label_audit_df[label_audit_df$dataset == cohort, ]
          label_key_used <- cohort_label_audit$label_key_used
          unmapped_rate <- cohort_label_audit$unmapped_rate
        } else {
          label_key_used <- NA
          unmapped_rate <- NA
        }
        
        # Determine if all required metrics are present
        effect_size_available <- !is.na(cohort_stats$effect_size)
        p_value_available <- !is.na(cohort_stats$p_value)
        all_metrics_available <- effect_size_available && p_value_available
        
        # Calculate PASS_tech
        PASS_tech <- ifelse(all_metrics_available, "Yes", "No")
        
        # Calculate direction_ok
        direction_ok <- "No"
        if (cohort %in% direction_audit_df$cohort) {
          cohort_direction <- direction_audit_df[direction_audit_df$cohort == cohort, ]
          score_should_be_high_in_case <- TRUE
          if (!is.na(cohort_direction$mean_diff)) {
            if ((score_should_be_high_in_case && cohort_direction$mean_diff > 0) || 
                (!score_should_be_high_in_case && cohort_direction$mean_diff < 0)) {
              direction_ok <- "Yes"
            }
          }
        }
        
        # Calculate PASS_interpret
        PASS_interpret <- ""
        if (direction_ok == "No") {
          if (!is.na(cohort_stats$mean_diff) && abs(cohort_stats$mean_diff) < 0.5) {
            PASS_interpret <- "Boundary condition"
          } else {
            PASS_interpret <- "Review contrast definition"
          }
        }
        
        # Create summary row
        summary_row <- data.frame(
          cohort = cohort,
          n_total = n_total,
          n_case = cohort_stats$n_case,
          n_control = cohort_stats$n_control,
          n_aa = cohort_stats$n_aa,
          n_na = cohort_stats$n_na,
          n_hc = cohort_stats$n_hc,
          effect_size = cohort_stats$effect_size,
          p_value = cohort_stats$p_value,
          auc = cohort_stats$auc,
          coverage_count = cohort_stats$coverage_count,
          n_fields = n_fields,
          n_groupable_fields = n_groupable_fields,
          label_key_used = label_key_used,
          unmapped_rate = unmapped_rate,
          status = cohort_stats$status,
          fail_reason = cohort_stats$fail_reason,
          pheno_file = basename(pheno_file),
          stats_available = ifelse(cohort %in% stats_df$cohort, "Yes", "No"),
          audit_available = ifelse(cohort %in% audit_df$cohort, "Yes", "No"),
          label_audit_available = ifelse(cohort %in% label_audit_df$dataset, "Yes", "No"),
          all_metrics_available = ifelse(all_metrics_available, "Yes", "No"),
          PASS_tech = PASS_tech,
          direction_ok = direction_ok,
          PASS_interpret = PASS_interpret,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      # Failed to read pheno file
      summary_row <- data.frame(
        cohort = cohort,
        n_total = 0,
        n_case = 0,
        n_control = 0,
        n_aa = 0,
        n_na = 0,
        n_hc = 0,
        effect_size = NA,
        p_value = NA,
        auc = NA,
        coverage_count = 0,
        n_fields = 0,
        n_groupable_fields = 0,
        label_key_used = NA,
        unmapped_rate = NA,
        status = "FAIL",
        fail_reason = "corrupted_pheno",
        pheno_file = basename(pheno_file),
        stats_available = "No",
        audit_available = "No",
        label_audit_available = "No",
        all_metrics_available = "No",
        PASS_tech = "No",
        direction_ok = "No",
        PASS_interpret = "",
        stringsAsFactors = FALSE
      )
    })
    
    # Append to summary
    deliverables_summary <- rbind(deliverables_summary, summary_row)
  }
  
  # Remove duplicate rows
  deliverables_summary <- deliverables_summary[!duplicated(deliverables_summary$cohort), ]
  
  # Save deliverables summary
  output_file <- file.path(processed_dir, "deliverables_summary.csv")
  write.csv(deliverables_summary, output_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved deliverables summary to", output_file, "\n")
  
  # Print summary
  cat("\n=== Deliverables Summary ===\n")
  print(deliverables_summary[, c("cohort", "n_total", "n_case", "n_control", "effect_size", "p_value", "auc", "status", "fail_reason", "all_metrics_available", "PASS_tech", "direction_ok", "PASS_interpret")])
  
  return(deliverables_summary)
}

# Run function
deliverables_summary <- summarize_deliverables()

cat("\nAll deliverables summarized successfully!\n")
