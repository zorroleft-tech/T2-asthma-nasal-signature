#!/usr/bin/env Rscript

# =========================================================================
# Script: 06_generate_audit_report.R
# Purpose: Generate field audit reports for raw data sources
# Author: Zuo
# Date: 2026-02-28
#
# Inputs:
#   - data/raw/gse45111&GPL6104/GSE45111_series_matrix.txt.gz
#   - data/raw/gse45111&GPL6104/GSE45111_family.soft.gz
#   - data/raw/gse201955&GPL20301&GPL6791/GSE201955-GPL16791_series_matrix.txt.gz
#   - data/raw/gse201955&GPL20301&GPL6791/GSE201955_family.soft.gz
#   - data/raw/gse65204&GPL14550/GSE65204_series_matrix.txt.gz
#   - data/raw/gse65204&GPL14550/GSE65204_family.soft.gz
#   - data/raw/gse118761&GPL11154/GSE118761_series_matrix.txt.gz
#   - data/raw/gse118761&GPL11154/GSE118761_family.soft.gz
#   - data/raw/gse40888&GPL6244/GSE40888_series_matrix.txt.gz
#   - data/raw/gse40888&GPL6244/GSE40888_family.soft.gz
#
# Outputs:
#   - data/processed/phase3_outputs/field_audit_report.csv  # Field audit report
#   - data/processed/phase3_outputs/label_mapping_audit.csv  # Label mapping audit report
# =========================================================================

library(dplyr)
library(stringr)

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
raw_dir <- file.path(base_dir, "data", "raw")
processed_dir <- file.path(base_dir, "data", "processed", "phase3_outputs")

# Create output directory if it doesn't exist
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Define function: Extract field information from Series Matrix file
extract_fields_from_matrix <- function(matrix_file) {
  cat("Reading Series Matrix file for audit:", matrix_file, "\n")
  
  tryCatch({
    # Read file content directly
    con <- gzfile(matrix_file, "rt")
    lines <- readLines(con)
    close(con)
    
    # Extract all fields starting with !Sample_
    sample_lines <- grep("^!Sample_", lines, value = TRUE)
    
    # Extract field names
    fields <- unique(gsub("^!Sample_([^=]+) = .*", "\1", sample_lines))
    
    # Extract dataset ID
    dataset <- str_extract(matrix_file, "GSE\\d+")
    
    # Initialize field information list
    field_info <- list()
    
    # Analyze each field
    for (field in fields) {
      # Extract all values for this field
      field_lines <- grep(paste0("^!Sample_", field, " = "), lines, value = TRUE)
      field_values <- gsub(paste0("^!Sample_", field, " = "), "", field_lines)
      
      # Calculate non-empty rate
      non_empty_count <- sum(field_values != "")
      total_count <- length(field_values)
      non_empty_rate <- ifelse(total_count > 0, non_empty_count / total_count, 0)
      
      # Extract unique values
      unique_values <- unique(field_values[field_values != ""])
      
      # Determine if field can be used for grouping
      can_be_used_for_grouping <- length(unique_values) > 1 && length(unique_values) < total_count * 0.5
      
      # Save field information
      field_info[[field]] <- list(
        dataset = dataset,
        field_name = field,
        non_empty_rate = non_empty_rate,
        unique_values = unique_values,
        can_be_used_for_grouping = can_be_used_for_grouping,
        source = "matrix"
      )
    }
    
    return(field_info)
  }, error = function(e) {
    cat("Error reading matrix file for audit:", e$message, "\n")
    return(list())
  })
}

# Define function: Extract field information from SOFT file
extract_fields_from_soft <- function(soft_file) {
  cat("Reading SOFT file for audit:", soft_file, "\n")
  
  tryCatch({
    # Read file content directly
    con <- gzfile(soft_file, "rt")
    lines <- readLines(con)
    close(con)
    
    # Extract all fields starting with !Sample_
    sample_lines <- grep("^!Sample_", lines, value = TRUE)
    
    # Extract field names
    fields <- unique(gsub("^!Sample_([^=]+) = .*", "\1", sample_lines))
    
    # Extract dataset ID
    dataset <- str_extract(soft_file, "GSE\\d+")
    
    # Initialize field information list
    field_info <- list()
    
    # Analyze each field
    for (field in fields) {
      # Extract all values for this field
      field_lines <- grep(paste0("^!Sample_", field, " = "), lines, value = TRUE)
      field_values <- gsub(paste0("^!Sample_", field, " = "), "", field_lines)
      
      # Calculate non-empty rate
      non_empty_count <- sum(field_values != "")
      total_count <- length(field_values)
      non_empty_rate <- ifelse(total_count > 0, non_empty_count / total_count, 0)
      
      # Extract unique values
      unique_values <- unique(field_values[field_values != ""])
      
      # Determine if field can be used for grouping
      can_be_used_for_grouping <- length(unique_values) > 1 && length(unique_values) < total_count * 0.5
      
      # Save field information
      field_info[[field]] <- list(
        dataset = dataset,
        field_name = field,
        non_empty_rate = non_empty_rate,
        unique_values = unique_values,
        can_be_used_for_grouping = can_be_used_for_grouping,
        source = "soft"
      )
    }
    
    return(field_info)
  }, error = function(e) {
    cat("Error reading soft file for audit:", e$message, "\n")
    return(list())
  })
}

# Define function: Generate field audit report
generate_audit_report <- function() {
  # Define file paths for cohorts
  cohorts <- list(
    "GSE45111" = list(
      matrix = file.path(raw_dir, "gse45111&GPL6104", "GSE45111_series_matrix.txt.gz"),
      soft = file.path(raw_dir, "gse45111&GPL6104", "GSE45111_family.soft.gz")
    ),
    "GSE201955" = list(
      matrix = file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955-GPL16791_series_matrix.txt.gz"),
      soft = file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_family.soft.gz")
    ),
    "GSE65204" = list(
      matrix = file.path(raw_dir, "gse65204&GPL14550", "GSE65204_series_matrix.txt.gz"),
      soft = file.path(raw_dir, "gse65204&GPL14550", "GSE65204_family.soft.gz")
    ),
    "GSE118761" = list(
      matrix = file.path(raw_dir, "gse118761&GPL11154", "GSE118761_series_matrix.txt.gz"),
      soft = file.path(raw_dir, "gse118761&GPL11154", "GSE118761_family.soft.gz")
    ),
    "GSE40888" = list(
      matrix = file.path(raw_dir, "gse40888&GPL6244", "GSE40888_series_matrix.txt.gz"),
      soft = file.path(raw_dir, "gse40888&GPL6244", "GSE40888_family.soft.gz")
    )
  )
  
  # Initialize audit report list
  audit_report <- list()
  
  # Analyze each cohort
  for (cohort in names(cohorts)) {
    cat("\nGenerating audit report for", cohort, "...\n")
    
    # Extract field information
    matrix_fields <- extract_fields_from_matrix(cohorts[[cohort]]$matrix)
    soft_fields <- extract_fields_from_soft(cohorts[[cohort]]$soft)
    
    # Merge field information
    all_fields <- c(matrix_fields, soft_fields)
    
    # Deduplicate (keep soft fields if duplicate with matrix fields)
    unique_fields <- list()
    for (field_name in names(all_fields)) {
      if (!(field_name %in% names(unique_fields)) || all_fields[[field_name]]$source == "soft") {
        unique_fields[[field_name]] <- all_fields[[field_name]]
      }
    }
    
    # Save to audit report
    audit_report[[cohort]] <- unique_fields
  }
  
  # Generate CSV format audit report
  audit_df <- data.frame()
  
  for (cohort in names(audit_report)) {
    fields <- audit_report[[cohort]]
    for (field_name in names(fields)) {
      field_info <- fields[[field_name]]
      row <- data.frame(
        cohort = cohort,
        field_name = field_info$field_name,
        source = field_info$source,
        non_empty_rate = field_info$non_empty_rate,
        unique_values_count = length(field_info$unique_values),
        unique_values = paste(field_info$unique_values, collapse = "; "),
        can_be_used_for_grouping = field_info$can_be_used_for_grouping
      )
      audit_df <- rbind(audit_df, row)
    }
  }
  
  # Save audit report
  output_file <- file.path(processed_dir, "field_audit_report.csv")
  write.csv(audit_df, output_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat("\nSaved field audit report to", output_file, "\n")
  
  # Print summary
  cat("\n=== Field Audit Summary ===\n")
  print(table(audit_df$cohort, audit_df$can_be_used_for_grouping))
  
  return(audit_report)
}

# Define function: Generate label mapping audit report
generate_label_mapping_audit <- function() {
  # Read all pheno_raw files
  pheno_files <- list.files(processed_dir, pattern = "_pheno_raw\\.rds$", full.names = TRUE)
  
  # Initialize label mapping audit list
  label_audit_list <- list()
  
  # Analyze each cohort
  for (pheno_file in pheno_files) {
    # Extract cohort name
    cohort <- str_extract(basename(pheno_file), "GSE\\d+")
    cat("\nGenerating label mapping audit for", cohort, "...\n")
    
    # Read phenotype data
    pheno <- readRDS(pheno_file)
    
    # Extract label key used
    label_key_used <- NA
    if ("std_label_reference" %in% colnames(pheno)) {
      # Extract the key from the reference
      references <- pheno$std_label_reference[!is.na(pheno$std_label_reference)]
      if (length(references) > 0) {
        # Get the most common reference
        reference_counts <- table(references)
        most_common <- names(reference_counts)[which.max(reference_counts)]
        # Extract the key from the reference string
        label_key_used <- str_extract(most_common, "[^:]+$")
      }
    } else if ("std_group_label_source" %in% colnames(pheno)) {
      # Fallback to group label source
      sources <- pheno$std_group_label_source[!is.na(pheno$std_group_label_source)]
      if (length(sources) > 0) {
        # Get the most common source
        source_counts <- table(sources)
        most_common <- names(source_counts)[which.max(source_counts)]
        # Extract the key from the source string
        label_key_used <- str_extract(most_common, "[^:]+$")
      }
    }
    
    # Extract raw unique values top
    raw_unique_values_top <- ""
    if ("raw_characteristics_all" %in% colnames(pheno)) {
      # Extract values for the label key
      if (!is.na(label_key_used)) {
        # Extract values from raw_characteristics_all
        values <- sapply(pheno$raw_characteristics_all, function(x) {
          match_pattern <- paste0(label_key_used, ": ([^|]+)")
          match <- str_match(x, match_pattern)
          if (!is.na(match[1, 2])) {
            match[1, 2]
          } else {
            NA
          }
        })
        
        # Calculate value counts
        value_counts <- table(values[!is.na(values)])
        # Get top 5 values
        top_values <- head(sort(value_counts, decreasing = TRUE), 5)
        if (length(top_values) > 0) {
          raw_unique_values_top <- paste(names(top_values), top_values, sep = ": ", collapse = "; ")
        }
      }
    }
    
    # Calculate std_group_label counts
    std_group_label_counts <- ""
    if ("std_group_label" %in% colnames(pheno)) {
      label_counts <- table(pheno$std_group_label, useNA = "always")
      std_group_label_counts <- paste(names(label_counts), label_counts, sep = ": ", collapse = "; ")
    }
    
    # Calculate unmapped rate
    unmapped_rate <- 0
    if ("std_group_label" %in% colnames(pheno)) {
      n_total <- nrow(pheno)
      n_unmapped <- sum(is.na(pheno$std_group_label))
      unmapped_rate <- (n_unmapped / n_total) * 100
    }
    
    # Extract examples of unmapped samples
    examples_unmapped <- ""
    if ("std_group_label" %in% colnames(pheno) && "raw_characteristics_all" %in% colnames(pheno)) {
      unmapped_samples <- pheno[is.na(pheno$std_group_label), ]
      if (nrow(unmapped_samples) > 0) {
        # Take up to 3 examples
        n_examples <- min(3, nrow(unmapped_samples))
        sampled <- unmapped_samples[sample(nrow(unmapped_samples), n_examples), ]
        examples <- paste(sampled$std_sample_id, sampled$raw_characteristics_all, sep = ": ", collapse = " || ")
        examples_unmapped <- examples
      }
    }
    
    # Create audit entry with all possible columns
    audit_entry <- data.frame(
      dataset = cohort,
      label_key_used = label_key_used,
      raw_unique_values_top = raw_unique_values_top,
      std_group_label_counts = std_group_label_counts,
      unmapped_rate = unmapped_rate,
      examples_unmapped = examples_unmapped,
      unique_subject_count = NA,
      unique_site_count = NA,
      complete_pairs = NA,
      total_subjects = NA,
      pair_completeness = NA,
      aa_count = NA,
      na_count = NA,
      hc_count = NA,
      stringsAsFactors = FALSE
    )
    
    # Add direction cohort specific audit fields
    if (cohort == "GSE118761") {
      # Add subject_id/site unique counts and pair completeness
      if ("std_subject_id" %in% colnames(pheno) && "std_site" %in% colnames(pheno)) {
        # Calculate unique subject count
        unique_subject_count <- length(unique(pheno$std_subject_id))
        # Calculate unique site count
        unique_site_count <- length(unique(pheno$std_site))
        # Check pair completeness
        subject_site_counts <- table(pheno$std_subject_id, pheno$std_site)
        complete_pairs <- sum(rowSums(subject_site_counts > 0) == 2)
        total_subjects <- nrow(subject_site_counts)
        pair_completeness <- (complete_pairs / total_subjects) * 100
        
        # Add to audit entry
        audit_entry$unique_subject_count <- unique_subject_count
        audit_entry$unique_site_count <- unique_site_count
        audit_entry$complete_pairs <- complete_pairs
        audit_entry$total_subjects <- total_subjects
        audit_entry$pair_completeness <- pair_completeness
      }
    } else if (cohort == "GSE40888") {
      # Add AA/NA/HC specific fields
      if ("std_group_label" %in% colnames(pheno)) {
        # Calculate AA/NA/HC counts
        aa_count <- sum(pheno$std_group_label == "AA", na.rm = TRUE)
        na_count <- sum(pheno$std_group_label == "NA", na.rm = TRUE)
        hc_count <- sum(pheno$std_group_label == "HC", na.rm = TRUE)
        
        # Add to audit entry
        audit_entry$aa_count <- aa_count
        audit_entry$na_count <- na_count
        audit_entry$hc_count <- hc_count
      }
    }
    
    label_audit_list[[cohort]] <- audit_entry
  }
  
  # Combine into a single data frame
  label_audit_df <- do.call(rbind, label_audit_list)
  
  # Save label mapping audit report
  output_file <- file.path(processed_dir, "label_mapping_audit.csv")
  write.csv(label_audit_df, output_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat("\nSaved label mapping audit report to", output_file, "\n")
  
  # Print summary
  cat("\n=== Label Mapping Audit Summary ===\n")
  print(label_audit_df[, c("dataset", "label_key_used", "unmapped_rate")])
  
  return(label_audit_df)
}

# Run functions
generate_audit_report()
generate_label_mapping_audit()

cat("\nAll field audit reports generated successfully!\n")