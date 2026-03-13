#!/usr/bin/env Rscript

# =========================================================================
# Script: 04_extract_pheno_from_raw.R
# Purpose: Extract phenotype information from raw data sources using GEOquery
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
#
# Outputs:
#   - data/processed/phase3_outputs/*_pheno_samples.rds  # Sample-level phenotype table
#   - data/processed/phase3_outputs/*_pheno_kv_long.rds  # Key-value long table
#   - data/processed/phase3_outputs/*_pheno_raw.rds       # Final standardized phenotype table
# =========================================================================

# Load necessary packages
required_packages <- c("GEOquery", "dplyr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(GEOquery)
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

cat("Base directory:", base_dir, "\n")
cat("Raw data directory:", raw_dir, "\n")
cat("Processed directory:", processed_dir, "\n")

# =========================================================================
# Part 1: Core functions for extracting and parsing
# =========================================================================

# Normalize key names: lowercase, spaces→underscores, remove parentheses
normalize_key <- function(key) {
  key %>%
    tolower() %>%
    str_replace_all("[[:space:]]+", "_") %>%
    str_replace_all("[()]", "") %>%
    str_trim()
}

# Parse a single "key: value" string
parse_kv_string <- function(kv_string) {
  if (str_detect(kv_string, ": ")) {
    parts <- str_split(kv_string, ": ", n = 2)[[1]]
    key <- normalize_key(parts[1])
    value <- str_trim(parts[2])
  } else {
    key <- "free_text"
    value <- str_trim(kv_string)
  }
  return(list(key = key, value = value))
}

# Extract phenotype directly from Series Matrix file (handling multiple characteristics_ch1 lines)
extract_from_matrix_direct <- function(matrix_file, dataset) {
  cat("Reading Series Matrix directly:", matrix_file, "\n")
  
  tryCatch({
    # Read all lines
    con <- gzfile(matrix_file, "rt")
    lines <- readLines(con)
    close(con)
    
    # Step 1: Find sample IDs
    sample_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
    if (length(sample_line) == 0) {
      stop("No !Sample_geo_accession line found")
    }
    
    # Parse sample IDs
    sample_parts <- str_split(sample_line, "\t")[[1]]
    sample_ids <- sample_parts[-1]
    sample_ids <- gsub("\"", "", sample_ids)
    sample_ids <- sample_ids[sample_ids != ""]
    
    n_samples <- length(sample_ids)
    cat("Found", n_samples, "samples\n")
    
    # Step 2: Find all characteristics_ch1 lines
    char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
    
    # Initialize data structures
    pheno_samples <- data.frame(
      std_sample_id = sample_ids,
      std_dataset = rep(dataset, n_samples),
      raw_characteristics_all = character(n_samples),
      stringsAsFactors = FALSE
    )
    
    pheno_kv_long <- data.frame(
      std_sample_id = character(),
      key = character(),
      value = character(),
      source = character(),
      stringsAsFactors = FALSE
    )
    
    # For each sample, collect all characteristics
    sample_chars <- vector("list", n_samples)
    names(sample_chars) <- sample_ids
    
    # Process each characteristics_ch1 line
    for (char_line in char_lines) {
      parts <- str_split(char_line, "\t")[[1]]
      if (length(parts) < 2) next
      
      values <- parts[-1]
      values <- gsub("\"", "", values)
      
      if (length(values) != n_samples) {
        values <- rep(values, length.out = n_samples)
      }
      
      for (i in seq_len(n_samples)) {
        val <- values[i]
        if (val != "" && !is.na(val)) {
          kv <- parse_kv_string(val)
          if (kv$key != "") {
            if (is.null(sample_chars[[i]])) {
              sample_chars[[i]] <- list()
            }
            sample_chars[[i]][[kv$key]] <- kv$value
            
            pheno_kv_long <- rbind(pheno_kv_long, data.frame(
              std_sample_id = sample_ids[i],
              key = kv$key,
              value = kv$value,
              source = "matrix",
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    # Build raw_characteristics_all for each sample
    for (i in seq_len(n_samples)) {
      chars <- sample_chars[[i]]
      if (!is.null(chars) && length(chars) > 0) {
        char_strings <- paste(names(chars), ": ", unlist(chars), sep = "")
        pheno_samples$raw_characteristics_all[i] <- paste(char_strings, collapse = " | ")
      } else {
        pheno_samples$raw_characteristics_all[i] <- "no_characteristics_available"
      }
    }
    
    cat("Successfully extracted data for", n_samples, "samples\n")
    cat("Number of KV pairs:", nrow(pheno_kv_long), "\n")
    
    return(list(
      pheno_samples = pheno_samples,
      pheno_kv_long = pheno_kv_long
    ))
    
  }, error = function(e) {
    cat("Error reading matrix directly:", e$message, "\n")
    return(NULL)
  })
}

# =========================================================================
# Part 2: Cohort-specific mapping functions
# =========================================================================

# GSE65204 (IgE anchor) mapping
map_group_label_GSE65204 <- function(pheno_samples, pheno_kv_long) {
  dataset <- "GSE65204"
  cat("\n=== Mapping group labels for", dataset, "===\n")
  
  pheno_samples$std_group_label <- NA
  pheno_samples$std_group_label_source <- NA
  pheno_samples$std_extract_notes <- NA
  pheno_samples$std_label_definition <- NA
  pheno_samples$std_label_reference <- NA
  
  kv_subset <- pheno_kv_long[pheno_kv_long$std_sample_id %in% pheno_samples$std_sample_id, ]
  
  ige_keys <- grep("ige|atopy|allergy|sensit", kv_subset$key, value = TRUE, ignore.case = TRUE)
  
  if (length(ige_keys) > 0) {
    cat("Found candidate IgE keys:", unique(ige_keys), "\n")
    
    target_key <- "lnige"
    if (!(target_key %in% kv_subset$key)) {
      target_key <- ige_keys[1]
    }
    ige_values <- kv_subset[kv_subset$key == target_key, ]
    
    ige_values$numeric_value <- as.numeric(str_extract(ige_values$value, "[0-9.]+"))
    
    valid_numeric <- ige_values[!is.na(ige_values$numeric_value), ]
    if (nrow(valid_numeric) >= 2) {
      median_val <- median(valid_numeric$numeric_value, na.rm = TRUE)
      cat("Using median split at:", median_val, "\n")
      
      for (i in seq_len(nrow(pheno_samples))) {
        gsm_id <- pheno_samples$std_sample_id[i]
        sample_value <- ige_values$numeric_value[ige_values$std_sample_id == gsm_id]
        
        if (length(sample_value) > 0 && !is.na(sample_value)) {
          if (sample_value > median_val) {
            pheno_samples$std_group_label[i] <- "case"
          } else {
            pheno_samples$std_group_label[i] <- "control"
          }
          pheno_samples$std_group_label_source[i] <- paste0("matrix:characteristics_ch1:", target_key)
          pheno_samples$std_extract_notes[i] <- paste0("IgE median split for minimal closure (median=", round(median_val, 3), ")")
          pheno_samples$std_label_definition[i] <- "lnIgE median split (推进用)"
          pheno_samples$std_label_reference[i] <- paste0("matrix:characteristics_ch1:", target_key)
        }
      }
    }
  }
  
  cat("Group label distribution:\n")
  print(table(pheno_samples$std_group_label, useNA = "always"))
  
  return(pheno_samples)
}

# GSE201955 (FeNO anchor) mapping
map_group_label_GSE201955 <- function(pheno_samples, pheno_kv_long) {
  dataset <- "GSE201955"
  cat("\n=== Mapping group labels for", dataset, "===\n")
  
  pheno_samples$std_group_label <- NA
  pheno_samples$std_group_label_source <- NA
  pheno_samples$std_extract_notes <- NA
  pheno_samples$std_label_definition <- NA
  pheno_samples$std_label_reference <- NA
  
  kv_subset <- pheno_kv_long[pheno_kv_long$std_sample_id %in% pheno_samples$std_sample_id, ]
  
  target_key <- "asthma"
  if (target_key %in% kv_subset$key) {
    cat("Using", target_key, "field for mapping\n")
    asthma_values <- kv_subset[kv_subset$key == target_key, ]
    
    for (i in seq_len(nrow(pheno_samples))) {
      gsm_id <- pheno_samples$std_sample_id[i]
      sample_value <- asthma_values$value[asthma_values$std_sample_id == gsm_id]
      
      if (length(sample_value) > 0) {
        if (str_detect(tolower(sample_value), "asthma")) {
          pheno_samples$std_group_label[i] <- "case"
          pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: Asthma"
          pheno_samples$std_label_definition[i] <- "Direct categorical mapping: Asthma vs Control"
          pheno_samples$std_label_reference[i] <- paste0("matrix:characteristics_ch1:", target_key)
        } else if (str_detect(tolower(sample_value), "control")) {
          pheno_samples$std_group_label[i] <- "control"
          pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: Control"
          pheno_samples$std_label_definition[i] <- "Direct categorical mapping: Asthma vs Control"
          pheno_samples$std_label_reference[i] <- paste0("matrix:characteristics_ch1:", target_key)
        }
      }
    }
  } else {
    feno_keys <- grep("feno|no|fractional_exhaled", kv_subset$key, value = TRUE, ignore.case = TRUE)
    
    if (length(feno_keys) > 0) {
      cat("Found candidate FeNO keys:", unique(feno_keys), "\n")
      target_key <- feno_keys[1]
      feno_values <- kv_subset[kv_subset$key == target_key, ]
      
      feno_values$numeric_value <- as.numeric(str_extract(feno_values$value, "[0-9.]+"))
      
      valid_numeric <- feno_values[!is.na(feno_values$numeric_value), ]
      if (nrow(valid_numeric) >= 2) {
        median_val <- median(valid_numeric$numeric_value, na.rm = TRUE)
        cat("Using FeNO median split at:", median_val, "\n")
        
        for (i in seq_len(nrow(pheno_samples))) {
          gsm_id <- pheno_samples$std_sample_id[i]
          sample_value <- feno_values$numeric_value[feno_values$std_sample_id == gsm_id]
          
          if (length(sample_value) > 0 && !is.na(sample_value)) {
            if (sample_value > median_val) {
              pheno_samples$std_group_label[i] <- "case"
            } else {
              pheno_samples$std_group_label[i] <- "control"
            }
            pheno_samples$std_group_label_source[i] <- paste0("matrix:characteristics_ch1:", target_key)
            pheno_samples$std_extract_notes[i] <- paste0("FeNO median split for minimal closure (median=", round(median_val, 3), ")")
            pheno_samples$std_label_definition[i] <- "FeNO median split (推进用)"
            pheno_samples$std_label_reference[i] <- paste0("matrix:characteristics_ch1:", target_key)
          }
        }
      }
    }
  }
  
  cat("Group label distribution:\n")
  print(table(pheno_samples$std_group_label, useNA = "always"))
  
  return(pheno_samples)
}

# GSE45111 (sputum phenotype anchor) mapping
map_group_label_GSE45111 <- function(pheno_samples, pheno_kv_long) {
  dataset <- "GSE45111"
  cat("\n=== Mapping group labels for", dataset, "===\n")
  
  pheno_samples$std_group_label <- NA
  pheno_samples$std_group_label_source <- NA
  pheno_samples$std_extract_notes <- NA
  pheno_samples$std_label_definition <- NA
  pheno_samples$std_label_reference <- NA
  
  kv_subset <- pheno_kv_long[pheno_kv_long$std_sample_id %in% pheno_samples$std_sample_id, ]
  
  sputum_keys <- grep("sputum|eosin|neutro|phenotype|inflamm", kv_subset$key, value = TRUE, ignore.case = TRUE)
  
  if (length(sputum_keys) > 0) {
    cat("Found candidate sputum keys:", unique(sputum_keys), "\n")
    
    target_key <- NULL
    for (key in unique(sputum_keys)) {
      key_values <- kv_subset$value[kv_subset$key == key]
      if (any(str_detect(tolower(key_values), "eosinophilic|neutrophilic"))) {
        target_key <- key
        break
      }
    }
    
    if (!is.null(target_key)) {
      cat("Using categorical key:", target_key, "\n")
      eos_values <- kv_subset[kv_subset$key == target_key, ]
      
      for (i in seq_len(nrow(pheno_samples))) {
        gsm_id <- pheno_samples$std_sample_id[i]
        sample_value <- eos_values$value[eos_values$std_sample_id == gsm_id]
        
        if (length(sample_value) > 0) {
          if (str_detect(tolower(sample_value), "eosinophilic")) {
            pheno_samples$std_group_label[i] <- "case"
            pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: eosinophilic"
            pheno_samples$std_label_definition[i] <- "asthma_phenotype collapsed to EA vs non-EA"
            pheno_samples$std_label_reference[i] <- paste0("matrix:characteristics_ch1:", target_key)
          } else if (str_detect(tolower(sample_value), "neutrophilic|non-eosinophilic|paucigranulocytic")) {
            pheno_samples$std_group_label[i] <- "control"
            pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: neutrophilic/non-eosinophilic"
            pheno_samples$std_label_definition[i] <- "asthma_phenotype collapsed to EA vs non-EA"
            pheno_samples$std_label_reference[i] <- paste0("matrix:characteristics_ch1:", target_key)
          }
        }
      }
    }
  }
  
  cat("Group label distribution:\n")
  print(table(pheno_samples$std_group_label, useNA = "always"))
  
  return(pheno_samples)
}

# GSE118761 (paired nasal-tracheal) mapping
map_group_label_GSE118761 <- function(pheno_samples, pheno_kv_long) {
  dataset <- "GSE118761"
  cat("\n=== Mapping group labels for", dataset, "===\n")
  
  pheno_samples$std_group_label <- NA
  pheno_samples$std_group_label_source <- NA
  pheno_samples$std_extract_notes <- NA
  pheno_samples$std_label_definition <- NA
  pheno_samples$std_label_reference <- NA
  pheno_samples$std_subject_id <- NA
  pheno_samples$std_site <- NA
  pheno_samples$std_pair_id <- NA
  
  kv_subset <- pheno_kv_long[pheno_kv_long$std_sample_id %in% pheno_samples$std_sample_id, ]
  
  # Extract subject_id
  if ("subject_id" %in% kv_subset$key) {
    subject_values <- kv_subset[kv_subset$key == "subject_id", ]
    for (i in seq_len(nrow(pheno_samples))) {
      gsm_id <- pheno_samples$std_sample_id[i]
      sample_value <- subject_values$value[subject_values$std_sample_id == gsm_id]
      if (length(sample_value) > 0) {
        pheno_samples$std_subject_id[i] <- sample_value
        pheno_samples$std_pair_id[i] <- sample_value
      }
    }
  }
  
  # Extract tissue/site
  if ("tissue" %in% kv_subset$key) {
    tissue_values <- kv_subset[kv_subset$key == "tissue", ]
    for (i in seq_len(nrow(pheno_samples))) {
      gsm_id <- pheno_samples$std_sample_id[i]
      sample_value <- tissue_values$value[tissue_values$std_sample_id == gsm_id]
      if (length(sample_value) > 0) {
        if (str_detect(tolower(sample_value), "nasal")) {
          pheno_samples$std_site[i] <- "nasal"
          pheno_samples$std_group_label[i] <- "nasal"
        } else if (str_detect(tolower(sample_value), "bronchial|tracheal")) {
          pheno_samples$std_site[i] <- "tracheal"
          pheno_samples$std_group_label[i] <- "tracheal"
        }
      }
    }
  }
  
  cat("Subject ID distribution:\n")
  print(length(unique(pheno_samples$std_subject_id)))
  cat("Site distribution:\n")
  print(table(pheno_samples$std_site, useNA = "always"))
  cat("Group label distribution:\n")
  print(table(pheno_samples$std_group_label, useNA = "always"))
  
  return(pheno_samples)
}

# GSE40888 (AA/NA/HC) mapping
map_group_label_GSE40888 <- function(pheno_samples, pheno_kv_long) {
  dataset <- "GSE40888"
  cat("\n=== Mapping group labels for", dataset, "===\n")
  
  pheno_samples$std_group_label <- NA
  pheno_samples$std_group_label_source <- NA
  pheno_samples$std_extract_notes <- NA
  pheno_samples$std_label_definition <- NA
  pheno_samples$std_label_reference <- NA
  
  kv_subset <- pheno_kv_long[pheno_kv_long$std_sample_id %in% pheno_samples$std_sample_id, ]
  
  # Extract group
  if ("group" %in% kv_subset$key) {
    group_values <- kv_subset[kv_subset$key == "group", ]
    for (i in seq_len(nrow(pheno_samples))) {
      gsm_id <- pheno_samples$std_sample_id[i]
      sample_value <- group_values$value[group_values$std_sample_id == gsm_id]
      if (length(sample_value) > 0) {
        sample_value_lower <- tolower(sample_value)
        if (str_detect(sample_value_lower, "naa|non-allergic asthmatic")) {
          pheno_samples$std_group_label[i] <- "NA"
          pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: NA"
          pheno_samples$std_label_definition[i] <- "AA/NA/HC classification"
          pheno_samples$std_label_reference[i] <- "matrix:characteristics_ch1:group"
        } else if (str_detect(sample_value_lower, "aa|allergic asthmatic")) {
          pheno_samples$std_group_label[i] <- "AA"
          pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: AA"
          pheno_samples$std_label_definition[i] <- "AA/NA/HC classification"
          pheno_samples$std_label_reference[i] <- "matrix:characteristics_ch1:group"
        } else if (str_detect(sample_value_lower, "hc|healthy control|control")) {
          pheno_samples$std_group_label[i] <- "HC"
          pheno_samples$std_extract_notes[i] <- "Direct categorical mapping: HC"
          pheno_samples$std_label_definition[i] <- "AA/NA/HC classification"
          pheno_samples$std_label_reference[i] <- "matrix:characteristics_ch1:group"
        }
      }
    }
  }
  
  cat("Group label distribution:\n")
  print(table(pheno_samples$std_group_label, useNA = "always"))
  
  return(pheno_samples)
}

# Main processing function for a cohort
process_cohort <- function(dataset, matrix_file) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing cohort:", dataset, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  extraction_result <- extract_from_matrix_direct(matrix_file, dataset)
  
  if (is.null(extraction_result)) {
    cat("Failed to extract data for", dataset, "\n")
    return(NULL)
  }
  
  pheno_samples <- extraction_result$pheno_samples
  pheno_kv_long <- extraction_result$pheno_kv_long
  
  cat("\n=== Random sample characteristics validation (3 samples) ===\n")
  set.seed(42)
  random_indices <- sample(nrow(pheno_samples), min(3, nrow(pheno_samples)))
  for (i in random_indices) {
    cat("\nSample", pheno_samples$std_sample_id[i], ":\n")
    cat(pheno_samples$raw_characteristics_all[i], "\n")
  }
  
  cat("\n=== Top 20 keys from KV table ===\n")
  key_counts <- table(pheno_kv_long$key)
  print(head(sort(key_counts, decreasing = TRUE), 20))
  
  mapping_function <- paste0("map_group_label_", dataset)
  if (exists(mapping_function, mode = "function")) {
    pheno_samples <- do.call(mapping_function, list(pheno_samples, pheno_kv_long))
  } else {
    cat("Warning: No mapping function found for", dataset, "\n")
    pheno_samples$std_group_label <- NA
    pheno_samples$std_group_label_source <- NA
    pheno_samples$std_extract_notes <- "No mapping function available"
  }
  
  output_samples <- file.path(processed_dir, paste0(dataset, "_pheno_samples.rds"))
  output_kv <- file.path(processed_dir, paste0(dataset, "_pheno_kv_long.rds"))
  output_final <- file.path(processed_dir, paste0(dataset, "_pheno_raw.rds"))
  
  saveRDS(pheno_samples, output_samples)
  saveRDS(pheno_kv_long, output_kv)
  saveRDS(pheno_samples, output_final)
  
  cat("\nSaved outputs:\n")
  cat("  -", output_samples, "\n")
  cat("  -", output_kv, "\n")
  cat("  -", output_final, "\n")
  
  # Calculate mapping statistics
  n_total <- nrow(pheno_samples)
  n_mapped <- sum(!is.na(pheno_samples$std_group_label))
  n_unmapped <- n_total - n_mapped
  unmapped_rate <- (n_unmapped / n_total) * 100
  
  # Add mapping statistics to pheno_samples
  pheno_samples$n_total <- n_total
  pheno_samples$n_mapped <- n_mapped
  pheno_samples$n_unmapped <- n_unmapped
  pheno_samples$unmapped_rate <- unmapped_rate
  
  # Print mapping statistics
  cat("\n=== Mapping Statistics ===\n")
  cat("Total samples:", n_total, "\n")
  cat("Mapped samples:", n_mapped, "\n")
  cat("Unmapped samples:", n_unmapped, "\n")
  cat("Unmapped rate:", round(unmapped_rate, 2), "%\n")
  
  # Check if this is a truth cohort and validate unmapped rate
  truth_cohorts <- c("GSE65204", "GSE201955", "GSE45111")
  if (dataset %in% truth_cohorts) {
    if (unmapped_rate >= 5) {
      cat("FAIL: Truth cohort has unmapped rate >= 5%\n")
    } else {
      cat("PASS: Truth cohort has unmapped rate < 5%\n")
    }
  }
  
  cat("\n=== std_group_label_source (head 3) ===\n")
  print(head(pheno_samples$std_group_label_source, 3))
  
  return(list(
    pheno_samples = pheno_samples,
    pheno_kv_long = pheno_kv_long
  ))
}

# =========================================================================
# Main execution
# =========================================================================

cohorts <- list(
  "GSE65204" = file.path(raw_dir, "gse65204&GPL14550", "GSE65204_series_matrix.txt.gz"),
  "GSE201955" = file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955-GPL16791_series_matrix.txt.gz"),
  "GSE45111" = file.path(raw_dir, "gse45111&GPL6104", "GSE45111_series_matrix.txt.gz"),
  "GSE118761" = file.path(raw_dir, "gse118761&GPL11154", "GSE118761_series_matrix.txt.gz"),
  "GSE40888" = file.path(raw_dir, "gse40888&GPL6244", "GSE40888_series_matrix.txt.gz")
)

results <- list()
for (dataset in names(cohorts)) {
  matrix_file <- cohorts[[dataset]]
  
  if (file.exists(matrix_file)) {
    result <- process_cohort(dataset, matrix_file)
    if (!is.null(result)) {
      results[[dataset]] <- result
    }
  } else {
    cat("Warning: Matrix file not found for", dataset, ":", matrix_file, "\n")
  }
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("All cohorts processed!\n")
cat(rep("=", 60), "\n", sep = "")
