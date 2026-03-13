# =========================================================================
# Script: 13_tableS8_age_stratum_extract.R
# Purpose: Extract age strata from series matrix files for Table S8
# Author: Zuo 
# Date: 2026-03-01
#
# Inputs:
#   - data/raw/gse123750&GPL13158/GSE123750_series_matrix.txt.gz
#   - data/raw/gse43696&GPL6480/GSE43696_series_matrix.txt.gz
#
# Outputs:
#   - data/derived/TableS8_age_stratum_map.csv
#   - output/logs/13_tableS8_age_stratum_extract.log
# =========================================================================

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "analysis")) || !dir.exists(file.path(base_dir, "data_preparation"))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define paths using file.path()
raw_dir <- file.path(base_dir, "data", "raw")
derived_dir <- file.path(base_dir, "data", "derived")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "13_tableS8_age_stratum_extract.log")
sink(log_file, append = FALSE, split = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(derived_dir, "TableS8_age_stratum_map.csv")))
cat(sprintf("  - %s\n\n", log_file))

cat("=======================================================\n")
cat("=== Extracting Age Strata for Table S8 ===\n")
cat("=======================================================\n")

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(stringr)
})

# -------------------------------------------------------------------------# Real Analysis Logic Starts Here# -------------------------------------------------------------------------

# Function to extract age stratum from GSE123750
extract_age_stratum_gse123750 <- function(gse) {
  # Get pheno data
  pheno_data <- gse@phenoData@data
  
  # Print column names for debugging
  cat("GSE123750 pheno data columns:\n")
  print(colnames(pheno_data))
  
  # Get sample titles and geo accessions
  if ("title" %in% colnames(pheno_data)) {
    sample_titles <- pheno_data$title
  } else {
    sample_titles <- NA_character_
    cat("Warning: 'title' column not found in GSE123750 pheno data\n")
  }
  
  if ("geo_accession" %in% colnames(pheno_data)) {
    geo_accessions <- pheno_data$geo_accession
  } else if ("sample_id" %in% colnames(pheno_data)) {
    geo_accessions <- pheno_data$sample_id
  } else {
    geo_accessions <- NA_character_
    cat("Warning: No geo_accession column found in GSE123750 pheno data\n")
  }
  
  # Print first few sample titles for debugging
  cat("First 10 sample titles in GSE123750:\n")
  print(head(sample_titles, 10))
  
  # Extract age stratum from title using simple string matching
  age_stratum <- rep(NA_character_, length(sample_titles))
  
  # Use simple substring matching
  for (i in 1:length(sample_titles)) {
    title <- sample_titles[i]
    if (grepl("pre-school age", title, ignore.case = TRUE)) {
      age_stratum[i] <- "Pre-school age"
    } else if (grepl("school age", title, ignore.case = TRUE)) {
      age_stratum[i] <- "School age"
    }
  }
  
  # Print first few results for debugging
  cat("First 10 age_stratum results:\n")
  print(head(age_stratum, 10))
  
  # Determine source based on where age_stratum was extracted from
  age_stratum_source <- rep(NA_character_, length(sample_titles))
  
  for (i in 1:length(sample_titles)) {
    title <- sample_titles[i]
    if (grepl("pre-school age", title, ignore.case = TRUE) || grepl("school age", title, ignore.case = TRUE)) {
      age_stratum_source[i] <- "title"
    }
  }
  
  # Also check age group:ch1 column if available
  if ("age group:ch1" %in% colnames(pheno_data)) {
    age_group_ch1 <- pheno_data$`age group:ch1`
    # Update age_stratum where it's NA
    age_stratum <- case_when(
      !is.na(age_stratum) ~ age_stratum,
      str_detect(age_group_ch1, regex("pre.?school", ignore_case = TRUE)) ~ "Pre-school age",
      str_detect(age_group_ch1, regex("school age", ignore_case = TRUE)) ~ "School age",
      str_detect(age_group_ch1, regex("pre.?school", ignore_case = TRUE)) ~ "Pre-school age",
      str_detect(age_group_ch1, regex("school", ignore_case = TRUE)) ~ "School age",
      TRUE ~ NA_character_
    )
    # Update source
    age_stratum_source <- case_when(
      !is.na(age_stratum_source) ~ age_stratum_source,
      !is.na(age_stratum) ~ "age group:ch1",
      TRUE ~ NA_character_
    )
  }
  
  # Create dataframe
  result <- data.frame(
    `Dataset (GEO)` = "GSE123750",
    geo_accession = geo_accessions,
    age_years = NA_real_,
    `Age stratum` = age_stratum,
    age_stratum_source = age_stratum_source,
    sample_title = sample_titles,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Check for NA values
  na_count <- sum(is.na(result$`Age stratum`))
  if (na_count > 0) {
    cat(sprintf("\nERROR: Found %d NA values in GSE123750 age stratum extraction\n", na_count))
    cat("First 20 samples with NA age stratum:\n")
    na_samples <- result %>% filter(is.na(`Age stratum`)) %>% head(20)
    print(na_samples[, c("geo_accession", "sample_title")])
    # For debugging, don't stop yet
    # stop("NA values found in GSE123750 age stratum extraction")
  }
  
  return(result)
}

# Function to extract age stratum from GSE43696
extract_age_stratum_gse43696 <- function(gse) {
  # Get pheno data
  pheno_data <- gse@phenoData@data
  
  # Print column names for debugging
  cat("GSE43696 pheno data columns:\n")
  print(colnames(pheno_data))
  
  # Get geo accessions and sample titles
  if ("geo_accession" %in% colnames(pheno_data)) {
    geo_accessions <- pheno_data$geo_accession
  } else if ("sample_id" %in% colnames(pheno_data)) {
    geo_accessions <- pheno_data$sample_id
  } else {
    geo_accessions <- NA_character_
    cat("Warning: No geo_accession column found in GSE43696 pheno data\n")
  }
  
  if ("title" %in% colnames(pheno_data)) {
    sample_titles <- pheno_data$title
  } else {
    sample_titles <- NA_character_
    cat("Warning: 'title' column not found in GSE43696 pheno data\n")
  }
  
  # Extract age from age:ch1 column
  if ("age:ch1" %in% colnames(pheno_data)) {
    age_years <- pheno_data$`age:ch1`
    # Convert to numeric
    age_years <- as.numeric(age_years)
  } else {
    age_years <- NA_real_
    cat("Warning: 'age:ch1' column not found in GSE43696 pheno data\n")
  }
  
  # Print first few age values for debugging
  cat("First 10 age values in GSE43696:\n")
  print(head(age_years, 10))
  
  # Create age stratum
  age_stratum <- rep(NA_character_, length(age_years))
  
  for (i in 1:length(age_years)) {
    age <- age_years[i]
    if (!is.na(age)) {
      if (age < 40) {
        age_stratum[i] <- "Young adult"
      } else {
        age_stratum[i] <- "Older adult"
      }
    }
  }
  
  # Determine source
  age_stratum_source <- rep(NA_character_, length(age_years))
  
  for (i in 1:length(age_years)) {
    if (!is.na(age_years[i])) {
      age_stratum_source[i] <- "age:ch1"
    }
  }
  
  # Create dataframe
  result <- data.frame(
    `Dataset (GEO)` = "GSE43696",
    geo_accession = geo_accessions,
    age_years = age_years,
    `Age stratum` = age_stratum,
    age_stratum_source = age_stratum_source,
    sample_title = sample_titles,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Check for NA values
  na_count <- sum(is.na(result$`Age stratum`))
  if (na_count > 0) {
    cat(sprintf("\nERROR: Found %d NA values in GSE43696 age stratum extraction\n", na_count))
    cat("First 20 samples with NA age stratum:\n")
    na_samples <- result %>% filter(is.na(`Age stratum`)) %>% head(20)
    print(na_samples[, c("geo_accession", "age_years", "Age stratum")])
    # For debugging, don't stop yet
    # stop("NA values found in GSE43696 age stratum extraction")
  }
  
  return(result)
}

# Process GSE123750
cat("\nProcessing GSE123750...\n")
gse123750_path <- file.path(raw_dir, "gse123750&GPL13158", "GSE123750_series_matrix.txt.gz")
if (!file.exists(gse123750_path)) {
  stop("GSE123750 series matrix file not found!")
}

gse123750 <- getGEO(filename = gse123750_path, GSEMatrix = TRUE, getGPL = FALSE)
gse123750_data <- extract_age_stratum_gse123750(gse123750)

# Print验收信息 for GSE123750
cat(sprintf("GSE123750 - Total samples: %d\n", nrow(gse123750_data)))

# Add debug info
cat("First 5 rows of gse123750_data:\n")
print(head(gse123750_data, 5))

cat("Checking Age stratum column:\n")
print(gse123750_data$`Age stratum`[1:10])

cat("Checking NA values:\n")
print(is.na(gse123750_data$`Age stratum`[1:10]))

success_count_123750 <- sum(!is.na(gse123750_data$`Age stratum`))
cat(sprintf("GSE123750 - Successfully parsed age strata: %d\n", success_count_123750))
cat("GSE123750 - Age stratum counts:\n")
print(table(gse123750_data$`Age stratum`))

# Check for NA values and stop if any
if (success_count_123750 < nrow(gse123750_data)) {
  stop("NA values found in GSE123750 age stratum extraction")
}

# Process GSE43696
cat("\nProcessing GSE43696...\n")
gse43696_path <- file.path(raw_dir, "gse43696&GPL6480", "GSE43696_series_matrix.txt.gz")
if (!file.exists(gse43696_path)) {
  stop("GSE43696 series matrix file not found!")
}

gse43696 <- getGEO(filename = gse43696_path, GSEMatrix = TRUE, getGPL = FALSE)
gse43696_data <- extract_age_stratum_gse43696(gse43696)

# Print验收信息 for GSE43696
cat(sprintf("GSE43696 - Total samples: %d\n", nrow(gse43696_data)))
success_count_43696 <- sum(!is.na(gse43696_data$`Age stratum`))
cat(sprintf("GSE43696 - Successfully parsed age strata: %d\n", success_count_43696))
cat("GSE43696 - Age stratum counts:\n")
print(table(gse43696_data$`Age stratum`))

# Check for NA values and stop if any
if (success_count_43696 < nrow(gse43696_data)) {
  stop("NA values found in GSE43696 age stratum extraction")
}

# Combine results
combined_data <- bind_rows(gse123750_data, gse43696_data)

# Save to file
output_file <- file.path(derived_dir, "TableS8_age_stratum_map.csv")
write_csv(combined_data, output_file)
cat(sprintf("\nSaved TableS8_age_stratum_map.csv: %s\n", output_file))

# Memory cleanup
rm(list = ls(pattern = "gse", all.names = TRUE))
gc()

# Close logging
sink()
