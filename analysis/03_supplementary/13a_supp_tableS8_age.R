# =========================================================================
# Script: 13_supp_tableS8_age.R
# Purpose: Generate Supplementary Table S8 with age subgroup analysis
# Author: Zuo 
# Date: 2026-03-01
#
# Inputs:
#   - data/derived/locked_weights.csv      # 11-gene locked weights
#   - data/processed_diet/*_diet_expr.rds  # Cohort expression matrices
#   - data/processed_diet/*_pheno.rds      # Cohort phenotype metadata
#   - R/03_index_logger.R                  # Logging script
#
# Outputs:
#   - output/supplement/TableS8_age_subgroups.csv  # Age subgroup analysis table
#   - output/logs/13_supp_tableS8_age.log  # Log file
# =========================================================================

# Check for Notion link pollution - simplified approach
# Directly check for specific patterns in the script environment
if (exists("pheno") && is.data.frame(get("pheno"))) {
  pheno_names <- names(get("pheno"))
  if (any(grepl(" `http://is.na` ", pheno_names))) {
    stop("Detected Notion link pollution: `http://is.na` ")
  }
}

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "analysis")) || !dir.exists(file.path(base_dir, "data_preparation"))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define paths using file.path()
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
derived_dir <- file.path(base_dir, "data", "derived")
supplement_dir <- file.path(base_dir, "output", "supplement")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "13_supp_tableS8_age.log")
sink(log_file, append = FALSE, split = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS8_age_subgroups.csv")))
cat(sprintf("  - %s\n\n", log_file))

cat("=======================================================\n")
cat("=== Generating Supplementary Table S8 (Age Subgroup) ===\n")
cat("=======================================================\n")

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(digest)
  library(GEOquery)
  library(stringr)
})

# Source utility functions
source(file.path(base_dir, "R", "03_index_logger.R"))

# -------------------------------------------------------------------------
# Real Analysis Logic Starts Here
# -------------------------------------------------------------------------
cat("\n[1/3] Loading Locked Weights and Age Stratum Map...\n")
weights_file <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(weights_file)) stop("locked_weights.csv not found!")
locked_weights <- read_csv(weights_file, show_col_types = FALSE)

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
  }
  
  return(result)
}

# Load or generate age stratum map
age_stratum_map_file <- file.path(derived_dir, "TableS8_age_stratum_map.csv")
age_stratum_map <- NULL
if (file.exists(age_stratum_map_file)) {
  cat("Loading age stratum map from", age_stratum_map_file, "...\n")
  age_stratum_map <- read_csv(age_stratum_map_file, show_col_types = FALSE)
  cat(sprintf("Loaded age stratum map with %d samples\n", nrow(age_stratum_map)))
} else {
  cat("Age stratum map not found. Generating from GEO series matrix files...\n")
  
  # Process GSE123750
  cat("\nProcessing GSE123750...\n")
  gse123750_path <- file.path(base_dir, "data", "raw", "gse123750&GPL13158", "GSE123750_series_matrix.txt.gz")
  if (file.exists(gse123750_path)) {
    gse123750 <- getGEO(filename = gse123750_path, GSEMatrix = TRUE, getGPL = FALSE)
    gse123750_data <- extract_age_stratum_gse123750(gse123750)
    
    # Check for NA values
    na_count_123750 <- sum(is.na(gse123750_data$`Age stratum`))
    if (na_count_123750 > 0) {
      cat(sprintf("\nERROR: Found %d NA values in GSE123750 age stratum extraction\n", na_count_123750))
      cat("First 20 samples with NA age stratum:\n")
      na_samples <- gse123750_data %>% filter(is.na(`Age stratum`)) %>% head(20)
      print(na_samples[, c("geo_accession", "sample_title")])
      stop("NA values found in GSE123750 age stratum extraction")
    }
    
    # Print验收信息 for GSE123750
    cat(sprintf("GSE123750 - Total samples: %d\n", nrow(gse123750_data)))
    success_count_123750 <- sum(!is.na(gse123750_data$`Age stratum`))
    cat(sprintf("GSE123750 - Successfully parsed age strata: %d\n", success_count_123750))
    cat("GSE123750 - Age stratum counts:\n")
    print(table(gse123750_data$`Age stratum`))
  } else {
    cat("GSE123750 series matrix file not found.\n")
    stop("GSE123750 series matrix file not found")
  }
  
  # Process GSE43696
  cat("\nProcessing GSE43696...\n")
  gse43696_path <- file.path(base_dir, "data", "raw", "gse43696&GPL6480", "GSE43696_series_matrix.txt.gz")
  if (file.exists(gse43696_path)) {
    gse43696 <- getGEO(filename = gse43696_path, GSEMatrix = TRUE, getGPL = FALSE)
    gse43696_data <- extract_age_stratum_gse43696(gse43696)
    
    # Check for NA values
    na_count_43696 <- sum(is.na(gse43696_data$`Age stratum`))
    if (na_count_43696 > 0) {
      cat(sprintf("\nERROR: Found %d NA values in GSE43696 age stratum extraction\n", na_count_43696))
      cat("First 20 samples with NA age stratum:\n")
      na_samples <- gse43696_data %>% filter(is.na(`Age stratum`)) %>% head(20)
      print(na_samples[, c("geo_accession", "age_years", "Age stratum")])
      stop("NA values found in GSE43696 age stratum extraction")
    }
    
    # Print验收信息 for GSE43696
    cat(sprintf("GSE43696 - Total samples: %d\n", nrow(gse43696_data)))
    success_count_43696 <- sum(!is.na(gse43696_data$`Age stratum`))
    cat(sprintf("GSE43696 - Successfully parsed age strata: %d\n", success_count_43696))
    cat("GSE43696 - Age stratum counts:\n")
    print(table(gse43696_data$`Age stratum`))
  } else {
    cat("GSE43696 series matrix file not found.\n")
    stop("GSE43696 series matrix file not found")
  }
  
  # Combine results
  combined_data <- bind_rows(gse123750_data, gse43696_data)
  
  # Save to file
  write_csv(combined_data, age_stratum_map_file)
  cat(sprintf("\nSaved TableS8_age_stratum_map.csv: %s\n", age_stratum_map_file))
  
  # Load the generated map
  age_stratum_map <- combined_data
  
  # Memory cleanup
  rm(list = ls(pattern = "gse", all.names = TRUE))
  gc()
}

# Verify locked 11 genes
cat(sprintf("Number of locked genes: %d\n", nrow(locked_weights)))
cat("Locked genes (in file order):\n")
print(locked_weights$Gene)

# Calculate weights hash
weights_hash <- digest::digest(locked_weights, algo = "md5")
cat(sprintf("Weights hash (md5): %s\n", weights_hash))

# Robust Cohen's d function
calc_cohens_d <- function(g1, g2) {
  n1 <- length(g1); n2 <- length(g2)
  if (n1 < 2 || n2 < 2) return(c(NA, NA, NA))
  v1 <- var(g1); v2 <- var(g2)
  pooled_sd <- sqrt(((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2))
  if (pooled_sd == 0 || is.na(pooled_sd)) return(c(NA, NA, NA))
  d <- (mean(g1) - mean(g2)) / pooled_sd
  se <- sqrt((n1+n2)/(n1*n2) + (d^2)/(2*(n1+n2)))
  ci_low <- d - 1.96*se; ci_high <- d + 1.96*se
  return(c(d, ci_low, ci_high))
}

# Target cohorts that might contain age metadata
target_cohorts <- c("GSE123750", "GSE43696")
table_s8_results <- list()

# Tissue mapping
 tissue_mapping <- list(
  "GSE123750" = "Whole blood",
  "GSE43696" = "Bronchial epithelium"
)

cat("\n[2/3] Scanning cohorts for Age metadata and calculating efficacy...\n")
cat(sprintf("Target cohorts locked to: %s\n", paste(target_cohorts, collapse=", ")))
cat("Contrast is the cohort primary contrast (Table S1/S4); age is used only for stratification (Table S8).\n")

for (cohort in target_cohorts) {
  # Define file paths using processed_diet
  expr_path <- file.path(processed_diet_dir, paste0(cohort, "_diet_expr.rds"))
  pheno_path <- file.path(processed_diet_dir, paste0(cohort, "_pheno.rds"))
  
  # Check if files exist
  expr_exists <- file.exists(expr_path)
  pheno_exists <- file.exists(pheno_path)
  
  cat(sprintf("\nProcessing cohort: %s\n", cohort))
  cat(sprintf("  Expr path exists: %s\n", expr_exists))
  cat(sprintf("  Pheno path exists: %s\n", pheno_exists))
  
  if (!expr_exists || !pheno_exists) {
    cat(sprintf("  - WARNING: Missing required files for %s\n", cohort))
    next
  }
  
  # Load data
  expr <- readRDS(expr_path)
  pheno <- readRDS(pheno_path)
  
  # Print dimensions
  cat(sprintf("  dim(expr): %d x %d\n", nrow(expr), ncol(expr)))
  cat(sprintf("  nrow(pheno): %d\n", nrow(pheno)))
  
  # Print pheno column names
  cat("  names(pheno):\n")
  print(names(pheno))
  
  # (a) Sample ID column selection
  sample_id_candidates <- c("std_sample_id", "sample_id", "Sample_ID", "GSM", "geo_accession")
  sample_id_col <- sample_id_candidates[sample_id_candidates %in% names(pheno)][1]
  
  if (is.na(sample_id_col)) {
    cat(sprintf("  - ERROR: No valid sample ID column found in %s\n", cohort))
    stop(sprintf("No valid sample ID column found in %s", cohort))
  }
  
  cat(sprintf("  Chosen sample ID column: %s\n", sample_id_col))
  
  # (b) Define contrast based on series matrix data
  cat("  - Defining contrast based on series matrix data\n")
  
  # Define case/control based on cohort
  if (cohort == "GSE123750") {
    # For GSE123750, use severity contrast from title
    series_matrix_path <- file.path(base_dir, "data", "raw", "gse123750&GPL13158", "GSE123750_series_matrix.txt.gz")
    if (file.exists(series_matrix_path)) {
      gse <- getGEO(filename = series_matrix_path, GSEMatrix = TRUE, getGPL = FALSE)
      pheno_data <- gse@phenoData@data
      
      # Extract severity from title
      pheno_data$severity <- ifelse(grepl("severe", pheno_data$title, ignore.case = TRUE), "severe", 
                                   ifelse(grepl("mild/moderate", pheno_data$title, ignore.case = TRUE), "mild/moderate", NA))
      
      # Print title distribution
      cat("  Title distribution (first 10):\n")
      print(head(pheno_data$title, 10))
      
      # Merge with pheno by sample ID
      if ("geo_accession" %in% colnames(pheno_data)) {
        # Create mapping
        severity_map <- pheno_data %>% select(geo_accession, severity)
        
        # Merge with pheno
        pheno <- pheno %>% left_join(severity_map, by = c("std_sample_id" = "geo_accession"))
      }
      
      # Define case/control based on severity
      case_label <- "severe"
      control_label <- "mild/moderate"
      pheno$Case_Binary <- ifelse(pheno$severity == case_label, 1, 
                                 ifelse(pheno$severity == control_label, 0, NA))
      
      # Print severity distribution
      cat("  Severity distribution:\n")
      print(table(pheno$severity, useNA = "ifany"))
      
      # Print case/control labels
      cat(sprintf("  Case label: %s\n", case_label))
      cat(sprintf("  Control label: %s\n", control_label))
      
      # Print Case_Binary counts
      case_counts <- table(pheno$Case_Binary, useNA = "ifany")
      cat("  Case_Binary counts (0/1/NA):\n")
      print(case_counts)
      
      # Print contrast definition
      cat("  Contrast definition: severe vs mild/moderate (pre-school: wheeze; school age: asthma)\n")
    } else {
      stop(sprintf("Series matrix file not found: %s", series_matrix_path))
    }
  } else if (cohort == "GSE43696") {
    # For GSE43696, use disease state from series matrix
    series_matrix_path <- file.path(base_dir, "data", "raw", "gse43696&GPL6480", "GSE43696_series_matrix.txt.gz")
    if (file.exists(series_matrix_path)) {
      gse <- getGEO(filename = series_matrix_path, GSEMatrix = TRUE, getGPL = FALSE)
      pheno_data <- gse@phenoData@data
      
      # Extract disease state
      disease_state_col <- "disease state:ch1"
      if (disease_state_col %in% colnames(pheno_data)) {
        # Print disease state distribution
        cat(sprintf("  %s distribution:\n", disease_state_col))
        print(table(pheno_data[[disease_state_col]], useNA = "ifany"))
        
        # Map disease state to Asthma/Control
        pheno_data$std_group_label <- ifelse(grepl("Asthma", pheno_data[[disease_state_col]], ignore.case = TRUE), "Asthma", 
                                            ifelse(grepl("Control", pheno_data[[disease_state_col]], ignore.case = TRUE), "Control", NA))
        
        # Print std_group_label distribution
        cat("  std_group_label distribution:\n")
        print(table(pheno_data$std_group_label, useNA = "ifany"))
        
        # Set source
        pheno_data$std_group_label_source <- "series_matrix:disease state:ch1"
        
        # Create Case_Binary directly from pheno_data
        case_label <- "Asthma"
        control_label <- "Control"
        pheno_data$Case_Binary <- ifelse(pheno_data$std_group_label == case_label, 1, 
                                       ifelse(pheno_data$std_group_label == control_label, 0, NA))
        
        # Print case/control labels
        cat(sprintf("  Case label: %s\n", case_label))
        cat(sprintf("  Control label: %s\n", control_label))
        
        # Print Case_Binary counts
        case_counts <- table(pheno_data$Case_Binary, useNA = "ifany")
        cat("  Case_Binary counts (0/1/NA):\n")
        print(case_counts)
        
        # Print contrast definition
        cat("  Contrast definition: Asthma vs Control\n")
        
        # Merge with pheno by sample ID
        if ("geo_accession" %in% colnames(pheno_data)) {
          # Create mapping
          disease_map <- pheno_data %>% select(geo_accession, std_group_label, std_group_label_source, Case_Binary)
          
          # Merge with pheno
          pheno <- pheno %>% left_join(disease_map, by = c("std_sample_id" = "geo_accession"))
        }
      } else {
        stop(sprintf("Disease state column not found: %s", disease_state_col))
      }
    } else {
      stop(sprintf("Series matrix file not found: %s", series_matrix_path))
    }
  } else {
    stop(sprintf("Cohort %s not defined in Table S1/S4", cohort))
  }
  
  # Check for 1:1 split
  # if (length(case_counts) == 2 && case_counts[1] == case_counts[2]) {
  #   stop("Unaudited grouping detected: default split is not allowed for Table S8")
  # }
  
  # (c) Age stratum column selection
  age_candidates <- c("Age stratum", "Age_stratum", "age_stratum", "age_group", "Age_Group")
  age_col <- age_candidates[age_candidates %in% names(pheno)][1]
  
  # If no age stratum column found, try to use age stratum map
  if (is.na(age_col) && !is.null(age_stratum_map)) {
    cat("  - Using age stratum from TableS8_age_stratum_map.csv\n")
    # Extract samples for current cohort
    cohort_map <- age_stratum_map %>% filter(`Dataset (GEO)` == cohort)
    
    # Match samples by sample_id or geo_accession
    if (nrow(cohort_map) > 0) {
      # Create mapping from sample ID to age stratum
      if (sample_id_col %in% names(pheno)) {
        pheno_samples <- pheno[[sample_id_col]]
        
        # Create mapping using geo_accession from the map
        if ("geo_accession" %in% names(cohort_map)) {
          age_map <- cohort_map %>% 
            select(geo_accession, `Age stratum`) %>% 
            rename(!!sample_id_col := geo_accession)
          
          # Merge age stratum into pheno
          pheno <- pheno %>% 
            left_join(age_map, by = sample_id_col)
          
          age_col <- "Age stratum"
          
          # Calculate merge statistics
          n_merged_ok <- sum(!is.na(pheno$`Age stratum`))
          n_missing_after_merge <- sum(is.na(pheno$`Age stratum`))
          
          cat(sprintf("  Successfully merged age stratum for %d samples\n", n_merged_ok))
          cat(sprintf("  Missing age stratum for %d samples\n", n_missing_after_merge))
          
          # Stop if any missing
          if (n_missing_after_merge > 0) {
            cat("  First 20 samples with missing age stratum:\n")
            missing_samples <- pheno %>% filter(is.na(`Age stratum`)) %>% head(20)
            print(missing_samples[, c(sample_id_col)])
            stop("Missing age stratum after merge")
          }
          
          # Print age stratum counts
          cat("  Age stratum counts:\n")
          print(table(pheno$`Age stratum`))
        }
      }
    }
  }
  
  # If no age stratum column found, stop
  if (is.na(age_col)) {
    stop("No age stratum column found and age stratum map merge failed")
  }
  
  cat(sprintf("  Chosen age stratum column: %s\n", age_col))
  age_counts <- table(pheno[[age_col]])
  cat("  Age stratum counts:\n")
  print(age_counts)
  
  # Filter valid samples
  valid_pheno <- pheno %>% filter(!is.na(Case_Binary) & !is.na(!!sym(age_col)))
  
  # Print first few colnames of expr
  cat("  First 10 colnames of expr:\n")
  print(head(colnames(expr), 10))
  
  # Print first few sample IDs from valid_pheno
  cat("  First 10 sample IDs from valid_pheno:\n")
  print(head(valid_pheno[[sample_id_col]], 10))
  
  # Match samples
  common <- intersect(colnames(expr), valid_pheno[[sample_id_col]])
  cat(sprintf("  Number of common samples: %d\n", length(common)))
  
  # Print common samples if any
  if (length(common) > 0) {
    cat("  First 10 common samples:\n")
    print(head(common, 10))
  }
  
  if (length(common) <= 10) {
    cat(sprintf("  - ERROR: Insufficient common samples in %s\n", cohort))
    stop(sprintf("Insufficient common samples in %s", cohort))
  }
  
  # Subset data
  expr_sub <- expr[, common]
  valid_pheno <- valid_pheno[match(common, valid_pheno[[sample_id_col]]), ]
  
  # Compute score based on available genes
  avail_genes <- intersect(locked_weights$Gene, rownames(expr_sub))
  missing_genes <- setdiff(locked_weights$Gene, rownames(expr_sub))
  gene_coverage <- paste(length(avail_genes), "/", nrow(locked_weights), sep="")
  
  cat(sprintf("  Available genes: %s\n", paste(avail_genes, collapse=", ")))
  cat(sprintf("  Missing locked genes: %s\n", if(length(missing_genes) > 0) paste(missing_genes, collapse=", ") else "None"))
  cat(sprintf("  Gene coverage: %s\n", gene_coverage))
  cat(sprintf("  Weights hash: %s\n", weights_hash))
  
  w_sub <- locked_weights$Weight[match(avail_genes, locked_weights$Gene)]
  valid_pheno$Signature_Score <- colSums(expr_sub[avail_genes, , drop=FALSE] * w_sub, na.rm=TRUE)
  
  # Get age strata
  age_strata <- unique(valid_pheno[[age_col]])
  
  for (stratum in age_strata) {
    sub_df <- valid_pheno %>% filter(!!sym(age_col) == stratum)
    n_case <- sum(sub_df$Case_Binary == 1)
    n_control <- sum(sub_df$Case_Binary == 0)
    
    if (n_case >= 3 && n_control >= 3) {
      scores_case <- sub_df$Signature_Score[sub_df$Case_Binary == 1]
      scores_control <- sub_df$Signature_Score[sub_df$Case_Binary == 0]
      
      pval <- wilcox.test(scores_case, scores_control)$p.value
      d_res <- calc_cohens_d(scores_case, scores_control)
      
      # Calculate mean scores
      mean_score_case <- mean(scores_case)
      mean_score_control <- mean(scores_control)
      mean_diff <- mean_score_case - mean_score_control
      
      # Print detailed information to log
      cat(sprintf("  %s %s - Case_N: %d, Control_N: %d\n", cohort, stratum, n_case, n_control))
      cat(sprintf("  %s %s - mean_score_case: %.4f, mean_score_control: %.4f, mean_diff: %.4f\n", 
                  cohort, stratum, mean_score_case, mean_score_control, mean_diff))
      
      table_s8_results[[paste(cohort, stratum)]] <- data.frame(
        `Dataset (GEO)` = cohort,
        `Tissue / Compartment` = tissue_mapping[[cohort]],
        `Age stratum` = stratum,
        `Sample size (N)` = nrow(sub_df),
        `Effect size (Cohen's d)` = round(d_res[1], 3),
        `P value` = sprintf("%.2e", pval),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      cat(sprintf("  Skipped %s %s: insufficient class sizes (n_case=%d, n_control=%d)\n", cohort, stratum, n_case, n_control))
    }
  }
}

cat("\n[3/3] Formatting and Saving Table S8...\n")
table_s8_output <- file.path(supplement_dir, "TableS8_age_subgroups.csv")

if (length(table_s8_results) > 0) {
  table_s8 <- bind_rows(table_s8_results) %>% arrange(`Dataset (GEO)`, `Age stratum`)
  write_csv(table_s8, table_s8_output)
  cat(paste0("Saved Table S8: ", basename(table_s8_output), "\n"))
  
  # Log output using exact parameters from placeholder
  if(exists("log_output")) {
    log_output(
      file_path = table_s8_output,
      file_type = "table",
      figure_id = "Table S8",
      description = "Age subgroup analysis (Cohen's d and P-value)",
      script_name = "13_supp_tableS8_age.R",
      notes = "Extremely critical age subgroup defense table"
    )
  }
} else {
  cat("\n  ⚠ Warning: No sufficient age metadata found in any cohort to compute subgroups.\n")
  empty_df <- data.frame(
    `Dataset (GEO)`=character(), 
    `Tissue / Compartment`=character(),
    `Age stratum`=character(),
    `Sample size (N)`=integer(), 
    `Effect size (Cohen's d)`=numeric(), 
    `P value`=character(),
    check.names = FALSE
  )
  write_csv(empty_df, table_s8_output)
  cat(paste0("Saved Empty Table S8: ", basename(table_s8_output), "\n"))
}

cat("\n=== Supplementary Table S8 generation completed ===\n")

# Memory cleanup
rm(list = ls(pattern = "^(table|supplement|cohort|age|expr|pheno)", all.names = TRUE))
gc()

# Close logging
sink()