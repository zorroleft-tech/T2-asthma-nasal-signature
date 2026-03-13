# =========================================================================
# Script: 00c_prep_pheno_anchors_65204_201955_45111.R
# Purpose: Prepare standardized phenotype data for three anchor cohorts (GSE65204, GSE201955, GSE45111)
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/raw/*/series_matrix.txt.gz  # Series matrix files for each GSE
#
# Outputs:
#   - data/processed_full/pheno_GSE65204__FULL.rds  # Phenotype data for GSE65204
#   - data/processed_full/pheno_GSE201955__FULL.rds  # Phenotype data for GSE201955
#   - data/processed_full/pheno_GSE45111__FULL.rds  # Phenotype data for GSE45111
#   - output/logs/00c_prep_pheno_anchors.log  # Log file
# =========================================================================

# Load required libraries
# library(GEOquery)

# Disable GEOquery network requests (local file mode)
# options('download.file.method.GEOquery' = 'auto')
# options('GEOquery.inmemory.gpl' = FALSE)

# Define functions

# Normalize column names
norm_names <- function(names) {
  names <- tolower(names)
  names <- gsub('[^a-z0-9_]', '_', names)
  names <- gsub('^_|_$', '', names)
  names <- gsub('__+', '_', names)
  return(names)
}

# Select matching column from candidate list
pick_col <- function(df, candidates) {
  col_names <- colnames(df)
  for (candidate in candidates) {
    matches <- grep(paste0('^', candidate, '$'), col_names, ignore.case = TRUE)
    if (length(matches) > 0) {
      return(col_names[matches[1]])
    }
  }
  return(NULL)
}

# Safe conversion to numeric
as_numeric_safe <- function(x) {
  if (is.character(x)) {
    x <- gsub('[^0-9.-]', '', x)
  }
  return(as.numeric(x))
}

# Parse series_matrix.txt.gz file
parse_series_matrix <- function(file_path) {
  # Read the file
  con <- gzfile(file_path, "r")
  lines <- readLines(con)
  close(con)
  
  # Extract sample IDs
  sample_id_line <- grep("^!Sample_geo_accession", lines)
  if (length(sample_id_line) == 0) {
    sample_id_line <- grep("^!Sample_title", lines)
  }
  sample_ids <- strsplit(lines[sample_id_line], "\t")[[1]]
  sample_ids <- sample_ids[2:length(sample_ids)]  # Remove the first element which is the key
  
  # Extract characteristics
  char_lines <- grep("^!Sample_characteristics", lines)
  
  # Create data frame
  pheno_data <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  
  # Track column names to avoid overwriting
  col_name_count <- list()
  
  # Process each characteristic line
  for (line in lines[char_lines]) {
    # Split the line by tabs
    parts <- strsplit(line, "\t")[[1]]
    # Extract values for each sample
    values <- parts[2:length(parts)]
    # Clean values and extract characteristic name from the first value
    if (length(values) > 0) {
      first_value <- gsub("\"", "", values[1])
      original_char_name <- sub(":.*", "", first_value)
      original_char_name <- trimws(original_char_name)
      # Normalize characteristic name
      char_name <- norm_names(original_char_name)
      
      # Handle duplicate column names
      if (char_name %in% names(col_name_count)) {
        col_name_count[[char_name]] <- col_name_count[[char_name]] + 1
        char_name <- paste0(char_name, "_", sprintf("%02d", col_name_count[[char_name]]))
      } else {
        col_name_count[[char_name]] <- 1
      }
      
      # Clean all values
      values <- gsub("\"", "", values)
      # Use original char name to clean values
      values <- sub(paste0(original_char_name, ":\\s*"), "", values)
      values <- trimws(values)
      # Add to data frame
      pheno_data[[char_name]] <- values
    }
  }
  
  return(pheno_data)
}

# Extract values from pheno data
extract_from_pheno <- function(pheno_data, pattern) {
  # Look for columns matching the pattern
  matching_cols <- grep(pattern, colnames(pheno_data), ignore.case = TRUE)
  
  if (length(matching_cols) == 0) return(rep(NA, nrow(pheno_data)))
  
  # Use the first matching column
  values <- as.character(pheno_data[[matching_cols[1]]])
  return(values)
}

# Main function
main <- function() {
  # Define paths
  # Get base directory (project root)
  base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  # Verify we're running from project root
  if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation")))) {
    stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
  }
  
  # Define directories
  raw_dir <- file.path(base_dir, "data", "raw")
  output_dir <- file.path(base_dir, "data", "processed_full")
  logs_dir <- file.path(base_dir, "output", "logs")
  
  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set up logging
  log_file <- file.path(logs_dir, "00c_prep_pheno_anchors.log")
  sink(log_file, append = TRUE)
  
  # Log base directory and paths
  cat("\n=== Starting phenotype preparation for anchor cohorts ===\n")
  cat(paste0("Execution time: ", Sys.time(), "\n"))
  cat(paste0("Base directory: ", base_dir, "\n"))
  cat(paste0("Raw directory: ", raw_dir, "\n"))
  cat(paste0("Output directory: ", output_dir, "\n"))
  cat(paste0("Logs directory: ", logs_dir, "\n"))
  
  # Define GSE list
  gse_list <- c('GSE65204', 'GSE201955', 'GSE45111')
  
  # Process each GSE
  for (gse_id in gse_list) {
    cat(paste0('\nProcessing ', gse_id, '...\n'))
    
    # Search for pheno files
    files <- list.files(raw_dir, pattern = gse_id, recursive = TRUE, full.names = TRUE)
    
    # Prioritize series_matrix.txt.gz files
    series_matrix_files <- grep('series_matrix\\.txt\\.gz', files, value = TRUE)
    
    if (length(series_matrix_files) > 0) {
      cat(paste0('Found series_matrix file: ', series_matrix_files[1], '\n'))
      # Parse series_matrix file
      cat('Parsing series_matrix file...\n')
      pheno_data <- parse_series_matrix(series_matrix_files[1])
      cat('Parsed pheno data dimensions: ', nrow(pheno_data), ' x ', ncol(pheno_data), '\n')
      
      # Extract sample_id
      sample_id_col <- pick_col(pheno_data, c('sample_id', 'geo_accession', 'gsm', 'title'))
      if (is.null(sample_id_col)) {
        stop(paste0(gse_id, ': sample_id column not found'))
      }
      sample_id <- pheno_data[[sample_id_col]]
      
      # Check if sample_id is unique
      if (any(duplicated(sample_id))) {
        stop(paste0(gse_id, ': sample_id contains duplicate values'))
      }
      
      # Initialize result data frame
      result <- data.frame(
        sample_id = sample_id,
        dataset = gse_id,
        age_years = NA,
        sex = NA,
        asthma_status = NA,
        site = NA,
        race_ethnicity = NA,
        stringsAsFactors = FALSE
      )
      
      # Extract common fields
      # age_years
      age_col <- pick_col(pheno_data, c('age', 'age_years', 'age_yrs'))
      if (!is.null(age_col)) {
        result$age_years <- as_numeric_safe(pheno_data[[age_col]])
      } else {
        # Try to extract from pheno data
        age_values <- extract_from_pheno(pheno_data, 'age')
        result$age_years <- as_numeric_safe(age_values)
      }
      
      # sex
      # For GSE201955, use broader extraction and normalization
      if (gse_id == 'GSE201955') {
        # Look for columns containing sex or gender
        sex_cols <- grep('sex|gender', colnames(pheno_data), ignore.case = TRUE)
        if (length(sex_cols) > 0) {
          for (col in sex_cols) {
            sex_values <- as.character(pheno_data[[col]])
            # Clean and normalize values
            x_chr <- trimws(tolower(sex_values))
            
            # Map text values first
            result$sex[x_chr %in% c('female', 'f')] <- 'Female'
            result$sex[x_chr %in% c('male', 'm')] <- 'Male'
            
            # Only handle numeric values if they look like numbers
            numeric_mask <- grepl('^[0-9]+(\\.[0-9]+)?$', x_chr)
            if (any(numeric_mask)) {
              num_values <- as.numeric(x_chr[numeric_mask])
              result$sex[numeric_mask & num_values == 2] <- 'Female'
              result$sex[numeric_mask & num_values == 1] <- 'Male'
            }
          }
        } else {
          # Try to extract from pheno data
          sex_values <- as.character(extract_from_pheno(pheno_data, 'sex|gender'))
          # Clean and normalize values
          x_chr <- trimws(tolower(sex_values))
          
          # Map text values first
          result$sex[x_chr %in% c('female', 'f')] <- 'Female'
          result$sex[x_chr %in% c('male', 'm')] <- 'Male'
          
          # Only handle numeric values if they look like numbers
          numeric_mask <- grepl('^[0-9]+(\\.[0-9]+)?$', x_chr)
          if (any(numeric_mask)) {
            num_values <- as.numeric(x_chr[numeric_mask])
            result$sex[numeric_mask & num_values == 2] <- 'Female'
            result$sex[numeric_mask & num_values == 1] <- 'Male'
          }
        }
        # Print unique sex values for debugging
        if (length(sex_cols) > 0) {
          unique_sex <- unique(pheno_data[[sex_cols[1]]])
          cat(paste0('GSE201955 unique sex values: ', paste(unique_sex, collapse = ', '), '\n'))
        }
      } else {
        # For other GSEs, use original logic
        sex_col <- pick_col(pheno_data, c('sex', 'gender'))
        if (!is.null(sex_col)) {
          sex_values <- tolower(pheno_data[[sex_col]])
          result$sex[sex_values %in% c('female', 'f')] <- 'Female'
          result$sex[sex_values %in% c('male', 'm')] <- 'Male'
        } else {
          # Try to extract from pheno data
          sex_values <- tolower(extract_from_pheno(pheno_data, 'sex|gender'))
          result$sex[sex_values %in% c('female', 'f')] <- 'Female'
          result$sex[sex_values %in% c('male', 'm')] <- 'Male'
        }
      }
      
      # asthma_status
      # Try to find asthma status from various possible columns
      asthma_candidates <- c('asthma', 'asthma_status', 'disease_status', 'case_control', 'status', 'disease')
      
      # For GSE65204, prioritize extracting from characteristics with keyword matching
      if (gse_id == 'GSE65204') {
        # Extract asthma status from characteristics
        for (pattern in asthma_candidates) {
          # Look for columns containing the pattern
          matching_cols <- grep(pattern, colnames(pheno_data), ignore.case = TRUE)
          if (length(matching_cols) > 0) {
            for (col in matching_cols) {
              asthma_values <- tolower(pheno_data[[col]])
              # Map values to standard categories
              result$asthma_status[asthma_values %in% c('asthma', 'yes', '1', 'true', 'case')] <- 'Asthma'
              result$asthma_status[asthma_values %in% c('control', 'no', '0', 'false', 'control')] <- 'Control'
            }
            break
          }
        }
      } else {
        # For other GSEs, use original logic
        asthma_col <- pick_col(pheno_data, asthma_candidates)
        
        if (!is.null(asthma_col)) {
          asthma_values <- tolower(pheno_data[[asthma_col]])
          # Map values to standard categories
          result$asthma_status[asthma_values %in% c('asthma', 'yes', '1', 'true', 'case')] <- 'Asthma'
          result$asthma_status[asthma_values %in% c('control', 'no', '0', 'false', 'control')] <- 'Control'
        } else {
          # Try to extract from pheno data using multiple patterns
          for (pattern in asthma_candidates) {
            asthma_values <- tolower(extract_from_pheno(pheno_data, pattern))
            if (any(!is.na(asthma_values))) {
              result$asthma_status[asthma_values %in% c('asthma', 'yes', '1', 'true', 'case')] <- 'Asthma'
              result$asthma_status[asthma_values %in% c('control', 'no', '0', 'false', 'control')] <- 'Control'
              break
            }
          }
        }
      }
      
      # Special handling for GSE45111 (if only asthma samples)
      if (gse_id == 'GSE45111') {
        # If no asthma status found, assume all are Asthma
        if (all(is.na(result$asthma_status))) {
          result$asthma_status <- 'Asthma'
        }
      }
      
      # site
      site_col <- pick_col(pheno_data, c('site', 'center', 'location'))
      if (!is.null(site_col)) {
        result$site <- as.character(pheno_data[[site_col]])
      } else {
        # Try to extract from pheno data
        result$site <- extract_from_pheno(pheno_data, 'site|center|location')
      }
      
      # race_ethnicity
      race_col <- pick_col(pheno_data, c('race', 'ethnicity', 'race_ethnicity', 'self_reported_race_ethnicity'))
      if (!is.null(race_col)) {
        result$race_ethnicity <- as.character(pheno_data[[race_col]])
      } else {
        # Try to extract from pheno data
        result$race_ethnicity <- extract_from_pheno(pheno_data, 'race|ethnicity')
      }
      
      # Process specific GSE anchor fields
      if (gse_id == 'GSE65204') {
        # lnige
        ige_col <- pick_col(pheno_data, c('lnige', 'ige', 'log_ige'))
        if (!is.null(ige_col)) {
          result$lnige <- as_numeric_safe(pheno_data[[ige_col]])
          result$lnige_raw <- pheno_data[[ige_col]]  # Keep original column
        } else {
          # Try to extract from pheno data
          ige_values <- extract_from_pheno(pheno_data, 'ige|lnige')
          result$lnige <- as_numeric_safe(ige_values)
          result$lnige_raw <- ige_values  # Keep original column
        }
        
        # Check if lnige exists
        if (all(is.na(result$lnige))) {
          stop('GSE65204: lnige field not found')
        }
        
        # Check if asthma_status is not all missing
        stopifnot(mean(is.na(result$asthma_status)) < 0.01)
        
      } else if (gse_id == 'GSE201955') {
        # feno
        feno_col <- pick_col(pheno_data, c('feno', 'fractional_exhaled_nitric_oxide'))
        if (!is.null(feno_col)) {
          result$feno <- as_numeric_safe(pheno_data[[feno_col]])
          result$feno_raw <- pheno_data[[feno_col]]  # Keep original column
        } else {
          # Try to extract from pheno data
          feno_values <- extract_from_pheno(pheno_data, 'feno|fractional_exhaled_nitric_oxide')
          result$feno <- as_numeric_safe(feno_values)
          result$feno_raw <- feno_values  # Keep original column
        }
        
        # Check if feno exists
        if (all(is.na(result$feno))) {
          stop('GSE201955: feno field not found')
        }
        
        # Check if sex is not mostly missing
        stopifnot(mean(is.na(result$sex)) < 0.05)
        
        # Extract optional covariates
        # bmi
        bmi_col <- pick_col(pheno_data, c('bmi', 'body_mass_index'))
        if (!is.null(bmi_col)) {
          result$bmi <- as_numeric_safe(pheno_data[[bmi_col]])
        } else {
          bmi_values <- extract_from_pheno(pheno_data, 'bmi|body_mass_index')
          result$bmi <- as_numeric_safe(bmi_values)
        }
        
        # fev1_abs, fvc_abs, fev1_pct
        fev1_abs_col <- pick_col(pheno_data, c('fev1', 'fev1_abs', 'forced_expiratory_volume_in_1_second'))
        if (!is.null(fev1_abs_col)) {
          result$fev1_abs <- as_numeric_safe(pheno_data[[fev1_abs_col]])
        } else {
          fev1_abs_values <- extract_from_pheno(pheno_data, 'fev1|forced_expiratory_volume_in_1_second')
          result$fev1_abs <- as_numeric_safe(fev1_abs_values)
        }
        
        fvc_abs_col <- pick_col(pheno_data, c('fvc', 'fvc_abs', 'forced_vital_capacity'))
        if (!is.null(fvc_abs_col)) {
          result$fvc_abs <- as_numeric_safe(pheno_data[[fvc_abs_col]])
        } else {
          fvc_abs_values <- extract_from_pheno(pheno_data, 'fvc|forced_vital_capacity')
          result$fvc_abs <- as_numeric_safe(fvc_abs_values)
        }
        
        fev1_pct_col <- pick_col(pheno_data, c('fev1_pct', 'fev1_percent', 'fev1_predicted'))
        if (!is.null(fev1_pct_col)) {
          result$fev1_pct <- as_numeric_safe(pheno_data[[fev1_pct_col]])
        } else {
          fev1_pct_values <- extract_from_pheno(pheno_data, 'fev1_pct|fev1_percent|fev1_predicted')
          result$fev1_pct <- as_numeric_safe(fev1_pct_values)
        }
        
        # bal_neut, bal_eo
        bal_neut_col <- pick_col(pheno_data, c('bal_neut', 'bal_neutrophils', 'bronchoalveolar_lavage_neutrophils'))
        if (!is.null(bal_neut_col)) {
          result$bal_neut <- as_numeric_safe(pheno_data[[bal_neut_col]])
        } else {
          bal_neut_values <- extract_from_pheno(pheno_data, 'bal_neut|bal_neutrophils|bronchoalveolar_lavage_neutrophils')
          result$bal_neut <- as_numeric_safe(bal_neut_values)
        }
        
        bal_eo_col <- pick_col(pheno_data, c('bal_eo', 'bal_eosinophils', 'bronchoalveolar_lavage_eosinophils'))
        if (!is.null(bal_eo_col)) {
          result$bal_eo <- as_numeric_safe(pheno_data[[bal_eo_col]])
        } else {
          bal_eo_values <- extract_from_pheno(pheno_data, 'bal_eo|bal_eosinophils|bronchoalveolar_lavage_eosinophils')
          result$bal_eo <- as_numeric_safe(bal_eo_values)
        }
        
        # blood_eo_abs
        blood_eo_col <- pick_col(pheno_data, c('blood_eo', 'blood_eosinophils', 'blood_eo_abs'))
        if (!is.null(blood_eo_col)) {
          result$blood_eo_abs <- as_numeric_safe(pheno_data[[blood_eo_col]])
        } else {
          blood_eo_values <- extract_from_pheno(pheno_data, 'blood_eo|blood_eosinophils')
          result$blood_eo_abs <- as_numeric_safe(blood_eo_values)
        }
        
        # current_smok
        current_smok_col <- pick_col(pheno_data, c('current_smok', 'smoking', 'smoker'))
        if (!is.null(current_smok_col)) {
          current_smok_values <- tolower(pheno_data[[current_smok_col]])
          result$current_smok[current_smok_values %in% c('yes', 'y', '1', 'current')] <- 'Y'
          result$current_smok[current_smok_values %in% c('no', 'n', '0', 'never', 'former')] <- 'N'
        } else {
          current_smok_values <- tolower(extract_from_pheno(pheno_data, 'smok|smoking|smoker'))
          result$current_smok[current_smok_values %in% c('yes', 'y', '1', 'current')] <- 'Y'
          result$current_smok[current_smok_values %in% c('no', 'n', '0', 'never', 'former')] <- 'N'
        }
        
        # ancestry_pc1, ancestry_pc2
        pc1_col <- pick_col(pheno_data, c('ancestry_pc1', 'pc1', 'principal_component_1'))
        if (!is.null(pc1_col)) {
          result$ancestry_pc1 <- as_numeric_safe(pheno_data[[pc1_col]])
        } else {
          pc1_values <- extract_from_pheno(pheno_data, 'pc1|principal_component_1')
          result$ancestry_pc1 <- as_numeric_safe(pc1_values)
        }
        
        pc2_col <- pick_col(pheno_data, c('ancestry_pc2', 'pc2', 'principal_component_2'))
        if (!is.null(pc2_col)) {
          result$ancestry_pc2 <- as_numeric_safe(pheno_data[[pc2_col]])
        } else {
          pc2_values <- extract_from_pheno(pheno_data, 'pc2|principal_component_2')
          result$ancestry_pc2 <- as_numeric_safe(pc2_values)
        }
        
        # ige (raw value)
        ige_raw_col <- pick_col(pheno_data, c('ige', 'immunoglobulin_e'))
        if (!is.null(ige_raw_col)) {
          result$ige <- as_numeric_safe(pheno_data[[ige_raw_col]])
        } else {
          ige_raw_values <- extract_from_pheno(pheno_data, 'ige|immunoglobulin_e')
          result$ige <- as_numeric_safe(ige_raw_values)
        }
        
        # maternal_asthma
        maternal_asthma_col <- pick_col(pheno_data, c('maternal_asthma', 'mother_asthma'))
        if (!is.null(maternal_asthma_col)) {
          result$maternal_asthma <- as.character(pheno_data[[maternal_asthma_col]])
        } else {
          result$maternal_asthma <- extract_from_pheno(pheno_data, 'maternal_asthma|mother_asthma')
        }
        
        # recruitment_source
        recruitment_source_col <- pick_col(pheno_data, c('recruitment_source', 'recruitment', 'source'))
        if (!is.null(recruitment_source_col)) {
          result$recruitment_source <- as.character(pheno_data[[recruitment_source_col]])
        } else {
          result$recruitment_source <- extract_from_pheno(pheno_data, 'recruitment_source|recruitment|source')
        }
        
        # step
        step_col <- pick_col(pheno_data, c('step', 'study_step'))
        if (!is.null(step_col)) {
          result$step <- as_numeric_safe(pheno_data[[step_col]])
        } else {
          step_values <- extract_from_pheno(pheno_data, 'step|study_step')
          result$step <- as_numeric_safe(step_values)
        }
        
        # tissue (fixed value)
        result$tissue <- 'primary bronchial epithelial cells'
        
      } else if (gse_id == 'GSE45111') {
        # asthma_phenotype
        phenotype_col <- pick_col(pheno_data, c('asthma_phenotype', 'phenotype', 'disease_phenotype'))
        if (!is.null(phenotype_col)) {
          result$asthma_phenotype <- as.character(pheno_data[[phenotype_col]])
          result$asthma_phenotype_raw <- pheno_data[[phenotype_col]]  # Keep original column
        } else {
          # Try to extract from pheno data
          phenotype_values <- extract_from_pheno(pheno_data, 'phenotype|asthma_phenotype')
          result$asthma_phenotype <- phenotype_values
          result$asthma_phenotype_raw <- phenotype_values  # Keep original column
        }
        
        # Check if asthma_phenotype exists
        if (all(is.na(result$asthma_phenotype))) {
          stop('GSE45111: asthma_phenotype field not found')
        }
        
        # Try to extract sputum_eos or sputum_eos_pct
        sputum_eos_col <- pick_col(pheno_data, c('sputum_eos', 'sputum_eosinophils', 'sputum_eos_pct'))
        if (!is.null(sputum_eos_col)) {
          result$sputum_eos <- as_numeric_safe(pheno_data[[sputum_eos_col]])
          result$sputum_eos_raw <- pheno_data[[sputum_eos_col]]  # Keep original column
        } else {
          # Try to extract from pheno data
          sputum_eos_values <- extract_from_pheno(pheno_data, 'sputum_eos|sputum_eosinophils')
          result$sputum_eos <- as_numeric_safe(sputum_eos_values)
          result$sputum_eos_raw <- sputum_eos_values  # Keep original column
        }
        
        # Extract optional covariates
        # smoking_pack_years
        smoking_col <- pick_col(pheno_data, c('smoking_pack_years', 'pack_years', 'smoking'))
        if (!is.null(smoking_col)) {
          result$smoking_pack_years <- as_numeric_safe(pheno_data[[smoking_col]])
        } else {
          smoking_values <- extract_from_pheno(pheno_data, 'smoking_pack_years|pack_years|smoking')
          result$smoking_pack_years <- as_numeric_safe(smoking_values)
        }
        
        # sample_type (fixed value)
        result$sample_type <- 'sputum'
        
        # disease_state (fixed value)
        result$disease_state <- 'asthma'
        
      }
      
      # Calculate missing rates for key columns
      key_columns <- c('sample_id', 'age_years', 'sex', 'asthma_status')
      if (gse_id == 'GSE65204') key_columns <- c(key_columns, 'lnige')
      if (gse_id == 'GSE201955') key_columns <- c(key_columns, 'feno')
      if (gse_id == 'GSE45111') key_columns <- c(key_columns, 'asthma_phenotype')
      
      missing_rates <- sapply(result[, key_columns], function(x) mean(is.na(x)) * 100)
      
      # Print log
      cat(paste0('\n', gse_id, ' processing results:\n'))
      cat(paste0('Sample count: ', nrow(result), '\n'))
      cat('Missing rates for key columns:\n')
      print(missing_rates)
      
      # Write output file
      output_file <- file.path(output_dir, paste0('pheno_', gse_id, '__FULL.rds'))
      saveRDS(result, file = output_file)
      cat(paste0('\nOutput file written to: ', output_file, '\n'))
      
    } else {
      stop(paste0(gse_id, ': series_matrix.txt.gz file not found'))
    }
  }
  
  cat('\nAll GSEs processed!\n')
}

# Execute main function
main()
