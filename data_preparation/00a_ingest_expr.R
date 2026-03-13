#!/usr/bin/env Rscript

# =========================================================================
# Script: 00a_ingest_expr.R
# Purpose: Find expression files in each cohort, read into matrix, determine ID type, write raw expression and diagnostic logs
# Author: Zuo
# Date: 2026-02-27
#
# Inputs:
#   - data/raw/              # Cohort directories with expression files
#
# Outputs:
#   - data/processed_full/expr_*__RAW.rds  # Raw expression matrices
#   - output/logs/00a_ingest_expr.log      # Log file
#   - output/logs/00_matrix_id_type.tsv    # Matrix ID type table
#   - output/logs/00_cohort_files.tsv      # Cohort files table
# =========================================================================

# Load required packages
library(tidyverse)
library(data.table)

# Set paths
# Get base directory from project root
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

raw_dir <- file.path(base_dir, "data", "raw")
processed_full_dir <- file.path(base_dir, "data", "processed_full")
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(processed_full_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_diet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00a_ingest_expr.log")
sink(log_file, append = TRUE)

cat("\n=== Starting expression data ingestion ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Raw directory: ", raw_dir, "\n"))
cat(paste0("Processed full directory: ", processed_full_dir, "\n"))
cat(paste0("Logs directory: ", logs_dir, "\n"))

# Function to check if series_matrix file is empty
check_series_matrix_empty <- function(file_path) {
  if (!grepl("series_matrix", file_path, ignore.case = TRUE)) {
    return(FALSE)
  }
  
  con <- gzfile(file_path, "rt")
  on.exit(close(con), add=TRUE)
  
  data_row_count <- 0
  in_data_section <- FALSE
  
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    
    # Check for Sample_data_row_count
    if (grepl("!Sample_data_row_count", line)) {
      # Extract the first numeric value after the key
      count_str <- sub("!Sample_data_row_count.*?([0-9]+).*", "\\1", line)
      count <- as.integer(count_str)
      if (!is.na(count) && count == 0) {
        return(TRUE)
      }
    }
    
    # Check for start of data section
    if (grepl("!series_matrix_table_begin", line)) {
      in_data_section <- TRUE
      next
    }
    
    # Count data rows after header (exclude header and end marker)
    if (in_data_section && nchar(line) > 0 && !grepl("^\"ID_REF\"", line) && !grepl("!series_matrix_table_end", line)) {
      data_row_count <- data_row_count + 1
    }
  }
  
  # If no Sample_data_row_count found, check actual data rows
  return(data_row_count == 0)
}

# Function to find expression files in cohort directory (prioritizing supplementary counts)
pick_expr_file_in_cohort <- function(cohort_dir) {
  files <- list.files(cohort_dir, recursive = TRUE, full.names = TRUE)
  bn <- basename(files)
  
  # Filter out unwanted files
  unwanted_patterns <- c("^9606\\.protein", "^KEGG_", "^pheno", "report", "profile", "annotation", "mapping", "id_symbol", "stats")
  for (pattern in unwanted_patterns) {
    files <- files[!grepl(pattern, bn, ignore.case = TRUE)]
    bn <- basename(files)
  }
  
  # Special case for gse152004&GPL11154: prioritize expr_norm
  if (basename(cohort_dir) == "gse152004&GPL11154") {
    expr_norm_files <- files[grepl("expr_norm.*\\.txt\\.gz$", bn, ignore.case = TRUE)]
    if (length(expr_norm_files) > 0) {
      return(list(file = expr_norm_files[1], type = "EXPR_NORM"))
    }
  }
  
  # First priority: supplementary counts files
  counts_patterns <- c(
    "genecounts.*\\.csv\\.gz",
    "raw_counts.*\\.txt\\.gz",
    "counts.*\\.csv\\.gz",
    "counts.*\\.txt\\.gz"
  )
  
  for (pattern in counts_patterns) {
    matches <- files[grepl(pattern, files, ignore.case = TRUE)]
    if (length(matches) > 0) {
      return(list(file = matches[1], type = "COUNTS"))
    }
  }
  
  # Second priority: processeddata files
  processeddata_patterns <- c(
    "RNAseq_.*processeddata\\.txt\\.gz$",
    "processeddata\\.txt\\.gz$"
  )
  
  for (pattern in processeddata_patterns) {
    matches <- files[grepl(pattern, files, ignore.case = TRUE)]
    if (length(matches) > 0) {
      return(list(file = matches[1], type = "PROCESSEDDATA"))
    }
  }
  
  # Third priority: series_matrix files (only if not empty)
  series_matrix_files <- files[grepl("_series_matrix\\.txt\\.gz$", files, ignore.case = TRUE) |
                               grepl("_series_matrix\\.txt$", files, ignore.case = TRUE)]
  
  for (file in series_matrix_files) {
    if (!check_series_matrix_empty(file)) {
      return(list(file = file, type = "SERIES_MATRIX"))
    }
  }
  
  # Third priority: other expression files
  other_expr_files <- files[grepl("_full\\.rds$", bn, ignore.case = TRUE) |
                           grepl("_expr\\.rds$", bn, ignore.case = TRUE) |
                           grepl("\\.txt\\.gz$", bn, ignore.case = TRUE) |
                           grepl("\\.tsv\\.gz$", bn, ignore.case = TRUE) |
                           grepl("\\.csv\\.gz$", bn, ignore.case = TRUE) |
                           grepl("\\.txt$", bn, ignore.case = TRUE) |
                           grepl("\\.tsv$", bn, ignore.case = TRUE) |
                           grepl("\\.csv$", bn, ignore.case = TRUE)]
  
  if (length(other_expr_files) > 0) {
    return(list(file = other_expr_files[1], type = "OTHER"))
  }
  
  return(list(file = NA_character_, type = NA_character_))
}

# Function to sniff separator from file
 sniff_sep <- function(file_path) {
   con <- if (grepl("\\.gz$", file_path)) gzfile(file_path, "rt") else file(file_path, "rt")
   on.exit(close(con), add=TRUE)
   # Skip lines starting with ! (GEO metadata)
   line <- ""
   while (length(line)==0 || grepl("^!", line)) {
     line <- readLines(con, n=1, warn=FALSE)
     if (length(line)==0) return("\t")
   }
   if (stringr::str_count(line, ",") > stringr::str_count(line, "\t")) "," else "\t"
 }

# Function to read expression file
read_expr_file <- function(file_path) {
  tryCatch({
    cat(paste0("Reading file: ", basename(file_path), "\n"))
    
    if (grepl("\\.rds$", file_path)) {
      # Read RDS file
      expr <- readRDS(file_path)
      if (is.matrix(expr)) {
        return(expr)
      } else if (is.data.frame(expr)) {
        return(as.matrix(expr))
      } else if (class(expr)[1] == "ExpressionSet") {
        # Extract expression matrix from ExpressionSet
        library(Biobase)
        return(exprs(expr))
      } else if (is.list(expr)) {
        # Try to find a matrix or data.frame in the list
        for (item in expr) {
          if (is.matrix(item) || is.data.frame(item)) {
            return(as.matrix(item))
          }
        }
      }
    } else if (grepl("GSE115770", file_path)) {
      # Special handling for GSE115770 (comma-separated, no header for first column)
      cat("Reading GSE115770 file...\n")
      
      # Use fread for faster reading
      library(data.table)
      df <- fread(file_path, sep=",", header=TRUE, stringsAsFactors=FALSE)
      
      # Get gene IDs from first column
      gene_ids <- as.character(df[[1]])
      
      # Extract expression data
      expr_data <- df[, -1, with=FALSE]
      
      # Convert to matrix
      mat <- as.matrix(expr_data)
      mode(mat) <- "numeric"
      
      # Set row and column names
      rownames(mat) <- gene_ids
      colnames(mat) <- colnames(expr_data)
      
      # Handle duplicate gene IDs by summing
      if (any(duplicated(gene_ids))) {
        cat("Found duplicate gene IDs, aggregating with sum\n")
        mat <- rowsum(mat, group = gene_ids, reorder = FALSE)
      }
      
      # Check if sample column names follow expected format
      if (ncol(mat) > 0) {
        sample_names <- colnames(mat)
        # Remove quotes from sample names
        sample_names <- gsub("^\"|\"$", "", sample_names)
        # Update column names
        colnames(mat) <- sample_names
        # Check if sample names are valid (not empty and contain alphanumeric characters)
        valid_sample_names <- sapply(sample_names, function(x) {
          nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
        })
        if (any(!valid_sample_names)) {
          invalid_cols <- which(!valid_sample_names)
          cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
          mat <- mat[, valid_sample_names, drop=FALSE]
        }
      }
      
      cat(paste0("Successfully read GSE115770 matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
      return(mat)
    } else if (basename(file_path) == "GSE118761_genecounts.csv.gz") {
      # Special handling for GSE118761_genecounts.csv.gz (white list special case)
      cat("Special handling for GSE118761_genecounts.csv.gz (white list)\n")
      
      # Use fread with gunzip to avoid column name issues
      library(data.table)
      library(R.utils)
      
      # Create temporary file for unzipped data
      temp_file <- tempfile(fileext = ".csv")
      on.exit(unlink(temp_file), add=TRUE)
      
      # Unzip file to temporary location (Windows compatible)
      cat("Unzipping file...\n")
      gunzip(file_path, destname = temp_file, remove = FALSE)
      
      # Read with fread, data.table = FALSE to get data.frame, header=TRUE to preserve column names
      cat("Reading unzipped file...\n")
      df <- fread(temp_file, data.table = FALSE)
      
      # Check file structure
      cat(paste0("File has ", nrow(df), " rows and ", ncol(df), " columns\n"))
      cat("First few rows:\n")
      print(head(df))
      
      # Identify gene ID column (last 3 columns are EnsemblID, Symbol, Description)
      cat("Identifying gene ID column\n")
      # Find EnsemblID column
      ensembl_col <- which(tolower(colnames(df)) == "ensemblid")
      if (length(ensembl_col) == 0) {
        # If EnsemblID column not found, use last column before description
        ensembl_col <- ncol(df) - 2
      }
      
      # Extract gene_id from EnsemblID column
      cat(paste0("Extracting gene_id from column: ", colnames(df)[ensembl_col], "\n"))
      gene_id <- df[[ensembl_col]]
      
      # Extract expression data (all columns except last 3 which are gene info)
      mat <- as.matrix(df[, 1:(ncol(df) - 3)])
      mode(mat) <- "numeric"
      
      # Clean gene_id conservatively
      cat("Cleaning gene_id...\n")
      valid_gene_id <- !is.na(gene_id) & gene_id != "" & gene_id != "0"
      cat(paste0("Removed ", sum(!valid_gene_id), " invalid gene_id (empty, NA, or '0')\n"))
      
      # Filter valid gene_id
      gene_id <- gene_id[valid_gene_id]
      mat <- mat[valid_gene_id, , drop = FALSE]
      
      # Set row names
      rownames(mat) <- gene_id
      
      # Handle duplicate gene IDs by summing
      if (any(duplicated(gene_id))) {
        cat("Found duplicate gene IDs, aggregating with sum\n")
        mat <- rowsum(mat, group = gene_id, reorder = FALSE)
        cat(paste0("Aggregated to ", nrow(mat), " unique genes\n"))
      }
      
      # Check if sample column names follow expected format
      if (ncol(mat) > 0) {
        sample_names <- colnames(mat)
        # Remove quotes from sample names
        sample_names <- gsub("^\"|\"$", "", sample_names)
        # Update column names
        colnames(mat) <- sample_names
        # Check if sample names are valid (not empty and contain alphanumeric characters)
        valid_sample_names <- sapply(sample_names, function(x) {
          nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
        })
        if (any(!valid_sample_names)) {
          invalid_cols <- which(!valid_sample_names)
          cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
          mat <- mat[, valid_sample_names, drop=FALSE]
        }
      }
      
      # Validate matrix size
      if (nrow(mat) < 20000) {
        stop("GSE118761 whitelist parsing failed: insufficient genes after processing (", nrow(mat), " < 20000)")
      }
      
      cat(paste0("Successfully read GSE118761 matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
      return(mat)
    } else if (grepl("processeddata", file_path, ignore.case = TRUE)) {
      # Special handling for processeddata files (e.g., GSE201955)
      cat("Reading processeddata file...\n")
      
      # Determine separator
      sep <- sniff_sep(file_path)
      
      # Read file
      if (grepl("\\.gz$", file_path)) {
        con <- gzfile(file_path, "rt")
        on.exit(close(con), add=TRUE)
        df <- read.table(con, header = TRUE, sep=sep, comment.char="!", quote="\"", stringsAsFactors = FALSE)
      } else {
        df <- read.table(file_path, header = TRUE, sep=sep, comment.char="!", quote="\"", stringsAsFactors = FALSE)
      }
      
      # Check if data frame is valid
      if (ncol(df) < 2) {
        cat("❌ File has insufficient columns\n")
        return(NULL)
      }
      
      # Identify gene ID column (ID_REF)
      cat("Identifying gene ID column\n")
      gene_id_col <- which(tolower(colnames(df)) == "id_ref")
      if (length(gene_id_col) == 0) {
        # If ID_REF column not found, use first column
        gene_id_col <- 1
      }
      
      # Extract gene IDs
      gene_id <- as.character(df[[gene_id_col]])
      
      # Extract expression data (all columns except gene ID column)
      mat <- as.matrix(df[, -gene_id_col, drop=FALSE])
      mode(mat) <- "numeric"
      
      # Set row names
      rownames(mat) <- gene_id
      
      # Handle duplicate gene IDs by summing
      if (any(duplicated(gene_id))) {
        cat("Found duplicate gene IDs, aggregating with sum\n")
        mat <- rowsum(mat, group = gene_id, reorder = FALSE)
      }
      
      # Check if sample column names follow expected format
      if (ncol(mat) > 0) {
        sample_names <- colnames(mat)
        # Remove quotes from sample names
        sample_names <- gsub("^\"|\"$", "", sample_names)
        # Update column names
        colnames(mat) <- sample_names
        # Check if sample names are valid (not empty and contain alphanumeric characters)
        valid_sample_names <- sapply(sample_names, function(x) {
          nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
        })
        if (any(!valid_sample_names)) {
          invalid_cols <- which(!valid_sample_names)
          cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
          mat <- mat[, valid_sample_names, drop=FALSE]
        }
      }
      
      # Sanity checks
      if (ncol(mat) < 2) {
        cat("❌ Matrix has insufficient columns after processing\n")
        return(NULL)
      }
      
      if (nrow(mat) < 1000) {
        cat("❌ Matrix has insufficient rows after processing\n")
        return(NULL)
      }
      
      na_percentage <- sum(is.na(mat)) / length(mat)
      if (na_percentage > 0.05) {
        cat(paste0("❌ High NA percentage: ", round(na_percentage * 100, 2), "%\n"))
        return(NULL)
      }
      
      cat(paste0("Successfully read processeddata matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
      return(mat)
    } else if (grepl("GSE118761", file_path)) {
      # General handling for other GSE118761 files
      con <- gzfile(file_path, "rt")
      on.exit(close(con), add=TRUE)
      
      # Read all lines
      lines <- readLines(con)
      
      # First line is header
      header <- strsplit(lines[1], ",")[[1]]
      
      # Find where gene info starts
      gene_info_start <- which(grepl("EnsemblID", header))[1]
      if (is.na(gene_info_start)) {
        gene_info_start <- length(header) - 2  # Last 3 columns
      }
      
      # Extract sample columns and gene info columns
      sample_cols <- header[1:(gene_info_start - 1)]
      
      # Read data lines
      data_lines <- lines[-1]
      
      # Parse each line
      data_list <- lapply(data_lines, function(line) {
        parts <- strsplit(line, ",")[[1]]
        counts <- as.numeric(parts[1:(gene_info_start - 1)])
        gene_id <- parts[gene_info_start]  # EnsemblID
        c(gene_id, counts)
      })
      
      # Convert to matrix
      data_matrix <- do.call(rbind, data_list)
      gene_ids <- data_matrix[, 1]
      mat <- matrix(as.numeric(data_matrix[, -1]), nrow = length(gene_ids), ncol = length(sample_cols))
      colnames(mat) <- sample_cols
      rownames(mat) <- gene_ids
      
      # Handle duplicate gene IDs by summing
      if (any(duplicated(gene_ids))) {
        cat("Found duplicate gene IDs, aggregating with sum\n")
        mat <- rowsum(mat, group = gene_ids, reorder = FALSE)
      }
      
      # Check for non-numeric columns (potential non-sample columns)
      non_numeric_cols <- c()
      for (i in 1:ncol(mat)) {
        col <- mat[, i]
        # Check if all values are numeric
        if (!all(!is.na(as.numeric(col)))) {
          non_numeric_cols <- c(non_numeric_cols, i)
        }
      }
      if (length(non_numeric_cols) > 0) {
        cat(paste0("Found ", length(non_numeric_cols), " non-numeric columns, removing them\n"))
        mat <- mat[, -non_numeric_cols, drop=FALSE]
      }
      
      # Check if sample column names follow expected format
      if (ncol(mat) > 0) {
        sample_names <- colnames(mat)
        # Remove quotes from sample names
        sample_names <- gsub("^\"|\"$", "", sample_names)
        # Update column names
        colnames(mat) <- sample_names
        # Check if sample names are valid (not empty and contain alphanumeric characters)
        valid_sample_names <- sapply(sample_names, function(x) {
          nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
        })
        if (any(!valid_sample_names)) {
          invalid_cols <- which(!valid_sample_names)
          cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
          mat <- mat[, valid_sample_names, drop=FALSE]
        }
      }
      
      cat(paste0("Successfully read GSE118761 matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
      return(mat)
    } else if (grepl("gse230048", file_path, ignore.case = TRUE)) {
      # Special handling for gse230048
      cat("Special handling for gse230048\n")
      
      # Use fread for faster and more reliable reading of CSV files
      library(data.table)
      cat("Reading file with fread...\n")
      df <- fread(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
      
      # Check if data frame is valid
      cat(paste0("File read successfully with ", nrow(df), " rows and ", ncol(df), " columns\n"))
      cat("Column names: ", paste(head(colnames(df)), collapse=", "), "\n")
      
      if (ncol(df) < 2) {
        cat("❌ File has insufficient columns\n")
        return(NULL)
      }
      
      # Check if data frame has any columns
      if (is.null(colnames(df))) {
        cat("❌ Data frame has no column names\n")
        return(NULL)
      }
      
      # Identify gene ID column
      cat("Identifying gene ID column\n")
      tryCatch({
        gene_id_col <- which(tolower(colnames(df)) == "gene_id")
        if (length(gene_id_col) == 0) {
          # If gene_id column not found, use first column
          gene_id_col <- 1
        }
        cat(paste0("Gene ID column identified: ", gene_id_col, " (", colnames(df)[gene_id_col], ")\n"))
      }, error = function(e) {
        cat(paste0("❌ Error identifying gene ID column: ", conditionMessage(e), "\n"))
        # Use first column as fallback
        gene_id_col <- 1
        cat(paste0("Using first column as gene ID column: ", gene_id_col, "\n"))
      })
      
      # Extract gene IDs
      gene_id <- as.character(df[[gene_id_col]])
      cat(paste0("Extracted ", length(gene_id), " gene IDs\n"))
      
      # Extract expression data (all columns except gene ID column)
      # Use data.table syntax to select all columns except gene_id_col
      expr_data <- df[, -c(gene_id_col), with=FALSE]
      mat <- as.matrix(expr_data)
      mode(mat) <- "numeric"
      cat(paste0("Expression matrix dimensions: ", nrow(mat), "x", ncol(mat), "\n"))
      
      # Check if gene_id length matches matrix rows
      if (length(gene_id) != nrow(mat)) {
        cat(paste0("❌ Mismatch: gene_id length (", length(gene_id), ") != matrix rows (", nrow(mat), ")\n"))
        # Use row numbers as gene IDs temporarily
        gene_id <- as.character(1:nrow(mat))
      }
      
      # Set row names
      rownames(mat) <- gene_id
      cat("Row names set successfully\n")
      
      # Handle duplicate gene IDs by summing
      if (any(duplicated(gene_id))) {
        cat("Found duplicate gene IDs, aggregating with sum\n")
        mat <- rowsum(mat, group = gene_id, reorder = FALSE)
        cat(paste0("Aggregated matrix dimensions: ", nrow(mat), "x", ncol(mat), "\n"))
      }
      
      # Check if sample column names follow expected format
      if (ncol(mat) > 0) {
        sample_names <- colnames(mat)
        # Remove quotes from sample names
        sample_names <- gsub("^\"|\"$", "", sample_names)
        # Update column names
        colnames(mat) <- sample_names
        # Check if sample names are valid (not empty and contain alphanumeric characters)
        valid_sample_names <- sapply(sample_names, function(x) {
          nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
        })
        if (any(!valid_sample_names)) {
          invalid_cols <- which(!valid_sample_names)
          cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
          mat <- mat[, valid_sample_names, drop=FALSE]
        }
      }
      
      # Sanity checks
      if (ncol(mat) < 2) {
        cat("❌ Matrix has insufficient columns after processing\n")
        return(NULL)
      }
      
      if (nrow(mat) < 1000) {
        cat("❌ Matrix has insufficient rows after processing\n")
        return(NULL)
      }
      
      na_percentage <- sum(is.na(mat)) / length(mat)
      if (na_percentage > 0.05) {
        cat(paste0("❌ High NA percentage: ", round(na_percentage * 100, 2), "%\n"))
        return(NULL)
      }
      
      cat(paste0("Successfully read gse230048 matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
      return(mat)
    } else if (basename(file_path) == "GSE152004_695_expr_norm.txt.gz") {
    # Special handling for GSE152004_695_expr_norm.txt.gz - hard rule to use first column as gene_id
    cat("Special handling for GSE152004_695_expr_norm.txt.gz - hard rule to use first column as gene_id\n")
    
    # Use data.table::fread() to read the file
    library(data.table)
    dt <- fread(file_path, header = TRUE)
    
    # Force gene_id_col = 1
    gene_id_col <- 1
    
    # Extract gene_id from first column
    gene_id <- as.character(dt[[gene_id_col]])
    
    # Remove first column from data table
    dt[, (gene_id_col) := NULL]
    
    # Convert to matrix and set to numeric
    mat <- as.matrix(dt)
    mode(mat) <- "numeric"
    
    # Set row names to uppercase trimmed gene symbols
    rownames(mat) <- toupper(trimws(gene_id))
    
    # Handle duplicate gene IDs
    if (any(duplicated(rownames(mat)))) {
      cat(paste0("Found duplicate gene IDs before aggregation: ", sum(duplicated(rownames(mat))), " duplicates\n"))
      cat("Aggregating duplicate gene IDs with sum\n")
      original_rows <- nrow(mat)
      mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
      cat(paste0("Aggregation complete: rows before=", original_rows, ", rows after=", nrow(mat), "\n"))
    }
    
    # Check if sample column names follow expected format
    if (ncol(mat) > 0) {
      sample_names <- colnames(mat)
      # Remove quotes from sample names
      sample_names <- gsub("^\"|\"$", "", sample_names)
      # Update column names
      colnames(mat) <- sample_names
      # Check if sample names are valid (not empty and contain alphanumeric characters)
      valid_sample_names <- sapply(sample_names, function(x) {
        nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
      })
      if (any(!valid_sample_names)) {
        invalid_cols <- which(!valid_sample_names)
        cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
        mat <- mat[, valid_sample_names, drop=FALSE]
      }
    }
    
    # Sanity checks
    if (ncol(mat) < 2) {
      cat("❌ Matrix has insufficient columns after processing\n")
      return(NULL)
    }
    
    if (nrow(mat) < 1000) {
      cat("❌ Matrix has insufficient rows after processing\n")
      return(NULL)
    }
    
    na_percentage <- sum(is.na(mat)) / length(mat)
    if (na_percentage > 0.05) {
      cat(paste0("❌ High NA percentage: ", round(na_percentage * 100, 2), "%\n"))
      return(NULL)
    }
    
    cat(paste0("Successfully read GSE152004_695_expr_norm.txt.gz matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
    return(mat)
  } else {
    # Determine separator by sniffing
    sep <- sniff_sep(file_path)
    cat(paste0("Detected separator: '", sep, "'\n"))
    
    # Read file
    if (grepl("\\.gz$", file_path)) {
      con <- gzfile(file_path, "rt")
      on.exit(close(con), add=TRUE)
      # For gse230048, explicitly use comma separator since it's a CSV file
      if (grepl("gse230048", file_path, ignore.case = TRUE)) {
        cat("Forcing comma separator for gse230048\n")
        df <- read.table(con, header = TRUE, sep=",", comment.char="!", quote="\"", stringsAsFactors = FALSE)
      } else {
        df <- read.table(con, header = TRUE, sep=sep, comment.char="!", quote="\"", stringsAsFactors = FALSE)
      }
    } else {
      # For gse230048, explicitly use comma separator since it's a CSV file
      if (grepl("gse230048", file_path, ignore.case = TRUE)) {
        cat("Forcing comma separator for gse230048\n")
        df <- read.table(file_path, header = TRUE, sep=",", comment.char="!", quote="\"", stringsAsFactors = FALSE)
      } else {
        df <- read.table(file_path, header = TRUE, sep=sep, comment.char="!", quote="\"", stringsAsFactors = FALSE)
      }
    }
    
    # Check if data frame is valid
    cat(paste0("File read successfully with ", nrow(df), " rows and ", ncol(df), " columns\n"))
    cat("Column names: ", paste(head(colnames(df)), collapse=", "), "\n")
    if (ncol(df) < 2) {
      cat("❌ File has insufficient columns\n")
      return(NULL)
    }
      
      # Extract gene IDs and expression data
      gene_id <- as.character(df[[1]])
      mat <- as.matrix(df[, -1, drop=FALSE])
      
      # Convert to numeric
      mode(mat) <- "numeric"
      
      # Handle duplicate gene IDs by summing
      if (any(duplicated(gene_id))) {
        cat("Found duplicate gene IDs, aggregating with sum\n")
        mat <- rowsum(mat, group = gene_id, reorder = FALSE)
      } else {
        rownames(mat) <- gene_id
      }
      
      # Check for non-numeric columns (potential non-sample columns)
      non_numeric_cols <- c()
      for (i in 1:ncol(mat)) {
        col <- mat[, i]
        # Check if all values are numeric
        if (!all(!is.na(as.numeric(col)))) {
          non_numeric_cols <- c(non_numeric_cols, i)
        }
      }
      if (length(non_numeric_cols) > 0) {
        cat(paste0("Found ", length(non_numeric_cols), " non-numeric columns, removing them\n"))
        mat <- mat[, -non_numeric_cols, drop=FALSE]
      }
      
      # Check if sample column names follow expected format
      if (ncol(mat) > 0) {
        sample_names <- colnames(mat)
        # Remove quotes from sample names
        sample_names <- gsub("^\"|\"$", "", sample_names)
        # Update column names
        colnames(mat) <- sample_names
        # Check if sample names are valid (not empty and contain alphanumeric characters)
        valid_sample_names <- sapply(sample_names, function(x) {
          nchar(x) > 0 && grepl("[A-Za-z0-9]", x)
        })
        if (any(!valid_sample_names)) {
          invalid_cols <- which(!valid_sample_names)
          cat(paste0("Found ", length(invalid_cols), " columns with invalid sample names, removing them\n"))
          mat <- mat[, valid_sample_names, drop=FALSE]
        }
      }
      
      # Sanity checks
      if (ncol(mat) < 2) {
        cat("❌ Matrix has insufficient columns after processing\n")
        return(NULL)
      }
      
      if (nrow(mat) < 1000) {
        cat("❌ Matrix has insufficient rows after processing\n")
        return(NULL)
      }
      
      na_percentage <- sum(is.na(mat)) / length(mat)
      if (na_percentage > 0.05) {
        cat(paste0("❌ High NA percentage: ", round(na_percentage * 100, 2), "%\n"))
        return(NULL)
      }
      
      cat(paste0("Successfully read matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
      return(mat)
    }
    
    cat("❌ Unsupported file format\n")
    return(NULL)
  }, error = function(e) {
    cat("❌ Error reading file: ", conditionMessage(e), "\n")
    return(NULL)
  }
)
}

# Function to score how gene-like an ID is
gene_like_score <- function(ids) {
  score <- 0
  
  is_probe <- function(x) {
    grepl("^ILMN_", x) ||
    grepl("^A_\\d+_", x) ||                 # Agilent
    grepl("_at$|_s_at$|_x_at$|_a_at$", x) || # Affy
    grepl("^\\d{6,}$", x)                   # 纯数字长串（常见探针/探针集）
  }
  
  is_entrez <- function(x) {
    grepl("^[0-9]+$", x)  # 纯数字
  }
  
  is_ensg <- function(x) {
    grepl("^ENSG", x)  # ENSG IDs
  }
  
  is_symbol <- function(x) {
    # 基因符号模式：以大写字母开头，后面跟着字母、数字或连字符，但不是GSM格式
    grepl("^[A-Z][A-Z0-9-]{1,}$", x) && !grepl("^GSM", x, ignore.case = TRUE)  # 基因符号
  }
  
  for (id in ids) {
    id <- as.character(id)
    if (is_probe(id) || is_entrez(id) || is_ensg(id) || is_symbol(id)) {
      score <- score + 1
    }
  }
  
  return(score / length(ids))
}

# Function to score how sample-like an ID is
sample_like_score <- function(ids) {
  score <- 0
  
  is_sample <- function(x) {
    grepl("^GSM\\d+", x, ignore.case = TRUE) ||  # GEO sample IDs
    grepl("^AF\\d+$", x) ||  # AF001 style
    grepl("^lib\\d+", x, ignore.case = TRUE) ||  # lib11081 style
    grepl("^sample", x, ignore.case = TRUE) ||  # sample prefix
    grepl("^s\\d+$", x, ignore.case = TRUE) ||  # s1, s2 style
    grepl("^QU\\d+$", x, ignore.case = TRUE)  # QU01 style
  }
  
  for (id in ids) {
    id <- as.character(id)
    if (is_sample(id)) {
      score <- score + 1
    }
  }
  
  return(score / length(ids))
}

# Function to check if a column name is sample-like
is_sample_like_column <- function(col_name) {
  # Check if column name matches sample patterns
  return(grepl("^GSM\\d+$", col_name, ignore.case = TRUE) ||
         grepl("^AF\\d+$", col_name) ||
         grepl("^lib\\d+", col_name, ignore.case = TRUE) ||
         grepl("^QU\\d+$", col_name, ignore.case = TRUE) ||
         grepl("^HR\\d+$", col_name, ignore.case = TRUE) ||  # HR sample IDs (GSE152004)
         grepl("^SJ\\d+$", col_name, ignore.case = TRUE))  # SJ sample IDs (GSE152004)
}

# Function to identify gene ID column in data frame
identify_gene_id_column <- function(df) {
  # Separate columns into sample-like and candidate columns
  sample_like_cols <- c()
  candidate_cols <- c()
  
  for (i in 1:ncol(df)) {
    if (is_sample_like_column(colnames(df)[i])) {
      sample_like_cols <- c(sample_like_cols, i)
    } else {
      candidate_cols <- c(candidate_cols, i)
    }
  }
  
  cat(paste0("Found ", length(sample_like_cols), " sample-like columns and ", length(candidate_cols), " candidate columns\n"))
  
  # If no candidate columns, gene_id must be in row names or first column (special case)
  if (length(candidate_cols) == 0) {
    cat("No candidate columns found, checking if gene_id is in row names or first column\n")
    # Check if row names are gene-like
    if (!is.null(rownames(df)) && is_gene_like_rownames(rownames(df))) {
      cat("Gene-like row names found, assuming gene_id is in row names\n")
      # Return 0 to indicate gene_id is in row names
      return(0)
    } else {
      # Check if first column could be gene_id despite sample-like name
      # This is a special case for files like GSE118761 where first column is gene_id but has sample-like name
      cat("No gene-like row names found, checking first column as potential gene_id\n")
      # Calculate gene-like score for first column values
      first_col_values <- as.character(df[[1]])[1:min(1000, nrow(df))]
      first_col_score <- gene_like_score(first_col_values)
      
      if (first_col_score >= 0.1) {
        cat("First column has gene-like values despite sample-like name, using it as gene_id column\n")
        return(1)
      } else {
        # No gene_id found, stop with error
        stop("No candidate gene_id column found. All columns appear to be sample columns. Please manually specify the gene_id column.")
      }
    }
  }
  
  # Priority A: Explicit gene ID column names in all columns
  gene_id_cols <- c("", "gene", "geneid", "gene_id", "id_ref", "id", "feature_id",
                    "ensg", "ensembl", "ensembl_gene_id", 
                    "entrez", "entrezid", "symbol", "gene_symbol")
  
  colnames_lower <- tolower(colnames(df))
  for (col in gene_id_cols) {
    for (i in 1:ncol(df)) {
      if (col == colnames_lower[i]) {
        cat(paste0("Found explicit gene ID column: ", colnames(df)[i], "\n"))
        return(i)
      }
    }
  }
  
  # Priority B: Content-based detection among candidate columns
  max_gene_score <- 0
  best_col <- candidate_cols[1]  # Start with first candidate column
  
  for (i in candidate_cols) {
    # Get first 1000 values
    values <- as.character(df[[i]])[1:min(1000, nrow(df))]
    # Calculate gene-like score
    score <- gene_like_score(values)
    if (score > max_gene_score) {
      max_gene_score <- score
      best_col <- i
    }
  }
  
  cat(paste0("Selected gene ID column based on content: ", colnames(df)[best_col], " with score: ", round(max_gene_score, 3), "\n"))
  return(best_col)
}

# Function to check if row names are gene-like
is_gene_like_rownames <- function(rownames) {
  if (is.null(rownames)) return(FALSE)
  # Check first 200 row names
  sample_rownames <- rownames[1:min(200, length(rownames))]
  # Calculate gene-like score
  score <- gene_like_score(sample_rownames)
  # Return TRUE if score is high enough (≥ 0.1, which means at least 20 out of 200 are gene-like)
  return(score >= 0.1)
}

# Function to check if row names are just row numbers
are_rownames_row_numbers <- function(rownames) {
  if (is.null(rownames)) return(FALSE)
  # Check if all row names are numeric and match row indices
  numeric_rownames <- suppressWarnings(as.numeric(rownames))
  if (any(is.na(numeric_rownames))) return(FALSE)
  # Check if they match 1:n
  return(all(numeric_rownames == 1:length(numeric_rownames)))
}

# Function to standardize expression matrix shape
standardize_expr_shape <- function(input) {
  cat("=== Standardizing expression matrix shape ===\n")
  
  # Handle different input types
  if (is.matrix(input)) {
    # Check if row names are already gene-like
    if (!is.null(rownames(input)) && is_gene_like_rownames(rownames(input))) {
      # Input is already a matrix with gene-like row names, use them as gene IDs
      cat("Input is a matrix with gene-like row names, using them as gene IDs\n")
      gene_ids <- rownames(input)
      mat <- input
    } else if (is.null(rownames(input)) || are_rownames_row_numbers(rownames(input))) {
      # Row names are empty or just row numbers, need to identify gene ID column from data frame
      cat("Matrix row names are empty or row numbers, converting to data frame to find gene ID column\n")
      # Convert to data frame
      df <- as.data.frame(input)
      # Identify gene ID column
      gene_id_col <- identify_gene_id_column(df)
      
      if (gene_id_col == 0) {
        # Gene ID is in row names
        cat("Gene ID is in row names\n")
        gene_ids <- rownames(df)
        mat <- as.matrix(df)
      } else {
        cat(paste0("Identified gene ID column: ", colnames(df)[gene_id_col], "\n"))
        # Extract gene IDs and expression data
        gene_ids <- as.character(df[[gene_id_col]])
        expr_data <- df[, -gene_id_col, drop = FALSE]
        # Convert back to matrix
        mat <- as.matrix(expr_data)
      }
    } else {
      # Input is a matrix but row names are not gene-like, use them anyway
      cat("Input is a matrix, using row names as gene IDs\n")
      gene_ids <- rownames(input)
      mat <- input
    }
  } else if (is.data.frame(input)) {
    # Check if we need to identify gene ID column
    if (is.null(rownames(input)) || are_rownames_row_numbers(rownames(input))) {
      # Input is a data frame, identify gene ID column
      gene_id_col <- identify_gene_id_column(input)
      
      if (gene_id_col == 0) {
        # Gene ID is in row names
        cat("Gene ID is in row names\n")
        gene_ids <- rownames(input)
        # Convert to matrix
        mat <- as.matrix(input)
      } else {
        cat(paste0("Identified gene ID column: ", colnames(input)[gene_id_col], "\n"))
        # Extract gene IDs and expression data
        gene_ids <- as.character(input[[gene_id_col]])
        expr_data <- input[, -gene_id_col, drop = FALSE]
        # Convert to matrix
        mat <- as.matrix(expr_data)
      }
    } else if (is_gene_like_rownames(rownames(input))) {
      # Data frame has gene-like row names, use them as gene IDs
      cat("Data frame has gene-like row names, using them as gene IDs\n")
      gene_ids <- rownames(input)
      # Convert to matrix
      mat <- as.matrix(input)
    } else {
      # Data frame has row names but they're not gene-like, identify gene ID column
      gene_id_col <- identify_gene_id_column(input)
      cat(paste0("Identified gene ID column: ", colnames(input)[gene_id_col], "\n"))
      # Extract gene IDs and expression data
      gene_ids <- as.character(input[[gene_id_col]])
      expr_data <- input[, -gene_id_col, drop = FALSE]
      # Convert to matrix
      mat <- as.matrix(expr_data)
    }
  } else {
    stop("Input must be a matrix or data frame")
  }
  
  # Convert to numeric if not already
  if (mode(mat) != "numeric") {
    mode(mat) <- "numeric"
  }
  
  # Step 4: Check NA percentage (skip for large matrices to save time)
  if (nrow(mat) < 50000) {
    na_count <- sum(is.na(mat))
    na_rate <- na_count / length(mat)
    cat(paste0("NA rate after numeric conversion: ", round(na_rate * 100, 2), "%\n"))
    
    if (na_rate > 0.05) {
      # Find first 5 columns with NA values
      na_cols <- which(colSums(is.na(mat)) > 0)[1:5]
      if (length(na_cols) > 0) {
        cat("First 5 columns with NA values: ", paste(colnames(mat)[na_cols], collapse=", "), "\n")
      }
      stop("High NA rate (>5%) after numeric conversion")
    }
  } else {
    cat("Skipping NA rate check for large matrix to save time\n")
  }
  
  # Step 5: Set row names
  rownames(mat) <- gene_ids
  
  # Step 6: Check and fix transposition (skip for large matrices to save time)
  if (nrow(mat) < 1000) {
    if (nrow(mat) < 100 && ncol(mat) > 1000) {
      # Likely transposed
      cat("Detected potential transposition: rows < 100 and columns > 1000\n")
      mat <- t(mat)
      cat(paste0("Transposed matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
    } else {
      # Score gene-like vs sample-like for rows and columns
      row_score <- gene_like_score(rownames(mat)[1:min(200, nrow(mat))])
      col_score <- gene_like_score(colnames(mat)[1:min(200, ncol(mat))])
      
      cat(paste0("Row gene-like score: ", round(row_score, 3), "\n"))
      cat(paste0("Column gene-like score: ", round(col_score, 3), "\n"))
      
      if (col_score > row_score * 1.5) {
        # Columns are more gene-like, transpose
        cat("Columns are more gene-like than rows, transposing\n")
        mat <- t(mat)
        cat(paste0("Transposed matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
        
        # Recheck after transposition
        new_row_score <- gene_like_score(rownames(mat)[1:min(200, nrow(mat))])
        new_col_score <- gene_like_score(colnames(mat)[1:min(200, ncol(mat))])
        
        if (new_col_score > new_row_score) {
          stop("Failed to correct transposition: columns still more gene-like than rows")
        }
      }
    }
  } else {
    cat("Skipping transposition check for large matrix to save time\n")
  }
  
  # Step 7: Handle duplicate gene IDs
  if (any(duplicated(rownames(mat)))) {
    cat("Found duplicate gene IDs, aggregating with sum\n")
    # For large matrices, use a more efficient approach
    if (nrow(mat) > 50000) {
      cat("Using efficient aggregation for large matrix\n")
      # Convert to data.table for faster aggregation
      library(data.table)
      dt <- data.table(mat)
      dt[, gene_id := rownames(mat)]
      dt_agg <- dt[, lapply(.SD, sum), by = gene_id]
      rownames(dt_agg) <- dt_agg$gene_id
      mat <- as.matrix(dt_agg[, -"gene_id", with = FALSE])
    } else {
      mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
    }
  }
  
  # Step 8: Validate matrix
  if (nrow(mat) < 2000) {
    stop(paste0("Insufficient rows: ", nrow(mat), " < 2000"))
  }
  
  if (ncol(mat) < 5) {
    stop(paste0("Insufficient columns: ", ncol(mat), " < 5"))
  }
  
  if (any(rownames(mat) == "")) {
    stop("Empty row names found")
  }
  
  if (any(duplicated(rownames(mat)))) {
    stop("Duplicate row names found after aggregation")
  }
  
  if (any(colnames(mat) == "")) {
    stop("Empty column names found")
  }
  
  if (any(duplicated(colnames(mat)))) {
    stop("Duplicate column names found")
  }
  
  cat(paste0("Standardized matrix: dim=", nrow(mat), "x", ncol(mat), "\n"))
  return(mat)
}

# Function to determine ID type
determine_id_type <- function(expr) {
  if (!is.matrix(expr)) {
    return("UNKNOWN")
  }
  
  # Get first 20 row names to determine ID type
  ids <- rownames(expr)[1:min(20, nrow(expr))]
  
  # Count gene-like vs probe-like IDs
  gene_count <- 0
  entrez_count <- 0
  ensg_count <- 0
  probe_count <- 0
  
  is_probe <- function(x) {
    grepl("^ILMN_", x) ||
    grepl("^A_\\d+_", x) ||                 # Agilent
    grepl("_at$|_s_at$|_x_at$|_a_at$", x) || # Affy
    grepl("^\\d{6,}$", x)                   # 纯数字长串（常见探针/探针集）
  }
  
  is_entrez <- function(x) {
    grepl("^[0-9]+$", x)  # 纯数字
  }
  
  is_ensg <- function(x) {
    grepl("^ENSG", x)  # ENSG IDs
  }
  
  is_symbol <- function(x) {
    # 基因符号模式：以大写字母开头，后面跟着字母、数字或连字符，但不是GSM格式
    grepl("^[A-Z][A-Z0-9-]{1,}$", x) && !grepl("^GSM", x, ignore.case = TRUE)  # 基因符号
  }
  
  for (id in ids) {
    id <- as.character(id)
    if (is_probe(id)) {
      probe_count <- probe_count + 1
    } else if (is_entrez(id)) {
      entrez_count <- entrez_count + 1
    } else if (is_ensg(id)) {
      ensg_count <- ensg_count + 1
    } else if (is_symbol(id)) {
      gene_count <- gene_count + 1
    }
  }
  
  # Determine the dominant ID type
  if (entrez_count > 0 && entrez_count >= max(ensg_count, gene_count, probe_count)) {
    return("ENTREZ")
  } else if (ensg_count > 0 && ensg_count >= max(entrez_count, gene_count, probe_count)) {
    return("ENSG")
  } else if (gene_count > 0 && gene_count >= max(entrez_count, ensg_count, probe_count)) {
    return("SYMBOL")
  } else if (probe_count > 0) {
    return("PROBE")
  } else {
    return("UNKNOWN")
  }
}

# Function to extract GSE/GPL from path
extract_gse_gpl_from_path <- function(path) {
  cohort <- basename(dirname(path))  # e.g. "gse65204&GPL14550"
  gse <- regmatches(cohort, regexpr("gse[0-9]+", cohort, ignore.case=TRUE))
  gpl <- regmatches(cohort, regexpr("GPL[0-9]+", cohort))
  gse <- if (length(gse)==0) NA_character_ else toupper(gse)
  gpl <- if (length(gpl)==0) NA_character_ else gpl
  list(gse=gse, gpl=gpl, cohort=cohort)
}

# Get cohort directories
cohort_dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
cohort_dirs <- cohort_dirs[grepl("gse", basename(cohort_dirs), ignore.case = TRUE)]

# Process each cohort
matrix_id_type_df <- data.frame()
cohort_files_df <- data.frame()

for (cohort_dir in cohort_dirs) {
  cat(paste0("\n=== Processing cohort: ", basename(cohort_dir), " ===\n"))
  
  # Pick expression file
  expr_file_info <- pick_expr_file_in_cohort(cohort_dir)
  if (is.na(expr_file_info$file)) {
    cat("❌ No expression file found\n")
    next
  }
  
  # Read expression file
  expr <- read_expr_file(expr_file_info$file)
  if (is.null(expr)) {
    cat("❌ Failed to read expression file\n")
    next
  }
  
  # Extract GSE/GPL info
  path_info <- extract_gse_gpl_from_path(expr_file_info$file)
  
  # Standardize expression matrix shape
  standardize_success <- TRUE
  tryCatch({
    expr <- standardize_expr_shape(as.data.frame(expr))
  }, error = function(e) {
    cat("❌ Error standardizing expression shape: ", conditionMessage(e), "\n")
    standardize_success <- FALSE
  })
  
  if (!standardize_success) {
    next
  }
  
  # Determine ID type
  id_type <- determine_id_type(expr)
  
  # Print and log RAW matrix information
  cat(paste0("=== RAW Matrix Information ===\n"))
  cat(paste0("Cohort: ", path_info$cohort, "\n"))
  cat(paste0("Dimensions: ", nrow(expr), "x", ncol(expr), "\n"))
  cat(paste0("Head of rownames: ", paste(head(rownames(expr)), collapse=", "), "\n"))
  cat(paste0("Head of colnames: ", paste(head(colnames(expr)), collapse=", "), "\n"))
  
  # Write raw expression matrix to processed_full
  raw_output_file <- file.path(processed_full_dir, paste0("expr_", path_info$cohort, "__RAW.rds"))
  saveRDS(expr, raw_output_file)
  cat("Saved raw expression: ", basename(raw_output_file), " | dim=", nrow(expr), "x", ncol(expr), "\n")
  
  # Update matrix ID type table
  matrix_id_type_df <- rbind(matrix_id_type_df, data.frame(
    gse = path_info$gse,
    gpl = path_info$gpl,
    cohort = path_info$cohort,
    id_type = id_type,
    file_type = expr_file_info$type,
    stringsAsFactors = FALSE
  ))
  
  # Update cohort files table
  cohort_files_df <- rbind(cohort_files_df, data.frame(
    cohort = path_info$cohort,
    gse = path_info$gse,
    gpl = path_info$gpl,
    file_path = expr_file_info$file,
    file_type = expr_file_info$type,
    stringsAsFactors = FALSE
  ))
  
  # Log FULL expr source information
  cat(paste0("=== FULL expr source information ===\n"))
  cat(paste0("Cohort: ", path_info$cohort, "\n"))
  cat(paste0("Source file: ", expr_file_info$file, "\n"))
  cat(paste0("File type: ", expr_file_info$type, "\n"))
  cat(paste0("ID type: ", id_type, "\n"))
  cat(paste0("Sample count: ", ncol(expr), "\n"))
  
  # Determine data type (counts/TPM/FPKM/normalized)
  data_type <- "UNKNOWN"
  file_name <- basename(expr_file_info$file)
  if (grepl("count", file_name, ignore.case = TRUE)) {
    data_type <- "COUNTS"
  } else if (grepl("tpm", file_name, ignore.case = TRUE)) {
    data_type <- "TPM"
  } else if (grepl("fpkm", file_name, ignore.case = TRUE)) {
    data_type <- "FPKM"
  } else if (grepl("norm", file_name, ignore.case = TRUE) || grepl("processed", file_name, ignore.case = TRUE)) {
    data_type <- "NORMALIZED"
  }
  cat(paste0("Data type: ", data_type, "\n"))
  
  # Check if multiple files were concatenated (based on cohort name)
  is_concatenated <- grepl("&", path_info$cohort)
  cat(paste0("Is concatenated: ", ifelse(is_concatenated, "YES", "NO"), "\n"))
  if (is_concatenated) {
    cat(paste0("Platforms: ", gsub("gse[0-9]+&", "", path_info$cohort, ignore.case = TRUE), "\n"))
  }
  cat(paste0("=================================\n"))
  
  # Memory cleanup
  rm(expr)
  gc()
  Sys.sleep(1)
}

# Save diagnostic tables
matrix_id_type_output <- file.path(logs_dir, "00_matrix_id_type.tsv")
write.table(matrix_id_type_df, matrix_id_type_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved matrix ID type table: ", basename(matrix_id_type_output), "\n")

cohort_files_output <- file.path(logs_dir, "00_cohort_files.tsv")
write.table(cohort_files_df, cohort_files_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved cohort files table: ", basename(cohort_files_output), "\n")

cat("\n=== Expression data ingestion completed ===\n")
sink()

# Memory cleanup
rm(list = ls(pattern = "^(matrix|expr|data|probe|gpl|cohort)", all.names = TRUE))
gc()