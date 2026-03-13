#!/usr/bin/env Rscript

# =========================================================================
# Script: 00b_id_standardize_and_map.R
# Purpose: Standardize IDs and map probes to genes based on platform strategy, write standardized expression and mapping files
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/processed_full/expr_*.rds  # Raw expression files for each cohort
#   - data_preparation/platform_strategy.tsv  # Platform strategy table
#   - data_preparation/gse_strategy.tsv  # GSE strategy table
#
# Outputs:
#   - data/processed_full/expr_*__FULL.rds  # Standardized expression files
#   - data/processed_full/id_map_*.tsv.gz  # ID mapping tables
#   - data/processed_full/probe_mapping_*.rds  # Probe mapping files
#   - output/logs/00b_id_standardize_and_map.log  # Log file
# =========================================================================

# Load required packages
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
processed_full_dir <- file.path(base_dir, "data", "processed_full")
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
logs_dir <- file.path(base_dir, "output", "logs")
platform_strategy_file <- file.path(base_dir, "data_preparation", "platform_strategy.tsv")

# Create output directories if not exist
dir.create(processed_full_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_diet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00b_id_standardize_and_map.log")
sink(log_file, append = TRUE)

# Log base directory and paths
cat("\n=== Starting ID standardization and mapping ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Input files: \n"))
cat(paste0("  - Platform strategy: ", platform_strategy_file, "\n"))
cat(paste0("  - GSE strategy: ", file.path(base_dir, "data_preparation", "gse_strategy.tsv"), "\n"))
cat(paste0("Output directories: \n"))
cat(paste0("  - Processed full: ", processed_full_dir, "\n"))
cat(paste0("  - Processed diet: ", processed_diet_dir, "\n"))
cat(paste0("  - Logs: ", logs_dir, "\n"))

# Read platform strategy table
platform_strategy <- read_tsv(platform_strategy_file)
cat("Loaded platform strategy table:\n")
print(platform_strategy)

# Read GSE strategy table
gse_strategy_file <- file.path(base_dir, "data_preparation", "gse_strategy.tsv")
gse_strategy <- read_tsv(gse_strategy_file)
cat("Loaded GSE strategy table:\n")
print(gse_strategy)

# Function to read raw expression file
read_raw_expr <- function(cohort) {
  raw_file <- file.path(processed_full_dir, paste0("expr_", cohort, "__RAW.rds"))
  if (!file.exists(raw_file)) {
    cat(paste0("❌ Raw expression file not found for cohort: ", cohort, "\n"))
    stop("Input file not found: ", raw_file)
  }
  readRDS(raw_file)
}

# Function to clean ENSG IDs (remove version)
clean_ensg <- function(ensg_ids) {
  sub("\\..*$", "", ensg_ids)
}

# Function to clean Entrez IDs
clean_entrez <- function(entrez_ids) {
  # Remove non-numeric characters
  gsub("[^0-9]", "", entrez_ids)
}

# Function to calculate HGNC symbol hit rate using org.Hs.eg.db
calculate_hgnc_hit_rate <- function(ids) {
  if (length(ids) == 0) {
    return(0)
  }
  
  # Get a representative sample
  sample_size <- min(2000, length(ids))
  sample_ids <- sample(ids, sample_size)
  sample_ids <- as.character(sample_ids)
  
  suppressPackageStartupMessages(library(AnnotationDbi))
  
  tryCatch({
    # Use mapIds to check if IDs are valid symbols, handle multiple mappings
    result <- mapIds(org.Hs.eg.db, keys = sample_ids, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    valid_symbols <- names(result[!is.na(result)])
    hit_rate <- length(valid_symbols) / length(sample_ids)
    return(hit_rate)
  }, error = function(e) {
    # If AnnotationDbi fails, return 0
    return(0)
  })
}

# Function to determine ID type
determine_id_type <- function(ids) {
  if (length(ids) == 0) {
    return("UNKNOWN")
  }
  
  # Get a representative sample of IDs
  sample_size <- min(2000, length(ids))
  sample_ids <- sample(ids, sample_size)
  sample_ids <- as.character(sample_ids)
  
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
    grepl("^[0-9]+$", x)  # 纯数字，全匹配
  }
  
  is_ensg <- function(x) {
    grepl("^ENSG", x)  # ENSG IDs
  }
  
  for (id in sample_ids) {
    id <- as.character(id)
    if (is_entrez(id)) {
      entrez_count <- entrez_count + 1
    } else if (is_probe(id)) {
      probe_count <- probe_count + 1
    } else if (is_ensg(id)) {
      ensg_count <- ensg_count + 1
    } else {
      gene_count <- gene_count + 1
    }
  }
  
  # Calculate percentages
  total_count <- length(sample_ids)
  gene_pct <- gene_count / total_count
  entrez_pct <- entrez_count / total_count
  ensg_pct <- ensg_count / total_count
  probe_pct <- probe_count / total_count
  
  # Determine the dominant ID type with threshold
  if (entrez_pct > 0.5) {
    return("ENTREZ")
  } else if (ensg_pct > 0.5) {
    return("ENSG")
  } else if (probe_pct > 0.5) {
    return("PROBE")
  } else if (gene_pct > 0.3) {
    # Use org.Hs.eg.db to validate if these are real HGNC symbols
    suppressPackageStartupMessages(library(AnnotationDbi))
    
    tryCatch({
        # Use mapIds to check if IDs are valid symbols, handle multiple mappings
        result <- mapIds(org.Hs.eg.db, keys = sample_ids, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
        valid_symbols <- names(result[!is.na(result)])
        hit_rate <- length(valid_symbols) / length(sample_ids)
        
        if (hit_rate > 0.3) {
          return("SYMBOL")
        } else {
          return("PROBE")
        }
      }, error = function(e) {
        # If AnnotationDbi fails, fall back to probe
        return("PROBE")
      })
  } else {
    # If no type dominates, use the original logic
    max_count <- max(ensg_count, gene_count, probe_count, entrez_count)
    if (entrez_count == max_count && entrez_count > 0) {
      return("ENTREZ")
    } else if (ensg_count == max_count && ensg_count > 0) {
      return("ENSG")
    } else if (gene_count == max_count && gene_count > 0) {
      # Validate gene symbols with org.Hs.eg.db
      suppressPackageStartupMessages(library(AnnotationDbi))
      
      tryCatch({
        result <- mapIds(org.Hs.eg.db, keys = sample_ids, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
        valid_symbols <- names(result[!is.na(result)])
        hit_rate <- length(valid_symbols) / length(sample_ids)
        
        if (hit_rate > 0.3) {
          return("SYMBOL")
        } else {
          return("PROBE")
        }
      }, error = function(e) {
        return("PROBE")
      })
    } else if (probe_count > 0) {
      return("PROBE")
    } else {
      return("UNKNOWN")
    }
  }
}

# Function to build probe mapping for GPL platforms
build_probe_mapping <- function(gpl_id, key_col, symbol_col) {
  mapping_output <- file.path(processed_full_dir, paste0("probe_mapping_", gpl_id, ".rds"))
  
  # Check if mapping already exists (cache priority)
  if (file.exists(mapping_output)) {
    existing_map <- readRDS(mapping_output)
    if (nrow(existing_map) > 100) {
      # For GPL23961, check if it has the required columns
      if (gpl_id == "GPL23961") {
        if ("spot_id" %in% tolower(colnames(existing_map)) && "orf" %in% tolower(colnames(existing_map))) {
          cat(paste0("✓ Mapping already exists for GPL", gpl_id, " (with required columns)\n"))
          return(existing_map)
        } else {
          cat(paste0("⚠️  Existing mapping for GPL", gpl_id, " missing required columns, rebuilding\n"))
        }
      } else {
        cat(paste0("✓ Mapping already exists for GPL", gpl_id, "\n"))
        return(existing_map)
      }
    }
  }
  
  # Offline parsing for GPL6104
  if (gpl_id == "GPL6104") {
    annot_file <- file.path(base_dir, "data", "raw", "GPL", "GPL6104.annot.gz")
    if (file.exists(annot_file)) {
      cat("Special handling for GPL6104: Using offline GPL6104.annot.gz\n")
      
      tryCatch({
        # Read the file and find the platform table
        con <- gzfile(annot_file, "r")
        lines <- readLines(con)
        close(con)
        
        # Remove empty lines
        lines <- lines[lines != ""]
        
        # Find the platform table begin
        table_start <- which(grepl("!platform_table_begin", lines))[1]
        if (!is.na(table_start)) {
          # Find the header line (first line after table start that doesn't start with !)
          header_idx <- table_start + 1
          while (header_idx <= length(lines) && grepl("^!", lines[header_idx])) {
            header_idx <- header_idx + 1
          }
          
          if (header_idx <= length(lines)) {
            # Get table data
            table_data <- lines[header_idx:length(lines)]
            
            # Write to a temporary file to handle line endings properly
            temp_file <- tempfile(fileext = ".txt")
            writeLines(table_data, temp_file)
            
            # Read the table with fill=TRUE to handle inconsistent line lengths
            tbl <- read.table(temp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
            unlink(temp_file)
            
            # Use ID and Gene symbol columns
            if ("ID" %in% colnames(tbl) && "Gene.symbol" %in% colnames(tbl)) {
              # Use Gene.symbol column (note the dot instead of space)
              probe_map <- data.frame(
                ID = as.character(tbl$ID),
                Gene_Symbol = as.character(tbl$Gene.symbol),
                stringsAsFactors = FALSE
              )
            } else if ("ID" %in% colnames(tbl) && "Gene symbol" %in% colnames(tbl)) {
              # Use Gene symbol column (with space)
              probe_map <- data.frame(
                ID = as.character(tbl$ID),
                Gene_Symbol = as.character(tbl$`Gene symbol`),
                stringsAsFactors = FALSE
              )
            } else {
              cat("❌ Required columns not found in GPL6104.annot.gz\n")
              return(NULL)
            }
            
            # Filter and deduplicate
            probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
            probe_map <- probe_map[!is.na(probe_map$ID) & probe_map$ID != "", ]
            probe_map <- probe_map[!duplicated(probe_map$ID), ]
            
            # Save mapping
            if (nrow(probe_map) > 100) {
              saveRDS(probe_map, mapping_output)
              cat(paste0("Saved GPL6104 mapping from annot.gz: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
              return(probe_map)
            } else {
              cat("⚠️  GPL6104 mapping has insufficient rows\n")
              return(NULL)
            }
          } else {
            cat("❌ No table data found in GPL6104.annot.gz\n")
            return(NULL)
          }
        } else {
          cat("❌ !platform_table_begin not found in GPL6104.annot.gz\n")
          return(NULL)
        }
      }, error = function(e) {
        cat(paste0("❌ Error parsing GPL6104.annot.gz: ", conditionMessage(e), "\n"))
        # Continue to other methods
      })
    }
  }
  
  # Offline parsing for GPL6244
  if (gpl_id == "GPL6244") {
    txt_file <- file.path(base_dir, "data", "raw", "GPL", "GPL6244-17930.txt")
    if (file.exists(txt_file)) {
      cat("Special handling for GPL6244: Using offline GPL6244-17930.txt\n")
      
      tryCatch({
        # Read the file, skipping comment lines
        lines <- readLines(txt_file)
        lines <- lines[!grepl("^#", lines)]
        lines <- lines[lines != ""]
        
        # Get header line
        header_line <- lines[1]
        header_cols <- strsplit(header_line, "\t")[[1]]
        
        # Find indices of ID and gene_assignment columns
        id_col_idx <- which(header_cols == "ID")
        gene_assignment_col_idx <- which(header_cols == "gene_assignment")
        
        if (length(id_col_idx) == 0 || length(gene_assignment_col_idx) == 0) {
          cat("❌ Required columns not found in GPL6244-17930.txt\n")
          return(NULL)
        }
        
        # Function to extract gene symbol from gene_assignment
        extract_symbol <- function(assignment) {
          if (is.na(assignment) || assignment == "" || assignment == "---") {
            return(NA)
          }
          # Split by ///
          entries <- strsplit(assignment, "///", fixed = TRUE)[[1]]
          for (entry in entries) {
            # Split by // 
            parts <- strsplit(entry, " // ", fixed = TRUE)[[1]]
            if (length(parts) >= 2) {
              symbol <- trimws(parts[2])
              if (symbol != "" && symbol != "---") {
                return(symbol)
              }
            }
          }
          return(NA)
        }
        
        # Process each line individually to handle inconsistent spacing
        ids <- character()
        gene_symbols <- character()
        
        for (line in lines[-1]) {  # Skip header
          # Split by tabs
          parts <- strsplit(line, "\t")[[1]]
          
          # Ensure we have enough columns
          if (length(parts) >= max(id_col_idx, gene_assignment_col_idx)) {
            id <- parts[id_col_idx]
            gene_assignment <- parts[gene_assignment_col_idx]
            
            # Extract symbol
            symbol <- extract_symbol(gene_assignment)
            
            ids <- c(ids, id)
            gene_symbols <- c(gene_symbols, symbol)
          }
        }
        
        # Create mapping
        probe_map <- data.frame(
          ID = as.character(ids),
          Gene_Symbol = as.character(gene_symbols),
          stringsAsFactors = FALSE
        )
        
        # Filter and deduplicate
        probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
        probe_map <- probe_map[!is.na(probe_map$ID) & probe_map$ID != "", ]
        probe_map <- probe_map[!duplicated(probe_map$ID), ]
        
        # Save mapping
        if (nrow(probe_map) > 100) {
          saveRDS(probe_map, mapping_output)
          cat(paste0("Saved GPL6244 mapping from txt: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
          return(probe_map)
        } else {
          cat("⚠️  GPL6244 mapping has insufficient rows\n")
          return(NULL)
        }
      }, error = function(e) {
        cat(paste0("❌ Error parsing GPL6244-17930.txt: ", conditionMessage(e), "\n"))
        # Continue to other methods
      })
    }
  }
  
  # Special handling for other GPL platforms
  if (gpl_id == "GPL23961") {
    cat("Special handling for GPL23961: Using direct GEO query and extracting only SPOT_ID and ORF columns\n")
    suppressPackageStartupMessages(library(GEOquery))
    
    tryCatch({
      # Get GPL23961 directly from GEO (smaller than local family.soft.gz)
      gpl_obj <- getGEO("GPL23961")
      tbl <- Table(gpl_obj)
      
      # Extract only the needed columns to minimize memory usage
      probe_map <- data.frame(
        spot_id = as.character(tbl$SPOT_ID),
        orf = as.character(tbl$ORF),
        stringsAsFactors = FALSE
      )
      
      # Clean up large objects immediately
      rm(gpl_obj, tbl)
      gc()
      
      # Filter and deduplicate
      probe_map <- probe_map[!is.na(probe_map$orf) & probe_map$orf != "" & probe_map$orf != "---", ]
      probe_map <- probe_map[!duplicated(probe_map$spot_id), ]
      
      # Save mapping
      if (nrow(probe_map) > 100) {
        saveRDS(probe_map, mapping_output)
        cat(paste0("Saved GPL23961 mapping: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
        return(probe_map)
      } else {
        cat("⚠️  GPL23961 mapping has insufficient rows\n")
        return(NULL)
      }
    }, error = function(e) {
      cat(paste0("❌ Error building GPL23961 mapping: ", conditionMessage(e), "\n"))
      return(NULL)
    })
  } else if (gpl_id == "GPL13158") {
    # First check for offline annotation file (GPL13158-5065.txt)
    txt_file <- file.path(base_dir, "data", "raw", "GPL", "GPL13158-5065.txt")
    if (file.exists(txt_file)) {
      cat("Special handling for GPL13158: Using offline GPL13158-5065.txt\n")
      
      tryCatch({
        # Read the file, skipping comment lines
        lines <- readLines(txt_file)
        lines <- lines[!grepl("^#", lines)]
        lines <- lines[lines != ""]
        
        # Get header line
        header_line <- lines[1]
        header_cols <- strsplit(header_line, "\t")[[1]]
        
        # Find indices of ID and Gene Symbol columns
        id_col_idx <- which(header_cols == "ID")
        gene_symbol_col_idx <- which(header_cols == "Gene Symbol")
        
        if (length(id_col_idx) == 0 || length(gene_symbol_col_idx) == 0) {
          cat("❌ Required columns not found in GPL13158-5065.txt\n")
          return(NULL)
        }
        
        # Process each line individually to handle inconsistent spacing
        ids <- character()
        gene_symbols <- character()
        
        for (line in lines[-1]) {  # Skip header
          # Split by tabs
          parts <- strsplit(line, "\t")[[1]]
          
          # Ensure we have enough columns
          if (length(parts) >= max(id_col_idx, gene_symbol_col_idx)) {
            id <- parts[id_col_idx]
            gene_symbol <- parts[gene_symbol_col_idx]
            
            ids <- c(ids, id)
            gene_symbols <- c(gene_symbols, gene_symbol)
          }
        }
        
        # Create mapping
        probe_map <- data.frame(
          ID = as.character(ids),
          Gene_Symbol = as.character(gene_symbols),
          stringsAsFactors = FALSE
        )
        
        # Filter and deduplicate
        probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
        probe_map <- probe_map[!is.na(probe_map$ID) & probe_map$ID != "", ]
        probe_map <- probe_map[!duplicated(probe_map$ID), ]
        
        # Save mapping
        if (nrow(probe_map) > 100) {
          saveRDS(probe_map, mapping_output)
          cat(paste0("Saved GPL13158 mapping from txt: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
          return(probe_map)
        } else {
          cat("⚠️  GPL13158 mapping has insufficient rows\n")
          return(NULL)
        }
      }, error = function(e) {
        cat(paste0("❌ Error parsing GPL13158-5065.txt: ", conditionMessage(e), "\n"))
        # Continue to other methods
      })
    }
    
    # Fallback to annot.gz file
    annot_file <- file.path(base_dir, "data", "raw", "GPL", "GPL13158.annot.gz")
    if (file.exists(annot_file)) {
      cat("Special handling for GPL13158: Using offline GPL13158.annot.gz\n")
      
      tryCatch({
        # Read the file and find the platform table
        con <- gzfile(annot_file, "r")
        lines <- readLines(con)
        close(con)
        
        # Remove empty lines
        lines <- lines[lines != ""]
        
        # Find the platform table begin
        table_start <- which(grepl("!platform_table_begin", lines))[1]
        if (!is.na(table_start)) {
          # Find the header line (first line after table start that doesn't start with !)
          header_idx <- table_start + 1
          while (header_idx <= length(lines) && grepl("^!", lines[header_idx])) {
            header_idx <- header_idx + 1
          }
          
          if (header_idx <= length(lines)) {
            # Get table data
            table_data <- lines[header_idx:length(lines)]
            
            # Write to a temporary file to handle line endings properly
            temp_file <- tempfile(fileext = ".txt")
            writeLines(table_data, temp_file)
            
            # Read the table with fill=TRUE to handle inconsistent line lengths
            tbl <- read.table(temp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
            unlink(temp_file)
            
            # Use ID and Gene symbol columns
            if ("ID" %in% colnames(tbl) && "Gene.symbol" %in% colnames(tbl)) {
              # Use Gene.symbol column (note the dot instead of space)
              probe_map <- data.frame(
                ID = as.character(tbl$ID),
                Gene_Symbol = as.character(tbl$Gene.symbol),
                stringsAsFactors = FALSE
              )
            } else if ("ID" %in% colnames(tbl) && "Gene symbol" %in% colnames(tbl)) {
              # Use Gene symbol column (with space)
              probe_map <- data.frame(
                ID = as.character(tbl$ID),
                Gene_Symbol = as.character(tbl$`Gene symbol`),
                stringsAsFactors = FALSE
              )
            } else {
              cat("❌ Required columns not found in GPL13158.annot.gz\n")
              return(NULL)
            }
            
            # Filter and deduplicate
            probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
            probe_map <- probe_map[!is.na(probe_map$ID) & probe_map$ID != "", ]
            probe_map <- probe_map[!duplicated(probe_map$ID), ]
            
            # Save mapping
            if (nrow(probe_map) > 100) {
              saveRDS(probe_map, mapping_output)
              cat(paste0("Saved GPL13158 mapping from annot.gz: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
              return(probe_map)
            } else {
              cat("⚠️  GPL13158 mapping has insufficient rows\n")
              return(NULL)
            }
          } else {
            cat("❌ No table data found in GPL13158.annot.gz\n")
            return(NULL)
          }
        } else {
          cat("❌ !platform_table_begin not found in GPL13158.annot.gz\n")
          return(NULL)
        }
      }, error = function(e) {
        cat(paste0("❌ Error parsing GPL13158.annot.gz: ", conditionMessage(e), "\n"))
        # Continue to other methods
      })
    }
  }
  
  # Special handling for GPL6244
  if (gpl_id == "GPL6244") {
    cat("Special handling for GPL6244: Using direct GEO query and printing column names\n")
    suppressPackageStartupMessages(library(GEOquery))
    
    tryCatch({
      # Get GPL6244 directly from GEO
      gpl_obj <- getGEO("GPL6244")
      tbl <- Table(gpl_obj)
      
      # Print column names for debugging
      cat("GPL6244 column names: ", paste(colnames(tbl), collapse = ", "), "\n")
      
      # Try to find appropriate columns for ID and gene symbol
      colnames_lower <- tolower(colnames(tbl))
      
      # Find key column with extended candidates
      key_candidates <- c("id", "id_ref", "probe", "ilmnid", "illuminaid", "probesetid", "probe_id", "probeset")
      key_col_idx <- which(colnames_lower %in% key_candidates)[1]
      
      # Find symbol column with extended candidates
      symbol_candidates <- c("gene_symbol", "symbol", "genesymbol", "orf", "gene", "gb_acc", "gene symbol")
      symbol_col_idx <- which(colnames_lower %in% symbol_candidates)[1]
      
      if (!is.na(key_col_idx) && !is.na(symbol_col_idx)) {
        key_col_name <- colnames(tbl)[key_col_idx]
        symbol_col_name <- colnames(tbl)[symbol_col_idx]
        
        cat(paste0("Using columns for GPL6244: ", key_col_name, " (key) and ", symbol_col_name, " (symbol)\n"))
        
        # Extract only the needed columns to minimize memory usage
        probe_map <- data.frame(
          ID = as.character(tbl[[key_col_name]]),
          Gene_Symbol = as.character(tbl[[symbol_col_name]]),
          stringsAsFactors = FALSE
        )
        
        # Clean up large objects immediately
        rm(gpl_obj, tbl)
        gc()
        
        # Filter and deduplicate using base R to reduce memory usage
        probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
        probe_map <- probe_map[!duplicated(probe_map$ID), ]
        
        # Save mapping
        if (nrow(probe_map) > 100) {
          saveRDS(probe_map, mapping_output)
          cat(paste0("Saved GPL6244 mapping: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
          return(probe_map)
        } else {
          cat("⚠️  GPL6244 mapping has insufficient rows\n")
          return(NULL)
        }
      } else {
        cat("❌ Could not find suitable columns for GPL6244\n")
        return(NULL)
      }
    }, error = function(e) {
      cat(paste0("❌ Error building GPL6244 mapping: ", conditionMessage(e), "\n"))
      return(NULL)
    })
  } else if (gpl_id == "GPL6104") {
    cat("Special handling for GPL6104: Using direct GEO query with caching\n")
    suppressPackageStartupMessages(library(GEOquery))
    
    tryCatch({
      # Get GPL6104 directly from GEO
      gpl_obj <- getGEO("GPL6104")
      tbl <- Table(gpl_obj)
      
      # Try to find appropriate columns for ID and gene symbol
      colnames_lower <- tolower(colnames(tbl))
      
      # Find key column
      key_candidates <- c("id", "id_ref", "probe", "ilmnid", "illuminaid", "probesetid")
      key_col_idx <- which(colnames_lower %in% key_candidates)[1]
      
      # Find symbol column
      symbol_candidates <- c("gene_symbol", "symbol", "genesymbol", "orf", "gene", "gb_acc")
      symbol_col_idx <- which(colnames_lower %in% symbol_candidates)[1]
      
      if (!is.na(key_col_idx) && !is.na(symbol_col_idx)) {
        key_col_name <- colnames(tbl)[key_col_idx]
        symbol_col_name <- colnames(tbl)[symbol_col_idx]
        
        cat(paste0("Using columns for GPL6104: ", key_col_name, " (key) and ", symbol_col_name, " (symbol)\n"))
        
        # Extract only the needed columns to minimize memory usage
        probe_map <- data.frame(
          ID = as.character(tbl[[key_col_name]]),
          Gene_Symbol = as.character(tbl[[symbol_col_name]]),
          stringsAsFactors = FALSE
        )
        
        # Clean up large objects immediately
        rm(gpl_obj, tbl)
        gc()
        
        # Filter and deduplicate using base R to reduce memory usage
        probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
        probe_map <- probe_map[!duplicated(probe_map$ID), ]
        
        # Save mapping for offline reuse
        if (nrow(probe_map) > 100) {
          saveRDS(probe_map, mapping_output)
          cat(paste0("Saved GPL6104 mapping: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
          return(probe_map)
        } else {
          cat("⚠️  GPL6104 mapping has insufficient rows\n")
          return(NULL)
        }
      } else {
        cat("❌ Could not find suitable columns for GPL6104\n")
        return(NULL)
      }
    }, error = function(e) {
      cat(paste0("❌ Error building GPL6104 mapping: ", conditionMessage(e), "\n"))
      cat("⚠️  Marking GPL6104 as FAILED and skipping\n")
      return(NULL)
    })
  } else if (gpl_id == "GPL14550") {
    # First check for the new text file
    txt_file <- file.path(base_dir, "data", "raw", "GPL", "GPL14550-9757.txt")
    if (file.exists(txt_file)) {
      cat("Special handling for GPL14550: Using offline GPL14550-9757.txt\n")
      
      tryCatch({
        # Read the file, skipping comment lines
        lines <- readLines(txt_file)
        lines <- lines[!grepl("^#", lines)]
        lines <- lines[lines != ""]
        
        # Get header line
        header_line <- lines[1]
        header_cols <- strsplit(header_line, "\t")[[1]]
        
        # Find indices of ID, GENE_SYMBOL, GENE, and REFSEQ columns
        id_col_idx <- which(header_cols == "ID")
        gene_symbol_col_idx <- which(header_cols == "GENE_SYMBOL")
        gene_col_idx <- which(header_cols == "GENE")
        refseq_col_idx <- which(header_cols == "REFSEQ")
        
        if (length(id_col_idx) == 0) {
          cat("❌ Required ID column not found in GPL14550-9757.txt\n")
          return(NULL)
        }
        
        # Process each line individually to handle inconsistent spacing
        ids <- character()
        gene_symbols <- character()
        
        for (line in lines[-1]) {  # Skip header
          # Split by tabs
          parts <- strsplit(line, "\t")[[1]]
          
          # Ensure we have enough columns for ID
          if (length(parts) >= id_col_idx) {
            id <- parts[id_col_idx]
            gene_symbol <- NA
            
            # Try to get gene symbol from GENE_SYMBOL column first
            if (length(gene_symbol_col_idx) > 0 && length(parts) >= gene_symbol_col_idx) {
              gene_symbol <- parts[gene_symbol_col_idx]
            }
            
            # If no gene symbol, try to get from GENE column (Entrez ID)
            if (is.na(gene_symbol) || gene_symbol == "" || gene_symbol == "---") {
              if (length(gene_col_idx) > 0 && length(parts) >= gene_col_idx) {
                gene_symbol <- parts[gene_col_idx]
              }
            }
            
            # If still no gene symbol, try to get from REFSEQ column
            if (is.na(gene_symbol) || gene_symbol == "" || gene_symbol == "---") {
              if (length(refseq_col_idx) > 0 && length(parts) >= refseq_col_idx) {
                gene_symbol <- parts[refseq_col_idx]
              }
            }
            
            ids <- c(ids, id)
            gene_symbols <- c(gene_symbols, gene_symbol)
          }
        }
        
        # Create mapping
        probe_map <- data.frame(
          ID = as.character(ids),
          Gene_Symbol = as.character(gene_symbols),
          stringsAsFactors = FALSE
        )
        
        # Filter and deduplicate
        probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
        probe_map <- probe_map[!is.na(probe_map$ID) & probe_map$ID != "", ]
        probe_map <- probe_map[!duplicated(probe_map$ID), ]
        
        # Convert Entrez IDs to HGNC symbols if needed
        suppressPackageStartupMessages(library(AnnotationDbi))
        suppressPackageStartupMessages(library(org.Hs.eg.db))
        
        # Identify Entrez IDs (numeric values)
        entrez_ids <- probe_map$Gene_Symbol[grepl("^[0-9]+", probe_map$Gene_Symbol)]
        if (length(entrez_ids) > 0) {
          cat("Converting Entrez IDs to HGNC symbols for GPL14550\n")
          # Convert Entrez IDs to symbols
          symbol_map <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
          # Update probe_map with symbols
          probe_map$Gene_Symbol[grepl("^[0-9]+", probe_map$Gene_Symbol)] <- as.character(symbol_map[entrez_ids])
          # Remove rows with NA symbols
          probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "", ]
        }
        
        # Also handle RefSeq IDs if present
        refseq_ids <- probe_map$Gene_Symbol[grepl("^NM_|^NR_|^XR_", probe_map$Gene_Symbol)]
        if (length(refseq_ids) > 0) {
          cat("Converting RefSeq IDs to HGNC symbols for GPL14550\n")
          # Convert RefSeq IDs to symbols
          symbol_map <- mapIds(org.Hs.eg.db, keys = refseq_ids, column = "SYMBOL", keytype = "REFSEQ", multiVals = "first")
          # Update probe_map with symbols
          probe_map$Gene_Symbol[grepl("^NM_|^NR_|^XR_", probe_map$Gene_Symbol)] <- as.character(symbol_map[refseq_ids])
          # Remove rows with NA symbols
          probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "", ]
        }
        
        # Remove rows with non-standard gene symbols (like XLOC_*)
        probe_map <- probe_map[!grepl("^XLOC_", probe_map$Gene_Symbol), ]
        
        # Save mapping
        if (nrow(probe_map) > 100) {
          saveRDS(probe_map, mapping_output)
          cat(paste0("Saved GPL14550 mapping from txt: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
          return(probe_map)
        } else {
          cat("⚠️  GPL14550 mapping has insufficient rows\n")
          return(NULL)
        }
      }, error = function(e) {
        cat(paste0("❌ Error parsing GPL14550-9757.txt: ", conditionMessage(e), "\n"))
        # Continue to other methods
      })
    }
    
    # Fallback to family.soft.gz file
    family_file <- file.path(base_dir, "data", "raw", "GPL", "GPL14550_family.soft.gz")
    if (file.exists(family_file)) {
      cat("Special handling for GPL14550: Using offline family.soft.gz file\n")
      
      tryCatch({
        # Read the file and find the platform table
        con <- gzfile(family_file, "r")
        lines <- readLines(con)
        close(con)
        
        # Remove empty lines
        lines <- lines[lines != ""]
        
        # Find the platform table begin
        table_start <- which(grepl("!platform_table_begin", lines))[1]
        if (!is.na(table_start)) {
          # Find the header line (first line after table start that doesn't start with !)
          header_idx <- table_start + 1
          while (header_idx <= length(lines) && grepl("^!", lines[header_idx])) {
            header_idx <- header_idx + 1
          }
          
          if (header_idx <= length(lines)) {
            # Get table data
            table_data <- lines[header_idx:length(lines)]
            
            # Write to a temporary file to handle line endings properly
            temp_file <- tempfile(fileext = ".txt")
            writeLines(table_data, temp_file)
            
            # Read the table with fill=TRUE to handle inconsistent line lengths
            tbl <- read.table(temp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
            unlink(temp_file)
            
            # Use ID and Gene symbol columns
            if ("ID" %in% colnames(tbl) && "Gene.symbol" %in% colnames(tbl)) {
              # Use Gene.symbol column (note the dot instead of space)
              probe_map <- data.frame(
                ID = as.character(tbl$ID),
                Gene_Symbol = as.character(tbl$Gene.symbol),
                stringsAsFactors = FALSE
              )
            } else if ("ID" %in% colnames(tbl) && "Gene symbol" %in% colnames(tbl)) {
              # Use Gene symbol column (with space)
              probe_map <- data.frame(
                ID = as.character(tbl$ID),
                Gene_Symbol = as.character(tbl$`Gene symbol`),
                stringsAsFactors = FALSE
              )
            } else {
              cat("❌ Required columns not found in GPL14550_family.soft.gz\n")
              return(NULL)
            }
            
            # Filter and deduplicate
            probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
            probe_map <- probe_map[!is.na(probe_map$ID) & probe_map$ID != "", ]
            probe_map <- probe_map[!duplicated(probe_map$ID), ]
            
            # Save mapping
            if (nrow(probe_map) > 100) {
              saveRDS(probe_map, mapping_output)
              cat(paste0("Saved GPL14550 mapping from family.soft.gz: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
              return(probe_map)
            } else {
              cat("⚠️  GPL14550 mapping has insufficient rows\n")
              return(NULL)
            }
          } else {
            cat("❌ No table data found in GPL14550_family.soft.gz\n")
            return(NULL)
          }
        } else {
          cat("❌ !platform_table_begin not found in GPL14550_family.soft.gz\n")
          return(NULL)
        }
      }, error = function(e) {
        cat(paste0("❌ Error parsing GPL14550_family.soft.gz: ", conditionMessage(e), "\n"))
        return(NULL)
      })
    }
    
    cat("❌ No suitable GPL14550 files found\n")
    return(NULL)
  }
  
  # Try to find GPL files
  gpl_dir <- file.path(base_dir, "data", "raw", "GPL")
  gpl_files <- list.files(gpl_dir, recursive = TRUE, full.names = TRUE)
  gpl_files <- gpl_files[grepl(gpl_id, basename(gpl_files))]
  
  if (length(gpl_files) == 0) {
    cat(paste0("❌ No GPL files found for ", gpl_id, "\n"))
    return(NULL)
  }
  
  # For GPL14550, use family.soft.gz file
  if (gpl_id == "GPL14550") {
    family_file <- gpl_files[grepl("family\\.soft\\.gz$", gpl_files)]
    if (length(family_file) > 0) {
      cat("Special handling for GPL14550: Using offline family.soft.gz file\n")
      gpl_files <- family_file
    } else {
      cat("❌ No GPL14550 family.soft.gz file found\n")
      return(NULL)
    }
  } else {
    # Filter out family.soft.gz files, .bgx.gz files, and large files (>1GB) to avoid memory issues
    gpl_files <- gpl_files[!grepl("family\\.soft\\.gz$", gpl_files)]
    gpl_files <- gpl_files[!grepl("\\.bgx\\.gz$", gpl_files)]  # Skip bgx.gz files
    gpl_files <- gpl_files[file.size(gpl_files) < 1073741824]  # < 1GB
  }
  
  if (length(gpl_files) == 0) {
    cat(paste0("❌ No suitable GPL files found for ", gpl_id, " (all filtered out)\n"))
    # Try online getGEO as fallback
    cat(paste0("Trying online getGEO for ", gpl_id, "\n"))
    suppressPackageStartupMessages(library(GEOquery))
    
    tryCatch({
      gpl_obj <- getGEO(gpl_id)
      tbl <- Table(gpl_obj)
      
      # Try to find appropriate columns for ID and gene symbol
      colnames_lower <- tolower(colnames(tbl))
      
      # Find key column
      key_candidates <- c("id", "id_ref", "probe", "ilmnid", "illuminaid", "probesetid")
      key_col_idx <- which(colnames_lower %in% key_candidates)[1]
      
      # Find symbol column
      symbol_candidates <- c("gene_symbol", "symbol", "genesymbol", "orf", "gene", "gb_acc")
      symbol_col_idx <- which(colnames_lower %in% symbol_candidates)[1]
      
      if (!is.na(key_col_idx) && !is.na(symbol_col_idx)) {
        key_col_name <- colnames(tbl)[key_col_idx]
        symbol_col_name <- colnames(tbl)[symbol_col_idx]
        
        cat(paste0("Using columns for ", gpl_id, ": ", key_col_name, " (key) and ", symbol_col_name, " (symbol)\n"))
        
        # Extract only the needed columns to minimize memory usage
        probe_map <- data.frame(
          ID = as.character(tbl[[key_col_name]]),
          Gene_Symbol = as.character(tbl[[symbol_col_name]]),
          stringsAsFactors = FALSE
        )
        
        # Clean up large objects immediately
        rm(gpl_obj, tbl)
        gc()
        
        # Filter and deduplicate using base R to reduce memory usage
        probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
        probe_map <- probe_map[!duplicated(probe_map$ID), ]
        
        # Save mapping
        if (nrow(probe_map) > 100) {
          saveRDS(probe_map, mapping_output)
          cat(paste0("Saved ", gpl_id, " mapping: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
          return(probe_map)
        } else {
          cat(paste0("⚠️  ", gpl_id, " mapping has insufficient rows\n"))
          return(NULL)
        }
      } else {
        cat(paste0("❌ Could not find suitable columns for ", gpl_id, "\n"))
        return(NULL)
      }
    }, error = function(e) {
      cat(paste0("❌ Error building ", gpl_id, " mapping: ", conditionMessage(e), "\n"))
      return(NULL)
    })
  }
  
  # Try to read GPL file using GEOquery
  suppressPackageStartupMessages(library(GEOquery))
  
  # Prioritize .soft.gz files over other file types
  soft_files <- gpl_files[grepl("\\.soft\\.gz$", gpl_files)]
  if (length(soft_files) > 0) {
    # Pick the smallest soft file to avoid memory issues
    soft_files <- soft_files[order(file.size(soft_files))]
    gpl_file <- soft_files[1]
  } else {
    # If no soft files, pick the smallest file
    gpl_files <- gpl_files[order(file.size(gpl_files))]
    gpl_file <- gpl_files[1]
  }
  
  # Fix bug: ensure gpl_file is not NA
  if (is.na(gpl_file)) {
    cat(paste0("❌ No valid GPL file found for ", gpl_id, "\n"))
    return(NULL)
  }
  
  cat(paste0("Building mapping for ", gpl_id, " using file: ", basename(gpl_file), "\n"))
  
  tryCatch({
    gpl_obj <- getGEO(filename = gpl_file)
    gpl_tbl <- Table(gpl_obj)
    
    if (is.null(gpl_tbl) || nrow(gpl_tbl) == 0) {
      cat(paste0("❌ Failed to read GPL table for ", gpl_id, "\n"))
      return(NULL)
    }
    
    # Convert column names to lowercase
    colnames(gpl_tbl) <- tolower(colnames(gpl_tbl))
    
    # Try to find appropriate columns for ID and gene symbol
    colnames_lower <- colnames(gpl_tbl)
    
    # First try the provided columns (case-insensitive)
    key_col_lower <- tolower(key_col)
    symbol_col_lower <- tolower(symbol_col)
    
    # If provided columns not found, try candidate columns
    if (!(key_col_lower %in% colnames_lower) || !(symbol_col_lower %in% colnames_lower)) {
      cat(paste0("⚠️  Provided columns not found: ", key_col, " and ", symbol_col, "\n"))
      cat("Trying candidate columns...\n")
      
      # Find key column from candidates
      key_candidates <- c("id", "id_ref", "probe", "ilmnid", "illuminaid", "probesetid")
      key_col_idx <- which(colnames_lower %in% key_candidates)[1]
      
      # Find symbol column from candidates
      symbol_candidates <- c("gene_symbol", "symbol", "genesymbol", "orf", "gene", "gb_acc")
      symbol_col_idx <- which(colnames_lower %in% symbol_candidates)[1]
      
      if (!is.na(key_col_idx) && !is.na(symbol_col_idx)) {
        key_col_lower <- colnames_lower[key_col_idx]
        symbol_col_lower <- colnames_lower[symbol_col_idx]
        cat(paste0("Using candidate columns: ", key_col_lower, " (key) and ", symbol_col_lower, " (symbol)\n"))
      } else {
        cat("❌ Could not find suitable columns in GPL table\n")
        return(NULL)
      }
    }
    
    # Build mapping
    probe_map <- data.frame(
      ID = as.character(gpl_tbl[[key_col_lower]]),
      Gene_Symbol = as.character(gpl_tbl[[symbol_col_lower]]),
      stringsAsFactors = FALSE
    )
    
    # Filter and deduplicate using base R to reduce memory usage
    probe_map <- probe_map[!is.na(probe_map$Gene_Symbol) & probe_map$Gene_Symbol != "" & probe_map$Gene_Symbol != "---", ]
    probe_map <- probe_map[!duplicated(probe_map$ID), ]
    
    # Write threshold check
    min_mapping_rows <- 100
    if (nrow(probe_map) < min_mapping_rows) {
      cat(paste0("⛔ Mapping has insufficient rows (", nrow(probe_map), "), not saving\n"))
      return(NULL)
    }
    
    # Save mapping
    saveRDS(probe_map, mapping_output)
    cat(paste0("Saved mapping file: ", basename(mapping_output), " | n=", nrow(probe_map), "\n"))
    
    return(probe_map)
  }, error = function(e) {
    cat(paste0("❌ Error building mapping for ", gpl_id, ": ", conditionMessage(e), "\n"))
    return(NULL)
  }
)
}

# Function to map IDs to symbols using org.Hs.eg.db with batch processing and deduplication
map_ids_to_symbols <- function(ids, id_type, batch_size = 10000) {
  ids <- as.character(ids)
  ids_unique <- unique(ids)
  
  keytype <- if (id_type == "ENSG") "ENSEMBL" else if (id_type == "ENTREZ") "ENTREZID" else return(rep(NA_character_, length(ids)))
  out_unique <- rep(NA_character_, length(ids_unique))
  names(out_unique) <- ids_unique
  
  idx <- split(seq_along(ids_unique), ceiling(seq_along(ids_unique)/batch_size))
  for (ii in idx) {
    keys <- ids_unique[ii]
    
    # Filter keys to avoid invalid keys causing hard stop
    keys <- unique(keys)
    keys <- keys[!is.na(keys) & keys != ""]
    
    if (id_type == "ENSG") {
      keys <- keys[grepl("^ENSG\\d+$", keys)]
    } else if (id_type == "ENTREZ") {
      keys <- keys[grepl("^[0-9]+$", keys)]
    }
    
    if (length(keys) == 0) {
      # No valid keys, return NA for this batch
      out_unique[ii] <- NA_character_
    } else {
      m <- suppressMessages(
        mapIds(org.Hs.eg.db, keys = keys, column = "SYMBOL", keytype = keytype, multiVals = "first")
      )
      out_unique[ii] <- as.character(m)
      rm(m); gc()
    }
  }
  out_unique[ids]
}

# Function to process HTS counts (GPL11154, GPL16791)
process_hts_counts <- function(expr, id_type, gpl, gse, cohort) {
  # Standardize IDs based on type
  feature_id_raw <- rownames(expr)
  feature_id_clean <- feature_id_raw
  
  if (id_type == "ENSG") {
    # Clean ENSG IDs (remove version)
    feature_id_clean <- clean_ensg(feature_id_clean)
  } else if (id_type == "ENTREZ") {
    # Clean Entrez IDs
    feature_id_clean <- clean_entrez(feature_id_clean)
  }
  
  # Map IDs to symbols
  gene_symbols <- map_ids_to_symbols(feature_id_clean, id_type)
  
  # Create ID mapping table
  id_map <- create_id_map(
    feature_id_raw = feature_id_raw,
    feature_id_clean = feature_id_clean,
    gene_symbol = gene_symbols,
    mapping_source = "HTS_ANNOT",
    gpl = gpl,
    gse = gse,
    cohort = cohort
  )
  
  # Remove features with no gene mapping
  valid_features <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  
  if (sum(valid_features) == 0) {
    cat("⚠️  No valid features with gene mapping\n")
    # Save empty ID map
    id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
    write_tsv(id_map, id_map_output)
    cat("Saved ID mapping table: ", basename(id_map_output), "\n")
    return(list(expr = NULL, id_map = id_map))
  }
  
  # Subset expression matrix and ID map
  expr <- expr[valid_features, ]
  gene_symbols <- gene_symbols[valid_features]
  id_map <- id_map[valid_features, ]
  
  # Set row names to gene symbols
  rownames(expr) <- gene_symbols
  
  # Deduplicate by gene symbol
  gene_symbols <- as.character(gene_symbols)  # Ensure character type to avoid copying
  original_probe_count <- length(gene_symbols)
  unique_gene_count <- length(unique(gene_symbols))
  duplicate_count <- original_probe_count - unique_gene_count
  
  if (duplicate_count > 0) {
    cat(paste0("Original probe count: ", original_probe_count, "\n"))
    cat(paste0("Mapped unique GeneSymbol count: ", unique_gene_count, "\n"))
    cat(paste0("Duplicate genes aggregated: ", duplicate_count, "\n"))
    
    # Use max of mean expression as aggregation strategy
    cat("Aggregating duplicate gene symbols with max(mean expression)\n")
    
    # Calculate mean expression for each probe
    probe_means <- rowMeans(expr)
    
    # Create a data frame with gene symbols, mean expression, and original indices
    probe_df <- data.frame(
      gene_symbol = gene_symbols,
      mean_expr = probe_means,
      original_idx = 1:length(gene_symbols),
      stringsAsFactors = FALSE
    )
    
    # For each gene symbol, keep the probe with the highest mean expression
    # Using base R for compatibility with Bioconductor packages
    max_mean_idx <- integer(0)
    unique_genes <- unique(gene_symbols)
    
    for (gene in unique_genes) {
      gene_indices <- which(gene_symbols == gene)
      gene_means <- probe_means[gene_indices]
      max_idx <- gene_indices[which.max(gene_means)][1]  # Keep first occurrence in case of ties
      max_mean_idx <- c(max_mean_idx, max_idx)
    }
    
    # Subset expression matrix and ID map
    expr <- expr[max_mean_idx, ]
    id_map <- id_map[max_mean_idx, ]
    rownames(expr) <- gene_symbols[max_mean_idx]
    
    # Double-check for any remaining duplicates
    if (any(duplicated(rownames(expr)))) {
      cat("⚠️  Found remaining duplicates after aggregation, removing them\n")
      unique_idx <- !duplicated(rownames(expr))
      expr <- expr[unique_idx, ]
      id_map <- id_map[unique_idx, ]
    }
  }
  
  return(list(expr = expr, id_map = id_map))
}

# Function to process GPL23961 with SPOT_ID mapping
process_gpl23961 <- function(expr, mapping, gpl, gse, cohort) {
  # Convert probe IDs to SPOT_IDs by removing _at suffix only
  probe_ids <- sub("_at$", "", rownames(expr))  # Remove _at suffix
  
  # Create probe to gene mapping - FORCE SPOT_ID and ORF columns
  if ("spot_id" %in% tolower(colnames(mapping)) && "orf" %in% tolower(colnames(mapping))) {
    # Use SPOT_ID and ORF columns
    spot_id_col <- which(tolower(colnames(mapping)) == "spot_id")
    symbol_col <- which(tolower(colnames(mapping)) == "orf")
    cat("Using SPOT_ID and ORF columns for mapping\n")
  } else {
    cat("❌ Required columns (SPOT_ID, ORF) not found in mapping\n")
    # Create empty ID map
    id_map <- create_id_map(
      feature_id_raw = rownames(expr),
      feature_id_clean = probe_ids,
      gene_symbol = rep(NA, length(probe_ids)),
      mapping_source = "GPL_TABLE_SPOTID",
      gpl = gpl,
      gse = gse,
      cohort = cohort
    )
    # Save empty ID map
    id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
    write_tsv(id_map, id_map_output)
    cat("Saved ID mapping table: ", basename(id_map_output), "\n")
    return(list(expr = NULL, id_map = id_map))
  }
  
  # Create mapping
  probe_to_gene <- setNames(as.character(mapping[[symbol_col]]), as.character(mapping[[spot_id_col]]))
  
  # Map probes to genes
  gene_symbols <- probe_to_gene[probe_ids]
  
  # Create ID mapping table
  id_map <- create_id_map(
    feature_id_raw = rownames(expr),
    feature_id_clean = probe_ids,
    gene_symbol = gene_symbols,
    mapping_source = "GPL_TABLE_SPOTID",
    gpl = gpl,
    gse = gse,
    cohort = cohort
  )
  
  # Remove probes with no gene mapping
  valid_probes <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  
  if (sum(valid_probes) == 0) {
    cat("⚠️  No valid probes with gene mapping\n")
    # Save ID map
    id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
    write_tsv(id_map, id_map_output)
    cat("Saved ID mapping table: ", basename(id_map_output), "\n")
    return(list(expr = NULL, id_map = id_map))
  }
  
  # Subset expression matrix and ID map
  expr <- expr[valid_probes, ]
  gene_symbols <- gene_symbols[valid_probes]
  id_map <- id_map[valid_probes, ]
  
  # Set row names to gene symbols
  rownames(expr) <- gene_symbols
  
  # Deduplicate by gene symbol
  gene_symbols <- as.character(gene_symbols)  # Ensure character type to avoid copying
  original_probe_count <- length(gene_symbols)
  unique_gene_count <- length(unique(gene_symbols))
  duplicate_count <- original_probe_count - unique_gene_count
  
  if (duplicate_count > 0) {
    cat(paste0("Original probe count: ", original_probe_count, "\n"))
    cat(paste0("Mapped unique GeneSymbol count: ", unique_gene_count, "\n"))
    cat(paste0("Duplicate genes aggregated: ", duplicate_count, "\n"))
    
    # Use max of mean expression as aggregation strategy
    cat("Aggregating duplicate gene symbols with max(mean expression)\n")
    
    # Calculate mean expression for each probe
    probe_means <- rowMeans(expr)
    
    # Create a data frame with gene symbols, mean expression, and original indices
    probe_df <- data.frame(
      gene_symbol = gene_symbols,
      mean_expr = probe_means,
      original_idx = 1:length(gene_symbols),
      stringsAsFactors = FALSE
    )
    
    # For each gene symbol, keep the probe with the highest mean expression
    # Using base R for compatibility with Bioconductor packages
    max_mean_idx <- integer(0)
    unique_genes <- unique(gene_symbols)
    
    for (gene in unique_genes) {
      gene_indices <- which(gene_symbols == gene)
      gene_means <- probe_means[gene_indices]
      max_idx <- gene_indices[which.max(gene_means)][1]  # Keep first occurrence in case of ties
      max_mean_idx <- c(max_mean_idx, max_idx)
    }
    
    # Subset expression matrix and ID map
    expr <- expr[max_mean_idx, ]
    id_map <- id_map[max_mean_idx, ]
    rownames(expr) <- gene_symbols[max_mean_idx]
    
    # Double-check for any remaining duplicates
    if (any(duplicated(rownames(expr)))) {
      cat("⚠️  Found remaining duplicates after aggregation, removing them\n")
      unique_idx <- !duplicated(rownames(expr))
      expr <- expr[unique_idx, ]
      id_map <- id_map[unique_idx, ]
    }
  }
  
  return(list(expr = expr, id_map = id_map))
}

# Function to create ID mapping table
create_id_map <- function(feature_id_raw, feature_id_clean, gene_symbol, mapping_source, gpl, gse, cohort) {
  data.frame(
    feature_id_raw = feature_id_raw,
    feature_id_clean = feature_id_clean,
    gene_symbol = gene_symbol,
    mapping_source = mapping_source,
    gpl = gpl,
    gse = gse,
    cohort = cohort,
    stringsAsFactors = FALSE
  )
}

# Function to process general probe platforms
process_general_probe <- function(expr, mapping, gpl, gse, cohort) {
  # Create probe to gene mapping
  probe_to_gene <- setNames(mapping$Gene_Symbol, mapping$ID)
  
  # Map probes to genes
  gene_symbols <- probe_to_gene[rownames(expr)]
  
  # Create ID mapping table
  id_map <- create_id_map(
    feature_id_raw = rownames(expr),
    feature_id_clean = rownames(expr),  # For probes, use as-is
    gene_symbol = gene_symbols,
    mapping_source = "GPL_TABLE",
    gpl = gpl,
    gse = gse,
    cohort = cohort
  )
  
  # Remove probes with no gene mapping
  valid_probes <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  
  if (sum(valid_probes) == 0) {
    cat("⚠️  No valid probes with gene mapping\n")
    # Save empty ID map
    id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
    write_tsv(id_map, id_map_output)
    cat("Saved ID mapping table: ", basename(id_map_output), "\n")
    return(list(expr = NULL, id_map = id_map))
  }
  
  # Subset expression matrix and ID map
  expr <- expr[valid_probes, ]
  gene_symbols <- gene_symbols[valid_probes]
  id_map <- id_map[valid_probes, ]
  
  # Set row names to gene symbols
  rownames(expr) <- gene_symbols
  
  # Deduplicate by gene symbol
  gene_symbols <- as.character(gene_symbols)  # Ensure character type to avoid copying
  original_probe_count <- length(gene_symbols)
  unique_gene_count <- length(unique(gene_symbols))
  duplicate_count <- original_probe_count - unique_gene_count
  
  if (duplicate_count > 0) {
    cat(paste0("Original probe count: ", original_probe_count, "\n"))
    cat(paste0("Mapped unique GeneSymbol count: ", unique_gene_count, "\n"))
    cat(paste0("Duplicate genes aggregated: ", duplicate_count, "\n"))
    
    # Use max of mean expression as aggregation strategy
    cat("Aggregating duplicate gene symbols with max(mean expression)\n")
    
    # Calculate mean expression for each probe
    probe_means <- rowMeans(expr)
    
    # Create a data frame with gene symbols, mean expression, and original indices
    probe_df <- data.frame(
      gene_symbol = gene_symbols,
      mean_expr = probe_means,
      original_idx = 1:length(gene_symbols),
      stringsAsFactors = FALSE
    )
    
    # For each gene symbol, keep the probe with the highest mean expression
    # Using base R for compatibility with Bioconductor packages
    max_mean_idx <- integer(0)
    unique_genes <- unique(gene_symbols)
    
    for (gene in unique_genes) {
      gene_indices <- which(gene_symbols == gene)
      gene_means <- probe_means[gene_indices]
      max_idx <- gene_indices[which.max(gene_means)][1]  # Keep first occurrence in case of ties
      max_mean_idx <- c(max_mean_idx, max_idx)
    }
    
    # Subset expression matrix and ID map
    expr <- expr[max_mean_idx, ]
    id_map <- id_map[max_mean_idx, ]
    rownames(expr) <- gene_symbols[max_mean_idx]
    
    # Double-check for any remaining duplicates
    if (any(duplicated(rownames(expr)))) {
      cat("⚠️  Found remaining duplicates after aggregation, removing them\n")
      unique_idx <- !duplicated(rownames(expr))
      expr <- expr[unique_idx, ]
      id_map <- id_map[unique_idx, ]
    }
  }
  
  return(list(expr = expr, id_map = id_map))
}

# Read matrix ID type table
matrix_id_type_file <- file.path(logs_dir, "00_matrix_id_type.tsv")
if (!file.exists(matrix_id_type_file)) {
  cat("❌ Matrix ID type table not found. Run 00a_ingest_expr.R first.\n")
  sink()
  stop("Matrix ID type table not found")
}

matrix_id_type_df <- read_tsv(matrix_id_type_file)

# Process each cohort
gpl_status_df <- data.frame()

for (i in 1:nrow(matrix_id_type_df)) {
  row <- matrix_id_type_df[i,]
  cohort <- row$cohort
  gpl <- row$gpl
  id_type <- row$id_type
  gse <- row$gse
  
  cat(paste0("\n=== Processing cohort: ", cohort, " (", gpl, ") ===\n"))
  
  # Read raw expression
  expr <- read_raw_expr(cohort)
  if (is.null(expr)) {
    next
  }
  
  # Check if there's a strategy for this GSE
  current_gse <- tolower(gse)
  cat(paste0("Processing GSE: ", gse, " (lowercase: ", current_gse, ")\n"))
  gse_strategy_row <- gse_strategy %>% filter(tolower(gse) == current_gse)
  cat(paste0("Found ", nrow(gse_strategy_row), " strategy rows for ", current_gse, "\n"))
  if (nrow(gse_strategy_row) > 0) {
    # GSE has a specific strategy
    strategy <- gse_strategy_row$strategy[1]
    cat(paste0("Using GSE-specific strategy: ", strategy, "\n"))
    
    if (strategy == "GENE_LEVEL_SYMBOL") {
      # Check if column names are numeric (potential pollution)
      if (mean(grepl("^[0-9.]+$", colnames(expr))) > 0.1) {
        cat("Detected numeric column names - attempting to fix\n")
        
        correct_colnames <- NULL
        
        # Priority 1: Try to get sample IDs from expr_norm file
        expr_norm_file <- file.path(base_dir, "data", "raw", "GSE152004_695_expr_norm.txt.gz")
        if (file.exists(expr_norm_file)) {
          cat("Found expr_norm file: ", expr_norm_file, "\n")
          # Read just the header to get sample names
          con <- gzfile(expr_norm_file, "r")
          header_line <- readLines(con, n = 1)
          close(con)
          
          # Parse header
          header <- strsplit(header_line, "\t")[[1]]
          # First element is gene ID column, rest are sample IDs
          sample_ids <- header[-1]
          cat("Found ", length(sample_ids), " sample IDs from expr_norm\n")
          
          if (length(sample_ids) == ncol(expr)) {
            correct_colnames <- sample_ids
            cat("Updated column names from expr_norm\n")
          } else {
            cat("❌ Column count mismatch: expr has ", ncol(expr), " columns, expr_norm has ", length(sample_ids), " columns\n")
          }
        } else {
          cat("❌ expr_norm file not found: ", expr_norm_file, "\n")
        }
        
        # Priority 2: Fall back to diet file if expr_norm failed
        if (is.null(correct_colnames)) {
          diet_file <- file.path(processed_diet_dir, paste0(gsub("&.*$", "", cohort), "_diet_expr.rds"))
          if (file.exists(diet_file)) {
            diet_expr <- readRDS(diet_file)
            diet_colnames <- colnames(diet_expr)
            cat("Found diet file with ", length(diet_colnames), " sample IDs\n")
            
            if (ncol(expr) == length(diet_colnames)) {
              correct_colnames <- diet_colnames
              cat("Updated column names from diet file\n")
            } else {
              cat("❌ Column count mismatch: expr has ", ncol(expr), " columns, diet has ", length(diet_colnames), " columns\n")
            }
          } else {
            cat("❌ Diet file not found: ", diet_file, "\n")
          }
        }
        
        # Apply correct column names if found
        if (!is.null(correct_colnames)) {
          colnames(expr) <- correct_colnames
        }
      }
      
      # Already gene-level symbol data, just clean and deduplicate
      feature_id_raw <- rownames(expr)
      feature_id_clean <- feature_id_raw
      gene_symbols <- feature_id_raw
      
      # Create ID mapping table
      id_map <- create_id_map(
        feature_id_raw = feature_id_raw,
        feature_id_clean = feature_id_clean,
        gene_symbol = gene_symbols,
        mapping_source = "GENE_LEVEL_SYMBOL",
        gpl = gpl,
        gse = gse,
        cohort = cohort
      )
      
      # Clean symbols: remove empty, "---", and NA
      valid_features <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
      
      if (sum(valid_features) == 0) {
        cat("⚠️  No valid features\n")
        # Save ID map
        id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
        write_tsv(id_map, id_map_output)
        cat("Saved ID mapping table: ", basename(id_map_output), "\n")
        expr_processed <- NULL
        status <- "FAILED"
        reason <- "No valid gene symbols"
      } else {
        # Subset expression matrix and ID map
        expr_processed <- expr[valid_features, ]
        gene_symbols <- gene_symbols[valid_features]
        id_map <- id_map[valid_features, ]
        
        # Set row names to gene symbols
        rownames(expr_processed) <- gene_symbols
        
        # Deduplicate by gene symbol - only if duplicate rate is significant
        gene_symbols <- as.character(gene_symbols)  # Ensure character type to avoid copying
        dup_rate <- mean(duplicated(gene_symbols))
        
        if (dup_rate > 0.01) {  # Only aggregate if >1% duplicates
          cat("Found duplicate gene symbols (rate: ", round(dup_rate, 3), "), aggregating with sum\n")
          # Use rowsum to aggregate duplicates
          aggregated <- rowsum(expr_processed, group = gene_symbols)
          
          # Update ID map for duplicates
          # Keep first occurrence of each gene symbol in ID map
          id_map <- id_map[!duplicated(id_map$gene_symbol), ]
          
          expr_processed <- aggregated
        } else if (any(duplicated(gene_symbols))) {
          # Low duplicate rate, just keep first occurrence without aggregation
          cat("Found low rate of duplicate gene symbols (rate: ", round(dup_rate, 3), "), keeping first occurrence\n")
          unique_idx <- !duplicated(gene_symbols)
          expr_processed <- expr_processed[unique_idx, ]
          id_map <- id_map[unique_idx, ]
          rownames(expr_processed) <- gene_symbols[unique_idx]
        }
        
        status <- "GENE_LEVEL_SYMBOL"
        reason <- "Processed as gene-level symbol data"
      }
    } else if (strategy == "MSTRG_SYMBOL_EXTRACTION") {
      # Special handling for MSTRG IDs with |SYMBOL format
      cat("Special handling for MSTRG IDs with |SYMBOL format\n")
      
      # Extract gene symbols from MSTRG IDs (format: MSTRG.XXX|SYMBOL)
      feature_id_raw <- rownames(expr)
      
      # Extract symbols from | separator
      gene_symbols <- sapply(feature_id_raw, function(x) {
        parts <- strsplit(x, "|", fixed = TRUE)[[1]]
        if (length(parts) > 1) {
          parts[2]
        } else {
          NA
        }
      })
      
      # Create ID mapping table
      id_map <- create_id_map(
        feature_id_raw = feature_id_raw,
        feature_id_clean = feature_id_raw,
        gene_symbol = gene_symbols,
        mapping_source = "MSTRG_SYMBOL_EXTRACTION",
        gpl = gpl,
        gse = gse,
        cohort = cohort
      )
      
      # Clean symbols: remove empty, "---", and NA
      valid_features <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
      
      # Calculate discard rate
      total_features <- length(feature_id_raw)
      valid_count <- sum(valid_features)
      discard_rate <- 1 - (valid_count / total_features)
      cat(paste0("Discard rate: ", round(discard_rate * 100, 2), "% (", total_features - valid_count, "/", total_features, ")\n"))
      
      if (sum(valid_features) == 0) {
        cat("⚠️  No valid features\n")
        # Save ID map
        id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
        write_tsv(id_map, id_map_output)
        cat("Saved ID mapping table: ", basename(id_map_output), "\n")
        expr_processed <- NULL
        status <- "FAILED"
        reason <- "No valid gene symbols"
      } else {
        # Subset expression matrix and ID map
        expr_processed <- expr[valid_features, ]
        gene_symbols <- gene_symbols[valid_features]
        id_map <- id_map[valid_features, ]
        
        # Set row names to gene symbols
        rownames(expr_processed) <- gene_symbols
        
        # Deduplicate by gene symbol - only if duplicate rate is significant
        gene_symbols <- as.character(gene_symbols)  # Ensure character type to avoid copying
        dup_rate <- mean(duplicated(gene_symbols))
        
        if (dup_rate > 0.01) {  # Only aggregate if >1% duplicates
          cat("Found duplicate gene symbols (rate: ", round(dup_rate, 3), "), aggregating with sum\n")
          # Use rowsum to aggregate duplicates
          aggregated <- rowsum(expr_processed, group = gene_symbols)
          
          # Update ID map for duplicates
          # Keep first occurrence of each gene symbol in ID map
          id_map <- id_map[!duplicated(id_map$gene_symbol), ]
          
          expr_processed <- aggregated
        } else if (any(duplicated(gene_symbols))) {
          # Low duplicate rate, just keep first occurrence without aggregation
          cat("Found low rate of duplicate gene symbols (rate: ", round(dup_rate, 3), "), keeping first occurrence\n")
          unique_idx <- !duplicated(gene_symbols)
          expr_processed <- expr_processed[unique_idx, ]
          id_map <- id_map[unique_idx, ]
          rownames(expr_processed) <- gene_symbols[unique_idx]
        }
        
        status <- "MSTRG_SYMBOL_EXTRACTION"
        reason <- "Processed MSTRG IDs with |SYMBOL extraction"
      }
    }
    
    # Write processed expression and ID mapping table for GSE-specific strategy
    if (!is.null(expr_processed)) {
      # Calculate mapping rate
      raw_rows <- nrow(expr)
      mapped_rows <- nrow(expr_processed)
      mapped_rate <- mapped_rows / raw_rows
      
      # Threshold gating
      if (nrow(expr_processed) >= 2000 && mapped_rate >= 0.2) {
        # Calculate NA percentage
        na_percentage <- sum(is.na(expr_processed)) / (nrow(expr_processed) * ncol(expr_processed)) * 100
        if (!is.na(na_percentage) && na_percentage < 50) {
          # Hard assertion: ensure final row names are SYMBOL
          stopifnot(is.matrix(expr_processed))
          final_type <- determine_id_type(rownames(expr_processed))
          if (final_type != "SYMBOL") {
            cat("❌ Final rownames not SYMBOL (", final_type, ") for cohort=", cohort, "\n")
            status <- "FAILED_FINAL_NOT_SYMBOL"
            reason <- paste("Final rownames not SYMBOL (", final_type, ")")
          } else {
            # Calculate HGNC symbol hit rate for QC
            hgnc_hit_rate <- calculate_hgnc_hit_rate(rownames(expr_processed))
            cat(paste0("HGNC symbol hit rate: ", round(hgnc_hit_rate, 3), "\n"))
            
            if (hgnc_hit_rate < 0.3) {
              cat("❌ Low HGNC symbol hit rate (", round(hgnc_hit_rate, 3), "), skipping FULL file\n")
              status <- "FAILED_MAPPING_QC"
              reason <- "Row names not recognized as HGNC symbols (low org.Hs.eg.db hit rate)"
            } else {
              # Assertions for output matrix
              stopifnot(!is.null(rownames(expr_processed)), !is.null(colnames(expr_processed)))
              stopifnot(length(unique(rownames(expr_processed))) == nrow(expr_processed))
              
              # Save processed expression with compression to reduce disk usage and memory pressure
              full_output_file <- file.path(processed_full_dir, paste0("expr_", cohort, "__FULL.rds"))
              saveRDS(expr_processed, full_output_file, compress = "xz")  # xz compression for better space savings
              cat("Saved processed expression (compressed): ", basename(full_output_file), " | dim=", nrow(expr_processed), "x", ncol(expr_processed), "\n")
              
              # Save ID mapping table
              id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
              write_tsv(id_map, id_map_output)
              cat("Saved ID mapping table: ", basename(id_map_output), "\n")
              
              cat(paste0("Mapping rate: ", round(mapped_rate, 3), " (", mapped_rows, "/", raw_rows, ")\n"))
              
              # For GSE65204, create a specific mapping report
              if (grepl("GSE65204", gse)) {
                gse65204_log_file <- file.path(logs_dir, "00b_map_gse65204.log")
                sink(gse65204_log_file, append = TRUE)
                cat("\n=== GSE65204 Mapping Report ===\n")
                cat(paste0("Execution time: ", Sys.time(), "\n"))
                cat(paste0("Original probe count: ", raw_rows, "\n"))
                cat(paste0("Mapped unique GeneSymbol count: ", mapped_rows, "\n"))
                cat(paste0("Mapping rate: ", round(mapped_rate, 3), "\n"))
                cat(paste0("HGNC symbol hit rate: ", round(hgnc_hit_rate, 3), "\n"))
                
                # Calculate missing probe percentage
                missing_probe_count <- raw_rows - mapped_rows
                missing_probe_percentage <- (missing_probe_count / raw_rows) * 100
                cat(paste0("Missing probes: ", missing_probe_count, " (", round(missing_probe_percentage, 2), "%)\n"))
                
                # Identify and list unmapped probes
                unmapped_probes <- id_map[is.na(id_map$gene_symbol) | id_map$gene_symbol == "" | id_map$gene_symbol == "---", "feature_id_raw"]
                if (length(unmapped_probes) > 0) {
                  cat(paste0("Unmapped probe count: ", length(unmapped_probes), "\n"))
                  if (length(unmapped_probes) <= 10) {
                    cat("Unmapped probes: ", paste(unmapped_probes, collapse = ", "), "\n")
                  } else {
                    cat("First 10 unmapped probes: ", paste(head(unmapped_probes, 10), collapse = ", "), "...\n")
                  }
                }
                
                cat("=== GSE65204 Mapping Report End ===\n")
                sink()
                cat(paste0("Saved GSE65204 mapping report: ", basename(gse65204_log_file), "\n"))
              }
            }
          }
        } else {
          cat("❌ High NA percentage (", round(na_percentage, 2), "%), skipping\n")
          status <- "FAILED"
          reason <- paste("High NA percentage:", round(na_percentage, 2), "%")
        }
      } else {
        cat("❌ Low mapping quality:\n")
        cat(paste0("  - Mapped rows: ", mapped_rows, " (required: >= 2000)\n"))
        cat(paste0("  - Mapping rate: ", round(mapped_rate, 3), " (required: >= 0.2)\n"))
        status <- "FAILED_LOW_MAPPING"
        reason <- paste("Low mapping quality: rows=", mapped_rows, ", rate=", round(mapped_rate, 3))
      }
    } else {
      # Save ID mapping table even if expression processing failed
      id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
      write_tsv(id_map, id_map_output)
      cat("Saved ID mapping table: ", basename(id_map_output), "\n")
      
      # For GSE65204, create a specific mapping report even if failed
      if (grepl("GSE65204", gse)) {
        gse65204_log_file <- file.path(logs_dir, "00b_map_gse65204.log")
        sink(gse65204_log_file, append = TRUE)
        cat("\n=== GSE65204 Mapping Report ===\n")
        cat(paste0("Execution time: ", Sys.time(), "\n"))
        cat(paste0("Mapping failed: ", reason, "\n"))
        cat("=== GSE65204 Mapping Report End ===\n")
        sink()
        cat(paste0("Saved GSE65204 mapping report: ", basename(gse65204_log_file), "\n"))
      }
    }
    
    # Calculate mapping statistics
    n_row_raw <- if (is.null(expr)) 0 else nrow(expr)
    n_col_raw <- if (is.null(expr)) 0 else ncol(expr)
    n_row_full <- if (is.null(expr_processed)) 0 else nrow(expr_processed)
    mapped_rate <- if (n_row_raw > 0) n_row_full / n_row_raw else 0
    
    # Determine final ID type
    final_id_type <- if (is.null(expr_processed)) "UNKNOWN" else determine_id_type(rownames(expr_processed))
    
    # Determine key used
    key_used <- if (gpl == "GPL23961") "SPOT_ID" else "ID"
    
    # Update GPL status table
    gpl_status_df <- rbind(gpl_status_df, data.frame(
      gpl_id = gpl,
      cohort = cohort,
      status = status,
      reason = reason,
      n_row_raw = n_row_raw,
      n_col_raw = n_col_raw,
      n_row_full = n_row_full,
      mapped_rate = mapped_rate,
      raw_id_type = id_type,
      final_id_type = final_id_type,
      key_used = key_used,
      stringsAsFactors = FALSE
    ))
    
    # Memory cleanup and move to next cohort
    rm(expr)
    if (exists("expr_processed")) {
      rm(expr_processed)
    }
    if (exists("id_map")) {
      rm(id_map)
    }
    gc()
    Sys.sleep(1)
    next
  }
  
  # Get platform strategy
  strategy_row <- platform_strategy %>% filter(gpl == !!gpl)
  
  # Initialize variables
  expr_processed <- NULL
  id_map <- NULL
  status <- "FAILED"
  reason <- "No strategy applied"
  
  if (nrow(strategy_row) > 0) {
    # Platform has a specific strategy
    strategy <- strategy_row$strategy[1]
    key_col <- strategy_row$key[1]
    symbol_col <- strategy_row$symbol_col[1]
    
    cat(paste0("Using strategy: ", strategy, "\n"))
    
    if (strategy == "HTS_COUNTS") {
      # Process HTS counts
      result <- process_hts_counts(expr, id_type, gpl, gse, cohort)
      expr_processed <- result$expr
      id_map <- result$id_map
      status <- "HTS_COUNTS"
      reason <- "Processed as HTS counts with symbol mapping"
    } else if (strategy == "GPL_TABLE_SPOTID" && gpl == "GPL23961") {
      # Process GPL23961 with SPOT_ID mapping
      mapping <- build_probe_mapping(gpl, key_col, symbol_col)
      if (!is.null(mapping)) {
        result <- process_gpl23961(expr, mapping, gpl, gse, cohort)
        expr_processed <- result$expr
        id_map <- result$id_map
        if (!is.null(expr_processed)) {
          status <- "GPL23961_SPOTID"
          reason <- "Processed with SPOT_ID mapping"
        } else {
          cat("❌ Failed to process GPL23961 mapping\n")
          expr_processed <- NULL
          status <- "FAILED_GPL23961_MAPPING_COLS"
          reason <- "Failed GPL23961 mapping: missing SPOT_ID or ORF columns"
        }
      } else {
        cat("❌ Failed to build mapping for GPL23961\n")
        expr_processed <- NULL
        id_map <- create_id_map(
          feature_id_raw = rownames(expr),
          feature_id_clean = rownames(expr),
          gene_symbol = rep(NA, nrow(expr)),
          mapping_source = "GPL_TABLE_SPOTID",
          gpl = gpl,
          gse = gse,
          cohort = cohort
        )
        # Save ID map
        id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
        write_tsv(id_map, id_map_output)
        cat("Saved ID mapping table: ", basename(id_map_output), "\n")
        status <- "FAILED"
        reason <- "Failed to build mapping"
      }
    } else if (strategy == "GPL_TABLE") {
      # Process general probe platform with specified columns
      mapping <- build_probe_mapping(gpl, key_col, symbol_col)
      if (!is.null(mapping)) {
        result <- process_general_probe(expr, mapping, gpl, gse, cohort)
        expr_processed <- result$expr
        id_map <- result$id_map
        status <- "GPL_TABLE"
        reason <- "Processed with GPL table mapping"
      } else {
        cat("❌ Failed to build mapping for GPL", gpl, "\n")
        expr_processed <- NULL
        id_map <- create_id_map(
          feature_id_raw = rownames(expr),
          feature_id_clean = rownames(expr),
          gene_symbol = rep(NA, nrow(expr)),
          mapping_source = "GPL_TABLE",
          gpl = gpl,
          gse = gse,
          cohort = cohort
        )
        # Save ID map
        id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
        write_tsv(id_map, id_map_output)
        cat("Saved ID mapping table: ", basename(id_map_output), "\n")
        status <- "FAILED"
        reason <- "Failed to build mapping"
      }
    } else {
      # Unknown strategy
      cat("❌ Unknown strategy: ", strategy, "\n")
      expr_processed <- NULL
      id_map <- create_id_map(
        feature_id_raw = rownames(expr),
        feature_id_clean = rownames(expr),
        gene_symbol = rep(NA, nrow(expr)),
        mapping_source = "UNKNOWN_STRATEGY",
        gpl = gpl,
        gse = gse,
        cohort = cohort
      )
      # Save ID map
      id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
      write_tsv(id_map, id_map_output)
      cat("Saved ID mapping table: ", basename(id_map_output), "\n")
      status <- "FAILED"
      reason <- "Unknown strategy"
    }
  } else {
    # General probe platform
    if (id_type == "PROBE") {
      # Try to build mapping with default columns
      mapping <- build_probe_mapping(gpl, "id", "gene_symbol")
      if (!is.null(mapping)) {
        result <- process_general_probe(expr, mapping, gpl, gse, cohort)
        expr_processed <- result$expr
        id_map <- result$id_map
        status <- "GPL_PROBE"
        reason <- "Processed with general probe mapping"
      } else {
        cat("❌ Failed to build mapping for GPL", gpl, "\n")
        expr_processed <- NULL
        id_map <- create_id_map(
          feature_id_raw = rownames(expr),
          feature_id_clean = rownames(expr),
          gene_symbol = rep(NA, nrow(expr)),
          mapping_source = "GPL_TABLE",
          gpl = gpl,
          gse = gse,
          cohort = cohort
        )
        # Save ID map
        id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
        write_tsv(id_map, id_map_output)
        cat("Saved ID mapping table: ", basename(id_map_output), "\n")
        status <- "FAILED"
        reason <- "Failed to build mapping"
      }
    } else {
      # Already gene-level data, map to symbols if needed
      feature_id_raw <- rownames(expr)
      feature_id_clean <- feature_id_raw
      
      if (id_type == "ENSG") {
        feature_id_clean <- clean_ensg(feature_id_clean)
      } else if (id_type == "ENTREZ") {
        feature_id_clean <- clean_entrez(feature_id_clean)
      }
      
      # Map to symbols
      gene_symbols <- map_ids_to_symbols(feature_id_clean, id_type)
      
      # Create ID mapping table
      id_map <- create_id_map(
        feature_id_raw = feature_id_raw,
        feature_id_clean = feature_id_clean,
        gene_symbol = gene_symbols,
        mapping_source = "GENE_LEVEL",
        gpl = gpl,
        gse = gse,
        cohort = cohort
      )
      
      # Remove features with no gene mapping
      valid_features <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
      
      if (sum(valid_features) == 0) {
        cat("⚠️  No valid features with gene mapping\n")
        # Save ID map
        id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
        write_tsv(id_map, id_map_output)
        cat("Saved ID mapping table: ", basename(id_map_output), "\n")
        expr_processed <- NULL
        status <- "FAILED"
        reason <- "No valid gene symbols"
      } else {
        # Subset expression matrix and ID map
        expr_processed <- expr[valid_features, ]
        gene_symbols <- gene_symbols[valid_features]
        id_map <- id_map[valid_features, ]
        
        # Set row names to gene symbols
        rownames(expr_processed) <- gene_symbols
        
        # Deduplicate by gene symbol - only if duplicate rate is significant
        gene_symbols <- as.character(gene_symbols)  # Ensure character type to avoid copying
        dup_rate <- mean(duplicated(gene_symbols))
        
        if (dup_rate > 0.01) {  # Only aggregate if >1% duplicates
          cat("Found duplicate gene symbols (rate: ", round(dup_rate, 3), "), aggregating with sum\n")
          # Use rowsum to aggregate duplicates
          aggregated <- rowsum(expr_processed, group = gene_symbols)
          expr_processed <- aggregated
          
          # Update ID map for duplicates
          id_map <- id_map[!duplicated(id_map$gene_symbol), ]
        } else if (any(duplicated(gene_symbols))) {
          # Low duplicate rate, just keep first occurrence without aggregation
          cat("Found low rate of duplicate gene symbols (rate: ", round(dup_rate, 3), "), keeping first occurrence\n")
          unique_idx <- !duplicated(gene_symbols)
          expr_processed <- expr_processed[unique_idx, ]
          id_map <- id_map[unique_idx, ]
          rownames(expr_processed) <- gene_symbols[unique_idx]
        }
        
        status <- "GENE_LEVEL"
        reason <- "Processed as gene-level data with symbol mapping"
      }
    }
  }
  
  # Write processed expression and ID mapping table
  if (!is.null(expr_processed)) {
    # Calculate mapping rate
    raw_rows <- nrow(expr)
    mapped_rows <- nrow(expr_processed)
    mapped_rate <- mapped_rows / raw_rows
    
    # Threshold gating
    if (nrow(expr_processed) >= 2000 && mapped_rate >= 0.2) {
      # Calculate NA percentage
      na_percentage <- sum(is.na(expr_processed)) / (nrow(expr_processed) * ncol(expr_processed)) * 100
      if (!is.na(na_percentage) && na_percentage < 50) {
        # Hard assertion: ensure final row names are SYMBOL
        stopifnot(is.matrix(expr_processed))
        final_type <- determine_id_type(rownames(expr_processed))
        if (final_type != "SYMBOL") {
          cat("❌ Final rownames not SYMBOL (", final_type, ") for cohort=", cohort, "\n")
          status <- "FAILED_FINAL_NOT_SYMBOL"
          reason <- paste("Final rownames not SYMBOL (", final_type, ")")
        } else {
          # Calculate HGNC symbol hit rate for QC
          hgnc_hit_rate <- calculate_hgnc_hit_rate(rownames(expr_processed))
          cat(paste0("HGNC symbol hit rate: ", round(hgnc_hit_rate, 3), "\n"))
          
          if (hgnc_hit_rate < 0.3) {
            cat("❌ Low HGNC symbol hit rate (", round(hgnc_hit_rate, 3), "), skipping FULL file\n")
            status <- "FAILED_MAPPING_QC"
            reason <- "Row names not recognized as HGNC symbols (low org.Hs.eg.db hit rate)"
          } else {
            # Check for numeric column names before saving
            if (mean(grepl("^[0-9.]+$", colnames(expr_processed))) > 0.1) {
              cat("❌ Invalid colnames: look like numeric expression values, header parse failure\n")
              cat("First 10 column names: ", head(colnames(expr_processed), 10), "\n")
              cat("First 10 row names: ", head(rownames(expr_processed), 10), "\n")
              stop("Invalid colnames: look like numeric expression values, header parse failure")
            }
            
            # Additional assertions for GSE152004
            if (tolower(gse) == "gse152004") {
              # Assert that column names are not numeric
              if (mean(grepl("^[0-9.]+$", colnames(expr_processed))) > 0.1) {
                stop("Invalid colnames: look like numeric expression values, header parse failure")
              }
              # Assert that column count is 695
              if (ncol(expr_processed) != 695) {
                stop(paste("Invalid column count: expected 695, got", ncol(expr_processed)))
              }
              cat("✅ GSE152004 assertions passed: column names valid and count=695\n")
            }
            
            # Assertions for output matrix
            stopifnot(!is.null(rownames(expr_processed)), !is.null(colnames(expr_processed)))
            stopifnot(length(unique(rownames(expr_processed))) == nrow(expr_processed))
            
            # Assertions for output matrix
            stopifnot(!is.null(rownames(expr_processed)), !is.null(colnames(expr_processed)))
            stopifnot(length(unique(rownames(expr_processed))) == nrow(expr_processed))
            
            # Assertions for output matrix
            stopifnot(!is.null(rownames(expr_processed)), !is.null(colnames(expr_processed)))
            stopifnot(length(unique(rownames(expr_processed))) == nrow(expr_processed))
            
            # Save processed expression with compression to reduce disk usage and memory pressure
            full_output_file <- file.path(processed_full_dir, paste0("expr_", cohort, "__FULL.rds"))
            saveRDS(expr_processed, full_output_file, compress = "xz")  # xz compression for better space savings
            cat("Saved processed expression (compressed): ", basename(full_output_file), " | dim=", nrow(expr_processed), "x", ncol(expr_processed), "\n")
            
            # Save ID mapping table
            id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
            write_tsv(id_map, id_map_output)
            cat("Saved ID mapping table: ", basename(id_map_output), "\n")
            
            cat(paste0("Mapping rate: ", round(mapped_rate, 3), " (", mapped_rows, "/", raw_rows, ")\n"))
            
            # For GSE65204, create a specific mapping report
            if (grepl("GSE65204", gse)) {
              gse65204_log_file <- file.path(logs_dir, "00b_map_gse65204.log")
              sink(gse65204_log_file, append = TRUE)
              cat("\n=== GSE65204 Mapping Report ===\n")
              cat(paste0("Execution time: ", Sys.time(), "\n"))
              cat(paste0("Original probe count: ", raw_rows, "\n"))
              cat(paste0("Mapped unique GeneSymbol count: ", mapped_rows, "\n"))
              cat(paste0("Mapping rate: ", round(mapped_rate, 3), "\n"))
              cat(paste0("HGNC symbol hit rate: ", round(hgnc_hit_rate, 3), "\n"))
              
              # Calculate missing probe percentage
              missing_probe_count <- raw_rows - mapped_rows
              missing_probe_percentage <- (missing_probe_count / raw_rows) * 100
              cat(paste0("Missing probes: ", missing_probe_count, " (", round(missing_probe_percentage, 2), "%)\n"))
              
              # Identify and list unmapped probes
              unmapped_probes <- id_map[is.na(id_map$gene_symbol) | id_map$gene_symbol == "" | id_map$gene_symbol == "---", "feature_id_raw"]
              if (length(unmapped_probes) > 0) {
                cat(paste0("Unmapped probe count: ", length(unmapped_probes), "\n"))
                if (length(unmapped_probes) <= 10) {
                  cat("Unmapped probes: ", paste(unmapped_probes, collapse = ", "), "\n")
                } else {
                  cat("First 10 unmapped probes: ", paste(head(unmapped_probes, 10), collapse = ", "), "...\n")
                }
              }
              
              cat("=== GSE65204 Mapping Report End ===\n")
              sink()
              cat(paste0("Saved GSE65204 mapping report: ", basename(gse65204_log_file), "\n"))
            }
          }
        }
      } else {
        cat("❌ High NA percentage (", round(na_percentage, 2), "%), skipping\n")
        status <- "FAILED"
        reason <- paste("High NA percentage:", round(na_percentage, 2), "%")
      }
    } else {
      cat("❌ Low mapping quality:\n")
      cat(paste0("  - Mapped rows: ", mapped_rows, " (required: >= 2000)\n"))
      cat(paste0("  - Mapping rate: ", round(mapped_rate, 3), " (required: >= 0.2)\n"))
      status <- "FAILED_LOW_MAPPING"
      reason <- paste("Low mapping quality: rows=", mapped_rows, ", rate=", round(mapped_rate, 3))
    }
  } else {
    # Save ID mapping table even if expression processing failed
    id_map_output <- file.path(processed_full_dir, paste0("id_map_", cohort, ".tsv.gz"))
    write_tsv(id_map, id_map_output)
    cat("Saved ID mapping table: ", basename(id_map_output), "\n")
  }
  
  # Calculate mapping statistics
  n_row_raw <- if (is.null(expr)) 0 else nrow(expr)
  n_col_raw <- if (is.null(expr)) 0 else ncol(expr)
  n_row_full <- if (is.null(expr_processed)) 0 else nrow(expr_processed)
  mapped_rate <- if (n_row_raw > 0) n_row_full / n_row_raw else 0
  
  # Determine final ID type
  final_id_type <- if (is.null(expr_processed)) "UNKNOWN" else determine_id_type(rownames(expr_processed))
  
  # Determine key used
  key_used <- if (gpl == "GPL23961") "SPOT_ID" else "ID"
  
  # Update GPL status table
  gpl_status_df <- rbind(gpl_status_df, data.frame(
    gpl_id = gpl,
    cohort = cohort,
    status = status,
    reason = reason,
    n_row_raw = n_row_raw,
    n_col_raw = n_col_raw,
    n_row_full = n_row_full,
    mapped_rate = mapped_rate,
    raw_id_type = id_type,
    final_id_type = final_id_type,
    key_used = key_used,
    stringsAsFactors = FALSE
  ))
  
  # Memory cleanup
  rm(expr)
  if (exists("expr_processed")) {
    rm(expr_processed)
  }
  if (exists("mapping")) {
    rm(mapping)
  }
  gc()
  Sys.sleep(1)
}

# Save GPL status table
gpl_status_output <- file.path(logs_dir, "00_gpl_status.tsv")
write.table(gpl_status_df, gpl_status_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved GPL status table: ", basename(gpl_status_output), "\n")

cat("\n=== ID standardization and mapping completed ===\n")
sink()

# Function to validate output files (memory-efficient version)
validate_outputs <- function() {
  cat("\n=== Validating output files ===\n")
  
  # Get all FULL files
  full_files <- list.files(processed_full_dir, pattern = "__FULL\\.rds$", full.names = TRUE)
  
  all_valid <- TRUE
  
  # Limit validation to first 10 files to reduce memory usage
  validation_files <- head(full_files, 10)
  if (length(full_files) > 10) {
    cat(paste0("⚠️  Validating first 10 of ", length(full_files), " files to reduce memory usage\n"))
  }
  
  for (full_file in validation_files) {
    cat(paste0("\nValidating: ", basename(full_file), "\n"))
    
    # Read the file - use a tryCatch to handle errors
    expr <- tryCatch({
      readRDS(full_file)
    }, error = function(e) {
      cat("❌ Failed to read file: ", conditionMessage(e), "\n")
      all_valid <<- FALSE
      return(NULL)
    })
    
    if (is.null(expr)) next
    
    # Check if it's a matrix
    if (!is.matrix(expr)) {
      cat("❌ Not a matrix\n")
      all_valid <- FALSE
    } else {
      # Check dimensions (lightweight check)
      cat(paste0("  Dimensions: ", nrow(expr), "x", ncol(expr), "\n"))
      
      # Check number of rows
      if (nrow(expr) < 2000) {
        cat(paste0("❌ Insufficient rows: ", nrow(expr), " < 2000\n"))
        all_valid <- FALSE
      }
      
      # Check row names (lightweight check - only first 100)
      row_names <- head(rownames(expr), 100)
      final_id_type <- determine_id_type(row_names)
      if (final_id_type != "SYMBOL") {
        cat(paste0("❌ Final ID type is not SYMBOL: ", final_id_type, "\n"))
        all_valid <- FALSE
      }
      
      # Check for duplicate row names in sample
      if (any(duplicated(row_names))) {
        cat("❌ Duplicate row names in sample\n")
        all_valid <- FALSE
      }
      
      # Check for empty or '---' row names in sample
      if (any(row_names == "")) {
        cat("❌ Empty row names in sample\n")
        all_valid <- FALSE
      }
      if (any(row_names == "---")) {
        cat("❌ Row names contain '---' in sample\n")
        all_valid <- FALSE
      }
    }
    
    # Clean up immediately to free memory
    rm(expr)
    gc()
    
    if (all_valid) {
      cat("✓ Valid\n")
    }
  }
  
  if (all_valid) {
    cat("\n=== All validated output files are valid ===\n")
  } else {
    cat("\n=== Some output files are invalid ===\n")
    stop("Validation failed")
  }
}

# Validate outputs
validate_outputs()

# Memory cleanup
rm(list = ls(pattern = "^(matrix|expr|data|probe|gpl|cohort|strategy)", all.names = TRUE))
gc()
