#!/usr/bin/env Rscript

# Script 03: Logging functions to automatically record all output figures to OUTPUTS_INDEX.csv
# Purpose: Record information about each generated figure, including file name, path, type, generation time, etc.

# Load necessary packages
library(tidyverse)

# Initialize output index file
init_output_index <- function(index_file = NULL) {
  if (is.null(index_file)) {
    # Default path
    base_dir <- getwd()
    index_file <- file.path(base_dir, "output", "OUTPUTS_INDEX.csv")
  }
  
  # Check if directory exists
  index_dir <- dirname(index_file)
  if (!dir.exists(index_dir)) {
    dir.create(index_dir, recursive = TRUE)
  }
  
  # Check if file exists, create if it doesn't
  if (!file.exists(index_file)) {
    index_df <- data.frame(
      file_name = character(),
      file_path = character(),
      file_type = character(),
      figure_id = character(),
      description = character(),
      script_name = character(),
      generated_at = as.POSIXct(character()),
      resolution = character(),
      file_size = numeric(),
      gene_coverage = character(),
      weights_hash = character(),
      notes = character(),
      stringsAsFactors = FALSE
    )
    write_csv(index_df, index_file)
    cat(paste0("Initializing output index file: ", index_file, "\n"))
  } else {
    cat(paste0("Output index file already exists: ", index_file, "\n"))
  }
  
  return(index_file)
}

# Log output file
log_output <- function(
  file_path,
  file_type = c("figure", "table", "other"),
  figure_id = "",
  description = "",
  script_name = "",
  resolution = "",
  gene_coverage = "",
  weights_hash = "",
  notes = "",
  index_file = NULL
) {
  # Verify file exists
  if (!file.exists(file_path)) {
    warning(paste0("File does not exist: ", file_path))
    return(NULL)
  }
  
  # Determine index file path
  if (is.null(index_file)) {
    base_dir <- getwd()
    index_file <- file.path(base_dir, "output", "OUTPUTS_INDEX.csv")
  }
  
  # Initialize index file if it doesn't exist
  if (!file.exists(index_file)) {
    init_output_index(index_file)
  }
  
  # Read existing index
  index_df <- read_csv(index_file)
  
  # Extract file name
  file_name <- basename(file_path)
  
  # Calculate file size (bytes)
  file_size <- file.size(file_path)
  
  # Get current time
  generated_at <- Sys.time()
  
  # Create new record
  new_record <- data.frame(
    file_name = file_name,
    file_path = file_path,
    file_type = match.arg(file_type),
    figure_id = figure_id,
    description = description,
    script_name = script_name,
    generated_at = generated_at,
    resolution = resolution,
    file_size = file_size,
    gene_coverage = gene_coverage,
    weights_hash = weights_hash,
    notes = notes,
    stringsAsFactors = FALSE
  )
  
  # Add new record to index
  if (nrow(index_df) > 0) {
    # Convert columns to correct types if needed
    if (is.character(index_df$generated_at)) {
      index_df$generated_at <- as.POSIXct(index_df$generated_at)
    }
    if (is.character(index_df$file_size)) {
      index_df$file_size <- as.numeric(index_df$file_size)
    }
    index_df <- bind_rows(index_df, new_record)
  } else {
    # If index is empty, use new_record as the first row
    index_df <- new_record
  }
  
  # Save updated index
  write_csv(index_df, index_file)
  
  cat(paste0("Output file logged: ", file_name, "\n"))
  
  return(new_record)
}

# Batch log output files
batch_log_output <- function(
  file_paths,
  file_type = c("figure", "table", "other"),
  figure_ids = NULL,
  descriptions = NULL,
  script_name = "",
  resolutions = NULL,
  gene_coverages = NULL,
  weights_hashes = NULL,
  notes = "",
  index_file = NULL
) {
  # Validate input lengths
  n_files <- length(file_paths)
  if (!is.null(figure_ids) && length(figure_ids) != n_files) {
    stop("figure_ids length must match file_paths length")
  }
  if (!is.null(descriptions) && length(descriptions) != n_files) {
    stop("descriptions length must match file_paths length")
  }
  if (!is.null(resolutions) && length(resolutions) != n_files) {
    stop("resolutions length must match file_paths length")
  }
  if (!is.null(gene_coverages) && length(gene_coverages) != n_files) {
    stop("gene_coverages length must match file_paths length")
  }
  if (!is.null(weights_hashes) && length(weights_hashes) != n_files) {
    stop("weights_hashes length must match file_paths length")
  }
  
  # Log each file
  records <- list()
  for (i in 1:n_files) {
    record <- log_output(
      file_path = file_paths[i],
      file_type = file_type,
      figure_id = ifelse(is.null(figure_ids), "", figure_ids[i]),
      description = ifelse(is.null(descriptions), "", descriptions[i]),
      script_name = script_name,
      resolution = ifelse(is.null(resolutions), "", resolutions[i]),
      gene_coverage = ifelse(is.null(gene_coverages), "", gene_coverages[i]),
      weights_hash = ifelse(is.null(weights_hashes), "", weights_hashes[i]),
      notes = notes,
      index_file = index_file
    )
    records[[i]] <- record
  }
  
  return(records)
}

# Example usage
# source("R/03_index_logger.R")
# 
# # Initialize index file
# init_output_index()
# 
# # Log single output file
# log_output(
#   file_path = "output/figures_main/fig1b_density.pdf",
#   file_type = "figure",
#   figure_id = "Fig 1B",
#   description = "Z-score density plot",
#   script_name = "07_fig1_discovery.R",
#   resolution = "300 dpi",
#   notes = "Discovery cohort analysis"
# )
# 
# # Batch log output files
# file_paths <- c(
#   "output/figures_main/fig2a_log2fc.pdf",
#   "output/figures_main/fig2b_heatmap.pdf"
# )
# figure_ids <- c("Fig 2A", "Fig 2B")
# descriptions <- c("log2FC plot", "Heatmap")
# batch_log_output(
#   file_paths = file_paths,
#   file_type = "figure",
#   figure_ids = figure_ids,
#   descriptions = descriptions,
#   script_name = "08_fig2_replication.R"
# )
