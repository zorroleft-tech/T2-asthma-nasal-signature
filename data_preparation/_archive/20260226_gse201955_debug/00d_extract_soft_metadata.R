#!/usr/bin/env Rscript

# Script 00d: Extract sample metadata from family.soft.gz file
# Function: Parse GSE201955_family.soft.gz and extract sample metadata

# Load required packages
library(tidyverse)
library(data.table)

# Set paths
# Check if we're in the project root or data_preparation directory
current_dir <- getwd()
if (basename(current_dir) == "data_preparation") {
  # Running from data_preparation directory
  base_dir <- dirname(current_dir)
} else {
  # Running from project root directory
  base_dir <- current_dir
}
raw_dir <- file.path(base_dir, "data", "raw")
processed_full_dir <- file.path(base_dir, "data", "processed_full")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(processed_full_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00d_extract_soft_metadata.log")
sink(log_file, append = TRUE)

cat("\n=== Starting soft metadata extraction ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))

# Define GSE201955 soft file path
soft_file <- file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_family.soft.gz")

if (!file.exists(soft_file)) {
  cat(paste0("❌ Soft file not found: ", soft_file, "\n"))
  sink()
  stop("Soft file not found")
}

cat(paste0("Reading soft file: ", basename(soft_file), "\n"))

# Read soft file line by line
con <- gzfile(soft_file, "rt")
lines <- readLines(con)
close(con)

# Initialize variables
samples <- list()
current_sample <- list()
in_sample <- FALSE

# Parse soft file
for (line in lines) {
  # Check for start of sample section
  if (grepl("^!Sample_geo_accession", line)) {
    # If we were in a sample, save it
    if (in_sample && length(current_sample) > 0) {
      samples <- c(samples, list(current_sample))
    }
    # Start new sample
    current_sample <- list()
    in_sample <- TRUE
    # Extract GSM
    gsm <- sub("^!Sample_geo_accession[[:space:]]+", "", line)
    current_sample$GSM <- gsm
  } else if (in_sample) {
    # Extract title
    if (grepl("^!Sample_title", line)) {
      title <- sub("^!Sample_title[[:space:]]+", "", line)
      current_sample$title <- title
    } 
    # Extract source_name_ch1
    else if (grepl("^!Sample_source_name_ch1", line)) {
      source_name <- sub("^!Sample_source_name_ch1[[:space:]]+", "", line)
      current_sample$source_name_ch1 <- source_name
    }
    # Extract characteristics_ch1 (including multiple values)
    else if (grepl("^!Sample_characteristics_ch1", line)) {
      characteristic <- sub("^!Sample_characteristics_ch1[[:space:]]+", "", line)
      # If characteristics already exists, append to it
      if ("characteristics_ch1" %in% names(current_sample)) {
        current_sample$characteristics_ch1 <- c(current_sample$characteristics_ch1, characteristic)
      } else {
        current_sample$characteristics_ch1 <- characteristic
      }
    }
    # Extract other relevant fields if needed
    else if (grepl("^!Sample_", line) && !grepl("^!Sample_geo_accession", line)) {
      # Skip other sample fields for now
    }
    # Check for end of sample section (start of new section or end of file)
    else if (grepl("^!", line) && !grepl("^!Sample_", line)) {
      # Save current sample
      if (length(current_sample) > 0) {
        samples <- c(samples, list(current_sample))
      }
      in_sample <- FALSE
    }
  }
}

# Save the last sample if any
if (in_sample && length(current_sample) > 0) {
  samples <- c(samples, list(current_sample))
}

# Convert to data frame
cat(paste0("Found ", length(samples), " samples in soft file\n"))

# Initialize vectors for each column
gsm_vec <- character(length(samples))
title_vec <- character(length(samples))
source_name_vec <- character(length(samples))
characteristics_vec <- character(length(samples))

# Populate vectors
for (i in 1:length(samples)) {
  sample <- samples[[i]]
  gsm_vec[i] <- ifelse(is.null(sample$GSM), NA, sample$GSM)
  title_vec[i] <- ifelse(is.null(sample$title), NA, sample$title)
  source_name_vec[i] <- ifelse(is.null(sample$source_name_ch1), NA, sample$source_name_ch1)
  
  # Handle characteristics_ch1
  if (is.list(sample$characteristics_ch1)) {
    characteristics_vec[i] <- paste(sample$characteristics_ch1, collapse = "; ")
  } else if (is.null(sample$characteristics_ch1)) {
    characteristics_vec[i] <- NA
  } else {
    characteristics_vec[i] <- sample$characteristics_ch1
  }
}

# Create data frame
pheno_df <- data.frame(
  GSM = gsm_vec,
  title = title_vec,
  source_name_ch1 = source_name_vec,
  characteristics_ch1 = characteristics_vec,
  stringsAsFactors = FALSE
)

# Remove NA rows
pheno_df <- pheno_df[!is.na(pheno_df$GSM), ]

cat(paste0("Extracted ", nrow(pheno_df), " valid samples\n"))

# Save output files
output_csv <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.csv")
output_rds <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.rds")

write.csv(pheno_df, output_csv, row.names = FALSE)
cat(paste0("Saved CSV: ", basename(output_csv), "\n"))

saveRDS(pheno_df, output_rds)
cat(paste0("Saved RDS: ", basename(output_rds), "\n"))

# Show first few rows
cat("\nFirst few rows of extracted metadata:\n")
print(head(pheno_df))

# Summary
cat("\n=== Extraction Summary ===\n")
cat(paste0("Total samples in soft file: ", length(samples), "\n"))
cat(paste0("Valid samples extracted: ", nrow(pheno_df), "\n"))
cat(paste0("Output files created: ", basename(output_csv), " and ", basename(output_rds), "\n"))

cat("\n=== Soft metadata extraction completed ===\n")
sink()

# Memory cleanup
rm(list = ls())
gc()
