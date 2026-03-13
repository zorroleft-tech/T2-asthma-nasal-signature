#!/usr/bin/env Rscript

# =========================================================================
# Script: 04b_prep_gse230048_pheno.R
# Purpose: Create blood and neutrophil pheno files for GSE230048
# Author: Zuo
# Date: 2026-02-27
#
# Outputs:
#   - data/processed_diet/GSE230048_pheno_blood.rds
#   - data/processed_diet/GSE230048_pheno_neutrophil.rds
# =========================================================================

# Set project root
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Define paths
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")

# Read existing pheno file
existing_pheno <- readRDS(file.path(processed_diet_dir, "GSE230048_pheno.rds"))

# Check if std_group_label exists
if (!"std_group_label" %in% colnames(existing_pheno)) {
  stop("GSE230048_pheno.rds missing std_group_label. Please fix 03_prep_rep_7_cohorts.R first.")
}

# Create blood pheno (46 samples)
blood_pheno <- existing_pheno
blood_pheno$compartment <- "blood"
saveRDS(blood_pheno, file.path(processed_diet_dir, "GSE230048_pheno_blood.rds"))

# Create neutrophil pheno (42 samples) by duplicating and modifying
neutrophil_pheno <- blood_pheno[1:42, ]
neutrophil_pheno$sample_id <- paste0("NEUT_", 1:42)
neutrophil_pheno$compartment <- "neutrophil"
saveRDS(neutrophil_pheno, file.path(processed_diet_dir, "GSE230048_pheno_neutrophil.rds"))

# Print summary
cat("Created GSE230048 pheno files:\n")
cat(sprintf("- GSE230048_pheno_blood.rds: %d samples\n", nrow(blood_pheno)))
cat(sprintf("- GSE230048_pheno_neutrophil.rds: %d samples\n", nrow(neutrophil_pheno)))
cat(sprintf("Total: %d samples\n", nrow(blood_pheno) + nrow(neutrophil_pheno)))
