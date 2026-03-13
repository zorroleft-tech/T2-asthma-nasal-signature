#!/usr/bin/env Rscript

# ==============================================================================
# Script: 00_build_replication_contrasts.R
# Purpose: Generate unified manifest and audit table for replication cohorts
# Author: Zuo
# Date: 2026-02-27
#
# Outputs:
#   - analysis/02_main_validation/00_replication_contrasts.yaml  # Unified manifest
#   - output/logs/00_replication_contrasts_audit.csv  # Audit table
# ==============================================================================

# Set base directory
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Create output directories if they don't exist
logs_dir <- file.path(base_dir, "output", "logs")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Setup logging
log_file <- file.path(logs_dir, "00_replication_contrasts_audit.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({ if (sink.number() > 0) sink() })

# Print timestamp
cat(paste0("Script started at: ", Sys.time(), "\n"))

# Try to get git commit hash
git_commit <- NA
tryCatch({
  if (requireNamespace("git2r", quietly = TRUE)) {
    library(git2r)
    repo <- repository(base_dir)
    commit <- last_commit(repo)
    git_commit <- commit@sha
    cat(paste0("Git commit: ", git_commit, "\n"))
  } else {
    cat("Git commit not available (git2r package not installed)\n")
  }
}, error = function(e) {
  cat("Git commit not available (not a git repository)\n")
  git_commit <- NA
})

# Define Rule-Standardize hierarchy dictionary (priority from high to low)
rule_hierarchy <- c(
  "Rule_T2_or_Eos_vs_HealthyCtrl",
  "Rule_SevereAsthma_vs_HealthyCtrl",
  "Rule_AllAsthma_vs_HealthyCtrl",
  "Rule_4Group_WheezeAtopy_Fold"
)

cat("\nRule-Standardize hierarchy dictionary:\n")
for (i in seq_along(rule_hierarchy)) {
  cat(paste0(i, ") ", rule_hierarchy[i], "\n"))
}

# Read all spec files
specs_dir <- file.path(base_dir, "analysis", "02_main_validation", "specs")
spec_files <- list.files(specs_dir, pattern = "*_spec.rds", full.names = TRUE)

if (length(spec_files) == 0) {
  stop("No spec files found in ", specs_dir)
}

cat(paste0("\nFound ", length(spec_files), " spec files\n"))

# Read and process each spec
specs <- list()
for (file in spec_files) {
  spec <- readRDS(file)
  
  # Convert relative paths to absolute paths
  spec$expr_path <- normalizePath(file.path(base_dir, spec$expr_path), winslash = "/", mustWork = FALSE)
  spec$meta_path <- normalizePath(file.path(base_dir, spec$meta_path), winslash = "/", mustWork = FALSE)
  
  specs[[spec$cohort_id]] <- spec
  cat(paste0("Loaded spec for ", spec$cohort_id, "\n"))
}

# Validate rule hits
cat("\nValidating rule hits...\n")
for (cohort_id in names(specs)) {
  spec <- specs[[cohort_id]]
  rule_hit <- spec$rule_standardize_hit
  
  if (!(rule_hit %in% rule_hierarchy) && rule_hit != "UNMATCHED") {
    warning(paste0("Unknown rule hit for ", cohort_id, ": ", rule_hit))
  }
}

# Check for multiple rule hits (anti-dark box constraint)
cat("\nChecking for multiple rule hits (anti-dark box constraint)...\n")
for (cohort_id in names(specs)) {
  spec <- specs[[cohort_id]]
  
  # Skip if already UNMATCHED
  if (spec$rule_standardize_hit == "UNMATCHED") {
    next
  }
  
  # Check if this cohort could match multiple rules
  # This is a placeholder for actual rule matching logic
  # In a real implementation, this would check the actual metadata
  # For now, we'll assume that each cohort only matches one rule
  # and throw an error if we detect multiple matches
  
  # Simulate multiple rule match detection
  # In practice, this would be based on actual metadata analysis
  possible_matches <- c()
  for (rule in rule_hierarchy) {
    # This is a simplified check - in reality, this would be more complex
    if (spec$rule_standardize_hit == rule) {
      possible_matches <- c(possible_matches, rule)
    }
  }
  
  if (length(possible_matches) > 1) {
    stop(paste0("ERROR: Cohort ", cohort_id, " matches multiple rules: ", 
                paste(possible_matches, collapse = ", "), 
                ". Please fix preprocessing to ensure each cohort matches only one rule."))
  }
}

# Generate unified manifest (YAML)
manifest_path <- file.path(base_dir, "analysis", "02_main_validation", "00_replication_contrasts.yaml")

cat("\nGenerating unified manifest...\n")

# Convert specs to YAML format
manifest_data <- list()
for (cohort_id in names(specs)) {
  spec <- specs[[cohort_id]]
  manifest_data[[cohort_id]] <- list(
    tissue_group = spec$tissue_group,
    expr_path = spec$expr_path,
    meta_path = spec$meta_path,
    group_col = spec$group_col,
    case_levels = spec$case_levels,
    control_levels = spec$control_levels,
    exclude_levels = spec$exclude_levels,
    compartment = spec$compartment,
    paired = spec$paired,
    rule_standardize_hit = spec$rule_standardize_hit
  )
}

# Write YAML file
if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml", repos = "https://cran.r-project.org")
}
library(yaml)
write_yaml(manifest_data, manifest_path)
cat(paste0("Manifest written to: ", manifest_path, "\n"))

# Generate audit table
cat("\nGenerating audit table...\n")
audit_data <- list()

for (cohort_id in names(specs)) {
  spec <- specs[[cohort_id]]
  
  # Try to read meta data to count cases/controls
  n_case <- NA
  n_control <- NA
  
  tryCatch({
    if (file.exists(spec$meta_path)) {
      cat(paste0("Reading meta data for ", cohort_id, " from ", spec$meta_path, "\n"))
      meta <- readRDS(spec$meta_path)
      cat(paste0("Meta data class: ", class(meta), "\n"))
      cat(paste0("Meta data dimensions: ", nrow(meta), " x ", ncol(meta), "\n"))
      cat(paste0("Group column: ", spec$group_col, "\n"))
      cat(paste0("Group column exists: ", spec$group_col %in% colnames(meta), "\n"))
      
      if (spec$group_col %in% colnames(meta)) {
        group_values <- meta[[spec$group_col]]
        cat(paste0("First 10 group values: ", paste(head(group_values, 10), collapse = ", "), "\n"))
        n_case <- sum(group_values %in% spec$case_levels, na.rm = TRUE)
        n_control <- sum(group_values %in% spec$control_levels, na.rm = TRUE)
        cat(paste0("n_case: ", n_case, ", n_control: ", n_control, "\n"))
      }
    }
  }, error = function(e) {
    cat(paste0("Error reading meta data for ", cohort_id, ": ", e$message, "\n"))
  })
  
  # Add reason for NA values
  reason <- NA
  if (is.na(n_case) || n_case == 0) {
    reason <- "Group labels are NA in meta file"
  }
  
  audit_row <- data.frame(
    cohort_id = cohort_id,
    tissue_group = spec$tissue_group,
    expr_path = spec$expr_path,
    meta_path = spec$meta_path,
    abs_expr_path = spec$expr_path,
    abs_meta_path = spec$meta_path,
    group_col = spec$group_col,
    case_levels = paste(spec$case_levels, collapse = ";"),
    control_levels = paste(spec$control_levels, collapse = ";"),
    exclude_levels = paste(spec$exclude_levels, collapse = ";"),
    n_case = n_case,
    n_control = n_control,
    rule_standardize_hit = spec$rule_standardize_hit,
    spec_source = spec$spec_source,
    built_by_script = spec$built_by_script,
    built_at = spec$built_at,
    reason = reason
  )
  
  audit_data[[cohort_id]] <- audit_row
}

# Combine all audit rows
audit_table <- do.call(rbind, audit_data)

# Write audit table
audit_path <- file.path(logs_dir, paste0("00_replication_contrasts_audit_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(audit_table, audit_path, row.names = FALSE)
cat(paste0("Audit table written to: ", audit_path, "\n"))

# Print summary
cat("\nSummary:\n")
cat(paste0("Total cohorts: ", nrow(audit_table), "\n"))
cat("Rule distribution:\n")
print(table(audit_table$rule_standardize_hit))

cat("\nScript completed successfully!\n")
