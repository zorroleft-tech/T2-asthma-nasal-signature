#!/usr/bin/env Rscript

# ==============================================================================
# Script: 00_build_replication_specs_all.R
# Purpose: Generate spec files for all replication-only cohorts
# Author: Zuo
# Date: 2026-02-27
#
# Outputs:
#   - analysis/02_main_validation/specs/<GSE>_spec.rds  # Spec files for each cohort
# ==============================================================================

# Set base directory
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Create specs directory if it doesn't exist
specs_dir <- file.path(base_dir, "analysis", "02_main_validation", "specs")
dir.create(specs_dir, recursive = TRUE, showWarnings = FALSE)

# Function to create spec for a cohort
create_cohort_spec <- function(gse_id, tissue_group, expr_path, meta_path, 
                              group_col, case_levels, control_levels, 
                              exclude_levels = character(0), compartment = NULL, 
                              paired = FALSE, rule_standardize_hit, spec_source) {
  spec <- list(
    cohort_id = gse_id,
    tissue_group = tissue_group,
    expr_path = expr_path,
    meta_path = meta_path,
    group_col = group_col,
    case_levels = case_levels,
    control_levels = control_levels,
    exclude_levels = exclude_levels,
    compartment = compartment,
    paired = paired,
    rule_standardize_hit = rule_standardize_hit,
    spec_source = spec_source,
    built_by_script = "00_build_replication_specs_all.R",
    built_at = Sys.time()
  )
  
  # Save spec to RDS file
  spec_file <- file.path(specs_dir, paste0(gse_id, "_spec.rds"))
  saveRDS(spec, spec_file)
  
  cat(paste0("Created spec file for ", gse_id, " at ", spec_file, "\n"))
  
  return(spec)
}

# Define replication-only cohorts
replication_cohorts <- c(
  "GSE103166",
  "GSE115770",
  "GSE118761",
  "GSE43696",
  "GSE123750",
  "GSE40888",
  "GSE230048"
)

# Create spec files for each cohort

# GSE103166
create_cohort_spec(
  gse_id = "GSE103166",
  tissue_group = "Airway",
  expr_path = "data/processed_diet/GSE103166_diet_expr.rds",
  meta_path = "data/processed_diet/GSE103166_pheno.rds",
  group_col = "std_group_label",
  case_levels = c("Asthma"),
  control_levels = c("Healthy"),
  rule_standardize_hit = "Rule_AllAsthma_vs_HealthyCtrl",
  spec_source = "GEO metadata + preprocessing script"
)

# GSE115770
create_cohort_spec(
  gse_id = "GSE115770",
  tissue_group = "Blood_PBMC",
  expr_path = "data/processed_diet/GSE115770_diet_expr.rds",
  meta_path = "data/processed_diet/GSE115770_pheno.rds",
  group_col = "std_group_label",
  case_levels = c("Severe Asthma"),
  control_levels = c("Healthy"),
  rule_standardize_hit = "Rule_SevereAsthma_vs_HealthyCtrl",
  spec_source = "GEO metadata + preprocessing script"
)

# GSE118761
create_cohort_spec(
  gse_id = "GSE118761",
  tissue_group = "Airway",
  expr_path = "data/processed_diet/GSE118761_diet_expr.rds",
  meta_path = "data/processed_diet/GSE118761_pheno.rds",
  group_col = "std_group_label",
  case_levels = c("Asthma"),
  control_levels = c("Healthy"),
  paired = TRUE,
  rule_standardize_hit = "Rule_AllAsthma_vs_HealthyCtrl",
  spec_source = "GEO metadata + preprocessing script"
)

# GSE43696
create_cohort_spec(
  gse_id = "GSE43696",
  tissue_group = "Blood_PBMC",
  expr_path = "data/processed_diet/GSE43696_diet_expr.rds",
  meta_path = "data/processed_diet/GSE43696_pheno.rds",
  group_col = "std_group_label",
  case_levels = c("Asthma"),
  control_levels = c("Healthy"),
  rule_standardize_hit = "Rule_AllAsthma_vs_HealthyCtrl",
  spec_source = "GEO metadata + preprocessing script"
)

# GSE123750
create_cohort_spec(
  gse_id = "GSE123750",
  tissue_group = "Nasal",
  expr_path = "data/processed_diet/GSE123750_diet_expr.rds",
  meta_path = "data/processed_diet/GSE123750_pheno.rds",
  group_col = "std_group_label",
  case_levels = c("Asthma"),
  control_levels = c("Healthy"),
  rule_standardize_hit = "Rule_AllAsthma_vs_HealthyCtrl",
  spec_source = "GEO metadata + preprocessing script"
)

# GSE40888
create_cohort_spec(
  gse_id = "GSE40888",
  tissue_group = "Blood_PBMC",
  expr_path = "data/processed_diet/GSE40888_diet_expr.rds",
  meta_path = "data/processed_diet/GSE40888_pheno.rds",
  group_col = "std_group_label",
  case_levels = c("Atopic wheeze"),
  control_levels = c("Healthy non-wheeze"),
  rule_standardize_hit = "Rule_4Group_WheezeAtopy_Fold",
  spec_source = "GEO metadata + preprocessing script (folded from 4 original groups)"
)

# GSE230048
create_cohort_spec(
  gse_id = "GSE230048",
  tissue_group = "Blood_PBMC",
  expr_path = "data/processed_diet/GSE230048_diet_expr.rds",
  meta_path = "data/processed_diet/GSE230048_pheno_blood.rds",
  group_col = "std_group_label",
  case_levels = c("Asthma"),
  control_levels = c("Healthy"),
  compartment = "blood",
  rule_standardize_hit = "Rule_AllAsthma_vs_HealthyCtrl",
  spec_source = "GEO metadata + preprocessing script"
)

cat("\nAll spec files created successfully!\n")
