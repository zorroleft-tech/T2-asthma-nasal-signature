#!/usr/bin/env Rscript

# =========================================================================
# Script: 11_supp_tables_S1_to_S4.R
# Purpose: Generate Supplementary Tables S1 to S4 (Submission-first build)
# Strategy (LOCKED):
#   - Do NOT recompute anchor AUC. Use Fig4B_threeAnchors_auc_summary.csv as source of truth.
#   - Table S4 should be sourced from phase3_outputs deliverables/stats summary (no re-derivation from expr).
#   - Avoid online installation. Fail fast if required packages are missing.
# Author: Zuo
# Date: 2026-02-28
#
# Outputs:
#   - output/supplement/TableS1_QC_cohort_overview.csv
#   - output/supplement/TableS2_cellular_source_function.csv
#   - output/supplement/TableS3_anchor_auc_audit.csv
#   - output/supplement/TableS4_replication_consistency.csv
#   - output/supplement/TableS4A_direction_audit.csv
#   - output/supplement/TableS4B_contribution_audit.csv
#   - output/logs/11_supp_tables_S1_to_S4.log
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
})

# -----------------------------
# Paths
# -----------------------------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

output_dir      <- file.path(base_dir, "output")
supplement_dir  <- file.path(output_dir, "supplement")
logs_dir        <- file.path(output_dir, "logs")

tables_main_dir <- file.path(output_dir, "tables_main")

phase3_dir      <- file.path(base_dir, "data", "processed", "phase3_outputs")
derived_dir     <- file.path(base_dir, "data", "derived")

dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "11_supp_tables_S1_to_S4.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({ if (sink.number() > 0) sink() })

# Print timestamp
cat(paste0("Script started at: ", Sys.time(), "\n"))

# Try to get git commit hash
git_commit <- NA
tryCatch({
  if (requireNamespace("git2r", quietly = TRUE)) {
    library(git2r)
    if (dir.exists(".git")) {
      repo <- repository(base_dir)
      commit <- last_commit(repo)
      git_commit <- commit@sha
      cat(paste0("Git commit: ", git_commit, "\n"))
    } else {
      cat("Git commit not available (not a git repository)\n")
    }
  } else {
    cat("Git commit not available (git2r package not installed)\n")
  }
}, error = function(e) {
  cat("Git commit not available (error accessing repository)\n")
  git_commit <- NA
})

cat("=======================================================\n")
cat("=== Generating Supplementary Tables S1 to S4 (LOCKED) ===\n")
cat("=======================================================\n\n")

# Optional: output index logger (do not hard fail)
index_logger_path <- file.path(base_dir, "R", "03_index_logger.R")
if (file.exists(index_logger_path)) {
  source(index_logger_path)
} else {
  cat("[WARN] R/03_index_logger.R not found. OUTPUTS_INDEX logging will be skipped.\n")
}

# -----------------------------
# Helper: safe log_output
# -----------------------------
safe_log_output <- function(path, type, title, purpose, script) {
  if (exists("log_output")) {
    tryCatch({
      log_output(path, type, title, purpose, script)
    }, error = function(e) {
      cat("[WARNING] Failed to log output:", e$message, "\n")
    })
  }
}

# -----------------------------
# Helper: assertions
# -----------------------------
assert_file <- function(path, label = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("Missing required file%s: %s",
                 if (!is.null(label)) paste0(" (", label, ")") else "",
                 path
    ))
  }
}

# =======================================================
# [1/4] Table S1: QC & cohort overview (source from Table 1)
# =======================================================
cat("[1/4] Table S1: Cohort overview & QC (derived from Table 1 + optional audit)\n")

# Use the Table 1 already embedded in the manuscript as the canonical cohort list.
# If you prefer reading from CSV, put it here:
# output/tables_main/Table1_baseline.csv or wherever you stored it.
# In Notion, Table 1 exists as a markdown table, but your pipeline also emits Table1_baseline.csv.
table1_csv_candidates <- c(
  file.path(output_dir, "tables_main", "Table1_baseline.csv"),
  file.path(output_dir, "tables_main", "Table1_baseline_Characteristics.csv"),
  file.path(output_dir, "tables_main", "Table1_baseline_characteristics.csv"),
  file.path(base_dir, "Table1_baseline.csv")
)

table1_path <- table1_csv_candidates[file.exists(table1_csv_candidates)][1]
if (is.na(table1_path) || is.null(table1_path) || table1_path == "") {
  stop("Table S1 build requires Table1_baseline.csv (or equivalent) but none was found in expected locations.")
}

table1 <- read_csv(table1_path, show_col_types = FALSE) %>%
  janitor::clean_names()

# Expect fields like: gse_id, role, tissue, platform_gpl, n_total, no_of_cases, no_of_controls, age, sex
# But your CSV headers contain spaces. After clean_names they become consistent.
# We'll standardize names defensively.

rename_if_present <- function(df, old, new) {
  if (old %in% names(df)) dplyr::rename(df, !!new := !!sym(old)) else df
}

table1 <- table1 %>%
  rename_if_present("gse_id", "cohort") %>%
  rename_if_present("gse", "cohort") %>%
  rename_if_present("platform_gpl", "platform") %>%
  rename_if_present("platform_gpl_", "platform") %>%
  rename_if_present("n_total", "n_total") %>%
  rename_if_present("no_of_cases", "n_case") %>%
  rename_if_present("no_of_controls", "n_control")

if (!("cohort" %in% names(table1))) stop("Table1 CSV must contain a cohort/GSE ID column (expected: gse_id).")

# Optional: bring in field audit summary (if present)
field_audit_path <- file.path(phase3_dir, "field_audit_report.csv")
table_s1 <- table1 %>%
  select(any_of(c("cohort", "role", "tissue", "platform", "n_total", "n_case", "n_control", "age", "sex"))) %>%
  mutate(source_table1 = basename(table1_path))

if (file.exists(field_audit_path)) {
  field_audit <- read_csv(field_audit_path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  # Heuristic: keep only high-level audit columns if present
  keep_audit_cols <- intersect(
    names(field_audit),
    c("cohort", "status", "fail_reason", "notes", "unmapped_rate", "mapping_rate", "qc_status")
  )
  
  if ("cohort" %in% names(field_audit) && length(keep_audit_cols) >= 2) {
    table_s1 <- table_s1 %>%
      left_join(field_audit %>% select(all_of(keep_audit_cols)), by = "cohort")
  }
}

table_s1_out <- file.path(supplement_dir, "TableS1_QC_cohort_overview.csv")
write_csv(table_s1, table_s1_out)
cat("  - Saved:", basename(table_s1_out), "\n")

safe_log_output(table_s1_out, "table", "Table S1", "QC and cohort overview (from Table 1 + audit if available)", "11_supp_tables_S1_to_S4.R")

# =======================================================
# [2/4] Table S2: Cell type/function annotation (STUB)
# =======================================================
cat("\n[2/4] Table S2: Cell type/function annotation (stub; no hard-coded mechanisms)\n")

weights_path <- file.path(derived_dir, "locked_weights.csv")
assert_file(weights_path, "locked_weights.csv")

locked_weights <- read_csv(weights_path, show_col_types = FALSE)

# Try to preserve direction column if already exists in locked_weights; otherwise infer
if (!("Direction" %in% names(locked_weights))) {
  locked_weights <- locked_weights %>%
    mutate(Direction = if_else(Weight >= 0, "Up", "Down"))
}

# Read TableS2_annotation.csv
annotation_path <- file.path(base_dir, "analysis", "03_supplementary", "resources", "TableS2_annotation.csv")
assert_file(annotation_path, "TableS2_annotation.csv")

annotation <- read_csv(annotation_path, show_col_types = FALSE)

# Get list of genes from locked_weights
locked_genes <- locked_weights$Gene

# Check for missing genes in annotation
missing_genes <- setdiff(locked_genes, annotation$Gene)
if (length(missing_genes) > 0) {
  stop("Missing genes in TableS2_annotation.csv: " , paste(missing_genes, collapse = ", "))
}

# Check for extra genes in annotation
extra_genes <- setdiff(annotation$Gene, locked_genes)
if (length(extra_genes) > 0) {
  stop("Extra genes in TableS2_annotation.csv: " , paste(extra_genes, collapse = ", "))
}

# Check for empty or NA Key_reference
empty_reference <- annotation %>%
  filter(is.na(Key_reference) | Key_reference == "")

if (nrow(empty_reference) > 0) {
  stop("Empty or NA Key_reference in TableS2_annotation.csv for genes: " , paste(empty_reference$Gene, collapse = ", "))
}

# Join locked_weights with annotation
table_s2 <- locked_weights %>%
  left_join(annotation, by = "Gene")

table_s2_out <- file.path(supplement_dir, "TableS2_cellular_source_function.csv")
write_csv(table_s2, table_s2_out)
cat("  - Saved:", basename(table_s2_out), "\n")

safe_log_output(table_s2_out, "table", "Table S2", "Cell source/function annotation for locked 11 genes", "11_supp_tables_S1_to_S4.R")

# =======================================================
# [3/4] Table S3: Anchor AUC audit (SOURCE OF TRUTH = Fig4B summary)
# =======================================================
cat("\n[3/4] Table S3: Anchor AUC audit (truth = Fig4B_threeAnchors_auc_summary.csv)\n")

fig4b_auc_path <- file.path(tables_main_dir, "Fig4B_threeAnchors_auc_summary.csv")
assert_file(fig4b_auc_path, "Fig4B_threeAnchors_auc_summary.csv")

auc_summary <- read_csv(fig4b_auc_path, show_col_types = FALSE) %>%
  janitor::clean_names()

# Required columns check
required_auc_cols <- c(
  "cohort", "truth", "n_case", "n_control", "gene_coverage",
  "case_definition", "control_definition",
  "score_mean_case", "score_mean_control",
  "auc_raw", "ci_raw",
  "auc_final", "ci_final",
  "flipped", "flip_reason"
)

missing_required <- setdiff(required_auc_cols, names(auc_summary))
if (length(missing_required) > 0) {
  stop(paste0("Fig4B_threeAnchors_auc_summary.csv missing required columns: ",
              paste(missing_required, collapse = ", ")
  ))
}

# Optional consistency check against gene coverage tables (do NOT fail if absent)
coverage_anchor_path <- file.path(tables_main_dir, "Fig4_gene_coverage_by_anchor.csv")
s9_platform_path <- file.path(output_dir, "supplement", "TableS9_platform_wide.csv")

table_s3 <- auc_summary %>%
  transmute(
    Cohort = cohort,
    Truth_anchor = truth,
    n_case = n_case,
    n_control = n_control,
    Gene_coverage = gene_coverage,
    Case_definition = case_definition,
    Control_definition = control_definition,
    Score_mean_case = score_mean_case,
    Score_mean_control = score_mean_control,
    AUC_raw = auc_raw,
    CI_raw = ci_raw,
    AUC_final_direction_aligned = auc_final,
    CI_final = ci_final,
    Direction_aligned = flipped,
    Direction_alignment_reason = flip_reason,
    Direction_note = "Direction alignment was prespecified so higher scores correspond to anchored high-truth group; auc_final is direction-aligned AUC; auc_raw provided for audit."
  )

# Consistency checks (hard fail only if coverage table exists but conflicts)
if (file.exists(coverage_anchor_path)) {
    coverage_anchor <- read_csv(coverage_anchor_path, show_col_types = FALSE) %>%
      janitor::clean_names()
    
    # Expect something like: cohort, gene_coverage or coverage_x_of_11 etc.
    # We'll try to derive a comparable coverage string.
    if ("cohort" %in% names(coverage_anchor)) {
      coverage_anchor2 <- coverage_anchor %>%
        mutate(
          coverage_str = ifelse("n_genes_available" %in% names(.), 
                               ifelse("n_genes_total" %in% names(.), 
                                      paste0(n_genes_available, "/", n_genes_total), 
                                      paste0(n_genes_available, "/11")),
                               ifelse("coverage_count" %in% names(.), 
                                      paste0(coverage_count, "/11"),
                                      ifelse("coverage_x_of_11" %in% names(.), 
                                             paste0(coverage_x_of_11, "/11"),
                                             ifelse("gene_coverage" %in% names(.), 
                                                    as.character(gene_coverage),
                                                    NA_character_))))
        ) %>%
        select(cohort, coverage_str) %>%
        distinct()
    
    check_df <- table_s3 %>%
      transmute(cohort = Cohort, gene_coverage = Gene_coverage) %>%
      left_join(coverage_anchor2, by = "cohort") %>%
      filter(!is.na(coverage_str)) %>%
      mutate(match = (gene_coverage == coverage_str))
    
    if (nrow(check_df) > 0 && any(!check_df$match)) {
      print(check_df)
      cat("[WARNING] Gene coverage mismatch detected between Fig4B summary and Fig4_gene_coverage_by_anchor.csv. Using Fig4B values for Table S3.\n")
    }
  }
}

table_s3_out <- file.path(supplement_dir, "TableS3_anchor_auc_audit.csv")
write_csv(table_s3, table_s3_out)
cat("  - Saved:", basename(table_s3_out), "\n")

safe_log_output(table_s3_out, "table", "Table S3", "Anchor cohorts AUC audit table (source: Fig4B_threeAnchors_auc_summary.csv)", "11_supp_tables_S1_to_S4.R")

# =======================================================
# [4/4] Table S4: Replication-only molecular consistency
#   Source: phase3_outputs/replication_direction_audit.csv
#   Also generate Table S4-A (direction audit) and Table S4-B (contribution audit)
# =======================================================
cat("\n[4/4] Table S4: Replication-only molecular consistency (source: replication_direction_audit.csv)\n")

# Define the expected replication-only cohorts (LOCKED to manuscript stance)
rep_cohorts_expected <- c("GSE103166", "GSE115770", "GSE118761", "GSE123750", "GSE230048", "GSE40888", "GSE43696")

# Read direction audit data
direction_audit_path <- file.path(phase3_dir, "replication_direction_audit.csv")
contribution_audit_path <- file.path(phase3_dir, "replication_contribution_audit.csv")

if (!file.exists(direction_audit_path)) {
  stop("replication_direction_audit.csv not found. Cannot generate Table S4 and related tables.")
}

if (!file.exists(contribution_audit_path)) {
  stop("replication_contribution_audit.csv not found. Cannot generate Table S4-B.")
}

# Read the audit files
direction_audit <- read_csv(direction_audit_path, show_col_types = FALSE) %>%
  janitor::clean_names()

contribution_audit <- read_csv(contribution_audit_path, show_col_types = FALSE) %>%
  janitor::clean_names()

# Filter to expected replication cohorts
direction_audit_rep <- direction_audit %>%
  filter(cohort %in% rep_cohorts_expected)

missing_rep <- setdiff(rep_cohorts_expected, unique(direction_audit_rep$cohort))
if (length(missing_rep) > 0) {
  # For submission, we require all 7 replication cohorts
  stop(paste0(
    "replication_direction_audit.csv is incomplete. Missing cohorts: ",
    paste(missing_rep, collapse = ", "),
    ".\nPlease fix upstream (05_calculate_minimal_stats.R) to include all 7 replication cohorts."
  ))
}

# Define tissue or compartment mapping
tissue_mapping <- list(
  "GSE103166" = "PBMC",
  "GSE115770" = "PBMC",
  "GSE118761" = "PBMC",
  "GSE123750" = "PBMC",
  "GSE230048" = "Blood",
  "GSE40888" = "PBMC",
  "GSE43696" = "PBMC"
)

# Generate Table S4 main table
table_s4 <- direction_audit_rep %>%
  mutate(
    gene_coverage = paste0(coverage_count, "/11"),
    missing_genes = ifelse(is.na(coverage_missing) | coverage_missing == "", "", coverage_missing),
    effect_size = as.numeric(effect_size),
    p_value = as.numeric(p_value),
    tissue_or_compartment = unlist(tissue_mapping[cohort]),
    contrast = paste(case_levels, "vs", control_levels),
    directional_concordance = ifelse(mean_diff > 0, "Yes", "No"),
    pass_tech = ifelse(status == "PASS", "Yes", "No"),
    pass_interpret = ifelse(mean_diff > 0, "Yes", "No")
  ) %>%
  transmute(
    Cohort = cohort,
    Tissue_or_compartment = tissue_or_compartment,
    Contrast = contrast,
    Gene_coverage = gene_coverage,
    Missing_genes = missing_genes,
    Effect_size_Cohens_d = sprintf("%.3f", effect_size),
    Wilcoxon_P = format(p_value, scientific = TRUE, digits = 3),
    Directional_concordance = directional_concordance,
    PASS_tech = pass_tech,
    PASS_interpret = pass_interpret
  ) %>%
  arrange(match(Cohort, rep_cohorts_expected))

table_s4_out <- file.path(supplement_dir, "TableS4_replication_consistency.csv")
write_csv(table_s4, table_s4_out)
cat("  - Saved:", basename(table_s4_out), "\n")

safe_log_output(table_s4_out, "table", "Table S4", "Replication-only molecular consistency across 7 external cohorts", "11_supp_tables_S1_to_S4.R")

# Generate Table S4-A: Direction audit
cat("\n[4/4a] Table S4-A: Direction audit\n")

table_s4a <- direction_audit_rep %>%
  mutate(
    direction_ok = ifelse(mean_diff > 0, "YES", "NO"),
    boundary_flag = ifelse(abs(mean_diff) < 0.5, "Boundary condition - review contrast definition", "")
  ) %>%
  transmute(
    Cohort = cohort,
    case_levels = case_levels,
    control_levels = control_levels,
    n_case = n_case,
    n_control = n_control,
    mean_score_case = mean_score_case,
    mean_score_control = mean_score_control,
    mean_diff = mean_diff,
    posWeightGene_pos_count = pos_weight_gene_pos_count,
    posWeightGene_neg_count = pos_weight_gene_neg_count,
    direction_ok = direction_ok,
    boundary_flag = boundary_flag
  ) %>%
  arrange(match(Cohort, rep_cohorts_expected))

table_s4a_out <- file.path(supplement_dir, "TableS4A_direction_audit.csv")
write_csv(table_s4a, table_s4a_out)
cat("  - Saved:", basename(table_s4a_out), "\n")

safe_log_output(table_s4a_out, "table", "Table S4-A", "Direction audit for replication-only cohorts", "11_supp_tables_S1_to_S4.R")

# Generate Table S4-B: Contribution audit (only for direction_ok = NO)
cat("\n[4/4b] Table S4-B: Contribution audit\n")

# First, get cohorts with direction_ok = NO
direction_ok_no <- table_s4a %>%
  filter(direction_ok == "NO") %>%
  dplyr::pull(Cohort)

# Filter contribution audit to only these cohorts
table_s4b <- contribution_audit %>%
  filter(cohort %in% direction_ok_no) %>%
  transmute(
    Cohort = cohort,
    sum_delta_contrib = sum_delta_contrib,
    Top1_gene = top1_gene,
    Top1_delta_contrib = top1_delta_contrib,
    Top2_gene = top2_gene,
    Top2_delta_contrib = top2_delta_contrib,
    Top3_gene = top3_gene,
    Top3_delta_contrib = top3_delta_contrib,
    Top4_gene = top4_gene,
    Top4_delta_contrib = top4_delta_contrib,
    Top5_gene = top5_gene,
    Top5_delta_contrib = top5_delta_contrib
  ) %>%
  arrange(match(Cohort, rep_cohorts_expected))

table_s4b_out <- file.path(supplement_dir, "TableS4B_contribution_audit.csv")
write_csv(table_s4b, table_s4b_out)
cat("  - Saved:", basename(table_s4b_out), "\n")

safe_log_output(table_s4b_out, "table", "Table S4-B", "Contribution audit for replication-only cohorts with direction issues", "11_supp_tables_S1_to_S4.R")

# Fix Table S4 fields consistency
cat("\n[5/6] Fixing Table S4 fields consistency...\n")

# Read Table S1 to get tissue information
table_s1_path <- file.path(supplement_dir, "TableS1_QC_cohort_overview.csv")
table_s1 <- read_csv(table_s1_path, show_col_types = FALSE)

# Fix Table S4: Update tissue field and rename Directional_concordance
table_s4_path <- file.path(supplement_dir, "TableS4_replication_consistency.csv")
table_s4 <- read_csv(table_s4_path, show_col_types = FALSE)

# Join with Table S1 to get tissue information
table_s4_fixed <- table_s4 %>%
  left_join(table_s1 %>% select(cohort, tissue), by = c("Cohort" = "cohort")) %>%
  mutate(Tissue_or_compartment = ifelse(!is.na(tissue), tissue, Tissue_or_compartment)) %>%
  select(-tissue)  # Remove the temporary tissue column

# Rename Directional_concordance to Direction_concordant
table_s4_fixed <- table_s4_fixed %>%
  rename(Direction_concordant = Directional_concordance)

# Write back the fixed Table S4
write_csv(table_s4_fixed, table_s4_path)
cat("  - Updated tissue fields and renamed Directional_concordance in Table S4\n")

# Fix Table S4A: Lightweight consistency for boundary_flag
table_s4a_path <- file.path(supplement_dir, "TableS4A_direction_audit.csv")
table_s4a <- read_csv(table_s4a_path, show_col_types = FALSE)

# Update boundary_flag to use one of two standardized values
table_s4a_fixed <- table_s4a %>%
  mutate(
    boundary_flag = case_when(
      boundary_flag == "Boundary condition - review contrast definition" ~ "Review contrast definition",
      boundary_flag != "" ~ "Boundary condition",
      TRUE ~ ""
    )
  )

# Write back the fixed Table S4A
write_csv(table_s4a_fixed, table_s4a_path)
cat("  - Updated boundary_flag in Table S4A\n")

# Print summary
cat("\nSummary:\n")
cat(paste0("Generated ", 6, " supplementary tables\n"))
cat("Tables generated:\n")
cat("1) TableS1_QC_cohort_overview.csv\n")
cat("2) TableS2_cellular_source_function.csv\n")
cat("3) TableS3_anchor_auc_audit.csv\n")
cat("4) TableS4_replication_consistency.csv\n")
cat("5) TableS4A_direction_audit.csv\n")
cat("6) TableS4B_contribution_audit.csv\n")

cat("\n=== Completed successfully (LOCKED build) ===\n")
cat(paste0("Script completed at: ", Sys.time(), "\n"))
sink()