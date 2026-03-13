#!/usr/bin/env Rscript

# ==============================================================================
# Script: 07a_fig1_discovery_common_pathways.R
# Purpose: Generate common pathways across replication cohorts for Figure 1C
# Author: Zuo
# Date: 2026-03-02
#
# Inputs:
#   - Various expression and phenotype files for replication cohorts (processed automatically)
#
# Outputs:
#   - output/tables_main/Fig1C_common_pathways_replication7.csv  # Common pathways across 7 replication cohorts
#   - output/tables_main/Fig1C_GSEA_<cohort_id>.csv  # GSEA results for each replication cohort
#   - output/logs/07a_fig1_discovery_common_pathways.log  # Log file
# ==============================================================================

# Load necessary packages
library(tidyverse)
library(limma)
library(fgsea)
library(rlang)

# ==============================================================================
# 1. Initialize base directory and setup logging
# ==============================================================================

# Set base directory
base_dir <- normalizePath(getwd(), winslash="/", mustWork=TRUE)

# Verify project structure
required_dirs <- c("data", "analysis", "data_preparation")
for (dir in required_dirs) {
  if (!dir.exists(file.path(base_dir, dir))) {
    stop(paste0("Required directory not found: ", file.path(base_dir, dir)))
  }
}

# Create output directories
output_dir <- file.path(base_dir, "output")
logs_dir <- file.path(output_dir, "logs")

for (dir in c(logs_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("Created directory: ", dir, "\n"))
  }
}

# Setup logging
log_file <- file.path(logs_dir, "07a_fig1_discovery_common_pathways.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({
  if (sink.number() > 0) {
    sink()
  }
})

cat("=== Starting common pathways generation ===\n")
cat(paste0("Base directory: ", base_dir, "\n"))

# ==============================================================================
# 2. Load theme and utility functions
# ==============================================================================

# Add file existence check to avoid source failure
theme_script <- file.path(base_dir, "R", "00_theme_setup.R")
logger_script <- file.path(base_dir, "R", "03_index_logger.R")

if (file.exists(theme_script)) {
  source(theme_script)
} else {
  warning(paste0("Theme script not found: ", theme_script, " - using default theme"))
  set_plot_theme <- function() { theme_set(theme_bw()) }
}

if (file.exists(logger_script)) {
  source(logger_script)
} else {
  warning(paste0("Logger script not found: ", logger_script, " - log_output function disabled"))
  log_output <- function(...) { cat("Log output disabled - script missing\n") }
}

set_plot_theme()

# ==============================================================================
# 3. Replication Cohorts GSEA Analysis
# ==============================================================================

cat("\n=== Generating Replication Cohorts GSEA Analysis ===\n")

# Define replication-only cohorts
replication_only_cohorts <- c("GSE103166", "GSE115770", "GSE118761", "GSE123750", "GSE230048", "GSE40888", "GSE43696")
cat("Replication-only cohorts: ", paste(replication_only_cohorts, collapse = ", "), "\n")

# Function to resolve expression and phenotype files for a cohort
resolve_replication_expr_pheno <- function(cohort_id, base_dir) {
  # Define potential file paths
  # First check processed_full directory (FULL version first, then RAW)
  processed_full_dir <- file.path(base_dir, "data", "processed_full")
  cohort_lower <- tolower(cohort_id)
  
  # List all files in processed_full
  all_files <- list.files(processed_full_dir, full.names = TRUE)
  
  # Find expression files in processed_full
  expr_files_full <- character(0)
  expr_files_raw <- character(0)
  pheno_files_full <- character(0)
  
  for (file in all_files) {
    if (grepl(paste0("expr_", cohort_lower), file) && grepl("__FULL.rds", file)) {
      expr_files_full <- c(expr_files_full, file)
    } else if (grepl(paste0("expr_", cohort_lower), file) && grepl("__RAW.rds", file)) {
      expr_files_raw <- c(expr_files_raw, file)
    } else if (grepl(paste0("pheno_", cohort_id), file) && grepl(".rds", file)) {
      pheno_files_full <- c(pheno_files_full, file)
    }
  }
  
  # Define fallback paths
  expr_paths <- c(
    if (length(expr_files_full) > 0) expr_files_full[1] else NULL,
    if (length(expr_files_raw) > 0) expr_files_raw[1] else NULL,
    file.path(base_dir, "data", "processed_diet", paste0(cohort_id, "_diet_expr.rds"))
  )
  
  pheno_paths <- c(
    if (length(pheno_files_full) > 0) pheno_files_full[1] else NULL,
    file.path(base_dir, "data", "processed", paste0(cohort_id, "_pheno.rds")),
    file.path(base_dir, "data", "processed_diet", paste0(cohort_id, "_pheno.rds"))
  )
  
  # Find existing expression file
  expr_file <- NULL
  for (path in expr_paths) {
    if (!is.null(path) && file.exists(path)) {
      expr_file <- path
      break
    }
  }
  
  # Find existing phenotype file
  pheno_file <- NULL
  for (path in pheno_paths) {
    if (!is.null(path) && file.exists(path)) {
      pheno_file <- path
      break
    }
  }
  
  return(list(expr_file = expr_file, pheno_file = pheno_file))
}

# Function to run KEGG GSEA for a single cohort
run_kegg_fgsea_for_cohort <- function(cohort_id, base_dir) {
  cat(paste0("\nRunning GSEA for cohort: ", cohort_id, "\n"))
  
  # Load required packages for gene ID mapping
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
  })
  
  # Resolve input files
  input_files <- resolve_replication_expr_pheno(cohort_id, base_dir)
  expr_file <- input_files$expr_file
  pheno_file <- input_files$pheno_file
  
  if (is.null(expr_file) || is.null(pheno_file)) {
    cat(paste0("  ERROR: Missing input files for cohort ", cohort_id, "\n"))
    if (is.null(expr_file)) cat("  - Expression file not found\n")
    if (is.null(pheno_file)) cat("  - Phenotype file not found\n")
    return(NULL)
  }
  
  cat(paste0("  Using expression file: ", expr_file, "\n"))
  cat(paste0("  Using phenotype file: ", pheno_file, "\n"))
  
  # Load expression data
  expr_data <- readRDS(expr_file)
  cat(paste0("  Expression data dimensions: ", nrow(expr_data), " genes x ", ncol(expr_data), " samples\n"))
  
  # Load phenotype data
  pheno_data <- readRDS(pheno_file)
  
  # Check for required columns in phenotype data
  if (!all(c("SampleID", "group") %in% colnames(pheno_data))) {
    # Try alternative sample ID column names
    if ("std_sample_id" %in% colnames(pheno_data)) {
      colnames(pheno_data)[colnames(pheno_data) == "std_sample_id"] <- "SampleID"
    } else if ("sample_id" %in% colnames(pheno_data)) {
      colnames(pheno_data)[colnames(pheno_data) == "sample_id"] <- "SampleID"
    } else if ("geo_accession" %in% colnames(pheno_data)) {
      colnames(pheno_data)[colnames(pheno_data) == "geo_accession"] <- "SampleID"
    }
    
    # Try alternative group column names
    if (!"group" %in% colnames(pheno_data)) {
      group_cols <- c("Group", "condition", "phenotype", "disease_status", "std_group_label", "asthma_status", "T2_Status")
      for (col in group_cols) {
        if (col %in% colnames(pheno_data)) {
          colnames(pheno_data)[colnames(pheno_data) == col] <- "group"
          break
        }
      }
    }
  }
  
  if (!all(c("SampleID", "group") %in% colnames(pheno_data))) {
    cat("  ERROR: Phenotype file missing required columns: SampleID and group\n")
    return(NULL)
  }
  
  # Ensure sample IDs match
  common_samples <- intersect(colnames(expr_data), pheno_data$SampleID)
  cat(paste0("  Common samples: ", length(common_samples), "\n"))
  
  if (length(common_samples) == 0) {
    cat("  ERROR: No common samples between expression and phenotype data\n")
    return(NULL)
  }
  
  # Subset data
  expr_data <- expr_data[, common_samples]
  pheno_data <- pheno_data %>% filter(SampleID %in% common_samples)
  
  # Ensure order matches
  expr_data <- expr_data[, pheno_data$SampleID]
  
  # Determine case and control labels
  groups <- unique(pheno_data$group)
  if (length(groups) != 2) {
    cat(paste0("  ERROR: Expected 2 groups, found ", length(groups), "\n"))
    return(NULL)
  }
  
  # Assume the first group is control and the second is case
  # This may need to be adjusted based on actual data
  control_label <- groups[1]
  case_label <- groups[2]
  
  cat(paste0("  Case label: ", case_label, "\n"))
  cat(paste0("  Control label: ", control_label, "\n"))
  
  # Create design matrix
  design <- model.matrix(~ 0 + pheno_data$group)
  colnames(design) <- c("control", "case")
  
  # Fit linear model
  fit <- lmFit(expr_data, design)
  
  # Create contrast matrix
  contrast.matrix <- makeContrasts(case - control, levels = design)
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Get DE results
  top_table <- topTable(fit2, n = Inf, adjust.method = "fdr")
  
  # Add gene symbols (assuming rownames are Entrez IDs or symbols)
  if (all(grepl("^[0-9]", rownames(top_table)))) {
    # Map Entrez IDs to symbols
    gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(top_table), column = "SYMBOL", keytype = "ENTREZID")
    top_table$GeneSymbol <- gene_symbols
  } else {
    # Assume rownames are already symbols
    top_table$GeneSymbol <- rownames(top_table)
  }
  
  # Create rank vector: rank_stat = log2FC
  rank_stat <- top_table$logFC
  names(rank_stat) <- top_table$GeneSymbol
  
  # Remove NA values
  rank_stat <- rank_stat[!is.na(names(rank_stat))]
  rank_stat <- rank_stat[!is.na(rank_stat)]
  
  # Map gene symbols to Entrez IDs for KEGG compatibility
  cat("  Mapping gene symbols to Entrez IDs...\n")
  tryCatch({
    # Use org.Hs.eg.db to map symbols to Entrez IDs
    entrez_ids <- mapIds(org.Hs.eg.db, keys = names(rank_stat), column = "ENTREZID", keytype = "SYMBOL")
    
    # Create a new rank_stat with Entrez IDs as names
    rank_stat_entrez <- rank_stat
    names(rank_stat_entrez) <- entrez_ids
    
    # Remove NA values
    rank_stat_entrez <- rank_stat_entrez[!is.na(names(rank_stat_entrez))]
    
    cat(paste0("  Rank vector length after mapping: ", length(rank_stat_entrez), "\n"))
    
    # Use the Entrez ID version for GSEA
    rank_stat <- rank_stat_entrez
  }, error = function(e) {
    cat("  WARNING: Gene ID mapping failed. Using gene symbols directly.\n")
    cat(paste0("  Error: ", e$message, "\n"))
  })
  
  cat(paste0("  Final rank vector length: ", length(rank_stat), "\n"))
  
  # Load KEGG gene sets
  cat("  Loading KEGG gene sets...\n")
  kegg_files <- c(
    file.path(base_dir, "data", "raw", "kegg", "KEGG_hsa_raw_list.rds"),
    file.path(base_dir, "01_data_raw", "kegg", "processed", "KEGG_hsa_raw_list.rds")
  )
  
  kegg_file_used <- NULL
  for (f in kegg_files) {
    if (file.exists(f)) {
      kegg_file_used <- f
      break
    }
  }
  
  if (is.null(kegg_file_used)) {
    cat("  ERROR: No KEGG gene sets found\n")
    return(NULL)
  }
  
  cat(paste0("  Using KEGG file: ", kegg_file_used, "\n"))
  
  # Load KEGG data
  kegg_data <- readRDS(kegg_file_used)
  
  # Process KEGG data
  if (all(c("KEGGPATHID2EXTID", "KEGGPATHID2NAME") %in% names(kegg_data))) {
    cat("  Processing KEGG mapping tables...\n")
    
    # Extract mappings
    pathway_to_gene <- kegg_data$KEGGPATHID2EXTID
    pathway_to_name <- kegg_data$KEGGPATHID2NAME
    
    # Convert to data frames
    pathway_to_gene_df <- as.data.frame(pathway_to_gene)
    pathway_to_name_df <- as.data.frame(pathway_to_name)
    
    # Rename columns for clarity
    colnames(pathway_to_gene_df) <- c("pathway_id", "gene_id")
    colnames(pathway_to_name_df) <- c("pathway_id", "pathway_name")
    
    # Build gene sets
    pathway_genes <- split(pathway_to_gene_df$gene_id, pathway_to_gene_df$pathway_id)
    
    # Create a data frame with pathway info
    pathway_info <- data.frame(
      pathway_id = names(pathway_genes),
      stringsAsFactors = FALSE
    )
    
    # Merge with pathway names
    pathway_info <- merge(pathway_info, pathway_to_name_df, by = "pathway_id", all.x = TRUE)
    
    # Create pathway names in the format "pathway_id - pathway_name"
    pathway_info$pathway <- paste0(pathway_info$pathway_id, " - ", pathway_info$pathway_name)
    
    # Build the gene sets list
    kegg_gene_sets <- pathway_genes
    names(kegg_gene_sets) <- pathway_info$pathway[match(names(kegg_gene_sets), pathway_info$pathway_id)]
    
    # Remove any pathways with NA names
    kegg_gene_sets <- kegg_gene_sets[!is.na(names(kegg_gene_sets))]
    
    cat(paste0("  KEGG gene sets count: ", length(kegg_gene_sets), "\n"))
  } else {
    # Assume it's already a list of gene sets
    kegg_gene_sets <- kegg_data
    cat(paste0("  KEGG gene sets count: ", length(kegg_gene_sets), "\n"))
  }
  
  # Check gene ID matching
  cat("  Checking gene ID matching...\n")
  # Extract all genes from KEGG sets
  kegg_genes <- unique(unlist(kegg_gene_sets))
  cat(paste0("  Total genes in KEGG sets: ", length(kegg_genes), "\n"))
  cat(paste0("  Total genes in rank_stat: ", length(rank_stat), "\n"))
  
  # Check overlap
  common_genes <- intersect(names(rank_stat), kegg_genes)
  cat(paste0("  Common genes: ", length(common_genes), "\n"))
  
  if (length(common_genes) == 0) {
    cat("  WARNING: No common genes between rank_stat and KEGG sets!\n")
    # Show sample genes from both
    cat("  Sample genes from rank_stat: ", paste(head(names(rank_stat), 5), collapse = ", "), "\n")
    cat("  Sample genes from KEGG sets: ", paste(head(kegg_genes, 5), collapse = ", "), "\n")
  }
  
  # Set seed for reproducibility
  set.seed(108)
  
  # Run fgsea
  cat("  Running fgsea...\n")
  fgsea_results <- fgsea(
    pathways = kegg_gene_sets,
    stats = rank_stat,
    minSize = 10,
    maxSize = 500,
    nperm = 10000
  )
  
  # Process results
  # Convert data.table to data frame first
  fgsea_results <- as.data.frame(fgsea_results)
  
  # Then convert to tibble
  fgsea_results <- as_tibble(fgsea_results)
  
  # Process results
  fgsea_results <- fgsea_results %>%
    dplyr::mutate(
      pathway_id = sub("^hsa", "hsa", pathway),
      pathway_name = gsub(" - Homo sapiens \\(human\\)", "", pathway),
      pvalue = pval  # Rename pval to pvalue for consistency
    ) %>%
    dplyr::select(
      pathway_id,
      pathway_name,
      NES,
      pvalue,
      padj,
      size,
      leadingEdge
    ) %>%
    dplyr::arrange(padj, pvalue)
  
  cat(paste0("  GSEA results rows: ", nrow(fgsea_results), "\n"))
  
  # Save GSEA results
  output_tables_dir <- file.path(base_dir, "output", "tables_main")
  if (!dir.exists(output_tables_dir)) {
    dir.create(output_tables_dir, recursive = TRUE)
  }
  
  gsea_output_file <- file.path(output_tables_dir, paste0("Fig1C_GSEA_", cohort_id, ".csv"))
  write_csv(fgsea_results, gsea_output_file)
  cat(paste0("  GSEA results saved to: ", gsea_output_file, "\n"))
  
  # Log the GSEA results file
  if (exists("log_output")) {
    log_output(
      file_path = gsea_output_file,
      file_type = "table",
      figure_id = paste0("Fig1C_GSEA_", cohort_id),
      description = paste0("GSEA results for Fig1C (KEGG pathways) using ", cohort_id, " DE log2FC"),
      script_name = "07a_fig1_discovery_common_pathways.R"
    )
  }
  
  return(fgsea_results)
}

# Function to build common pathways across replication cohorts
build_common_pathways_replication7 <- function(gsea_list) {
  cat("\n=== Building common pathways across 7 replication cohorts ===\n")
  
  # Combine all GSEA results
  all_results <- bind_rows(gsea_list, .id = "cohort_id")
  
  # Create audit table
  audit_table <- all_results %>%
    mutate(
      pass_fdr = padj <= 0.10
    )
  
  # Filter pathways with padj <= 0.10 for each cohort
  filtered_results <- audit_table %>%
    filter(pass_fdr)
  
  # Calculate cohort count and other metrics
  common_pathways <- filtered_results %>%
    group_by(pathway_id, pathway_name) %>%
    summarise(
      Cohort_Count = n(),
      Cohorts = paste(cohort_id, collapse = ","),
      Avg_NES = mean(NES),
      Avg_pvalue = mean(pvalue),
      Avg_FDR = mean(padj),
      .groups = "drop"
    ) %>%
    # Filter pathways present in at least 3 cohorts
    filter(Cohort_Count >= 3) %>%
    # Sort by Cohort_Count descending, then by absolute Avg_NES descending
    arrange(desc(Cohort_Count), desc(abs(Avg_NES)))
  
  # Add used_in_common flag to audit table
  audit_table <- audit_table %>%
    mutate(
      used_in_common = pathway_id %in% common_pathways$pathway_id
    )
  
  cat(paste0("  Common pathways found: ", nrow(common_pathways), "\n"))
  
  # Save common pathways
  output_tables_dir <- file.path(base_dir, "output", "tables_main")
  common_pathways_file <- file.path(output_tables_dir, "Fig1C_common_pathways_replication7.csv")
  write_csv(common_pathways, common_pathways_file)
  cat(paste0("  Common pathways saved to: ", common_pathways_file, "\n"))
  
  # Save audit table
  audit_file <- file.path(output_tables_dir, "Fig1C_common_pathways_replication7_audit.csv")
  write_csv(audit_table, audit_file)
  cat(paste0("  Audit table saved to: ", audit_file, "\n"))
  
  # Log the common pathways file
  if (exists("log_output")) {
    log_output(
      file_path = common_pathways_file,
      file_type = "table",
      figure_id = "Fig1C_common_pathways_replication7",
      description = "Common pathways across 7 replication cohorts for Fig1C",
      script_name = "07a_fig1_discovery_common_pathways.R"
    )
    
    # Log the audit table
    log_output(
      file_path = audit_file,
      file_type = "table",
      figure_id = "Fig1C_common_pathways_replication7_audit",
      description = "Audit table for common pathways across 7 replication cohorts",
      script_name = "07a_fig1_discovery_common_pathways.R"
    )
  }
  
  return(common_pathways)
}

# Run GSEA for all replication cohorts
cat("\n=== Running GSEA for all replication cohorts ===\n")
gsea_results_list <- list()

for (cohort_id in replication_only_cohorts) {
  gsea_result <- run_kegg_fgsea_for_cohort(cohort_id, base_dir)
  if (!is.null(gsea_result)) {
    gsea_results_list[[cohort_id]] <- gsea_result
  }
}

# Build common pathways
if (length(gsea_results_list) > 0) {
  common_pathways <- build_common_pathways_replication7(gsea_results_list)
} else {
  cat("\nERROR: No GSEA results found for any cohort\n")
}

# ==============================================================================
# 4. Fig1C_GSEA: Calculate DE log2FC and run KEGG GSEA for discovery cohort
# ==============================================================================

cat("\n=== Generating Fig1C_GSEA_results.csv for discovery cohort ===\n")

# Load required packages for DE and GSEA
suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(org.Hs.eg.db)
})

# Create output directory for tables
output_tables_dir <- file.path(output_dir, "tables_main")
if (!dir.exists(output_tables_dir)) {
  dir.create(output_tables_dir, recursive = TRUE)
  cat(paste0("Created directory: ", output_tables_dir, "\n"))
}

# ------------------------------------------------------------------------------
# 1) Input resolution: Expression matrix
# ------------------------------------------------------------------------------
cat("\n[Step 1] Resolving input files ...\n")

expr_file <- file.path(base_dir, "data", "processed", "expr_gse152004&GPL11154__FULL.rds")
cat(paste0("Checking expression file: ", expr_file, "\n"))
cat(paste0("File exists: ", file.exists(expr_file), "\n"))

if (!file.exists(expr_file)) {
  stop(paste0("Expression file not found: ", expr_file))
}

# Load expression data
expr_data <- readRDS(expr_file)
cat(paste0("Expression data dimensions: ", nrow(expr_data), " genes x ", ncol(expr_data), " samples\n"))

# ------------------------------------------------------------------------------
# 2) Input resolution: Phenotype data (T2 labels)
# ------------------------------------------------------------------------------
pheno_files <- c(
  file.path(base_dir, "data", "processed", "phase3_outputs", "GSE152004_T2_labels.csv"),
  file.path(base_dir, "data", "derived", "GSE152004_T2_labels.csv"),
  file.path(base_dir, "data", "processed", "GSE152004_T2_labels.csv")
)

pheno_file_used <- NULL
for (f in pheno_files) {
  if (file.exists(f)) {
    pheno_file_used <- f
    break
  }
}

if (is.null(pheno_file_used)) {
  stop("No T2 phenotype file found. Please check pheno files in data/processed or data/derived")
}

cat(paste0("Using phenotype file: ", pheno_file_used, "\n"))

# Load phenotype data
pheno_data <- read_csv(pheno_file_used, show_col_types = FALSE)

# Check for required columns
if (!all(c("Sample_ID", "T2_Status") %in% colnames(pheno_data))) {
  stop("Phenotype file missing required columns: Sample_ID or T2_Status")
}

# Check if column names are sample IDs (HRxxxx format)
colnames_are_sample_ids <- any(grepl("^HR", colnames(expr_data)))

if (colnames_are_sample_ids) {
  # Filter to relevant samples
  common_samples <- intersect(colnames(expr_data), pheno_data$Sample_ID)
  cat(paste0("Common samples: ", length(common_samples), "\n"))

  if (length(common_samples) == 0) {
    stop("No common samples between expression and phenotype data")
  }

  # Subset data
  expr_data <- expr_data[, common_samples]
  pheno_data <- pheno_data %>% filter(Sample_ID %in% common_samples)

  # Ensure order matches
  expr_data <- expr_data[, pheno_data$Sample_ID]
}

# ------------------------------------------------------------------------------
# 3) DE analysis: T2-high vs T2-low (exact replica of 02_discover_T2.R)
# ------------------------------------------------------------------------------
cat("\n[Step 2] DE analysis: T2-high vs T2-low ...\n")

# Create design matrix
design <- model.matrix(~ 0 + pheno_data$T2_Status)
colnames(design) <- c("T2low", "T2high")

# Fit linear model
fit <- lmFit(expr_data, design)

# Create contrast matrix
contrast.matrix <- makeContrasts(T2high - T2low, levels = design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get DE results
top_table <- topTable(fit2, n = Inf, adjust.method = "fdr")

# Add gene symbols (assuming rownames are Entrez IDs or symbols)
if (all(grepl("^[0-9]", rownames(top_table)))) {
  # Map Entrez IDs to symbols
  gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(top_table), column = "SYMBOL", keytype = "ENTREZID")
  top_table$GeneSymbol <- gene_symbols
} else {
  # Assume rownames are already symbols
  top_table$GeneSymbol <- rownames(top_table)
}

# Create rank vector: rank_stat = log2FC
rank_stat <- top_table$logFC
names(rank_stat) <- top_table$GeneSymbol

# Remove NA values
rank_stat <- rank_stat[!is.na(names(rank_stat))]
rank_stat <- rank_stat[!is.na(rank_stat)]

# Map gene symbols to Entrez IDs for KEGG compatibility
cat("  Mapping gene symbols to Entrez IDs...\n")
tryCatch({
  # Use org.Hs.eg.db to map symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = names(rank_stat), column = "ENTREZID", keytype = "SYMBOL")
  
  # Create a new rank_stat with Entrez IDs as names
  rank_stat_entrez <- rank_stat
  names(rank_stat_entrez) <- entrez_ids
  
  # Remove NA values
  rank_stat_entrez <- rank_stat_entrez[!is.na(names(rank_stat_entrez))]
  
  cat(paste0("  Rank vector length after mapping: ", length(rank_stat_entrez), "\n"))
  
  # Use the Entrez ID version for GSEA
  rank_stat <- rank_stat_entrez
}, error = function(e) {
  cat("  WARNING: Gene ID mapping failed. Using gene symbols directly.\n")
  cat(paste0("  Error: ", e$message, "\n"))
})

cat(paste0("  Final rank vector length: ", length(rank_stat), "\n"))

# Load KEGG gene sets
cat("  Loading KEGG gene sets...\n")
kegg_files <- c(
  file.path(base_dir, "data", "raw", "kegg", "KEGG_hsa_raw_list.rds"),
  file.path(base_dir, "01_data_raw", "kegg", "processed", "KEGG_hsa_raw_list.rds")
)

kegg_file_used <- NULL
for (f in kegg_files) {
  if (file.exists(f)) {
    kegg_file_used <- f
    break
  }
}

if (is.null(kegg_file_used)) {
  cat("  ERROR: No KEGG gene sets found\n")
  # Skip GSEA for discovery cohort
} else {
  cat(paste0("  Using KEGG file: ", kegg_file_used, "\n"))
  
  # Load KEGG data
  kegg_data <- readRDS(kegg_file_used)
  
  # Process KEGG data
  if (all(c("KEGGPATHID2EXTID", "KEGGPATHID2NAME") %in% names(kegg_data))) {
    cat("  Processing KEGG mapping tables...\n")
    
    # Extract mappings
    pathway_to_gene <- kegg_data$KEGGPATHID2EXTID
    pathway_to_name <- kegg_data$KEGGPATHID2NAME
    
    # Convert to data frames
    pathway_to_gene_df <- as.data.frame(pathway_to_gene)
    pathway_to_name_df <- as.data.frame(pathway_to_name)
    
    # Rename columns for clarity
    colnames(pathway_to_gene_df) <- c("pathway_id", "gene_id")
    colnames(pathway_to_name_df) <- c("pathway_id", "pathway_name")
    
    # Build gene sets
    pathway_genes <- split(pathway_to_gene_df$gene_id, pathway_to_gene_df$pathway_id)
    
    # Create a data frame with pathway info
    pathway_info <- data.frame(
      pathway_id = names(pathway_genes),
      stringsAsFactors = FALSE
    )
    
    # Merge with pathway names
    pathway_info <- merge(pathway_info, pathway_to_name_df, by = "pathway_id", all.x = TRUE)
    
    # Create pathway names in the format "pathway_id - pathway_name"
    pathway_info$pathway <- paste0(pathway_info$pathway_id, " - ", pathway_info$pathway_name)
    
    # Build the gene sets list
    kegg_gene_sets <- pathway_genes
    names(kegg_gene_sets) <- pathway_info$pathway[match(names(kegg_gene_sets), pathway_info$pathway_id)]
    
    # Remove any pathways with NA names
    kegg_gene_sets <- kegg_gene_sets[!is.na(names(kegg_gene_sets))]
    
    cat(paste0("  KEGG gene sets count: ", length(kegg_gene_sets), "\n"))
  } else {
    # Assume it's already a list of gene sets
    kegg_gene_sets <- kegg_data
    cat(paste0("  KEGG gene sets count: ", length(kegg_gene_sets), "\n"))
  }
  
  # Check gene ID matching
  cat("  Checking gene ID matching...\n")
  # Extract all genes from KEGG sets
  kegg_genes <- unique(unlist(kegg_gene_sets))
  cat(paste0("  Total genes in KEGG sets: ", length(kegg_genes), "\n"))
  cat(paste0("  Total genes in rank_stat: ", length(rank_stat), "\n"))
  
  # Check overlap
  common_genes <- intersect(names(rank_stat), kegg_genes)
  cat(paste0("  Common genes: ", length(common_genes), "\n"))
  
  if (length(common_genes) == 0) {
    cat("  WARNING: No common genes between rank_stat and KEGG sets!\n")
    # Show sample genes from both
    cat("  Sample genes from rank_stat: ", paste(head(names(rank_stat), 5), collapse = ", "), "\n")
    cat("  Sample genes from KEGG sets: ", paste(head(kegg_genes, 5), collapse = ", "), "\n")
  } else {
    # Set seed for reproducibility
    set.seed(108)
    
    # Run fgsea
    cat("  Running fgsea...\n")
    fgsea_results <- fgsea(
      pathways = kegg_gene_sets,
      stats = rank_stat,
      minSize = 10,
      maxSize = 500,
      nperm = 10000
    )
    
    # Process results
    # Convert data.table to data frame first
    fgsea_results <- as.data.frame(fgsea_results)
    
    # Then convert to tibble
    fgsea_results <- as_tibble(fgsea_results)
    
    # Process results
    # Convert leadingEdge from list to comma-separated string
    fgsea_results <- fgsea_results %>%
      dplyr::mutate(
        pathway_id = sub("^hsa", "hsa", pathway),
        pathway_name = gsub(" - Homo sapiens \\(human\\)", "", pathway),
        pvalue = pval,  # Rename pval to pvalue for consistency
        leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ", "), character(1))
      ) %>%
      dplyr::select(
        pathway_id,
        pathway_name,
        NES,
        pvalue,
        padj,
        size,
        leadingEdge
      ) %>%
      dplyr::arrange(padj, pvalue)
    
    cat(paste0("  GSEA results rows: ", nrow(fgsea_results), "\n"))
    
    # Save GSEA results for discovery cohort
    gsea_output_file <- file.path(output_tables_dir, "Fig1C_GSEA_results.csv")
    write_csv(fgsea_results, gsea_output_file)
    cat(paste0("  Discovery cohort GSEA results saved to: ", gsea_output_file, "\n"))
    
    # Log the GSEA results file
    if (exists("log_output")) {
      log_output(
        file_path = gsea_output_file,
        file_type = "table",
        figure_id = "Fig1C_GSEA_results",
        description = "GSEA results for Fig1C (KEGG pathways) using discovery cohort DE log2FC",
        script_name = "07a_fig1_discovery_common_pathways.R"
      )
    }
  }
}

cat("\n=== Common pathways generation completed ===\n")

# Close sink
if (sink.number() > 0) {
  sink()
}
