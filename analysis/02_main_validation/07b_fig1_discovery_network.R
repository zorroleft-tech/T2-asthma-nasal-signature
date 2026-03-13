#!/usr/bin/env Rscript

# ==============================================================================
# Script: 07b_fig1_discovery_network.R
# Purpose: Generate Figure 1B (density plot) and Figure 1C (KEGG network)
# Author: Zuo
# Date: 2026-03-07
#
# Inputs:
#   - data/derived/GSE152004_T2_labels.csv  # T2 signature scores and labels
#   - output/tables_main/Fig1C_common_pathways_replication7.csv  # Common pathways from 07a
#   - data/raw/kegg/KEGG_hsa_raw_list.rds  # KEGG pathway-to-gene mapping
#
# Outputs:
#   - output/figures_main/Fig1B_density.pdf  # T2 signature score density plot
#   - output/figures_main/Fig1C_KEGG_network.pdf  # KEGG pathway network
#   - output/figures_main/Fig1C_input_pathways.csv  # QC: Input pathways for Fig1C
#   - output/figures_main/Fig1C_module_annotation.csv  # QC: Module annotations for Fig1C
#   - output/figures_main/Fig1C_module_summary.csv  # QC: Module summary for Fig1C
#   - output/figures_main/Fig1C_hub_pathways.csv  # QC: Hub pathways for Fig1C
#   - output/figures_main/Fig1C_network_stats.csv  # QC: Network statistics for Fig1C
#   - output/figures_main/Fig1C_removed_singleton_modules.csv  # QC: Removed singleton modules
#   - output/logs/07b_fig1_discovery_network.log  # Log file
# ==============================================================================

# Load necessary packages
library(tidyverse)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggforce)
library(ggrepel)
library(rlang)

# ==============================================================================
# 1. Initialize base directory and setup logging
# ==============================================================================

# Set base directory - use script location to determine project root
tryCatch({
  # Try to get script path from RStudio
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
}, error = function(e) {
  # If not in RStudio, try to get from sys.frame
  tryCatch({
    script_dir <- dirname(sys.frame(1)$ofile)
  }, error = function(e) {
    # Fallback to current working directory
    script_dir <- getwd()
  })
})
base_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash="/", mustWork=TRUE)

# Verify project structure
required_dirs <- c("data", "analysis", "data_preparation")
for (dir in required_dirs) {
  if (!dir.exists(file.path(base_dir, dir))) {
    stop(paste0("Required directory not found: ", file.path(base_dir, dir)))
  }
}

# Create output directories
output_dir <- file.path(base_dir, "output")
figures_dir <- file.path(output_dir, "figures_main")
logs_dir <- file.path(output_dir, "logs")

for (dir in c(figures_dir, logs_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("Created directory: ", dir, "\n"))
  }
}

# Setup logging
log_file <- file.path(logs_dir, "07b_fig1_discovery_network.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({
  if (sink.number() > 0) {
    sink()
  }
})

cat("=== Starting Figure 1C network generation ===\n")
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
# 3. Figure 1B: Reference T2 signature score density plot
# ==============================================================================

cat("\n=== Generating Figure 1B ===\n")

# Read T2 labels data
labels_file <- file.path(base_dir, "data", "derived", "GSE152004_T2_labels.csv")
cat(paste0("Reading T2 labels file: ", labels_file, "\n"))
cat(paste0("File exists: ", file.exists(labels_file), "\n"))

# Add file existence check
if (!file.exists(labels_file)) {
  stop(paste0("T2 labels file not found: ", labels_file))
}

t2_labels <- read_csv(labels_file, show_col_types = FALSE)

# Check for required columns
if (!all(c("T2_Signature_Score", "T2_Status") %in% colnames(t2_labels))) {
  stop("T2 labels file missing required columns: T2_Signature_Score or T2_Status")
}

# Calculate median cutoff
median_cutoff <- median(t2_labels$T2_Signature_Score, na.rm = TRUE)
cat(paste0("Median cutoff: ", median_cutoff, "\n"))

# Count groups
t2_high_count <- sum(t2_labels$T2_Status == "T2-high", na.rm = TRUE)
t2_low_count <- sum(t2_labels$T2_Status == "T2-low", na.rm = TRUE)
cat(paste0("T2-high count: ", t2_high_count, "\n"))
cat(paste0("T2-low count: ", t2_low_count, "\n"))

# Create density plot
fig1b_plot <- ggplot(t2_labels, aes(x = T2_Signature_Score, fill = T2_Status)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = median_cutoff, linetype = "dashed") +
  scale_fill_manual(values = c("T2-high" = "#E64B35", "T2-low" = "#4DBBD5")) +
  labs(
    title = "Reference T2 Signature Score Distribution (Woodruff 2009)",
    x = "Reference T2 Signature Score",
    y = "Density",
    fill = "T2 Status"
  ) +
  theme_bw() +
  annotate("text", x = max(t2_labels$T2_Signature_Score, na.rm = TRUE) * 0.8, 
           y = max(density(t2_labels$T2_Signature_Score, na.rm = TRUE)$y) * 0.9,
           label = paste0("T2-high (n=", t2_high_count, ")\nT2-low (n=", t2_low_count, ")"),
           hjust = 1, vjust = 1)

# Save Figure 1B
fig1b_path <- file.path(figures_dir, "Fig1B_density.pdf")
ggsave(
  filename = fig1b_path,
  plot = fig1b_plot,
  width = 8,
  height = 6,
  dpi = 300,
  device = "pdf"
)

cat(paste0("Figure 1B saved to: ", fig1b_path, "\n"))
if (file.exists(fig1b_path)) {
  cat(paste0("File size: ", round(file.size(fig1b_path) / 1024, 2), " KB\n"))
} else {
  warning("Figure 1B not saved successfully")
}

# ==============================================================================
# 4. Figure 1C: KEGG GSEA pathway similarity network
# ==============================================================================

cat("\n=== Generating Figure 1C ===\n")

# Input files (use relative paths based on base_dir)
in_common_pathways <- file.path(base_dir, "output", "tables_main", "Fig1C_common_pathways_replication7.csv")
kegg_rds <- file.path(base_dir, "data", "raw", "kegg", "KEGG_hsa_raw_list.rds")

cat(paste0("Reading common pathways file: ", in_common_pathways, "\n"))
cat(paste0("File exists: ", file.exists(in_common_pathways), "\n"))
cat(paste0("Reading KEGG RDS file: ", kegg_rds, "\n"))
cat(paste0("File exists: ", file.exists(kegg_rds), "\n"))

# Check if input files exist
if (!file.exists(in_common_pathways) || !file.exists(kegg_rds)) {
  cat("⚠️ Fig1C inputs not found. Skipping Figure 1C (KEGG network).\n")
  cat("   Missing file(s):\n")
  if (!file.exists(in_common_pathways)) cat("   - ", in_common_pathways, "\n", sep="")
  if (!file.exists(kegg_rds)) cat("   - KEGG RDS file: ", kegg_rds, "\n", sep="")
  fig1c_path <- NULL
} else {
  # Load required packages
  suppressPackageStartupMessages({
    library(ggraph)
    library(ggforce)
    library(ggrepel)
  })

  # Read input data
  common_pathways <- read_csv(in_common_pathways, show_col_types = FALSE)
  
  # Build ptg2 from KEGG RDS file
  build_ptg2_from_kegg_rds <- function(kegg_rds_path) {
    kegg_obj <- readRDS(kegg_rds_path)

    # Case 1: clusterProfiler-style env/list with KEGGPATHID2EXTID
    if (is.list(kegg_obj) && all(c("KEGGPATHID2EXTID") %in% names(kegg_obj))) {
      ptg <- as.data.frame(kegg_obj$KEGGPATHID2EXTID)
      colnames(ptg) <- c("pathway_id", "gene_id")
      ptg$pathway_id <- sub("^path:", "", as.character(ptg$pathway_id))
      ptg$gene_id    <- as.character(ptg$gene_id)
      return(dplyr::as_tibble(ptg))
    }

    # Case 2: named list: each element is a gene vector; names are pathway ids or "hsaXXXX - name"
    if (is.list(kegg_obj)) {
      nm <- names(kegg_obj)

      pathway_id <- sub(" -.*$", "", nm)        # keep "hsa04060"
      pathway_id <- sub("^path:", "", pathway_id)

      ptg2 <- tibble::tibble(
        pathway_id = rep(pathway_id, lengths(kegg_obj)),
        gene_id    = as.character(unlist(kegg_obj, use.names = FALSE))
      ) |> 
        dplyr::filter(!is.na(pathway_id), !is.na(gene_id), nzchar(pathway_id), nzchar(gene_id)) |> 
        dplyr::distinct()

      return(ptg2)
    }

    stop("Unrecognized KEGG_hsa_raw_list.rds format.")
  }
  
  ptg2 <- build_ptg2_from_kegg_rds(kegg_rds)
  
  # Sanity check
  cat("\n=== Sanity check for ptg2 ===\n")
  cat("First few rows of ptg2:\n")
  print(head(ptg2))
  cat("Number of rows in ptg2:", nrow(ptg2), "\n")
  cat("Number of unique pathways:", length(unique(ptg2$pathway_id)), "\n")
  cat("Number of unique genes:", length(unique(ptg2$gene_id)), "\n")

  # Figure 1C parameters (locked)
  set.seed(108)
  min_cohort_count <- 3
  max_fdr <- 0.10
  max_nodes <- 35
  k_modules <- 6
  edge_jaccard_cut <- 0.15
  force_keep_one_edge <- TRUE
  drop_singleton_modules <- FALSE  # Set to FALSE to keep all pathways
  n_hubs_to_label <- 10

  cat("\nFigure 1C parameters:\n")
  cat(paste0("set.seed: 108\n"))
  cat(paste0("min_cohort_count: ", min_cohort_count, "\n"))
  cat(paste0("max_fdr: ", max_fdr, "\n"))
  cat(paste0("max_nodes: ", max_nodes, "\n"))
  cat(paste0("k_modules: ", k_modules, "\n"))
  cat(paste0("edge_jaccard_cut: ", edge_jaccard_cut, "\n"))
  cat(paste0("force_keep_one_edge: ", force_keep_one_edge, "\n"))
  cat(paste0("drop_singleton_modules: ", drop_singleton_modules, "\n"))
  cat(paste0("n_hubs_to_label: ", n_hubs_to_label, "\n"))

  # -----------------------------------------------------------------------------  
  # 1) Standardize fields (compatible with different column names)
  # -----------------------------------------------------------------------------  
  # ID
  if (!"ID" %in% names(common_pathways)) {
    if ("pathway_id" %in% names(common_pathways)) {
      # Extract just the pathway ID (e.g., "hsa04060" from "hsa04060 - ...")
      common_pathways$ID <- sub(" -.*$", "", common_pathways$pathway_id)
    } else {
      stop("common_pathways.csv missing column: ID or pathway_id")
    }
  }

  # Description
  if (!"Description" %in% names(common_pathways)) {
    if ("pathway_name" %in% names(common_pathways)) {
      common_pathways$Description <- common_pathways$pathway_name
    } else {
      stop("common_pathways.csv missing column: Description or pathway_name")
    }
  }

  # Cohort_Count
  if (!"Cohort_Count" %in% names(common_pathways)) {
    if ("n_cohorts" %in% names(common_pathways)) {
      common_pathways$Cohort_Count <- common_pathways$n_cohorts
    } else {
      stop("common_pathways.csv missing column: Cohort_Count or n_cohorts")
    }
  }

  # NES
  if (!"NES" %in% names(common_pathways)) {
    if ("Avg_NES" %in% names(common_pathways)) {
      common_pathways$NES <- common_pathways$Avg_NES
    } else {
      stop("common_pathways.csv missing column: NES or Avg_NES")
    }
  }

  # FDR
  if (!"FDR" %in% names(common_pathways)) {
    if ("Avg_FDR" %in% names(common_pathways)) {
      common_pathways$FDR <- common_pathways$Avg_FDR
    } else if ("p.adjust" %in% names(common_pathways)) {
      common_pathways$FDR <- common_pathways$p.adjust
    } else if ("Avg_pvalue" %in% names(common_pathways)) {
      common_pathways$FDR <- common_pathways$Avg_pvalue
    } else {
      stop("common_pathways.csv missing column: FDR / Avg_FDR / p.adjust / Avg_pvalue")
    }
  }

  # -----------------------------------------------------------------------------  
  # 2) Pathway filtering (Surgical Fix: Remove Digestion/Viral)
  # -----------------------------------------------------------------------------  
  cat("\n[Step 3] Filter pathways (Surgical Fix: Remove Digestion/Viral) ...\n")

  # Calculate gene set sizes (from pathway_to_gene.csv)
  pathway_gene_counts <- ptg2 %>%
    group_by(pathway_id) %>%
    summarise(setSize = n_distinct(gene_id), .groups = "drop")

  # 1. Absolute blacklist (The Kill List) - Add Digestion and Viral
  blacklist_kw <- c(
    # --- Added: Targeted removal ---
    "digestion", "absorption", "bile", "vitamin", "fatty acid",
    "viral", "virus", "infection", "influenza", "herpes", "measles", "hepatitis",
    # -----------------------
    # General impurities
    "cancer", "carcinoma", "disease",
    "metabolism", "biosynthesis", "degradation", "longevity", "ampk", "hippo",
    "mtor", "pi3k", "mapk", "ras", "rap1", "calcium", "camp", "cgmp", "foxo",
    "autophagy", "mitophagy", "apoptosis", "necroptosis", "ferroptosis",
    "cysteine", "methionine", "amino acid",
    "adherens", "gap junction", "tight junction", "focal adhesion",
    "platelet", "fluid shear stress", "endocrine", "thermogenesis",
    "cardiac", "muscle", "synaptic", "gaba", "cholinergic"
  )

  # 2. Absolute whitelist (The Core List) - Remove "ige" to prevent false match with "digestion"
  whitelist_kw <- c(
    "asthma", "airway", "lung", "bronchial",             # Disease/location
    "interleukin", "cytokine", "chemokine",              # Signaling molecules
    "immune", "inflammatory", "response",                # Processes
    "t cell", "b cell", "eosinophil", "mast cell",       # Cells
    "jak", "stat", "nf-kappa", "toll-like", "nod-like",  # Core pathways
    "mucin", "fc epsilon", "graft", "rejection",         # Specific features (Fc Epsilon covers IgE receptor)
    "iga", "intestinal immune", "antigen"                # Mucosal immunity (IgA is safe)
    # Note: Removed "ige" because it matches "digestion"
  )

  # 3. Prepare data
  processed_gsea <- common_pathways %>%
    mutate(
      ID = as.character(ID),
      Description = as.character(Description),
      Cohort_Count = as.numeric(Cohort_Count),
      NES = as.numeric(NES),
      FDR = as.numeric(FDR),
      Description_lower = tolower(Description)
    ) %>%
    filter(!is.na(ID), !is.na(Description)) %>%
    left_join(pathway_gene_counts, by = c("ID" = "pathway_id"))

  # 4. Perform filtering: statistical threshold first, then light blacklist
  filtered <- processed_gsea %>%
    filter(Cohort_Count >= min_cohort_count, FDR <= max_fdr) %>%
    filter(!purrr::map_lgl(Description_lower, function(x) any(sapply(blacklist_kw, function(k) grepl(k, x, fixed=TRUE))))) %>%
    arrange(desc(Cohort_Count), desc(abs(NES))) %>%
    slice_head(n = 35)

  # Output log
  cat("  - Blacklist keywords:", length(blacklist_kw), "\n")
  cat("  - Whitelist keywords:", length(whitelist_kw), "\n")
  cat("  - Final Selected Pathways:", nrow(filtered), "\n")
  print(filtered$Description)

  # Output input node list (for traceability)
  fig1c_input_pathways <- file.path(figures_dir, "Fig1C_input_pathways.csv")
  write_csv(
    filtered %>% dplyr::select(ID, Description, Cohort_Count, NES, FDR, setSize),
    fig1c_input_pathways
  )
  cat(paste0("Input pathways saved to: ", fig1c_input_pathways, "\n"))

  # -----------------------------------------------------------------------------  
  # 3) Prepare gene sets (from pathway_to_gene.csv) and compute Jaccard matrix
  # -----------------------------------------------------------------------------  
  cat("\n[Step 4] Build gene sets and compute Jaccard matrix ...\n")

  path_ids <- filtered$ID

  match_rate <- mean(path_ids %in% ptg2$pathway_id)
  cat("  - Pathway ID match rate:", round(match_rate, 3), "\n")
  if (match_rate < 0.8) {
    warning("pathway_id match rate < 0.8, please check if common_pathways.csv IDs are in hsa**** format")
  }

  # Build gene sets
  gene_sets <- lapply(path_ids, function(pid) unique(ptg2$gene_id[ptg2$pathway_id == pid]))
  names(gene_sets) <- path_ids

  # Remove pathways with no genes (should not happen in theory)
  non_empty <- sapply(gene_sets, function(x) length(x) > 0)
  if (any(!non_empty)) {
    cat("  - Remove empty gene sets:", sum(!non_empty), "\n")
    keep_ids <- names(gene_sets)[non_empty]
    filtered <- filtered %>% filter(ID %in% keep_ids)
    gene_sets <- gene_sets[keep_ids]
    path_ids <- filtered$ID
  }

  # Compute Jaccard matrix
  n <- length(path_ids)
  jacc <- matrix(0, nrow = n, ncol = n)
  rownames(jacc) <- path_ids
  colnames(jacc) <- path_ids
  diag(jacc) <- 1

  for (i in 1:(n - 1)) {
    gi <- gene_sets[[i]]
    for (j in (i + 1):n) {
      gj <- gene_sets[[j]]
      inter <- length(intersect(gi, gj))
      uni   <- length(union(gi, gj))
      sim   <- ifelse(uni == 0, 0, inter / uni)
      jacc[i, j] <- sim
      jacc[j, i] <- sim
    }
  }

  nz <- jacc[jacc > 0 & jacc < 1]
  cat("  - Jaccard >0 count:", length(nz), "\n")
  if (length(nz) > 0) {
    cat("  - Jaccard range (non-diagonal, >0): [",
        round(min(nz), 3), ", ", round(max(nz), 3), "]\n", sep = "")
  }

  # -----------------------------------------------------------------------------  
  # 4) Module detection: hierarchical clustering on (1 - Jaccard) distance + cutree(k=5~7)
  # -----------------------------------------------------------------------------  
  cat("\n[Step 5] Define modules by hierarchical clustering ...\n")

  dist_mat <- as.dist(1 - jacc)
  hc <- hclust(dist_mat, method = "average")
  
  # Adjust k_modules based on the number of pathways
  n_pathways <- length(path_ids)
  k_use <- min(k_modules, n_pathways)
  if (k_use < 2) k_use <- 1
  
  cat(paste0("  - Adjusted k_modules: ", k_use, " (original: ", k_modules, ")\n"))
  modules <- cutree(hc, k = k_use)

  filtered$module <- modules[match(filtered$ID, names(modules))]

  module_counts <- sort(table(filtered$module), decreasing = TRUE)
  cat("  - Modules:", length(module_counts), "\n")
  print(module_counts)

  # ---- Optional: drop singleton modules for main figure readability ----
  if (drop_singleton_modules) {
    module_sizes <- table(filtered$module)
    singleton_modules <- names(module_sizes[module_sizes == 1])
    
    if (length(singleton_modules) > 0) {
      removed_singletons <- filtered %>% filter(module %in% singleton_modules)
      
      # Record to QC (strongly recommended for traceability)
      fig1c_removed_singletons <- file.path(figures_dir, "Fig1C_removed_singleton_modules.csv")
      write_csv(
        removed_singletons %>% dplyr::select(ID, Description, module, Cohort_Count, NES, FDR, setSize),
        fig1c_removed_singletons
      )
      cat(paste0("Removed singleton modules saved to: ", fig1c_removed_singletons, "\n"))
      
      # Remove from main figure data
      filtered <- filtered %>% filter(!module %in% singleton_modules)
      
      cat("  - Dropped singleton modules:", paste(singleton_modules, collapse = ","), "\n")
      cat("  - Dropped pathways (singleton):", nrow(removed_singletons), "\n")
    }
  }
  
  # Check if there are any pathways left after filtering
  if (nrow(filtered) == 0) {
    cat("⚠️ No pathways left after filtering. Skipping Figure 1C (KEGG network).\n")
    fig1c_path <- NULL
  } else {

  # -----------------------------------------------------------------------------  
  # 5) Module naming (Manual Mapping & Validation)
  # -----------------------------------------------------------------------------  
  cat("\n[Step 6] Module Naming (Manual Mapping & Validation) ...\n")

  # 1. First print pathways in each module for manual mapping reference
  # ----------------------------------------------------------
  cat("  - Diagnostic: Content of each module:\n")
  module_info <- list()
  for (m in sort(unique(filtered$module))) {
    descs <- filtered$Description[filtered$module == m]
    # Simplify printing
    short_descs <- gsub(" - Homo sapiens \\(human\\)", "", descs)
    cat(sprintf("    [Module %d] (n=%d): %s\n", m, length(descs), paste(head(short_descs, 3), collapse = ", ")))
    
    # Save module info for manual mapping
    module_info[[as.character(m)]] <- list(
      descs = descs,
      short_descs = short_descs,
      n = length(descs),
      top3 = head(short_descs, 3)
    )
  }

  # 2. Manual module name mapping (based on diagnostic results)
  # ----------------------------------------------------------
  # Manual mapping based on module content
  manual_module_mapping <- function(module_id, descs) {
    txt <- paste(tolower(descs), collapse = " ")
    
    # M2: Neuro-Immune - Highest priority to avoid misclassification
    # Only match truly neuro-related pathways, not TRP-containing pathways that are mainly immune-related
    if (grepl("neuroactive ligand", txt)) {
      return("Neuro-Immune")
    }
    
    # M1: Adaptive Immunity & FcεRI (Core immunity + Fc receptors)
    if (grepl("t cell|b cell|fc epsilon", txt)) {
      return("Adaptive Immunity & FcεRI")
    }
    
    # M6: Innate Immunity (TLR/NOD-NF-kappaB) - ASCII-friendly version
    if (grepl("toll-like|nod-like|nf-kappa", txt)) {
      return("Innate Immunity (TLR/NOD-NF-kappaB)")
    }
    
    # M5: Cytokine-JAK/STAT (Cytokine signaling) - ASCII-friendly version
    if (grepl("cytokine|jak-stat", txt)) {
      return("Cytokine-JAK/STAT")
    }
    
    # M3: Antigen Presentation & Asthma / Mucosal (Antigen presentation + Asthma + Mucosal)
    if (grepl("antigen|asthma|allograft|intestinal immune|iga", txt)) {
      return("Antigen Presentation & Asthma / Mucosal")
    }
    
    # Default
    return("Other Pathways")
  }

  # Apply manual mapping
  # First create an empty module name vector
  module_names <- character(nrow(filtered))

  # Process each unique module
  unique_modules <- sort(unique(filtered$module))
  for (m in unique_modules) {
    # Get pathway descriptions for this module
    descs <- filtered$Description[filtered$module == m]
    module_size <- length(descs)
    
    # Generate module name based on module size
    if (module_size == 1) {
      # Rule 2: Module size = 1, use the node label as module name
      # Extract label from Description (part after "- ")
      if (grepl("- ", descs[1])) {
        name <- sub("^.*- ", "", descs[1])
      } else {
        name <- filtered$ID[filtered$module == m][1]
      }
    } else {
      # Rule 2: Module size ≥ 2, use thematic name
      name <- manual_module_mapping(m, descs)
    }
    
    # Add module number
    full_name <- paste0(name, " (M", m, ")")
    # Assign to all rows in this module
    module_names[filtered$module == m] <- full_name
  }

  # Assign to filtered dataframe
  filtered$module_name <- module_names

  # 3. Output validation information (ensure names are correctly assigned to modules)
  # ----------------------------------------------------------
  cat("\n  - Final Module Names (Validation):\n")
  validation_data <- list()
  for (m in sort(unique(filtered$module))) {
    module_name <- unique(filtered$module_name[filtered$module == m])[1]
    module_size <- sum(filtered$module == m)
    top3_paths <- head(filtered$Description[filtered$module == m], 3)
    short_top3 <- gsub(" - Homo sapiens \\(human\\)", "", top3_paths)
    
    cat(sprintf("    M%d: %s (n=%d) | Top3: %s\n", 
              m, module_name, module_size, paste(short_top3, collapse = ", ")))
    
    # Save validation data
    validation_data[[as.character(m)]] <- list(
      module_id = m,
      module_name = module_name,
      n = module_size,
      top3_pathways = paste(short_top3, collapse = " | ")
    )
  }

  # 4. Output module annotation table
  # ----------------------------------------------------------
  fig1c_module_annotation <- file.path(figures_dir, "Fig1C_module_annotation.csv")
  write_csv(
    filtered %>% dplyr::select(ID, Description, module, module_name, Cohort_Count, NES, FDR),
    fig1c_module_annotation
  )
  cat(paste0("Module annotation saved to: ", fig1c_module_annotation, "\n"))

  # 5. Generate module summary table (for validation)
  # ----------------------------------------------------------
  module_summary <- filtered %>%
    group_by(module, module_name) %>%
    summarise(
      n_pathways = n(),
      avg_abs_NES = mean(abs(NES), na.rm = TRUE),
      top_pathways = paste(head(Description[order(abs(NES), decreasing = TRUE)], 3), collapse = " | "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_pathways), desc(avg_abs_NES))

  fig1c_module_summary <- file.path(figures_dir, "Fig1C_module_summary.csv")
  write_csv(module_summary, fig1c_module_summary)
  cat(paste0("Module summary saved to: ", fig1c_module_summary, "\n"))

  # -----------------------------------------------------------------------------  
  # 6) Plotting edges: Jaccard > edge_jaccard_cut
  #    Optional: Force each node to keep at least one strongest edge (avoid losing nodes)
  # -----------------------------------------------------------------------------  
  cat("\n[Step 7] Build plotting graph ...\n")

  # Update path_ids to pathways after removing singletons
  path_ids <- filtered$ID

  # Recalculate jacc matrix (only includes remaining nodes)
  n <- length(path_ids)
  jacc_filtered <- matrix(0, nrow = n, ncol = n)
  rownames(jacc_filtered) <- path_ids
  colnames(jacc_filtered) <- path_ids
  diag(jacc_filtered) <- 1

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      jacc_filtered[i, j] <- jacc[path_ids[i], path_ids[j]]
      jacc_filtered[j, i] <- jacc[path_ids[i], path_ids[j]]
    }
  }

  # Get threshold edges
  idx <- which(jacc_filtered > edge_jaccard_cut & jacc_filtered < 1, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

  edge_list <- tibble(
    from = rownames(jacc_filtered)[idx[, 1]],
    to   = colnames(jacc_filtered)[idx[, 2]],
    weight = jacc_filtered[idx]
  )

  # Force each node to keep at least one strongest edge (only when nodes would be lost)
  if (force_keep_one_edge) {
    for (pid in path_ids) {
      # Check if this pid currently has any edges
      has_edge <- any(edge_list$from == pid | edge_list$to == pid)
      if (!has_edge) {
        # Find strongest neighbor (exclude self)
        v <- jacc_filtered[pid, ]
        v[names(v) == pid] <- 0
        if (max(v, na.rm = TRUE) > 0) {
          nb <- names(which.max(v))
          w  <- max(v, na.rm = TRUE)
          # Standardize from/to order
          a <- sort(c(pid, nb))[1]
          b <- sort(c(pid, nb))[2]
          edge_list <- bind_rows(edge_list, tibble(from = a, to = b, weight = w))
        }
      }
    }
    edge_list <- edge_list %>% distinct(from, to, .keep_all = TRUE)
  }

  cat("  - Edges plotted (Jaccard >", edge_jaccard_cut, "):", nrow(edge_list), "\n")

  # vertices
  vertices <- filtered %>%
    transmute(
      name = ID,
      # Add a column for short plot labels (simplified version)
      # Rule 1: If Description contains "- ", use the part after it
      # Rule 2: Otherwise, use pathway_id as fallback
      plot_label = ifelse(grepl("- ", Description),
                         stringr::str_wrap(sub("^.*- ", "", Description), width = 15),
                         stringr::str_wrap(ID, width = 15)),
      Description = Description,
      module = as.factor(module),
      module_name = module_name,
      NES = NES,
      Cohort_Count = Cohort_Count,
      FDR = FDR
    )

  g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = vertices)
  V(g)$degree <- degree(g)

  cat("  - Nodes in graph:", vcount(g), " (from selected:", nrow(filtered), ")\n")
  cat("  - Avg degree:", round(mean(V(g)$degree), 2), "\n")
  cat("  - Density:", round(edge_density(g), 3), "\n")

  # -----------------------------------------------------------------------------  
  # 7) Identify hubs (only label these)
  # -----------------------------------------------------------------------------  
  hub_nodes <- tibble(
    ID = V(g)$name,
    degree = V(g)$degree
  ) %>%
    arrange(desc(degree), ID) %>%
    slice_head(n = min(n_hubs_to_label, nrow(.)))

  hub_ids <- hub_nodes$ID

  # Save hub pathways
  fig1c_hub_pathways <- file.path(figures_dir, "Fig1C_hub_pathways.csv")
  hub_table <- filtered %>%
    filter(ID %in% hub_ids) %>%
    left_join(hub_nodes, by = c("ID")) %>%
    dplyr::select(ID, Description, module, module_name, NES, FDR, degree) %>%
    arrange(desc(degree), desc(abs(NES)))
  write_csv(hub_table, fig1c_hub_pathways)
  cat(paste0("Hub pathways saved to: ", fig1c_hub_pathways, "\n"))

  # -----------------------------------------------------------------------------  
  # 8) Plot (Final Fix: Leader Lines & Legend)
  # -----------------------------------------------------------------------------  
  cat("\n[Step 8] Plot (Final Fix: Leader Lines & Legend) ...\n")

  # Module colors (colorblind-friendly Nature style)
  module_levels <- sort(unique(vertices$module))
  palette <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
  module_colors <- setNames(palette[seq_along(module_levels)], module_levels)

  layout <- create_layout(g, layout = "fr", niter = 10000) # Increase iterations for more舒展 layout

  p <- ggraph(layout) +
    # 1. Background clusters (Hull) - Add module labels and label lines
    geom_mark_hull(
      aes(x, y, group = module, fill = module, label = module_name), # Add module name labels
      concavity = 4,
      expand = unit(3, "mm"),
      alpha = 0.15,
      # Label style
      label.fontsize = 9,     # Font size
      label.fontface = "bold", # Bold font
      label.colour = "grey30",  # Font color
      label.buffer = unit(8, "mm"), # Distance from label to cluster
      con.cap = 0,             # Connector end style
      con.colour = "grey40",   # Connector color
      con.size = 0.5,          # Connector thickness
      show.legend = FALSE      # No legend for clusters
    ) +
    
    # 2. Edges
    geom_edge_link(aes(width = weight, alpha = weight),
                   color = "grey50",
                   show.legend = FALSE) +
    scale_edge_width(range = c(0.7, 1.5)) + # Slightly thicker
    scale_edge_alpha(range = c(0.4, 0.8)) +
    
    # 3. Nodes (solid points)
    geom_node_point(aes(size = abs(NES), fill = module),
                    shape = 21, color = "white", stroke = 1.2) +
    
    # 4. Color legend settings
    scale_fill_manual(
      values = module_colors,
      name = "Functional Module",
      guide = guide_legend(override.aes = list(size = 5), order = 2) # Larger legend points, placed later
    ) +
    
    # 5. Size legend settings
    scale_size_continuous(
      range = c(4, 12), # Slightly larger nodes for more impact
      name = "|NES| (absolute NES)",
      guide = guide_legend(override.aes = list(fill = "grey50", color = "white"), order = 1)
    ) +
    
    # 6. Node text labels (with leader lines!)
    geom_node_text(
      data = subset(layout, name %in% hub_ids),
      aes(label = plot_label),
      size = 3.5,                 # Slightly larger font
      fontface = "bold",
      repel = TRUE,
      bg.color = "white",
      bg.r = 0.15,
      
      # --- [Key Fix] Force show leader lines ---
      segment.color = "grey30",   # Darker line color for visibility
      segment.size = 0.6,         # Thicker lines
      min.segment.length = 0,     # Force drawing lines regardless of distance
      force = 2,                  # Increase repulsion to push labels away and make space for lines
      # -------------------------------
      
      box.padding = unit(0.6, "lines"),   
      point.padding = unit(0.6, "lines"), # Increase space around points
      max.overlaps = Inf
    ) +
    
    # 7. Theme settings (keep large right margin)
    theme_void() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
      # Keep 120 on the right to ensure legend is never cut
      plot.margin = margin(30, 120, 30, 30),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    # Remove title and subtitle for cleaner figure
    labs()

  # Save Figure 1C
  fig1c_path <- file.path(figures_dir, "Fig1C_KEGG_network.pdf")
  ggsave(fig1c_path, plot = p, width = 15, height = 11, device = "pdf")

  cat(paste0("Figure 1C saved to: ", fig1c_path, "\n"))
  if (file.exists(fig1c_path)) {
    cat(paste0("File size: ", round(file.size(fig1c_path) / 1024, 2), " KB\n"))
  } else {
    warning("Figure 1C not saved successfully")
  }

  # -----------------------------------------------------------------------------  
  # 9) Export network statistics (for Results section & reproducibility)
  # -----------------------------------------------------------------------------  
  cat("\n[Step 9] Export network statistics ...\n")

  network_stats <- tibble(
    Metric = c(
      "Input pathways (after filtering)",
      "Modules (k)",
      "Edges plotted (thresholded + keep-one-edge)",
      "Nodes plotted (in graph)",
      "Average degree",
      "Network density",
      "Edge jaccard cutoff"
    ),
    Value = c(
      nrow(filtered),
      k_use,
      nrow(edge_list),
      vcount(g),
      round(mean(V(g)$degree), 2),
      round(edge_density(g), 3),
      edge_jaccard_cut
    )
  )

  fig1c_network_stats <- file.path(figures_dir, "Fig1C_network_stats.csv")
  write_csv(network_stats, fig1c_network_stats)
  cat(paste0("Network stats saved to: ", fig1c_network_stats, "\n"))

  # Log network stats
  cat("\nNetwork statistics:\n")
  cat(paste0("Nodes: ", vcount(g), "\n"))
  cat(paste0("Edges: ", ecount(g), "\n"))
  cat(paste0("Density: ", round(edge_density(g), 4), "\n"))
  cat(paste0("Average degree: ", round(mean(V(g)$degree), 2), "\n"))
  cat(paste0("Module size distribution:\n"))
  print(sort(table(filtered$module), decreasing = TRUE))
  }
}

# ==============================================================================
# 5. Log outputs to index
# ==============================================================================

cat("\n=== Logging outputs ===\n")

# Log Figure 1B
log_output(
  file_path = fig1b_path,
  file_type = "figure",
  figure_id = "Fig 1B",
  description = "Reference T2 signature score density in discovery cohort (GSE152004) with locked median cutoff",
  script_name = "07b_fig1_discovery_network.R"
)

# Log Figure 1C if it exists and was generated
if (!is.null(fig1c_path) && exists("fig1c_path") && file.exists(fig1c_path)) {
  log_output(
    file_path = fig1c_path,
    file_type = "figure",
    figure_id = "Fig 1C",
    description = "KEGG GSEA pathway similarity network (seed=108; Cohort_Count>=3; FDR<=0.10; max 35 pathways; k=6 modules)",
    script_name = "07b_fig1_discovery_network.R"
  )
}

cat("\n=== Figure 1B and 1C generation completed ===\n")

# Close sink
if (sink.number() > 0) {
  sink()
}
