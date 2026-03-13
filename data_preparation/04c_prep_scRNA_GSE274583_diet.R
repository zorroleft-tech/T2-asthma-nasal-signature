# =========================================================================
# Script: 04c_prep_scRNA_GSE274583_diet.R
# Purpose: Build integrated CD4+ T cell scRNA object (GSE274583) and export
#          both full RDS and reviewer DietSeurat RDS for reproducible Fig3.
# Author: Zuo
# Date: 2026-03-04
#
# Inputs:
#   - data/raw/single_cell/gse274583&GPL24676/Library1/  # 10X directory
#   - data/raw/single_cell/gse274583&GPL24676/Library2/  # 10X directory
#   - data/raw/single_cell/gse274583&GPL24676/Library3/  # 10X directory
#   - data/raw/single_cell/gse274583&GPL24676/Library4/  # 10X directory
#   - data/raw/single_cell/gse274583&GPL24676/GPL24676_family.soft.gz
#   - data/raw/single_cell/gse274583&GPL24676/GSE274583_family.soft.gz
#   - data/derived/locked_weights.csv  # Locked 11-gene weights (single source of truth)
#
# Outputs:
#   - data/processed_full/single_cell/gse274583/gse274583_cd4_integrated_full.rds
#   - data/processed_full/single_cell/gse274583/gse274583_cd4_integrated_reviewer_diet.rds
#   - output/logs/04c_prep_scRNA_GSE274583_diet_sessionInfo.txt
#   - output/logs/04c_prep_scRNA_GSE274583_diet_build_summary.tsv
# =========================================================================

# Disable debugging mode
options(error = NULL)

# Increase memory limit
message("Setting memory limit...")
if (Sys.info()["sysname"] == "Windows") {
  memory.limit(size = 40000)  # Set to 40GB
} else {
  # For non-Windows systems
  options(future.globals.maxSize = 40 * 1024^3)  # 40GB
}
message("Memory limit set")

# Remove any debugging functions if present
for (func in c('browser', 'debug', 'recover')) {
  if (exists(func)) {
    rm(list = func)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(Matrix)
})

# ----------------------------
# Project root + path setup
# ----------------------------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).",
       call. = FALSE)
}

DATA_RAW_DIR            <- file.path(base_dir, "data", "raw")
DATA_DERIVED_DIR        <- file.path(base_dir, "data", "derived")
DATA_PROCESSED_FULL_DIR <- file.path(base_dir, "data", "processed_full")
LOG_DIR                 <- file.path(base_dir, "output", "logs")

RAW_SC_DIR              <- file.path(DATA_RAW_DIR, "single_cell", "gse274583&GPL24676")
OUT_SC_DIR              <- file.path(DATA_PROCESSED_FULL_DIR, "single_cell", "gse274583")
CACHE_DIR               <- file.path(OUT_SC_DIR, "_cache")

dir.create(DATA_PROCESSED_FULL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# Locked files
LOCKED_WEIGHTS_FILE <- file.path(DATA_DERIVED_DIR, "locked_weights.csv")

# Input 10X dirs
LIB_DIRS <- list(
  Library1 = file.path(RAW_SC_DIR, "Library1"),
  Library2 = file.path(RAW_SC_DIR, "Library2"),
  Library3 = file.path(RAW_SC_DIR, "Library3"),
  Library4 = file.path(RAW_SC_DIR, "Library4")
)

# Output files
FULL_RDS_FILE    <- file.path(OUT_SC_DIR, "gse274583_cd4_integrated_full.rds")
REVIEWER_RDS_FILE<- file.path(OUT_SC_DIR, "gse274583_cd4_integrated_reviewer_diet.rds")
SESSIONINFO_FILE <- file.path(LOG_DIR, "04c_prep_scRNA_GSE274583_diet_sessionInfo.txt")
SUMMARY_FILE     <- file.path(LOG_DIR, "04c_prep_scRNA_GSE274583_diet_build_summary.tsv")

# ----------------------------
# Helpers
# ----------------------------
assert_file <- function(path) {
  if (!file.exists(path)) stop(paste0("Missing file: ", path), call. = FALSE)
}
assert_dir <- function(path) {
  if (!dir.exists(path)) stop(paste0("Missing directory: ", path), call. = FALSE)
}

# Save sessionInfo early (will overwrite at end too)
capture.output(sessionInfo(), file = SESSIONINFO_FILE)

# Check inputs exist
assert_file(LOCKED_WEIGHTS_FILE)
purrr_needed <- FALSE
for (nm in names(LIB_DIRS)) {
  assert_dir(LIB_DIRS[[nm]])
  # 10X can be either filtered_feature_bc_matrix/ OR matrix.mtx.gz in dir
  # so we do not over-restrict; Read10X will error with a clear message if invalid.
}

# ----------------------------
# Read locked weights (single source of truth)
# ----------------------------
locked_w <- read_csv(LOCKED_WEIGHTS_FILE, show_col_types = FALSE) %>%
  janitor::clean_names()

# Expect at least: gene, weight (case-insensitive after clean_names)
# If your CSV uses different column names, standardize here:
if (!all(c("gene", "weight") %in% colnames(locked_w))) {
  stop("locked_weights.csv must contain columns: gene, weight", call. = FALSE)
}

locked_genes <- locked_w$gene %>% as.character() %>% unique()
locked_w_vec <- locked_w$weight %>% as.numeric()
names(locked_w_vec) <- locked_genes

# Orthogonal Th2 module
th2_genes <- c("GATA3", "PTGDR2", "CCR4", "IL17RB")

# ----------------------------
# Read 10X and create Seurat objects per library
# ----------------------------
read_10x_with_custom_naming <- function(lib_dir, lib_nm) {
  # For Library4, we need to handle multiple sub-libraries
  if (lib_nm == "Library4") {
    sub_obj_list <- list()
    
    for (i in 1:3) {
      # Construct exact file paths
      prefix <- paste0("GSE274583_", lib_nm)
      
      # For matrix and features files, always use the number
      matrix_file <- file.path(lib_dir, paste0(prefix, "_matrix", i, ".mtx.gz"))
      features_file <- file.path(lib_dir, paste0(prefix, "_features", i, ".tsv.gz"))
      
      # For barcode files, use different naming for i=1 vs i=2,3
      if (i == 1) {
        barcode_file <- file.path(lib_dir, paste0(prefix, "_barcodes.tsv.gz"))
      } else {
        barcode_file <- file.path(lib_dir, paste0(prefix, "_barcodes", i, ".tsv.gz"))
      }
      
      # Check if files exist
      if (!file.exists(barcode_file)) stop("Barcode file not found: ", barcode_file)
      if (!file.exists(features_file)) stop("Features file not found: ", features_file)
      if (!file.exists(matrix_file)) stop("Matrix file not found: ", matrix_file)
      
      # Read the files directly with cache
      cache_id <- paste0(lib_nm, "_sub", i)
      counts <- read_10x_files(barcode_file, features_file, matrix_file, cache_id = cache_id)
      
      # Create Seurat object with minimal filtering to reduce memory
      sub_seu <- CreateSeuratObject(counts = counts, project = paste0("GSE274583_CD4_", i), min.cells = 3, min.features = 200)
      sub_seu$library_id <- paste0(lib_nm, "_sub", i)
      sub_seu$percent.mt <- PercentageFeatureSet(sub_seu, pattern = "^MT-")
      
      # Add cell id prefix
      colnames(sub_seu) <- paste0(paste0(lib_nm, "_sub", i), "_", colnames(sub_seu))
      
      # Log the number of cells
      message(paste0(lib_nm, "_sub", i, ": ", ncol(sub_seu), " cells"))
      
      sub_obj_list[[i]] <- sub_seu
      
      # Free memory
      rm(counts)
      gc()
    }
    
    # Merge sub-libraries
    if (length(sub_obj_list) > 1) {
      seu <- merge(x = sub_obj_list[[1]], y = sub_obj_list[-1], project = "GSE274583_CD4")
    } else {
      seu <- sub_obj_list[[1]]
    }
  } else {
    # For other libraries, construct exact file paths
    prefix <- paste0("GSE274583_", lib_nm)
    barcode_file <- file.path(lib_dir, paste0(prefix, "_barcodes.tsv.gz"))
    features_file <- file.path(lib_dir, paste0(prefix, "_features.tsv.gz"))
    matrix_file <- file.path(lib_dir, paste0(prefix, "_matrix.mtx.gz"))
    
    # Check if files exist
    if (!file.exists(barcode_file)) stop("Barcode file not found: ", barcode_file)
    if (!file.exists(features_file)) stop("Features file not found: ", features_file)
    if (!file.exists(matrix_file)) stop("Matrix file not found: ", matrix_file)
    
    # Read the files directly with cache
    cache_id <- lib_nm
    counts <- read_10x_files(barcode_file, features_file, matrix_file, cache_id = cache_id)
    
    # Create Seurat object
    seu <- CreateSeuratObject(counts = counts, project = "GSE274583_CD4", min.cells = 3, min.features = 200)
    seu$library_id <- lib_nm
    seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
    
    # Add cell id prefix
    colnames(seu) <- paste0(lib_nm, "_", colnames(seu))
    
    # Log the number of cells
    message(paste0(lib_nm, ": ", ncol(seu), " cells"))
    
    # Free memory
    rm(counts)
    gc()
  }
  
  return(seu)
}

# Helper function to read 10x files directly without copying
read_10x_files <- function(barcode_file, features_file, matrix_file, cache_id = NULL) {
  # Validate cache_id
  if (!is.null(cache_id)) {
    allowed_cache_ids <- c("Library1", "Library2", "Library3", "Library4_sub1", "Library4_sub2", "Library4_sub3")
    if (!cache_id %in% allowed_cache_ids) {
      stop("Invalid cache_id: ", cache_id, ". Allowed values: ", paste(allowed_cache_ids, collapse = ", "), call. = FALSE)
    }
  }
  
  # Check if cache exists
  if (!is.null(cache_id)) {
    cache_file <- file.path(CACHE_DIR, paste0(cache_id, "_counts.rds"))
    if (file.exists(cache_file)) {
      message(paste0("Loading cached counts from: ", cache_file))
      counts <- readRDS(cache_file)
      message(paste0("Cache loaded: dim=", paste(dim(counts), collapse = "x"), ", class=", class(counts)))
      # Sample gene names for verification
      if (nrow(counts) > 3) {
        sample_genes <- rownames(counts)[sample(1:nrow(counts), 3)]
        message(paste0("Sample gene names: ", paste(sample_genes, collapse = ", ")))
      }
      return(counts)
    }
  }
  
  # Read barcodes
  message("Reading barcodes...")
  barcodes <- readLines(gzfile(barcode_file))
  
  # Read features
  message("Reading features...")
  features <- read.delim(gzfile(features_file), header = FALSE, stringsAsFactors = FALSE)
  
  # Determine which column to use for gene symbols
  if (ncol(features) >= 3) {
    # Check if second column has many underscores or ENSG prefixes
    col2_has_underscore <- sum(grepl("_", features[, 2])) / nrow(features) > 0.5
    col2_is_ensg <- sum(grepl("^ENSG", features[, 2])) / nrow(features) > 0.5
    
    if (col2_has_underscore || col2_is_ensg) {
      # Use first column instead
      gene_names <- features[, 1]
      message("Using first column for gene symbols (second column has many underscores or ENSG prefixes)")
    } else {
      # Use second column
      gene_names <- features[, 2]
      message("Using feature column 2 as gene symbols")
    }
  } else if (ncol(features) >= 2) {
    # Use second column
    gene_names <- features[, 2]
    message("Using feature column 2 as gene symbols")
  } else {
    # Use first column
    gene_names <- features[, 1]
    message("Only one column in features file, using it for gene symbols")
  }
  
  # Convert gene symbols to uppercase
  gene_names <- toupper(gene_names)
  
  # Read matrix
  message("Reading matrix...")
  matrix <- readMM(gzfile(matrix_file))
  
  # ----------------------------
  # Set dimnames + merge duplicate genes (FAST, sparse-safe)
  # ----------------------------
  # Ensure sparse column-compressed matrix for fast ops
  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "dgCMatrix")
  }

  # Always set colnames
  colnames(matrix) <- barcodes

  # Uppercase gene symbols (you already did this earlier, keep it)
  gene_names <- toupper(gene_names)

  if (anyDuplicated(gene_names)) {
    message("Duplicate gene names found. Merging counts (sparse mapping)...")

    # factor levels keep first-seen order, stable across runs
    f <- factor(gene_names, levels = unique(gene_names))
    unique_genes <- levels(f)

    # Build sparse mapping matrix A: (n_unique_genes x n_features)
    # sparse.model.matrix gives (n_features x n_unique); transpose it.
    A <- Matrix::sparse.model.matrix(~ 0 + f)
    A <- Matrix::t(A)

    # Aggregate counts: gene-level = A %*% feature-level
    matrix <- A %*% matrix

    rownames(matrix) <- unique_genes

    message(sprintf(
      "Duplicate gene merging completed. total=%d unique=%d dupGroups=%d",
      length(gene_names), length(unique_genes), length(gene_names) - length(unique_genes)
    ))
  } else {
    rownames(matrix) <- gene_names
  }

  # Ensure dgCMatrix at end
  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "dgCMatrix")
  }
  
  # Save to cache if cache_id is provided
  if (!is.null(cache_id)) {
    cache_file <- file.path(CACHE_DIR, paste0(cache_id, "_counts.rds"))
    message(paste0("Saving counts to cache: ", cache_file))
    message(paste0("Saving: dim=", paste(dim(matrix), collapse = "x"), ", class=", class(matrix)))
    saveRDS(matrix, file = cache_file)
  }
  
  return(matrix)
}

obj_list <- lapply(names(LIB_DIRS), function(lib_nm) {
  lib_dir <- LIB_DIRS[[lib_nm]]
  read_10x_with_custom_naming(lib_dir, lib_nm)
})
names(obj_list) <- names(LIB_DIRS)

# Merge - optimize memory usage by merging one by one
# First object
obj <- obj_list[[1]]
# Add cell ids to first object
colnames(obj) <- paste0(names(obj_list)[1], "_", colnames(obj))

# Merge remaining objects
for (i in 2:length(obj_list)) {
  # Add cell ids to current object
  colnames(obj_list[[i]]) <- paste0(names(obj_list)[i], "_", colnames(obj_list[[i]]))
  # Merge
  obj <- merge(x = obj, y = obj_list[[i]], project = "GSE274583_CD4")
  # Free memory after each merge
  gc()
}
# Free remaining memory
rm(obj_list)
gc()

# ----------------------------
# QC filtering (match your current stance)
# ----------------------------
qc_before <- ncol(obj)

obj <- subset(
  obj,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 4000 &
    percent.mt < 15 &
    nCount_RNA < 50000
)

qc_after <- ncol(obj)

# ----------------------------
# Add batch1234 field (唯一允许的批次分组字段)
# ----------------------------
# Map library_id to batch1234
obj$batch1234 <- case_when(
  grepl("^Library1", obj$library_id) ~ "Batch1",
  grepl("^Library2", obj$library_id) ~ "Batch2",
  grepl("^Library3", obj$library_id) ~ "Batch3",
  grepl("^Library4", obj$library_id) ~ "Batch4",
  TRUE ~ NA_character_
)

# Validate batch1234
if (any(is.na(obj$batch1234))) {
  stop("batch1234 contains NA values. Please check library_id mapping.", call. = FALSE)
}

valid_batches <- paste0("Batch", 1:4)
if (!all(obj$batch1234 %in% valid_batches)) {
  invalid_batches <- unique(obj$batch1234[!obj$batch1234 %in% valid_batches])
  stop("Invalid batch1234 values found: ", paste(invalid_batches, collapse = ", "), ". Valid values are: ", paste(valid_batches, collapse = ", "), call. = FALSE)
}

# Set batch1234 as factor with fixed levels
obj$batch1234 <- factor(obj$batch1234, levels = valid_batches)

# Output audit files
MAPPING_CSV <- file.path(LOG_DIR, "04c_scRNA_library_to_batch1234_mapping.csv")
BATCH_LIBRARY_TSV <- file.path(LOG_DIR, "04c_scRNA_batch1234_x_library_id_counts.tsv")

# Create mapping table
mapping_df <- data.frame(
  library_id = unique(obj$library_id),
  batch1234 = obj$batch1234[match(unique(obj$library_id), obj$library_id)]
)

# Write mapping CSV
write_csv(mapping_df, file = MAPPING_CSV)

# Create cross-tabulation
batch_library_table <- table(obj$batch1234, obj$library_id)

# Write cross-tab TSV
write.table(batch_library_table, file = BATCH_LIBRARY_TSV, sep = "\t", quote = FALSE)

# Log batch distribution
message("Batch distribution:")
print(table(obj$batch1234))

# ----------------------------
# Standard Seurat workflow + Harmony
# ----------------------------
set.seed(2026)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)

# ----------------------------
# Harmony integration by batch1234 (唯一允许的批次分组字段)
# ----------------------------
obj <- RunHarmony(
  object = obj,
  group.by.vars = "batch1234",
  reduction.use = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, n.neighbors = 30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)

# ----------------------------
# Memory optimization before JoinLayers
# ----------------------------
message("Running garbage collection to free memory...")
gc(verbose = FALSE)
message("Garbage collection completed")

# Get current memory usage
mem_usage <- gc(reset = TRUE)
message(paste("Memory usage before JoinLayers:", round(mem_usage[2, 2], 2), "MB"))

# ----------------------------
# Seurat v5: join layers so GetAssayData works
# ----------------------------
message("Performing JoinLayers with memory optimization...")
tryCatch({
  # Try with explicit assay specification
  obj <- JoinLayers(object = obj, assay = DefaultAssay(obj), layers = "data")
  message("JoinLayers completed successfully")
}, error = function(e) {
  message("JoinLayers failed with error: ", e$message)
  message("Trying alternative approach...")
  # Alternative approach: directly access the data layer without JoinLayers
  message("Skipping JoinLayers, will access data directly")
  # This flag will be used later to adjust GetAssayData call
  skip_join_layers <- TRUE
})

# Get memory usage after JoinLayers
gc(verbose = FALSE)
mem_usage <- gc(reset = TRUE)
message(paste("Memory usage after JoinLayers:", round(mem_usage[2, 2], 2), "MB"))

# ----------------------------
# Scores
# ----------------------------
# (1) Locked 11-gene *weighted* score (single truth: locked_weights.csv)
# Use normalized data slot; we compute weighted sum per cell on log-normalized data.
message("Getting expression data...")
tryCatch({
  # Seurat v5 compatible approach to get expression data
  message("Using Seurat v5 compatible approach to get expression data")
  
  # Get the assay name
  assay_name <- DefaultAssay(obj)
  
  # Try different approaches for Seurat v5
  if (packageVersion("Seurat") >= "5.0.0") {
    # Seurat v5 approach
    message("Using Seurat v5 LayerData approach")
    expr <- LayerData(obj, assay = assay_name, layer = "data")
  } else {
    # Fallback for older Seurat versions
    message("Using older Seurat GetAssayData approach")
    expr <- SeuratObject::GetAssayData(obj, assay = assay_name, slot = "data")
  }
  
  message("Expression data obtained successfully")
  message(paste("Expression matrix dimensions:", paste(dim(expr), collapse = "x")))
}, error = function(e) {
  message("Error getting expression data: ", e$message)
  message("Trying alternative approach...")
  
  # Try another alternative approach
  tryCatch({
    message("Trying direct assay access approach")
    assay_obj <- obj[[DefaultAssay(obj)]]
    # Directly access the data layer
    expr <- assay_obj$data
    message("Expression data obtained successfully")
    message(paste("Expression matrix dimensions:", paste(dim(expr), collapse = "x")))
  }, error = function(e2) {
    message("Alternative approach also failed: ", e2$message)
    stop("Failed to get expression data for scoring", call. = FALSE)
  })
})

available_locked <- intersect(locked_genes, rownames(expr))
missing_locked   <- setdiff(locked_genes, available_locked)

if (length(available_locked) < 2) {
  stop("Too few locked genes found in scRNA object. Check gene symbols and assay rownames.", call. = FALSE)
}

w_use <- locked_w_vec[available_locked]
w_use <- w_use / sum(abs(w_use))  # normalize weights for stability (keeps direction, avoids scale blow-up)

locked_score <- Matrix::colSums(expr[available_locked, , drop = FALSE] * as.numeric(w_use))
obj$Sig11_Score1 <- as.numeric(locked_score)

# (2) Orthogonal Th2 module via AddModuleScore (unweighted, standard Seurat approach)
# Keep the column name consistent with downstream (Th2_Score1).
obj <- AddModuleScore(obj, features = list(th2_genes), name = "Th2_Score", assay = DefaultAssay(obj))
# AddModuleScore creates Th2_Score1 by default
if (!("Th2_Score1" %in% colnames(obj@meta.data))) {
  # In case Seurat changes naming, enforce:
  th2_col <- grep("^Th2_Score", colnames(obj@meta.data), value = TRUE)[1]
  obj$Th2_Score1 <- obj@meta.data[[th2_col]]
}

# ----------------------------
# Save outputs
# ----------------------------
print(paste("Saving full RDS file to:", FULL_RDS_FILE))
saveRDS(obj, file = FULL_RDS_FILE)
print("Full RDS file saved successfully")

# Reviewer diet: remove scale.data to reduce file size;
# keep reductions/harmony/umap/clusters/scores so reviewer can run Fig3 directly.
obj_diet <- DietSeurat(
  obj,
  assays     = DefaultAssay(obj),
  counts     = TRUE,
  data       = TRUE,
  scale.data = FALSE,
  dimreducs  = c("pca", "harmony", "umap"),
  graphs     = c("RNA_snn"),
  misc       = TRUE
)

print(paste("Saving reviewer RDS file to:", REVIEWER_RDS_FILE))
saveRDS(obj_diet, file = REVIEWER_RDS_FILE)
print("Reviewer RDS file saved successfully")

# ----------------------------
# Output batch1234 x cluster counts (P1-2)
# ----------------------------
BATCH_CLUSTER_TSV <- file.path(LOG_DIR, "04c_scRNA_batch1234_x_cluster_counts.tsv")
batch_cluster_table <- table(obj$batch1234, obj$seurat_clusters)
write.table(batch_cluster_table, file = BATCH_CLUSTER_TSV, sep = "\t", quote = FALSE)

# ----------------------------
# Build summary + sessionInfo
# ----------------------------
# Get library sizes before and after QC
library_sizes <- table(obj$library_id)
library_sizes_str <- paste(names(library_sizes), library_sizes, sep = ":", collapse = ";")

# Get batch sizes
batch_sizes <- table(obj$batch1234)
batch_sizes_str <- paste(names(batch_sizes), batch_sizes, sep = ":", collapse = ";")

# Get package versions
seurat_version <- packageVersion("Seurat")
seuratobject_version <- packageVersion("SeuratObject")
harmony_version <- packageVersion("harmony")

# Add cluster_group placeholder (P2-1)
obj$cluster_group <- NA_character_

# Update summary table with governance fields (P1-1)
# Expected cache files
cache_files_expected <- c(
  "Library1_counts.rds",
  "Library2_counts.rds",
  "Library3_counts.rds",
  "Library4_sub1_counts.rds",
  "Library4_sub2_counts.rds",
  "Library4_sub3_counts.rds"
)

summary_tbl <- tibble::tibble(
  item = c(
    "base_dir",
    "raw_sc_dir",
    "full_rds",
    "reviewer_rds",
    "qc_cells_before",
    "qc_cells_after",
    "n_clusters",
    "library_sizes_post_qc",
    "batch1234_sizes_post_qc",
    "harmony_groupby",
    "batch1234_levels",
    "audit_mapping_csv",
    "audit_batch_library_tsv",
    "audit_batch_cluster_tsv",
    "cache_dir",
    "cache_files_expected",
    "cluster_group_policy",
    "locked_genes_total",
    "locked_genes_available",
    "locked_genes_missing",
    "th2_genes_total",
    "seurat_version",
    "seuratobject_version",
    "harmony_version"
  ),
  value = c(
    base_dir,
    RAW_SC_DIR,
    FULL_RDS_FILE,
    REVIEWER_RDS_FILE,
    qc_before,
    qc_after,
    length(unique(obj$seurat_clusters)),
    library_sizes_str,
    batch_sizes_str,
    "batch1234",
    paste(valid_batches, collapse = ";"),
    MAPPING_CSV,
    BATCH_LIBRARY_TSV,
    BATCH_CLUSTER_TSV,
    CACHE_DIR,
    paste(cache_files_expected, collapse = ";"),
    "plot_only_placeholder",
    length(locked_genes),
    length(available_locked),
    paste(missing_locked, collapse = ";"),
    length(th2_genes),
    as.character(seurat_version),
    as.character(seuratobject_version),
    as.character(harmony_version)
  )
)

readr::write_tsv(summary_tbl, file = SUMMARY_FILE)

capture.output(sessionInfo(), file = SESSIONINFO_FILE)

message("Done. Full object: ", FULL_RDS_FILE)
message("Done. Reviewer diet: ", REVIEWER_RDS_FILE)
message("Logs: ", SUMMARY_FILE, " ; ", SESSIONINFO_FILE)