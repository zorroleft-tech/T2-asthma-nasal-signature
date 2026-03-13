# WARNING: Do not run in routine pipeline. Requires manual curation.
# =========================================================================
# Script: 00b_candidate_annotation.R
# Purpose: Add biological annotations and perform quality check
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/raw/01_candidate_genes_passed.csv  # Genes passing PubMed threshold
#   - data/raw/01_candidate_genes_failed.csv  # Genes below PubMed threshold
#
# Outputs:
#   - data/raw/candidate_genes_annotated_frozen.csv  # Annotated candidate genes
#   - data/raw/02_Supplementary_Table_S1_Candidate_Genes.csv  # Supplementary table
# =========================================================================

library(dplyr)
library(readr)
library(tidyr)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
DATA_RAW_DIR <- file.path(base_dir, "data", "raw")

# Load candidates from Step 1
candidates <- read_csv(file.path(DATA_RAW_DIR, "01_candidate_genes_passed.csv"))

# Log base directory and paths
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Raw directory: ", DATA_RAW_DIR, "\n\n"))

cat(sprintf("Annotating %d candidate genes...\n", nrow(candidates)))
cat("(If actual count differs, check Script 1 output)\n\n")

# ============================================================================
# Complete Gene Annotation Template
# ============================================================================

gene_annotations <- tribble(
  ~Gene,       ~Primary_Pathway,       ~Secondary_Pathway,         ~Drug_Target,              ~Protein_Level_Biomarker,
  #------------|---------------------|-------------------------|--------------------------|--------------------------|
  
  # ===== IL-4/IL-13 axis (8 genes) =====
  "IL4",       "IL-4/13 axis",        "Th2 cytokine",           "Dupilumab",              "Serum",
  "IL13",      "IL-4/13 axis",        "Th2 cytokine",           "Lebrikizumab",           "Serum",
  "IL4R",      "IL-4/13 axis",        "Receptor",               "Dupilumab",              "Tissue",
  "IL13RA1",   "IL-4/13 axis",        "Receptor",               "None",                   "Tissue",
  "IL13RA2",   "IL-4/13 axis",        "Decoy receptor",         "None",                   "Tissue",
  "STAT6",     "IL-4/13 axis",        "Transcription factor",   "None",                   "Tissue",
  "CCL11",     "IL-4/13 axis",        "Eosinophil chemokine",   "Bertilimumab (disc.)",   "Serum",
  "CCL26",     "IL-4/13 axis",        "Eosinophil chemokine",   "None",                   "Tissue",
  
  # ===== IL-5/eosinophil (10 genes) =====
  "IL5",       "IL-5/Eosinophil",     "Th2 cytokine",           "Mepolizumab",            "Serum",
  "IL5RA",     "IL-5/Eosinophil",     "Receptor",               "Benralizumab",           "Tissue",
  "CLC",       "IL-5/Eosinophil",     "Eosinophil marker",      "None",                   "Tissue",
  "EPX",       "IL-5/Eosinophil",     "Eosinophil peroxidase",  "None",                   "Tissue",
  "RNASE2",    "IL-5/Eosinophil",     "Eosinophil-derived",     "None",                   "Serum (EDN)",
  "RNASE3",    "IL-5/Eosinophil",     "Eosinophil-derived",     "None",                   "Serum (ECP)",
  "PRG2",      "IL-5/Eosinophil",     "Major basic protein",    "None",                   "Tissue",
  "PRG3",      "IL-5/Eosinophil",     "Major basic protein",    "None",                   "Tissue",
  "SIGLEC8",   "IL-5/Eosinophil",     "Eosinophil receptor",    "Lirentelimab",           "Tissue",
  "CCR3",      "IL-5/Eosinophil",     "Chemokine receptor",     "None",                   "Tissue",
  
  # ===== IL-33/TSLP/ILC2 (6 genes) =====
  "IL33",      "IL-33/TSLP/ILC2",     "Alarmin",                "Etokimab (dev.)",        "Serum",
  "TSLP",      "IL-33/TSLP/ILC2",     "Epithelial cytokine",    "Tezepelumab",            "Serum",
  "IL1RL1",    "IL-33/TSLP/ILC2",     "ST2 receptor",           "Astegolimab (dev.)",     "Serum (sST2)",
  "AREG",      "IL-33/TSLP/ILC2",     "Epithelial repair",      "None",                   "Tissue",
  "IL17RB",    "IL-33/TSLP/ILC2",     "IL-25 receptor",         "None",                   "Tissue",
  "IL25",      "IL-33/TSLP/ILC2",     "Epithelial cytokine",    "None",                   "Serum",
  
  # ===== Transcription factors (5 genes) =====
  "GATA3",     "Transcription",       "Th2 master regulator",   "None",                   "Tissue",
  "GATA1",     "Transcription",       "Eosinophil development", "None",                   "Tissue",
  "STAT5A",    "Transcription",       "IL-5 signaling",         "None",                   "Tissue",
  "STAT5B",    "Transcription",       "IL-5 signaling",         "None",                   "Tissue",
  "IRF4",      "Transcription",       "Th2 differentiation",    "None",                   "Tissue",
  
  # ===== Mucus/epithelial (7 genes) =====
  "CLCA1",     "Mucus/Epithelial",    "Chloride channel",       "None",                   "Tissue",
  "MUC5AC",    "Mucus/Epithelial",    "Mucin",                  "None",                   "Sputum",
  "MUC5B",     "Mucus/Epithelial",    "Mucin",                  "None",                   "Sputum",
  "CST1",      "Mucus/Epithelial",    "Cystatin (protease inh.)","None",                  "Tissue",
  "CST2",      "Mucus/Epithelial",    "Cystatin (protease inh.)","None",                  "Tissue",
  "AGR2",      "Mucus/Epithelial",    "ER stress response",     "None",                   "Tissue",
  "FOXA3",     "Mucus/Epithelial",    "Mucus transcription",    "None",                   "Tissue",
  
  # ===== Chemokines (4 genes) =====
  "CCL17",     "Chemokine",           "Th2 recruitment (TARC)", "None",                   "Serum",
  "CCL22",     "Chemokine",           "Th2 recruitment (MDC)",  "None",                   "Serum",
  "CCL24",     "Chemokine",           "Eosinophil recruitment", "None",                   "Tissue",
  "CXCL8",     "Chemokine",           "Neutrophil chemokine",   "None",                   "Serum",
  
  # ===== Tissue remodeling (5 genes) =====
  "POSTN",     "Tissue remodeling",   "ECM protein",            "None",                   "Serum",
  "TGM2",      "Tissue remodeling",   "Transglutaminase",       "None",                   "Tissue",
  "SERPINB2",  "Tissue remodeling",   "PAI-2 (protease inh.)",  "None",                   "Tissue",
  "CHI3L1",    "Tissue remodeling",   "YKL-40 (chitinase-like)","None",                   "Serum",
  "SERPINB10", "Tissue remodeling",   "Bomapin (protease inh.)","None",                   "Tissue",
  
  # ===== Prostaglandin pathway (4 genes) =====
  "PTGDR2",    "PGD2 pathway",        "CRTH2 receptor",         "Fevipiprant (disc.)",    "Tissue",
  "PTGS2",     "PGD2 pathway",        "COX-2 (PG synthesis)",   "NSAIDs",                 "Tissue",
  "HPGDS",     "PGD2 pathway",        "PGD2 synthase",          "None",                   "Tissue",
  "ALOX15",    "PGD2 pathway",        "15-lipoxygenase",        "None",                   "Tissue",
  
  # ===== Histamine (4 genes) =====
  "HDC",       "Histamine",           "Histamine synthesis",    "None",                   "Tissue",
  "HRH1",      "Histamine",           "H1 receptor",            "Antihistamines",         "Tissue",
  "HRH2",      "Histamine",           "H2 receptor",            "H2 blockers",            "Tissue",
  "HRH4",      "Histamine",           "H4 receptor",            "None",                   "Tissue",
  
  # ===== Mast cell/basophil (6 genes) =====
  "CPA3",      "Mast cell",           "Carboxypeptidase A3",    "None",                   "Tissue",
  "TPSAB1",    "Mast cell",           "Tryptase alpha/beta",    "None",                   "Serum",
  "TPSB2",     "Mast cell",           "Tryptase beta",          "None",                   "Serum",
  "MS4A2",     "Mast cell",           "FcεRI beta chain",       "None",                   "Tissue",
  "FCER1A",    "Mast cell",           "FcεRI alpha chain",      "Omalizumab (IgE)",       "Tissue",
  "FCER1G",    "Mast cell",           "FcεRI gamma chain",      "None",                   "Tissue",
  
  # ===== IgE pathway (1 gene) =====
  "FCER2",     "IgE pathway",         "CD23 (low-affinity)",    "Quilizumab (disc.)",     "Serum (sCD23)"
)

cat(sprintf("✓ Annotation template loaded: %d genes defined\n", nrow(gene_annotations)))

# ============================================================================
# Merge annotations with PubMed counts
# ============================================================================

candidates_annotated <- candidates %>%
  left_join(gene_annotations, by = "Gene")

unannotated <- candidates_annotated %>%
  filter(is.na(Primary_Pathway))

if (nrow(unannotated) > 0) {
  cat("\n⚠ WARNING: The following genes passed Step 1 but lack annotations:\n")
  print(unannotated$Gene)
  cat("→ Please add them to gene_annotations tribble\n\n")
}

# ============================================================================
# Quality Checks
# ============================================================================

cat("\n[Quality Checks]\n")
cat("=" , rep("-", 60), "\n", sep = "")

# Check 1: Pathway coverage
pathway_coverage <- candidates_annotated %>%
  filter(!is.na(Primary_Pathway)) %>%
  group_by(Primary_Pathway) %>%
  summarize(N_genes = n()) %>%
  arrange(desc(N_genes))

cat("\n✓ Check 1: Pathway coverage\n")
print(pathway_coverage)

required_pathways <- c(
  "IL-4/13 axis",
  "IL-5/Eosinophil",
  "IL-33/TSLP/ILC2",
  "Mucus/Epithelial",
  "Tissue remodeling"
)

missing_pathways <- setdiff(required_pathways, pathway_coverage$Primary_Pathway)
if (length(missing_pathways) > 0) {
  cat("\n⚠ WARNING: Missing essential pathways:\n")
  cat(paste("  -", missing_pathways, collapse = "\n"), "\n")
} else {
  cat("  → All 5 essential pathways covered ✓\n")
}

# Check 2: Drug target representation
drug_targets <- candidates_annotated %>%
  filter(!is.na(Drug_Target) & Drug_Target != "None") %>%
  dplyr::select(Gene, Primary_Pathway, Drug_Target) %>%
  arrange(Primary_Pathway)

cat(sprintf("\n✓ Check 2: Drug targets identified: %d genes\n", nrow(drug_targets)))
if (nrow(drug_targets) > 0) {
  print(drug_targets)
}

# Check 3: Biomarker measurability
protein_biomarkers <- candidates_annotated %>%
  filter(Protein_Level_Biomarker %in% c("Serum", "Sputum")) %>%
  dplyr::select(Gene, Protein_Level_Biomarker)

cat(sprintf("\n✓ Check 3: Serum/sputum-measurable biomarkers: %d genes\n",
            nrow(protein_biomarkers)))

# Check 4: Platform compatibility estimation
common_platform_genes <- c(
  "IL13", "IL4", "IL5", "TSLP", "IL1RL1", "CLCA1", "POSTN", 
  "CLC", "CCL11", "CCL17", "CCL26", "PTGDR2", "SERPINB2", 
  "CPA3", "CST1", "GATA3", "STAT6"
)

platform_compatible <- candidates_annotated %>%
  filter(Gene %in% common_platform_genes) %>%
  dplyr::select(Gene)

cat(sprintf("\n✓ Check 4: Platform-compatible genes (GPL570/GPL6480/GPL96): %d genes\n",
            nrow(platform_compatible)))
cat("  (These genes have good probe coverage on major platforms)\n")

# Check 5: Cell type diversity
cell_type_markers <- list(
  Epithelial = c("CLCA1", "MUC5AC", "MUC5B", "CST1", "CST2", "AGR2", "FOXA3", "TSLP", "IL33", "IL25"),
  Th2_ILC2 = c("IL4", "IL13", "IL5", "GATA3", "IL1RL1", "PTGDR2", "CCL17"),
  Eosinophil = c("CLC", "EPX", "RNASE2", "RNASE3", "PRG2", "PRG3", "SIGLEC8", "CCR3", "IL5RA"),
  Mast_Basophil = c("CPA3", "TPSAB1", "TPSB2", "MS4A2", "FCER1A", "HDC"),
  Tissue_ECM = c("POSTN", "TGM2", "CHI3L1", "SERPINB2", "SERPINB10")
)

cat("\n✓ Check 5: Cell type marker coverage\n")
for (ct in names(cell_type_markers)) {
  markers_found <- sum(candidates_annotated$Gene %in% cell_type_markers[[ct]])
  total_markers <- length(cell_type_markers[[ct]])
  cat(sprintf("  - %-15s: %2d/%2d genes (%.0f%%)\n", 
              ct, markers_found, total_markers, 
              100 * markers_found / total_markers))
}

# ============================================================================
# Additional genes from pathway gaps (if needed)
# ============================================================================

borderline <- read_csv(file.path(DATA_RAW_DIR, "01_candidate_genes_failed.csv")) %>%
  filter(PubMed_Count >= 8)

if (nrow(borderline) > 0 & length(missing_pathways) > 0) {
  cat("\n[Optional] Borderline genes for pathway coverage:\n")
  print(borderline)
  cat("\nManually review if any belong to missing pathways\n")
}

# ============================================================================
# Export annotated candidate list
# ============================================================================

cat("\n[Export Results]\n")
cat("=" , rep("-", 60), "\n", sep = "")

write_csv(
  candidates_annotated,
  file.path(DATA_RAW_DIR, "candidate_genes_annotated_frozen.csv")
)

cat("\n✓ Annotated candidate list exported\n")
cat(sprintf("  → File: %s/candidate_genes_annotated_frozen.csv\n", DATA_RAW_DIR))
cat(sprintf("  → Final candidate pool: %d genes\n", nrow(candidates_annotated)))

# ============================================================================
# Generate summary table for paper
# ============================================================================

summary_table <- candidates_annotated %>%
  filter(!is.na(Primary_Pathway)) %>%
  dplyr::select(
    Gene,
    Primary_Pathway,
    Secondary_Pathway,
    PubMed_Count,
    Drug_Target,
    Protein_Level_Biomarker
  ) %>%
  arrange(Primary_Pathway, desc(PubMed_Count))

write_csv(
  summary_table,
  file.path(DATA_RAW_DIR, "02_Supplementary_Table_S1_Candidate_Genes.csv")
)

cat("\n✓ Supplementary Table S1 generated\n")
cat(sprintf("  → File: %s/02_Supplementary_Table_S1_Candidate_Genes.csv\n", DATA_RAW_DIR))
cat("  (This will be included in the paper)\n\n")

# ============================================================================
# Summary Statistics for Methods Section
# ============================================================================

cat("\n[Summary Statistics for Methods]\n")
cat("=" , rep("-", 60), "\n", sep = "")

cat(sprintf("\n- Total candidate genes: %d\n", nrow(candidates_annotated)))
cat(sprintf("- Pathways represented: %d\n", nrow(pathway_coverage)))
cat(sprintf("- Drug-targetable genes: %d (%.1f%%)\n", 
            nrow(drug_targets),
            100 * nrow(drug_targets) / nrow(candidates_annotated)))
cat(sprintf("- Serum/sputum biomarkers: %d (%.1f%%)\n",
            nrow(protein_biomarkers),
            100 * nrow(protein_biomarkers) / nrow(candidates_annotated)))

# ============================================================================
# Add Detailed Drug Target Information
# ============================================================================

cat("\n[Adding Detailed Drug Target Information]\n")
cat("=" , rep("-", 60), "\n", sep = "")

# Create drug target database
drug_targets <- data.frame(
  Gene = c("IL4R", "SIGLEC8", "CCR3", "STAT6"),
  Drug_Name = c(
    "Dupilumab",
    "Lirentelimab",
    "Not in clinical trials",
    "Not in clinical trials"
  ),
  Drug_Status = c(
    "FDA approved (2017)",
    "Phase 3 clinical trial",
    "Research stage",
    "Research stage"
  ),
  Clinical_Indication = c(
    "Moderate-to-severe asthma, atopic dermatitis",
    "Eosinophilic disorders",
    "N/A",
    "N/A"
  ),
  stringsAsFactors = FALSE
)

# Merge to annotation file
candidates_annotated <- candidates_annotated %>%
  left_join(drug_targets, by = "Gene") %>%
  mutate(
    Drug_Target_Detailed = ifelse(!is.na(Drug_Name) & Drug_Status != "Research stage",
                                 paste0(Drug_Name, " (", Drug_Status, ")"),
                                 Drug_Target)
  )

# Update Drug_Target column with detailed information
candidates_annotated <- candidates_annotated %>%
  mutate(Drug_Target = Drug_Target_Detailed) %>%
  dplyr::select(-Drug_Target_Detailed)

# Re-export annotated candidate list
write_csv(
  candidates_annotated,
  file.path(DATA_RAW_DIR, "candidate_genes_annotated_frozen.csv")
)

cat(sprintf("\n✓ Drug target information updated: %d genes with detailed drug data\n",
            nrow(drug_targets)))

cat("\n" , rep("=", 70), "\n", sep = "")
cat("✓ Script 2 Complete: Candidate pool annotation finished\n")
cat("→ Next: Proceed to Phase 2 (differential expression analysis)\n\n")
