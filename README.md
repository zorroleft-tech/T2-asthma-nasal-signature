# Biomarker-Anchored Validation of a Locked Nasal Immune Signature

## Project Overview
This repository contains the complete analytical R scripts, auditable run-all pipeline, and intermediate matrices supporting our manuscript submitted to the ***Journal of Translational Medicine***.

The project aims to provide a transparent, *in silico* validation workspace with "industrial-grade reproducibility", "absolute data leakage prevention", and "extremely low memory consumption." By enforcing mathematically locked weights and prespecified reporting boundaries, this pipeline serves as a reproducible blueprint for non-invasive T2-high asthma stratification, specifically aligning with data-driven clinical decision processes.

## Core Research Stance (Audit-Ready Methodology)
- **Anchored ROC/AUC Only**: We report diagnostic discrimination (ROC/AUC) *exclusively* in cohorts with independent clinical biomarker anchors (IgE, FeNO, sputum phenotypes).
- **Replication-only Boundary**: In external cohorts lacking independent clinical labels, we enforce a strict "Replication-only" strategy. We report only molecular consistency (Cohen's d, Wilcoxon P) and strictly prohibit calculating or reporting AUC to completely eliminate circular validation and cohort-specific refitting bias.
- **Locked Weights**: All external validations are forced to read a single `locked_weights.csv` file. No algorithmic refitting or dynamic thresholding is allowed across distinct platforms.
- **Prespecified Signature-Size Rule**: 70% Elastic Net stability selection threshold; retain all genes if 10–18; if 19–25 retain top 18 by frequency; if <10 stepwise lower threshold to ensure ≥8; all decisions log-documented.

## Directory Structure

```
asthma11gene_allergy/
├── README.md                  # Core instructions: running guide, environment requirements, and figure correspondences
├── run_all.R                  # Only master script: sequentially call each analysis module using independent subprocesses
├── renv.lock                  # Environment lock file: lock R package versions
│
├── 📂 data_preparation/       # [Phase 0: Raw material processing plant] (extract minimal matrix, prevent memory explosion, output logs)
│   ├── 00a_ingest_expr.R                    # Read expression matrices under raw/ and perform initial data ingestion
│   ├── 00b_id_standardize_and_map.R         # Standardize identifiers and complete probe to Gene Symbol mapping
│   ├── 00c_prep_pheno_anchors_65204_201955_45111.R  # Prepare phenotype data for anchor cohorts
│   ├── 00d_gse201955_build_sample_map.R     # Build sample mapping for GSE201955
│   ├── 01_prep_discovery_GSE152004.R         # Extract full discovery cohort data
│   ├── 02_prep_anchors_65204_201955_45111.R  # Extract three truth cohorts (retain continuous variables like FeNO/Eos for subsequent splitting)
│   ├── 03_prep_rep_7_cohorts.R               # Extract 7 directional validation cohorts
│   ├── 04a_extract_pheno_from_raw.R         # Extract phenotype info from raw Series Matrix files using GEOquery and generate standardized tables
│   ├── 04b_prep_gse230048_pheno.R           # Prepare phenotype data for GSE230048
│   ├── 04c_prep_scRNA_GSE274583_diet.R      # Single-cell extreme slimming (only keep model genes + dozens of Canonical Markers + UMAP)
│   ├── 05_calculate_minimal_stats.R          # Calculate minimal statistical closure for each cohort
│   ├── 06_generate_audit_report.R           # Generate field audit reports for raw data sources
│   ├── 07_summarize_deliverables.R          # Summarize all cohort results and prepare deliverables
│   ├── check_rds_files.R                      # Check generated RDS files
│   ├── check_generated_files.R               # Check generated files for validation
│   ├── gse_strategy.tsv                      # GSE processing strategy table
│   ├── platform_strategy.tsv                # Platform processing strategy table
│   └── validate_outputs_only.R                # Validate output files
│
├── 📂 data/                   # [Data center: strictly follow one-way flow from raw to processed to derived]
│   ├── 📂 raw/                # (Static input: 14 raw cohort files, GPL files, 60 candidate gene table, pathway/drug databases, etc.)
│   ├── 📂 processed_full/     # (Storage for full transcriptome matrices)
│   ├── 📂 processed_diet/     # (Storage for minimal model-gene matrices)
│   ├── 📂 processed/          # (Storage for processed data including phenotype tables and stats)
│   │   └── 📂 phase3_outputs/ # (Outputs from Phase 3 phenotype processing workflow)
│   └── 📂 derived/            # [Core moat: feature selection results]
│       ├── gene_panel_final.csv  # Final gene panel
│       └── locked_weights.csv # Locked coefficients for model genes, external validation must and only read this file!
│
├── 📂 R/                      # [Underlying function arsenal]
│   ├── 00_theme_setup.R       # Global unified plotting theme, color palette (ensure publication-quality figures)
│   ├── 01_scoring_logic.R     # Strict scoring function based on locked_weights
│   ├── 02_stats_funcs.R       # Cohen's d / Wilcoxon calculation encapsulation (AUC calculation module needs physical isolation)
│   ├── 03_index_logger.R      # Log function that automatically records all output figures to OUTPUTS_INDEX.csv
│   ├── 04_auc_anchor_only.R   # AUC calculation functions - strictly for anchor cohorts only
│   └── 📂 figures/            # [Figure generation functions]
│       ├── fig2A_logFC_concordance.R  # Generate Fig2A - log2FC concordance plot
│       ├── fig2B_circular_heatmap.R   # Generate Fig2B - circular heatmap
│       └── fig2C_attenuation.R        # Generate Fig2C - tissue gradient attenuation plot
│
├── 📂 analysis/               # [Core analysis area: globally sequential numbered 01~20, highly mapped to main paper]
│   │
│   ├── 📂 01_feature_selection/        # [Model training area: discovery cohort only, not run by default]
│   │   ├── 00a_pubmed_literature_mining.R  # PubMed literature mining (API-dependent)
│   │   ├── 00b_candidate_annotation.R      # Candidate gene annotation (static)
│   │   ├── 00c_literature_sensitivity.R    # Literature sensitivity analysis
│   │   ├── 01_t2_labeling.R            # Define T2-high/low (Woodruff median split)
│   │   ├── 02_limma_filtering.R        # Candidate genes logFC > 1, FDR < 0.001 filtering
│   │   ├── 03_stability_selection.R    # 1000 Bootstrap + ElasticNet (>70% threshold)
│   │   └── 04_finalize_model.R         # Lock final model genes, export to data/derived/
│   │
│   ├── 📂 02_main_validation/          # [Main paper figure area: no parameter adjustment allowed, forced to read locked weights]
│   │   ├── 05_table1_baseline.R        # Output Table 1 (integrate cohort baseline data)
│   │   ├── 06_table2_model_coef.R      # Output Table 2 (locked-signature weight table)
│   │   ├── 07a_fig1_discovery_common_pathways.R  # Output Fig 1B (Z-score density plot), generate common pathways
│   │   ├── 07b_fig1_discovery_network.R           # Output Fig 1C (KEGG network)
│   │   ├── 08_fig2_replication.R       # Output Fig 2A(log2FC), 2B(heatmap), 2C(attenuation effect), 2D(118761 paired correlation), 2E(40888 PBMC boxplot)
│   │   ├── 09_fig3_scRNA.R             # Output Fig 3A(UMAP), 3B(model-gene localization), 3C(Th2 module), 3D(violin plot), 3E(bubble plot)
│   │   └── 10_fig4_anchors_trans.R     # Output Fig 4A(IgE scatter), 4B(side-by-side three truth ROC puzzle), 4C(DCA), 4D(DGIdb interaction)
│   │
│   └── 📂 03_supplementary/            # [Supplementary figure and table area: detailed supporting materials]
│       ├── 11_supp_tables_S1_to_S4.R   # Output Table S1(QC), S2(cell/function), S3(65204 detailed table), S4(7-cohort replication consistency)
│       ├── 12_supp_tables_S5_to_S7.R   # Output Table S5(NRI/IDI), S6(GSEA full table), S7(drug interaction full table)
│       ├── 13a_supp_tableS8_age.R      # Output Table S8 (extremely critical age subgroup defense table: only Cohen's d)
│       ├── 13b_tableS8_age_stratum_extract.R  # Extract age strata from series matrix files
│       ├── 14_supp_tableS9_platform.R  # Output Table S9 (platform probe coverage matrix)
│       ├── 15_supp_figS1_qc.R          # Output Fig S1 (discovery cohort volcano plot, PCA)
│       ├── 15_supp_figS1_qc_fixed.R    # Fixed version of Fig S1
│       ├── 16_supp_figS2_landscape.R   # Output Fig S2 (8-cohort extended heatmap, coverage tile plot)
│       ├── 17_supp_figS3_scRNA_qc.R    # Output Fig S3 (single-cell Library QC and Harmony correction)
│       ├── 18_supp_figS4_architecture.R# Output Fig S4 (IgE correlation Bootstrap forest plot, cross-cohort module consistency)
│       ├── 19_supp_figS5_robustness.R  # Output Fig S5 (dose-response, calibration curves, permutation test distribution)
│       └── 20_supp_figS6_mechanistic.R # Output Fig S6 (5 representative genes IgE-stratified expression boxplots)
│
└── 📂 output/                 # [Finished product delivery area]
    ├── OUTPUTS_INDEX.csv      # Core tracking table: automatically record each generated figure information by 03_index_logger.R
    ├── 📂 logs/               # Store log files from all scripts
    ├── 📂 figures_main/       # (Store Fig 1 to Fig 4)
    ├── 📂 tables_main/        # (Store Table 1 and Table 2)
    └── 📂 supplement/         # (Store S1-S9 tables, S1-S6 supplementary figures)
```

## Running Guide

### Environment Requirements
- Recommended / tested: R 4.4.3
- Minimum: R ≥ 4.2.0
- Required R packages: tidyverse, limma, glmnet, Seurat, patchwork, cowplot, ggpubr, dplyr, ggplot2
- Use `renv` to manage package versions to ensure reproducibility

### Running Process
1. **Data Preparation**: First run scripts in the `data_preparation` directory to generate:
   - Full transcriptome matrices in `data/processed_full/`
   - Minimal model-gene matrices in `data/processed_diet/`
   - Standardized phenotype data in `data/processed/phase3_outputs/`
2. **Main Analysis**: Run the `run_all.R` script, which will sequentially call each analysis module
3. **Result Viewing**: All outputs will be saved in the `output` directory

### Script Execution Requirements
- **Working Directory**: All scripts must be run from the project root directory
- **Path Validation**: Scripts will automatically validate the working directory structure
- **Log Files**: All scripts generate log files in `output/logs/` directory
- **Progress Output**: Scripts provide step-by-step progress updates via `cat()` statements
- **Error Handling**: Scripts include file existence checks and appropriate error messages

### Running Options
- **Reviewer mode (default)**: `RUN_PREP = FALSE`, `RUN_TRAINING = FALSE`, `ALLOW_NET = FALSE` - No internet access required, no GPL downloads needed, directly reads processed data from `data/processed_diet/` directory (model-gene minimal matrices) (and scRNA diet products if applicable)
- **Full/Training mode (optional)**: Can read from `data/processed_full/` for training, troubleshooting, and supplementary analyses, but this is not part of the default reviewer input
- **From-scratch mode (optional)**: Set `RUN_PREP = TRUE` or use `--run-prep` command line argument - Requires internet access/proxy, downloads and processes GPL platform annotation files
- **Training mode (optional)**: Set `RUN_TRAINING = TRUE` or use `--run-training` command line argument - Executes feature selection and model retraining
- **Network access (optional)**: Set `ALLOW_NET = TRUE` or use `--allow-net` command line argument - Enables network access for API calls

## New Phenotype Processing Workflow

### Phase 3: Raw Phenotype Extraction and Standardization

This new workflow ensures reproducible, traceable phenotype extraction from raw data sources with strict validation gates, using GEOquery for reliable data access.

#### Key Scripts and Their Functions

1. **04a_extract_pheno_from_raw.R**
   - **Purpose**: Extract phenotype information from raw data sources using GEOquery
   - **Features**:
     - Direct parsing of `!Sample_characteristics_ch1` lines from Series Matrix files
     - Per-sample characteristics aggregation (each GSM gets its own `raw_characteristics_all`)
     - Key-value long table output (`*_pheno_kv_long.rds`)
     - Cohort-specific mapping functions for anchor cohorts (GSE65204, GSE201955, GSE45111)
     - Normalization of key names for consistent parsing
   - **Outputs**:
     - `data/processed/phase3_outputs/*_pheno_samples.rds` (sample-level table)
     - `data/processed/phase3_outputs/*_pheno_kv_long.rds` (key-value long table)
     - `data/processed/phase3_outputs/*_pheno_raw.rds` (final standardized table)

2. **04b_prep_gse230048_pheno.R**
   - **Purpose**: Prepare phenotype data for GSE230048
   - **Features**:
     - Specialized processing for GSE230048 cohort
     - Phenotype data standardization and validation
   - **Outputs**:
     - Standardized phenotype data for GSE230048

3. **04c_prep_scRNA_GSE274583_diet.R**
   - **Purpose**: Single-cell extreme slimming
   - **Features**:
     - Only keep model genes + dozens of Canonical Markers + UMAP
     - Memory optimization for single-cell data
   - **Outputs**:
     - Slimmed single-cell data for GSE274583

4. **05_calculate_minimal_stats.R**
   - **Purpose**: Calculate minimal statistical closure for each cohort, including sample counts, effect sizes, and p-values
   - **Features**:
     - Enforces validation gates: n_case > 0 && n_control > 0 for binary cohorts
     - n_AA > 0 && n_NA > 0 && n_HC > 0 for GSE40888
     - Calculates sample distributions and basic statistics
     - Placeholder for effect size and p-value calculation (requires expression data)
   - **Outputs**:
     - `data/processed/phase3_outputs/minimal_stats_report.csv` (summary of all cohort stats)

5. **06_generate_audit_report.R**
   - **Purpose**: Generate field audit reports for raw data sources
   - **Features**:
     - Field-level audit: non-empty rate, unique values, grouping potential
     - Analysis of both Series Matrix and SOFT files
     - Merges field information from multiple sources with priority to SOFT files
     - Identifies fields that can be used for grouping
   - **Outputs**:
     - `data/processed/phase3_outputs/field_audit_report.csv` (comprehensive field audit report)

6. **07_summarize_deliverables.R**
   - **Purpose**: Summarize all cohort results and prepare deliverables
   - **Features**:
     - Compiles information from phenotype tables, stats reports, and audit reports
     - Provides comprehensive cohort summary including sample counts, group distributions, and available fields
     - Identifies cohorts with complete vs. missing information
   - **Outputs**:
     - `data/processed/phase3_outputs/deliverables_summary.csv` (summary of all deliverables)

### Anchor Cohort Validation

All three anchor cohorts have passed validation with proper case/control mapping:

- **GSE65204 (IgE anchor)**: case=34, control=35 (median split on lnige)
- **GSE201955 (FeNO anchor)**: case=48, control=26 (direct mapping from asthma field)
- **GSE45111 (Sputum anchor)**: case=17, control=30 (phenotype-based mapping)

### Validation Gates

- **Binary cohorts**: Must have n_case > 0 AND n_control > 0
- **GSE40888 (AA/NA/HC)**: Must have all three groups present
- **Traceability**: All group labels must have a clear source in `std_group_label_source`
- **Reproducibility**: All parsing is based on raw files, no manual curation
- **Consistency**: Group labels must match predefined case/control levels exactly

### Execution Flow

1. **Extract Phenotypes**: Run `04_extract_pheno_from_raw.R` to extract and standardize phenotype data from raw files
2. **Calculate Stats**: Run `05_calculate_minimal_stats.R` to compute sample distributions and basic statistics
3. **Generate Audit Reports**: Run `06_generate_audit_report.R` to analyze field quality and grouping potential
4. **Summarize Deliverables**: Run `07_summarize_deliverables.R` to compile all results into a comprehensive summary

### Key Improvements

- **Reliable Data Access**: Uses GEOquery for consistent data retrieval
- **Traceability**: All group labels include clear source information
- **Standardization**: Normalized key names and consistent field parsing
- **Validation**: Enforces validation gates to ensure data quality
- **Audit Trail**: Comprehensive field audit reports for transparency
- **Reproducibility**: All parsing based on raw files with no manual curation


## Figure Correspondences

### Main Paper Figures
- **Fig 1**: Discovery cohort analysis
  - Fig 1B: Z-score density plot (generated by 07a_fig1_discovery_common_pathways.R)
  - Fig 1C: KEGG network (generated by 07b_fig1_discovery_network.R)
- **Fig 2**: External cohort replication
  - Fig 2A: log2FC plot
  - Fig 2B: Heatmap
  - Fig 2C: Attenuation effect
  - Fig 2D: GSE118761 paired correlation plot
  - Fig 2E: GSE40888 PBMC boxplot
- **Fig 3**: Single-cell analysis
  - Fig 3A: UMAP
  - Fig 3B: Model-gene localization
  - Fig 3C: Th2 module
  - Fig 3D: Violin plot
  - Fig 3E: Bubble plot
- **Fig 4**: Clinical anchor validation
  - Fig 4A: IgE scatter plot
  - Fig 4B: Side-by-side three truth ROC puzzle
  - Fig 4C: DCA
  - Fig 4D: DGIdb interaction

### Supplementary Figures
- **Table S1**: QC table
- **Table S2**: Cell/function table
- **Table S3**: GSE65204 detailed table
- **Table S4**: 7-cohort replication consistency table
- **Table S5**: NRI/IDI table
- **Table S6**: GSEA full table
- **Table S7**: Drug interaction full table
- **Table S8**: Age subgroup defense table
- **Table S9**: Platform probe coverage matrix
- **Fig S1**: Discovery cohort volcano plot and PCA
- **Fig S2**: 8-cohort extended heatmap and coverage tile plot
- **Fig S3**: Single-cell Library QC and Harmony correction
- **Fig S4**: IgE correlation Bootstrap forest plot and cross-cohort module consistency
- **Fig S5**: Dose-response, calibration curves, and permutation test distribution
- **Fig S6**: 5 representative genes IgE-stratified expression boxplots

## Code Disciplines

1. **Feature Selection Isolation and Default Skip**: run_all.R must not execute scripts in the analysis/01_feature_selection/ directory by default. The main process must directly and only read data/derived/locked_weights.csv.
2. **"Zero AUC" Red Line for Non-anchor Cohorts**: When processing non-anchor cohorts, strictly prohibit importing or using packages like pROC, ROCR, and strictly prohibit calculating any AUC, Sensitivity, or Specificity.
   - AUC/ROC is only allowed in scripts for anchor cohorts and only for Fig4-related outputs
   - Explicit whitelist of anchor cohorts:
     - GSE65204 (IgE anchor)
     - GSE201955 (FeNO anchor)
     - GSE45111 (sputum phenotype anchor)
   - For all other cohorts: only Cohen's d / Wilcoxon P / directional consistency and other replication-only metrics are allowed
3. **Mandatory Statement for Fig 4B Three Truth ROC**: When plotting ROC curves for GSE65204 (IgE), GSE201955 (FeNO), GSE45111 (Sputum eos), the code comments and figure title must clearly state that this is a "Parallel display of distinct clinical anchors", and strictly prohibit calculating a pooled AUC.
4. **Strict Positioning for Fig 2D and 2E**: Fig 2D must plot a paired correlation scatter plot or paired difference plot of the same subject's nose-bronchus; Fig 2E plots a boxplot of AA/NA/HC groups in PBMC samples, labeling Cohen's d and Wilcoxon P, strictly prohibiting adding diagnostic threshold lines.
5. **Extreme Data Slimming and Memory Control**: The submission/reproduction default pipeline (reviewer mode) strictly uses diet data (model-gene minimal matrices) and does not require full transcriptome matrices. Analysis scripts must include garbage collection code at the end of their logic.
   - `processed_full/` is allowed as internal training/troubleshooting input but is not part of the default reviewer pipeline
6. **Forced Independent Process Execution**: run_all.R must use system2("Rscript", "script_path") or callr::r_script() to call analysis scripts, ensuring each run executes in a clean R environment.
7. **Path Standardization**: All scripts must use consistent path handling:
   - Base directory detection: `normalizePath(getwd(), winslash = "/", mustWork = TRUE)`
   - Path construction: `file.path()` for cross-platform compatibility
   - Output directory creation: `dir.create()` with `recursive = TRUE`
   - Log files: Stored in `output/logs/<script_name>.log`
   - No hardcoded absolute paths allowed
   - No automatic `setwd()` calls
   - No scattered path logic

8. **Script Header Standardization**: All scripts must include a standardized header with the following sections:
   - Script: Name and brief description
   - Purpose: Detailed description of the script's purpose
   - Author: Script author
   - Date: Creation or last modification date
   - Inputs: List of input files and parameters
   - Outputs: List of output files and results

9. **Execution Flow Standardization**: All scripts must follow a consistent execution flow:
   - Use `cat()` for step-by-step progress output
   - Output dimension information after loading files
   - Log key parameters and settings
   - Perform file existence checks before processing
   - Provide clear completion messages
   - Use `sink()` to properly close log files at the end

## Data Security
- All raw data must be stored in the data/raw/ directory, and only necessary information should be extracted during processing
- Processed data must be slimmed, only retaining the minimal dataset required for analysis
- Strictly prohibit hardcoding any sensitive information in code
- All output files must be checked to ensure they do not contain raw data or sensitive information
- `processed_full/` directory should not be included in reviewer deliverables

## Memory Optimization
- Use streaming processing during data processing to avoid loading all data at once
- Immediately execute garbage collection after analysis to release memory
- Use DietSeurat() function for single-cell analysis to reduce memory usage
- Use chunk processing for large matrix operations

## Contact
For any questions, please contact the project's chief bioinformatics research engineer.
