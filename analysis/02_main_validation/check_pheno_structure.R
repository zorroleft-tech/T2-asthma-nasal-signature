# Script to check pheno file structure

# Set project root
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Read a pheno file
pheno_path <- file.path(base_dir, "data", "processed_diet", "GSE103166_pheno.rds")
pheno_data <- readRDS(pheno_path)

# Print structure
cat("Structure of pheno data:\n")
str(pheno_data)

# Print column names
cat("\nColumn names:\n")
print(colnames(pheno_data))

# Print first few rows
cat("\nFirst few rows:\n")
print(head(pheno_data))
