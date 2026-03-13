#!/usr/bin/env Rscript

library(dplyr)

base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(base_dir, "data", "processed", "phase3_outputs")

cat("Checking GSE65204...\n")
kv_65204 <- readRDS(file.path(processed_dir, "GSE65204_pheno_kv_long.rds"))
cat("\nGSE65204 keys:\n")
print(table(kv_65204$key))
cat("\nGSE65204 first 30 rows:\n")
print(head(kv_65204, 30))

cat("\n\nChecking GSE201955...\n")
kv_201955 <- readRDS(file.path(processed_dir, "GSE201955_pheno_kv_long.rds"))
cat("\nGSE201955 keys:\n")
print(table(kv_201955$key))
cat("\nGSE201955 first 30 rows:\n")
print(head(kv_201955, 30))

cat("\n\nChecking GSE45111...\n")
kv_45111 <- readRDS(file.path(processed_dir, "GSE45111_pheno_kv_long.rds"))
cat("\nGSE45111 keys:\n")
print(table(kv_45111$key))

cat("\n\nChecking GSE45111 pheno_samples...\n")
pheno_45111 <- readRDS(file.path(processed_dir, "GSE45111_pheno_samples.rds"))
print(table(pheno_45111$std_group_label, useNA = "always"))
cat("\nstd_group_label_source:\n")
print(head(pheno_45111$std_group_label_source, 10))
