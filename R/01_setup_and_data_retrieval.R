# 01_setup_and_data_retrieval.R
# Purpose: install/load required packages and retrieve TCGA KIRC data from cBioPortal.

# ---- Package setup ----

# CRAN packages used for data handling, plotting, and survival analysis
cran_packages <- c(
  "BiocManager",
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "survival",
  "survminer",
  "glmnet"
)

# Install missing packages if not already installed
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor package for cBioPortal access
if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
  BiocManager::install("cBioPortalData")
}

# Load packages
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(survival)
library(survminer)
library(glmnet)
library(cBioPortalData)

# ---- Project settings ----

study_id <- "kirc_tcga_pan_can_atlas_2018"

message("Package setup complete.")
message("Working study ID: ", study_id)