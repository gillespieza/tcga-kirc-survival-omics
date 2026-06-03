# Project package setup and configuration -------------------------------------
#
# Bootstraps the R environment for this TCGA KIRC survival analysis mini ML 
# project. Checks for and installs any missing CRAN packages, loads the core 
# libraries required by downstream scripts, and declares project-level constants 
# for the study identifier and raw data file paths.
#
# Packages installed/loaded:
#   here,     - Reliable project-relative file paths
#   tidyverse - Data manipulation and visualisation (dplyr, ggplot2, stringr, etc.)
#   ggpubr    - Publication-ready plot helpers (dependency of survminer)
#   car       - Statistical utilities (dependency of ggpubr)
#   markdown  - Markdown rendering (dependency of plotting/report helpers)
#   survival  - Core survival-analysis models (Surv, coxph, etc.)
#   survminer - Kaplan-Meier plots and survival-curve visualisation
#   glmnet    - Penalised regression, including LASSO Cox models
#   broom     - Tidy model summaries for Cox model results
#   msigdbr   - MSigDB gene sets for pathway analysis
#
# Global constants defined:
#   study_id              - cBioPortal study identifier
#   clinical_patient_file - Path to patient-level clinical data
#   clinical_sample_file  - Path to sample-level clinical data
#   mutation_file         - Path to somatic mutation data
#   cna_file              - Path to copy-number alteration data
#   rppa_file             - Path to RPPA protein expression z-scores
#
# Note: data files must be pre-downloaded from cBioPortal and placed in data/
# so the analysis can be reproduced without a live API connection. A global
# random seed is set here to ensure reproducibility of stochastic steps
# (e.g. glmnet cross-validation) across all downstream scripts.
#
# Source: https://www.cbioportal.org/study/summary?id=kirc_tcga_pan_can_atlas_2018
# Usage:  source("00_setup.R")  # call before any other project script


# R version guard -------------------------------------------------------------

if (getRversion() < "4.1.0") stop("R >= 4.1.0 is required.")


# Package setup ---------------------------------------------------------------
# CRAN packages used for data handling, plotting, and survival analysis.
# All packages are explicitly loaded below, including those that are also
# required as dependencies (ggpubr, car, markdown), so that missing or
# broken installations surface with a clear error at startup.

cran_packages <- c(
   "here",      # Reliable project-relative file paths
   "tidyverse", # Data manipulation and visualisation
   "ggpubr",    # Publication-ready plot helpers (dependency of survminer)
   "car",       # Statistical utilities (dependency of ggpubr)
   "markdown",  # Markdown rendering (dependency of plotting/report helpers)
   "survival",  # Core survival-analysis models
   "survminer", # Kaplan-Meier plots and survival-curve visualisation
   "glmnet",    # Penalised regression, including LASSO Cox models
   "broom",     # Tidy model summaries for Cox model results
   # "ellmer",    # LLM
   "msigdbr"    # MSigDB gene sets for pathway analysis
)

for (pkg in cran_packages) {
   if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
   }
}

library(here)
library(tidyverse)
library(ggpubr)
library(car)
library(markdown)
library(survival)
library(survminer)
library(glmnet)
library(broom)
library(msigdbr)
# library(ellmer)

message("Packages loaded successfully.")


# Reproducibility -------------------------------------------------------------
# Set a global seed to ensure stochastic steps (e.g. glmnet cross-validation)
# produce identical results across runs.

set.seed(42)


# Project settings ------------------------------------------------------------
# Kidney Renal Clear Cell Carcinoma, TCGA PanCancer Atlas

study_id <- "kirc_tcga_pan_can_atlas_2018"

# These files were downloaded from cBioPortal and saved locally so that the
# analysis can be rerun even if the cBioPortal API is unavailable.
# here::here() ensures paths resolve correctly regardless of working directory.

clinical_patient_file <- here::here("data", "data_clinical_patient.txt")
clinical_sample_file  <- here::here("data", "data_clinical_sample.txt")
mutation_file         <- here::here("data", "data_mutations.txt")
cna_file              <- here::here("data", "data_cna.txt")
rppa_file             <- here::here("data", "data_rppa_zscores.txt")
rnaseq_file           <- here::here("data", "data_mrna_seq_v2_rsem.txt")


# File existence check --------------------------------------------------------
# Validate all data files are present before any downstream script attempts to
# read them, so failures are caught early with an informative error message.

data_files <- c(
   clinical_patient = clinical_patient_file,
   clinical_sample  = clinical_sample_file,
   mutation         = mutation_file,
   cna              = cna_file,
   rppa             = rppa_file,
   rnaseq           = rnaseq_file
)

missing_files <- data_files[!file.exists(data_files)]

if (length(missing_files) > 0) {
   stop(
      "The following data files were not found:\n",
      paste0("  [", names(missing_files), "] ", missing_files, collapse = "\n"),
      "\nDownload them from cBioPortal and place them in the data/ directory."
   )
}

message("All data files found.")
message("Package setup complete.")
message("Working study ID: ", study_id)