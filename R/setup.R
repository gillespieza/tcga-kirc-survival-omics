# Project package setup and configuration -------------------------------------
#
# Bootstraps the R environment for this TCGA KIRC survival analysis mini ML
# project. Checks for and installs any missing CRAN packages, loads the core
# libraries required by downstream scripts, and declares project-level constants
# for the study identifier and raw data file paths.
#
# Packages installed/loaded:
#   here      - Reliable project-relative file paths
#   tidyverse - Data manipulation/visualisation (dplyr, ggplot2, stringr, etc.)
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
#   rnaseq_file           - Path to RNA-seq expression data
#
# Note: data files must be pre-downloaded from cBioPortal and placed in data/
# so the analysis can be reproduced without a live API connection. A global
# random seed is set here to ensure reproducibility of stochastic steps
# (e.g. glmnet cross-validation) across all downstream scripts.
#
# src: https://www.cbioportal.org/study/summary?id=kirc_tcga_pan_can_atlas_2018
#
# Usage: source("run_analysis.R") from the project root
# This script is intended to be sourced by run_analysis.R, not run directly.


# Project Root Validation ------------------------------------------------------
# Confirm that the script is being run from the project root. This helps avoid
# confusing path errors later when here::here() is used to build file paths.

if (!dir.exists("R") || !dir.exists("data")) {
   stop(
      "Project directories 'R' and/or 'data' were not found.\n",
      "Run run_analysis.R from the project root.",
      call. = FALSE
   )
}


# R Version Guard --------------------------------------------------------------
# Stop early if the installed R version is too old for this project.

if (getRversion() < "4.1.0") {
   stop(
      "R >= 4.1.0 is required.",
      call. = FALSE
   )
}


# Package Configuration --------------------------------------------------------
# Define the CRAN packages used across the project. Missing packages are
# installed automatically so the pipeline is easier to run while learning.

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
   "progress",  # Progress bars for long-running steps 
   "msigdbr"    # MSigDB gene sets for pathway analysis
)


# Dependency Management Helpers ------------------------------------------------

# Install any package that is not already available in the current R library.
install_missing_packages <- function(packages) {
   # Identify packages not already present in the library
   is_installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
   missing_packages <- packages[!is_installed]
   
   if (length(missing_packages) > 0L) {
      install.packages(
         missing_packages,
         repos = "https://cloud.r-project.org"
      )
   }
}

# Load each required package so downstream scripts can use them.
load_required_packages <- function(packages) {
   for (pkg in packages) {
      library(
         pkg,
         character.only = TRUE
      )
   }
}

# Execute installation and loading passes
install_missing_packages(cran_packages)
load_required_packages(cran_packages)


# Package Version Assertions ---------------------------------------------------
# Ensure package versions match the validated baseline analysis environment
# to protect the pipeline against breaking API changes.

current_msigdbr_version <- utils::packageVersion("msigdbr")

# Stop execution if the version is older than what your analysis expects
if (current_msigdbr_version < "7.5.1") {
   stop(
      "msigdbr >= 7.5.1 is required for gene set compatibility.\n",
      "Currently installed: ", current_msigdbr_version,
      call. = FALSE
   )
}

# message("msigdbr version verified: ", as.character(current_msigdbr_version))
message("Packages loaded successfully.")


# Deterministic Reproducibility ------------------------------------------------
# Set a global seed to ensure stochastic steps (e.g. glmnet cross-validation)
# produce identical results across runs.

set.seed(42)


# Global Theme and Typography Settings -----------------------------------------

ggplot2::theme_set(
   ggplot2::theme_classic(
      base_size   = 13,
      base_family = "sans"
   )
)


# Study Cohort Configurations --------------------------------------------------
# Kidney Renal Clear Cell Carcinoma, TCGA PanCancer Atlas

study_id <- "kirc_tcga_pan_can_atlas_2018"


# Multi-Omics Data Layout Paths ------------------------------------------------
# Build project-relative paths to the local cBioPortal data files.

clinical_patient_file <- here::here("data", "data_clinical_patient.txt")
clinical_sample_file  <- here::here("data", "data_clinical_sample.txt")
mutation_file         <- here::here("data", "data_mutations.txt")
cna_file              <- here::here("data", "data_cna.txt")
rppa_file             <- here::here("data", "data_rppa_zscores.txt")
rnaseq_file           <- here::here("data", "data_mrna_seq_v2_rsem.txt")


# Raw Workspace Verification ---------------------------------------------------
# Check that all required raw data files exist before any downstream script
# tries to read them.

data_files <- c(
   clinical_patient = clinical_patient_file,
   clinical_sample  = clinical_sample_file,
   mutation         = mutation_file,
   cna              = cna_file,
   rppa             = rppa_file,
   rnaseq           = rnaseq_file
)

missing_files <- data_files[!file.exists(data_files)]

if (length(missing_files) > 0L) {
   stop(
      "The following data files were not found:\n",
      paste0(
         "  - ", names(missing_files), ": ", missing_files,
         collapse = "\n"
      ),
      "\nDownload them from cBioPortal and place them in the data/ directory.",
      call. = FALSE
   )
}

message("All data files found.")
message("Package setup complete.")
#message("Working study ID: ", study_id)