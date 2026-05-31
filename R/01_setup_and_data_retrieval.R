# 01_setup_and_data_retrieval.R
# Purpose: install/load required packages and retrieve TCGA KIRC data from cBioPortal.

# Package Setup -------------------------------------------------------------------

# CRAN packages used for data handling, plotting, and survival analysis
cran_packages <- c(
  "BiocManager", # Installs packages from Bioconductor
  "tidyverse",   # Loads dplyr, tidyr, readr, ggplot2, stringr, and related tools
  "ggpubr",      # Provides plotting helpers used by survminer
  "car",         # Provides statistical helpers required by ggpubr
  "survival",    # Fits survival models
  "survminer",   # Creates Kaplan-Meier and survival plots
  "glmnet"       # Fits penalised models, including LASSO Cox regression
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
library(tidyverse)
library(survival)
library(survminer)
library(glmnet)
library(cBioPortalData)

# Project settings --------------------------------------------------------

# Kidney Renal Clear Cell Carcinoma, TCGA PanCancer Atlas
study_id <- "kirc_tcga_pan_can_atlas_2018"

# GENOMICS/MUTATION profile - can use it to check for alterations in known ccRCC
# driver genes such as VHL, PBRM1, BAP1, SETD2, KDM5C, MTOR, PTEN, TSC1
mutation_profile_id <- "kirc_tcga_pan_can_atlas_2018_mutations"

# GENOMICS/COPY NUMBER ALTERATION: profile from GISTIC
cna_profile_id <- "kirc_tcga_pan_can_atlas_2018_gistic"

# PROTEOMICS/RPPA protein-level RPPA profile as Z-scores - Z-scores are
# standardised, make comparing proteins easier in survival modelling than raw abundance
rppa_profile_id <- "kirc_tcga_pan_can_atlas_2018_rppa_Zscores"

# GENOMICS/RNA expression - may want this later. Chose median Z-score version
# rather than raw mRNA because it is already centred relative to the cohort
rna_profile_id <- "kirc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores"

# Sample lists define which samples are available for each data type
all_sample_list_id <- "kirc_tcga_pan_can_atlas_2018_all"             # All samples in the study
rppa_sample_list_id <- "kirc_tcga_pan_can_atlas_2018_rppa"           # Samples with RPPA protein data
mutation_sample_list_id <- "kirc_tcga_pan_can_atlas_2018_sequenced"  # Samples with mutation data
cna_sample_list_id <- "kirc_tcga_pan_can_atlas_2018_cna"             # Samples with copy-number data
complete_sample_list_id <- "kirc_tcga_pan_can_atlas_2018_3way_complete" # Samples with mutation, CNA, and mRNA data

message("Package setup complete.")
message("Working study ID: ", study_id)

# Connect to cBioPortal ---------------------------------------------------

cbio <- cBioPortal()

# Check that the selected study is available
studies <- getStudies(cbio)

kirc_study <- studies %>%
  filter(studyId == study_id)

# Uncomment to confirm the selected study metadata
# print(kirc_study)

# List available molecular profiles for this study
molecular_profiles <- molecularProfiles(cbio, studyId = study_id)

# Uncomment to inspect available molecular profiles
# molecular_profiles %>%
#   select(molecularProfileId, name, molecularAlterationType, datatype) %>%
#   print(n = Inf)

# Uncomment to print the full molecular profile IDs without truncation
# molecular_profiles %>%
#   select(molecularProfileId, name, molecularAlterationType, datatype) %>%
#   as.data.frame()

# Inspect sample lists ----------------------------------------------------

sample_lists <- sampleLists(cbio, studyId = study_id)

# Uncomment to inspect the sample lists
# sample_lists %>%
#   select(sampleListId, name, category) %>%
#   print(n = Inf)