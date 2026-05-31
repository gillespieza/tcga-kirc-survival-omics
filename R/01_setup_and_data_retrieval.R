# 01_setup_and_data_retrieval.R
# Purpose: install/load required packages and retrieve TCGA KIRC data from cBioPortal.

# Package Setup -------------------------------------------------------------------

# CRAN packages used for data handling, plotting, and survival analysis
cran_packages <- c(
  "BiocManager", # Installs packages from Bioconductor
  "tidyverse",   # Loads dplyr, tidyr, readr, ggplot2, stringr, and related tools
  "ggpubr",      # Provides plotting helpers used by survminer
  "car",         # Provides statistical helpers required by ggpubr
  "markdown",    # Supports markdown rendering used by plotting/report helpers
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
  # print(n = Inf)


# Retrieve clinical data ---------------------------------------------------

clinical_data <- clinicalData(cbio, studyId = study_id)

# Uncomment to find all the column names
# Print all clinical column names, one per line
# cat(names(clinical_data), sep = "\n")


# Prepare clinical survival table -----------------------------------------

clinical_survival <- clinical_data %>%
   select(
      patientId,
      sampleId,
      OS_MONTHS,
      OS_STATUS,
      AGE,
      SEX,
      AJCC_PATHOLOGIC_TUMOR_STAGE,
      GRADE
   ) %>%
   mutate(
      # Convert overall survival time from text to numeric months
      os_months = as.numeric(OS_MONTHS),
      
      # Convert overall survival status into a Cox model event indicator:
      # 1 = death event, 0 = censored/alive
      os_event = if_else(str_detect(OS_STATUS, "DECEASED"), 1, 0),
      
      # Convert age from text to numeric years
      age = as.numeric(AGE),
      
      # Convert clinical categories into factors for modelling
      sex = as.factor(SEX),
      stage = as.factor(AJCC_PATHOLOGIC_TUMOR_STAGE),
      grade = as.factor(GRADE)
   ) %>%
   select(
      patientId,
      sampleId,
      os_months,
      os_event,
      age,
      sex,
      stage,
      grade
   ) %>%
   # Keep only patients with usable overall survival time and status
   filter(!is.na(os_months), !is.na(os_event))

# Summarise available survival data
clinical_survival %>%
   summarise(
      patients = n_distinct(patientId),
      samples = n_distinct(sampleId),
      events = sum(os_event),
      median_follow_up_months = median(os_months, na.rm = TRUE)
   ) %>%
   print()


# Quick overall survival check --------------------------------------------

# Create a survival object from the raw clinical columns.
# time = follow-up or survival time in months
# event = TRUE if the patient died, FALSE if the patient was censored/alive
surv_obj <- Surv(
   time = as.numeric(clinical_data$OS_MONTHS),
   event = clinical_data$OS_STATUS == "1:DECEASED"
)

# Fit a Kaplan-Meier curve for all patients together
km_fit <- survfit(surv_obj ~ 1)

# Print a short summary of the Kaplan-Meier fit
summary(km_fit)

# Plot the overall survival curve
km_plot <- ggsurvplot(
   km_fit,
   data = clinical_data,
   risk.table = TRUE,
   xlab = "Time (months)",
   ylab = "Overall survival probability"
)

print(km_plot)
