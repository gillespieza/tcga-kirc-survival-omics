# 01_setup_and_data_retrieval.R
# Purpose: load required packages and read local TCGA KIRC cBioPortal data files.

# Package Setup ----------------------------------------------------------------

# CRAN packages used for data handling, plotting, and survival analysis
cran_packages <- c(
   "tidyverse", # Loads dplyr, tidyr, readr, ggplot2, stringr, and related tools
   "ggpubr",    # Provides plotting helpers used by survminer
   "car",       # Provides statistical helpers required by ggpubr
   "markdown",  # Supports markdown rendering used by plotting/report helpers
   "survival",  # Fits survival models
   "survminer", # Creates Kaplan-Meier and survival plots
   "glmnet"     # Fits penalised models, including LASSO Cox regression
)

# Install missing packages if not already installed
for (pkg in cran_packages) {
   if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
   }
}

# Load packages
library(tidyverse)
library(survival)
library(survminer)
library(glmnet)

# Project settings --------------------------------------------------------

# Kidney Renal Clear Cell Carcinoma, TCGA PanCancer Atlas
study_id <- "kirc_tcga_pan_can_atlas_2018"

# These files were downloaded from cBioPortal and saved locally so that the
# analysis can be rerun even if the cBioPortal API is unavailable.
clinical_patient_file <- "data/data_clinical_patient.txt"
clinical_sample_file <- "data/data_clinical_sample.txt"
mutation_file <- "data/data_mutations.txt"
cna_file <- "data/data_cna.txt"
rppa_file <- "data/data_rppa_zscores.txt"

message("Package setup complete.")
message("Working study ID: ", study_id)

# Read local cBioPortal data files ----------------------------------------

# Patient-level clinical data contain survival and patient-level variables.
clinical_patient <- read_tsv(
   clinical_patient_file,
   comment = "#",
   show_col_types = FALSE
)

# Sample-level clinical data contain sample IDs and tumour/sample variables.
clinical_sample <- read_tsv(
   clinical_sample_file,
   comment = "#",
   show_col_types = FALSE
)

# Mutation, copy-number, and RPPA files provide the ~omics data layers.
mutation_data <- read_tsv(
   mutation_file,
   comment = "#",
   show_col_types = FALSE
)

cna_data <- read_tsv(
   cna_file,
   comment = "#",
   show_col_types = FALSE
)

rppa_data <- read_tsv(
   rppa_file,
   comment = "#",
   show_col_types = FALSE
)

# Combine patient-level and sample-level clinical data.
# This creates one table containing survival variables and sample identifiers.
clinical_data <- clinical_sample %>%
   left_join(clinical_patient, by = "PATIENT_ID")

# Quick checks on imported data -------------------------------------------

cat("Clinical patient data:", nrow(clinical_patient), "rows x", ncol(clinical_patient), "columns\n")
cat("Clinical sample data:", nrow(clinical_sample), "rows x", ncol(clinical_sample), "columns\n")
cat("Combined clinical data:", nrow(clinical_data), "rows x", ncol(clinical_data), "columns\n")
cat("Mutation data:", nrow(mutation_data), "rows x", ncol(mutation_data), "columns\n")
cat("CNA data:", nrow(cna_data), "rows x", ncol(cna_data), "columns\n")
cat("RPPA data:", nrow(rppa_data), "rows x", ncol(rppa_data), "columns\n")

# Uncomment to inspect column names if needed
# cat(names(clinical_data), sep = "\n")
# cat(names(mutation_data), sep = "\n")
# cat(names(cna_data), sep = "\n")
# cat(names(rppa_data), sep = "\n")

# Prepare clinical survival table -----------------------------------------

clinical_survival <- clinical_data %>%
   dplyr::transmute(
      # Keep identifiers and rename them to tidy names
      patient_id = .data$PATIENT_ID,
      sample_id = .data$SAMPLE_ID,
      
      # Convert overall survival time from text to numeric months
      os_months = as.numeric(.data$OS_MONTHS),
      
      # Convert overall survival status into a Cox model event indicator:
      # 1 = death event, 0 = censored/alive
      os_event = dplyr::if_else(
         stringr::str_detect(.data$OS_STATUS, "DECEASED"),
         1,
         0
      ),
      
      # Convert age from text to numeric years
      age = as.numeric(.data$AGE),
      
      # Convert clinical categories into factors for modelling
      sex = as.factor(.data$SEX),
      stage = as.factor(.data$AJCC_PATHOLOGIC_TUMOR_STAGE),
      grade = as.factor(.data$GRADE)
   ) %>%
   # Keep only patients with usable overall survival time and status
   dplyr::filter(!is.na(.data$os_months), !is.na(.data$os_event))

# Summarise available survival data
clinical_survival_summary <- clinical_survival %>%
   summarise(
      patients = n_distinct(patient_id),
      samples = n_distinct(sample_id),
      events = sum(os_event),
      median_follow_up_months = median(os_months, na.rm = TRUE)
   )

print(clinical_survival_summary)

# Quick overall survival check --------------------------------------------

# Create a survival object from the cleaned clinical survival table.
# time = follow-up or survival time in months
# event = 1 if the patient died, 0 if the patient was censored/alive
surv_obj <- Surv(
   time = clinical_survival$os_months,
   event = clinical_survival$os_event
)

# Fit a Kaplan-Meier curve for all patients together
km_fit <- survfit(surv_obj ~ 1)

# Print a compact summary of the Kaplan-Meier fit
print(km_fit)

# Plot the overall survival curve
suppressMessages(
   suppressWarnings({
      km_plot <- ggsurvplot(
         km_fit,
         data = clinical_survival,
         risk.table = TRUE,
         legend = "none",
         xlab = "Time (months)",
         ylab = "Overall survival probability"
      )
      
      print(km_plot)
   })
)