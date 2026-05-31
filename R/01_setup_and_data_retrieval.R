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
# print(km_fit)
km_median <- as.numeric(km_fit$table[["median"]])

cat(
   "Overall survival:",
   km_fit$n,
   "patients,",
   km_fit$events,
   "events, median survival",
   if_else(is.na(km_median), "not reached", as.character(round(km_median, 1))),
   "months\n"
)

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


# Inspect RPPA protein feature names --------------------------------------

# RPPA files usually have one row per antibody/protein feature and one column
# per tumour sample. The first column identifies the protein feature.
# rppa_data %>%
#   dplyr::select(1:6) %>%
#   head(10) %>%
#   print()


# Prepare RPPA proteomics table -------------------------------------------

# Convert RPPA from wide feature-by-sample into sample-by-feature
# This gives one row/tumour sample and one column/protein feature
rppa_proteomics <- rppa_data %>%
   pivot_longer(
      cols = -Composite.Element.REF,
      names_to = "sample_id",
      values_to = "rppa_zscore"
   ) %>%
   separate(
      Composite.Element.REF,
      into = c("gene_symbol", "protein_feature"),
      sep = "\\|",
      remove = FALSE,
      extra = "merge",
      fill = "right"
   ) %>%
   mutate(
      # Make protein feature names safe to use as column names later
      protein_feature_clean = Composite.Element.REF %>%
         str_replace_all("[^A-Za-z0-9]+", "_") %>%
         str_replace_all("_$", "")
   ) %>%
   select(sample_id, protein_feature_clean, rppa_zscore) %>%
   pivot_wider(
      names_from = protein_feature_clean,
      values_from = rppa_zscore
   )

# Inspect model-ready RPPA table
cat("Model-ready RPPA data:", nrow(rppa_proteomics), "rows x", ncol(rppa_proteomics), "columns\n")

# rppa_proteomics %>%
#   select(1:6) %>%
#   head() %>%
#   print()

# Integrate clinical and RPPA proteomics data -----------------------------

# inner_join will keep only rows with both data
clinical_rppa <- clinical_survival %>%
   inner_join(rppa_proteomics, by = "sample_id")

cat("Clinical + RPPA data:", nrow(clinical_rppa), "rows x", ncol(clinical_rppa), "columns\n")

# clinical_rppa %>%
#   select(1:10) %>%
#   head() %>%
#   print()

# Inspect mutation data ----------------------------------------------------

# cat("Mutation data column names:\n")
# cat(names(mutation_data), sep = "\n")

# Prepare mutation features ------------------------------------------------

# Commonly altered ccRCC genes selected from known kidney cancer biology.
# from doi:10.1038/nature12222 and 10.1158/1078-0432.CCR-15-2631
driver_genes <- c(
   "VHL",   # Core ccRCC tumour suppressor; VHL loss drives HIF/hypoxia and angiogenesis biology
   "PBRM1", # Chromatin-remodelling tumour suppressor frequently mutated in ccRCC
   "BAP1",  # Tumour suppressor associated with more aggressive ccRCC and poorer prognosis
   "SETD2", # Chromatin/histone methyltransferase gene altered in ccRCC and linked to genomic regulation
   "KDM5C", # Chromatin-regulation gene recurrently altered in ccRCC
   "MTOR",  # Kinase pathway gene; links directly to PI3K/AKT/mTOR signalling and RCC targeted therapy
   "PTEN",  # Negative regulator of PI3K/AKT signalling; recurrently altered in TCGA ccRCC
   "TSC1",  # mTOR pathway regulator; TSC1/TSC2/MTOR mutations are linked to rapalog response in metastatic RCC
   "TSC2"   # mTOR pathway regulator that functions with TSC1 to suppress mTORC1 signalling
)

# Convert the long mutation table into sample-level binary mutation features.

# Keep only sample ID and gene symbol, then restrict to selected genes.
mutation_long <- mutation_data %>%
   dplyr::transmute(
      sample_id = .data$Tumor_Sample_Barcode,
      gene_symbol = .data$Hugo_Symbol
   ) %>%
   dplyr::filter(.data$gene_symbol %in% driver_genes) %>%
   dplyr::distinct(.data$sample_id, .data$gene_symbol)

# Create a sample-by-gene mutation matrix.
# table() counts whether each selected gene appears in each sample.
mutation_matrix <- table(
   mutation_long$sample_id,
   mutation_long$gene_symbol
)
   

# Convert the matrix to a data frame and make values binary:
# 1 = mutated, 0 = not mutated.
mutation_features <- as.data.frame.matrix(mutation_matrix) %>%
   tibble::rownames_to_column("sample_id") %>%
   dplyr::mutate(
      dplyr::across(
         -sample_id,
         ~ as.integer(.x > 0)
      )
   ) %>%
   dplyr::rename_with(
      ~ paste0("mut_", .x),
      -sample_id
   )

# Add any selected driver genes that were not observed in the mutation data.
# This keeps the feature table consistent even if a gene has zero mutations.
missing_mutation_cols <- paste0("mut_", driver_genes)[
   !paste0("mut_", driver_genes) %in% names(mutation_features)
]

for (col in missing_mutation_cols) {
   mutation_features[[col]] <- 0
}

# Keep sample_id first, then mutation columns in the same order as driver_genes.
mutation_features <- mutation_features %>%
   dplyr::select(
      sample_id,
      dplyr::all_of(paste0("mut_", driver_genes))
   )

cat("Mutation feature data:", nrow(mutation_features), "rows x", ncol(mutation_features), "columns\n")

mutation_features %>%
   dplyr::summarise(dplyr::across(dplyr::starts_with("mut_"), sum)) %>%
   print()