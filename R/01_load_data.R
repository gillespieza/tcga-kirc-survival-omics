#' @file        01_load_data.R
#' @title       Load Raw cBioPortal Data Files
#' @description Reads all five local cBioPortal data files into memory as
#'              tibbles. File paths are expected to have been declared in
#'              00_setup.R. Comment lines (prefixed with #) are skipped
#'              automatically by read_tsv.
#'
#' @details
#'   Objects created in the global environment:
#'   \describe{
#'     \item{clinical_patient}{Patient-level clinical data. Contains survival
#'           outcomes and patient-level covariates (e.g. age, sex, stage).}
#'     \item{clinical_sample}{Sample-level clinical data. Contains sample IDs
#'           and tumour/sample-level variables (e.g. tumour site, purity).}
#'     \item{mutation_data}{Somatic mutation calls, one row per
#'           variant/sample combination.}
#'     \item{cna_data}{Gene-level copy-number alterations, typically in a
#'           wide gene x sample matrix format.}
#'     \item{rppa_data}{RPPA protein expression z-scores, typically in a
#'           wide protein x sample matrix format.}
#'   }
#'
#' @pre  00_setup.R must be sourced first so that the following path variables
#'       are defined: \code{clinical_patient_file}, \code{clinical_sample_file},
#'       \code{mutation_file}, \code{cna_file}, \code{rppa_file}.
#'
#' @post Five tibbles are available in the global environment for use by
#'       downstream preprocessing and analysis scripts.
#'
#' @usage source("01_load_data.R")

# Helper ----------------------------------------------------------------------

#' Read a cBioPortal TSV, skipping header comment lines
#'
#' @param path Character. Path to a cBioPortal-format TSV file.
#' @return A tibble.
read_cbio <- function(path) {
   read_tsv(path, comment = "#", show_col_types = FALSE)
}


# Load data -------------------------------------------------------------------

# Patient-level clinical data: survival outcomes and patient-level covariates.
clinical_patient <- read_cbio(clinical_patient_file)
stopifnot(
   "clinical_patient must be non-empty"          = nrow(clinical_patient) > 0,
   "clinical_patient must contain PATIENT_ID"    = "PATIENT_ID"  %in% names(clinical_patient),
   "clinical_patient must contain OS_STATUS"     = "OS_STATUS"   %in% names(clinical_patient),
   "clinical_patient must contain OS_MONTHS"     = "OS_MONTHS"   %in% names(clinical_patient)
)
message("Loaded clinical_patient:  ", nrow(clinical_patient), " rows x ",
        ncol(clinical_patient), " cols")

# Sample-level clinical data: sample IDs and tumour/sample-level variables.
clinical_sample <- read_cbio(clinical_sample_file)
stopifnot(
   "clinical_sample must be non-empty"           = nrow(clinical_sample) > 0,
   "clinical_sample must contain PATIENT_ID"     = "PATIENT_ID" %in% names(clinical_sample),
   "clinical_sample must contain SAMPLE_ID"      = "SAMPLE_ID"  %in% names(clinical_sample)
)
message("Loaded clinical_sample:   ", nrow(clinical_sample), " rows x ",
        ncol(clinical_sample), " cols")

# Somatic mutation calls: one row per variant/sample combination.
mutation_data <- read_cbio(mutation_file)
stopifnot(
   "mutation_data must be non-empty"             = nrow(mutation_data) > 0,
   "mutation_data must contain Hugo_Symbol"      = "Hugo_Symbol"         %in% names(mutation_data),
   "mutation_data must contain Tumor_Sample_Barcode" = "Tumor_Sample_Barcode" %in% names(mutation_data)
)
message("Loaded mutation_data:     ", nrow(mutation_data), " rows x ",
        ncol(mutation_data), " cols")

# Gene-level copy-number alterations: wide gene x sample matrix.
cna_data <- read_cbio(cna_file)
stopifnot(
   "cna_data must be non-empty"                  = nrow(cna_data) > 0,
   "cna_data must contain Hugo_Symbol"           = "Hugo_Symbol" %in% names(cna_data)
)
message("Loaded cna_data:          ", nrow(cna_data), " rows x ",
        ncol(cna_data), " cols")


# RPPA protein expression z-scores: wide protein x sample matrix.
rppa_data <- read_cbio(rppa_file)
stopifnot(
   "rppa_data must be non-empty"                 = nrow(rppa_data) > 0
)
message("Loaded rppa_data:         ", nrow(rppa_data), " rows x ",
        ncol(rppa_data), " cols")


message("All data files loaded successfully.")
