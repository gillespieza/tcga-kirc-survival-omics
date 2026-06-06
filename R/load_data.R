# Load raw cBioPortal data files ----------------------------------------------
#
# Reads all local cBioPortal data files into memory as tibbles. File paths are
# expected to have been declared in setup.R. Comment lines (prefixed with #)
# are skipped automatically by read_tsv.
#
# Requires: setup.R to have been sourced by run_analysis.R so that the
# following path variables are defined: clinical_patient_file,
# clinical_sample_file, mutation_file, cna_file, rppa_file, rnaseq_file.
#
# Produces:
# clinical_patient - Patient-level clinical data. Contains survival outcomes
#                    and patient-level covariates (e.g. age, sex, stage).
# clinical_sample  - Sample-level clinical data. Contains sample IDs and
#                    tumour/sample-level variables (e.g. tumour site, purity).
# mutation_data    - Somatic mutation calls, one row per variant/sample.
# cna_data         - Gene-level copy-number alterations, wide gene x sample
#                    matrix.
# rppa_data        - RPPA protein expression z-scores, wide protein x sample
#                    matrix.
# rnaseq_data      - RNA-seq expression data, wide gene x sample matrix.
#
# Note: this script is intended to be sourced by run_analysis.R as part of
# the full pipeline, not run directly.


# Required Objects Check -------------------------------------------------------
# Confirm that setup.R has already created the file path objects needed here.

required_path_objects <- c(
  "clinical_patient_file",
  "clinical_sample_file",
  "mutation_file",
  "cna_file",
  "rppa_file",
  "rnaseq_file"
)

# Create a vector of any missing path objects.
missing_path_objects <- required_path_objects[
  !vapply(
    required_path_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

# If any are missing, stop with an informative error message.
if (length(missing_path_objects) > 0L) {
  stop(
    "The following required file path objects are missing:\n",
    paste0("  - ", missing_path_objects, collapse = "\n"),
    "\nSource setup.R via run_analysis.R before running load_data.R.",
    call. = FALSE
  )
}


# File IO Helper ---------------------------------------------------------------
# Read a cBioPortal-format TSV file, skipping comment lines that begin with #.

read_cbio <- function(path) {
  readr::read_tsv(
    file           = path,
    comment        = "#",
    show_col_types = FALSE
  )
}


# Validation Helper ------------------------------------------------------------
# Load one dataset, check that it is non-empty, check for required columns,
# print a short summary message, and return the loaded tibble.

load_and_validate_data <- function(object_name,
                                   file_path,
                                   required_columns,
                                   description) {
  data <- read_cbio(file_path)

  if (nrow(data) == 0L) {
    stop(
      object_name, " is empty after loading from:\n",
      file_path,
      call. = FALSE
    )
  }

  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      object_name, " is missing required columns:\n",
      paste0("  - ", missing_columns, collapse = "\n"),
      "\nFile: ", file_path,
      call. = FALSE
    )
  }

  message(
    "Loaded ", object_name, " (", description, "): ",
    nrow(data), " rows x ",
    ncol(data), " cols"
  )

  data
}


# Dataset Definitions ----------------------------------------------------------
# Define each dataset once so the loading loop below is shorter and easier to
# maintain.

dataset_specs <- tibble::tibble(
  object_name = c(
    "clinical_patient",
    "clinical_sample",
    "mutation_data",
    "cna_data",
    "rppa_data",
    "rnaseq_data"
  ),
  path_object = c(
    "clinical_patient_file",
    "clinical_sample_file",
    "mutation_file",
    "cna_file",
    "rppa_file",
    "rnaseq_file"
  ),
  description = c(
    "patient-level clinical data",
    "sample-level clinical data",
    "somatic mutation calls",
    "copy-number alteration matrix",
    "RPPA protein expression z-scores",
    "RNA-seq expression matrix"
  ),
  required_columns = list(
    c("PATIENT_ID", "OS_STATUS", "OS_MONTHS"),
    c("PATIENT_ID", "SAMPLE_ID"),
    c("Hugo_Symbol", "Tumor_Sample_Barcode"),
    c("Hugo_Symbol"),
    c("Composite.Element.REF"),
    c("Hugo_Symbol")
  )
)


# Load Data Execution ----------------------------------------------------------
# Iterate over the dataset definitions, load each file, validate it, and save
# the resulting tibble into the current environment.

for (i in seq_len(nrow(dataset_specs))) {
  object_name <- dataset_specs$object_name[[i]]
  path_object <- dataset_specs$path_object[[i]]
  description <- dataset_specs$description[[i]]
  required_columns <- dataset_specs$required_columns[[i]]

  file_path <- get(path_object, inherits = TRUE)

  loaded_data <- load_and_validate_data(
    object_name      = object_name,
    file_path        = file_path,
    required_columns = required_columns,
    description      = description
  )

  assign(
    x     = object_name,
    value = loaded_data,
    envir = .GlobalEnv
  )
}

message("All data files loaded successfully.")
