# Integrate clinical, RPPA, RNA-seq, and CNA data ----------------------------
#
# Combines the prepared clinical survival, RPPA proteomics, RNA-seq, and
# copy-number alteration tables into a single analysis-ready tibble.
# Standardises all sample keys to ensure perfect multi-omics cross-matching.
#
# Requires: prepare_clinical.R, prepare_rppa.R, prepare_rnaseq.R, and
#           prepare_cna.R to have been sourced.
#
# Produces:
#   clinical_rppa_rna_cna - Fully integrated analysis tibble used
#                           downstream by modelling scripts.
#
# Usage: this script is intended to be sourced by run_analysis.R.


# Validate Inputs -------------------------------------------------------------

required_tables <- c(
  "clinical_survival",
  "rppa_proteomics",
  "rnaseq_expression",
  "cna_features"
)

# Enforce defensive checks to ensure all upstream tables are present in memory
check_required_objects(required_tables)

# Guarantee every table uses identical key formatting
clinical_survival <- clinical_survival |>
  dplyr::mutate(sample_id = standardise_sample_id(.data$sample_id))

rppa_proteomics <- rppa_proteomics |>
  dplyr::mutate(sample_id = standardise_sample_id(.data$sample_id))

rnaseq_expression <- rnaseq_expression |>
  dplyr::mutate(sample_id = standardise_sample_id(.data$sample_id))

cna_features <- cna_features |>
  dplyr::mutate(sample_id = standardise_sample_id(.data$sample_id))

for (tbl in required_tables) {
  check_has_sample_id(tbl)
}


# 1. Integrate Clinical and RPPA ----------------------------------------------

clinical_rppa <- clinical_survival |>
  dplyr::inner_join(
    rppa_proteomics,
    by = "sample_id"
  )


# 2. Integrate RNA-seq --------------------------------------------------------

clinical_rppa_rna <- clinical_rppa |>
  dplyr::inner_join(
    rnaseq_expression,
    by = "sample_id"
  )

n_dropped_rna <- length(setdiff(
  clinical_rppa$sample_id,
  clinical_rppa_rna$sample_id
))

if (n_dropped_rna > 0L) {
  message(
    n_dropped_rna,
    " sample(s) with clinical + RPPA data had no matching RNA-seq ",
    "and were dropped by the RNA-seq inner join."
  )
}


# 3. Integrate Copy Number Alteration (CNA) Features -------------------------
# Missing values are only treated as unaltered (0L) if the sample is confirmed
# to have successfully undergone copy-number profiling.

cna_sample_universe <- standardise_sample_id(setdiff(
  names(cna_data),
  c("Hugo_Symbol", "Entrez_Gene_Id")
))

clinical_rppa_rna_cna <- clinical_rppa_rna |>
  dplyr::left_join(
    cna_features,
    by = "sample_id"
  )

# Pre-calculate logical alignment vector to avoid across() scoping bugs
is_cna_profiled <- clinical_rppa_rna_cna$sample_id %in%
  cna_sample_universe

clinical_rppa_rna_cna <- clinical_rppa_rna_cna |>
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("cna_", ignore.case = FALSE),
      ~ dplyr::if_else(is_cna_profiled & is.na(.x), 0L, .x)
    )
  )

abort_if_false(
  nrow(clinical_rppa_rna_cna) == nrow(clinical_rppa_rna),
  paste(
    "Row count changed after CNA join. Check cna_features for duplicate",
    "sample_id values."
  )
)


# Validate Final Output -------------------------------------------------------

abort_if_false(
  nrow(clinical_rppa_rna_cna) > 0L,
  "clinical_rppa_rna_cna must be non-empty."
)

check_has_columns(
  "clinical_rppa_rna_cna",
  c("sample_id", "os_months", "os_event")
)


# Missingness Diagnostic Audit ------------------------------------------------
# Check if any numeric multi-omics features are carrying unintended NAs

na_audit <- clinical_rppa_rna_cna |>
  dplyr::summarise(
    total_rows      = dplyr::n(),
    clinical_age_na = sum(is.na(.data$age)),
    cna_na          = sum(is.na(dplyr::pick(dplyr::starts_with("cna_"))))
  )

message("Multi-omics integration quality audit:")
print(na_audit)

message(
  "Integration complete. ",
  nrow(clinical_rppa_rna_cna),
  " samples fully combined into multi-omics workspace."
)
