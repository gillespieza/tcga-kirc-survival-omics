# Integrate clinical, RPPA, RNA-seq, mutations, and CNA data -----------------
#
# Combines the prepared clinical survival, RPPA proteomics, RNA-seq, binary
# mutation, and copy-number alteration tables into a single analysis-ready tibble.
# Standardises all sample keys to ensure perfect multi-omics cross-matching.
#
# Requires: prepare_clinical.R, prepare_rppa.R, prepare_rnaseq.R,
#           prepare_mutations.R, and prepare_cna.R to have been sourced.
#
# Produces:
#   clinical_rppa_rna_mutation_cna - Fully integrated analysis tibble used
#                                    downstream by modelling scripts.
#
# Usage: this script is intended to be sourced by run_analysis.R.


# Validate inputs -------------------------------------------------------------

required_tables <- c(
   "clinical_survival",
   "rppa_proteomics",
   "rnaseq_expression",
   "mutation_features",
   "cna_features"
)

check_required_objects(required_tables)

for (tbl in required_tables) {
   check_has_sample_id(tbl)
}


# 1. Integrate clinical and RPPA ----------------------------------------------

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

if (n_dropped_rna > 0) {
   message(
      n_dropped_rna,
      " sample(s) with clinical + RPPA data had no matching RNA-seq ",
      "and were dropped by the RNA-seq inner join."
   )
}


# 3. Integrate mutation features ----------------------------------------------
# Missing values are only treated as wild-type (0) if the sample is confirmed 
# to have undergone genomic sequencing.

sequenced_sample_universe <- standardise_sample_id(unique(mutation_data$Tumor_Sample_Barcode))

clinical_rppa_rna_mutation <- clinical_rppa_rna |>
   dplyr::left_join(
      mutation_features,
      by = "sample_id"
   )

# Pre-calculate logical alignment vector to avoid across() scoping bugs
is_sequenced <- clinical_rppa_rna_mutation$sample_id %in% sequenced_sample_universe

clinical_rppa_rna_mutation <- clinical_rppa_rna_mutation |>
   dplyr::mutate(
      dplyr::across(
         dplyr::starts_with("mut_", ignore.case = FALSE),
         ~ dplyr::if_else(is_sequenced & is.na(.x), 0L, .x)
      )
   )


# 4. Integrate Copy Number Alteration (CNA) features -------------------------
# Missing values are only treated as unaltered (0) if the sample is confirmed
# to have successfully undergone copy-number profiling.

cna_sample_universe <- standardise_sample_id(setdiff(
   names(cna_data),
   c("Hugo_Symbol", "Entrez_Gene_Id")
))

clinical_rppa_rna_mutation_cna <- clinical_rppa_rna_mutation |>
   dplyr::left_join(
      cna_features,
      by = "sample_id"
   )

# Pre-calculate logical alignment vector to avoid across() scoping bugs
is_cna_profiled <- clinical_rppa_rna_mutation_cna$sample_id %in% cna_sample_universe

clinical_rppa_rna_mutation_cna <- clinical_rppa_rna_mutation_cna |>
   dplyr::mutate(
      dplyr::across(
         dplyr::starts_with("cna_", ignore.case = FALSE),
         ~ dplyr::if_else(is_cna_profiled & is.na(.x), 0L, .x)
      )
   )

abort_if_false(
   nrow(clinical_rppa_rna_mutation_cna) == nrow(clinical_rppa_rna),
   "Row count changed after CNA join. Check cna_features for duplicate sample_id values."
)

# Rename to the standard master name expected by downstream files
clinical_rppa_rna_mutation <- clinical_rppa_rna_mutation_cna


# Validate final output -------------------------------------------------------

abort_if_false(
   nrow(clinical_rppa_rna_mutation) > 0,
   "clinical_rppa_rna_mutation must be non-empty."
)

check_has_columns(
   "clinical_rppa_rna_mutation",
   c("sample_id", "os_months", "os_event")
)


# Summarise integrated dataset ------------------------------------------------

mutation_summary <- summarise_mutations(clinical_rppa_rna_mutation)

message(
   "Integration complete. ",
   nrow(clinical_rppa_rna_mutation),
   " samples fully combined into multi-omics workspace with zeroed wild-types."
)

print(mutation_summary)