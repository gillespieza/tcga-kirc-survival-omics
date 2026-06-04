# Integrate clinical, RPPA, RNA-seq, mutations, and CNA data -----------------
#
# Combines the prepared clinical survival, RPPA proteomics, RNA-seq, binary
# mutation, and copy-number alteration tables into a single analysis-ready tibble.
#
# Order of integration:
#   1. Inner join clinical_survival with rppa_proteomics on sample_id.
#   2. Inner join that result with rnaseq_expression on sample_id.
#   3. Left join mutation_features on sample_id (handle missing entries safely).
#   4. Left join cna_features on sample_id (handle missing entries safely).
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
# RNA-seq is high-dimensional and may have slightly different sample coverage.
# We restrict to the intersection of clinical + RPPA + RNA-seq samples.

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

message(
   "Clinical + RPPA + RNA-seq: ",
   nrow(clinical_rppa_rna), " samples x ",
   ncol(clinical_rppa_rna), " features."
)


# 3. Integrate mutation features ----------------------------------------------
# Missing values are only treated as wild-type (0) if the sample is confirmed 
# to have undergone genomic sequencing.

sequenced_sample_universe <- unique(mutation_data$Tumor_Sample_Barcode)

clinical_rppa_rna_mutation <- clinical_rppa_rna |>
   dplyr::left_join(
      mutation_features,
      by = "sample_id"
   ) |>
   dplyr::mutate(
      dplyr::across(
         dplyr::starts_with("mut_", ignore.case = FALSE),
         function(x) {
            dplyr::case_when(
               !is.na(x) ~ x,
               is.na(x) & .data$sample_id %in% sequenced_sample_universe ~ 0L,
               TRUE ~ NA_integer_
            )
         }
      )
   )


# 4. Integrate Copy Number Alteration (CNA) features -------------------------
# Missing values are only treated as unaltered (0) if the sample is confirmed
# to have successfully undergone copy-number profiling.

cna_sample_universe <- setdiff(names(cna_data), "Hugo_Symbol")

clinical_rppa_rna_mutation_cna <- clinical_rppa_rna_mutation |>
   dplyr::left_join(
      cna_features,
      by = "sample_id"
   ) |>
   dplyr::mutate(
      dplyr::across(
         dplyr::starts_with("cna_", ignore.case = FALSE),
         function(x) {
            dplyr::case_when(
               !is.na(x) ~ x,
               is.na(x) & .data$sample_id %in% cna_sample_universe ~ 0L,
               TRUE ~ NA_integer_
            )
         }
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
   " samples fully combined into multi-omics workspace."
)

print(mutation_summary)