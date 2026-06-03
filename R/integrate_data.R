# Integrate clinical, RPPA, RNA-seq, and mutation data ------------------------
#
# Combines the prepared clinical survival, RPPA proteomics, RNA-seq, and binary
# mutation feature tables into a single analysis-ready tibble.
#
# Order of integration:
#   1. Inner join clinical_survival with rppa_proteomics on sample_id.
#   2. Inner join that result with rnaseq_expression on sample_id.
#   3. Left join mutation_features on sample_id (so samples without mutation
#      calls are retained and coded as wild-type = 0).
#
# Requires: prepare_clinical.R, prepare_rppa.R, prepare_rnaseq.R, and
#           prepare_mutations.R to have been sourced so that clinical_survival,
#           rppa_proteomics, rnaseq_expression, and mutation_features are
#           available in the global environment.
#
# Produces:
#   clinical_rppa            - Clinical + RPPA table.
#   clinical_rppa_rna          - Clinical + RPPA + RNA-seq table.
#   clinical_rppa_rna_mutation - Fully integrated analysis tibble used
#                                downstream by modelling scripts.
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.


# Validate inputs -------------------------------------------------------------

required_tables <- c(
   "clinical_survival",
   "rppa_proteomics",
   "rnaseq_expression",
   "mutation_features"
)

check_required_objects(required_tables)

for (tbl in required_tables) {
   check_has_sample_id(tbl)
}

# 1. Integrate clinical and RPPA ----------------------------------------------
# Inner join retains only samples present in both tables.

common_samples_clin_rppa <- intersect(
   clinical_survival$sample_id,
   rppa_proteomics$sample_id
)

clinical_rppa <- clinical_survival |>
   dplyr::inner_join(
      rppa_proteomics,
      by = "sample_id"
   )

n_dropped_rppa <- length(setdiff(
   clinical_survival$sample_id,
   clinical_rppa$sample_id
))

if (n_dropped_rppa > 0) {
   message(
      n_dropped_rppa,
      " sample(s) in clinical_survival had no matching RPPA data ",
      "and were dropped by the inner join."
   )
}

message(
   "Clinical + RPPA: ",
   nrow(clinical_rppa), " samples x ",
   ncol(clinical_rppa), " features."
)


# 2. Integrate RNA-seq --------------------------------------------------------
# RNA-seq is high-dimensional and may have slightly different sample coverage.
# We restrict to the intersection of clinical + RPPA + RNA-seq samples.

common_samples_all_omics <- intersect(
   clinical_rppa$sample_id,
   rnaseq_expression$sample_id
)

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
# Left join retains all samples with full omics profile.
# Missing mutation values are treated as wild-type (0).

clinical_rppa_rna_mutation <- clinical_rppa_rna |>
   dplyr::left_join(
      mutation_features,
      by = "sample_id"
   ) |>
   dplyr::mutate(
      dplyr::across(
         dplyr::starts_with("mut_"),
         ~ tidyr::replace_na(.x, 0L)
      )
   )

abort_if_false(
   nrow(clinical_rppa_rna_mutation) == nrow(clinical_rppa_rna),
   paste(
      "Row count changed after mutation join:",
      nrow(clinical_rppa_rna), "->",
      nrow(clinical_rppa_rna_mutation),
      ". Check mutation_features for duplicate sample_id values."
   )
)

message(
   "Clinical + RPPA + RNA-seq + mutation: ",
   nrow(clinical_rppa_rna_mutation), " samples x ",
   ncol(clinical_rppa_rna_mutation), " features."
)

# Validate output -------------------------------------------------------------

abort_if_false(
   nrow(clinical_rppa_rna_mutation) > 0,
   "clinical_rppa_rna_mutation must be non-empty."
)

abort_if_false(
   "sample_id" %in% names(clinical_rppa_rna_mutation),
   "clinical_rppa_rna_mutation must contain sample_id."
)

abort_if_false(
   "os_months" %in% names(clinical_rppa_rna_mutation),
   "clinical_rppa_rna_mutation must contain os_months."
)

abort_if_false(
   "os_event" %in% names(clinical_rppa_rna_mutation),
   "clinical_rppa_rna_mutation must contain os_event."
)


# Summarise integrated dataset ------------------------------------------------
# Mutation counts are reported as a ranked long table so you can compare with
# the standalone mutation summary.

n_wt_samples <- sum(
   rowSums(
      dplyr::select(
         clinical_rppa_rna_mutation,
         dplyr::starts_with("mut_")
      )
   ) == 0
)

mutation_summary <- clinical_rppa_rna_mutation |>
   dplyr::summarise(
      dplyr::across(
         dplyr::starts_with("mut_"),
         sum
      )
   ) |>
   tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "gene",
      values_to = "n_mutated"
   ) |>
   dplyr::mutate(
      gene        = stringr::str_remove(.data$gene, "^mut_"),
      pct_mutated = round(
         100 * .data$n_mutated / nrow(clinical_rppa_rna_mutation),
         1
      )
   ) |>
   dplyr::arrange(
      dplyr::desc(.data$n_mutated)
   )

message(
   "Integration complete. ",
   nrow(clinical_rppa_rna_mutation),
   " samples ready for modelling. ",
   n_wt_samples,
   " sample(s) are wild-type for all selected driver genes."
)

print(mutation_summary)