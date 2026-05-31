#' @file        05_integrate_data.R
#' @title       Integrate Clinical, RPPA, and Mutation Data
#' @description Combines the prepared clinical survival, RPPA proteomics, and
#'              binary mutation feature tables into a single analysis-ready
#'              tibble. Clinical and RPPA tables are inner-joined on sample ID
#'              to retain only samples with both data types. Mutation features
#'              are left-joined so that samples absent from the mutation table
#'              (i.e. wild-type for all selected driver genes) are retained and
#'              coded as 0.
#'
#' @details
#'   Inputs:
#'   \describe{
#'     \item{clinical_survival}{Analysis-ready survival tibble. One row per
#'           sample. Produced by 02_prepare_clinical.R.}
#'     \item{rppa_proteomics}{Sample-by-feature RPPA z-score tibble. One row
#'           per sample. Produced by 03_prepare_rppa.R.}
#'     \item{mutation_features}{Sample-by-gene binary mutation tibble. One row
#'           per sample with at least one driver gene mutation. Produced by
#'           04_prepare_mutations.R.}
#'   }
#'
#'   Objects created in the global environment:
#'   \describe{
#'     \item{clinical_rppa}{Tibble combining clinical survival variables and
#'           RPPA protein z-score features. One row per sample present in
#'           both \code{clinical_survival} and \code{rppa_proteomics}
#'           (inner join).}
#'     \item{clinical_rppa_mutation}{Fully integrated analysis tibble.
#'           Extends \code{clinical_rppa} with binary driver gene mutation
#'           features. Samples absent from \code{mutation_features} are
#'           retained with their mutation columns set to 0 (not mutated).}
#'   }
#'
#' @pre  02_prepare_clinical.R, 03_prepare_rppa.R, and
#'       04_prepare_mutations.R must be sourced first so that
#'       \code{clinical_survival}, \code{rppa_proteomics}, and
#'       \code{mutation_features} are available in the global environment.
#'
#' @post \code{clinical_rppa_mutation} is the primary analysis dataset for
#'       downstream modelling scripts.
#'
#' @usage source("05_integrate_data.R")


# Validate inputs -------------------------------------------------------------

required_tables <- c("clinical_survival", "rppa_proteomics", "mutation_features")
missing_tables  <- required_tables[!sapply(required_tables, exists)]

if (length(missing_tables) > 0) {
   stop(
      "The following required table(s) are missing from the environment: ",
      paste(missing_tables, collapse = ", "),
      ". Source the relevant prepare scripts first."
   )
}

if (!"sample_id" %in% names(clinical_survival)) {
   stop("clinical_survival is missing required join key: sample_id")
}
if (!"sample_id" %in% names(rppa_proteomics)) {
   stop("rppa_proteomics is missing required join key: sample_id")
}
if (!"sample_id" %in% names(mutation_features)) {
   stop("mutation_features is missing required join key: sample_id")
}


# Integrate clinical and RPPA proteomics data ---------------------------------
# Inner join retains only samples present in both tables. Samples in
# clinical_survival with no RPPA data (and vice versa) are dropped.

clinical_rppa <- clinical_survival %>%
   inner_join(rppa_proteomics, by = "sample_id")

n_dropped_rppa <- nrow(clinical_survival) - nrow(clinical_rppa)
if (n_dropped_rppa > 0) {
   message(
      n_dropped_rppa, " sample(s) in clinical_survival had no matching ",
      "RPPA data and were dropped by the inner join."
   )
}

message("Clinical + RPPA: ", nrow(clinical_rppa), " rows x ",
        ncol(clinical_rppa), " cols.")


# Integrate mutation features -------------------------------------------------
# Left join retains all rows from clinical_rppa. Samples absent from
# mutation_features had no observed mutation in any selected driver gene
# and are valid wild-type samples; their mutation columns are set to 0.

clinical_rppa_mutation <- clinical_rppa %>%
   left_join(mutation_features, by = "sample_id") %>%
   mutate(
      across(
         starts_with("mut_"),
         ~ replace_na(.x, 0L)
      )
   )

# The left join should never change the row count. Warn if it does, since
# this would indicate duplicate sample IDs in mutation_features.
if (nrow(clinical_rppa_mutation) != nrow(clinical_rppa)) {
   warning(
      "Row count changed after mutation left join: ",
      nrow(clinical_rppa), " -> ", nrow(clinical_rppa_mutation), " rows. ",
      "mutation_features may contain duplicate sample IDs."
   )
}

message("Clinical + RPPA + mutation: ", nrow(clinical_rppa_mutation), " rows x ",
        ncol(clinical_rppa_mutation), " cols.")


# Validate output -------------------------------------------------------------

stopifnot(
   "clinical_rppa_mutation must be non-empty"      = nrow(clinical_rppa_mutation) > 0,
   "clinical_rppa_mutation must contain sample_id" = "sample_id" %in% names(clinical_rppa_mutation),
   "clinical_rppa_mutation must contain os_months" = "os_months" %in% names(clinical_rppa_mutation),
   "clinical_rppa_mutation must contain os_event"  = "os_event"  %in% names(clinical_rppa_mutation)
)


# Summarise integrated dataset ------------------------------------------------
# Mutation counts are reported as a ranked long table (consistent with
# 04_prepare_mutations.R). Wild-type sample count is reported separately.

n_wt_samples <- sum(
   rowSums(dplyr::select(clinical_rppa_mutation, starts_with("mut_"))) == 0
)

mutation_summary <- clinical_rppa_mutation %>%
   dplyr::summarise(dplyr::across(dplyr::starts_with("mut_"), sum)) %>%
   tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "gene",
      values_to = "n_mutated"
   ) %>%
   dplyr::mutate(
      gene        = stringr::str_remove(.data$gene, "^mut_"),
      pct_mutated = round(100 * .data$n_mutated / nrow(clinical_rppa_mutation), 1)
   ) %>%
   dplyr::arrange(dplyr::desc(.data$n_mutated))

message(
   "Integration complete. ",
   nrow(clinical_rppa_mutation), " samples ready for modelling. ",
   n_wt_samples, " sample(s) are wild-type for all selected driver genes."
)

print(mutation_summary)