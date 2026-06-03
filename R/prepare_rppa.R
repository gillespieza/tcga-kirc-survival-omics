# Prepare RPPA proteomics data ------------------------------------------------
#
# Reshapes the raw RPPA data from a wide feature-by-sample matrix into a tidy
# sample-by-feature tibble ready for modelling. Parses the composite feature
# identifier into gene symbol and protein feature components, sanitises feature
# names for use as column names, and pivots to one row per tumour sample.
#
# Requires: load_data.R to have been sourced so that rppa_data is available in
#           the global environment.
#
# Produces:
#   rppa_proteomics - Sample-by-feature tibble. One row per tumour sample
#                     (sample_id) and one column per sanitised protein feature,
#                     containing RPPA z-scores.
#
# Transformation steps:
#   1. Validate input dimensions and feature name uniqueness.
#   2. Pivot to long format (one row per sample-feature combination).
#   3. Parse Composite.Element.REF into gene_symbol and protein_feature on
#      the | separator.
#   4. Sanitise feature names: replace non-alphanumeric characters with
#      underscores and strip trailing underscores.
#   5. Pivot back to wide format (one row per sample).
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.

# Validate input --------------------------------------------------------------

if (!"Composite.Element.REF" %in% names(rppa_data)) {
   stop(
      "rppa_data is missing expected column: Composite.Element.REF",
      call. = FALSE
   )
}

n_features <- nrow(rppa_data)
n_samples  <- ncol(rppa_data) - 1L

message(
   "Input RPPA data: ",
   n_features, " features x ",
   n_samples, " samples."
)

# Check for feature name collisions before pivoting ---------------------------
# Sanitise names here first so that any collisions are caught and reported
# clearly before they silently produce list columns in pivot_wider.

feature_names_clean <- rppa_data$Composite.Element.REF |>
   stringr::str_replace_all("[^A-Za-z0-9]+", "_") |>
   stringr::str_remove("_+$")

duplicated_names <- feature_names_clean[duplicated(feature_names_clean)]

if (length(duplicated_names) > 0) {
   stop(
      length(duplicated_names), " duplicate sanitised feature name(s) detected. ",
      "pivot_wider would silently produce list columns. ",
      "Affected name(s): ", paste(unique(duplicated_names), collapse = ", "),
      call. = FALSE
   )
}

# Prepare RPPA proteomics table -----------------------------------------------
# Convert RPPA from wide feature-by-sample to sample-by-feature:
#   1. Pivot to long (one row per sample-feature combination).
#   2. Parse Composite.Element.REF into gene_symbol | protein_feature.
#   3. Sanitise feature names for use as column names.
#   4. Pivot back to wide (one row per sample).

rppa_proteomics <- rppa_data |>
   tidyr::pivot_longer(
      cols      = -Composite.Element.REF,
      names_to  = "sample_id",
      values_to = "rppa_zscore"
   ) |>
   tidyr::separate(
      Composite.Element.REF,
      into   = c("gene_symbol", "protein_feature"),
      sep    = "\\|",
      remove = FALSE,
      extra  = "merge",
      fill   = "right"
   ) |>
   dplyr::mutate(
      protein_feature_clean = Composite.Element.REF |>
         stringr::str_replace_all("[^A-Za-z0-9]+", "_") |>
         stringr::str_remove("_+$")
   ) |>
   dplyr::select(sample_id, protein_feature_clean, rppa_zscore) |>
   tidyr::pivot_wider(
      names_from  = protein_feature_clean,
      values_from = rppa_zscore
   )

# Validate output -------------------------------------------------------------

if (nrow(rppa_proteomics) == 0) {
   stop("rppa_proteomics must be non-empty", call. = FALSE)
}

if (!"sample_id" %in% names(rppa_proteomics)) {
   stop("rppa_proteomics must contain sample_id", call. = FALSE)
}

# Report NA coverage so sparse or missing data surfaces early.
n_na    <- sum(is.na(dplyr::select(rppa_proteomics, -sample_id)))
n_total <- nrow(rppa_proteomics) * (ncol(rppa_proteomics) - 1L)
pct_na  <- round(100 * n_na / n_total, 1)

if (pct_na > 10) {
   warning(
      pct_na, "% of RPPA values are NA. ",
      "Check for sample coverage issues.",
      call. = FALSE
   )
}

message(
   "RPPA proteomics table prepared: ",
   nrow(rppa_proteomics), " samples x ",
   ncol(rppa_proteomics) - 1L, " features. ",
   pct_na, "% missing values."
)