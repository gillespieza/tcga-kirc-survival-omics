#' @file        03_prepare_rppa.R
#' @title       Prepare RPPA Proteomics Data
#' @description Reshapes the raw RPPA data from a wide feature-by-sample matrix
#'              into a tidy sample-by-feature tibble ready for modelling.
#'              Parses the composite feature identifier into gene symbol and
#'              protein feature components, sanitises feature names for use as
#'              column names, and pivots to one row per tumour sample.
#'
#' @details
#'   Input:
#'   \describe{
#'     \item{rppa_data}{Wide tibble with one row per antibody/protein feature
#'           and one column per tumour sample. The first column,
#'           \code{Composite.Element.REF}, identifies each feature in the
#'           format \code{GENE|Antibody-Name}.}
#'   }
#'
#'   Objects created in the global environment:
#'   \describe{
#'     \item{rppa_proteomics}{Sample-by-feature tibble. One row per tumour
#'           sample (\code{sample_id}) and one column per sanitised protein
#'           feature, containing RPPA z-scores.}
#'   }
#'
#'   Transformation steps:
#'   \enumerate{
#'     \item Validate input dimensions and feature name uniqueness.
#'     \item Pivot to long format (one row per sample-feature combination).
#'     \item Parse \code{Composite.Element.REF} into \code{gene_symbol} and
#'           \code{protein_feature} on the \code{|} separator.
#'     \item Sanitise feature names: replace non-alphanumeric characters with
#'           underscores and strip trailing underscores.
#'     \item Pivot back to wide format (one row per sample).
#'   }
#'
#' @pre  01_load_data.R must be sourced first so that \code{rppa_data} is
#'       available in the global environment.
#'
#' @post \code{rppa_proteomics} is available for use by downstream feature
#'       selection and modelling scripts.
#'
#' @usage source("03_prepare_rppa.R")


# Validate input --------------------------------------------------------------

if (!"Composite.Element.REF" %in% names(rppa_data)) {
   stop("rppa_data is missing expected column: Composite.Element.REF")
}

n_features <- nrow(rppa_data)
n_samples  <- ncol(rppa_data) - 1L   # exclude the Composite.Element.REF column

message("Input RPPA data: ", n_features, " features x ", n_samples, " samples.")


# Check for feature name collisions before pivoting ---------------------------
# Sanitise names here first so that any collisions are caught and reported
# clearly before they silently produce list columns in pivot_wider.

feature_names_clean <- rppa_data$Composite.Element.REF %>%
   stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
   stringr::str_remove("_+$")

duplicated_names <- feature_names_clean[duplicated(feature_names_clean)]

if (length(duplicated_names) > 0) {
   stop(
      length(duplicated_names), " duplicate sanitised feature name(s) detected. ",
      "pivot_wider would silently produce list columns. ",
      "Affected name(s): ", paste(unique(duplicated_names), collapse = ", ")
   )
}


# Prepare RPPA proteomics table -----------------------------------------------

# Convert RPPA from wide feature-by-sample to sample-by-feature:
#   1. Pivot to long (one row per sample-feature combination)
#   2. Parse Composite.Element.REF into gene_symbol | protein_feature
#   3. Sanitise feature names for use as column names
#   4. Pivot back to wide (one row per sample)

rppa_proteomics <- rppa_data %>%
   pivot_longer(
      cols      = -Composite.Element.REF,
      names_to  = "sample_id",
      values_to = "rppa_zscore"
   ) %>%
   separate(
      Composite.Element.REF,
      into   = c("gene_symbol", "protein_feature"),
      sep    = "\\|",
      remove = FALSE,
      extra  = "merge",
      fill   = "right"
   ) %>%
   mutate(
      # Sanitise protein feature names, make safe to use as column names later
      protein_feature_clean = Composite.Element.REF %>%
         stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
         stringr::str_remove("_+$")
   ) %>%
   select(sample_id, protein_feature_clean, rppa_zscore) %>%
   pivot_wider(
      names_from  = protein_feature_clean,
      values_from = rppa_zscore
   )


# Validate output -------------------------------------------------------------

stopifnot(
   "rppa_proteomics must be non-empty"      = nrow(rppa_proteomics) > 0,
   "rppa_proteomics must contain sample_id" = "sample_id" %in% names(rppa_proteomics)
)

# Report NA coverage so sparse or missing data surfaces early.
n_na    <- sum(is.na(select(rppa_proteomics, -sample_id)))
n_total <- nrow(rppa_proteomics) * (ncol(rppa_proteomics) - 1L)
pct_na  <- round(100 * n_na / n_total, 1)

if (pct_na > 10) {
   warning(pct_na, "% of RPPA values are NA. Check for sample coverage issues.")
}

message(
   "RPPA proteomics table prepared: ",
   nrow(rppa_proteomics), " samples x ",
   ncol(rppa_proteomics) - 1L, " features. ",
   pct_na, "% missing values."
)