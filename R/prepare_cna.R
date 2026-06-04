# Prepare binary Copy Number Alteration (CNA) features -----------------------
#
# Derives a sample-by-feature binary CNA table from the raw gene-level matrix.
# Restricts analysis to the curated set of established ccRCC driver genes,
# isolates high-impact focal events (deep homozygous deletions [-2] and
# high-level amplifications [2]), and applies a prevalence filter to keep
# the downstream model parameter space compact and stable.
#
# Requires: load_data.R to have been sourced so that cna_data is available
#           in the global environment.
#
# Produces:
#   cna_features - Sample-by-feature binary tibble. One row per tumour sample
#                  (sample_id) and columns representing high-impact copy
#                  number changes (prefixed with cna_loss_ or cna_gain_).
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.


# Validate input --------------------------------------------------------------

check_required_objects("cna_data")

check_has_columns("cna_data", "Hugo_Symbol")

abort_if_false(
  nrow(cna_data) > 0,
  "cna_data has zero rows after loading. Check data_cna.txt."
)

# Number of samples is total columns minus the Hugo_Symbol column
n_samples_raw <- ncol(cna_data) - 1L

message(
  "Input CNA data: ",
  nrow(cna_data), " genes x ",
  n_samples_raw, " samples."
)


# Define driver genes ---------------------------------------------------------
# Aligned exactly with the mutation layer to capture core ccRCC biology

driver_genes <- c(
  "VHL",
  "PBRM1",
  "BAP1",
  "SETD2",
  "KDM5C",
  "MTOR",
  "PTEN",
  "TSC1",
  "TSC2"
)


# Reshape matrix and extract focal alterations --------------------------------
# Reshapes wide gene-by-sample matrix into a clean, long format tibble,
# safely ignoring any Entrez ID metadata columns.

cna_long <- cna_data |>
   dplyr::filter(.data$Hugo_Symbol %in% driver_genes) |>
   dplyr::select(-dplyr::any_of("Entrez_Gene_Id")) |>
   tidyr::pivot_longer(
      cols      = -dplyr::all_of("Hugo_Symbol"),
      names_to  = "sample_id",
      values_to = "cna_status"
   ) |>
   tidyr::drop_na(.data$cna_status)


# Derive binary features ------------------------------------------------------
# Isolate high-impact events: Deep Loss (-2) and High Gain (2)

cna_binary_loss <- cna_long |>
  dplyr::transmute(
    sample_id   = .data$sample_id,
    feature     = paste0("cna_loss_", .data$Hugo_Symbol),
    altered_val = as.integer(.data$cna_status == -2)
  )

cna_binary_gain <- cna_long |>
  dplyr::transmute(
    sample_id   = .data$sample_id,
    feature     = paste0("cna_gain_", .data$Hugo_Symbol),
    altered_val = as.integer(.data$cna_status == 2)
  )

# Combine losses and gains into a single table, then pivot wide
cna_features_all <- dplyr::bind_rows(cna_binary_loss, cna_binary_gain) |>
  tidyr::pivot_wider(
    names_from  = "feature",
    values_from = "altered_val",
    values_fill = 0L
  )


# Apply prevalence filter -----------------------------------------------------
# Drops alterations with low prevalence to safeguard downstream models
# against severe overfitting (EPV deficit)

cna_feature_cols <- setdiff(names(cna_features_all), "sample_id")

passed_prevalence_cols <- cna_feature_cols[
  vapply(
    cna_features_all[cna_feature_cols],
    function(x) mean(x, na.rm = TRUE) >= 0.02,
    logical(1)
  )
]

if (length(passed_prevalence_cols) == 0) {
  warning(
    "No focal CNA features passed the 2% prevalence filter. ",
    "Creating an empty placeholder feature table to maintain pipeline integrity.",
    call. = FALSE
  )
}

cna_features <- cna_features_all |>
  dplyr::select(
    dplyr::all_of("sample_id"),
    dplyr::all_of(passed_prevalence_cols)
  )


# Validate output -------------------------------------------------------------

abort_if_false(
  nrow(cna_features) > 0,
  "cna_features table must be non-empty."
)

check_has_sample_id("cna_features")


# Summarise CNA frequencies --------------------------------------------------

cna_summary <- cna_features |>
  dplyr::summarise(
    dplyr::across(
      -dplyr::all_of("sample_id"),
      sum
    )
  ) |>
  tidyr::pivot_longer(
    cols      = dplyr::everything(),
    names_to  = "alteration",
    values_to = "n_altered"
  ) |>
  dplyr::mutate(
    pct_altered = round(100 * .data$n_altered / nrow(cna_features), 1)
  ) |>
  dplyr::arrange(dplyr::desc(.data$n_altered))

message(
  "CNA feature table prepared: ",
  nrow(cna_features), " samples x ",
  ncol(cna_features) - 1L, " features passed screening."
)

print(cna_summary)
