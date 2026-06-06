# Prepare Copy Number Alteration (CNA) features -------------------------------
#
# Extracts high-impact focal copy number alterations (deep deletions and high-
# level amplifications) for your curated panel of clear cell renal cell
# carcinoma driver genes. Applies a tiered prevalence filter to safeguard
# downstream survival cross-validation models from statistical separation
# and empty-matrix crashes.
#
# Requires: load_data.R to have been sourced, with cna_data available in
#           global memory.
#
# Produces:
#   cna_features - Wide binary matrix containing sample_id and prefixed columns
#                  (cna_loss_* and cna_gain_*).
#   cna_summary  - Tidy summary reporting cohort frequencies and prevalence.
#
# Outputs:
#   results/cna_prevalence_summary.csv
#   results/prepared_cna_features.csv
#
# Usage: this script is intended to be sourced by run_analysis.R.


# Validate Inputs -------------------------------------------------------------

check_required_objects(c("cna_data", "driver_genes"))

if (!all(c("Hugo_Symbol", "Entrez_Gene_Id") %in% names(cna_data))) {
  stop(
    "cna_data is missing required structural metadata columns.",
    call. = FALSE
  )
}


# 1. Reshape and Tidy Raw CNA Matrix -------------------------------------------

message("Processing and standardising raw Copy Number Alteration matrix...")

# Isolate columns that represent sample barcodes
sample_cols <- setdiff(names(cna_data), c("Hugo_Symbol", "Entrez_Gene_Id"))

cna_long_clean <- cna_data |>
  # Restrict rows strictly to your curated panel of core ccRCC driver genes
  dplyr::filter(.data$Hugo_Symbol %in% driver_genes) |>
  # Pivot wide patients into a clean long format vertical stream
  tidyr::pivot_longer(
    cols      = dplyr::all_of(sample_cols),
    names_to  = "sample_id",
    values_to = "gistic_score"
  ) |>
  # Standardise the barcodes to uniform character structures
  dplyr::mutate(
    sample_id    = standardise_sample_id(.data$sample_id),
    gistic_score = as.integer(.data$gistic_score)
  ) |>
  # Create distinct binary features for high-impact biological events
  dplyr::mutate(
    cna_loss = dplyr::if_else(
      !is.na(.data$gistic_score) & .data$gistic_score == -2L, 1L, 0L
    ),
    cna_gain = dplyr::if_else(
      !is.na(.data$gistic_score) & .data$gistic_score == 2L, 1L, 0L
    )
  ) |>
  # Pivot the discrete events into explicit, individual feature rows
  tidyr::pivot_longer(
    cols      = c("cna_loss", "cna_gain"),
    names_to  = "alteration_type",
    values_to = "status"
  ) |>
  # Construct clean, descriptive column headers
  dplyr::mutate(
    feature = paste0(.data$alteration_type, "_", .data$Hugo_Symbol)
  ) |>
  dplyr::select(sample_id, feature, status)


# 2. Execute Tiered Prevalence Filtering Engine --------------------------------

message("Evaluating cohort prevalence across available alterations...")

# Pre-calculate baseline statistics for all available driver features
cna_stats_base <- cna_long_clean |>
  dplyr::group_by(.data$feature) |>
  dplyr::summarise(
    n_altered = sum(.data$status == 1L, na.rm = TRUE),
    n_profiled = sum(!is.na(.data$status)),
    pct_altered = sum(.data$status == 1L, na.rm = TRUE) /
      sum(!is.na(.data$status)),
    .groups = "drop"
  )

# Tier 1: Target the standard robust 5% threshold
target_threshold <- 0.05
cna_features_filtered <- cna_stats_base |>
  dplyr::filter(.data$pct_altered >= target_threshold)

# Tier 2 Fallback: Relax to 2% if the strict threshold wipes out all features
if (nrow(cna_features_filtered) == 0L) {
  warning(
    "\u26a0\ufe0f No focal CNA features passed the strict 5% prevalence filter",
    call. = FALSE
  )
  message(
    "--> Engaging Tier 2 Fallback: Relaxing threshold to >= 2% to ",
    "retain available genomic signal."
  )

  target_threshold <- 0.02
  cna_features_filtered <- cna_stats_base |>
    dplyr::filter(.data$pct_altered >= target_threshold)
}

# Tier 3 Ultimate Fallback: Retain the single most frequent alteration
if (nrow(cna_features_filtered) == 0L) {
  warning(
    "\u26a0\ufe0f No focal CNA features passed relaxed 2% prevalence filter",
    call. = FALSE
  )
  message(
    "--> Engaging Tier 3 Ultimate Fallback: Retaining the single most ",
    "frequent alteration to maintain matrix structure."
  )

  cna_features_filtered <- cna_stats_base |>
    dplyr::arrange(dplyr::desc(.data$n_altered)) |>
    dplyr::slice_head(n = 1L)
}

# Pull the final vector of validated feature names
surviving_cna_features <- cna_features_filtered |>
  dplyr::pull(feature)


# 3. Create Wide Analysis Matrix & Summary Reports -----------------------------

# Construct the wide sample-by-feature matrix using the non-empty vector
cna_features <- cna_long_clean |>
  dplyr::filter(.data$feature %in% surviving_cna_features) |>
  tidyr::pivot_wider(
    names_from  = "feature",
    values_from = "status",
    values_fn   = max # Resolves duplicate rows if present defensively
  )

# Compile a clean summary reporting frame for your results chapter
cna_summary <- cna_stats_base |>
  dplyr::filter(.data$feature %in% surviving_cna_features) |>
  dplyr::transmute(
    feature    = .data$feature,
    n_altered  = .data$n_altered,
    pct_preval = round(.data$pct_altered * 100, 1)
  ) |>
  dplyr::arrange(dplyr::desc(.data$n_altered))


# Save and Export Clean Data Layers --------------------------------------------

readr::write_csv(cna_summary, "results/cna_prevalence_summary.csv")
readr::write_csv(cna_features, "results/prepared_cna_features.csv")

message(
  "CNA preparation complete. Isolated ",
  length(surviving_cna_features), " stable focal feature(s)."
)
print(cna_summary)
