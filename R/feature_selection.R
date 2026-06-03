# RPPA feature selection -------------------------------------------------------
#
# Screens RPPA protein expression features with univariable Cox models and keeps
# a compact, interpretable set of proteins for downstream survival modelling.
# Mutation features are carried forward from the curated ccRCC driver-gene set
# prepared in prepare_mutations.R, after filtering for minimum prevalence.
#
# Requires: quick_survival_check.R to have been sourced so that survival_data is
#           available in the global environment.
#
# Produces:
#   rppa_feature_cols          - Candidate RPPA feature column names.
#   mutation_feature_cols      - Curated binary mutation feature column names.
#   rppa_univariable_results   - Ranked univariable Cox screening results.
#   selected_rppa_features     - Top RPPA features used in model comparison.
#   selected_mutation_features - Mutation features used in model comparison.
#
# Outputs:
#   results/rppa_univariable_cox_results.csv
#   results/selected_rppa_features.csv
#   results/selected_mutation_features.csv
#
# Note: this script is intended to be sourced by run_analysis.R as part of
#       the full pipeline, not run directly.


# Validate inputs -------------------------------------------------------------

check_required_objects("survival_data")

required_survival_cols <- c(
   "patient_id",
   "sample_id",
   "os_months",
   "os_event",
   "age",
   "sex",
   "stage",
   "grade"
)

check_has_columns(
   "survival_data",
   required_survival_cols
)


# Column classification --------------------------------------------------------
# Identify clinical, mutation, RNA, and RPPA feature columns. Clinical, mutation,
# and RNA columns are excluded from RPPA screening so this script focuses on
# protein (RPPA) features.

all_col_names <- names(survival_data)

mutation_feature_cols <- all_col_names[
   stringr::str_starts(all_col_names, "mut_")
]

rna_feature_cols <- all_col_names[
   stringr::str_starts(all_col_names, "rna_")
]

rppa_feature_cols <- setdiff(
   all_col_names,
   c(clinical_cols, mutation_feature_cols, rna_feature_cols)
)

# Filter RPPA candidates -------------------------------------------------------
# Keep numeric RPPA columns with enough observed values and variation to support
# a Cox model.

rppa_feature_cols <- rppa_feature_cols[
   vapply(
      survival_data[rppa_feature_cols],
      function(x) {
         is.numeric(x) &&
            sum(!is.na(x)) >= 30 &&
            stats::sd(x, na.rm = TRUE) > 0
      },
      logical(1)
   )
]

abort_if_false(
   length(rppa_feature_cols) > 0,
   "No usable numeric RPPA features were found for Cox screening."
)

# Drop mutation features with prevalence below 5 % to avoid overfitting
# downstream on nearly-absent alterations.

mutation_feature_cols <- mutation_feature_cols[
   vapply(
      survival_data[mutation_feature_cols],
      function(x) mean(x, na.rm = TRUE) >= 0.05,
      logical(1)
   )
]

if (length(mutation_feature_cols) == 0) {
   warning(
      "No mutation feature columns passed the prevalence filter; ",
      "mutation models will be skipped.",
      call. = FALSE
   )
}

# Univariable Cox screen ------------------------------------------------------
# Rename the feature column inside the model data to avoid backtick quoting,
# which is fragile when feature names contain special characters.

fit_univariable_rppa <- function(feature_name) {
   model_data <- survival_data |>
      dplyr::select(
         os_months,
         os_event,
         feature = dplyr::all_of(feature_name)
      ) |>
      tidyr::drop_na()
   
   if (nrow(model_data) < 30 || length(unique(model_data$feature)) < 2) {
      return(NULL)
   }
   
   fit <- tryCatch(
      survival::coxph(
         survival::Surv(os_months, os_event) ~ feature,
         data = model_data,
         ties = "efron"
      ),
      error = function(e) NULL
   )
   
   if (is.null(fit)) {
      return(NULL)
   }
   
   broom::tidy(
      fit,
      conf.int     = TRUE,
      exponentiate = TRUE
   ) |>
      dplyr::slice(1) |>
      dplyr::transmute(
         feature      = feature_name,
         n            = nrow(model_data),
         events       = sum(model_data$os_event),  # <- use model_data here
         hazard_ratio = estimate,
         conf_low     = conf.low,
         conf_high    = conf.high,
         p_value      = p.value
      )
}

rppa_univariable_results <- purrr::map_dfr(
   rppa_feature_cols,
   fit_univariable_rppa
) |>
   dplyr::mutate(
      p_adjust_bh = stats::p.adjust(p_value, method = "BH"),
      abs_log_hr  = abs(log(hazard_ratio))
   ) |>
   dplyr::arrange(
      p_adjust_bh,
      dplyr::desc(abs_log_hr)
   )

abort_if_false(
   nrow(rppa_univariable_results) > 0,
   "RPPA univariable screening did not produce any valid Cox models."
)

# Select top features that survive FDR correction (BH-adjusted p < 0.05),
# capped at 10. Selection is based on adjusted p-value so the threshold is
# applied before ranking, not after.

n_selected_rppa <- 10L

selected_rppa_features <- rppa_univariable_results |>
   dplyr::filter(p_adjust_bh < 0.05) |>
   dplyr::slice_head(n = n_selected_rppa) |>
   dplyr::pull(feature)

if (length(selected_rppa_features) == 0) {
   warning(
      "No RPPA features passed the FDR threshold (BH-adjusted p < 0.05). ",
      "Consider relaxing the threshold or reviewing data quality.",
      call. = FALSE
   )
} else if (length(selected_rppa_features) < n_selected_rppa) {
   message(
      "Fewer than ", n_selected_rppa,
      " RPPA features passed the FDR threshold; ",
      length(selected_rppa_features), " selected."
   )
}

selected_mutation_features <- mutation_feature_cols


# Write outputs ---------------------------------------------------------------

readr::write_csv(
   rppa_univariable_results,
   "results/rppa_univariable_cox_results.csv"
)

readr::write_csv(
   rppa_univariable_results |>
      dplyr::filter(feature %in% selected_rppa_features) |>
      dplyr::select(
         feature,
         hazard_ratio,
         conf_low,
         conf_high,
         p_value,
         p_adjust_bh
      ),
   "results/selected_rppa_features.csv"
)

readr::write_csv(
   tibble::tibble(feature = selected_mutation_features),
   "results/selected_mutation_features.csv"
)

message(
   "RPPA feature selection complete: selected ",
   length(selected_rppa_features), " of ",
   length(rppa_feature_cols), " usable RPPA features."
)

message("Selected RPPA features:")
print(selected_rppa_features)