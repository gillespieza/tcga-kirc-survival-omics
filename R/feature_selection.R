# RPPA feature selection -------------------------------------------------------
#
# Screens RPPA protein expression features with univariable Cox models to
# support pathway correlation analysis, then applies a multivariable LASSO
# Cox model on a promising subset to isolate stable prognostic proteins for
# survival modelling. Mutation features are carried forward from the curated
# ccRCC driver-gene set prepared in prepare_mutations.R, after filtering for
# minimum prevalence.
#
# Requires: prepare_mutations.R and prepare_clinical.R to have been sourced,
#            and survival_data to be available in the global environment.
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
#   results/rppa_lasso_coefficients.csv
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.


# Validate Inputs ------------------------------------------------

check_required_objects(
  c(
    "survival_data",
    "selected_clinical_features",
    "clinical_candidate_cols"
  )
)

# Validate outcome columns only — clinical predictors are dynamic.
required_survival_cols <- c(
  "patient_id",
  "sample_id",
  "os_months",
  "os_event"
)

check_has_columns(
  "survival_data",
  c(required_survival_cols, selected_clinical_features)
)

# Column Classification --------------------------------------------------------
# Separate the individual feature blocks cleanly by parsing column name prefixes

all_col_names <- names(survival_data)

mutation_feature_cols <- all_col_names[
  stringr::str_starts(all_col_names, "mut_")
]

rna_feature_cols <- all_col_names[
  stringr::str_starts(all_col_names, "rna_")
]

# Exclude outcome/ID columns AND all clinical candidate columns from the RPPA
# candidate set. clinical_candidate_cols (the full list from prepare_clinical.R)
# is used here rather than selected_clinical_features so that numeric clinical
# variables such as aneuploidy_score and tmb_nonsynonymous are not mistakenly
# screened as RPPA proteins.
non_rppa_cols <- c(
  required_survival_cols,
  clinical_candidate_cols,
  mutation_feature_cols,
  rna_feature_cols
)

rppa_feature_cols <- setdiff(all_col_names, non_rppa_cols)


# Filter RPPA Candidates -------------------------------------------------------
# Retain numeric proteomic variables that carry adequate sample coverage

rppa_feature_cols <- rppa_feature_cols[
  vapply(
    survival_data[rppa_feature_cols],
    function(x) {
      is.numeric(x) &&
        sum(!is.na(x)) >= 30L &&
        stats::sd(x, na.rm = TRUE) > 0L
    },
    logical(1L)
  )
]

abort_if_false(
  length(rppa_feature_cols) > 0L,
  "No usable numeric RPPA features were found for Cox screening."
)

# Drop rare somatic mutations with low prevalence (< 5%) to protect power
mutation_feature_cols <- mutation_feature_cols[
  vapply(
    survival_data[mutation_feature_cols],
    function(x) mean(x, na.rm = TRUE) >= 0.05,
    logical(1L)
  )
]

if (length(mutation_feature_cols) == 0L) {
  warning(
    "\u26a0\ufe0f No mutation feature columns passed the prevalence filter; ",
    "mutation models will be skipped.",
    call. = FALSE
  )
}


# 1. Univariable Cox Screening Loop --------------------------------------------
# Evaluates baseline prognostic value for each protein to support pathway scores

message("Running univariable screening across all valid proteomic markers ...")

rppa_results_list <- list()

for (feat in rppa_feature_cols) {
  # Build a narrow complete case subset for this single protein to maximise
  # sample retention
  df_temp <- survival_data |>
    dplyr::select(dplyr::all_of(c("os_months", "os_event", feat))) |>
    tidyr::drop_na()

  if (nrow(df_temp) >= 30L) {
    formula_uni <- stats::as.formula(
      paste("survival::Surv(os_months, os_event) ~", feat)
    )

    fit_uni <- tryCatch(
      survival::coxph(formula_uni, data = df_temp, ties = "efron"),
      error = function(e) NULL
    )

    if (!is.null(fit_uni)) {
      summary_uni <- summary(fit_uni)
      coef_info <- summary_uni$coefficients

      rppa_results_list[[feat]] <- tibble::tibble(
        feature    = feat,
        log_hr     = coef_info[1L, "coef"],
        abs_log_hr = abs(coef_info[1L, "coef"]),
        p_value    = coef_info[1L, "Pr(>|z|)"]
      )
    }
  }
}

rppa_univariable_results <- dplyr::bind_rows(rppa_results_list) |>
  dplyr::mutate(
    p_adjust_bh = stats::p.adjust(.data$p_value, method = "BH")
  ) |>
  dplyr::arrange(.data$p_adjust_bh, dplyr::desc(.data$abs_log_hr))

# Export the univariable reference matrix required by pathway_scores.R
readr::write_csv(
  rppa_univariable_results,
  "results/rppa_univariable_cox_results.csv"
)


# 2. Multivariable LASSO Feature Selection -------------------------------------
# Restrict the model space to the top 30 univariable proteins to avoid a
# sample size crash

top_uni_proteins <- rppa_univariable_results |>
  dplyr::slice_head(n = 30L) |>
  dplyr::pull(feature)

# Use the empirically selected clinical features rather than a hardcoded list.
clinical_vars <- selected_clinical_features

modelling_vars <- c(clinical_vars, top_uni_proteins)

lasso_data <- survival_data |>
  dplyr::select(dplyr::all_of(c("os_months", "os_event", modelling_vars))) |>
  tidyr::drop_na()

message(
  "Number of complete cases for multivariable regularisation: ",
  nrow(lasso_data)
)

# Construct design matrix including dummy expansion strings for factor variables
x_matrix <- stats::model.matrix(
  stats::as.formula(paste("~", paste(modelling_vars, collapse = " + "))),
  data = lasso_data
)[, -1L]

x_names <- colnames(x_matrix)

# Map clinical vs proteomic column position indices
clinical_indices <- which(!x_names %in% top_uni_proteins)

y_surv <- survival::Surv(
  time  = lasso_data$os_months,
  event = lasso_data$os_event
)

# Set penalty factors: 0.001 acts as a ridge stabiliser for clinical baseline
# indicators
penalty_vector <- rep(1L, ncol(x_matrix))
penalty_vector[clinical_indices] <- 0.001

message("Fitting cross-validated multivariable LASSO Cox model ...")

lasso_cv_fit <- glmnet::cv.glmnet(
  x              = x_matrix,
  y              = y_surv,
  family         = "cox",
  penalty.factor = penalty_vector,
  nfolds         = 10L,
  cox.ties       = "efron"
)


# 3. Tiered Feature Extraction & Fail-Safes ------------------------------------

# Tier A: Attempt parsimonious extraction via standard lambda.1se
lasso_coefs <- stats::coef(lasso_cv_fit, s = "lambda.1se") |> as.matrix()

lasso_results_df <- tibble::tibble(
  feature     = rownames(lasso_coefs),
  coefficient = as.numeric(lasso_coefs)
) |>
  dplyr::filter(
    .data$coefficient != 0L,
    !.data$feature %in% x_names[clinical_indices]
  ) |>
  dplyr::arrange(dplyr::desc(abs(.data$coefficient)))

selected_rppa_features <- lasso_results_df |> dplyr::pull(feature)

# Tier B Fallback: Revert to lambda.min if 1se is completely shrunk
if (length(selected_rppa_features) == 0L) {
  message(
    "LASSO shrunk all protein coefficients to zero at lambda.1se. ",
    "Reverting to lambda.min ..."
  )

  lasso_coefs_min <- stats::coef(lasso_cv_fit, s = "lambda.min") |> as.matrix()

  lasso_results_df <- tibble::tibble(
    feature     = rownames(lasso_coefs_min),
    coefficient = as.numeric(lasso_coefs_min)
  ) |>
    dplyr::filter(
      .data$coefficient != 0L,
      !.data$feature %in% x_names[clinical_indices]
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$coefficient)))

  selected_rppa_features <- lasso_results_df |> dplyr::pull(feature)
}

# Tier C Fallback: Revert to top 5 univariable features if clinical dominance
# is absolute
if (length(selected_rppa_features) == 0L) {
  warning(
    "\u26a0\ufe0f LASSO shrunk all protein coefficients to zero at lambda.min ",
    "due to clinical dominance. Falling back to the top 5 most robust ",
    "univariable prognostic markers to guarantee pipeline continuity.",
    call. = FALSE
  )

  lasso_results_df <- rppa_univariable_results |>
    dplyr::slice_head(n = 5L) |>
    dplyr::transmute(
      feature     = .data$feature,
      coefficient = .data$log_hr
    )

  selected_rppa_features <- lasso_results_df |> dplyr::pull(feature)
}

# Apply an upper cap to eliminate selection bias and structural EPV deficits
if (length(selected_rppa_features) > 10L) {
  message(
    "Selection isolated more than 10 proteins. Capping at top 10 ",
    "to protect model power."
  )
  selected_rppa_features <- selected_rppa_features[1L:10L]
  lasso_results_df <- lasso_results_df |>
    dplyr::filter(.data$feature %in% selected_rppa_features)
}

selected_mutation_features <- mutation_feature_cols


# Write Outputs ---------------------------------------------------------------

readr::write_csv(
  lasso_results_df,
  "results/rppa_lasso_coefficients.csv"
)

readr::write_csv(
  tibble::tibble(feature = selected_rppa_features),
  "results/selected_rppa_features.csv"
)

readr::write_csv(
  tibble::tibble(feature = selected_mutation_features),
  "results/selected_mutation_features.csv"
)

message(
  "RPPA feature selection complete: Isolated ",
  length(selected_rppa_features),
  " robust prognostic protein markers."
)

print(selected_rppa_features)
