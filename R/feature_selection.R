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
# Identify clinical, mutation, RNA, & RPPA feature columns. Clinical, mutation,
# and RNA columns are excluded from RPPA screening so this script focuses on
# protein (RPPA) features.

all_col_names <- names(survival_data)

clinical_cols <- required_survival_cols

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

# Multivariable LASSO Feature Selection ---------------------------------------
#
# Instead of biased univariable screening, we use a multivariable LASSO Cox 
# model with internal 10-fold cross-validation. Clinical variables are given 
# a tiny penalty factor of 0.001 to act as a ridge stabilizer, preventing 
# C++ convergence failures due to separation while ensuring they are never 
# regularised out of the model space.

# Prepare design matrix (X) and survival outcome (Y)
clinical_vars <- c("age", "sex", "stage", "grade")
modelling_vars <- c(clinical_vars, rppa_feature_cols)

lasso_data <- survival_data |>
   dplyr::select(dplyr::all_of(c("os_months", "os_event", modelling_vars))) |>
   tidyr::drop_na()

# Create the numeric design matrix including dummy codes for clinical factors
x_matrix <- stats::model.matrix(
   stats::as.formula(paste("~", paste(modelling_vars, collapse = " + "))),
   data = lasso_data
)[, -1]

x_names <- colnames(x_matrix)

# Map which columns belong to the proteomic features vs clinical variables
rppa_indices <- which(x_names %in% rppa_feature_cols)
clinical_indices <- which(!x_names %in% rppa_feature_cols)

y_surv <- survival::Surv(
   time  = lasso_data$os_months,
   event = lasso_data$os_event
)

# Set up penalty factors: 0.001 to stabilize baseline clinical features, 1 for proteins
penalty_vector <- rep(1, ncol(x_matrix))
penalty_vector[clinical_indices] <- 0.001

message("Fitting cross-validated multivariable LASSO Cox model...")

lasso_cv_fit <- glmnet::cv.glmnet(
   x              = x_matrix,
   y              = y_surv,
   family         = "cox",
   penalty.factor = penalty_vector,
   nfolds         = 10,
   cox.ties       = "efron"
)

# Extract coefficients at the optimal lambda value 
# 'lambda.1se' provides the most regularised, parsimonious model within 1 SE of the minimum
lasso_coefs <- stats::coef(
   lasso_cv_fit,
   s = "lambda.1se"
) |>
   as.matrix()

# Convert to a tidy data frame for filtering
lasso_results_df <- tibble::tibble(
   feature     = rownames(lasso_coefs),
   coefficient = as.numeric(lasso_coefs)
) |>
   dplyr::filter(
      .data$coefficient != 0,
      !.data$feature %in% x_names[clinical_indices]
   ) |>
   dplyr::arrange(dplyr::desc(abs(.data$coefficient)))

# Extract the selected protein feature names for downstream scripts
selected_rppa_features <- lasso_results_df |>
   dplyr::pull(.data$feature)

# Cap the maximum number of features to safeguard against the EPV deficit
n_selected_rppa <- length(selected_rppa_features)

if (n_selected_rppa > 10) {
   message("LASSO selected ", n_selected_rppa, " proteins. Capping at top 10 to protect model power.")
   selected_rppa_features <- lasso_results_df |>
      dplyr::slice_head(n = 10) |>
      dplyr::pull(.data$feature)
} else if (n_selected_rppa == 0) {
   warning("LASSO shrunk all protein coefficients to zero. Reverting to minimum lambda target.")
   
   lasso_coefs_min <- stats::coef(lasso_cv_fit, s = "lambda.min") |>
      as.matrix()
   
   selected_rppa_features <- tibble::tibble(
      feature = rownames(lasso_coefs_min),
      coefficient = as.numeric(lasso_coefs_min)
   ) |>
      dplyr::filter(
         .data$coefficient != 0,
         !.data$feature %in% x_names[clinical_indices]
      ) |>
      dplyr::slice_head(n = 5) |>
      dplyr::pull(.data$feature)
}

selected_mutation_features <- mutation_feature_cols


# Write outputs ---------------------------------------------------------------

# Save the non-zero regularisation coefficients
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
   "RPPA feature selection complete: LASSO isolated ",
   length(selected_rppa_features),
   " robust prognostic protein markers."
)

print(selected_rppa_features)
