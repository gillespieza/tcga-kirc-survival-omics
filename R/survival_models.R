# Cross-validated multi-omics survival modelling ------------------------------
#
# Evaluates and compares the out-of-fold predictive performance of seven
# distinct survival models using stratified 5-fold cross-validation with
# fold-level invariance guards to prevent statistical separation crashes:
#   1. Clinical Baseline (Age, Sex, Stage, Grade)
#   2. Clinical + Mutations
#   3. Clinical + CNA
#   4. Clinical + RPPA Proteomics (Top 5 LASSO-retained markers)
#   5. Clinical + RNA Pathways (All 8 independent MSigDB signatures)
#   6. Clinical + RNA Data-Driven Score (Top 20 univariable p-value genes)
#   7. Unified Multi-Omics Integrated Model (Clinical + Mutations + CNA +
#      RPPA + RNA Pathways)
#
# Requires: integrate_data.R, feature_selection.R, and pathway_scores.R
#            to have been sourced, with survival_data available in global
#            memory.
#
# Produces:
#   cv_cindex_results    - Comparative data frame of cross-validated C-indices
#   final_integrated_fit - Master Cox object trained on full cohort
#
# Outputs:
#   results/cross_validated_cindex_comparison.csv
#   results/final_integrated_cox_coefficients.csv
#
# Usage: this script is intended to be sourced by run_analysis.R.


# Validate inputs --------------------------------------------------------------

check_required_objects(c(
  "survival_data",
  "selected_clinical_features",
  "selected_rppa_features",
  "selected_mutation_features"
))

required_base_cols <- c("os_months", "os_event", selected_clinical_features)
check_has_columns("survival_data", required_base_cols)


# Feature space discovery ------------------------------------------------------
# Dynamically extract available multi-omics predictors using prefixes

clinical_vars <- selected_clinical_features

mutation_vars <- selected_mutation_features

cna_vars <- names(survival_data)[stringr::str_starts(
  names(survival_data),
  "cna_"
)]

rppa_vars <- selected_rppa_features

# Pull all 8 curated individual PC1 pathway scores (excluding means and
# data-driven scores)
rna_vars <- names(survival_data)[stringr::str_starts(
  names(survival_data),
  "score_"
) &
  !stringr::str_ends(names(survival_data), "_mean") &
  names(survival_data) != "score_rna_datadriven"]

# Isolate the standalone data-driven PC1 score
rna_datadriven_vars <- "score_rna_datadriven"

message("=== Modelling Feature Dimensions ===")
message("Clinical Covariates    : ", length(clinical_vars))
message("Somatic Mutations      : ", length(mutation_vars))
message("Copy Number Alterations: ", length(cna_vars))
message("RPPA Proteins          : ", length(rppa_vars))
message("RNA Curated Pathways   : ", length(rna_vars), " (8 signatures)")
message(
  "RNA Data-Driven Score  : ",
  length(rna_datadriven_vars),
  " (1 axis)"
)
message("====================================")


# Build explicit model formulae ------------------------------------------------

model_specs <- list(
  Clinical = clinical_vars,
  Mutations = c(clinical_vars, mutation_vars),
  CNA = c(clinical_vars, cna_vars),
  RPPA = c(clinical_vars, rppa_vars),
  RNA_Path = c(clinical_vars, rna_vars),
  RNA_DataDriven = c(clinical_vars, rna_datadriven_vars),
  Integrated = c(clinical_vars, mutation_vars, cna_vars, rppa_vars, rna_vars)
)


# Create stratified cross-validation folds -------------------------------------
# We set a static seed to ensure reproducible, deterministic data partitioning

set.seed(42)

# Subset to complete cases across all variables to ensure an identical base
all_modelling_vars <- unique(unlist(model_specs))
cv_data <- survival_data |>
  dplyr::select(dplyr::all_of(c(
    "os_months", "os_event", all_modelling_vars
  ))) |>
  tidyr::drop_na()

message(
  "Total complete-case samples available for 5-fold cross-validation: ",
  nrow(cv_data)
)

n_folds <- 5L
# Stratify folds by event status to distribute deaths equally across splits
cv_data <- cv_data |>
  dplyr::group_by(.data$os_event) |>
  dplyr::mutate(fold = sample(rep(1L:n_folds, length.out = dplyr::n()))) |>
  dplyr::ungroup()


# Execute the cross-validation loop --------------------------------------------

# Initialise data frames to hold out-of-fold linear risk predictors
predictions_df <- cv_data |>
  dplyr::select(os_months, os_event, fold)

for (model_name in names(model_specs)) {
  predictions_df[[paste0("lp_", model_name)]] <- NA_real_
}

message("Executing stratified 5-fold cross-validation loop + invariance guards")

for (f in seq_len(n_folds)) {
  train_fold <- cv_data |> dplyr::filter(.data$fold != f)
  test_fold <- cv_data |> dplyr::filter(.data$fold == f)

  for (model_name in names(model_specs)) {
    vars <- model_specs[[model_name]]

    # INVARIANCE GUARD: Filter out any binary variables that lack variation or
    # events inside this specific fold training split
    stable_vars <- c()
    for (v in vars) {
      if (v %in% selected_clinical_features ||
          (is.numeric(train_fold[[v]]) &&
           !all(train_fold[[v]] == train_fold[[v]][1L]))) {
        # Additional validation check for binary markers: ensure there are at
        # least 2 events in each group to prevent quasi-separation.
        # Clinical features are exempted from this check since they are not
        # binary — selected_clinical_features replaces the old hardcoded list.
        if (!v %in% selected_clinical_features &&
            all(range(train_fold[[v]], na.rm = TRUE) == c(0L, 1L))) {
          events_in_mutated <- sum(
            train_fold$os_event[train_fold[[v]] == 1L],
            na.rm = TRUE
          )
          events_in_wildtype <- sum(
            train_fold$os_event[train_fold[[v]] == 0L],
            na.rm = TRUE
          )
          if (events_in_mutated < 2L || events_in_wildtype < 2L) {
            next # Drop feature for this specific fold iteration dynamically
          }
        }
        stable_vars <- c(stable_vars, v)
      }
    }

    if (length(stable_vars) == 0L) {
      next
    }

    # Construct formula string dynamically
    formula_str <- paste(
      "survival::Surv(os_months, os_event) ~",
      paste(stable_vars, collapse = " + ")
    )
    form <- stats::as.formula(formula_str)

    # Fit Cox model with Efron tie handling and suppression of separation notes
    fit_cv <- tryCatch(
      suppressWarnings(survival::coxph(
        form,
        data = train_fold, ties = "efron"
      )),
      error = function(e) {
        NULL
      }
    )

    if (!is.null(fit_cv)) {
      # Predict risk scores onto the independent held-out test fold
      lp <- stats::predict(fit_cv, newdata = test_fold, type = "lp")
      predictions_df[[paste0(
        "lp_",
        model_name
      )]][predictions_df$fold == f] <- lp
    }
  }
}


# Compute out-of-fold cross-validated c-indices --------------------------------

cv_cindex_list <- list()

for (model_name in names(model_specs)) {
  lp_col <- paste0("lp_", model_name)

  # Check if predictions were generated successfully
  if (all(is.na(predictions_df[[lp_col]]))) {
    warning("No valid out-of-fold predictions found for model: ",
      model_name,
      call. = FALSE
    )
    next
  }

  # Evaluate true out-of-fold Concordance via Harrell's method
  c_fit <- survival::concordance(
    survival::Surv(os_months, os_event) ~ predictions_df[[lp_col]],
    data    = predictions_df,
    reverse = TRUE
  )

  # Calculate standard errors to establish rigorous confidence intervals
  c_index <- c_fit$concordance
  c_se <- sqrt(c_fit$var)

  cv_cindex_list[[model_name]] <- tibble::tibble(
    model = model_name,
    cv_concordance = c_index,
    conf_low = c_index - (1.96 * c_se),
    conf_high = c_index + (1.96 * c_se),
    incremental_vs_clinical = if (model_name == "Clinical") {
      0.0
    } else {
      c_index - cv_cindex_list[["Clinical"]]$cv_concordance
    }
  )
}

cv_cindex_results <- dplyr::bind_rows(cv_cindex_list) |>
  dplyr::arrange(dplyr::desc(.data$cv_concordance))

readr::write_csv(
  cv_cindex_results,
  "results/cross_validated_cindex_comparison.csv"
)


# Fit final integrated model on full cohort -----------------------------------
# Trains the ultimate unified model on full data to extract coefficients

message("\n=== Cross-Validated Model Comparison Summary ===")
print(cv_cindex_results)
message("================================================")

message("Fitting final integrated model on full complete-case cohort...")
final_formula <- stats::as.formula(paste(
  "survival::Surv(os_months, os_event) ~",
  paste(model_specs$Integrated, collapse = " + ")
))

final_integrated_fit <- survival::coxph(
  final_formula,
  data = cv_data,
  ties = "efron"
)

final_coefs_df <- broom::tidy(final_integrated_fit,
  exponentiate = TRUE,
  conf.int = TRUE
) |>
  dplyr::transmute(
    feature      = .data$term,
    hazard_ratio = .data$estimate,
    conf_low     = .data$conf.low,
    conf_high    = .data$conf.high,
    p_value      = .data$p.value
  ) |>
  dplyr::arrange(.data$p_value)

readr::write_csv(
  final_coefs_df,
  "results/final_integrated_cox_coefficients.csv"
)

message("Survival modelling complete. Model objects held in memory + exported.")
