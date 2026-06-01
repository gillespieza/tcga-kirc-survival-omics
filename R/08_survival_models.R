# Survival model comparison ---------------------------------------------------
#
# Fits and compares clinical-only, mutation-only, RPPA-only, integrated Cox, and
# optional penalised Cox models for TCGA KIRC overall survival. The goal is a
# compact assignment-ready comparison of whether omics features add useful
# survival signal beyond clinical covariates.
#
# Requires: 07_feature_selection.R to have been sourced.
#
# Produces:
#   model_comparison_results - C-index and sample/event counts for each model
#   fitted_survival_models   - Named list of fitted Cox/glmnet model objects
#   integrated_cox_results   - Hazard-ratio table for the integrated Cox model
#
# Outputs:
#   results/model_comparison_results.csv
#   results/integrated_cox_results.csv
#   results/lasso_selected_features.csv, when LASSO is fitted


# Validate inputs -------------------------------------------------------------

required_objects <- c("survival_data", "selected_rppa_features", "selected_mutation_features")
missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
   stop(
      "Missing required object(s): ", paste(missing_objects, collapse = ", "),
      ". Source 07_feature_selection.R first."
   )
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

clinical_model_terms <- c("age", "sex", "stage", "grade")
selected_rppa_features <- intersect(selected_rppa_features, names(survival_data))
selected_mutation_features <- intersect(selected_mutation_features, names(survival_data))

if (length(selected_rppa_features) == 0) {
   stop("No selected RPPA features are present in survival_data.")
}


# Helpers ---------------------------------------------------------------------

build_model_data <- function(feature_cols) {
   survival_data %>%
      dplyr::select(os_months, os_event, dplyr::all_of(feature_cols)) %>%
      tidyr::drop_na()
}

cox_formula_from_terms <- function(feature_cols) {
   stats::as.formula(
      paste(
         "survival::Surv(os_months, os_event) ~",
         paste(sprintf("`%s`", feature_cols), collapse = " + ")
      )
   )
}

fit_cox_model <- function(feature_cols, model_name) {
   model_data <- build_model_data(feature_cols)
   
   if (nrow(model_data) < 30 || sum(model_data$os_event) < 5) {
      warning(model_name, " skipped because it has too few complete cases/events.")
      return(NULL)
   }
   
   fit <- survival::coxph(
      cox_formula_from_terms(feature_cols),
      data = model_data,
      ties = "efron",
      x = TRUE
   )
   
   fit_summary <- summary(fit)
   tibble::tibble(
      model       = model_name,
      n           = nrow(model_data),
      events      = sum(model_data$os_event),
      predictors  = length(feature_cols),
      c_index     = unname(fit_summary$concordance[1]),
      c_index_se  = unname(fit_summary$concordance[2]),
      partial_aic = stats::extractAIC(fit)[2]
   ) %>%
      list(fit = fit, data = model_data, summary = .)
}

extract_cox_results <- function(fit) {
   fit_summary <- summary(fit)
   as.data.frame(fit_summary$coefficients) %>%
      tibble::rownames_to_column("term") %>%
      dplyr::as_tibble() %>%
      dplyr::transmute(
         term,
         hazard_ratio = exp(.data$coef),
         conf_low     = exp(.data$coef - 1.96 * .data$`se(coef)`),
         conf_high    = exp(.data$coef + 1.96 * .data$`se(coef)`),
         p_value      = .data$`Pr(>|z|)`
      ) %>%
      dplyr::arrange(.data$p_value)
}


# Fit interpretable Cox models ------------------------------------------------

model_terms <- list(
   clinical_only = clinical_model_terms,
   mutation_only = selected_mutation_features,
   rppa_only     = selected_rppa_features,
   integrated    = c(clinical_model_terms, selected_rppa_features, selected_mutation_features)
)

# Mutation-only is skipped automatically if no mutation features are available.
model_terms <- model_terms[lengths(model_terms) > 0]

fitted_survival_models <- purrr::imap(model_terms, fit_cox_model)
fitted_survival_models <- fitted_survival_models[!vapply(fitted_survival_models, is.null, logical(1))]

if (length(fitted_survival_models) == 0) {
   stop("No Cox models could be fitted.")
}

model_comparison_results <- purrr::map_dfr(
   fitted_survival_models,
   ~ .x$summary
) %>%
   dplyr::arrange(dplyr::desc(.data$c_index))

readr::write_csv(model_comparison_results, "results/model_comparison_results.csv")

if ("integrated" %in% names(fitted_survival_models)) {
   integrated_cox_results <- extract_cox_results(fitted_survival_models$integrated$fit)
   readr::write_csv(integrated_cox_results, "results/integrated_cox_results.csv")
} else {
   integrated_cox_results <- tibble::tibble()
}

message("Model comparison results:")
print(model_comparison_results)


# Optional penalised Cox model ------------------------------------------------
# LASSO is useful when the integrated feature count is large relative to events.
# Here it is fitted when enough events exist for cross-validation.

lasso_candidate_features <- c(
   clinical_model_terms,
   selected_rppa_features,
   selected_mutation_features
)
lasso_data <- build_model_data(lasso_candidate_features)

lasso_selected_features <- tibble::tibble()

if (nrow(lasso_data) >= 50 && sum(lasso_data$os_event) >= 20) {
   x_lasso <- stats::model.matrix(
      cox_formula_from_terms(lasso_candidate_features),
      data = lasso_data
   )[, -1, drop = FALSE]
   y_lasso <- survival::Surv(lasso_data$os_months, lasso_data$os_event)
   
   set.seed(42)
   lasso_cv_fit <- glmnet::cv.glmnet(
      x = x_lasso,
      y = y_lasso,
      family = "cox",
      alpha = 1,
      nfolds = min(10, sum(lasso_data$os_event)),
      standardize = TRUE
   )
   
   lasso_coef <- as.matrix(stats::coef(lasso_cv_fit, s = "lambda.1se"))
   lasso_selected_features <- tibble::tibble(
      feature = rownames(lasso_coef),
      coefficient = as.numeric(lasso_coef[, 1])
   ) %>%
      dplyr::filter(.data$coefficient != 0) %>%
      dplyr::arrange(dplyr::desc(abs(.data$coefficient)))
   
   lasso_lp <- as.numeric(stats::predict(lasso_cv_fit, newx = x_lasso, s = "lambda.1se", type = "link"))
   lasso_concordance <- survival::concordance(y_lasso ~ I(-lasso_lp))
   
   model_comparison_results <- model_comparison_results %>%
      dplyr::bind_rows(
         tibble::tibble(
            model       = "lasso_integrated",
            n           = nrow(lasso_data),
            events      = sum(lasso_data$os_event),
            predictors  = ncol(x_lasso),
            c_index     = unname(lasso_concordance$concordance),
            c_index_se  = sqrt(unname(lasso_concordance$var)),
            partial_aic = NA_real_
         )
      ) %>%
      dplyr::arrange(dplyr::desc(.data$c_index))
   
   fitted_survival_models$lasso_integrated <- list(
      fit = lasso_cv_fit,
      data = lasso_data,
      selected_features = lasso_selected_features
   )
   
   readr::write_csv(model_comparison_results, "results/model_comparison_results.csv")
   readr::write_csv(lasso_selected_features, "results/lasso_selected_features.csv")
   
   message("LASSO Cox model fitted with lambda.1se.")
   message("LASSO-selected features:")
   print(lasso_selected_features)
} else {
   readr::write_csv(lasso_selected_features, "results/lasso_selected_features.csv")
   message("LASSO Cox model skipped: fewer than 50 complete cases or 20 events.")
}
