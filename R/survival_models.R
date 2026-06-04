# Cross-validated multi-omics survival modelling ------------------------------
#
# Evaluates and compares the out-of-fold predictive performance of six distinct
# survival models using stratified 5-fold cross-validation with an invariance guard.
#
# Requires: integrate_data.R, feature_selection.R, and pathway_scores.R
#           to have been sourced, with survival_data available in global memory.
#
# Usage: this script is intended to be sourced by run_analysis.R.

# Validate inputs -------------------------------------------------------------
check_required_objects(c("survival_data", "selected_rppa_features", "selected_mutation_features"))
required_base_cols <- c("os_months", "os_event", "age", "sex", "stage", "grade")
check_has_columns("survival_data", required_base_cols)

# 1. Feature Space Discovery ---------------------------------------------------
clinical_vars <- c("age", "sex", "stage", "grade")
mutation_vars <- selected_mutation_features
cna_vars      <- names(survival_data)[stringr::str_starts(names(survival_data), "cna_")]
rppa_vars     <- selected_rppa_features
rna_vars      <- names(survival_data)[stringr::str_starts(names(survival_data), "score_") & !stringr::str_ends(names(survival_data), "_mean")]

message("=== Modelling Feature Dimensions ===")
message("Clinical Covariates : ", length(clinical_vars))
message("Somatic Mutations   : ", length(mutation_vars))
message("Copy Number Alterat.: ", length(cna_vars))
message("RPPA Proteins       : ", length(rppa_vars))
message("RNA Pathway Scores  : ", length(rna_vars))
message("====================================")

# 2. Build Explicit Model Formulae ---------------------------------------------
model_specs <- list(
   Clinical   = clinical_vars,
   Mutations  = c(clinical_vars, mutation_vars),
   CNA        = c(clinical_vars, cna_vars),
   RPPA       = c(clinical_vars, rppa_vars),
   RNA_Path   = c(clinical_vars, rna_vars),
   Integrated = c(clinical_vars, mutation_vars, cna_vars, rppa_vars, rna_vars)
)

# 3. Create Stratified Cross-Validation Folds ----------------------------------
set.seed(42)
all_modelling_vars <- unique(unlist(model_specs))
cv_data <- survival_data |>
   dplyr::select(dplyr::all_of(c("os_months", "os_event", all_modelling_vars))) |>
   tidyr::drop_na()

message("Total complete-case samples available for 5-fold cross-validation: ", nrow(cv_data))

n_folds <- 5L
cv_data <- cv_data |>
   dplyr::group_by(.data$os_event) |>
   dplyr::mutate(fold = sample(rep(1:n_folds, length.out = dplyr::n()))) |>
   dplyr::ungroup()

# 4. Execute the Cross-Validation Loop -----------------------------------------
predictions_df <- cv_data |> dplyr::select(os_months, os_event, fold)
for (model_name in names(model_specs)) {
   predictions_df[[paste0("lp_", model_name)]] <- NA_real_
}

message("Executing stratified 5-fold cross-validation loop with invariance guards...")

for (f in seq_len(n_folds)) {
   train_fold <- cv_data |> dplyr::filter(.data$fold != f)
   test_fold  <- cv_data |> dplyr::filter(.data$fold == f)
   
   for (model_name in names(model_specs)) {
      vars <- model_specs[[model_name]]
      
      # INVARIANCE GUARD: Filter out any binary variables that lack variation in this specific training fold
      stable_vars <- c()
      for (v in vars) {
         if (v %in% c("age", "sex", "stage", "grade") || is.numeric(train_fold[[v]]) && !all(train_fold[[v]] == train_fold[[v]][1])) {
            # Keep continuous variables or binary features with active variation
            
            # Additional safety check for binary variables: ensure there's at least 2 events in the minor group
            if (!v %in% c("age", "sex", "stage", "grade") && all(range(train_fold[[v]]) == c(0, 1))) {
               events_in_mutated <- sum(train_fold$os_event[train_fold[[v]] == 1L])
               events_in_wildtype <- sum(train_fold$os_event[train_fold[[v]] == 0L])
               if (events_in_mutated < 2L || events_in_wildtype < 2L) {
                  next # Drop feature for this fold to eliminate complete separation warnings
               }
            }
            stable_vars <- c(stable_vars, v)
         }
      }
      
      if (length(stable_vars) == 0) next
      
      formula_str <- paste("survival::Surv(os_months, os_event) ~", paste(stable_vars, collapse = " + "))
      form <- stats::as.formula(formula_str)
      
      fit_cv <- tryCatch(
         survival::coxph(form, data = train_fold, ties = "efron"),
         error = function(e) NULL
      )
      
      if (!is.null(fit_cv)) {
         lp <- stats::predict(fit_cv, newdata = test_fold, type = "lp")
         predictions_df[[paste0("lp_", model_name)]][predictions_df$fold == f] <- lp
      }
   }
}

# 5. Compute Out-of-Fold Cross-Validated C-Indices -----------------------------
cv_cindex_list <- list()
for (model_name in names(model_specs)) {
   lp_col <- paste0("lp_", model_name)
   
   c_fit <- survival::concordance(
      survival::Surv(os_months, os_event) ~ predictions_df[[lp_col]],
      data    = predictions_df,
      reverse = TRUE
   )
   
   c_index <- c_fit$concordance
   c_se    <- sqrt(c_fit$var)
   
   cv_cindex_list[[model_name]] <- tibble::tibble(
      model                  = model_name,
      cv_concordance         = c_index,
      conf_low               = c_index - (1.96 * c_se),
      conf_high              = c_index + (1.96 * c_se),
      incremental_vs_clinical = if (model_name == "Clinical") 0.0 else c_index - cv_cindex_list[["Clinical"]]$cv_concordance
   )
}

cv_cindex_results <- dplyr::bind_rows(cv_cindex_list) |>
   dplyr::arrange(dplyr::desc(.data$cv_concordance))

readr::write_csv(cv_cindex_results, "results/cross_validated_cindex_comparison.csv")

message("\n=== Cross-Validated Model Comparison Summary ===")
print(cv_cindex_results)
message("================================================")

# 6. Fit Final Model on Full Cohort -------------------------------------------
message("Fitting final integrated model on full complete-case cohort...")
final_formula <- stats::as.formula(
   paste("survival::Surv(os_months, os_event) ~", paste(model_specs$Integrated, collapse = " + "))
)
final_integrated_fit <- survival::coxph(final_formula, data = cv_data, ties = "efron")

final_coefs_df <- broom::tidy(final_integrated_fit, exponentiate = TRUE, conf.int = TRUE) |>
   dplyr::transmute(
      feature      = .data$term,
      hazard_ratio = .data$estimate,
      conf_low     = .data$conf.low,
      conf_high    = .data$conf.high,
      p_value      = .data$p.value
   ) |>
   dplyr::arrange(.data$p_value)

readr::write_csv(final_coefs_df, "results/final_integrated_cox_coefficients.csv")
message("Survival modelling complete. Model objects held in memory and exported to /results.")