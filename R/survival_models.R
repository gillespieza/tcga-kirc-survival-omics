# Compare survival models via cross-validation -------------------------------
#
# Fits and compares six Cox proportional hazards models using a strict 5-fold
# cross-validation loop to completely eliminate selection bias and overfitting.
# Feature selection (LASSO) is embedded entirely within the training folds.
#
# Models evaluated:
#   - Model 1: Clinical baseline (age + sex + stage + grade)
#   - Model 2: Clinical + LASSO-selected RPPA features
#   - Model 3: Clinical + RNA-seq pathway scores (PC1 scores)
#   - Model 4: Clinical + Curated binary mutation features
#   - Model 5: Clinical + Copy Number Alterations (CNA)
#   - Model 6: Full Multi-omics Integration (Clinical + RPPA + RNA + Mut + CNA)
#
# Outputs:
#   results/model_comparison_metrics.csv
#   results/integrated_model_hazard_ratios.csv
#
# Note: this script is intended to be sourced by run_analysis.R.


# Validate inputs -------------------------------------------------------------

check_required_objects(c(
   "survival_data",
   "mutation_feature_cols",
   "rppa_feature_cols"
))


# Define feature groups --------------------------------------------------------

clinical_vars <- c("age", "sex", "stage", "grade")

rna_vars <- c("score_neutrophil_deg", "score_ecm_deg", "score_ptk")

mutation_vars <- mutation_feature_cols

cna_vars <- names(survival_data)[stringr::str_starts(names(survival_data), "cna_")]


# Establish unified complete-case analysis cohort ------------------------------

all_possible_vars <- c(
   "os_months",
   "os_event",
   clinical_vars,
   rppa_feature_cols,
   rna_vars,
   mutation_vars,
   cna_vars
)

modelling_data <- enforce_complete_cases(
   data      = survival_data,
   variables = all_possible_vars
)

abort_if_false(
   nrow(modelling_data) >= 50,
   "Insufficient samples for robust 5-fold cross-validation."
)


# Construct design matrix for regularised models -------------------------------
# Converts factors to dummy variables globally to maintain matching row profiles

x_master <- stats::model.matrix(
   stats::as.formula(paste("~", paste(c(clinical_vars, rppa_feature_cols, rna_vars, mutation_vars, cna_vars), collapse = " + "))),
   data = modelling_data
)[, -1]

y_master <- survival::Surv(
   time  = modelling_data$os_months,
   event = modelling_data$os_event
)

# Identify the actual dummy-expanded clinical columns present in the model matrix
clinical_cols_in_x <- setdiff(
   colnames(x_master),
   c(rppa_feature_cols, rna_vars, mutation_vars, cna_vars)
)


# Setup cross-validation loops -------------------------------------------------

set.seed(42)

n_folds <- 5

folds <- sample(rep(seq_len(n_folds), length.out = nrow(modelling_data)))

# Data frame to collect unbiased, independent out-of-fold risk predictions
oof_predictions <- tibble::tibble(
   clinical   = rep(NA_real_, nrow(modelling_data)),
   rppa       = rep(NA_real_, nrow(modelling_data)),
   rna        = rep(NA_real_, nrow(modelling_data)),
   mutation   = rep(NA_real_, nrow(modelling_data)),
   cna        = rep(NA_real_, nrow(modelling_data)),
   integrated = rep(NA_real_, nrow(modelling_data))
)

message("Executing 5-fold cross-validation loop with embedded feature selection...")


# Execute cross-validation ----------------------------------------------------

for (f in seq_len(n_folds)) {
   message("  Processing fold ", f, " of ", n_folds, "...")
   
   train_idx <- which(folds != f)
   test_idx  <- which(folds == f)
   
   train_df <- modelling_data[train_idx, ]
   test_df  <- modelling_data[test_idx, ]
   
   # Model 1: Clinical Baseline
   fit_clin <- survival::coxph(
      survival::Surv(os_months, os_event) ~ age + sex + stage + grade,
      data = train_df,
      ties = "efron"
   )
   oof_predictions$clinical[test_idx] <- stats::predict(fit_clin, newdata = test_df, type = "lp")
   
   # Model 2: Clinical + Proteomics (Embedded LASSO selection)
   rppa_model_cols <- c(clinical_cols_in_x, rppa_feature_cols)
   x_train_rppa    <- x_master[train_idx, rppa_model_cols, drop = FALSE]
   x_test_rppa     <- x_master[test_idx, rppa_model_cols, drop = FALSE]
   
   p_fac_rppa      <- rep(1, ncol(x_train_rppa))
   p_fac_rppa[seq_along(clinical_cols_in_x)] <- 0.001
   
   cv_rppa <- glmnet::cv.glmnet(
      x              = x_train_rppa,
      y              = y_master[train_idx],
      family         = "cox",
      penalty.factor = p_fac_rppa,
      cox.ties       = "efron"
   )
   oof_predictions$rppa[test_idx] <- as.numeric(stats::predict(cv_rppa, newx = x_test_rppa, s = "lambda.1se", type = "link"))
   
   # Model 3: Clinical + RNA-seq Pathways
   fit_rna <- survival::coxph(
      survival::Surv(os_months, os_event) ~ age + sex + stage + grade + score_neutrophil_deg + score_ecm_deg + score_ptk,
      data = train_df,
      ties = "efron"
   )
   oof_predictions$rna[test_idx] <- stats::predict(fit_rna, newdata = test_df, type = "lp")
   
   # Model 4: Clinical + Driver Mutations
   mut_formula <- stats::as.formula(paste("survival::Surv(os_months, os_event) ~ age + sex + stage + grade +", paste(mutation_vars, collapse = " + ")))
   fit_mut <- survival::coxph(mut_formula, data = train_df, ties = "efron")
   oof_predictions$mutation[test_idx] <- stats::predict(fit_mut, newdata = test_df, type = "lp")
   
   # Model 5: Clinical + CNA Alterations
   cna_formula <- stats::as.formula(paste("survival::Surv(os_months, os_event) ~ age + sex + stage + grade +", paste(cna_vars, collapse = " + ")))
   fit_cna <- survival::coxph(cna_formula, data = train_df, ties = "efron")
   oof_predictions$cna[test_idx] <- stats::predict(fit_cna, newdata = test_df, type = "lp")
   
   # Model 6: Full Multi-omics Integration (Embedded LASSO selection across all layers)
   p_fac_int <- rep(1, ncol(x_master))
   p_fac_int[seq_along(clinical_cols_in_x)] <- 0.001
   
   cv_int <- glmnet::cv.glmnet(
      x              = x_master[train_idx, , drop = FALSE],
      y              = y_master[train_idx],
      family         = "cox",
      penalty.factor = p_fac_int,
      cox.ties       = "efron"
   )
   oof_predictions$integrated[test_idx] <- as.numeric(stats::predict(cv_int, newx = x_master[test_idx, , drop = FALSE], s = "lambda.1se", type = "link"))
}


# Calculate Unbiased Cross-Validated Performance Metrics ----------------------

calculate_cv_cindex <- function(predictions, model_name) {
   conc_obj <- survival::concordance(y_master ~ predictions)
   
   c_index <- conc_obj$concordance
   c_se    <- sqrt(conc_obj$var)
   
   tibble::tibble(
      model       = model_name,
      c_index     = c_index,
      c_conf_low  = max(0, c_index - (1.96 * c_se)),
      c_conf_high = min(1, c_index + (1.96 * c_se))
   )
}

model_names <- c("clinical", "rppa", "rna", "mutation", "cna", "integrated")

metrics_list <- purrr::map2(oof_predictions, model_names, calculate_cv_cindex)

model_comparison_df <- dplyr::bind_rows(metrics_list)


# Fit Full Cohort Master Model for Hazard Ratio Reporting ---------------------

message("Fitting final integrated model on full cohort for hazard ratio extraction...")

final_integrated_fit <- survival::coxph(
   stats::as.formula(paste("survival::Surv(os_months, os_event) ~ age + sex + stage + grade +", 
                           paste(c(selected_rppa_features, rna_vars, mutation_vars, cna_vars), collapse = " + "))),
   data = modelling_data,
   ties = "efron"
)

integrated_hr_df <- broom::tidy(
   final_integrated_fit,
   conf.int     = TRUE,
   exponentiate = TRUE
) |>
   dplyr::transmute(
      term         = .data$term,
      hazard_ratio = .data$estimate,
      conf_low     = .data$conf.low,
      conf_high    = .data$conf.high,
      p_value      = .data$p.value
   ) |>
   dplyr::arrange(.data$p_value)


# Save results ----------------------------------------------------------------

readr::write_csv(
   model_comparison_df,
   "results/model_comparison_metrics.csv"
)

readr::write_csv(
   integrated_hr_df,
   "results/integrated_model_hazard_ratios.csv"
)

message("\nCross-validated model comparison summary (leak-free C-indices):")
print(model_comparison_df)