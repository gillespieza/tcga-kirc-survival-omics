# Compare survival models -----------------------------------------------------
#
# Fits and compares multiple multivariable Cox proportional hazards models to
# evaluate the incremental prognostic value of adding different omics layers to
# standard clinical risk factors:
#   - Model 1: Clinical baseline (age + sex + stage + grade)
#   - Model 2: Clinical + Top LASSO-selected RPPA features
#   - Model 3: Clinical + RNA-seq pathway scores (PC1 scores)
#   - Model 4: Clinical + Curated binary mutation features
#   - Model 5: Clinical + Copy Number Alterations (CNA)
#   - Model 6: Full Multi-omics Integration (Clinical + RPPA + RNA + Mut + CNA)
#
# Computes concordance indices (C-index) with Wald 95 % confidence intervals,
# compares nested models with Likelihood Ratio Tests (LRT), and calculates
# Akaike Information Criterion (AIC) for non-nested comparison.
#
# Requires: pathway_scores.R and feature_selection.R to have been sourced 
#           so that survival_data and selected features are available.
#
# Produces:
#   model_fits           - Named list containing the 6 coxph model objects.
#   model_comparison_df  - Combined statistical summary table of all models.
#   integrated_hr_df     - Extracted hazard ratios from the full integrated model.
#
# Outputs:
#   results/model_comparison_metrics.csv
#   results/integrated_model_hazard_ratios.csv
#
# Note: this script is intended to be sourced by run_analysis.R.


# Validate inputs -------------------------------------------------------------

check_required_objects(c(
   "survival_data",
   "selected_rppa_features",
   "selected_mutation_features"
))


# Build model formulae dynamically --------------------------------------------

clinical_vars <- c("age", "sex", "stage", "grade")

# Define the feature sets explicitly from the previous pipeline phases
rppa_vars <- selected_rppa_features

rna_vars <- c("score_neutrophil_deg", "score_ecm_deg", "score_ptk")

mutation_vars <- selected_mutation_features

# Dynamically identify which CNA columns successfully passed the prevalence screen
cna_vars <- names(survival_data)[stringr::str_starts(names(survival_data), "cna_")]

# Combine variable groups into specific, incremental model definitions
formula_specs <- list(
   clinical   = clinical_vars,
   rppa       = c(clinical_vars, rppa_vars),
   rna        = c(clinical_vars, rna_vars),
   mutation   = c(clinical_vars, mutation_vars),
   cna        = c(clinical_vars, cna_vars),
   integrated = c(clinical_vars, rppa_vars, rna_vars, mutation_vars, cna_vars)
)

# Build formula objects cleanly by compiling string components and converting
# them explicitly with as.formula to avoid syntax breakdowns.
model_formulae <- purrr::map(
   formula_specs,
   function(vars) {
      formula_str <- paste0(
         "survival::Surv(os_months, os_event) ~ ",
         paste(vars, collapse = " + ")
      )
      stats::as.formula(formula_str)
   }
)


# Complete-case filter across all variables -----------------------------------
# To make Likelihood Ratio Tests and AIC metrics strictly comparable, all 
# models must be evaluated on identical samples.

all_modelling_vars <- c("os_months", "os_event", formula_specs$integrated)

modelling_data <- enforce_complete_cases(
   data      = survival_data,
   variables = all_modelling_vars
)

abort_if_false(
   nrow(modelling_data) >= 30,
   paste(
      "Insufficient data for multivariable modeling after complete-case filter.",
      "Required: >= 30 samples. Available: ", nrow(modelling_data)
   )
)

message(
   "Unified complete-case analysis cohort: ",
   nrow(modelling_data), " samples, ",
   sum(modelling_data$os_event), " events."
)


# Fit models ------------------------------------------------------------------

fit_cox_model <- function(formula_obj, model_name) {
   fit <- tryCatch(
      survival::coxph(
         formula_obj,
         data = modelling_data,
         ties = "efron"
      ),
      error = function(e) {
         stop(
            "Cox model fitting failed for [", model_name, "]: ", e$message,
            call. = FALSE
         )
      }
   )
   fit
}

model_fits <- purrr::imap(model_formulae, fit_cox_model)


# Extract summary metrics -----------------------------------------------------

extract_metrics <- function(fit, name) {
   sum_fit <- summary(fit)
   
   # Wald test for overall model significance
   wald_p <- sum_fit$waldtest[["pvalue"]]
   
   # Concordance index and standard error extracted safely by numeric position
   # to handle variations in survival package label names.
   c_index <- unname(sum_fit$concordance[1])
   c_se <- unname(sum_fit$concordance[2])
   
   # Compute 95 % Wald confidence limits for C-index
   c_low <- max(0, c_index - (1.96 * c_se))
   c_high <- min(1, c_index + (1.96 * c_se))
   
   tibble::tibble(
      model        = name,
      n_coef       = length(fit$coefficients),
      log_lik      = fit$loglik[2],
      aic          = stats::AIC(fit),
      c_index      = c_index,
      c_conf_low   = c_low,
      c_conf_high  = c_high,
      wald_p_value = wald_p
   )
}

metrics_df <- purrr::imap_dfr(model_fits, extract_metrics)


# Statistical nested model comparisons (LRT) ----------------------------------
# Likelihood Ratio Test evaluates whether an omics layer provides a
# statistically significant improvement over the clinical baseline.

compute_lrt_vs_clinical <- function(target_fit) {
   lrt <- stats::anova(model_fits$clinical, target_fit)
   # Extract p-value from the second row using the exact column mapping P(>|Chi|)
   lrt[["P(>|Chi|)"]][2]
}

lrt_p_values <- c(
   clinical   = NA_real_,
   rppa       = compute_lrt_vs_clinical(model_fits$rppa),
   rna        = compute_lrt_vs_clinical(model_fits$rna),
   mutation   = compute_lrt_vs_clinical(model_fits$mutation),
   cna        = compute_lrt_vs_clinical(model_fits$cna),
   integrated = compute_lrt_vs_clinical(model_fits$integrated)
)

model_comparison_df <- metrics_df |>
   dplyr::mutate(
      lrt_p_value_vs_clinical = lrt_p_values[match(
         .data$model,
         names(lrt_p_values)
      )]
   ) |>
   dplyr::arrange(.data$aic)


# Extract and save hazard ratios for plotting ---------------------------------
# Tidy and format the full integrated model coefficients to feed the forest plot.

integrated_hr_df <- broom::tidy(
   model_fits$integrated,
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


# Print summary ---------------------------------------------------------------

message("\nModel comparison summary (ordered by AIC):")
print(model_comparison_df)

best_model <- model_comparison_df$model[1]
message(
   "\nModel comparison complete. Best performing model by AIC: [",
   best_model, "]"
)