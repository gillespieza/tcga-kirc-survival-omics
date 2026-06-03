# Quick survival checks -------------------------------------------------------
#
# Creates the first report-ready survival outputs from the integrated TCGA KIRC
# dataset: cohort summary, overall and stage-stratified Kaplan-Meier curves,
# a clinical baseline Cox model using age, sex, pathological stage, and grade,
# and a proportional hazards assumption check.
#
# Requires: integrate_data.R to have been sourced so that 
#   clinical_rppa_rna_mutation is available in the global environment.
#   broom must be installed (it is installed alongside tidyverse but not
#   attached by library("tidyverse")).
#
# Produces:
#   survival_data            - Integrated table with complete survival data,
#                              stage and grade re-levelled for Cox modelling
#   overall_surv_obj         - Surv object for overall survival
#   overall_km_fit           - Kaplan-Meier fit for all samples
#   stage_km_fit             - Kaplan-Meier fit stratified by AJCC stage
#   overall_survival_summary - One-row cohort survival summary
#   clinical_cox_fit         - Cox model using clinical covariates
#   clinical_cox_results     - Tidy hazard-ratio table for clinical model
#   cox_zph                  - Schoenfeld residual test for proportional
#                              hazards assumption
#
# Outputs:
#   results/overall_survival_summary.csv
#   results/clinical_cox_results.csv
#   figures/overall_kaplan_meier.png
#   figures/stage_kaplan_meier.png
#
# Note: this script is intended to be sourced by run_analysis.R as part of
#       the full pipeline, not run directly.


# Validate inputs -------------------------------------------------------------

check_required_objects("clinical_rppa_rna_mutation")

required_survival_cols <- c(
   "sample_id",
   "patient_id",
   "os_months",
   "os_event",
   "age",
   "sex",
   "stage",
   "grade"
)

check_has_columns(
   "clinical_rppa_rna_mutation",
   required_survival_cols
)

# Prepare survival analysis table ---------------------------------------------

survival_data <- clinical_rppa_rna_mutation |>
   dplyr::filter(
      !is.na(.data$os_months),
      !is.na(.data$os_event)
   ) |>
   dplyr::mutate(
      os_event = as.integer(.data$os_event),
      sex      = droplevels(factor(.data$sex)),
      # Treat ordered clinical categories as reference-coded factors in Cox
      # models. This gives interpretable contrasts such as STAGE IV vs STAGE I
      # instead of polynomial terms like stage.L and stage.Q.
      stage    = stats::relevel(
         droplevels(factor(.data$stage, ordered = FALSE)),
         ref = "STAGE I"
      ),
      grade    = stats::relevel(
         droplevels(factor(.data$grade, ordered = FALSE)),
         ref = "G1"
      )
   )

abort_if_false(
   nrow(survival_data) > 0,
   "No rows with complete survival time and event status are available."
)

abort_if_false(
   all(survival_data$os_event %in% c(0L, 1L)),
   "os_event must be binary 0/1."
)


# Overall Kaplan-Meier summary ------------------------------------------------

overall_surv_obj <- survival::Surv(
   time  = survival_data$os_months,
   event = survival_data$os_event
)

overall_km_fit <- survival::survfit(
   overall_surv_obj ~ 1,
   data = survival_data
)

# Extract median survival using quantile() rather than indexing $table
# directly, which can silently return NaN when the KM curve does not cross 0.5
# (i.e. fewer than half of patients have experienced the event).
median_survival_months <- tryCatch(
   unname(stats::quantile(overall_km_fit, probs = 0.5)$quantile),
   error = function(e) NA_real_
)

overall_survival_summary <- survival_data |>
   dplyr::summarise(
      samples                 = dplyr::n(),
      patients                = dplyr::n_distinct(.data$patient_id),
      events                  = sum(.data$os_event),
      censored                = dplyr::n() - sum(.data$os_event),
      median_follow_up_months = stats::median(.data$os_months, na.rm = TRUE),
      median_survival_months  = median_survival_months
   )

readr::write_csv(
   overall_survival_summary,
   "results/overall_survival_summary.csv"
)

message("Overall survival summary:")
print(overall_survival_summary)


# Clinical baseline Cox model -------------------------------------------------
# Use complete cases for the clinical covariates so coxph uses a transparent
# and reproducible analysis set.

clinical_model_data <- survival_data |>
   dplyr::select(
      .data$os_months,
      .data$os_event,
      .data$age,
      .data$sex,
      .data$stage,
      .data$grade
   ) |>
   tidyr::drop_na()

n_dropped_cc <- nrow(survival_data) - nrow(clinical_model_data)

if (n_dropped_cc > 0) {
   message(
      n_dropped_cc,
      " row(s) dropped for complete case Cox analysis."
   )
}

if (nrow(clinical_model_data) < 20) {
   warning(
      "Clinical Cox model has fewer than 20 complete cases.",
      call. = FALSE
   )
}

clinical_cox_fit <- survival::coxph(
   survival::Surv(os_months, os_event) ~ age + sex + stage + grade,
   data = clinical_model_data,
   ties = "efron"
)

cox_summary <- summary(clinical_cox_fit)
c_index     <- cox_summary$concordance[["C"]]

message("Clinical Cox C-index: ", round(c_index, 3))

clinical_cox_results <- broom::tidy(
   clinical_cox_fit,
   conf.int     = TRUE,
   exponentiate = TRUE
) |>
   dplyr::transmute(
      term,
      hazard_ratio = .data$estimate,
      conf_low     = .data$conf.low,
      conf_high    = .data$conf.high,
      p_value      = .data$p.value
   ) |>
   dplyr::arrange(.data$p_value)

readr::write_csv(
   clinical_cox_results,
   "results/clinical_cox_results.csv"
)

message("Clinical Cox model complete cases: ", nrow(clinical_model_data))
message("Clinical Cox results:")
print(clinical_cox_results)


# Proportional hazards assumption test ----------------------------------------
# Schoenfeld residuals test. A significant p-value for a covariate indicates
# its hazard ratio changes over time, violating the PH assumption. A
# significant global p-value suggests the model as a whole should be revisited
# — for example by adding time-varying coefficients or stratifying on the
# offending covariate.

cox_zph <- survival::cox.zph(clinical_cox_fit)

message("Proportional hazards test (Schoenfeld residuals):")
print(cox_zph)

if (any(cox_zph$table[, "p"] < 0.05)) {
   warning(
      paste(
         "One or more covariates show evidence of non-proportional hazards.",
         "Consider time-varying coefficients or a stratified Cox model."
      ),
      call. = FALSE
   )
}


# Kaplan-Meier plots ----------------------------------------------------------

# Overall KM curve
overall_km_plot <- survminer::ggsurvplot(
   overall_km_fit,
   data       = survival_data,
   risk.table = TRUE,
   legend     = "none",
   xlab       = "Time (months)",
   ylab       = "Overall survival probability",
   title      = "TCGA KIRC: overall survival"
)

grDevices::png(
   filename = "figures/overall_kaplan_meier.png",
   width    = 1800,
   height   = 1600,
   res      = 220
)
print(overall_km_plot)
grDevices::dev.off()

message("Saved overall KM plot to figures/overall_kaplan_meier.png")

# Stage-stratified KM curve
# Log-rank test p-value and method label are included to indicate whether
# stage groups differ significantly in survival.
stage_km_fit <- survival::survfit(
   overall_surv_obj ~ stage,
   data = survival_data
)

stage_km_plot <- survminer::ggsurvplot(
   stage_km_fit,
   data         = survival_data,
   risk.table   = TRUE,
   pval         = TRUE,
   pval.method  = TRUE,
   legend.title = "Stage",
   legend.labs  = levels(survival_data$stage),
   xlab         = "Time (months)",
   ylab         = "Overall survival probability",
   title        = "TCGA KIRC: overall survival by AJCC stage"
)

grDevices::png(
   filename = "figures/stage_kaplan_meier.png",
   width    = 1800,
   height   = 1800,
   res      = 220
)
print(stage_km_plot)
grDevices::dev.off()

message("Saved stage-stratified KM plot to figures/stage_kaplan_meier.png")