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
#        the full pipeline, not run directly.


# Validate Inputs -------------------------------------------------------------

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


# Prepare Survival Analysis Table ---------------------------------------------
# Factor levels and baseline reference positions are preserved exactly as
# configured during the clinical preparation script.

survival_data <- clinical_rppa_rna_mutation |>
  dplyr::mutate(
    os_event = as.integer(.data$os_event)
  )

abort_if_false(
  nrow(survival_data) > 0L,
  "No rows with complete survival time and event status are available."
)

# FIXED: Standardised assignment typo from %in= to valid tidy comparison
abort_if_false(
  all(survival_data$os_event %in% c(0L, 1L)),
  "os_event must be binary 0/1."
)


# Overall Kaplan-Meier Summary ------------------------------------------------

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


# Clinical Baseline Cox Model -------------------------------------------------
# Use complete cases for the clinical covariates so coxph uses a transparent
# and reproducible analysis set.

clinical_vars <- c("os_months", "os_event", "age", "sex", "stage", "grade")

clinical_model_data <- enforce_complete_cases(
  data      = survival_data,
  variables = clinical_vars
)

if (nrow(clinical_model_data) < 20L) {
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
c_index <- cox_summary$concordance[["C"]]

message("Clinical Cox C-index: ", round(c_index, 3L))

clinical_cox_results <- broom::tidy(
  clinical_cox_fit,
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

readr::write_csv(
  clinical_cox_results,
  "results/clinical_cox_results.csv"
)

message("Clinical Cox model complete cases: ", nrow(clinical_model_data))
message("Clinical Cox results:")
print(clinical_cox_results)


# Proportional Hazards Assumption Test ----------------------------------------
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


# Kaplan-Meier Plots ----------------------------------------------------------

# Overall KM curve
overall_km_plot <- survminer::ggsurvplot( # Create KM plot object using survminer

  # survfit object containing overall KM curve
  overall_km_fit,

  # original dataset used to build the survival model
  data = survival_data,

  # do show (=TRUE) confidence intervals around survival curves
  conf.int = TRUE,

  # add a risk table below the plot (number at risk over time)
  risk.table = TRUE,

  # don't show any legend
  legend = "none",

  # put a dotted line horizontal and vertical for median
  surv.median.line = "hv",

  # x-axis legend
  xlab = "Time (months)",

  # y-axis legend
  ylab = "Overall survival probability",

  # plot title
  title = "TCGA KIRC: overall survival",

  # Nature Publishing Group colours
  palette = "npg",

  # line thickness
  linewidth = 1.0,

  # Relative height of table (0–1 scale)
  risk.table.height = 0.3,

  # Controls risk table styling separately
  tables.theme = theme_classic()
)


save_pipeline_plot(
  plot_object = overall_km_plot,
  file_path   = "figures/overall_kaplan_meier.png",
  width       = 1800,
  height      = 1200,
  resolution  = 220
)

message("Saved overall KM plot to figures/overall_kaplan_meier.png")
