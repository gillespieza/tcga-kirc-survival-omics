# Quick survival checks --------------------------------------------------------
#
# Creates the first report-ready survival outputs from the integrated TCGA KIRC
# dataset: cohort summary, overall Kaplan-Meier curve, and a clinical baseline
# Cox model using age, sex, pathological stage, and grade.
#
# Requires: 05_integrate_data.R to have been sourced so that
#   clinical_rppa_mutation is available.
#
# Produces:
#   survival_data            - Integrated table with complete survival data
#   overall_surv_obj         - Surv object for overall survival
#   overall_km_fit           - Kaplan-Meier fit for all samples
#   overall_survival_summary - One-row cohort survival summary
#   clinical_cox_fit         - Cox model using clinical covariates
#   clinical_cox_results     - Tidy hazard-ratio table for clinical model
#
# Outputs:
#   results/overall_survival_summary.csv
#   results/clinical_cox_results.csv
#   figures/overall_kaplan_meier.png


# Validate inputs -------------------------------------------------------------

if (!exists("clinical_rppa_mutation")) {
   stop("clinical_rppa_mutation is missing. Source 05_integrate_data.R first.")
}

required_survival_cols <- c(
   "sample_id", "patient_id", "os_months", "os_event",
   "age", "sex", "stage", "grade"
)
missing_survival_cols <- setdiff(required_survival_cols, names(clinical_rppa_mutation))

if (length(missing_survival_cols) > 0) {
   stop(
      "Integrated dataset is missing survival/clinical column(s): ",
      paste(missing_survival_cols, collapse = ", ")
   )
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)


# Prepare survival analysis table --------------------------------------------

survival_data <- clinical_rppa_mutation %>%
   dplyr::filter(!is.na(.data$os_months), !is.na(.data$os_event)) %>%
   dplyr::mutate(
      os_event = as.integer(.data$os_event),
      sex      = droplevels(factor(.data$sex)),
      # Treat ordered clinical categories as reference-coded factors in Cox
      # models. This gives interpretable contrasts such as STAGE IV vs STAGE I
      # instead of polynomial terms like stage.L and stage.Q.
      stage    = stats::relevel(droplevels(factor(.data$stage, ordered = FALSE)), ref = "STAGE I"),
      grade    = stats::relevel(droplevels(factor(.data$grade, ordered = FALSE)), ref = "G1")
   )

if (nrow(survival_data) == 0) {
   stop("No rows with complete survival time and event status are available.")
}

stopifnot(
   "os_event must be binary 0/1" = all(survival_data$os_event %in% c(0L, 1L))
)


# Overall Kaplan-Meier summary ------------------------------------------------

overall_surv_obj <- survival::Surv(
   time  = survival_data$os_months,
   event = survival_data$os_event
)

overall_km_fit <- survival::survfit(overall_surv_obj ~ 1, data = survival_data)
km_table <- overall_km_fit$table

median_survival_months <- unname(as.numeric(km_table[["median"]]))
if (length(median_survival_months) == 0 || is.nan(median_survival_months)) {
   median_survival_months <- NA_real_
}

overall_survival_summary <- survival_data %>%
   dplyr::summarise(
      samples                = dplyr::n(),
      patients               = dplyr::n_distinct(.data$patient_id),
      events                 = sum(.data$os_event),
      censored               = dplyr::n() - sum(.data$os_event),
      median_follow_up_months = median(.data$os_months, na.rm = TRUE),
      median_survival_months = median_survival_months
   )

readr::write_csv(overall_survival_summary, "results/overall_survival_summary.csv")

message("Overall survival summary:")
print(overall_survival_summary)


# Clinical baseline Cox model -------------------------------------------------
# Use complete cases for the clinical covariates so coxph uses a transparent
# and reproducible analysis set.

clinical_model_data <- survival_data %>%
   dplyr::select(os_months, os_event, age, sex, stage, grade) %>%
   tidyr::drop_na()

if (nrow(clinical_model_data) < 20) {
   warning("Clinical Cox model has fewer than 20 complete cases.")
}

clinical_cox_fit <- survival::coxph(
   survival::Surv(os_months, os_event) ~ age + sex + stage + grade,
   data = clinical_model_data,
   ties = "efron"
)

clinical_cox_results <- broom::tidy(
   clinical_cox_fit,
   conf.int    = TRUE,
   exponentiate = TRUE
) %>%
   dplyr::transmute(
      term,
      hazard_ratio = .data$estimate,
      conf_low     = .data$conf.low,
      conf_high    = .data$conf.high,
      p_value      = .data$p.value
   ) %>%
   dplyr::arrange(.data$p_value)

readr::write_csv(clinical_cox_results, "results/clinical_cox_results.csv")

message("Clinical Cox model complete cases: ", nrow(clinical_model_data))
message("Clinical Cox results:")
print(clinical_cox_results)


# Save Kaplan-Meier plot ------------------------------------------------------

km_plot <- survminer::ggsurvplot(
   overall_km_fit,
   data = survival_data,
   risk.table = TRUE,
   legend = "none",
   xlab = "Time (months)",
   ylab = "Overall survival probability",
   title = "TCGA KIRC overall survival"
)

png("figures/overall_kaplan_meier.png", width = 1800, height = 1600, res = 220)
print(km_plot)
dev.off()

message("Saved overall Kaplan-Meier plot to figures/overall_kaplan_meier.png")
