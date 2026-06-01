# RPPA feature selection -------------------------------------------------------
#
# Screens RPPA protein expression features with univariable Cox models and keeps
# a compact, interpretable set of proteins for downstream survival modelling.
# Mutation features are carried forward from the curated ccRCC driver-gene set
# prepared in 04_prepare_mutations.R.
#
# Requires: 06_quick_survival_check.R to have been sourced so that survival_data
#   is available.
#
# Produces:
#   rppa_feature_cols          - Candidate RPPA feature column names
#   mutation_feature_cols      - Curated binary mutation feature column names
#   rppa_univariable_results   - Ranked univariable Cox screening results
#   selected_rppa_features     - Top RPPA features used in model comparison
#   selected_mutation_features - Mutation features used in model comparison
#
# Outputs:
#   results/rppa_univariable_cox_results.csv
#   results/selected_rppa_features.csv
#   results/selected_mutation_features.csv


# Validate inputs -------------------------------------------------------------

if (!exists("survival_data")) {
   stop("survival_data is missing. Source 06_quick_survival_check.R first.")
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

clinical_cols <- c(
   "patient_id", "sample_id", "os_months", "os_event",
   "age", "sex", "stage", "grade"
)
mutation_feature_cols <- names(survival_data)[stringr::str_starts(names(survival_data), "mut_")]
rppa_feature_cols <- setdiff(names(survival_data), c(clinical_cols, mutation_feature_cols))

# Keep numeric RPPA columns with enough observed values and variation to support
# a Cox model. This avoids errors from all-missing or constant features.
rppa_feature_cols <- rppa_feature_cols[
   vapply(survival_data[rppa_feature_cols], is.numeric, logical(1))
]
rppa_feature_cols <- rppa_feature_cols[
   vapply(
      survival_data[rppa_feature_cols],
      function(x) sum(!is.na(x)) >= 30 && stats::sd(x, na.rm = TRUE) > 0,
      logical(1)
   )
]

if (length(rppa_feature_cols) == 0) {
   stop("No usable numeric RPPA features were found for Cox screening.")
}

if (length(mutation_feature_cols) == 0) {
   warning("No mutation feature columns were found; mutation models will be skipped.")
}


# Univariable Cox screen ------------------------------------------------------

fit_univariable_rppa <- function(feature_name) {
   model_data <- survival_data %>%
      dplyr::select(os_months, os_event, dplyr::all_of(feature_name)) %>%
      tidyr::drop_na()
   
   if (nrow(model_data) < 30 || length(unique(model_data[[feature_name]])) < 2) {
      return(NULL)
   }
   
   fit <- tryCatch(
      survival::coxph(
         stats::as.formula(paste0("survival::Surv(os_months, os_event) ~ `", feature_name, "`")),
         data = model_data,
         ties = "efron"
      ),
      error = function(e) NULL
   )
   
   if (is.null(fit)) return(NULL)
   
   fit_summary <- summary(fit)
   coef_row <- as.data.frame(fit_summary$coefficients)[1, , drop = FALSE]
   
   tibble::tibble(
      feature      = feature_name,
      n            = nrow(model_data),
      events       = sum(model_data$os_event),
      coef         = coef_row$coef,
      hazard_ratio = exp(coef_row$coef),
      conf_low     = exp(coef_row$coef - 1.96 * coef_row$`se(coef)`),
      conf_high    = exp(coef_row$coef + 1.96 * coef_row$`se(coef)`),
      p_value      = coef_row$`Pr(>|z|)`
   )
}

rppa_univariable_results <- purrr::map_dfr(rppa_feature_cols, fit_univariable_rppa) %>%
   dplyr::mutate(
      p_adjust_bh = p.adjust(.data$p_value, method = "BH"),
      abs_log_hr  = abs(log(.data$hazard_ratio))
   ) %>%
   dplyr::arrange(.data$p_value, dplyr::desc(.data$abs_log_hr))

if (nrow(rppa_univariable_results) == 0) {
   stop("RPPA univariable screening did not produce any valid Cox models.")
}

n_selected_rppa <- min(10L, nrow(rppa_univariable_results))
selected_rppa_features <- rppa_univariable_results %>%
   dplyr::slice_head(n = n_selected_rppa) %>%
   dplyr::pull(.data$feature)

selected_mutation_features <- mutation_feature_cols

readr::write_csv(rppa_univariable_results, "results/rppa_univariable_cox_results.csv")
readr::write_csv(
   rppa_univariable_results %>%
      dplyr::filter(.data$feature %in% selected_rppa_features) %>%
      dplyr::select(feature, hazard_ratio, conf_low, conf_high, p_value, p_adjust_bh),
   "results/selected_rppa_features.csv"
)
readr::write_csv(
   tibble::tibble(feature = selected_mutation_features),
   "results/selected_mutation_features.csv"
)

message(
   "RPPA feature selection complete: selected ",
   length(selected_rppa_features), " of ", length(rppa_feature_cols), " usable RPPA features."
)
message("Selected RPPA features:")
print(selected_rppa_features)
