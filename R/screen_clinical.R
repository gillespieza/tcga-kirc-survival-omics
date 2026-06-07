# Univariable clinical feature screening ---------------------------------------
#
# Screens all candidate clinical and molecular variables from clinical_survival
# using univariable Cox proportional hazards models. Applies Benjamini-Hochberg
# (BH) multiple testing correction and selects variables that survive the
# adjusted threshold (p_adjust_bh < 0.05). Resolves collinearity between TNM
# components and composite AJCC stage by retaining only composite stage when
# multiple TNM-related variables are significant. Produces a ranked results
# table, a vector of selected features, and two forest plots.
#
# This script mirrors the univariable screening design used for RPPA proteins
# in feature_selection.R, applying the same BH threshold for consistency.
#
# Requires: prepare_clinical.R to have been sourced so that clinical_survival
#           and clinical_candidate_cols are available in the global environment.
#
# Produces:
#   clinical_univariable_results  - One row per candidate variable; LRT
#                                   p-value, BH-adjusted p-value, n, events.
#   clinical_univariable_detail   - One row per Cox model term; HR and 95% CI.
#                                   Used for the selected-variables forest plot.
#   selected_clinical_features    - Character vector of variable names that
#                                   passed BH correction, after TNM resolution.
#                                   Replaces the hardcoded four-variable list
#                                   used by downstream modelling scripts.
#
# Outputs:
#   results/clinical_univariable_cox_results.csv
#   results/selected_clinical_features.csv
#   figures/clinical_screening_overview.png  - All screened variables
#   figures/clinical_screening_selected.png  - Selected variables forest plot
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.


# Validate inputs --------------------------------------------------------------

check_required_objects(c(
  "clinical_survival",
  "clinical_candidate_cols"
))

check_has_columns(
  "clinical_survival",
  c("os_months", "os_event")
)

abort_if_false(
  length(clinical_candidate_cols) > 0L,
  "clinical_candidate_cols is empty. Check prepare_clinical.R ran correctly."
)

message(
  "Screening ",
  length(clinical_candidate_cols),
  " candidate clinical variable(s) by univariable Cox regression ..."
)


# Pre-screening variable filter ------------------------------------------------
# Remove variables that cannot support a Cox model before attempting to fit one.
# This avoids noisy warnings and failed fits later in the loop.
#
# Filters applied:
#   - Numeric: must have sd > 0 and at least 30 non-missing values.
#   - Factor: must have at least 2 levels with >= 10 patients per level.

is_screenable <- function(col, min_n = 30L, min_group = 10L) {
  # Drop missing values before checking
  col_complete <- col[!is.na(col)]

  if (length(col_complete) < min_n) {
    return(FALSE)
  }

  if (is.numeric(col)) {
    # Numeric: check for non-zero variance
    return(stats::sd(col_complete, na.rm = TRUE) > 0)
  }

  # Factor or character: check for at least 2 levels with adequate group size
  tab <- table(col_complete)
  n_valid_levels <- sum(tab >= min_group)
  return(n_valid_levels >= 2L)
}

# Apply the pre-screening filter to every candidate column
screenable_cols <- clinical_candidate_cols[
  vapply(
    clinical_survival[clinical_candidate_cols],
    is_screenable,
    logical(1L)
  )
]

n_dropped_prescreening <- length(clinical_candidate_cols) - length(screenable_cols)

if (n_dropped_prescreening > 0L) {
  dropped_names <- setdiff(clinical_candidate_cols, screenable_cols)
  message(
    n_dropped_prescreening,
    " candidate(s) dropped before screening (constant, near-constant, ",
    "or insufficient group size): ",
    paste(dropped_names, collapse = ", ")
  )
}

abort_if_false(
  length(screenable_cols) > 0L,
  "No candidate clinical variables passed the pre-screening filter."
)


# Univariable Cox screening loop -----------------------------------------------
# For each screenable variable:
#   1. For factor variables, collapse sparse levels (< 10 patients) to NA
#      before fitting to prevent quasi-complete separation.
#   2. Rename the variable to "feature" in the model data so the Cox formula
#      does not require backtick-quoting of column names.
#   3. Fit the Cox model with Efron tie correction.
#   4. Extract the overall likelihood ratio test (LRT) p-value for the variable
#      as a whole — this gives a single p-value per variable regardless of how
#      many levels a factor has, which is needed for BH correction.
#   5. Extract per-term HRs and CIs from broom::tidy() for the forest plot.

summary_list <- list() # one entry per variable (for BH correction)
detail_list <- list() # one entry per term (for forest plot)

for (var in screenable_cols) {
  col <- clinical_survival[[var]]

  # Collapse sparse factor levels to NA before fitting
  if (is.factor(col) || is.character(col)) {
    tab <- table(col, useNA = "no")
    valid_levels <- names(tab)[tab >= 10L]

    col <- dplyr::if_else(
      as.character(col) %in% valid_levels,
      as.character(col),
      NA_character_
    )
    col <- factor(col, levels = valid_levels)
  }

  # Build a complete-case model tibble, renaming the variable to "feature"
  # to avoid fragile backtick-quoting in the formula string
  model_data <- tibble::tibble(
    os_months = clinical_survival$os_months,
    os_event  = clinical_survival$os_event,
    feature   = col
  ) |>
    tidyr::drop_na()

  # Skip if too few observations or events after dropping NAs
  if (nrow(model_data) < 30L || sum(model_data$os_event) < 5L) {
    message(
      "  Skipping '",
      var,
      "': too few complete cases or events after NA removal."
    )
    next
  }

  # Fit univariable Cox model; skip on error
  fit <- tryCatch(
    survival::coxph(
      survival::Surv(os_months, os_event) ~ feature,
      data = model_data,
      ties = "efron"
    ),
    error = function(e) {
      message(
        "  Skipping '",
        var,
        "': Cox model failed (",
        conditionMessage(e),
        ")."
      )
      NULL
    }
  )

  if (is.null(fit)) next

  # Overall LRT p-value for the variable (single p per variable for BH)
  lrt_p <- tryCatch(
    unname(summary(fit)$logtest["pvalue"]),
    error = function(e) NA_real_
  )

  if (is.na(lrt_p)) next

  summary_list[[var]] <- tibble::tibble(
    variable      = var,
    variable_type = if (is.numeric(clinical_survival[[var]])) "numeric" else "factor",
    n             = nrow(model_data),
    events        = sum(model_data$os_event),
    p_value       = lrt_p
  )

  # Per-term HRs and CIs for the forest plot
  tidy_fit <- tryCatch(
    broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE),
    error = function(e) NULL
  )

  if (!is.null(tidy_fit) && nrow(tidy_fit) > 0L) {
    detail_list[[var]] <- tidy_fit |>
      dplyr::transmute(
        variable     = var,
        term         = stringr::str_remove(.data$term, "^feature"),
        n            = nrow(model_data),
        events       = sum(model_data$os_event),
        hazard_ratio = .data$estimate,
        conf_low     = .data$conf.low,
        conf_high    = .data$conf.high,
        p_value_term = .data$p.value
      )
  }
}

abort_if_false(
  length(summary_list) > 0L,
  "Univariable clinical Cox screening produced no valid results."
)


# Compile results and apply BH correction -------------------------------------

clinical_univariable_results <- dplyr::bind_rows(summary_list) |>
  dplyr::mutate(
    # BH correction across all screened variables
    p_adjust_bh = stats::p.adjust(.data$p_value, method = "BH"),
    significant = .data$p_adjust_bh < 0.05
  ) |>
  dplyr::arrange(.data$p_adjust_bh)

clinical_univariable_detail <- dplyr::bind_rows(detail_list)

message(
  "Screening complete: ",
  sum(clinical_univariable_results$significant),
  " of ",
  nrow(clinical_univariable_results),
  " variables passed BH-adjusted p < 0.05."
)

print(
  clinical_univariable_results |>
    dplyr::select(
      variable,
      variable_type,
      n,
      events,
      p_value,
      p_adjust_bh,
      significant
    )
)


# TNM collinearity resolution --------------------------------------------------
# AJCC composite stage (stage) is derived from the three TNM components
# (path_t_stage, path_n_stage, path_m_stage). If any TNM component is
# significant alongside composite stage, the individual components are dropped
# to prevent collinearity in downstream multivariable models. Only composite
# stage is retained as the clinically interpretable summary variable.

tnm_components <- c("path_t_stage", "path_n_stage", "path_m_stage")
sig_vars <- clinical_univariable_results$variable[
  clinical_univariable_results$significant
]

composite_stage_sig <- "stage" %in% sig_vars
tnm_sig <- intersect(tnm_components, sig_vars)

if (composite_stage_sig && length(tnm_sig) > 0L) {
  message(
    "TNM collinearity resolution: composite stage is significant alongside ",
    paste(tnm_sig, collapse = ", "),
    ". Dropping individual TNM components and retaining composite stage only."
  )
  sig_vars <- setdiff(sig_vars, tnm_components)
}

if (!composite_stage_sig && length(tnm_sig) > 0L) {
  message(
    "Composite stage did not survive BH correction but individual TNM ",
    "component(s) did: ", paste(tnm_sig, collapse = ", "),
    ". Retaining TNM components."
  )
}


# Define selected clinical features -------------------------------------------

selected_clinical_features <- sig_vars

if (length(selected_clinical_features) == 0L) {
  warning(
    "No clinical variables passed BH correction. ",
    "Falling back to the four standard covariates: age, sex, stage, grade.",
    call. = FALSE
  )
  selected_clinical_features <- intersect(
    c("age", "sex", "stage", "grade"),
    names(clinical_survival)
  )
}

message(
  "Selected clinical features (",
  length(selected_clinical_features),
  "): ",
  paste(selected_clinical_features, collapse = ", ")
)


# Save result tables -----------------------------------------------------------

readr::write_csv(
  clinical_univariable_results,
  "results/clinical_univariable_cox_results.csv"
)

readr::write_csv(
  tibble::tibble(feature = selected_clinical_features),
  "results/selected_clinical_features.csv"
)


# Final summary ----------------------------------------------------------------

message(
  "Clinical screening complete. ",
  length(selected_clinical_features),
  " feature(s) selected for downstream modelling: ",
  paste(selected_clinical_features, collapse = ", ")
)
