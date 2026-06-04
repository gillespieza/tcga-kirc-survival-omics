# Prepare clinical survival data ----------------------------------------------
#
# Merges sample-level and patient-level clinical tables into a single
# analysis-ready survival tibble. Renames columns to tidy snake_case names,
# coerces types, derives the overall survival event indicator, encodes
# categorical variables as ordered or unordered factors with explicit levels,
# removes records with missing survival information, and prints a summary of
# the resulting cohort.
#
# Requires: load_data.R to have been sourced so that clinical_patient and
#           clinical_sample are available in the global environment.
#
# Produces:
# clinical_data             - Wide tibble from left-joining clinical_patient
#                             onto clinical_sample by PATIENT_ID. One row
#                             per sample.
# clinical_survival         - Analysis-ready survival tibble with tidy column
#                             names. Records with missing os_months or
#                             os_event are dropped. One row per sample.
#                             Columns:
#                               patient_id - TCGA patient barcode
#                               sample_id  - TCGA sample barcode
#                               os_months  - Overall survival time (months)
#                               os_event   - Event indicator (1 = death,
#                                            0 = censored)
#                               age        - Age at diagnosis (years)
#                               sex        - Factor
#                               stage      - Ordered factor
#                                            (STAGE I < II < III < IV)
#                               grade      - Ordered factor (G1 < G2 < G3 < G4);
#                                            GX (not assessable) recoded to NA
# clinical_survival_summary - Single-row tibble of patient/event counts and
#                             median follow-up duration.
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.

# Validate join keys -----------------------------------------------------------
# Only the join keys are checked here. All other required columns are validated
# against clinical_data after the join, since column provenance (patient vs
# sample table) varies by study.

if (!"PATIENT_ID" %in% names(clinical_patient)) {
  stop(
    "clinical_patient is missing required join key: PATIENT_ID",
    call. = FALSE
  )
}

if (!all(c("PATIENT_ID", "SAMPLE_ID") %in% names(clinical_sample))) {
  stop(
    "clinical_sample is missing required join key(s): ",
    "PATIENT_ID and/or SAMPLE_ID",
    call. = FALSE
  )
}

# Combine patient-level and sample-level clinical data -------------------------
# Keep clinical_sample on the left so that the joined output remains one row
# per sample. Patient-level variables, including survival outcomes, are added
# from clinical_patient via PATIENT_ID.

n_samples_before <- nrow(clinical_sample)

clinical_data <- clinical_sample |>
  dplyr::left_join(clinical_patient, by = "PATIENT_ID")

# Warn if the join changed the row count. That would suggest duplicate
# PATIENT_ID values in clinical_patient, because multiple samples per patient
# are already expected in clinical_sample and should not by themselves change
# the number of rows in this join.

if (nrow(clinical_data) != n_samples_before) {
  warning(
    "Row count changed after join: ",
    n_samples_before, " sample rows -> ",
    nrow(clinical_data), " joined rows. ",
    "This suggests PATIENT_ID may not be unique in clinical_patient.",
    call. = FALSE
  )
}

# Validate all columns needed for transmute against the joined table.
# Checked here rather than before the join since columns may come from either
# source table.

required_cols <- c(
  "PATIENT_ID",
  "SAMPLE_ID",
  "OS_MONTHS",
  "OS_STATUS",
  "AGE",
  "SEX",
  "AJCC_PATHOLOGIC_TUMOR_STAGE",
  "GRADE"
)

missing_cols <- setdiff(required_cols, names(clinical_data))

if (length(missing_cols) > 0) {
  stop(
    "Joined clinical_data is missing expected column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Joined clinical tables: ",
  nrow(clinical_data), " rows x ",
  ncol(clinical_data), " cols."
)

# Prepare clinical survival table ----------------------------------------------
# Expected factor levels follow TCGA / cBioPortal conventions for this study.
# GX (grade not assessable) is intentionally excluded from the ordered grade
# levels and recoded to NA, as it cannot be meaningfully ranked alongside G1-G4.

stage_levels <- c("STAGE I", "STAGE II", "STAGE III", "STAGE IV")
grade_levels <- c("G1", "G2", "G3", "G4")

# Check for unexpected raw values before transmute, so warnings reflect actual
# problematic strings rather than just a post-hoc count of NAs.

unexpected_stages <- setdiff(
  unique(stats::na.omit(clinical_data$AJCC_PATHOLOGIC_TUMOR_STAGE)),
  stage_levels
)

unexpected_grades <- setdiff(
  unique(stats::na.omit(clinical_data$GRADE)),
  c(grade_levels, "GX")
)

if (length(unexpected_stages) > 0) {
  warning(
    "Unexpected AJCC stage value(s) will be coerced to NA: ",
    paste(unexpected_stages, collapse = ", "),
    call. = FALSE
  )
}

if (length(unexpected_grades) > 0) {
  warning(
    "Unexpected GRADE value(s) will be coerced to NA: ",
    paste(unexpected_grades, collapse = ", "),
    call. = FALSE
  )
}

n_before_filter <- nrow(clinical_data)

clinical_survival <- clinical_data |>
  dplyr::transmute(
    patient_id = .data$PATIENT_ID,
    sample_id = .data$SAMPLE_ID,
    os_months = as.numeric(.data$OS_MONTHS),
    os_event = dplyr::case_when(
      is.na(.data$OS_STATUS) ~ NA_integer_,
      stringr::str_detect(
        .data$OS_STATUS,
        stringr::regex("DECEASED", ignore_case = TRUE)
      ) ~ 1L,
      stringr::str_detect(
        .data$OS_STATUS,
        stringr::regex("LIVING", ignore_case = TRUE)
      ) ~ 0L,
      TRUE ~ NA_integer_
    ),
    age = as.numeric(.data$AGE),
    sex = factor(.data$SEX),
    stage = factor(
      .data$AJCC_PATHOLOGIC_TUMOR_STAGE,
      levels  = stage_levels, # Enforces the natural order sequence
      ordered = FALSE # Prevents polynomial contrasts in Cox models
    ) |> droplevels(), # prevents missing/excluded categories (like GX)
    grade = factor(
      dplyr::if_else(.data$GRADE == "GX", NA_character_, .data$GRADE),
      levels = grade_levels,
      ordered = FALSE # Prevents polynomial contrasts in Cox models
    ) |> droplevels() # prevents missing/excluded categories (like GX)
  ) |>
  dplyr::filter(
    !is.na(.data$os_months),
    !is.na(.data$os_event)
  )

n_dropped <- n_before_filter - nrow(clinical_survival)

message(
  "Removed ",
  n_dropped,
  " record(s) with missing os_months or os_event."
)

# Summarise available survival data -------------------------------------------

clinical_survival_summary <- clinical_survival |>
  dplyr::summarise(
    patients                = dplyr::n_distinct(patient_id),
    samples                 = dplyr::n_distinct(sample_id),
    events                  = sum(os_event),
    median_follow_up_months = median(os_months, na.rm = TRUE)
  )

message("Clinical survival table prepared.")
print(clinical_survival_summary)
