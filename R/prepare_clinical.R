# Prepare clinical survival data ----------------------------------------------
#
# Merges sample-level and patient-level clinical tables into a single
# analysis-ready survival tibble. Renames all candidate clinical columns to
# tidy snake_case names, coerces types, derives the overall survival event
# indicator, encodes categorical variables as factors, and removes records
# with missing survival information.
#
# Unlike the previous fixed-column approach, this script retains all candidate
# clinical and molecular variables for downstream univariable Cox screening in
# screen_clinical.R. The fixed set of four covariates (age, sex, stage, grade)
# has been replaced by a broader set of clinical candidates.
#
# Excluded variables and reasons:
#   - Identifiers/admin:
#       PATIENT_ID, OTHER_PATIENT_ID, CANCER_TYPE_ACRONYM, AJCC_STAGING_EDITION,
#       FORM_COMPLETION_DATE, INFORMED_CONSENT_VERIFIED,
#       IN_PANCANPATHWAYS_FREEZE, ICD_10, ICD_O_3_HISTOLOGY, ICD_O_3_SITE,
#       ONCOTREE_CODE, CANCER_TYPE, CANCER_TYPE_DETAILED, TUMOR_TYPE,
#       TISSUE_SOURCE_SITE_CODE, TISSUE_SOURCE_SITE, SOMATIC_STATUS, SAMPLE_ID
#   - Outcome variables (circular):
#       OS_STATUS, OS_MONTHS, DSS_STATUS/MONTHS, DFS_STATUS/MONTHS,
#       PFS_STATUS/MONTHS, DAYS_LAST_FOLLOWUP
#   - Redundant with AGE:
#       DAYS_TO_BIRTH, DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS
#   - Post-baseline/data leakage:
#       PERSON_NEOPLASM_CANCER_STATUS, NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT
#   - Constant in TCGA-KIRC (zero variance):
#       TUMOR_TISSUE_SITE (all kidney), SAMPLE_TYPE (all primary tumour),
#       SUBTYPE (all KIRC)
#   - Near-constant in TCGA-KIRC:
#       RADIATION_THERAPY, HISTORY_NEOADJUVANT_TRTYN
#   - Socioeconomic/healthcare confounders:
#       ETHNICITY, RACE, GENETIC_ANCESTRY_LABEL
#   - Collection process artefacts (selection bias, not biology):
#       TISSUE_PROSPECTIVE_COLLECTION_INDICATOR,
#       TISSUE_RETROSPECTIVE_COLLECTION_INDICATOR
#   - Procedural variables (no direct biological mechanism on survival):
#       PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT — binary flag recording
#         whether lymph nodes were examined at surgery; univariable association
#         with survival is confounded entirely by stage and disappears in
#         multivariable models
#       PRIOR_DX — records a prior cancer diagnosis; lacks a plausible direct
#         biological effect on post-diagnosis ccRCC survival and is more likely
#         a proxy for healthcare access or surveillance intensity.
#
# NOTE on TNM vs composite stage: PATH_T_STAGE, PATH_N_STAGE, PATH_M_STAGE,
# and AJCC_PATHOLOGIC_TUMOR_STAGE are all retained as candidates. They are
# screened independently in screen_clinical.R; if multiple TNM components
# survive BH correction, only the composite stage is retained in the final
# model to avoid collinearity.
#
# Requires: load_data.R to have been sourced so that clinical_patient and
#           clinical_sample are available in the global environment.
#
# Produces:
#   clinical_survival         - Analysis-ready tibble. One row per sample.
#                               Outcome columns: patient_id, sample_id,
#                               os_months, os_event.
#                               Candidate predictor columns (snake_case):
#                               age, sex, stage, path_t_stage, path_n_stage,
#                               path_m_stage, weight, subtype, ethnicity,
#                               race, genetic_ancestry_label,
#                               radiation_therapy, neoadjuvant_treatment,
#                               prior_dx, lymph_node_assessment, grade,
#                               aneuploidy_score, msi_score_mantis,
#                               msi_sensor_score, tmb_nonsynonymous,
#                               tbl_score, tissue_prospective_collection,
#                               tissue_retrospective_collection,
#                               tumor_tissue_site, sample_type.
#   clinical_candidate_cols   - Character vector of candidate predictor column
#                               names present in clinical_survival. Used by
#                               screen_clinical.R and downstream modelling
#                               scripts to identify clinical columns to screen.
#   clinical_survival_summary - Single-row tibble of patient/event counts and
#                               median follow-up duration.
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.


# Validate Join Keys -----------------------------------------------------------

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


# Combine Patient-Level and Sample-Level Clinical Data -------------------------

# Record row count before joining to detect unexpected duplications.
n_samples_before <- nrow(clinical_sample)

# Left join keeps clinical_sample as the anchor (one row per sample).
# Patient-level variables including survival outcomes are added from
# clinical_patient via PATIENT_ID.
clinical_data <- clinical_sample |>
  dplyr::left_join(clinical_patient, by = "PATIENT_ID")

if (nrow(clinical_data) != n_samples_before) {
  warning(
    "Row count changed after join: ",
    n_samples_before, " sample rows -> ",
    nrow(clinical_data), " joined rows. ",
    "This suggests PATIENT_ID may not be unique in clinical_patient.",
    call. = FALSE
  )
}

# Hard-check outcome columns — the pipeline cannot proceed without these.
required_outcome_cols <- c("PATIENT_ID", "SAMPLE_ID", "OS_MONTHS", "OS_STATUS")

missing_outcome_cols <- setdiff(required_outcome_cols, names(clinical_data))

if (length(missing_outcome_cols) > 0L) {
  stop(
    "Joined clinical_data is missing required outcome column(s): ",
    paste(missing_outcome_cols, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Joined clinical tables: ",
  nrow(clinical_data), " rows x ",
  ncol(clinical_data), " cols."
)


# Define Candidate Column Map --------------------------------------------------
# Maps each raw cBioPortal column name to its tidy snake_case output name and
# expected type. Columns absent from the joined data are silently skipped so
# the script is portable across TCGA studies with different column coverage.

candidate_col_map <- tibble::tibble(
  raw_name = c(
    # --- Standard clinical covariates ---
    "AGE",
    "SEX",
    "WEIGHT", # body weight; obesity linked to RCC outcomes

    # --- Pathological staging ---
    "AJCC_PATHOLOGIC_TUMOR_STAGE", # composite AJCC stage (screened alongside TNM)
    "PATH_T_STAGE", # TNM T: primary tumour size and local invasion
    "PATH_N_STAGE", # TNM N: regional lymph node involvement
    "PATH_M_STAGE", # TNM M: distant metastasis
    "GRADE", # histological nuclear grade

    # --- Genomic instability and mutational burden ---
    "ANEUPLOIDY_SCORE", # chromosomal instability; somatic copy-number burden
    "TBL_SCORE", # tumour break load; structural variant burden
    "TMB_NONSYNONYMOUS", # tumour mutational burden; non-synonymous somatic variants

    # --- Microsatellite instability ---
    "MSI_SCORE_MANTIS", # MSI score from MANTIS algorithm
    "MSI_SENSOR_SCORE" # MSI score from MSIsensor algorithm
  ),
  output_name = c(
    "age",
    "sex",
    "weight",
    "stage",
    "path_t_stage",
    "path_n_stage",
    "path_m_stage",
    "grade",
    "aneuploidy_score",
    "tbl_score",
    "tmb_nonsynonymous",
    "msi_score_mantis",
    "msi_sensor_score"
  ),
  col_type = c(
    "numeric", # age
    "factor", # sex
    "numeric", # weight
    "factor_special", # stage — explicit levels applied below
    "factor", # path_t_stage
    "factor", # path_n_stage
    "factor", # path_m_stage
    "factor_special", # grade — G1/G2 collapse applied below
    "numeric", # aneuploidy_score
    "numeric", # tbl_score
    "numeric", # tmb_nonsynonymous
    "numeric", # msi_score_mantis
    "numeric" # msi_sensor_score
  )
)

# Restrict to columns actually present in the joined data.
available_candidates <- candidate_col_map |>
  dplyr::filter(.data$raw_name %in% names(clinical_data))

skipped_candidates <- setdiff(
  candidate_col_map$raw_name,
  names(clinical_data)
)

if (length(skipped_candidates) > 0L) {
  message(
    length(skipped_candidates),
    " candidate column(s) not found in the joined data and will be skipped: ",
    paste(skipped_candidates, collapse = ", ")
  )
}


# Prepare Clinical Survival Table ----------------------------------------------
# Build a named vector for dplyr::rename() in the format c(new = "old").

rename_map <- stats::setNames(
  available_candidates$raw_name,
  available_candidates$output_name
)

# Identify output names by type for use in across() calls below.
numeric_cols <- available_candidates$output_name[
  available_candidates$col_type == "numeric"
]
plain_factor_cols <- available_candidates$output_name[
  available_candidates$col_type == "factor"
]

# Factor levels for variables with a known reference group.
# GX (grade not assessable) is recoded to NA before factoring since it
# cannot be meaningfully ranked alongside G1-G4 in a survival model.
stage_levels <- c("STAGE I", "STAGE II", "STAGE III", "STAGE IV")
grade_levels <- c("G1/G2", "G3", "G4")

# Warn about unexpected raw values before coercion so the warning is
# informative rather than appearing as silent NA introduction later.
if ("AJCC_PATHOLOGIC_TUMOR_STAGE" %in% names(clinical_data)) {
  unexpected_stages <- setdiff(
    unique(stats::na.omit(clinical_data$AJCC_PATHOLOGIC_TUMOR_STAGE)),
    stage_levels
  )
  if (length(unexpected_stages) > 0L) {
    warning(
      "Unexpected AJCC stage value(s) will be coerced to NA: ",
      paste(unexpected_stages, collapse = ", "),
      call. = FALSE
    )
  }
}

if ("GRADE" %in% names(clinical_data)) {
  unexpected_grades <- setdiff(
    unique(stats::na.omit(clinical_data$GRADE)),
    c("G1", "G2", "G3", "G4", "GX")
  )
  if (length(unexpected_grades) > 0L) {
    warning(
      "Unexpected GRADE value(s) will be coerced to NA: ",
      paste(unexpected_grades, collapse = ", "),
      call. = FALSE
    )
  }
}

n_before_filter <- nrow(clinical_data)

clinical_survival <- clinical_data |>
  # Step 1: Derive outcome and standardised ID columns from raw values.
  dplyr::mutate(
    patient_id = .data$PATIENT_ID,
    sample_id = standardise_sample_id(.data$SAMPLE_ID),
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
    )
  ) |>
  # Step 2: Select outcome columns and all available raw candidate columns.
  dplyr::select(
    patient_id,
    sample_id,
    os_months,
    os_event,
    dplyr::all_of(available_candidates$raw_name)
  ) |>
  # Step 3: Rename raw candidate columns to tidy snake_case output names.
  dplyr::rename(dplyr::all_of(rename_map)) |>
  # Step 4: Coerce numeric candidate columns.
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(numeric_cols),
      as.numeric
    )
  ) |>
  # Step 5: Coerce plain factor columns (no pre-specified levels).
  # This covers TNM components, sex, race, ethnicity, etc.
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(plain_factor_cols),
      ~ factor(as.character(.x))
    )
  ) |>
  # Step 6: Apply explicit factor levels to composite AJCC stage.
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of("stage"),
      ~ factor(.x, levels = stage_levels, ordered = FALSE) |> droplevels()
    )
  ) |>
  # Step 7: Collapse G1/G2 and apply explicit factor levels to grade.
  # G1 and G2 are merged into a single reference group to prevent
  # quasi-complete separation in Cox models (very few G1 events).
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of("grade"),
      ~ dplyr::case_when(
        .x == "GX" ~ NA_character_,
        .x %in% c("G1", "G2") ~ "G1/G2",
        TRUE ~ as.character(.x)
      ) |>
        factor(levels = grade_levels, ordered = FALSE) |>
        droplevels()
    )
  ) |>
  # Step 8: Remove records where survival outcome is unknown.
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


# Define Clinical Candidate Predictor Vector -----------------------------------
# clinical_candidate_cols is a character vector of all candidate predictor
# column names that are present in clinical_survival. It is used by
# screen_clinical.R to know which columns to pass to univariable Cox models,
# and by downstream modelling scripts to correctly exclude non-predictor columns
# (patient_id, sample_id, os_months, os_event) from feature lists.

clinical_candidate_cols <- intersect(
  available_candidates$output_name,
  names(clinical_survival)
)

message(
  "Candidate clinical predictors available for screening: ",
  length(clinical_candidate_cols),
  " (",
  paste(clinical_candidate_cols, collapse = ", "),
  ")"
)


# Summarise Available Survival Data -------------------------------------------

clinical_survival_summary <- clinical_survival |>
  dplyr::summarise(
    patients                = dplyr::n_distinct(.data$patient_id),
    samples                 = dplyr::n_distinct(.data$sample_id),
    events                  = sum(.data$os_event),
    median_follow_up_months = stats::median(.data$os_months, na.rm = TRUE)
  )

message("Clinical survival table prepared.")
# print(clinical_survival_summary)
