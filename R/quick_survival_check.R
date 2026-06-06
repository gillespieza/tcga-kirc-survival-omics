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
# Factor levels and baseline reference positions are preserved exactly as
# configured during the clinical preparation script.

survival_data <- clinical_rppa_rna_mutation |>
  dplyr::mutate(
    os_event = as.integer(.data$os_event)
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

clinical_vars <- c("os_months", "os_event", "age", "sex", "stage", "grade")

clinical_model_data <- enforce_complete_cases(
  data      = survival_data,
  variables = clinical_vars
)

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
c_index <- cox_summary$concordance[["C"]]

message("Clinical Cox C-index: ", round(c_index, 3))

clinical_cox_results <- broom::tidy(
  clinical_cox_fit,
  conf.int     = TRUE,
  exponentiate = TRUE
) |>
  dplyr::transmute(
    term = .data$term,
    hazard_ratio = .data$estimate,
    conf_low = .data$conf.low,
    conf_high = .data$conf.high,
    p_value = .data$p.value
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
  legend     = "none",
  
  # put a dotted line horizontal and vertical for median
  surv.median.line = "hv",
  
  # x-axis legend
  xlab       = "Time (months)",
  
  # y-axis legend
  ylab       = "Overall survival probability",
  
  # plot title
  title      = "TCGA KIRC: overall survival",
  
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

# Stage-stratified KM curve
# Log-rank test p-value and method label are included to indicate whether
# stage groups differ significantly in survival.
stage_km_fit <- survival::survfit(
  overall_surv_obj ~ stage,
  data = survival_data
)


# see comments in overall_km_plot for explanations of the parameters used here, 
# which are mostly the same.
stage_km_plot <- survminer::ggsurvplot(
   stage_km_fit,
   data = survival_data,
   conf.int = FALSE,
   risk.table = TRUE,
   pval = TRUE,
   pval.method  = TRUE,
   legend.title = "Stage",                
   legend.labs = levels(survival_data$stage),
   xlab = "Time (months)",
   ylab = "Overall survival probability",
   title = "TCGA KIRC: overall survival by AJCC stage",
   risk.table.y.text.col = TRUE, 
   risk.table.y.text = FALSE,
   risk.table.height = 0.3,
   palette = "npg",
   linewidth = 1.0,
   tables.theme = theme_classic()
)

save_pipeline_plot(
  plot_object = stage_km_plot,
  file_path   = "figures/stage_kaplan_meier.png",
  width       = 1800,
  height      = 1800,
  resolution  = 220
)

message("Saved stage-stratified KM plot to figures/stage_kaplan_meier.png")

# Grade-stratified KM curve
# Log-rank test p-value and method label are included to indicate whether
# stage groups differ significantly in survival.
grade_km_fit <- survival::survfit(
  overall_surv_obj ~ grade,
  data = survival_data
)


# see comments in overall_km_plot for explanations of the parameters used here, 
# which are mostly the same.
grade_km_plot <- survminer::ggsurvplot(
   grade_km_fit,
   data = survival_data,
   conf.int = FALSE,
   risk.table = TRUE,
   pval = TRUE,
   pval.method  = TRUE,
   legend.title = "Grade",                
   legend.labs = levels(survival_data$grade),
   xlab = "Time (months)",
   ylab = "Overall survival probability",
   title = "TCGA KIRC: overall survival by histological grade",
   risk.table.y.text.col = TRUE, 
   risk.table.y.text = FALSE,
   risk.table.height = 0.3,
   palette = "npg",
   linewidth = 1.0,
   tables.theme = theme_classic()
)

save_pipeline_plot(
  plot_object = grade_km_plot,
  file_path   = "figures/grade_kaplan_meier.png",
  width       = 1800,
  height      = 1800,
  resolution  = 220
)

message("Saved grade-stratified KM plot to figures/grade_kaplan_meier.png")



# EXPLORATORY P-VALUE FEATURE SELECTION --------------------------------------------------
# filter omics features at a statistical significance level 0.1 or 0.2 (larger than 0.05 
# to reduce false negative identification of omics features in multivariate analysis). 

# Use the raw RNA-seq data source so this step runs over the full gene matrix.
check_required_objects("rnaseq_data")
check_has_columns("rnaseq_data", "Hugo_Symbol")

rnaseq_all_expression <- rnaseq_data |>
  dplyr::mutate(
    row_mean = rowMeans(
      as.matrix(dplyr::pick(-dplyr::all_of("Hugo_Symbol"))),
      na.rm = TRUE
    )
  ) |>
  dplyr::arrange(dplyr::desc(.data$row_mean)) |>
  dplyr::distinct(.data$Hugo_Symbol, .keep_all = TRUE) |>
  dplyr::select(-dplyr::all_of("row_mean")) |>
  tidyr::pivot_longer(
    cols = -dplyr::all_of("Hugo_Symbol"),
    names_to = "sample_id",
    values_to = "rsem"
  ) |>
  dplyr::mutate(
    sample_id = standardise_sample_id(.data$sample_id),
    rsem = tidyr::replace_na(.data$rsem, 0),
    log2_expr = log2(.data$rsem + 1)
  ) |>
  dplyr::select(.data$sample_id, .data$Hugo_Symbol, .data$log2_expr) |>
  tidyr::pivot_wider(
    names_from = "Hugo_Symbol",
    values_from = "log2_expr"
  ) |>
  dplyr::rename_with(~ paste0("rna_", .x), -dplyr::all_of("sample_id"))

cohort_rna <- clinical_survival |>
  inner_join(rnaseq_all_expression, by = "sample_id")

gene_names <- names(cohort_rna)[startsWith(names(cohort_rna), "rna_")]
n_genes <- length(gene_names)

rna_pvalues <- numeric(length(gene_names))

message("Now running p-value feature selection for ", n_genes, " genes...")

pb <- utils::txtProgressBar(min = 0, max = n_genes, style = 3)

for (i in seq_along(gene_names)) {
  gene <- gene_names[i]
  formula <- as.formula(paste0(
    "survival::Surv(os_months, os_event) ~ `",
    gene,
    "`"
  ))
  fit <- survival::coxph(formula, data = cohort_rna)
  rna_pvalues[i] <- summary(fit)$coefficients[, "Pr(>|z|)"]
  utils::setTxtProgressBar(pb, i)
}

close(pb)

gene_pvalues <- tibble::tibble(
  gene = gene_names,
  p_value = rna_pvalues
)

significant_genes <- gene_pvalues |>
  dplyr::filter(p_value < 0.1)

# significant_genes <- sub("^rna_", "", significant_genes)

message(
  "Found ",
  nrow(significant_genes),
  " genes with p < 0.1."
)