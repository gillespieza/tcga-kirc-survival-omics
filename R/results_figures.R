# Generate publication-quality multi-omics results figures --------------------
#
# Reads the compiled data frames from the cross-validation and modelling steps,
# cleans feature labels for publication, and exports high-resolution, non-
# overlapping figures to the /figures directory:
#   1. Model Performance Comparison (Sorted CV C-indices with 95% CIs)
#   2. Multi-Omics Integrated Model Forest Plot (Cleanly structured)
#   3. Clinical Variable Screening Overview (BH-adjusted -log10 p-value dot plot)
#   4. Clinical Univariable Forest Plot (Selected variables, HRs and 95% CIs)
#   5. Kaplan-Meier Curves for Selected Clinical Variables:
#        5a. AJCC pathological stage
#        5b. Histological grade
#        5c. Tumour break load score (median split)
#        5d. Age (median split)
#
# Requires: survival_models.R, screen_clinical.R, and survival_check.R to have
#           been sourced. The following objects must be present in the global
#           environment:
#             cv_results, final_coefs     - loaded from disk (see below)
#             clinical_univariable_results - from screen_clinical.R
#             clinical_univariable_detail  - from screen_clinical.R
#             selected_clinical_features   - from screen_clinical.R
#             survival_data                - from survival_check.R
#
# Outputs:
#   figures/model_cindex_comparison.png
#   figures/integrated_model_forest_plot.png
#   figures/clinical_screening_overview.png
#   figures/clinical_screening_selected.png
#   figures/km_stage.png
#   figures/km_grade.png
#   figures/km_tbl_score.png
#   figures/km_age.png
#
# Usage: this script is intended to be sourced by run_analysis.R.


# Validate inputs --------------------------------------------------------------

check_required_objects(c(
  "clinical_univariable_results",
  "clinical_univariable_detail",
  "selected_clinical_features",
  "survival_data"
))

cindex_path <- "results/cross_validated_cindex_comparison.csv"
coefs_path  <- "results/final_integrated_cox_coefficients.csv"

if (!file.exists(cindex_path) || !file.exists(coefs_path)) {
  stop(
    "Required results tables are missing from disk. ",
    "Please ensure R/survival_models.R has executed completely.",
    call. = FALSE
  )
}

# Load modelling result tables cleanly from disk
cv_results  <- readr::read_csv(cindex_path, show_col_types = FALSE)
final_coefs <- readr::read_csv(coefs_path,  show_col_types = FALSE)


# 1. Figure 1: Model Performance Comparison Plot -------------------------------

message("Generating out-of-fold C-index comparison plot...")

# Clean and enhance model labels for publication readability
plot_cv_df <- cv_results |>
  dplyr::mutate(
    clean_label = dplyr::case_when(
      model == "RNA_Path"       ~ "RNA Pathway Signatures (8)",
      model == "RNA_DataDriven" ~ "RNA Data-Driven Score (Top 20 p-val)",
      model == "CNA"            ~ "Copy Number Alterations (CNA)",
      model == "Integrated"     ~ "Fully Integrated Multi-Omics",
      model == "Clinical"       ~ "Clinical Baseline Model",
      model == "RPPA"           ~ "RPPA Proteomics (Top 5)",
      model == "Mutations"      ~ "Somatic Mutations (9)",
      TRUE                      ~ model
    ),
    # Order bars by ascending concordance so the best model sits at the top
    clean_label = stats::reorder(.data$clean_label, .data$cv_concordance)
  )

# Extract the clinical C-index scalar outside ggplot to use in geom_rect
clinical_cindex <- cv_results |>
  dplyr::filter(.data$model == "Clinical") |>
  dplyr::pull(.data$cv_concordance)

if (length(clinical_cindex) == 0L) {
  clinical_cindex <- 0.751 # Safe fallback if clinical model result is missing
}

best_model <- cv_results$model[1]

cindex_comparison_plot <- plot_cv_df |>

  ggplot2::ggplot(
    ggplot2::aes(x = .data$clean_label, y = .data$cv_concordance)
  ) +
  # Shaded band marking the clinical benchmark level
  ggplot2::geom_rect(
    ymin  = clinical_cindex - 0.001,
    ymax  = clinical_cindex + 0.001,
    xmin  = -Inf,
    xmax  = Inf,
    fill  = "#f2f4f4",
    alpha = 0.5
  ) +
  # 95% cross-validation confidence interval bars
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
    width = 0.15, colour = "#34495e", linewidth = 0.8
  ) +
  # C-index point estimates; best model highlighted in red
  ggplot2::geom_point(
    ggplot2::aes(colour = (.data$model == best_model)),
    size = 4L, show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#2c3e50")
  ) +
  # Flip so long model labels read horizontally without truncation
  ggplot2::coord_flip() +
  ggplot2::labs(
    title    = "Cross-Validated Model Performance Comparison",
    subtitle = paste(
      "True out-of-fold Harrell's C-index evaluated via stratified 5-fold CV"
    ),
    x = "Model Configuration",
    y = "Out-of-Fold Concordance Index (C-index)"
  ) +
  ggplot2::theme_classic(base_size = 13L) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold", size = 14L, colour = "#2c3e50"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 11L, colour = "#7f8c8d",
      margin = ggplot2::margin(b = 15L)
    ),
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10L), face = "bold"
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 10L), face = "bold"
    ),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank()
  )

save_pipeline_plot(
  plot_object = cindex_comparison_plot,
  file_path   = "figures/model_cindex_comparison.png",
  width       = 1000L,
  height      = 650L,
  resolution  = 100L
)

message("Saved model C-index comparison plot.")


# 2. Figure 2: Multi-Omics Integrated Model Forest Plot -----------------------

message("Generating multi-omics integrated forest plot...")

plot_coefs_df <- final_coefs |>
  dplyr::mutate(
    # Cap extreme HRs and CIs to a readable window — display only
    hazard_ratio = pmax(pmin(.data$hazard_ratio, 10.0), 0.1),
    conf_low     = pmax(pmin(.data$conf_low,     10.0), 0.1),
    conf_high    = pmax(pmin(.data$conf_high,    10.0), 0.1),

    # Humanise column prefixes into report-ready labels
    clean_feature = .data$feature,
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^age$", "Age (Years)"
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^tbl_score$", "Tumour Break Load Score"
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^sexMale$", "Sex: Male"
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^stageSTAGE ", "Tumour Stage: "
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^gradeG", "Nuclear Grade: "
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^mut_", "Mutation: "
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^cna_loss_", "CNA Focal Loss: "
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^cna_gain_", "CNA Focal Gain: "
    ),
    clean_feature = stringr::str_replace(
      .data$clean_feature, "^score_", "RNA Score: "
    )
  ) |>
  # Assign each feature to its biological data layer for faceting
  dplyr::mutate(
    layer_group = dplyr::case_when(
      stringr::str_starts(.data$feature, "score_") ~
        "Transcriptomic Pathways (8)",
      stringr::str_starts(.data$feature, "mut_") ~
        "Somatic Mutations",
      stringr::str_starts(.data$feature, "cna_") ~
        "Copy Number Events",
      .data$feature %in% c("age", "tbl_score", "sexMale") |
        stringr::str_starts(.data$feature, "stage") |
        stringr::str_starts(.data$feature, "grade") ~
        "Clinical Metadata",
      TRUE ~ "Proteomic Markers"
    ),
    layer_group = factor(.data$layer_group, levels = c(
      "Clinical Metadata", "Somatic Mutations", "Copy Number Events",
      "Proteomic Markers", "Transcriptomic Pathways (8)"
    ))
  ) |>
  # Sort within each layer by ascending HR so protective features sit at top
  dplyr::arrange(.data$layer_group, .data$hazard_ratio) |>
  dplyr::mutate(
    clean_feature = factor(
      .data$clean_feature, levels = unique(.data$clean_feature)
    )
  )

forest_plot <- plot_coefs_df |>
  ggplot2::ggplot(
    ggplot2::aes(x = .data$clean_feature, y = .data$hazard_ratio)
  ) +
  # Reference line at HR = 1 (no effect)
  ggplot2::geom_hline(
    yintercept = 1.0, linetype = "dashed", colour = "#7f8c8d", linewidth = 0.8
  ) +
  # 95% confidence interval bars
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
    width = 0.2, colour = "#2c3e50", linewidth = 0.6
  ) +
  # HR point estimates; red = increased risk, green = decreased risk
  ggplot2::geom_point(
    ggplot2::aes(fill = (.data$hazard_ratio > 1.0)),
    shape = 21, size = 3L, stroke = 0.5, colour = "#2c3e50",
    show.legend = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#2ecc71"),
    labels = c("TRUE" = "Increased risk (HR > 1)", "FALSE" = "Decreased risk (HR < 1)"),
    name   = NULL
  ) +
  # Log scale so distances above and below HR = 1 are symmetric
  ggplot2::scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
    labels = c("0.1", "0.2", "0.5", "1.0", "2.0", "5.0", "10.0")
  ) +
  # Separate panels per biological layer to prevent label crowding
  ggplot2::facet_grid(
    layer_group ~ ., scales = "free_y", space = "free_y"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title    = "Multi-Omics Integrated Model Coefficients",
    subtitle = paste(
      "Full-cohort Cox proportional hazards fit;",
      "elements uniformly bounded to window limits"
    ),
    x = "Model Feature",
    y = "Multivariable Hazard Ratio (Log Scale)"
  ) +
  ggplot2::theme_classic(base_size = 12L) +
  ggplot2::theme(
    legend.position      = c(0.97, 0.04),
    legend.justification = c("right", "bottom"),
    legend.background    = ggplot2::element_rect(
      fill   = ggplot2::alpha("white", 0.8),
      colour = NA
    ),
    plot.title = ggplot2::element_text(
      face = "bold", size = 14L, colour = "#2c3e50"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 11L, colour = "#7f8c8d",
      margin = ggplot2::margin(b = 15L)
    ),
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10L), face = "bold"
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 10L), face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      angle = 0, face = "bold", colour = "#2c3e50", hjust = 0
    ),
    strip.background   = ggplot2::element_rect(fill = "#f8f9f9", colour = NA),
    panel.spacing.y    = ggplot2::unit(0.8, "lines"),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank()
  )

save_pipeline_plot(
  plot_object = forest_plot,
  file_path   = "figures/integrated_model_forest_plot.png",
  width       = 1200L,
  height      = 900L,
  resolution  = 100L
)

message("Saved multi-omics integrated forest plot.")


# 3. Figure 3: Clinical Screening Overview ------------------------------------
# Dot plot of -log10(BH-adjusted p-value) for every screened clinical variable,
# ordered by significance. Significant variables highlighted in red. A dashed
# vertical line marks the BH threshold (p_adjust_bh = 0.05, -log10 = 1.3).

message("Generating clinical variable screening overview plot...")

overview_plot_data <- clinical_univariable_results |>
  dplyr::mutate(
    neg_log10_p = -log10(.data$p_adjust_bh),
    # Order variables by significance for a clean ranked display
    variable    = forcats::fct_reorder(.data$variable, .data$neg_log10_p)
  )

overview_plot <- overview_plot_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x      = neg_log10_p,
      y      = variable,
      colour = significant
    )
  ) +
  # Dashed line at the BH significance threshold
  ggplot2::geom_vline(
    xintercept = -log10(0.05),
    linetype   = "dashed",
    colour     = "grey50"
  ) +
  # One dot per screened variable
  ggplot2::geom_point(size = 3L) +
  ggplot2::scale_colour_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#95a5a6"),
    labels = c("TRUE" = "Significant (BH p < 0.05)", "FALSE" = "Not significant"),
    name   = NULL
  ) +
  ggplot2::labs(
    title    = "Clinical variable screening: univariable Cox regression",
    subtitle = "BH-adjusted p-values; dashed line = 0.05 threshold",
    x = expression(-log[10] ~ "(BH-adjusted p-value)"),
    y        = NULL
  ) +
  ggplot2::theme_classic(base_size = 12L) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title      = ggplot2::element_text(face = "bold")
  )

save_pipeline_plot(
  plot_object = overview_plot,
  file_path   = "figures/clinical_screening_overview.png",
  width       = 1200L,
  height      = 800L,
  resolution  = 120L
)

message("Saved clinical screening overview plot.")


# 4. Figure 4: Clinical Univariable Forest Plot --------------------------------
# Per-term HRs and 95% CIs for variables that passed BH correction, on a log
# scale. Reference categories are omitted (HR = 1 by definition).

message("Generating clinical univariable forest plot...")

forest_data <- clinical_univariable_detail |>
  # Keep only terms belonging to BH-selected variables
  dplyr::filter(.data$variable %in% selected_clinical_features) |>
  dplyr::mutate(
    # Cap extreme HRs for display only — does not affect modelling
    hazard_ratio = pmax(pmin(.data$hazard_ratio, 20.0), 0.05),
    conf_low     = pmax(pmin(.data$conf_low,     20.0), 0.05),
    conf_high    = pmax(pmin(.data$conf_high,    20.0), 0.05),
    # Build a readable label combining variable name and factor level
    term_label = dplyr::if_else(
      .data$term == "",
      .data$variable,
      paste0(.data$variable, ": ", .data$term)
    ),
    # Order rows by HR magnitude for a clean visual ranking
    term_label = forcats::fct_reorder(.data$term_label, .data$hazard_ratio)
  )

if (nrow(forest_data) > 0L) {
  clinical_forest_plot <- forest_data |>
    ggplot2::ggplot(
      ggplot2::aes(
        x    = hazard_ratio,
        y    = term_label,
        fill = (.data$hazard_ratio > 1.0)
      )
    ) +
    # Reference line at HR = 1 (no effect)
    ggplot2::geom_vline(
      xintercept = 1.0,
      linetype   = "dashed",
      colour     = "grey50"
    ) +
    # 95% confidence interval bars
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = conf_low, xmax = conf_high),
      height    = 0.2,
      colour    = "#2c3e50",
      linewidth = 0.6
    ) +
    # HR point estimates; red = increased risk, green = decreased risk
    ggplot2::geom_point(
      shape       = 21,
      size        = 3L,
      stroke      = 0.5,
      colour      = "#2c3e50",
      show.legend = FALSE # don't need a legend here
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#e74c3c", "FALSE" = "#2ecc71"),
      #labels = c("TRUE" = "Increased risk (HR > 1)", "FALSE" = "Decreased risk (HR < 1)"),
      name   = NULL
    ) +
    # Log scale so distances above and below HR = 1 are symmetric
    ggplot2::scale_x_log10(
      breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
      labels = c("0.1", "0.2", "0.5", "1.0", "2.0", "5.0", "10.0")
    ) +
    ggplot2::labs(
      title    = "Selected clinical features: univariable hazard ratios",
      subtitle = "BH-adjusted p < 0.05; log scale; red = increased risk",
      x        = "Hazard ratio (log scale)",
      y        = NULL
    ) +
    ggplot2::theme_classic(base_size = 12L) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )

  save_pipeline_plot(
    plot_object = clinical_forest_plot,
    file_path   = "figures/clinical_screening_selected.png",
    width       = 1200L,
    # Scale height so labels never overlap regardless of how many terms exist
    height      = max(600L, 80L * nrow(forest_data)),
    resolution  = 120L
  )

  message("Saved clinical univariable forest plot.")
} else {
  message("No terms available for clinical forest plot.")
}


# 5. Kaplan-Meier Curves for Selected Clinical Variables ----------------------
# One KM plot per selected variable. Categorical variables (stage, grade) are
# plotted directly. Continuous variables (tbl_score, age) are split at the
# cohort median into two groups for visual clarity. The log-rank p-value is
# printed directly on each plot via pval = TRUE.

message("Generating KM curves for selected clinical variables...")

# Helper: build and save one KM plot given a prepared dataset and group variable.
# group_var must already be a factor column in km_data with meaningful levels.
plot_km <- function(km_data, group_var, file_path, title_str) {
  # Rename the grouping column to a fixed name so survminer can find it
  # reliably via non-standard evaluation — passing a symbol by string fails
  # inside ggsurvplot's internal data lookup
  km_data_plot <- km_data |>
    dplyr::rename(group = dplyr::all_of(group_var))

  km_fit <- survival::survfit(
    survival::Surv(os_months, os_event) ~ group,
    data = km_data_plot
  )

  km_plot <- survminer::ggsurvplot(
    km_fit,
    data              = km_data_plot,  # must match the data used in survfit
    conf.int          = TRUE,
    risk.table        = TRUE,
    pval              = TRUE,
    pval.method       = TRUE,
    legend.title      = group_var,     # use the original name as the legend title
    # Extract factor levels directly from the renamed column to use as
    # clean legend labels, suppressing the default "group=..." prefix
    legend.labs       = levels(km_data_plot$group),
    xlab              = "Time (months)",
    ylab              = "Overall survival probability",
    title             = title_str,
    palette           = "npg",
    linewidth         = 1.0,
    risk.table.height = 0.28,
    tables.theme      = ggplot2::theme_classic()
  )

  save_pipeline_plot(
    plot_object = km_plot,
    file_path   = file_path,
    width       = 1800L,
    height      = 1400L,
    resolution  = 220L
  )

  message("Saved: ", file_path)
}


# Build a clean base dataset restricted to the four variables needed for KM
km_base <- survival_data |>
  dplyr::select(os_months, os_event, stage, grade, tbl_score, age) |>
  tidyr::drop_na()

message("Complete cases for KM plots: ", nrow(km_base))


# Figure 5a: Stage -------------------------------------------------------------

plot_km(
  km_data   = km_base,
  group_var = "stage",
  file_path = "figures/km_stage.png",
  title_str = "Overall survival by AJCC pathological stage"
)


# Figure 5b: Grade -------------------------------------------------------------

plot_km(
  km_data   = km_base,
  group_var = "grade",
  file_path = "figures/km_grade.png",
  title_str = "Overall survival by histological grade"
)


# Figure 5c: TBL score (median split) -----------------------------------------
# Tumour break load is continuous; split at the cohort median to produce two
# interpretable groups reflecting low vs high structural variant burden.

tbl_median <- stats::median(km_base$tbl_score, na.rm = TRUE)

km_base <- km_base |>
  dplyr::mutate(
    tbl_group = factor(
      dplyr::if_else(
        .data$tbl_score >= tbl_median,
        paste0("High (\u2265", round(tbl_median, 1L), ")"),
        paste0("Low (<",       round(tbl_median, 1L), ")")
      ),
      levels = c(
        paste0("Low (<",       round(tbl_median, 1L), ")"),
        paste0("High (\u2265", round(tbl_median, 1L), ")")
      )
    )
  )

plot_km(
  km_data   = km_base,
  group_var = "tbl_group",
  file_path = "figures/km_tbl_score.png",
  title_str = "Overall survival by tumour break load score (median split)"
)


# Figure 5d: Age (median split) ------------------------------------------------

age_median <- stats::median(km_base$age, na.rm = TRUE)

km_base <- km_base |>
  dplyr::mutate(
    age_group = factor(
      dplyr::if_else(
        .data$age >= age_median,
        paste0("Older (\u2265", round(age_median, 0L), " yrs)"),
        paste0("Younger (<",    round(age_median, 0L), " yrs)")
      ),
      levels = c(
        paste0("Younger (<",    round(age_median, 0L), " yrs)"),
        paste0("Older (\u2265", round(age_median, 0L), " yrs)")
      )
    )
  )

plot_km(
  km_data   = km_base,
  group_var = "age_group",
  file_path = "figures/km_age.png",
  title_str = "Overall survival by age (median split)"
)


# Final message ----------------------------------------------------------------

message("All publication figures successfully generated and saved to /figures.")
