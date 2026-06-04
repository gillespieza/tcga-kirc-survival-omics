# Create report-ready summary figures ------------------------------------------
#
# Generates final evaluation plots comparing performance across all six Cox
# survival models and visualises individual multi-omics hazard ratios. 
# Creates two clean, publication-quality figures:
#   1. A dot-and-whisker plot of Concordance Indices (C-index) with 95% CIs.
#   2. A forest plot of hazard ratios for all full integrated model predictors.
#
# Requires: survival_models.R to have been sourced so that model_comparison_df
#           and integrated_hr_df are available in the global environment.
#
# Produces:
#   c_index_comparison_plot  - ggplot2 object showing C-indices across models
#   hazard_ratio_forest_plot - ggplot2 object showing integrated model coefficients
#
# Outputs:
#   figures/model_c_index_comparison.png
#   figures/model_hazard_ratios_forest.png
#
# Note: this script is intended to be sourced by run_analysis.R.


# Validate inputs -------------------------------------------------------------

check_required_objects(c("model_comparison_df", "integrated_hr_df"))

required_metrics_cols <- c(
   "model",
   "c_index",
   "c_conf_low",
   "c_conf_high",
   "aic"
)

check_has_columns(
   "model_comparison_df",
   required_metrics_cols
)

check_has_columns(
   "integrated_hr_df",
   c("term", "hazard_ratio", "conf_low", "conf_high", "p_value")
)


# 1. Generate C-index Comparison Plot -----------------------------------------
# Re-level and capitalise model labels so they read nicely on a plot axis
# instead of appearing as raw code variable names.

plot_data <- model_comparison_df |>
   dplyr::mutate(
      model_clean = dplyr::case_when(
         .data$model == "clinical" ~ "Clinical Baseline",
         .data$model == "rppa" ~ "Clinical + Proteomics (RPPA)",
         .data$model == "rna" ~ "Clinical + RNA-seq Pathways",
         .data$model == "mutation" ~ "Clinical + Driver Mutations",
         .data$model == "cna" ~ "Clinical + CNA Alterations",
         .data$model == "integrated" ~ "Full Multi-omics Integration",
         TRUE ~ .data$model
      ),
      model_clean = factor(
         .data$model_clean,
         levels = c(
            "Clinical Baseline",
            "Clinical + Driver Mutations",
            "Clinical + CNA Alterations",
            "Clinical + RNA-seq Pathways",
            "Clinical + Proteomics (RPPA)",
            "Full Multi-omics Integration"
         )
      )
   )

c_index_comparison_plot <- ggplot2::ggplot(
   plot_data,
   ggplot2::aes(
      x    = .data$c_index,
      y    = .data$model_clean,
      xmin = .data$c_conf_low,
      xmax = .data$c_conf_high
   )
) +
   ggplot2::geom_vline(
      xintercept = 0.5,
      linetype   = "dashed",
      colour     = "grey50",
      linewidth  = 0.6
   ) +
   ggplot2::geom_errorbarh(
      height    = 0.15,
      colour    = "#2c3e50",
      linewidth = 0.8
   ) +
   ggplot2::geom_point(
      size   = 3.5,
      colour = "#e74c3c"
   ) +
   ggplot2::theme_minimal(base_size = 12) +
   ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y     = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 11)
   ) +
   ggplot2::labs(
      title    = "Survival Model Performance Comparison",
      subtitle = "TCGA-KIRC Cohort Concordance Indices (C-Index and 95% Confidence Interval)",
      x        = "Harrell's Concordance Index (C-Index)"
   ) +
   ggplot2::xlim(0.45, 0.9)


# 2. Generate Hazard Ratio Forest Plot ----------------------------------------
# Translate raw R model terms into clean, publication-ready descriptions.

hr_plot_data <- integrated_hr_df |>
   dplyr::mutate(
      term_clean = dplyr::case_when(
         .data$term == "age" ~ "Age (years)",
         .data$term == "sexMale" ~ "Sex (Male vs Female)",
         .data$term == "stageSTAGE II" ~ "Stage II vs Stage I",
         .data$term == "stageSTAGE III" ~ "Stage III vs Stage I",
         .data$term == "stageSTAGE IV" ~ "Stage IV vs Stage I",
         .data$term == "gradeG2" ~ "Grade 2 vs Grade 1",
         .data$term == "gradeG3" ~ "Grade 3 vs Grade 1",
         .data$term == "gradeG4" ~ "Grade 4 vs Grade 1",
         .data$term == "score_neutrophil_deg" ~ "RNA Neutrophil Degranulation",
         .data$term == "score_ecm_deg" ~ "RNA ECM Remodelling",
         .data$term == "score_ptk" ~ "RNA Protein Tyrosine Kinase",
         stringr::str_starts(.data$term, "mut_") ~ paste(stringr::str_remove(.data$term, "^mut_"), "Mutation"),
         stringr::str_starts(.data$term, "cna_loss_") ~ paste(stringr::str_remove(.data$term, "^cna_loss_"), "Deep Deletion"),
         stringr::str_starts(.data$term, "cna_gain_") ~ paste(stringr::str_remove(.data$term, "^cna_gain_"), "High Amplification"),
         TRUE ~ .data$term
      )
   )

hazard_ratio_forest_plot <- ggplot2::ggplot(
   hr_plot_data,
   ggplot2::aes(
      x    = .data$hazard_ratio,
      y    = stats::reorder(.data$term_clean, .data$hazard_ratio),
      xmin = .data$conf_low,
      xmax = .data$conf_high
   )
) +
   ggplot2::geom_vline(
      xintercept = 1.0,
      linetype   = "dashed",
      colour     = "grey50",
      linewidth  = 0.6
   ) +
   ggplot2::geom_errorbarh(
      height    = 0.2,
      colour    = "#2c3e50",
      linewidth = 0.8
   ) +
   ggplot2::geom_point(
      size   = 3.5,
      colour = "#2980b9"
   ) +
   ggplot2::scale_x_log10() +
   ggplot2::theme_minimal(base_size = 12) +
   ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y     = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 11)
   ) +
   ggplot2::labs(
      title    = "Multivariable Multi-omics Hazard Ratios",
      subtitle = "Full Integrated Cox Model Predictors (Hazard Ratio and 95% Confidence Interval)",
      x        = "Hazard Ratio (Log Scale)"
   )


# Save visualisations to disk and print to RStudio ---------------------------

save_pipeline_plot(
   plot_object = c_index_comparison_plot,
   file_path   = "figures/model_c_index_comparison.png",
   width       = 1600,
   height      = 1000,
   resolution  = 200
)

save_pipeline_plot(
   plot_object = hazard_ratio_forest_plot,
   file_path   = "figures/model_hazard_ratios_forest.png",
   width       = 1600,
   height      = 1300,
   resolution  = 200
)


# Final pipeline diagnostics --------------------------------------------------

message("Final performance and hazard ratio visualisations saved successfully.")
message("Outputs generated in figures/ folder.")
message("Pipeline figure output complete.")