# Create report-ready summary figures ------------------------------------------
#
# Generates final evaluation plots comparing performance across all five Cox
# survival models. Creates a clean, publication-quality dot-and-whisker plot
# of Concordance Indices (C-index) with their corresponding 95 % confidence
# intervals, and outputs a summary message.
#
# Requires: survival_models.R to have been sourced so that model_comparison_df
#           is available in the global environment.
#
# Produces:
#   c_index_comparison_plot - ggplot2 object showing C-indices across models
#
# Outputs:
#   figures/model_c_index_comparison.png
#
# Note: this script is intended to be sourced by run_analysis.R.


# Validate inputs -------------------------------------------------------------

check_required_objects("model_comparison_df")

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

# Standardise model display names ---------------------------------------------
# Re-level and capitalize model labels so they read nicely on a plot axis
# instead of appearing as raw code variable names.

plot_data <- model_comparison_df |>
  dplyr::mutate(
    model_clean = dplyr::case_when(
      .data$model == "clinical" ~ "Clinical Baseline",
      .data$model == "rppa" ~ "Clinical + Proteomics (RPPA)",
      .data$model == "rna" ~ "Clinical + RNA-seq Pathways",
      .data$model == "mutation" ~ "Clinical + Driver Mutations",
      .data$model == "integrated" ~ "Full Multi-omics Integration",
      TRUE ~ .data$model
    ),
    model_clean = factor(
      .data$model_clean,
      levels = c(
        "Clinical Baseline",
        "Clinical + Driver Mutations",
        "Clinical + RNA-seq Pathways",
        "Clinical + Proteomics (RPPA)",
        "Full Multi-omics Integration"
      )
    )
  )

# Generate C-index comparison plot --------------------------------------------

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
    height = 0.15,
    colour = "#2c3e50",
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
    subtitle = "TCGA-KIRC Cohort Concordance Indices (C-Index \u00b1 95% CI)",
    x        = "Harrell's Concordance Index (C-Index)"
  ) +
  ggplot2::xlim(0.45, 0.9)

# Save visualization to disk --------------------------------------------------

grDevices::png(
  filename = "figures/model_c_index_comparison.png",
  width    = 1600,
  height   = 1000,
  res      = 200
)

print(c_index_comparison_plot)

grDevices::dev.off()

# Final pipeline diagnostics --------------------------------------------------

message("Final performance visualization saved successfully.")
message("Plot location: figures/model_c_index_comparison.png")
message("Pipeline figure output complete.")
