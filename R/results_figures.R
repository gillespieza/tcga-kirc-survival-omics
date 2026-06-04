# Report-ready results and figures --------------------------------------------
#
# Creates final compact outputs for the assignment report from objects generated
# by the survival-analysis workflow.
#
# Intended to be sourced by the pipeline; not run directly.
#
# Requires: survival_models.R to have been sourced.
#
# Outputs:
#   figures/model_c_index_comparison.png
#   figures/top_rppa_features.png
#   figures/top_rppa_kaplan_meier.png
#   results/report_key_findings.csv


# Validate inputs -------------------------------------------------------------

check_required_objects(
   c(
      "survival_data",
      "rppa_univariable_results",
      "selected_rppa_features",
      "model_comparison_results"
   )
)

check_has_columns(
   "survival_data",
   c("os_months", "os_event")
)

check_has_columns(
   "rppa_univariable_results",
   c("feature", "hazard_ratio", "p_adjust_bh")
)

check_has_columns(
   "model_comparison_results",
   c("model", "c_index", "c_index_se")
)

abort_if_false(
   nrow(rppa_univariable_results) > 0,
   "rppa_univariable_results is empty."
)

abort_if_false(
   length(selected_rppa_features) > 0,
   "selected_rppa_features is empty."
)


# Model comparison figure -----------------------------------------------------
# model_comparison_results is sorted descending by c_index; factor levels are
# set in the same order so that after coord_flip() the best model appears at
# the top of the chart.

model_c_index_plot <- model_comparison_results |>
   dplyr::mutate(
      model = factor(.data$model, levels = .data$model)
   ) |>
   ggplot2::ggplot(
      ggplot2::aes(
         x = model,
         y = c_index
      )
   ) +
   ggplot2::geom_col(
      fill = "#2C7FB8",
      width = 0.7
   ) +
   ggplot2::geom_errorbar(
      ggplot2::aes(
         ymin = pmax(0, .data$c_index - 1.96 * .data$c_index_se),
         ymax = pmin(1, .data$c_index + 1.96 * .data$c_index_se)
      ),
      width = 0.2,
      colour = "grey30"
   ) +
   ggplot2::coord_flip() +
   ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
   ) +
   ggplot2::labs(
      x = NULL,
      y = "Concordance index",
      title = "Survival model comparison"
   ) +
   ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
   filename = "figures/model_c_index_comparison.png",
   plot = model_c_index_plot,
   width = 7,
   height = 4.5,
   dpi = 300
)


# Top RPPA feature figure -----------------------------------------------------
# Plot BH-adjusted p-values on the y-axis to be consistent with the FDR-gated
# feature selection applied upstream in feature_selection.R.

top_rppa_plot_data <- rppa_univariable_results |>
   dplyr::slice_head(
      n = min(10L, nrow(rppa_univariable_results))
   ) |>
   dplyr::mutate(
      feature = factor(.data$feature, levels = rev(.data$feature)),
      neg_log10_p_adj = -log10(.data$p_adjust_bh)
   )

top_rppa_plot <- top_rppa_plot_data |>
   ggplot2::ggplot(
      ggplot2::aes(
         x = feature,
         y = neg_log10_p_adj
      )
   ) +
   ggplot2::geom_col(
      fill = "#41AB5D",
      width = 0.7
   ) +
   ggplot2::coord_flip() +
   ggplot2::labs(
      x = NULL,
      y = "-log10 BH-adjusted p-value",
      title = "Top RPPA survival-associated proteins"
   ) +
   ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
   filename = "figures/top_rppa_features.png",
   plot = top_rppa_plot,
   width = 7,
   height = 4.8,
   dpi = 300
)


# Kaplan-Meier split for the top RPPA feature ---------------------------------
# Dichotomising a continuous feature loses information, so this is used only as
# an interpretable visualisation of the strongest screened RPPA association.
# selected_rppa_features is ordered by BH-adjusted p-value (ascending), so
# [[1]] is the most significant feature after FDR correction, consistent with
# the top row of rppa_univariable_results.

top_rppa_feature <- selected_rppa_features[[1]]

abort_if_false(
   top_rppa_feature %in% names(survival_data),
   paste0(
      "Top RPPA feature '",
      top_rppa_feature,
      "' is not present in survival_data."
   )
)

top_rppa_km_data <- survival_data |>
   dplyr::select(
      os_months,
      os_event,
      dplyr::all_of(top_rppa_feature)
   ) |>
   tidyr::drop_na() |>
   dplyr::mutate(
      rppa_group = dplyr::if_else(
         .data[[top_rppa_feature]] >= median(.data[[top_rppa_feature]], na.rm = TRUE),
         "High",
         "Low"
      ),
      rppa_group = factor(.data$rppa_group, levels = c("Low", "High"))
   )

abort_if_false(
   nrow(top_rppa_km_data) > 0,
   "No complete cases are available for the top RPPA Kaplan-Meier plot."
)

top_rppa_km_fit <- survival::survfit(
   survival::Surv(os_months, os_event) ~ rppa_group,
   data = top_rppa_km_data
)

top_rppa_km_plot <- survminer::ggsurvplot(
   top_rppa_km_fit,
   data = top_rppa_km_data,
   risk.table = TRUE,
   pval = TRUE,
   legend.title = top_rppa_feature,
   legend.labs = c("Low", "High"),
   xlab = "Time (months)",
   ylab = "Overall survival probability",
   title = paste("Overall survival by", top_rppa_feature, "expression")
)

grDevices::png(
   filename = "figures/top_rppa_kaplan_meier.png",
   width = 1800,
   height = 1600,
   res = 220
)
print(top_rppa_km_plot)
grDevices::dev.off()


# Compact key-findings table --------------------------------------------------

best_model <- model_comparison_results |>
   dplyr::slice_max(
      .data$c_index,
      n = 1,
      with_ties = FALSE
   )

if (
   nrow(model_comparison_results) > 1 &&
   sum(model_comparison_results$c_index == best_model$c_index) > 1
) {
   message(
      "Note: multiple models share the highest C-index (",
      round(best_model$c_index, 3),
      "); '",
      best_model$model,
      "' reported arbitrarily."
   )
}

top_rppa <- rppa_univariable_results |>
   dplyr::slice_head(n = 1)

report_key_findings <- tibble::tibble(
   item = c(
      "analysis_samples",
      "analysis_events",
      "top_rppa_feature",
      "top_rppa_hazard_ratio",
      "top_rppa_p_adjust_bh",
      "best_model",
      "best_model_c_index"
   ),
   value = as.character(
      c(
         nrow(survival_data),
         sum(survival_data$os_event),
         top_rppa$feature,
         round(top_rppa$hazard_ratio, 3),
         signif(top_rppa$p_adjust_bh, 3),
         best_model$model,
         round(best_model$c_index, 3)
      )
   )
)

readr::write_csv(
   report_key_findings,
   "results/report_key_findings.csv"
)

message("Report-ready figures and key findings saved.")
print(report_key_findings)