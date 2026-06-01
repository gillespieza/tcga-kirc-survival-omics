# Report-ready results and figures -------------------------------------------
#
# Creates final compact outputs for the assignment report from objects generated
# by the survival-analysis workflow.
#
# Requires: 08_survival_models.R to have been sourced.
#
# Outputs:
#   figures/model_c_index_comparison.png
#   figures/top_rppa_features.png
#   figures/top_rppa_kaplan_meier.png
#   results/report_key_findings.csv


# Validate inputs -------------------------------------------------------------

required_objects <- c(
   "survival_data", "rppa_univariable_results", "selected_rppa_features",
   "model_comparison_results"
)
missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
   stop(
      "Missing required object(s): ", paste(missing_objects, collapse = ", "),
      ". Source 08_survival_models.R first."
   )
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)


# Model comparison figure -----------------------------------------------------

model_c_index_plot <- model_comparison_results %>%
   dplyr::mutate(
      model = factor(.data$model, levels = rev(.data$model))
   ) %>%
   ggplot2::ggplot(ggplot2::aes(x = model, y = c_index)) +
   ggplot2::geom_col(fill = "#2C7FB8", width = 0.7) +
   ggplot2::geom_errorbar(
      ggplot2::aes(
         ymin = pmax(0, c_index - 1.96 * c_index_se),
         ymax = pmin(1, c_index + 1.96 * c_index_se)
      ),
      width = 0.2,
      colour = "grey30"
   ) +
   ggplot2::coord_flip() +
   ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
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

top_rppa_plot_data <- rppa_univariable_results %>%
   dplyr::slice_head(n = min(10L, dplyr::n())) %>%
   dplyr::mutate(
      feature = factor(.data$feature, levels = rev(.data$feature)),
      neg_log10_p = -log10(.data$p_value)
   )

top_rppa_plot <- top_rppa_plot_data %>%
   ggplot2::ggplot(ggplot2::aes(x = feature, y = neg_log10_p)) +
   ggplot2::geom_col(fill = "#41AB5D", width = 0.7) +
   ggplot2::coord_flip() +
   ggplot2::labs(
      x = NULL,
      y = "-log10 Cox p-value",
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

top_rppa_feature <- selected_rppa_features[[1]]

top_rppa_km_data <- survival_data %>%
   dplyr::select(os_months, os_event, dplyr::all_of(top_rppa_feature)) %>%
   tidyr::drop_na() %>%
   dplyr::mutate(
      rppa_group = dplyr::if_else(
         .data[[top_rppa_feature]] >= median(.data[[top_rppa_feature]], na.rm = TRUE),
         "High",
         "Low"
      ),
      rppa_group = factor(.data$rppa_group, levels = c("Low", "High"))
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

png("figures/top_rppa_kaplan_meier.png", width = 1800, height = 1600, res = 220)
print(top_rppa_km_plot)
dev.off()


# Compact key-findings table --------------------------------------------------

best_model <- model_comparison_results %>%
   dplyr::slice_max(.data$c_index, n = 1, with_ties = FALSE)

top_rppa <- rppa_univariable_results %>%
   dplyr::slice_head(n = 1)

report_key_findings <- tibble::tibble(
   item = c(
      "analysis_samples",
      "analysis_events",
      "top_rppa_feature",
      "top_rppa_hazard_ratio",
      "top_rppa_p_value",
      "best_model",
      "best_model_c_index"
   ),
   value = c(
      as.character(nrow(survival_data)),
      as.character(sum(survival_data$os_event)),
      top_rppa$feature,
      as.character(round(top_rppa$hazard_ratio, 3)),
      signif(top_rppa$p_value, 3),
      best_model$model,
      as.character(round(best_model$c_index, 3))
   )
)

readr::write_csv(report_key_findings, "results/report_key_findings.csv")

message("Report-ready figures and key findings saved.")
print(report_key_findings)
