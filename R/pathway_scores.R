# Pathway score computation from RNA-seq expression ---------------------------
#
# Computes per-sample pathway activity scores using:
#   1. Mean log2 expression (interpretability baseline)
#   2. PC1 score (primary modelling feature; captures co-expression structure)
#
# Also evaluates:
#   - gene coverage per pathway
#   - per-sample completeness
#   - correlation with RPPA features (collinearity screening)
#
# Requires:
#   rnaseq_expression          (from prepare_rnaseq.R; sample_id + rna_* columns)
#   rnaseq_gene_set_membership (from prepare_rnaseq.R; query, gs_name, gene_symbol)
#   survival_data              (from quick_survival_check.R)
#   rppa_univariable_results   (from feature_selection.R)
#
# Produces:
#   pathway_score_results      - Per-sample RNA pathway scores
#   coverage_df                - Pathway gene coverage summary
#   correlation_df             - RNA-pathway vs RPPA correlation table
#
# Outputs:
#   results/pathway_score_results.csv
#   results/pathway_coverage_report.csv
#   results/rna_rppa_correlation.csv
#
# Figures:
#   figures/pathway_coverage.png
#   figures/rna_rppa_correlations.png
#
# Note: this script is intended to be sourced by run_analysis.R.

# Validate inputs -------------------------------------------------------------

required_objects <- c(
   "rnaseq_expression",
   "rnaseq_gene_set_membership",
   "survival_data",
   "rppa_univariable_results"
)

check_required_objects(required_objects)
check_has_sample_id("rnaseq_expression")
check_has_columns("survival_data", "sample_id")

# Prepare RNA expression matrix -----------------------------------------------

rna_mat <- rnaseq_expression |>
   tibble::column_to_rownames("sample_id")

# Pathway definitions ---------------------------------------------------------

get_genes_for_query <- function(query_name) {
   rnaseq_gene_set_membership |>
      dplyr::filter(query == !!query_name) |>
      dplyr::pull(gene_symbol) |>
      unique()
}

pathways <- list(
   neutrophil_deg = list(
      genes = get_genes_for_query("neutrophil_degranulation")
   ),
   ecm_deg = list(
      genes = get_genes_for_query("ecm_remodelling")
   ),
   ptk = list(
      genes = unique(c(
         get_genes_for_query("kinase_activity"),
         get_genes_for_query("kinase_signalling")
      ))
   )
)

# Helper functions ------------------------------------------------------------

compute_mean_score <- function(x) {
   rowMeans(x, na.rm = TRUE)
}

compute_pc1_score <- function(x) {
   if (ncol(x) < 2L) {
      return(rep(NA_real_, nrow(x)))
   }
   
   stats::prcomp(x, center = TRUE, scale. = TRUE)$x[, 1]
}

# Compute pathway scores and coverage -----------------------------------------

score_list    <- list()
coverage_list <- list()

for (pw in names(pathways)) {
   gene_set <- pathways[[pw]]$genes
   
   # Map gene symbols -> rnaseq_expression column names (rna_<gene>)
   gene_cols    <- paste0("rna_", gene_set)
   genes_present <- intersect(gene_cols, colnames(rna_mat))
   
   mat_pw <- rna_mat[, genes_present, drop = FALSE]
   
   if (ncol(mat_pw) == 0L) {
      warning("No RNA genes found for pathway: ", pw, call. = FALSE)
      next
   }
   
   # Scores
   score_list[[paste0("score_", pw, "_mean")]] <- compute_mean_score(mat_pw)
   score_list[[paste0("score_", pw)]]          <- compute_pc1_score(mat_pw)
   
   # Coverage (in terms of original gene symbols)
   coverage_list[[pw]] <- tibble::tibble(
      pathway      = pw,
      n_present    = length(genes_present),
      n_total      = length(gene_set),
      pct_coverage = 100 * length(genes_present) / length(gene_set)
   )
}

# Pathway score table ---------------------------------------------------------

pathway_score_results <- tibble::tibble(
   sample_id = rownames(rna_mat)
)

for (nm in names(score_list)) {
   pathway_score_results[[nm]] <- score_list[[nm]]
}

# Coverage summary ------------------------------------------------------------

coverage_df <- dplyr::bind_rows(coverage_list) |>
   dplyr::arrange(pct_coverage)

# Save tables -----------------------------------------------------------------

readr::write_csv(
   pathway_score_results,
   "results/pathway_score_results.csv"
)

readr::write_csv(
   coverage_df,
   "results/pathway_coverage_report.csv"
)

# Coverage visualisation ------------------------------------------------------

grDevices::png(
   filename = "figures/pathway_coverage.png",
   width    = 800,
   height   = 500
)

coverage_df |>
   dplyr::mutate(
      pathway = factor(pathway, levels = pathway)
   ) |>
   ggplot2::ggplot(
      ggplot2::aes(x = pathway, y = pct_coverage)
   ) +
   ggplot2::geom_col() +
   ggplot2::geom_hline(
      yintercept = 50,
      linetype   = "dashed",
      colour     = "red"
   ) +
   ggplot2::coord_flip() +
   ggplot2::ylim(0, 100) +
   ggplot2::labs(
      title = "RNA pathway gene coverage",
      x     = "Pathway",
      y     = "% genes present"
   )

grDevices::dev.off()

# Merge pathway scores into survival_data -------------------------------------

survival_data <- survival_data |>
   dplyr::left_join(
      pathway_score_results,
      by = "sample_id"
   )

# RPPA correlation analysis ---------------------------------------------------

# Top 20 RPPA features by BH-adjusted p then |log HR|
top_rppa_features <- rppa_univariable_results |>
   dplyr::arrange(
      .data$p_adjust_bh,
      dplyr::desc(.data$abs_log_hr)
   ) |>
   dplyr::slice_head(n = 20L) |>
   dplyr::pull(.data$feature)

correlation_results <- list()

for (pw in c("neutrophil_deg", "ecm_deg", "ptk")) {
   pw_col <- paste0("score_", pw)
   
   if (!pw_col %in% names(survival_data)) {
      next
   }
   
   for (feat in top_rppa_features) {
      if (!feat %in% names(survival_data)) {
         next
      }
      
      merged <- survival_data |>
         dplyr::select(
            sample_id,
            dplyr::all_of(pw_col),
            dplyr::all_of(feat)
         ) |>
         tidyr::drop_na()
      
      if (nrow(merged) < 10L) {
         next
      }
      
      r_val <- suppressWarnings(
         stats::cor(
            merged[[pw_col]],
            merged[[feat]],
            method = "spearman",
            use    = "complete.obs"
         )
      )
      
      correlation_results[[length(correlation_results) + 1L]] <- tibble::tibble(
         pathway      = pw,
         rppa_feature = feat,
         spearman_r   = r_val
      )
   }
}

correlation_df <- dplyr::bind_rows(correlation_results) |>
   dplyr::arrange(
      dplyr::desc(abs(.data$spearman_r))
   )

readr::write_csv(
   correlation_df,
   "results/rna_rppa_correlation.csv"
)

# Correlation plot ------------------------------------------------------------

grDevices::png(
   filename = "figures/rna_rppa_correlations.png",
   width    = 900,
   height   = 600
)

correlation_df |>
   ggplot2::ggplot(
      ggplot2::aes(
         x = stats::reorder(rppa_feature, spearman_r),
         y = spearman_r
      )
   ) +
   ggplot2::geom_point() +
   ggplot2::coord_flip() +
   ggplot2::geom_hline(
      yintercept = c(-0.7, 0.7),
      linetype   = "dashed",
      colour     = "red"
   ) +
   ggplot2::labs(
      title = "RNA pathway vs RPPA feature correlations",
      x     = "RPPA feature",
      y     = "Spearman correlation"
   )

grDevices::dev.off()

# Collinearity check ----------------------------------------------------------

high_corr <- correlation_df |>
   dplyr::filter(
      abs(.data$spearman_r) > 0.7
   )

if (nrow(high_corr) > 0L) {
   warning(
      "High RNA–RPPA collinearity detected (|r| > 0.7)",
      call. = FALSE
   )
   print(high_corr)
}

# Final message ---------------------------------------------------------------

message("Pathway scoring complete.")
message("Results saved to results/.")
message("Figures saved to figures/.")
message("Pathway scores appended to survival_data.")