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
#   rnaseq_expression (from 03b_prepare_rnaseq.R)
#   survival_data
#   rppa_univariable_results (for correlation check)
#
# Outputs:
#   /results/pathway_score_results.csv
#   /results/pathway_coverage_report.csv
#   /results/rna_rppa_correlation.csv
#
# Figures:
#   /figures/pathway_coverage.png
#   /figures/rna_rppa_correlations.png


# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)


# ------------------------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------------------------

required_objects <- c(
   "rnaseq_expression",
   "survival_data",
   "rppa_univariable_results",
   "rnaseq_gene_sets"
)

missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
   stop("Missing objects: ", paste(missing_objects, collapse = ", "))
}

stopifnot("sample_id" %in% names(rnaseq_expression))


# ------------------------------------------------------------------------------
# Prepare RNA expression matrix
# ------------------------------------------------------------------------------

rna_mat <- rnaseq_expression %>%
   tibble::column_to_rownames("sample_id")


# ------------------------------------------------------------------------------
# Pathway definitions
# ------------------------------------------------------------------------------

pathways <- list(
   
   neutrophil_deg = list(
      genes = rnaseq_gene_sets$neutrophil_degranulation
   ),
   
   ecm_deg = list(
      genes = rnaseq_gene_sets$ecm_remodelling
   ),
   
   ptk = list(
      genes = unique(c(
         rnaseq_gene_sets$kinase_activity,
         rnaseq_gene_sets$kinase_signalling
      ))
   )
)


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

compute_mean_score <- function(x) {
   rowMeans(x, na.rm = TRUE)
}

compute_pc1_score <- function(x) {
   
   if (ncol(x) < 2) {
      return(rep(NA_real_, nrow(x)))
   }
   
   stats::prcomp(x, center = TRUE, scale. = TRUE)$x[, 1]
}


# ------------------------------------------------------------------------------
# Compute pathway scores + coverage
# ------------------------------------------------------------------------------

score_list <- list()
coverage_list <- list()

for (pw in names(pathways)) {
   
   gene_set <- pathways[[pw]]$genes
   genes_present <- intersect(gene_set, colnames(rna_mat))
   
   mat_pw <- rna_mat[, genes_present, drop = FALSE]
   
   if (ncol(mat_pw) == 0) {
      warning("No genes found for pathway: ", pw)
      next
   }
   
   # ---- Scores ----
   score_list[[paste0("score_", pw, "_mean")]] <- compute_mean_score(mat_pw)
   score_list[[paste0("score_", pw)]] <- compute_pc1_score(mat_pw)
   
   # ---- Coverage ----
   coverage_list[[pw]] <- tibble::tibble(
      pathway = pw,
      n_present = length(genes_present),
      n_total = length(gene_set),
      pct_coverage = 100 * length(genes_present) / length(gene_set)
   )
}


# ------------------------------------------------------------------------------
# Pathway score table
# ------------------------------------------------------------------------------

pathway_score_results <- tibble::tibble(
   sample_id = rownames(rna_mat)
)

for (nm in names(score_list)) {
   pathway_score_results[[nm]] <- score_list[[nm]]
}


# ------------------------------------------------------------------------------
# Coverage summary
# ------------------------------------------------------------------------------

coverage_df <- dplyr::bind_rows(coverage_list) %>%
   dplyr::arrange(pct_coverage)


# ------------------------------------------------------------------------------
# Save tables
# ------------------------------------------------------------------------------

readr::write_csv(
   pathway_score_results,
   "results/pathway_score_results.csv"
)

readr::write_csv(
   coverage_df,
   "results/pathway_coverage_report.csv"
)


# ------------------------------------------------------------------------------
# Coverage visualisation
# ------------------------------------------------------------------------------

png("figures/pathway_coverage.png", width = 800, height = 500)

coverage_df %>%
   dplyr::mutate(pathway = factor(pathway, levels = pathway)) %>%
   ggplot2::ggplot(ggplot2::aes(x = pathway, y = pct_coverage)) +
   ggplot2::geom_col() +
   ggplot2::geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
   ggplot2::coord_flip() +
   ggplot2::ylim(0, 100) +
   ggplot2::labs(
      title = "RNA Pathway Gene Coverage",
      x = "Pathway",
      y = "% genes present"
   )

dev.off()


# ------------------------------------------------------------------------------
# Merge with survival data
# ------------------------------------------------------------------------------

survival_data <- survival_data %>%
   dplyr::left_join(pathway_score_results, by = "sample_id")


# ------------------------------------------------------------------------------
# RPPA correlation analysis
# ------------------------------------------------------------------------------

top_rppa_features <- rppa_univariable_results %>%
   dplyr::arrange(desc(abs(coef))) %>%
   dplyr::slice_head(n = 20) %>%
   dplyr::pull(feature)


correlation_results <- list()

for (pw in c("neutrophil_deg", "ecm_deg", "ptk")) {
   
   pw_col <- paste0("score_", pw)
   
   if (!pw_col %in% names(pathway_score_results)) next
   
   for (feat in top_rppa_features) {
      
      if (!feat %in% colnames(rnaseq_expression)) next
      
      merged <- survival_data %>%
         dplyr::select(sample_id, dplyr::all_of(pw_col)) %>%
         dplyr::inner_join(
            rnaseq_expression %>%
               dplyr::select(sample_id, dplyr::all_of(feat)),
            by = "sample_id"
         )
      
      if (nrow(merged) < 10) next
      
      r_val <- suppressWarnings(
         cor(merged[[pw_col]], merged[[feat]], method = "spearman", use = "complete.obs")
      )
      
      correlation_results[[length(correlation_results) + 1]] <- tibble::tibble(
         pathway = pw,
         rppa_feature = feat,
         spearman_r = r_val
      )
   }
}


correlation_df <- dplyr::bind_rows(correlation_results) %>%
   dplyr::arrange(desc(abs(spearman_r)))


readr::write_csv(
   correlation_df,
   "results/rna_rppa_correlation.csv"
)


# ------------------------------------------------------------------------------
# Correlation plot
# ------------------------------------------------------------------------------

png("figures/rna_rppa_correlations.png", width = 900, height = 600)

correlation_df %>%
   ggplot2::ggplot(ggplot2::aes(x = reorder(rppa_feature, spearman_r), y = spearman_r)) +
   ggplot2::geom_point() +
   ggplot2::coord_flip() +
   ggplot2::geom_hline(yintercept = c(-0.7, 0.7), linetype = "dashed", colour = "red") +
   ggplot2::labs(
      title = "RNA Pathway vs RPPA Feature Correlations",
      x = "RPPA feature",
      y = "Spearman correlation"
   )

dev.off()


# ------------------------------------------------------------------------------
# Collinearity check
# ------------------------------------------------------------------------------

high_corr <- correlation_df %>%
   dplyr::filter(abs(spearman_r) > 0.7)

if (nrow(high_corr) > 0) {
   warning("High RNA–RPPA collinearity detected (|r| > 0.7)")
   print(high_corr)
}


# ------------------------------------------------------------------------------
# Final message
# ------------------------------------------------------------------------------

message("Pathway scoring complete")
message("Results saved to /results")
message("Figures saved to /figures")
message("Pathway scores appended to survival_data")