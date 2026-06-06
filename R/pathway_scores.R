# Pathway score computation from RNA-seq expression ---------------------------
#
# Computes per-sample pathway activity scores for all 8 individual MSigDB
# inflammatory and kinase signatures using:
#   1. Mean log2 expression (interpretability baseline)
#   2. PC1 score (primary modelling feature; captures co-expression structure)
#
# Also evaluates:
#   - gene coverage per pathway
#   - per-sample completeness
#   - correlation with selected RPPA features (collinearity screening)
#
# Requires:
#   rnaseq_expression          (prepare_rnaseq.R; sample_id + rna_* columns)
#   rnaseq_gene_set_membership (prepare_rnaseq.R; query, gs_name, gene_symbol)
#   selected_rppa_features     (feature_selection.R)
#   rppa_proteomics            (prepare_rppa.R)
#
# Produces:
#   pathway_score_results      - Per-sample RNA pathway scores for all 8 signatures
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
   "selected_rppa_features",
   "rppa_proteomics"
)

# Enforce defensive checks to ensure previous steps generated required data
check_required_objects(required_objects)
check_has_sample_id("rnaseq_expression")

# Prepare RNA expression matrix -----------------------------------------------

# Move sample_id to row names so the matrix contains purely numeric expression values
rna_mat <- rnaseq_expression |>
   tibble::column_to_rownames("sample_id")

# Dynamically extract all 8 signature queries present in the dataset
unique_queries <- unique(rnaseq_gene_set_membership$query)

# Helper functions ------------------------------------------------------------

# Computes the simple row-wise average of expressions for a given signature block
compute_mean_score <- function(x) {
   rowMeans(x, na.rm = TRUE)
}

# Extracts the first principal component (PC1) to capture maximum co-expression variance
compute_pc1_score <- function(x) {
   if (ncol(x) < 2L) {
      return(rep(NA_real_, nrow(x)))
   }
   stats::prcomp(x, center = TRUE, scale. = TRUE)$x[, 1]
}

# Compute pathway scores and coverage across all 8 signatures -----------------

score_list <- list()
coverage_list <- list()

for (pw in unique_queries) {
   # Isolate the curated gene symbols associated with the active signature query
   gene_set <- rnaseq_gene_set_membership |>
      dplyr::filter(.data$query == !!pw) |>
      dplyr::pull(.data$gene_symbol) |>
      unique()
   
   # Map gene symbols to matching rnaseq_expression column headers (prefixed with rna_)
   gene_cols <- paste0("rna_", gene_set)
   genes_present <- intersect(gene_cols, colnames(rna_mat))
   
   # Extract a narrow sub-matrix containing only the available genes for this signature
   mat_pw <- rna_mat[, genes_present, drop = FALSE]
   
   if (ncol(mat_pw) == 0L) {
      warning("No RNA genes found for signature pathway: ", pw, call. = FALSE)
      next
   }
   
   # Generate the unadjusted baseline average scores
   score_list[[paste0("score_", pw, "_mean")]] <- compute_mean_score(mat_pw)
   
   # Generate the primary PC1 multi-omics modelling features
   score_list[[paste0("score_", pw)]] <- compute_pc1_score(mat_pw)
   
   # Track data completeness and platform gene coverage metrics
   coverage_list[[pw]] <- tibble::tibble(
      pathway      = pw,
      n_present    = length(genes_present),
      n_total      = length(gene_set),
      pct_coverage = 100 * length(genes_present) / length(gene_set)
   )
}

# Pathway score table ---------------------------------------------------------

# Instantiate a wide table indexed by sample barcodes
pathway_score_results <- tibble::tibble(
   sample_id = rownames(rna_mat)
)

# Dynamically bind all calculated means and PC1 scores as distinct columns
for (nm in names(score_list)) {
   pathway_score_results[[nm]] <- score_list[[nm]]
}

# Coverage summary ------------------------------------------------------------

# Combine list elements into a single arranged reporting frame
coverage_df <- dplyr::bind_rows(coverage_list) |>
   dplyr::arrange(.data$pct_coverage)

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

coverage_plot <- coverage_df |>
   dplyr::mutate(
      pathway = factor(.data$pathway, levels = .data$pathway)
   ) |>
   ggplot2::ggplot(
      ggplot2::aes(x = .data$pathway, y = .data$pct_coverage)
   ) +
   ggplot2::geom_col(fill = "#2c3e50") +
   ggplot2::geom_hline(
      yintercept = 50,
      linetype   = "dashed",
      colour     = "#e74c3c"
   ) +
   ggplot2::coord_flip() +
   ggplot2::ylim(0, 100) +
   ggplot2::labs(
      title = "RNA Pathway Gene Coverage (8 Signatures)",
      x     = "MSigDB Signature Query",
      y     = "% Genes Present in Dataset"
   ) +
   ggplot2::theme_classic()

# Routing high-resolution export through our customized defensive graphics engine
save_pipeline_plot(
   plot_object = coverage_plot,
   file_path   = "figures/pathway_coverage.png",
   width       = 1000,
   height      = 600,
   resolution  = 100
)

# Merge pathway scores into survival_data -------------------------------------

# Safe inclusion to append scores when running steps interactively inside the console
if (exists("survival_data")) {
   survival_data <- survival_data |>
      dplyr::left_join(
         pathway_score_results,
         by = "sample_id"
      )
}

# RPPA correlation analysis ---------------------------------------------------

correlation_results <- list()

for (pw in unique_queries) {
   pw_col <- paste0("score_", pw)
   
   if (!pw_col %in% names(pathway_score_results)) {
      next
   }
   
   for (feat in selected_rppa_features) {
      if (!feat %in% names(rppa_proteomics)) {
         next
      }
      
      # Extract clean non-missing sample frames intersected across both layers
      merged <- rppa_proteomics |>
         dplyr::select(dplyr::all_of(c("sample_id", feat))) |>
         dplyr::inner_join(
            pathway_score_results |> dplyr::select(dplyr::all_of(c("sample_id", pw_col))),
            by = "sample_id"
         ) |>
         tidyr::drop_na()
      
      if (nrow(merged) < 10L) {
         next
      }
      
      # Evaluate non-parametric monotonic relationship via Spearman's Rho
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

correlation_df <- dplyr::bind_rows(correlation_results)

if (nrow(correlation_df) > 0) {
   correlation_df <- correlation_df |> 
      dplyr::arrange(dplyr::desc(abs(.data$spearman_r)))
} else {
   correlation_df <- tibble::tibble(
      pathway      = character(), 
      rppa_feature = character(), 
      spearman_r   = numeric()
   )
}

readr::write_csv(
   correlation_df,
   "results/rna_rppa_correlation.csv"
)

# Correlation plot ------------------------------------------------------------

if (nrow(correlation_df) > 0) {
   correlation_plot <- correlation_df |>
      ggplot2::ggplot(
         ggplot2::aes(
            x = stats::reorder(.data$rppa_feature, .data$spearman_r),
            y = .data$spearman_r
         )
      ) +
      ggplot2::geom_point(colour = "#2980b9", size = 2.5) +
      ggplot2::facet_wrap(~pathway, scales = "free_y") +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(
         yintercept = c(-0.7, 0.7),
         linetype   = "dashed",
         colour     = "#e74c3c"
      ) +
      ggplot2::labs(
         title = "Multi-Omics Collinearity: RNA Pathways vs Proteomic Markers",
         x     = "LASSO Selected RPPA Feature",
         y     = "Spearman Correlation Coefficient (r)"
      ) +
      ggplot2::theme_classic()
   
   save_pipeline_plot(
      plot_object = correlation_plot,
      file_path   = "figures/rna_rppa_correlations.png",
      width       = 1400,
      height      = 1000,
      resolution  = 100
   )
}

# Collinearity check ----------------------------------------------------------

high_corr <- correlation_df |>
   dplyr::filter(abs(.data$spearman_r) > 0.7)

if (nrow(high_corr) > 0L) {
   warning(
      "High RNA-RPPA multi-omics collinearity detected (|r| > 0.7)",
      call. = FALSE
   )
   print(high_corr)
}

message("Pathway scoring complete. All 8 distinct signatures processed successfully.")