# Pathway score computation & Data-Driven RNA Feature Selection ---------------
#
# Engine 1 (Curated): Computes PC1 scores for all 8 individual MSigDB
#                     signatures.
# Engine 2 (Data-Driven): Conducts an unbiased univariable Cox screening
#                         loop across ALL available RNA-seq genes, isolates
#                         the top 20 by p-value, and computes a summary
#                         PC1 score.
#
# Requires:
#   rnaseq_expression          (sample_id + rna_* columns)
#   rnaseq_gene_set_membership (query, gs_name, gene_symbol)
#   survival_data              (initialized by quick_survival_check.R)
#
# Produces:
#   pathway_score_results      - Combined table containing all 8 curated
#                                scores PLUS the new 'score_rna_datadriven'
#                                feature.
#
# Usage: this script is intended to be sourced by run_analysis.R.


# Validate inputs --------------------------------------------------------------

check_required_objects(
   c("rnaseq_expression", "rnaseq_gene_set_membership", "survival_data")
)


# Prepare data structures ------------------------------------------------------

rna_df         <- rnaseq_expression
unique_queries <- unique(rnaseq_gene_set_membership$query)


# Mathematical engine helpers --------------------------------------------------

compute_mean_score <- function(x) {
   rowMeans(x, na.rm = TRUE)
}

compute_pc1_score <- function(x) {
   if (ncol(x) < 2L) {
      return(rep(NA_real_, nrow(x)))
   }
   stats::prcomp(x, center = TRUE, scale. = TRUE)$x[, 1L]
}


# Engine 1: Curated knowledge-driven pathways (8 signatures) -------------------

message("Running Engine 1: Computing scores for 8 curated MSigDB signatures...")

score_list    <- list()
coverage_list <- list()

for (pw in unique_queries) {
   gene_set <- rnaseq_gene_set_membership |>
      dplyr::filter(.data$query == !!pw) |>
      dplyr::pull(gene_symbol) |>
      unique()
   
   gene_cols     <- paste0("rna_", gene_set)
   genes_present <- intersect(gene_cols, colnames(rna_df))
   
   if (length(genes_present) == 0L) {
      next
   }
   
   # Extract local continuous matrix only when required for numeric functions
   mat_pw <- rna_df |>
      dplyr::select(dplyr::all_of(genes_present)) |>
      as.matrix()
   
   score_list[[paste0("score_", pw, "_mean")]] <- compute_mean_score(mat_pw)
   score_list[[paste0("score_", pw)]]          <- compute_pc1_score(mat_pw)
   
   coverage_list[[pw]] <- tibble::tibble(
      pathway      = pw,
      n_present    = length(genes_present),
      n_total      = length(gene_set),
      pct_coverage = 100 * length(genes_present) / length(gene_set)
   )
}


# Engine 2: Unbiased data-driven feature selection by p-value ------------------

message(
   "Running Engine 2: Executing unbiased genome-wide univariable screening..."
)

# Align survival outcomes with the RNA expression matrix rows cleanly
outcome_df <- survival_data |>
   dplyr::select(sample_id, os_months, os_event) |>
   dplyr::mutate(sample_id = standardise_sample_id(.data$sample_id))

all_rna_cols            <- names(rna_df)[names(rna_df) != "sample_id"]
univariable_rna_results <- list()

# Run a univariable Cox regression for every single available gene column
for (g_col in all_rna_cols) {
   single_gene_df <- rna_df |>
      dplyr::select(sample_id, dplyr::all_of(g_col)) |>
      dplyr::rename(expression = dplyr::all_of(g_col)) |>
      dplyr::inner_join(outcome_df, by = "sample_id") |>
      tidyr::drop_na()
   
   if (nrow(single_gene_df) < 10L) {
      next
   }
   
   # Suppress single-gene separation warnings during high-throughput screening
   fit_gene <- tryCatch(
      suppressWarnings(
         survival::coxph(
            survival::Surv(os_months, os_event) ~ expression,
            data = single_gene_df,
            ties = "efron"
         )
      ),
      error = function(e) NULL
   )
   
   if (!is.null(fit_gene)) {
      coef_summary <- broom::tidy(fit_gene)
      if (nrow(coef_summary) > 0L) {
         univariable_rna_results[[g_col]] <- tibble::tibble(
            gene_feature = g_col,
            p_value      = coef_summary$p.value[1L]
         )
      }
   }
}

# Compile and rank the results by statistical significance
rna_screening_df <- dplyr::bind_rows(univariable_rna_results) |>
   dplyr::arrange(.data$p_value)

# Isolate the top 20 most statistically significant prognostic genes
top_20_datadriven_genes <- rna_screening_df |>
   dplyr::slice_head(n = 20L) |>
   dplyr::pull(gene_feature)

message(
   "Isolated the top 20 most significant data-driven genes. Top 3: ", 
   paste(
      head(stringr::str_remove(top_20_datadriven_genes, "^rna_"), 3L),
      collapse = ", "
   )
)

# Extract the continuous sub-matrix for PCA reduction via tidy select
mat_datadriven <- rna_df |>
   dplyr::select(dplyr::all_of(top_20_datadriven_genes)) |>
   as.matrix()

score_list[["score_rna_datadriven"]] <- compute_pc1_score(mat_datadriven)


# Compile and export cohort data -----------------------------------------------

pathway_score_results <- tibble::tibble(sample_id = rna_df$sample_id)
for (nm in names(score_list)) {
   pathway_score_results[[nm]] <- score_list[[nm]]
}

coverage_df <- dplyr::bind_rows(coverage_list) |>
   dplyr::arrange(.data$pct_coverage)

readr::write_csv(pathway_score_results, "results/pathway_score_results.csv")
readr::write_csv(coverage_df, "results/pathway_coverage_report.csv")
readr::write_csv(
   rna_screening_df, "results/rna_datadriven_screening_pvalues.csv"
)

# Append results back directly into the active console workspace
survival_data <- survival_data |>
   dplyr::select(
      -dplyr::any_of(
         names(pathway_score_results)[names(pathway_score_results) != "sample_id"]
      )
   ) |>
   dplyr::left_join(pathway_score_results, by = "sample_id")

message("Pathway scoring and data-driven RNA selection complete.")