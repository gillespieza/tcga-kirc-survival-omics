# Prepare RNA-seq gene expression features -----------------------------------
#
# Prepares a filtered RNA-seq feature table from the raw cBioPortal matrix.
# The script retrieves a small set of MSigDB gene sets, keeps genes present in
# the data, removes duplicates and low-quality genes, reshapes to one row per
# sample, and writes simple output tables for later modelling.
#
# Requires: load_data.R to have been sourced so that rnaseq_data is available in
#           the global environment.
#           Package msigdbr must be installed.
#
# Outputs:
#   results/rnaseq_gene_set_membership.csv
#   results/rnaseq_selected_genes.csv
#   rnaseq_expression (in memory)
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.

# Validate inputs -------------------------------------------------------------
abort_if_false <- function(condition, message_text) {
   if (!condition) {
      stop(message_text, call. = FALSE)
   }
}
abort_if_false(
   exists("rnaseq_data"),
   "rnaseq_data is missing. Source 01_load_data.R first."
)

abort_if_false(
   "Hugo_Symbol" %in% names(rnaseq_data),
   "rnaseq_data is missing expected column: Hugo_Symbol"
)

abort_if_false(
   requireNamespace("msigdbr", quietly = TRUE),
   "Package \"msigdbr\" is required. Run install.packages(\"msigdbr\")."
)

# Validate msigdbr version -------------------------------------------------------------
message("msigdbr version: ", as.character(utils::packageVersion("msigdbr")))

n_genes_raw   <- nrow(rnaseq_data)
n_samples_raw <- ncol(rnaseq_data) - 1L

message(
   "Input RNA-seq data: ",
   n_genes_raw, " genes x ",
   n_samples_raw, " samples."
)


# Gene set definitions --------------------------------------------------------
# These queries are intentionally small and targeted so the downstream gene
# universe stays interpretable.

gene_set_queries <- tibble::tibble(
   query = c(
      "neutrophil_degranulation",
      "neutrophil_activation",
      "inflammatory_response",
      "hallmark_inflammation",
      "ecm_remodelling",
      "metalloprotease_activity",
      "kinase_activity",
      "kinase_signalling"
   ),
   collection = c("C5", "C5", "C5", "H", "C5", "C5", "C5", "C5"),
   pattern = c(
      "NEUTROPHIL_DEGRANULATION",
      "NEUTROPHIL_(ACTIVATION|MEDIATED_IMMUNITY)",
      "INFLAMMATORY_RESPONSE",
      "INFLAMMATORY_RESPONSE",
      "EXTRACELLULAR_MATRIX|COLLAGEN|PROTEOLYSIS",
      "METALLO(ENDO|EXO)?PEPTIDASE_ACTIVITY",
      "PROTEIN_TYROSINE_KINASE_ACTIVITY",
      "TYROSINE_KINASE_SIGNALING|RECEPTOR_TYROSINE_KINASE"
   )
)

fetch_gene_set <- function(query_name, collection, pattern) {
   msig_tbl <- msigdbr::msigdbr(
      species = "Homo sapiens",
      category = collection
   ) |>
      dplyr::select(gs_name, gene_symbol)
   
   if (nrow(msig_tbl) == 0) {
      stop(
         "Empty MSigDB result for category: ",
         collection,
         call. = FALSE
      )
   }
   
   matched <- msig_tbl |>
      dplyr::filter(
         stringr::str_detect(
            gs_name,
            stringr::regex(pattern, ignore_case = TRUE)
         )
      ) |>
      dplyr::distinct(gs_name, gene_symbol)
   
   if (nrow(matched) == 0) {
      warning(
         "[",
         query_name,
         "] no gene sets matched pattern: ",
         pattern,
         call. = FALSE
      )
      return(
         tibble::tibble(
            query = character(),
            gs_name = character(),
            gene_symbol = character()
         )
      )
   }
   
   message(
      "[",
      query_name,
      "] ",
      dplyr::n_distinct(matched$gs_name),
      " gene sets; ",
      dplyr::n_distinct(matched$gene_symbol),
      " genes"
   )
   
   matched |>
      dplyr::mutate(query = query_name) |>
      dplyr::select(query, gs_name, gene_symbol)
}

# Retrieve MSigDB gene sets ---------------------------------------------------

message("Fetching MSigDB gene sets ...")

gene_set_long <- purrr::pmap_dfr(
   gene_set_queries,
   fetch_gene_set
)

if (nrow(gene_set_long) == 0) {
   stop(
      "No genes retrieved from MSigDB. Check queries and msigdbr version.",
      call. = FALSE
   )
}

rnaseq_gene_set_membership <- gene_set_long |>
   dplyr::distinct(query, gs_name, gene_symbol) |>
   dplyr::arrange(query, gs_name, gene_symbol)

# Gene universe ---------------------------------------------------------------

candidate_genes <- unique(gene_set_long$gene_symbol)
genes_in_data <- intersect(candidate_genes, rnaseq_data$Hugo_Symbol)

message(
   "Gene set universe: ",
   length(candidate_genes),
   " candidates; ",
   length(genes_in_data),
   " present in RNA-seq"
)

# Filter and deduplicate ------------------------------------------------------

rnaseq_filtered <- rnaseq_data |>
   dplyr::filter(Hugo_Symbol %in% genes_in_data) |>
   dplyr::mutate(
      row_mean = rowMeans(
         as.matrix(dplyr::pick(-Hugo_Symbol)),
         na.rm = TRUE
      )
   ) |>
   dplyr::arrange(dplyr::desc(row_mean)) |>
   dplyr::distinct(Hugo_Symbol, .keep_all = TRUE) |>
   dplyr::select(-row_mean)

message("Genes after deduplication: ", nrow(rnaseq_filtered))

# Long format and transform ---------------------------------------------------

rnaseq_long <- rnaseq_filtered |>
   tidyr::pivot_longer(
      cols = -Hugo_Symbol,
      names_to = "sample_id",
      values_to = "rsem"
   ) |>
   dplyr::mutate(
      rsem = tidyr::replace_na(rsem, 0),
      log2_expr = log2(rsem + 1)
   )
# QC --------------------------------------------------------------------------

gene_qc <- rnaseq_long |>
   dplyr::group_by(Hugo_Symbol) |>
   dplyr::summarise(
      pct_na = mean(is.na(log2_expr)),
      variance = stats::var(log2_expr, na.rm = TRUE),
      .groups = "drop"
   ) |>
   dplyr::mutate(
      low_variance = is.na(variance) | variance < 1e-6
   )

excluded_genes <- gene_qc |>
   dplyr::filter(pct_na > 0.20 | low_variance) |>
   dplyr::pull(Hugo_Symbol)

message("Excluded genes (QC): ", length(excluded_genes))

# Wide matrix -----------------------------------------------------------------

rnaseq_wide <- rnaseq_long |>
   dplyr::filter(!Hugo_Symbol %in% excluded_genes) |>
   dplyr::select(sample_id, Hugo_Symbol, log2_expr) |>
   tidyr::pivot_wider(
      names_from = Hugo_Symbol,
      values_from = log2_expr
   )

rnaseq_expression <- rnaseq_wide |>
   dplyr::rename_with(~ paste0("rna_", .x), -sample_id)

rnaseq_gene_union <- setdiff(genes_in_data, excluded_genes)

# Validation -----------------------------------------------------------------

if (nrow(rnaseq_expression) == 0) {
   stop("rnaseq_expression must be non-empty", call. = FALSE)
}

if (!"sample_id" %in% names(rnaseq_expression)) {
   stop("rnaseq_expression must contain sample_id", call. = FALSE)
}

message(
   "RNA-seq matrix: ",
   nrow(rnaseq_expression),
   " samples x ",
   ncol(rnaseq_expression) - 1L,
   " genes"
)

# Outputs --------------------------------------------------------------------

readr::write_csv(
   rnaseq_gene_set_membership,
   "results/rnaseq_gene_set_membership.csv"
)

selected_genes <- tibble::tibble(
   gene = rnaseq_gene_union,
   rna_col = paste0("rna_", rnaseq_gene_union)
)

selected_genes <- selected_genes |>
   dplyr::left_join(
      rnaseq_gene_set_membership |>
         dplyr::filter(gene_symbol %in% rnaseq_gene_union) |>
         dplyr::group_by(gene_symbol) |>
         dplyr::summarise(
            queries = paste(sort(unique(query)), collapse = "; "),
            gs_names = paste(sort(unique(gs_name)), collapse = "; "),
            .groups = "drop"
         ),
      by = c("gene" = "gene_symbol")
   )

readr::write_csv(
   selected_genes,
   "results/rnaseq_selected_genes.csv"
)

query_summary <- rnaseq_gene_set_membership |>
   dplyr::filter(gene_symbol %in% rnaseq_gene_union) |>
   dplyr::group_by(query) |>
   dplyr::summarise(
      n_gene_sets = dplyr::n_distinct(gs_name),
      n_genes = dplyr::n_distinct(gene_symbol),
      .groups = "drop"
   ) |>
   dplyr::arrange(query)

message("\nGenes available per query (after QC):")
print(query_summary)

message("\nTotal genes retained: ", length(rnaseq_gene_union))