# Prepare RNA-seq gene expression features ------------------------------------
#
# Prepares a filtered RNA-seq feature table from the raw cBioPortal matrix.
# The script retrieves a small set of MSigDB gene sets, keeps genes present in
# the data, removes duplicates and low-quality genes, reshapes to one row per
# sample, and writes simple output tables for later modelling.
#
# Requires: load_data.R to have been sourced so that rnaseq_data is available in
#            the global environment.
#
# Outputs:
#   results/rnaseq_gene_set_membership.csv
#   results/rnaseq_selected_genes.csv
#   rnaseq_expression (in memory)
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.
#
# Note: the raw data file from cBioPortal is already normalised (meta states
# mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2),
# datatype CONTINUOUS) so no need for DESeq2 or TMM (edgeR) normalization.


# Validate Inputs -------------------------------------------------------------

check_required_objects("rnaseq_data")
check_has_columns("rnaseq_data", "Hugo_Symbol")

n_genes_raw <- nrow(rnaseq_data)
n_samples_raw <- ncol(rnaseq_data) - 1L

message(
  "Input RNA-seq data: ",
  n_genes_raw, " genes x ",
  n_samples_raw, " samples."
)


# Knowledge Based Feature Selection --------------------------------------------
# Gene set definitions
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
  # Perform high-speed in-memory regex filtration over pre-loaded reference
  matched <- msig_reference |>
    dplyr::filter(
      .data$collection == !!collection,
      stringr::str_detect(
        .data$gs_name,
        stringr::regex(pattern, ignore_case = TRUE)
      )
    ) |>
    dplyr::distinct(gs_name, gene_symbol)

  if (nrow(matched) == 0L) {
    warning(
      "[", query_name, "] no gene sets matched pattern: ", pattern,
      call. = FALSE
    )
    return(
      tibble::tibble(
        query       = character(),
        gs_name     = character(),
        gene_symbol = character()
      )
    )
  }

  message(
    "[", query_name, "] ",
    dplyr::n_distinct(matched$gs_name), " gene sets; ",
    dplyr::n_distinct(matched$gene_symbol), " genes"
  )

  matched |>
    dplyr::mutate(query = query_name) |>
    dplyr::select(query, gs_name, gene_symbol)
}


# Retrieve MSigDB Gene Sets ----------------------------------------------------

message("Now running knowledge-based feature selection.")
message("Fetching MSigDB gene sets ...")

# Pre-load required MSigDB collections once to eliminate redundant disk reads
unique_colls <- unique(gene_set_queries$collection)

msig_reference <- purrr::map_dfr(unique_colls, function(coll) {
  msigdbr::msigdbr(species = "Homo sapiens", collection = coll) |>
    dplyr::select(gs_name, gene_symbol) |>
    dplyr::mutate(collection = coll)
})

gene_set_long <- purrr::pmap_dfr(
  gene_set_queries,
  fetch_gene_set
)

if (nrow(gene_set_long) == 0L) {
  stop(
    "No genes retrieved from MSigDB. Check queries and msigdbr version.",
    call. = FALSE
  )
}

rnaseq_gene_set_membership <- gene_set_long |>
  dplyr::distinct(.data$query, .data$gs_name, .data$gene_symbol) |>
  dplyr::arrange(.data$query, .data$gs_name, .data$gene_symbol)


# Gene Universe ----------------------------------------------------------------

candidate_genes <- unique(gene_set_long$gene_symbol)
genes_in_data <- intersect(candidate_genes, rnaseq_data$Hugo_Symbol)

message(
  "Gene set universe: ",
  length(candidate_genes), " candidates; ",
  length(genes_in_data), " present in RNA-seq"
)


# Full-Transcriptome Long Format and QC ---------------------------------------
# Pivots the complete raw expression matrix (all 20,531 genes) to long format,
# standardises sample IDs, log2-transforms, and applies the same QC filters
# used for the pathway-filtered matrix. Used exclusively by Engine 2.

rnaseq_long_full <- rnaseq_data |>
  # Deduplicate on the full gene set by mean expression before pivoting
  dplyr::mutate(
    row_mean = rowMeans(
      as.matrix(dplyr::pick(-Hugo_Symbol)),
      na.rm = TRUE
    )
  ) |>
  dplyr::arrange(dplyr::desc(.data$row_mean)) |>
  dplyr::distinct(.data$Hugo_Symbol, .keep_all = TRUE) |>
  dplyr::select(-row_mean) |>
  tidyr::pivot_longer(
    cols      = -Hugo_Symbol,
    names_to  = "sample_id",
    values_to = "rsem"
  ) |>
  dplyr::mutate(
    sample_id = standardise_sample_id(.data$sample_id),
    rsem      = tidyr::replace_na(.data$rsem, 0),
    log2_expr = log2(.data$rsem + 1)
  )

# Quality Control on the full transcriptome -----------------------------------
# Same variance and missingness thresholds as the pathway-filtered matrix.

gene_qc_full <- rnaseq_long_full |>
  dplyr::group_by(.data$Hugo_Symbol) |>
  dplyr::summarise(
    pct_na      = mean(is.na(.data$log2_expr)),
    variance    = stats::var(.data$log2_expr, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  dplyr::mutate(
    low_variance = is.na(.data$variance) | .data$variance < 1e-6
  )

excluded_genes_full <- gene_qc_full |>
  dplyr::filter(.data$pct_na > 0.20 | .data$low_variance) |>
  dplyr::pull(Hugo_Symbol)

message("Excluded genes (QC, full transcriptome): ", length(excluded_genes_full))

# Full-Transcriptome Wide Matrix (Engine 2 input) -----------------------------

rnaseq_expression_full <- rnaseq_long_full |>
  dplyr::filter(!.data$Hugo_Symbol %in% excluded_genes_full) |>
  dplyr::select(sample_id, Hugo_Symbol, log2_expr) |>
  tidyr::pivot_wider(
    names_from  = "Hugo_Symbol",
    values_from = "log2_expr"
  ) |>
  dplyr::rename_with(~ paste0("rna_", .x), -sample_id)

message(
  "Full-transcriptome matrix for Engine 2: ",
  nrow(rnaseq_expression_full), " samples x ",
  ncol(rnaseq_expression_full) - 1L, " genes"
)


# Filter and Deduplicate (pathway members only) --------------------------------
# From here, processing continues on the pathway-filtered gene universe only.

rnaseq_filtered <- rnaseq_data |>
  dplyr::filter(.data$Hugo_Symbol %in% genes_in_data) |>
  dplyr::mutate(
    row_mean = rowMeans(
      as.matrix(dplyr::pick(-Hugo_Symbol)),
      na.rm = TRUE
    )
  ) |>
  dplyr::arrange(dplyr::desc(.data$row_mean)) |>
  dplyr::distinct(.data$Hugo_Symbol, .keep_all = TRUE) |>
  dplyr::select(-row_mean)

message("Genes after deduplication: ", nrow(rnaseq_filtered))


# Long Format and Transform ----------------------------------------------------

rnaseq_long <- rnaseq_filtered |>
  tidyr::pivot_longer(
    cols      = -Hugo_Symbol,
    names_to  = "sample_id",
    values_to = "rsem"
  ) |>
  dplyr::mutate(
    sample_id = standardise_sample_id(.data$sample_id),
    rsem      = tidyr::replace_na(.data$rsem, 0),
    log2_expr = log2(.data$rsem + 1)
  )


# Quality Control (pathway members) -------------------------------------------

gene_qc <- rnaseq_long |>
  dplyr::group_by(.data$Hugo_Symbol) |>
  dplyr::summarise(
    pct_na      = mean(is.na(.data$log2_expr)),
    variance    = stats::var(.data$log2_expr, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  dplyr::mutate(
    low_variance = is.na(.data$variance) | .data$variance < 1e-6
  )

excluded_genes <- gene_qc |>
  dplyr::filter(.data$pct_na > 0.20 | .data$low_variance) |>
  dplyr::pull(Hugo_Symbol)

message("Excluded genes (QC): ", length(excluded_genes))

# Quality Control --------------------------------------------------------------

gene_qc <- rnaseq_long |>
  dplyr::group_by(.data$Hugo_Symbol) |>
  dplyr::summarise(
    pct_na       = mean(is.na(.data$log2_expr)),
    variance     = stats::var(.data$log2_expr, na.rm = TRUE),
    .groups      = "drop"
  ) |>
  dplyr::mutate(
    low_variance = is.na(.data$variance) | .data$variance < 1e-6
  )

excluded_genes <- gene_qc |>
  dplyr::filter(.data$pct_na > 0.20 | .data$low_variance) |>
  dplyr::pull(Hugo_Symbol)

message("Excluded genes (QC): ", length(excluded_genes))


# Wide Matrix ------------------------------------------------------------------

rnaseq_wide <- rnaseq_long |>
  dplyr::filter(!.data$Hugo_Symbol %in% excluded_genes) |>
  dplyr::select(sample_id, Hugo_Symbol, log2_expr) |>
  tidyr::pivot_wider(
    names_from  = "Hugo_Symbol",
    values_from = "log2_expr"
  )

rnaseq_expression <- rnaseq_wide |>
  dplyr::rename_with(~ paste0("rna_", .x), -sample_id)

rnaseq_gene_union <- setdiff(genes_in_data, excluded_genes)


# Validation -------------------------------------------------------------------

if (nrow(rnaseq_expression) == 0L) {
  stop("rnaseq_expression must be non-empty", call. = FALSE)
}

if (!"sample_id" %in% names(rnaseq_expression)) {
  stop("rnaseq_expression must contain sample_id", call. = FALSE)
}

message(
  "RNA-seq matrix: ",
  nrow(rnaseq_expression), " samples x ",
  ncol(rnaseq_expression) - 1L, " genes"
)


# Outputs ----------------------------------------------------------------------

readr::write_csv(
  rnaseq_gene_set_membership,
  "results/rnaseq_gene_set_membership.csv"
)

selected_genes <- tibble::tibble(
  gene    = rnaseq_gene_union,
  rna_col = paste0("rna_", rnaseq_gene_union)
)

selected_genes <- selected_genes |>
  dplyr::left_join(
    rnaseq_gene_set_membership |>
      dplyr::filter(.data$gene_symbol %in% rnaseq_gene_union) |>
      dplyr::group_by(.data$gene_symbol) |>
      dplyr::summarise(
        queries  = paste(sort(unique(.data$query)), collapse = "; "),
        gs_names = paste(sort(unique(.data$gs_name)), collapse = "; "),
        .groups  = "drop"
      ),
    by = c("gene" = "gene_symbol")
  )

readr::write_csv(
  selected_genes,
  "results/rnaseq_selected_genes.csv"
)

query_summary <- rnaseq_gene_set_membership |>
  dplyr::filter(.data$gene_symbol %in% rnaseq_gene_union) |>
  dplyr::group_by(.data$query) |>
  dplyr::summarise(
    n_gene_sets = dplyr::n_distinct(.data$gs_name),
    n_genes     = dplyr::n_distinct(.data$gene_symbol),
    .groups     = "drop"
  ) |>
  dplyr::arrange(.data$query)

message("Genes available per query (after QC):")
print(query_summary)

message("Total genes retained: ", length(rnaseq_gene_union))
