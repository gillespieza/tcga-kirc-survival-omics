# Prepare RNA-seq gene expression features -----------------------------------
#
# Refactored version:
# - strict MSigDB category filtering
# - deterministic gene set extraction
# - consistent QC ordering
# - explicit gene universe tracking
# - reproducible outputs
#
# Requires: 01_load_data.R (rnaseq_data in environment)
#           msigdbr installed
#
# Outputs:
#   results/rnaseq_gene_set_membership.csv
#   results/rnaseq_selected_genes.csv


# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------

if (!exists("rnaseq_data")) {
   stop("rnaseq_data is missing. Source 01_load_data.R first.")
}

if (!"Hugo_Symbol" %in% names(rnaseq_data)) {
   stop("rnaseq_data is missing expected column: Hugo_Symbol")
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
   stop("Package 'msigdbr' is required. Run install.packages('msigdbr').")
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

message("msigdbr version: ", packageVersion("msigdbr"))

n_genes_raw   <- nrow(rnaseq_data)
n_samples_raw <- ncol(rnaseq_data) - 1L
message("Input RNA-seq data: ", n_genes_raw, " genes x ", n_samples_raw, " samples.")


# ---------------------------------------------------------------------------
# Gene set definitions (strict + non-overlapping)
# ---------------------------------------------------------------------------

gene_set_queries <- list(
   
   neutrophil_degranulation = list(
      collection = "C5",
      pattern = "NEUTROPHIL_DEGRANULATION"
   ),
   
   neutrophil_activation = list(
      collection = "C5",
      pattern = "NEUTROPHIL_(ACTIVATION|MEDIATED_IMMUNITY)"
   ),
   
   inflammatory_response = list(
      collection = "C5",
      pattern = "INFLAMMATORY_RESPONSE"
   ),
   
   hallmark_inflammation = list(
      collection = "H",
      pattern = "INFLAMMATORY_RESPONSE"
   ),
   
   ecm_remodelling = list(
      collection = "C5",
      pattern = paste(
         "EXTRACELLULAR_MATRIX",
         "COLLAGEN",
         "PROTEOLYSIS",
         sep = "|"
      )
   ),
   
   metalloprotease_activity = list(
      collection = "C5",
      pattern = "METALLO(ENDO|EXO)?PEPTIDASE_ACTIVITY"
   ),
   
   kinase_activity = list(
      collection = "C5",
      pattern = "PROTEIN_TYROSINE_KINASE_ACTIVITY"
   ),
   
   kinase_signalling = list(
      collection = "C5",
      pattern = "TYROSINE_KINASE_SIGNALING|RECEPTOR_TYROSINE_KINASE"
   )
)


# ---------------------------------------------------------------------------
# MSigDB fetch (fixed, explicit, reproducible)
# ---------------------------------------------------------------------------

fetch_gene_set <- function(query_name, query_spec) {
   
   coll <- toupper(query_spec$collection)
   pattern <- query_spec$pattern
   
   msig_tbl <- msigdbr::msigdbr(
      species = "Homo sapiens",
      category = coll
   ) %>%
      dplyr::select(gs_name, gene_symbol)
   
   if (nrow(msig_tbl) == 0) {
      stop("Empty MSigDB result for category: ", coll)
   }
   
   matched <- msig_tbl %>%
      dplyr::filter(
         stringr::str_detect(
            gs_name,
            stringr::regex(pattern, ignore_case = TRUE)
         )
      ) %>%
      dplyr::distinct(gs_name, gene_symbol)
   
   if (nrow(matched) == 0) {
      warning("[", query_name, "] no gene sets matched pattern: ", pattern)
      
      return(tibble::tibble(
         query = character(),
         gs_name = character(),
         gene_symbol = character()
      ))
   }
   
   message(
      "[", query_name, "] ",
      dplyr::n_distinct(matched$gs_name), " gene sets; ",
      dplyr::n_distinct(matched$gene_symbol), " genes"
   )
   
   matched %>%
      dplyr::mutate(query = query_name) %>%
      dplyr::select(query, gs_name, gene_symbol)
}


# ---------------------------------------------------------------------------
# Retrieve MSigDB gene sets
# ---------------------------------------------------------------------------

message("\nFetching MSigDB gene sets ...")

gene_set_long <- purrr::map_dfr(
   names(gene_set_queries),
   function(qname) fetch_gene_set(qname, gene_set_queries[[qname]])
)

if (nrow(gene_set_long) == 0) {
   stop("No genes retrieved from MSigDB. Check queries and msigdbr version.")
}


# ---------------------------------------------------------------------------
# Gene set membership tables
# ---------------------------------------------------------------------------

rnaseq_gene_set_membership <- gene_set_long %>%
   dplyr::distinct(query, gs_name, gene_symbol) %>%
   dplyr::arrange(query, gs_name, gene_symbol)

rnaseq_gene_sets <- gene_set_long %>%
   dplyr::group_by(query) %>%
   dplyr::summarise(
      genes = list(unique(gene_symbol)),
      .groups = "drop"
   ) %>%
   tibble::deframe()


# ---------------------------------------------------------------------------
# Gene universe (strict + reproducible)
# ---------------------------------------------------------------------------

candidate_genes <- unique(gene_set_long$gene_symbol)
genes_in_data   <- intersect(candidate_genes, rnaseq_data$Hugo_Symbol)

message(
   "Gene set universe: ", length(candidate_genes),
   " candidates; ", length(genes_in_data), " present in RNA-seq"
)


# ---------------------------------------------------------------------------
# Filter RNA-seq matrix + deduplicate genes
# ---------------------------------------------------------------------------

rnaseq_filtered <- rnaseq_data %>%
   dplyr::filter(Hugo_Symbol %in% genes_in_data) %>%
   dplyr::mutate(
      row_mean = rowMeans(dplyr::select(., -Hugo_Symbol), na.rm = TRUE)
   ) %>%
   dplyr::arrange(dplyr::desc(row_mean)) %>%
   dplyr::distinct(Hugo_Symbol, .keep_all = TRUE) %>%
   dplyr::select(-row_mean)

message("Genes after deduplication: ", nrow(rnaseq_filtered))


# ---------------------------------------------------------------------------
# Long format + transform
# ---------------------------------------------------------------------------

rnaseq_long <- rnaseq_filtered %>%
   tidyr::pivot_longer(
      cols = -Hugo_Symbol,
      names_to = "sample_id",
      values_to = "rsem"
   ) %>%
   dplyr::mutate(
      rsem = tidyr::replace_na(rsem, 0),
      log2_expr = log2(rsem + 1)
   )


# ---------------------------------------------------------------------------
# QC (global gene-level)
# ---------------------------------------------------------------------------

gene_qc <- rnaseq_long %>%
   dplyr::group_by(Hugo_Symbol) %>%
   dplyr::summarise(
      pct_na = mean(is.na(log2_expr)),
      variance = stats::var(log2_expr, na.rm = TRUE),
      .groups = "drop"
   ) %>%
   dplyr::mutate(
      low_variance = is.na(variance) | variance < 1e-6
   )

excluded_genes <- gene_qc %>%
   dplyr::filter(pct_na > 0.20 | low_variance) %>%
   dplyr::pull(Hugo_Symbol)

message("Excluded genes (QC): ", length(excluded_genes))


# ---------------------------------------------------------------------------
# Wide matrix construction
# ---------------------------------------------------------------------------

rnaseq_wide <- rnaseq_long %>%
   dplyr::filter(!Hugo_Symbol %in% excluded_genes) %>%
   dplyr::select(sample_id, Hugo_Symbol, log2_expr) %>%
   tidyr::pivot_wider(
      names_from = Hugo_Symbol,
      values_from = log2_expr
   )

rnaseq_expression <- rnaseq_wide %>%
   dplyr::rename_with(~ paste0("rna_", .x), -sample_id)


# ---------------------------------------------------------------------------
# Final gene universe (post-QC)
# ---------------------------------------------------------------------------

rnaseq_gene_union <- setdiff(genes_in_data, excluded_genes)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

stopifnot(
   nrow(rnaseq_expression) > 0,
   "sample_id" %in% names(rnaseq_expression)
)

message(
   "RNA-seq matrix: ",
   nrow(rnaseq_expression), " samples x ",
   ncol(rnaseq_expression) - 1, " genes"
)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

readr::write_csv(
   rnaseq_gene_set_membership,
   "results/rnaseq_gene_set_membership.csv"
)

readr::write_csv(
   tibble::tibble(
      gene = rnaseq_gene_union,
      rna_col = paste0("rna_", rnaseq_gene_union)
   ) %>%
      dplyr::left_join(
         rnaseq_gene_set_membership %>%
            dplyr::filter(gene_symbol %in% rnaseq_gene_union) %>%
            dplyr::group_by(gene_symbol) %>%
            dplyr::summarise(
               queries = paste(sort(unique(query)), collapse = "; "),
               gs_names = paste(sort(unique(gs_name)), collapse = "; "),
               .groups = "drop"
            ),
         by = c("gene" = "gene_symbol")
      ),
   "results/rnaseq_selected_genes.csv"
)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

query_summary <- rnaseq_gene_set_membership %>%
   dplyr::filter(gene_symbol %in% rnaseq_gene_union) %>%
   dplyr::group_by(query) %>%
   dplyr::summarise(
      n_gene_sets = dplyr::n_distinct(gs_name),
      n_genes = dplyr::n_distinct(gene_symbol),
      .groups = "drop"
   ) %>%
   dplyr::arrange(query)

message("\nGenes available per query (after QC):")
print(query_summary)

message("\nTotal genes retained: ", length(rnaseq_gene_union))