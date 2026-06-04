# Prepare binary mutation features --------------------------------------------
#
# Derives a sample-by-gene binary mutation feature table from the raw somatic
# mutation data. Restricts to a curated set of established ccRCC driver genes,
# converts the long mutation table to a wide binary matrix (1 = mutated,
# 0 = not mutated), and ensures all driver genes are represented as columns
# even if absent from the data.
#
# Requires: load_data.R to have been sourced so that mutation_data is available
#           in the global environment.
#
# Produces:
#   driver_genes      - Character vector of curated ccRCC driver gene symbols.
#   mutation_features - Sample-by-gene binary tibble. One row per tumour
#                       sample (sample_id) and one column per driver gene
#                       (prefixed with mut_), ordered to match driver_genes.
#
# Driver gene selection is based on:
#   - Comprehensive molecular characterisation of clear cell RCC
#     (doi:10.1038/nature12222)
#   - Actionable mutations in metastatic RCC
#     (doi:10.1158/1078-0432.CCR-15-2631)
#
# Usage: this script is intended to be sourced by run_analysis.R as part of
#        the full pipeline, not run directly.


# Validate input --------------------------------------------------------------

check_required_objects("mutation_data")

check_has_columns(
  "mutation_data",
  c("Tumor_Sample_Barcode", "Hugo_Symbol")
)

abort_if_false(
  nrow(mutation_data) > 0,
  "mutation_data has zero rows after loading. Check data_mutations.txt."
)

message(
  "Input mutation data: ",
  nrow(mutation_data), " rows x ",
  ncol(mutation_data), " cols."
)


# Define driver genes ---------------------------------------------------------
# Commonly altered ccRCC genes selected from known kidney cancer biology.
# Sources: doi:10.1038/nature12222 and doi:10.1158/1078-0432.CCR-15-2631

driver_genes <- c(
  # Core ccRCC tumour suppressor; loss drives HIF/hypoxia, angiogenesis biology
  "VHL",

  # Chromatin-remodelling tumour suppressor frequently mutated in ccRCC
  "PBRM1",

  # Tumour suppressor associated with more aggressive ccRCC & poorer prognosis
  "BAP1",

  # Chromatin/histone methyltransferase gene altered in ccRCC
  "SETD2",

  # Chromatin-regulation gene recurrently altered in ccRCC
  "KDM5C",

  # Kinase pathway gene; links to PI3K/AKT/mTOR signalling and targeted therapy
  "MTOR",

  # Negative regulator of PI3K/AKT signalling; recurrently altered in ccRCC
  "PTEN",

  # mTOR pathway regulator; TSC1/TSC2/MTOR mutations linked to rapalog response
  "TSC1",

  # mTOR pathway regulator, functions with TSC1 to suppress mTORC1 signalling
  "TSC2"
)

# Convert long mutation table to binary sample-by-gene matrix ----------------
# Each driver gene that appears at least once in a sample is counted as
# mutated (1); all other driver genes are coded as not mutated (0).
mutation_long <- mutation_data |>
  dplyr::transmute(
    sample_id   = .data$Tumor_Sample_Barcode,
    gene_symbol = .data$Hugo_Symbol
  ) |>
  dplyr::filter(.data$gene_symbol %in% driver_genes) |>
  dplyr::distinct(sample_id, gene_symbol)

abort_if_false(
  nrow(mutation_long) > 0,
  paste(
    "No mutations in the selected driver genes were found in mutation_data.",
    "Check that the study and driver_genes are correct."
  )
)

# Warn early if any driver genes are entirely absent from the mutation data,
# before the matrix step silently drops them.
unobserved_genes <- setdiff(driver_genes, unique(mutation_long$gene_symbol))

if (length(unobserved_genes) > 0) {
  message(
    "The following driver gene(s) have no mutations in this dataset and ",
    "will be added as all-zero columns: ",
    paste(unobserved_genes, collapse = ", ")
  )
}

# Count at least one mutation per sample/gene, pivot to wide, and binarise.
mutation_features <- mutation_long |>
  dplyr::count(
    .data$sample_id,
    .data$gene_symbol,
    name = "mutated"
  ) |>
  tidyr::pivot_wider(
    names_from  = "gene_symbol",
    values_from = "mutated",
    values_fill = 0L
  ) |>
  dplyr::mutate(
    dplyr::across(
      -dplyr::all_of("sample_id"),
      ~ as.integer(.x > 0)
    )
  ) |>
  dplyr::rename_with(
    ~ paste0("mut_", .x),
    -dplyr::all_of("sample_id")
  )

# Add all-zero columns for driver genes not observed in the mutation data,
# so that the feature table is complete and consistently ordered regardless
# of which genes happen to appear in this dataset.

all_mut_cols <- paste0("mut_", driver_genes)
missing_mut_cols <- setdiff(all_mut_cols, names(mutation_features))

if (length(missing_mut_cols) > 0) {
  mutation_features[missing_mut_cols] <- 0L
}

mutation_features <- mutation_features |>
  dplyr::select(
    dplyr::all_of("sample_id"),
    dplyr::all_of(all_mut_cols)
  )


# Validate output -------------------------------------------------------------

abort_if_false(
  nrow(mutation_features) > 0,
  "mutation_features must be non-empty."
)

abort_if_false(
  "sample_id" %in% names(mutation_features),
  "mutation_features must contain sample_id."
)

abort_if_false(
  all(all_mut_cols %in% names(mutation_features)),
  "mutation_features must have all driver gene columns."
)


# Summarise mutation frequencies ----------------------------------------------
# Mutation counts are reported using the shared helper from utils_validation.R.

mutation_summary <- summarise_mutations(mutation_features)

message(
  "Mutation feature table prepared: ",
  nrow(mutation_features), " samples x ",
  length(all_mut_cols), " driver gene features."
)

print(mutation_summary)
