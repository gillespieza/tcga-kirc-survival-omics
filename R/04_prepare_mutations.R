#' @file        04_prepare_mutations.R
#' @title       Prepare Binary Mutation Features
#' @description Derives a sample-by-gene binary mutation feature table from the
#'              raw somatic mutation data. Restricts to a curated set of
#'              established ccRCC driver genes, converts the long mutation table
#'              to a wide binary matrix (1 = mutated, 0 = not mutated), and
#'              ensures all driver genes are represented as columns even if
#'              absent from the data.
#'
#' @details
#'   Input:
#'   \describe{
#'     \item{mutation_data}{Long tibble of somatic mutation calls. Must contain
#'           columns \code{Tumor_Sample_Barcode} and \code{Hugo_Symbol}.}
#'   }
#'
#'   Objects created in the global environment:
#'   \describe{
#'     \item{driver_genes}{Character vector of curated ccRCC driver gene
#'           symbols used to filter the mutation data.}
#'     \item{mutation_features}{Sample-by-gene binary tibble. One row per
#'           tumour sample (\code{sample_id}) and one column per driver gene
#'           (prefixed with \code{mut_}), ordered to match
#'           \code{driver_genes}.}
#'   }
#'
#'   Driver gene selection is based on:
#'   \itemize{
#'     \item Comprehensive molecular characterisation of clear cell RCC
#'           (doi:10.1038/nature12222)
#'     \item Actionable mutations in metastatic RCC
#'           (doi:10.1158/1078-0432.CCR-15-2631)
#'   }
#'
#' @pre  01_load_data.R must be sourced first so that \code{mutation_data} is
#'       available in the global environment.
#'
#' @post \code{mutation_features} is available for use by downstream feature
#'       integration and modelling scripts.
#'
#' @usage source("04_prepare_mutations.R")


# Validate input --------------------------------------------------------------

required_mutation_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol")
missing_mutation_cols  <- setdiff(required_mutation_cols, names(mutation_data))

if (length(missing_mutation_cols) > 0) {
   stop("mutation_data is missing expected column(s): ",
        paste(missing_mutation_cols, collapse = ", "))
}

message("Input mutation data: ", nrow(mutation_data), " rows x ",
        ncol(mutation_data), " cols.")


# Define driver genes ---------------------------------------------------------
# Commonly altered ccRCC genes selected from known kidney cancer biology.
# Sources: doi:10.1038/nature12222 and doi:10.1158/1078-0432.CCR-15-2631

driver_genes <- c(
   "VHL",   # Core ccRCC tumour suppressor; VHL loss drives HIF/hypoxia and angiogenesis biology
   "PBRM1", # Chromatin-remodelling tumour suppressor frequently mutated in ccRCC
   "BAP1",  # Tumour suppressor associated with more aggressive ccRCC and poorer prognosis
   "SETD2", # Chromatin/histone methyltransferase gene altered in ccRCC
   "KDM5C", # Chromatin-regulation gene recurrently altered in ccRCC
   "MTOR",  # Kinase pathway gene; links to PI3K/AKT/mTOR signalling and targeted therapy
   "PTEN",  # Negative regulator of PI3K/AKT signalling; recurrently altered in TCGA ccRCC
   "TSC1",  # mTOR pathway regulator; TSC1/TSC2/MTOR mutations linked to rapalog response
   "TSC2"   # mTOR pathway regulator that functions with TSC1 to suppress mTORC1 signalling
)


# Convert long mutation table to binary sample-by-gene matrix ----------------
# Each driver gene that appears at least once in a sample is counted as
# mutated (1); all other driver genes are coded as not mutated (0).

mutation_long <- mutation_data %>%
   dplyr::transmute(
      sample_id   = .data$Tumor_Sample_Barcode,
      gene_symbol = .data$Hugo_Symbol
   ) %>%
   dplyr::filter(.data$gene_symbol %in% driver_genes) %>%
   dplyr::distinct(.data$sample_id, .data$gene_symbol)

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

mutation_matrix <- table(
   mutation_long$sample_id,
   mutation_long$gene_symbol
)

mutation_features <- as.data.frame.matrix(mutation_matrix) %>%
   tibble::rownames_to_column("sample_id") %>%
   dplyr::mutate(
      dplyr::across(
         -sample_id,
         ~ as.integer(.x > 0)
      )
   ) %>%
   dplyr::rename_with(
      ~ paste0("mut_", .x),
      -sample_id
   )

# Add all-zero columns for driver genes not observed in the mutation data,
# so that the feature table is complete and consistently ordered regardless
# of which genes happen to appear in this dataset.
all_mut_cols     <- paste0("mut_", driver_genes)
missing_mut_cols <- setdiff(all_mut_cols, names(mutation_features))

if (length(missing_mut_cols) > 0) {
   mutation_features[missing_mut_cols] <- 0L
}

# Reorder columns: sample_id first, then driver genes in their defined order.
mutation_features <- mutation_features %>%
   dplyr::select(sample_id, dplyr::all_of(all_mut_cols))


# Validate output -------------------------------------------------------------

stopifnot(
   "mutation_features must be non-empty"      = nrow(mutation_features) > 0,
   "mutation_features must contain sample_id" = "sample_id" %in% names(mutation_features),
   "mutation_features must have all driver gene columns" =
      all(all_mut_cols %in% names(mutation_features))
)


# Summarise mutation frequencies ----------------------------------------------

mutation_summary <- mutation_features %>%
   dplyr::summarise(dplyr::across(dplyr::starts_with("mut_"), sum)) %>%
   tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "gene",
      values_to = "n_mutated"
   ) %>%
   dplyr::mutate(
      gene       = stringr::str_remove(.data$gene, "^mut_"),
      pct_mutated = round(100 * .data$n_mutated / nrow(mutation_features), 1)
   ) %>%
   dplyr::arrange(dplyr::desc(.data$n_mutated))

message(
   "Mutation feature table prepared: ",
   nrow(mutation_features), " samples x ",
   length(all_mut_cols), " driver gene features."
)
print(mutation_summary)