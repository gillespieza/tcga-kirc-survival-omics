# Utility validation helpers --------------------------------------------------
#
# Small shared helpers for common validation patterns and data streamlining
# used across the pipeline.
#
# Note: this script is intended to be sourced by run_analysis.R at the very
# beginning of the pipeline so these functions are available globally.


# Defensive Assertions --------------------------------------------------------

# Halts execution with a custom message if a required logical condition is FALSE
abort_if_false <- function(condition, message_text) {
   if (!condition) {
      stop(
         message_text,
         call. = FALSE # Suppresses internal function call in the error message
      )
   }
}

# Checks if a vector of object names exists in the specified environment
check_required_objects <- function(object_names, env = .GlobalEnv) {
   missing <- object_names[
      # vapply checks each object name and returns a true/false logical vector
      !vapply(object_names, exists, logical(1), envir = env, inherits = TRUE)
   ]
   
   # If any objects are missing, stop the pipeline with an informative error
   abort_if_false(
      length(missing) == 0L,
      paste(
         "The following required object(s) are missing from the environment:",
         paste(missing, collapse = ", "),
         "\nSource the appropriate setup/prepare script(s) via run_analysis.R."
      )
   )
}

# Validates that a specific data frame contains all required columns
check_has_columns <- function(object_name, columns, env = .GlobalEnv) {
   # Retrieve the actual data frame object using its string name
   obj <- get(object_name, envir = env)
   
   # Identify which required columns are missing from the data frame names
   missing_cols <- setdiff(columns, names(obj))
   
   abort_if_false(
      length(missing_cols) == 0L,
      paste0(
         object_name, " is missing required column(s): ",
         paste(missing_cols, collapse = ", ")
      )
   )
}

# Dedicated wrapper to verify that 'sample_id' is present in a dataset.
# Delegates to check_has_columns()
check_has_sample_id <- function(object_name, env = .GlobalEnv) {
   check_has_columns(
      object_name,
      "sample_id",
      env = env
   )
}


# Pipeline DRY Optimisations --------------------------------------------------

# Summarise Mutation Frequencies
# Takes a data frame containing binary mutation columns prefixed with "mut_"
# and returns a sorted, long-format summary table with sample percentages,
# accounting for missing profiles via non-missing denominators.
summarise_mutations <- function(data) {
   summary_table <- data |>
      # Isolate mutation columns cleanly
      dplyr::select(dplyr::starts_with("mut_", ignore.case = FALSE)) |>
      # Pivot into a long single vertical stream for uniform grouping
      tidyr::pivot_longer(
         cols      = dplyr::everything(),
         names_to  = "gene",
         values_to = "status"
      ) |>
      # Group by individual gene to calculate specific profile statistics
      dplyr::group_by(.data$gene) |>
      dplyr::summarise(
         n_mutated   = sum(.data$status == 1L, na.rm = TRUE),
         n_profiled  = sum(!is.na(.data$status)),
         pct_mutated = round(100 * sum(
            .data$status == 1L,
            na.rm = TRUE) / sum(!is.na(.data$status)),
            1
         ),
         .groups     = "drop"
      ) |>
      # Strip the prefix for clear report-ready tables
      dplyr::mutate(
         gene = stringr::str_remove(.data$gene, "^mut_")
      ) |>
      # Sort with the most frequently mutated driver genes at the top
      dplyr::arrange(
         dplyr::desc(.data$n_mutated)
      )
   
   summary_table
}

# Enforce Complete Cases across Specific Variables
# Subsets a data frame to chosen variables, drops any rows containing NA values,
# and cleanly logs the number of dropped samples to the console.
enforce_complete_cases <- function(data, variables) {
   # Keep only the subset of columns required for modeling
   filtered_data <- data |>
      dplyr::select(dplyr::all_of(variables)) |>
      tidyr::drop_na()
   
   # Calculate exactly how many rows were dropped by drop_na()
   rows_dropped <- nrow(data) - nrow(filtered_data)
   
   if (rows_dropped > 0L) {
      message(
         rows_dropped,
         " sample(s) dropped to establish a unified complete-case dataset."
      )
   }
   
   filtered_data
}


# Standardised Graphics Export -------------------------------------------------

# Formats and saves high-resolution graphics files to disk defensively
save_pipeline_plot <- function(
      plot_object,
      file_path,
      width,
      height,
      resolution = 72
   ) {
   # Open the PNG file device to capture the high-resolution plot on disk
   grDevices::png(
      filename = file_path,
      width    = width,
      height   = height,
      res      = resolution
   )
   
   # Handle survminer object wrapper separation safely
   if (inherits(plot_object, "ggsurvplot")) {
      print(plot_object$plot)
   } else {
      print(plot_object)
   }
   
   # Close the file device to finish writing and save the PNG file
   grDevices::dev.off()
   
   # Print the plot a second time to send it to the active R graphics canvas
   if (inherits(plot_object, "ggsurvplot")) {
      print(plot_object)
   } else {
      print(plot_object)
   }
}

# Standardise TCGA sample barcodes to a uniform 15-character hyphenated format
# (e.g., converts 'TCGA.KL.8323.01A' or 'TCGA-KL-8323-01A-11D' to 
# 'TCGA-KL-8323-01')
standardise_sample_id <- function(ids) {
   ids_clean <- as.character(ids)
   ids_clean <- stringr::str_replace_all(ids_clean, "\\.", "-")
   ids_clean <- stringr::str_sub(ids_clean, start = 1L, end = 15L)
   ids_clean
}