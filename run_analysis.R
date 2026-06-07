# Run the full TCGA KIRC survival analysis pipeline ---------------------------
#
# Overview:
#   An integrated, multi-omics machine learning and survival modeling pipeline
#   for Kidney Renal Clear Cell Carcinoma (TCGA-KIRC). The workflow ingests
#   curated clinical metadata, RPPA proteomics, RNA-seq expression, and
#   copy-number alterations (CNA). It executes automated data
#   cleansing, clinical grade pooling, multi-omics sample cross-matching, and
#   feature selection (multivariable LASSO for proteins and a comparative dual-
#   engine transcriptomic profiling of 8 curated MSigDB pathways versus an
#   unbiased data-driven univariable signature). Survival performance is
#   evaluated via stratified 5-fold cross-validation with invariance guards
#   across seven model configurations, culminating in the export of report-ready
#   tables and publication-quality forest and Kaplan-Meier plots.
#
# This script should be run from the project root.
# It sources each analysis script in order, prints clear progress messages,
# stops immediately if a step fails, and saves session information at the end
# for reproducibility.
#
# Note: must be run from the project root (the folder containing R/).
# An RStudio project (.Rproj) sets the working directory automatically.
# When using Rscript, cd to the project root first; do not rely on the
# shell's current directory matching the project root by coincidence.


# Configuration and Utilities --------------------------------------------------

# Ensure the 'cli' package is installed before loading it defensively
# For styled console output and consistent messaging
if (!requireNamespace("cli", quietly = TRUE)) {
  install.packages("cli")
}

library(cli) # Load cli early for consistent messaging


# Project Structure Validation -------------------------------------------------

if (!dir.exists("R")) {
  stop(
    "Directory 'R' not found.\n",
    "run_analysis.R must be run from the project root.\n",
    "Current working directory: ", getwd(),
    call. = FALSE
  )
}


# Output Directory Setup -------------------------------------------------------

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}


# Pipeline Log File Initialization --------------------------------------------

log_file_path <- file.path("results", "pipeline_log.txt")

# Initialize the pipeline log with a timestamped header.
# This overwrites any existing file so each run starts fresh.
writeLines(
  text = c(
    strrep("=", 60),
    paste("Pipeline run started:", Sys.time()),
    strrep("=", 60)
  ),
  con = log_file_path
)


# Logging Helpers --------------------------------------------------------------

# Send a message to the console and append the same text to the pipeline log.
log_message <- function(...) {
  message_text <- paste0(...)
  message(message_text)
  write(
    x = message_text,
    file = log_file_path,
    append = TRUE
  )
}

# Convert elapsed seconds into a human-readable string using ms for sub-second
# durations and seconds for longer durations.
format_elapsed_time <- function(time_seconds) {
  if (time_seconds < 1) {
    paste0(round(time_seconds * 1000, 1), " ms")
  } else {
    paste0(round(time_seconds, 2), " s")
  }
}


# Pipeline Step Definitions ----------------------------------------------------
# Re-sequenced to avoid circular scoping crashes and missing variables!

pipeline_steps <- tibble::tribble(
  ~path,                                    ~step_name,
  # Phase 1: Core Utilities & Workspace Initialization
  file.path("R", "utils_validation.R"),     "Validation error helpers",
  file.path("R", "setup.R"),                "Package setup and configuration",
  file.path("R", "load_data.R"),            "Load raw cBioPortal data",

  # Phase 2: Individual Layer Data Cleaning (Isolation Stage)
  file.path("R", "prepare_clinical.R"),     "Prepare clinical survival table",
  file.path("R", "screen_clinical.R"),      "Screen for clinical predictors",
  file.path("R", "prepare_rppa.R"),         "Prepare RPPA proteomics data",
  file.path("R", "prepare_rnaseq.R"),       "Prepare RNA-seq data",
  #file.path("R", "prepare_mutations.R"),    "Prepare binary mutation features",
  file.path("R", "prepare_cna.R"),          "Prepare binary CNA features",

  # Phase 3: Master Cohort Integration & Clinical Benchmarking
  file.path("R", "integrate_data.R"),       "Integrate multiomics data layers",
  file.path("R", "survival_check.R"),       "Baseline cohort survival check",

  # Phase 4: Downstream Multi-Omics Feature Extraction
  file.path("R", "feature_selection.R"),    "Screen RPPA: multivariable LASSO",
  file.path("R", "pathway_scores.R"),       "Calculate RNA pathway scores",

  # Phase 5: Final Cross-Validated Modelling & Publication Figures
  file.path("R", "survival_models.R"),      "Compare cross-validated models",
  file.path("R", "results_figures.R"),      "Create report-ready outputs"
)


# Verify Script Existence ------------------------------------------------------

# Create a tibble of pipeline steps whose script files are missing from disk.
missing_scripts <- pipeline_steps |>
  dplyr::filter(!file.exists(path))

if (nrow(missing_scripts) > 0) {
  stop(
    "The following pipeline scripts were not found:\n",
    paste(missing_scripts$path, collapse = "\n"),
    call. = FALSE
  )
}


# Helper function for one pipeline step ----------------------------------------

# Execute a single pipeline script while printing and logging progress,
# handling messages, warnings, and errors, and measuring elapsed time.
# This is all done using coloured console output for easier reading.
source_step <- function(path, step_name) {
  plain_line <- strrep("-", 60)
  plain_step <- paste0("Step: ", step_name, " (`", path, "`)")

  # Write the condensed step header to the console with clean formatting
  cli::cli_text(cli::col_black(plain_line))
  cli::cli_text(cli::style_bold(cli::col_black(plain_step)))
  cli::cli_text(cli::col_black(plain_line))

  # Append the identical condensed layout into your tracking file silently
  write(
    x      = c(plain_line, plain_step, plain_line),
    file   = log_file_path,
    append = TRUE
  )

  # Start a timer to measure how long the script is taking
  t_start <- proc.time()

  # Execute the script with automated error interception
  withCallingHandlers(
    tryCatch(
      source(path),
      error = function(e) {
        error_message <- paste0(
          "Pipeline failed at step: ", step_name, "\n",
          "Script: ", path, "\n",
          "Error: ", conditionMessage(e)
        )
        write(x = error_message, file = log_file_path, append = TRUE)
        stop(error_message, call. = FALSE)
      }
    ),
    message = function(m) {
      message_text <- conditionMessage(m)

      # Suppress harmless internal ggplot2 aesthetic note pollution cleanly
      if (grepl("Ignoring unknown labels", message_text)) {
        invokeRestart("muffleMessage")
      }

      cli::cli_text(paste0(
        cli::col_cyan("\u2139\ufe0f [MESSAGE] "),
        message_text
      ))
      write(
        x      = paste0("[MESSAGE] ", message_text),
        file   = log_file_path,
        append = TRUE
      )
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warning_text <- conditionMessage(w)
      cli::cli_text(paste0(
        cli::col_yellow("\u26a0\ufe0f [WARNING] "),
        cli::col_yellow(step_name), ": ", warning_text
      ))
      write(
        x      = paste0("[WARNING] ", step_name, ": ", warning_text),
        file   = log_file_path,
        append = TRUE
      )
      invokeRestart("muffleWarning")
    }
  )

  elapsed <- unname((proc.time() - t_start)[["elapsed"]])
  formatted_time <- format_elapsed_time(elapsed)

  cli::cli_text(cli::style_bold(cli::col_green(paste0(
    "\u2705 Completed ", step_name, " in ",
    formatted_time, "."
  ))))

  write(
    x      = paste0("Completed ", step_name, " in ", formatted_time, "."),
    file   = log_file_path,
    append = TRUE
  )
}


# Run Pipeline -----------------------------------------------------------------

# Create a start time reference point.
t_pipeline_start <- proc.time()

# Loop through each pipeline step and execute it using the source_step function.
for (i in seq_len(nrow(pipeline_steps))) {
  source_step(
    path = pipeline_steps$path[[i]],
    step_name = pipeline_steps$step_name[[i]]
  )
}

total_elapsed <- unname((proc.time() - t_pipeline_start)[["elapsed"]])
formatted_total_time <- format_elapsed_time(total_elapsed)


# Save Session Information -----------------------------------------------------

session_info_path <- file.path("results", "session_info.txt")

writeLines(
  text = capture.output(sessionInfo()),
  con = session_info_path
)

# Write a final completion message to the console and the pipeline log.
cli::cli_text(cli::col_green(strrep("=", 60)))
cli::cli_text(cli::style_bold(cli::col_green(
  paste0(
    "\ud83c\udf89 Pipeline complete. Total time: ",
    formatted_total_time, "."
  )
)))
cli::cli_text(cli::col_green(strrep("=", 60)))

write(
  x = c(
    strrep("=", 60),
    paste0("Pipeline complete. Total time: ", formatted_total_time, "."),
    strrep("=", 60)
  ),
  file = log_file_path,
  append = TRUE
)
