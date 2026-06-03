# Run the full TCGA KIRC survival analysis pipeline ----------------------------
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

# Project structure check ------------------------------------------------------
# Confirm that the R/ directory exists before trying to source any scripts.

if (!dir.exists("R")) {
   stop(
      "Directory 'R' not found.\n",
      "run_analysis.R must be run from the project root.\n",
      "Current working directory: ", getwd(),
      call. = FALSE
   )
}

# Output directory setup -------------------------------------------------------
# Create the results directory early so that log and session files can be saved.

if (!dir.exists("results")) {
   dir.create("results", recursive = TRUE)
}
if (!dir.exists("figures")) {
   dir.create("figures", recursive = TRUE)
}

# Pipeline log file ------------------------------------------------------------
# Save a plain-text log of the pipeline run in the results directory.

log_file_path <- file.path("results", "pipeline_log.txt")

writeLines(
   text = c(
      strrep("=", 60),
      paste("Pipeline run started:", Sys.time()),
      strrep("=", 60)
   ),
   con = log_file_path
)

# Logging helper ---------------------------------------------------------------
# Write wrapper messages to both the console and the pipeline log file.

log_message <- function(...) {
   message_text <- paste0(...)
   
   message(message_text)
   
   write(
      x = message_text,
      file = log_file_path,
      append = TRUE
   )
}

# Pipeline step definitions ----------------------------------------------------
# Store the pipeline order in one place so it is easy to edit later without
# renaming files.

pipeline_steps <- tibble::tribble(
   ~path,                                    ~step_name,
   file.path("R", "utils_validation.R"),     "Helper functions for validation errors",
   file.path("R", "setup.R"),                "Package setup and configuration",
   file.path("R", "load_data.R"),            "Load raw cBioPortal data",
   file.path("R", "prepare_clinical.R"),     "Prepare clinical survival table",
   file.path("R", "prepare_rppa.R"),         "Prepare RPPA proteomics data",
   file.path("R", "prepare_rnaseq.R"),       "Prepare RNA-seq data",
   file.path("R", "prepare_mutations.R"),    "Prepare binary mutation features",
   file.path("R", "integrate_data.R"),       "Integrate all data layers",
   file.path("R", "quick_survival_check.R"), "Quick survival check",
   file.path("R", "feature_selection.R"),    "Screen RPPA for features",
   file.path("R", "pathway_scores.R"),       "Calculate pathway scores",
   file.path("R", "survival_models.R"),      "Compare survival models",
   file.path("R", "results_figures.R"),      "Create report-ready outputs"
)

# Pipeline script existence check ----------------------------------------------
# Check that every expected script file exists before the pipeline starts.

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
# Source one script, print the step name, capture messages and warnings from
# the sourced script, and stop with a helpful error message if the step fails.

source_step <- function(path, step_name) {
   log_message(strrep("-", 60))
   log_message("Step: ", step_name)
   log_message("Script: ", path)
   log_message(strrep("-", 60))
   
   t_start <- proc.time()
   
   withCallingHandlers(
      tryCatch(
         source(path),
         error = function(e) {
            error_message <- paste0(
               "Pipeline failed at step: ", step_name, "\n",
               "Script: ", path, "\n",
               "Error: ", conditionMessage(e)
            )
            
            write(
               x = error_message,
               file = log_file_path,
               append = TRUE
            )
            
            stop(error_message, call. = FALSE)
         }
      ),
      message = function(m) {
         message_text <- paste0("[MESSAGE] ", conditionMessage(m))
         
         cat(message_text, "\n")
         
         write(
            x = message_text,
            file = log_file_path,
            append = TRUE
         )
         
         invokeRestart("muffleMessage")
      },
      warning = function(w) {
         warning_text <- paste0(
            "[WARNING] ",
            step_name,
            ": ",
            conditionMessage(w)
         )
         
         cat(warning_text, "\n")
         
         write(
            x = warning_text,
            file = log_file_path,
            append = TRUE
         )
         
         invokeRestart("muffleWarning")
      }
   )
   
   elapsed <- round((proc.time() - t_start)[["elapsed"]], 1)
   log_message("Completed ", step_name, " in ", elapsed, " s.")
}

# Run pipeline -----------------------------------------------------------------
# Source each script in order and measure the total runtime of the pipeline.

t_pipeline_start <- proc.time()

for (i in seq_len(nrow(pipeline_steps))) {
   source_step(
      path = pipeline_steps$path[[i]],
      step_name = pipeline_steps$step_name[[i]]
   )
}

total_elapsed <- round((proc.time() - t_pipeline_start)[["elapsed"]], 1)

# Save session information -----------------------------------------------------
# Save R version, platform, and loaded package versions for reproducibility.

session_info_path <- file.path("results", "session_info.txt")

writeLines(
   text = capture.output(sessionInfo()),
   con = session_info_path
)

log_message(strrep("=", 60))
log_message("Pipeline complete. Total time: ", total_elapsed, " s.")
log_message(strrep("=", 60))
log_message("Session info written to ", session_info_path)