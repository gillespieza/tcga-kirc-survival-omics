#' @file        run_analysis.R
#' @title       Run Full Analysis Workflow
#' @description Orchestrates the complete TCGA KIRC survival analysis pipeline
#'              by sourcing each numbered script in sequence. Each step is
#'              timed and wrapped in error handling so that failures surface
#'              with a clear step name and the underlying error message, rather
#'              than a bare R traceback.
#'
#' @details
#'   Pipeline steps executed in order:
#'   \enumerate{
#'     \item \code{00_setup.R}           — Install/load packages; define paths
#'     \item \code{01_load_data.R}       — Load raw cBioPortal data files
#'     \item \code{02_prepare_clinical.R} — Build clinical survival table
#'     \item \code{03_prepare_rppa.R}    — Reshape and clean RPPA data
#'     \item \code{04_prepare_mutations.R} — Build binary mutation feature table
#'     \item \code{05_integrate_data.R}  — Integrate all data layers
#'   }
#'
#' @note  This script must be run from the project root directory (the folder
#'        containing the \code{R/} subdirectory), or via an RStudio project
#'        (\code{.Rproj}) which sets the working directory automatically.
#'        Running from any other location will cause the \code{source()} calls
#'        to fail with a "cannot open file" error.
#'
#' @usage source("run_analysis.R")
#' @usage Rscript run_analysis.R


# Verify working directory ----------------------------------------------------
# Catch the most common failure mode (wrong working directory) before any
# script is sourced.

if (!dir.exists("R")) {
   stop(
      "Directory 'R/' not found. ",
      "run_analysis.R must be run from the project root. ",
      "Current working directory: ", getwd()
   )
}


# Pipeline helper -------------------------------------------------------------

#' Source a pipeline step with timing and error handling
#'
#' @param path      Character. Relative path to the R script to source.
#' @param step_name Character. Human-readable label used in log messages.
source_step <- function(path, step_name) {
   message("\n", strrep("-", 60))
   message("Step: ", step_name)
   message(strrep("-", 60))
   t_start <- proc.time()
   tryCatch(
      source(path),
      error = function(e) {
         stop(
            "Pipeline failed at step: ", step_name, "\n",
            "Script: ", path, "\n",
            "Error:  ", conditionMessage(e),
            call. = FALSE
         )
      }
   )
   elapsed <- round((proc.time() - t_start)[["elapsed"]], 1)
   message("Completed '", step_name, "' in ", elapsed, "s.")
}


# Run pipeline ----------------------------------------------------------------

t_pipeline_start <- proc.time()

source_step("R/00_setup.R",             "Package setup and configuration")
source_step("R/01_load_data.R",         "Load raw cBioPortal data")
source_step("R/02_prepare_clinical.R",  "Prepare clinical survival table")
source_step("R/03_prepare_rppa.R",      "Prepare RPPA proteomics data")
source_step("R/04_prepare_mutations.R", "Prepare binary mutation features")
source_step("R/05_integrate_data.R",    "Integrate all data layers")

total_elapsed <- round((proc.time() - t_pipeline_start)[["elapsed"]], 1)

message("\n", strrep("=", 60))
message("Pipeline complete. Total time: ", total_elapsed, "s.")
message(strrep("=", 60))