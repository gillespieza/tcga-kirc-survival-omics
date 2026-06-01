# Run full analysis workflow ---------------------------------------------------
#
# Orchestrates the complete TCGA KIRC survival analysis pipeline by sourcing
# each numbered script in sequence. Each step is timed and wrapped in error
# handling so that failures surface with a clear step name and error message
# rather than a bare R traceback. Warnings are also captured and re-emitted
# with the step name so they are not silently swallowed.
#
# Pipeline steps:
#   00_setup.R                - Install/load packages; define paths
#   01_load_data.R            - Load raw cBioPortal data files
#   02_prepare_clinical.R     - Build clinical survival table
#   03_prepare_rppa.R         - Reshape and clean RPPA data
#   04_prepare_mutations.R    - Build binary mutation feature table
#   05_integrate_data.R       - Integrate all data layers
#   06_quick_survival_check.R - Overall KM and clinical Cox checks
#   07_feature_selection.R    - RPPA univariable Cox feature selection
#   08_survival_models.R      - Clinical/omics/integrated model comparison
#   09_results_figures.R      - Save report-ready figures and tables
#
# Usage (from project root):
#   source("run_analysis.R")          # inside an RStudio project
#   cd /path/to/project && Rscript run_analysis.R   # from the shell
#
# Note: must be run from the project root (the folder containing R/).
# An RStudio project (.Rproj) sets the working directory automatically.
# When using Rscript, cd to the project root first; do not rely on the
# shell's current directory matching the project root by coincidence.


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
# Sources a script with per-step timing, error handling, and warning capture.
# Warnings are re-emitted with the step name prepended so they remain visible
# and attributable rather than being swallowed by tryCatch.
# Failures re-throw with the step name and script path prepended to the error
# message so the point of failure is immediately clear in the console output.

source_step <- function(path, step_name) {
   message("\n", strrep("-", 60))
   message("Step: ", step_name)
   message(strrep("-", 60))
   
   t_start <- proc.time()
   
   withCallingHandlers(
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
      ),
      warning = function(w) {
         message("Warning in '", step_name, "': ", conditionMessage(w))
         invokeRestart("muffleWarning")
      }
   )
   
   elapsed <- round((proc.time() - t_start)[["elapsed"]], 1)
   message("Completed '", step_name, "' in ", elapsed, "s.")
}


# Run pipeline ----------------------------------------------------------------

t_pipeline_start <- proc.time()

source_step("R/00_setup.R",                "Package setup and configuration")
source_step("R/01_load_data.R",            "Load raw cBioPortal data")
source_step("R/02_prepare_clinical.R",     "Prepare clinical survival table")
source_step("R/03_prepare_rppa.R",         "Prepare RPPA proteomics data")
source_step("R/04_prepare_mutations.R",    "Prepare binary mutation features")
source_step("R/05_integrate_data.R",       "Integrate all data layers")
source_step("R/06_quick_survival_check.R", "Quick survival check")
source_step("R/07_feature_selection.R",    "Screen RPPA for features")
source_step("R/08_survival_models.R",      "Compare survival models")
source_step("R/09_results_figures.R",      "Create report-ready outputs")

total_elapsed <- round((proc.time() - t_pipeline_start)[["elapsed"]], 1)

message("\n", strrep("=", 60))
message("Pipeline complete. Total time: ", total_elapsed, "s.")
message(strrep("=", 60))


# Session info ----------------------------------------------------------------
# Written to results/ on every successful run for reproducibility. Captures
# R version, platform, and all loaded package versions.

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
session_info_path <- "results/session_info.txt"
writeLines(capture.output(sessionInfo()), session_info_path)
message("Session info written to ", session_info_path)