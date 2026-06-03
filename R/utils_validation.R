# Utility validation helpers --------------------------------------------------
#
# Small shared helpers for common validation patterns used across the pipeline.
# Intended to be sourced by other scripts; not run directly.

abort_if_false <- function(condition, message_text) {
   if (!condition) {
      stop(
         message_text, 
         call. = FALSE
      )
   }
}

check_required_objects <- function(object_names, env = .GlobalEnv) {
   missing <- object_names[
      !vapply(object_names, exists, logical(1), envir = env, inherits = TRUE)
   ]
   
   abort_if_false(
      length(missing) == 0,
      paste(
         "The following required object(s) are missing from the environment:",
         paste(missing, collapse = ", "),
         "\nSource the appropriate setup/prepare script(s) via run_analysis.R."
      )
   )
}

check_has_columns <- function(object_name, columns, env = .GlobalEnv) {
   obj <- get(object_name, envir = env)
   missing_cols <- setdiff(columns, names(obj))
   
   abort_if_false(
      length(missing_cols) == 0,
      paste0(
         object_name, " is missing required column(s): ",
         paste(missing_cols, collapse = ", ")
      )
   )
}

check_has_sample_id <- function(object_name, env = .GlobalEnv) {
   check_has_columns(
      object_name, 
      "sample_id", 
      env = env
   )
   # delegates to check_has_columns() <--- error message is reported there.
}