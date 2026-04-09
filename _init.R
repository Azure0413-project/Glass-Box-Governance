# =============================================================================
# _init.R — Project Bootstrap
# IDH EBM Governance Reproducibility Pipeline
#
# Sets the project root directory and provides src() for cross-directory
# source() calls. All runner scripts source this file first.
#
# Usage:
#   source("_init.R")         # from project root
#   src("R/00_config.R")      # instead of src("R/00_config.R")
# =============================================================================

# Detect project root (the directory containing this file)
PROJ_ROOT <- tryCatch({
  # When sourced: use the directory of _init.R itself
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())

# Validate: R/ and pipeline/ must exist
if (!dir.exists(file.path(PROJ_ROOT, "R")) ||
    !dir.exists(file.path(PROJ_ROOT, "pipeline"))) {
  stop(sprintf(
    paste0("Project structure not found in '%s'.\n",
           "  Expected: R/ and pipeline/ subdirectories.\n",
           "  Fix: setwd() to the project root before source('_init.R')"),
    PROJ_ROOT), call. = FALSE)
}

#' Source a file relative to the project root
#'
#' @param rel_path  Path relative to PROJ_ROOT (e.g., "R/00_config.R")
#' @param ...       Additional arguments passed to source()
src <- function(rel_path, ...) {
  full_path <- file.path(PROJ_ROOT, rel_path)
  if (!file.exists(full_path)) {
    stop(sprintf("File not found: %s", full_path), call. = FALSE)
  }
  source(full_path, ...)
}

cat(sprintf("  [init] Project root: %s\n", PROJ_ROOT))
