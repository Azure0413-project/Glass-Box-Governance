# =============================================================================
# run_el.R — Master Orchestrator for E/L Temporal Split Pipeline
# IDH EBM Governance — Plan A: E/L Temporal Split with Protocol Lock
#
# v3.0 update: Pure E-only design
#   03_iecv_ebm.R v3.0 now performs E-only tuning/training/validation directly.
#   Phase 2 (04d) is deprecated — 03 output serves as E-only derivation.
#   Phase 4 (04f) reads E-only models from RUN_DIR and L-partition on-the-fly.
#
# Before running, ensure:
#   1. Paths in 00_config.R are properly configured
#   2. 03_iecv_ebm.R v3.0 has completed (E-only IECV artifacts in RUN_DIR)
#   3. 05_shapeqc.R has completed (on E-only models for tier assignments)
#   4. Python virtualenv is available
#   5. 00_utils_r.R includes apply_el_split() (Section 19)
#
# ┌──────────────────────────────────────────────────────────────┐
# │  Phase 1: Split Table + Frozen Analysis Plan   (04c)        │
# │  Phase 2: DEPRECATED — see 03_iecv_ebm.R v3.0  (04d shim)  │
# │  Phase 3: E-only ShapeQC + Agreement Check      (04e)        │
# │                                                              │
# │  ══════════  PROTOCOL LOCK  ══════════                       │
# │  (Review FAP, confirm tiers, then proceed)                   │
# │                                                              │
# │  Phase 4: L-only Confirmatory Scoring           (04f)        │
# │  Phase 5: Report Generation                     (04g)        │
# └──────────────────────────────────────────────────────────────┘
# =============================================================================
# Usage: EL_SCRIPT_DIR <- "/path/to/S2"; source(file.path(EL_SCRIPT_DIR, "run_el.R"))

EL_N_WORKERS <- 6L   # Used by 04d (PSOCK parallelism); default fallback is 3L

# --- Working directory: auto-detect from source() path, with fallback ---
if (!exists("EL_SCRIPT_DIR")) {
  EL_SCRIPT_DIR <- tryCatch({
    ofile <- sys.frame(1)$ofile
    if (!is.null(ofile)) dirname(normalizePath(ofile)) else getwd()
  }, error = function(e) getwd())
}

# Validate: _init.R must exist in the project root
if (!file.exists(file.path(EL_SCRIPT_DIR, "_init.R"))) {
  stop(sprintf(
    paste0("_init.R not found in '%s'.\n",
           "  Fix: setwd() to the project root before source('run_el.R'),\n",
           "  or:  EL_SCRIPT_DIR <- '/path/to/project' before source('run_el.R')"),
    EL_SCRIPT_DIR), call. = FALSE)
}

setwd(EL_SCRIPT_DIR)
source("_init.R")
cat(sprintf("  [run_el] Working directory: %s\n\n", EL_SCRIPT_DIR))

cat("\n")
cat(paste(rep("#", 70), collapse = ""), "\n")
cat("# IDH EBM Governance — E/L Temporal Split Pipeline (v3.0 Pure E-only)\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n")

t0 <- Sys.time()

# --- Foundation ---
src("R/00_config.R")
src("R/00_utils_r.R")
src("R/00_config_el.R")
# Note: 00_utils_python.R loaded on-demand by 04d and 04f

# Verify shared functions are available
if (!exists("temporal_filter", mode = "function"))
  stop("temporal_filter() not found. Check that 00_utils_r.R Section 18 is intact.")
if (!exists("apply_el_split", mode = "function"))
  stop("apply_el_split() not found. Check that 00_utils_r.R Section 19 is intact.")

# =============================================================================
# Phase 1: Split Table + Frozen Analysis Plan
# =============================================================================
cat("\n>>> Phase 1: E/L Split Table <<<\n")
src("pipeline/el/04c_el_split_table.R")

# =============================================================================
# Phase 2: E-only Derivation IECV (DEPRECATED — now handled by 03 v3.0)
# =============================================================================
cat("\n>>> Phase 2: E-only Derivation IECV (deprecated shim) <<<\n")
src("pipeline/el/04d_el_derivation_iecv.R")

# =============================================================================
# Phase 3: E-only ShapeQC + Agreement Check
# =============================================================================
cat("\n>>> Phase 3: E-only ShapeQC + Agreement <<<\n")
src("pipeline/el/04e_el_shapeqc_agreement.R")

# =============================================================================
# ══════════  PROTOCOL LOCK  ══════════
#
# Interactive checkpoint:
#   0 = LOCK  — stop here for manual review
#   1 = UNLOCK — proceed to Phase 4 + 5
#
# When unlocked, the pipeline automatically runs BOTH tier configurations:
#   Primary:     E-only tiers (clean derivation–confirmation separation)
#   Sensitivity: Full-data tiers (from 00_config.R)
# Results are saved to separate directories; no overwriting.
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  PROTOCOL LOCK CHECKPOINT\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("  Review before proceeding:\n")
cat(sprintf("    1. %s/Frozen_Analysis_Plan.txt\n", EL_DIR_SPLIT))
cat(sprintf("    2. %s/EL_Expert_Review_Decision.csv\n", EL_DIR_SHAPEQC))
cat(sprintf("    3. %s/EL_Eonly_SPAR_Map.csv\n", EL_DIR_SHAPEQC))
cat("\n  No L data should influence governance decisions.\n\n")

if (!exists("PROCEED_PAST_LOCK") || !isTRUE(PROCEED_PAST_LOCK)) {
  stop(
    paste0(
      "Protocol lock engaged. Review Phase 1–3 artifacts before continuing.\n",
      "To proceed intentionally, set PROCEED_PAST_LOCK <- TRUE and source('run_el.R') again."
    ),
    call. = FALSE
  )
}

# =============================================================================
# Phase 4 + 5: Dual-run (E-only primary + Full-data sensitivity)
#
# Each run writes to its own tagged subdirectory:
#   05_L_Only_Confirmatory_eonly/  +  06_Report_eonly/
#   05_L_Only_Confirmatory_full/   +  06_Report_full/
# =============================================================================

tier_runs <- c("eonly", "full")

for (tier_src in tier_runs) {
  EL_TIER_SOURCE <- tier_src
  run_label <- if (tier_src == "eonly") "PRIMARY (E-only tiers)" else "SENSITIVITY (full-data tiers)"

  cat("\n")
  cat(paste(rep("#", 70), collapse = ""), "\n")
  cat(sprintf("# Phase 4+5: %s\n", run_label))
  cat(paste(rep("#", 70), collapse = ""), "\n")

  cat("\n>>> Phase 4: L-only Confirmatory Scoring <<<\n")
  src("pipeline/el/04f_el_confirmatory.R")

  cat("\n>>> Phase 5: Report Generation <<<\n")
  src("pipeline/el/04g_el_report.R")
}

# =============================================================================
# Done
# =============================================================================

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("  E/L Pipeline complete. Total elapsed: %.1f minutes\n",
            as.numeric(elapsed)))
cat("  Results:\n")
cat(sprintf("    Primary:     %s/05_L_Only_Confirmatory_eonly/\n", EL_OUT_ROOT))
cat(sprintf("    Sensitivity: %s/05_L_Only_Confirmatory_full/\n", EL_OUT_ROOT))
cat(sprintf("    Comparison:  %s/07_Cross_Comparison/\n", EL_OUT_ROOT))
cat(paste(rep("=", 70), collapse = ""), "\n")
