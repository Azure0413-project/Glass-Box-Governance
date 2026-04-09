# =============================================================================
# run_all.R — Master Orchestrator
# IDH EBM Governance Reproducibility Pipeline
#
# Runs all analysis modules in manuscript-logical order.
# Each module can also be run independently (source 00_*.R first).
#
# Before running, ensure:
#   1. Paths in 00_config.R are configured for your local environment
#   2. Python venv has interpret==0.7.4, scikit-learn==1.8.0 installed
#   3. Required R packages are installed (see README.md)
# =============================================================================

# --- Working directory: auto-detect from source() path, with fallback ---
SCRIPT_DIR <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())

if (!file.exists(file.path(SCRIPT_DIR, "_init.R"))) {
  stop(sprintf(
    paste0("_init.R not found in '%s'.\n",
           "  Fix: setwd() to the project root before source('run_all.R')"),
    SCRIPT_DIR), call. = FALSE)
}
setwd(SCRIPT_DIR)
source("_init.R")

cat("\n")
cat(paste(rep("#", 70), collapse = ""), "\n")
cat("# IDH EBM Governance — Full Reproducibility Pipeline (E-only v3.0)\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n")

t0 <- Sys.time()

# --- Foundation ---
src("R/00_config.R")
src("R/00_utils_r.R")
# Note: 00_utils_python.R is loaded on-demand by modules 03-10

# --- Module 1: Data Preparation ---
cat("\n>>> Module 01: Data Preparation <<<\n")
src("pipeline/01_data_prep.R")           # → Table 1A, Supp S7

# --- Module 2: Cross-site Heterogeneity ---
cat("\n>>> Module 02: Cross-site Heterogeneity <<<\n")
src("pipeline/02_heterogeneity.R")       # → Table 1B, Fig 2, Supp Tables 1-2

# --- Module 3: IECV-EBM Training & External Validation ---
cat("\n>>> Module 03: IECV-EBM <<<\n")
src("pipeline/03_iecv_ebm.R")            # → Table 2, Supp S2 [SLOW: ~2-6 hours]
if (exists("OUT_ROOT") && dir.exists(OUT_ROOT)) {
  RUN_DIR <- OUT_ROOT
  cat(sprintf("  [run_all] RUN_DIR updated to current IECV output: %s\n", RUN_DIR))
}

# --- Module 4: XGBoost Benchmark ---
cat("\n>>> Module 04: XGBoost Benchmark <<<\n")
src("pipeline/04_iecv_xgboost.R")        # → Supp S1 [SLOW: ~1-3 hours]

# --- Module 4b: Temporal Hold-out ---
cat("\n>>> Module 04b: Temporal Hold-out <<<\n")
src("pipeline/04b_temporal_holdout.R")    # → Supp S6

# --- Module 5: ShapeQC v3.0 ---
cat("\n>>> Module 05: ShapeQC <<<\n")
src("pipeline/05_shapeqc.R")             # → Table 4A, Supp S3, S14

# --- Module 6: TV Decomposition ---
cat("\n>>> Module 06: TV Decomposition <<<\n")
src("pipeline/06_tv_decomposition.R")    # → Supp S5

# --- Module 7: Post-hoc Governance ---
# Requires: 03 (original EBMs) + 05 (SPAR map from ShapeQC)
cat("\n>>> Module 07: Post-hoc Governance <<<\n")
src("pipeline/07_posthoc_governance.R")   # → Table 4B

# --- Module 8: Sensitivity + LOFO ---
cat("\n>>> Module 08: Sensitivity & LOFO <<<\n")
src("pipeline/08_sensitivity_lofo.R")     # → Supp S8-S11

# --- Module 9: Zeroing Comparison ---
cat("\n>>> Module 09: Zeroing Comparison <<<\n")
src("pipeline/09_zeroing_comparison.R")   # → Supp S12, S15

# --- Module 10: Reclassification ---
cat("\n>>> Module 10: NRI / IDI <<<\n")
src("pipeline/10_reclassification.R")     # → Supp S13

# --- Module 11: Figures ---
cat("\n>>> Module 11: Figures <<<\n")
src("pipeline/11_figures.R")              # → Fig 1-5

# --- Done ---
elapsed <- difftime(Sys.time(), t0, units = "mins")
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("  Pipeline complete. Total elapsed: %.1f minutes\n", as.numeric(elapsed)))
cat(paste(rep("=", 70), collapse = ""), "\n")
