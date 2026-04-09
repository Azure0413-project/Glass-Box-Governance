# =============================================================================
# run_03_eonly.R — Standalone E-only IECV-EBM Hyperparameter Tuning
# IDH EBM Governance — v3.0 Pure E-only Design
#
# Usage:
#   1. Open RStudio or R console
#   2. src("run_03_eonly.R")
#      or
#      setwd("/path/to/project"); source("run_03_eonly.R")
#
# Description:
#   Within each IECV fold, uses only E-cohort (earlier 75%) for:
#     - Grid search (36 HP combos, 60% subsample, 5-fold grouped CV)
#     - Best HP refinement on full E-only dev pool
#     - Final model training (outer_bags=64)
#     - External validation on E-only held-out center
#   6 iterations (3 centers x 2 seeds) via PSOCK parallel
#
# Estimated runtime: 2-6 hours (depending on CPU cores)
#
# Outputs:
#   <BASE_DIR>/Run_<timestamp>/
#     ├── Iter1_External_TN_Seed1/
#     |     +-- *_tuning_results.csv      <- 36-combo HP OOF results
#     │     ├── *_Final_EBM.joblib        ← E-only trained model
#     │     ├── *_artifacts.rds           ← best HP + threshold
#     |     +-- *_EL_split_meta.csv       <- E/L split metadata
#     │     ├── *_external_predictions.fst
#     │     ├── *_boot_ci.csv
#     │     └── ...
#     ├── Iter2_External_TN_Seed2/ ...
#     +-- ... (6 iterations total)
#     +-- IECV_iteration_summary.csv      <- aggregate summary of all 6 iterations
#     └── IECV_Complete_Results.rds
#
# After completion:
#   1. Update RUN_DIR in 00_config.R to point to the new Run_* directory
#   2. Re-run 05_shapeqc.R to update tier assignments
#   3. Re-run run_el.R for L-only confirmatory evaluation
# =============================================================================

cat("\n")
cat(paste(rep("#", 70), collapse = ""), "\n")
cat("# E-only IECV-EBM Hyperparameter Tuning — Standalone Launcher\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n")

# --- Step 0: Set working directory -----------------------------------------------

SCRIPT_DIR <- "."  # <-- change to your project root directory path

if (!file.exists(file.path(SCRIPT_DIR, "_init.R"))) {
  stop(sprintf("_init.R not found in '%s'. Please update SCRIPT_DIR above.", SCRIPT_DIR),
       call. = FALSE)
}
setwd(SCRIPT_DIR)
source("_init.R")
cat(sprintf("  Working directory: %s\n\n", SCRIPT_DIR))

# ─── Step 1: Load dependencies ───────────────────────────────────────────────

cat("[Step 1] Loading configuration and utilities ...\n\n")

src("R/00_config.R")
src("R/00_utils_r.R")
src("R/00_utils_python.R")

# Verify E/L split config exists
stopifnot(
  "SESSION_DATE_COL not found — update 00_config.R" = exists("SESSION_DATE_COL"),
  "E_SPLIT_QUANTILE not found — update 00_config.R" = exists("E_SPLIT_QUANTILE"),
  "COHORT_COL not found — update 00_config.R"       = exists("COHORT_COL")
)

cat(sprintf("  SESSION_DATE_COL : %s\n", SESSION_DATE_COL))
cat(sprintf("  E_SPLIT_QUANTILE : %.0f%% E / %.0f%% L\n",
            E_SPLIT_QUANTILE * 100, (1 - E_SPLIT_QUANTILE) * 100))
cat(sprintf("  EBM Grid         : %d combinations\n", nrow(EBM_GRID)))
cat(sprintf("  Workers          : %d × %d threads = %d total\n",
            N_WORKERS, N_THREADS_PER_WORKER, N_WORKERS * N_THREADS_PER_WORKER))

# ─── Step 2: Preview E/L split (before committing to full run) ───────────────

cat("\n[Step 2] E/L split preview ...\n\n")

for (site in c("TN", "D6", "CY")) {
  dt_raw <- read_fst_dt(SITE_FILES[[site]])
  if (!(SESSION_DATE_COL %in% names(dt_raw))) {
    # Try after compat rename
    dt_raw <- apply_compat_rename(dt_raw, COMPAT_RENAME_MAP)
  }
  if (SESSION_DATE_COL %in% names(dt_raw)) {
    dates <- as.Date(dt_raw[[SESSION_DATE_COL]])
    cutpoint <- quantile(dates, probs = E_SPLIT_QUANTILE, type = 1, na.rm = TRUE)
    n_E <- sum(dates <= cutpoint, na.rm = TRUE)
    n_L <- sum(dates >  cutpoint, na.rm = TRUE)
    cat(sprintf("  %s: %s sessions → E=%s (≤%s) + L=%s\n",
                site,
                format(nrow(dt_raw), big.mark = ","),
                format(n_E, big.mark = ","),
                as.character(cutpoint),
                format(n_L, big.mark = ",")))
  } else {
    cat(sprintf("  %s: [WARNING] Column '%s' not found!\n", site, SESSION_DATE_COL))
  }
  rm(dt_raw); gc(verbose = FALSE)
}

# ─── Step 3: Run E-only IECV-EBM ─────────────────────────────────────────────

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Starting 03_iecv_ebm.R (E-only v3.0)\n")
cat(sprintf("  Output will be saved to: %s/Run_%s\n", BASE_DIR, RUN_TS))
cat(paste(rep("=", 70), collapse = ""), "\n\n")

run_start <- Sys.time()

src("pipeline/03_iecv_ebm.R")

run_elapsed <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))

# ─── Step 4: Post-run summary ────────────────────────────────────────────────

NEW_RUN_DIR <- file.path(BASE_DIR, paste0("Run_", RUN_TS))

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  E-only IECV-EBM COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("  Total elapsed    : %.1f minutes\n", run_elapsed))
cat(sprintf("  Output directory : %s\n\n", NEW_RUN_DIR))

# Check results
if (file.exists(file.path(NEW_RUN_DIR, "IECV_iteration_summary.csv"))) {
  summ <- data.table::fread(file.path(NEW_RUN_DIR, "IECV_iteration_summary.csv"))
  if (nrow(summ) > 0 && "external_AUROC" %in% names(summ)) {
    cat("  === E-only Results Summary ===\n")
    for (site in c("TN", "D6", "CY")) {
      row <- summ[external == site & seed_id == 1]
      if (nrow(row) > 0) {
        cat(sprintf("  %s: leaf=%d, smooth=%d, pos=%.1f → AUPRC=%.4f, AUROC=%.4f\n",
                    site, row$hp_min_samples_leaf, row$hp_smoothing_rounds,
                    row$hp_pos_multiplier, row$external_AUPRC, row$external_AUROC))
      }
    }
    cat(sprintf("\n  Mean: AUROC %.4f ± %.4f | AUPRC %.4f ± %.4f\n",
                mean(summ$external_AUROC, na.rm = TRUE), sd(summ$external_AUROC, na.rm = TRUE),
                mean(summ$external_AUPRC, na.rm = TRUE), sd(summ$external_AUPRC, na.rm = TRUE)))
  }
}

# ─── Step 5: Next steps guidance ─────────────────────────────────────────────

cat("\n")
cat(paste(rep("─", 70), collapse = ""), "\n")
cat("  NEXT STEPS:\n")
cat(paste(rep("─", 70), collapse = ""), "\n\n")

cat(sprintf('  1. Update RUN_DIR in 00_config.R:\n'))
cat(sprintf('     RUN_DIR <- "%s"\n\n', gsub("\\\\", "/", NEW_RUN_DIR)))

cat('  2. Re-run ShapeQC (on E-only models):\n')
cat('     src("pipeline/05_shapeqc.R")\n\n')

cat('  3. Update tier assignments in 00_config.R based on new 05 results\n\n')

cat('  4. Re-run L-only confirmatory:\n')
cat('     src("run_el.R")\n\n')

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Done.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
