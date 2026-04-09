# =============================================================================
# extract_ebm_importance.R  (v3 — weighted, safe Python indexing)
# Extract feature importance directly from trained EBM .joblib models
#
# Official EBM importance = weighted mean of |score|:
#   np.average(np.abs(term_scores_[i]), weights=bin_weights_[i])
#
# v3 fix: uses Python-side 0-indexed access to term_scores_ / bin_weights_
#          to avoid R/Python index confusion with reticulate convert=TRUE
#
# Usage:
#   src("R/00_config.R"); source("00_utils_r.R"); source("00_utils_python.R")
#   src("R/extract_ebm_importance.R")
#
# Standalone usage:
#   RUN_DIR <- "/path/to/Run_YYYY_MM_DD"
#   src("R/extract_ebm_importance.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
})

# --- Setup Python ---
if (!py_available(initialize = FALSE)) {
  use_virtualenv("~/.virtualenvs/r-interpret", required = FALSE)
}
joblib <- import("joblib", convert = TRUE)
np     <- import("numpy",  convert = TRUE)

# Helper: extract term_scores_[idx] and bin_weights_[idx] on Python side
# to avoid R/Python index confusion with reticulate convert=TRUE
py_run_string("
def _ebm_get_term(ebm, idx):
    '''Return (scores_array, weights_array) for term idx (0-based).'''
    import numpy as np
    scores = np.array(ebm.term_scores_[idx], dtype=np.float64).flatten()
    weights = np.array(ebm.bin_weights_[idx], dtype=np.float64).flatten()
    return scores, weights
")
py_get_term <- py$`_ebm_get_term`

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Extract EBM Feature Importance (weighted, official method)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

stopifnot("RUN_DIR must be set" = exists("RUN_DIR"))

iter_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
iter_dirs <- iter_dirs[grepl("^Iter[0-9]+_External_", basename(iter_dirs))]
if (length(iter_dirs) == 0) stop("No iteration directories found in RUN_DIR.")

cat(sprintf("  RUN_DIR: %s\n", RUN_DIR))
cat(sprintf("  Found %d iteration directories\n\n", length(iter_dirs)))

# =============================================================================
# Core extraction: matches InterpretML term_importances("avg_weight")
# =============================================================================

extract_ebm_importance <- function(model_path, verbose = TRUE) {
  ebm <- joblib$load(model_path)
  
  feature_names <- py_to_r(ebm$term_names_)
  n_terms       <- length(feature_names)
  
  results <- list()
  
  for (i in seq_len(n_terms)) {
    fname <- feature_names[i]
    
    # Skip interaction terms
    if (grepl(" x | & ", fname)) next
    
    # Use Python-side 0-indexed access (i-1 because R loop is 1-based)
    sw <- py_get_term(ebm, as.integer(i - 1L))
    scores_raw <- as.numeric(py_to_r(sw[[1]]))
    weights    <- as.numeric(py_to_r(sw[[2]]))
    
    abs_scores <- abs(scores_raw)
    
    # Official EBM importance: weighted mean of |score|
    total_w <- sum(weights)
    if (total_w > 0) {
      importance_weighted <- sum(abs_scores * weights) / total_w
    } else {
      importance_weighted <- mean(abs_scores)
    }
    
    max_abs   <- max(abs_scores)
    score_rng <- max(scores_raw) - min(scores_raw)
    
    results[[fname]] <- data.table(
      feature           = fname,
      importance        = importance_weighted,
      max_abs_score     = max_abs,
      score_range       = score_rng,
      n_bins            = length(scores_raw),
      n_bins_nonzero_wt = sum(weights > 0),
      total_samples     = total_w
    )
  }
  
  dt <- rbindlist(results)
  dt[, rank := frank(-importance, ties.method = "min")]
  setorder(dt, rank)
  
  if (verbose) {
    cat(sprintf("    %s: %d features, total_samples=%.0f\n",
                basename(model_path), nrow(dt), dt$total_samples[1]))
  }
  dt
}

# =============================================================================
# Cross-validation: compare with known-good CSV
# =============================================================================

verify_against_csv <- function(extracted_dt, csv_path) {
  if (!file.exists(csv_path)) {
    cat("  [Verify] Reference CSV not found; skipping verification.\n")
    return(invisible(NULL))
  }
  ref <- fread(csv_path)
  merged <- merge(extracted_dt[, .(feature, model, importance)],
                  ref[, .(feature, model, importance_ref = importance)],
                  by = c("feature", "model"), all = TRUE)
  
  if (nrow(merged) == 0) {
    cat("  [Verify] No overlapping rows; check column names.\n")
    return(invisible(NULL))
  }
  
  merged[, abs_diff := abs(importance - importance_ref)]
  merged[, rel_diff := abs_diff / pmax(importance_ref, 1e-10)]
  
  max_abs <- merged[, max(abs_diff, na.rm = TRUE)]
  max_rel <- merged[, max(rel_diff, na.rm = TRUE)]
  n_match <- merged[abs_diff < 1e-6, .N]
  
  cat(sprintf("  [Verify] Compared %d rows vs reference CSV\n", nrow(merged)))
  cat(sprintf("    Max absolute diff : %.2e\n", max_abs))
  cat(sprintf("    Max relative diff : %.2e\n", max_rel))
  cat(sprintf("    Exact match (<1e-6): %d / %d\n", n_match, nrow(merged)))
  
  if (max_abs > 0.001) {
    cat("    [WARN] Discrepancy > 0.001 detected!\n")
    worst <- merged[order(-abs_diff)][1:5]
    print(worst)
  } else {
    cat("    [OK] All values match within tolerance.\n")
  }
  invisible(merged)
}

# =============================================================================
# Main: extract from all iterations
# =============================================================================

all_imp <- list()

for (iter_dir in iter_dirs) {
  tag <- basename(iter_dir)
  
  model_files <- list.files(iter_dir, pattern = "_Final_EBM\\.joblib$",
                            full.names = TRUE, recursive = FALSE)
  model_files <- model_files[!grepl("POSTHOC|Governance", model_files)]
  
  if (length(model_files) == 0) {
    cat(sprintf("  [Skip] %s: no Final_EBM.joblib found\n", tag))
    next
  }
  
  cat(sprintf("  %s\n", tag))
  dt <- extract_ebm_importance(model_files[1])
  dt[, model := tag]
  all_imp[[tag]] <- dt
}

imp <- rbindlist(all_imp)

cat(sprintf("\n  Total: %d rows (%d models x %d features)\n",
            nrow(imp), length(all_imp), imp[, uniqueN(feature)]))

# =============================================================================
# Verify against existing CSV (if available)
# =============================================================================

ref_csv <- file.path(RUN_DIR, "GLOBAL_feature_importance.csv")
if (!file.exists(ref_csv) && exists("OUT_DIR")) {
  ref_csv <- file.path(OUT_DIR, "GLOBAL_feature_importance.csv")
}
cat("\n")
verify_against_csv(imp, ref_csv)

# =============================================================================
# Summary table
# =============================================================================

summary_dt <- imp[, .(
  mean_importance = mean(importance),
  sd_importance   = sd(importance),
  min_importance  = min(importance),
  max_importance  = max(importance),
  cv              = sd(importance) / mean(importance),
  mean_max_abs    = mean(max_abs_score),
  mean_range      = mean(score_range)
), by = feature]

# Normalized share (within-model, then averaged)
imp[, norm_share := importance / sum(importance), by = model]
norm_summary <- imp[, .(mean_norm_share = mean(norm_share)), by = feature]
summary_dt <- merge(summary_dt, norm_summary, by = "feature")

summary_dt[, consensus_rank := frank(-mean_importance, ties.method = "min")]
setorder(summary_dt, consensus_rank)

cat("\n  --- Feature Importance Summary (weighted, official) ---\n\n")
cat(sprintf("  %-25s %8s %8s %6s %8s %5s\n",
            "Feature", "Mean", "SD", "CV", "Share%", "Rank"))
cat(sprintf("  %s\n", paste(rep("-", 68), collapse = "")))
for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  cat(sprintf("  %-25s %8.4f %8.4f %6.2f %7.1f%% %5d\n",
              r$feature, r$mean_importance, r$sd_importance,
              r$cv, r$mean_norm_share * 100, r$consensus_rank))
}

top3_share <- summary_dt[consensus_rank <= 3, sum(mean_norm_share)] * 100
cat(sprintf("\n  Top-3 normalized share: %.1f%%\n", top3_share))

# =============================================================================
# Export
# =============================================================================

out_dir <- RUN_DIR
fwrite(imp,        file.path(out_dir, "GLOBAL_feature_importance_raw.csv"))
fwrite(summary_dt, file.path(out_dir, "GLOBAL_feature_importance_summary.csv"))

cat(sprintf("\n  Saved:\n"))
cat(sprintf("    %s/GLOBAL_feature_importance_raw.csv\n", out_dir))
cat(sprintf("    %s/GLOBAL_feature_importance_summary.csv\n", out_dir))
cat("\n  Done.\n\n")