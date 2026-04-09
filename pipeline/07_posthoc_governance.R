# =============================================================================
# 07_posthoc_governance.R — Post-hoc Governance: Yellow Smooth + Red Zero
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Table 4 Panel B — Aggregate performance change after editing (ΔAUPRC, ΔAUROC, ΔBrier)
#   Mean-preserve verification (6/6 pass)
#
# v3.0 update (Route B):
#   Loads original EBM from 03_iecv_ebm.R output, applies Yellow smoothing
#   (gold-standard smooth_scores_spline_fixed from 00_utils_r.R Section 20)
#   and Red zeroing in a single pass. Evaluates on E-only external cohort.
#   No longer requires a separate v2.0 smoothing step.
#
# Design:
#   (a) Load original EBM models from 03 output
#   (b) Smooth 3 Yellow features with J-adaptive spar (from 05_shapeqc.R)
#   (c) Zero out Red features with analytic intercept compensation
#   (d) Re-evaluate on E-only external validation sets (no retraining)
#   (e) Verify mean-preserve criterion
#
# Methods → "Three Yellow-tier features underwent targeted smoothing and
#  six Red-tier features were zeroed out"
# Methods → "Mean-preserve verification confirms population-average
#  predicted risk unchanged (discrepancy < 10^-6)"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 03_iecv_ebm.R and 05_shapeqc.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(fst); library(reticulate); library(precrec)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt") || !exists("load_site_eonly", mode = "function") ||
    !exists("load_spar_map", mode = "function") ||
    !exists("smooth_scores_spline_fixed", mode = "function")) src("R/00_utils_r.R")
if (!exists("init_python_full")) src("R/00_utils_python.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 07: Post-hoc Governance (Yellow Smooth + Red Zero)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

ANALYSIS_TS <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
init_python_full()
joblib <- reticulate::import("joblib", convert = TRUE)

# --- Smoothing parameters (gold-standard, matching 04f) ---
GOV_UPDATE_ONLY_PHYS_RANGE <- TRUE
GOV_UPDATE_MISSING_BIN     <- FALSE
GOV_UPDATE_UNKNOWN_BIN     <- FALSE
GOV_CLIP_TO_ORIGINAL_RANGE <- TRUE
GOV_CLIP_MARGIN            <- 0.25

# --- Load SPAR map from 05_shapeqc.R output ---
SPAR_MAP <- load_spar_map(RUN_DIR)

# --- Dynamic tier loading from 05_shapeqc.R output (overrides 00_config.R) ---
# Prevents stale hardcoded tiers when E-only ShapeQC differs from full-cohort.
tryCatch({
  qc_dirs_07 <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
  qc_dirs_07 <- sort(qc_dirs_07[grepl("ShapeQC_v3\\.0_", basename(qc_dirs_07))], decreasing = TRUE)
  if (length(qc_dirs_07) > 0) {
    dm_path <- file.path(qc_dirs_07[1], "Table_4A_Decision_Matrix.csv")
    if (file.exists(dm_path)) {
      dm <- fread(dm_path)
      if ("light" %in% names(dm)) setnames(dm, "light", "tier")
      GREEN_FEATURES  <- dm[tier == "Green",  feature]
      YELLOW_FEATURES <- dm[tier == "Yellow", feature]
      RED_FEATURES    <- dm[tier == "Red",    feature]
      GRAY_FEATURES   <- dm[tier == "Gray",   feature]
      ALL_EDIT_FEATURES <- c(YELLOW_FEATURES, RED_FEATURES)
      cat("  [INFO] Tiers loaded dynamically from ShapeQC Decision_Matrix\n")
    }
  }
}, error = function(e) {
  cat(sprintf("  [WARN] Dynamic tier load failed (%s); using 00_config.R hardcoded tiers\n", e$message))
})

get_spar <- function(feat) {
  val <- SPAR_MAP[feature == feat, spar]
  if (length(val) == 0 || is.na(val)) {
    warning(sprintf("[GOV] Feature '%s' not in SPAR map; using fallback 0.70", feat))
    return(0.70)
  }
  val[1]
}

# --- Output directory ---
OUT_ROOT <- file.path(RUN_DIR, paste0("Posthoc_v3.0_Governance_", ANALYSIS_TS))
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Output: %s\n", basename(OUT_ROOT)))
cat(sprintf("  Yellow (smooth): %s\n", paste(YELLOW_FEATURES, collapse = ", ")))
cat(sprintf("  Red    (zero):   %s\n\n", paste(RED_FEATURES, collapse = ", ")))

# =============================================================================
# Main Pipeline: Per-iteration Yellow smooth + Red zero + E-only evaluation
# =============================================================================

iter_tags <- vapply(ITERATIONS, function(it)
  sprintf("Iter%d_External_%s_Seed%d", it$iter, it$external, it$seed_id), character(1))

all_perf   <- data.table()
all_mp     <- data.table()
site_cache <- list()

for (i in seq_along(ITERATIONS)) {
  it       <- ITERATIONS[[i]]
  iter_tag <- iter_tags[i]
  ext_site <- it$external

  cat(sprintf("\n  === %s [External: %s] ===\n", iter_tag, ext_site))
  out_iter <- file.path(OUT_ROOT, iter_tag)
  dir.create(out_iter, recursive = TRUE, showWarnings = FALSE)

  # --- Load original EBM from 03 output ---
  orig_iter  <- file.path(RUN_DIR, iter_tag)
  model_path <- find_latest_file(orig_iter, "_Final_EBM\\.joblib$")
  if (is.na(model_path)) { cat("  [Skip] No original EBM found\n"); next }

  ebm_baseline <- tryCatch(joblib$load(model_path), error = function(e) NULL)
  ebm_new      <- tryCatch(joblib$load(model_path), error = function(e) NULL)
  if (is.null(ebm_baseline) || is.null(ebm_new)) { cat("  [Skip] Load failed\n"); next }

  shift_total <- 0.0
  n_smoothed  <- 0L
  n_zeroed    <- 0L

  # ── Phase 1: Smooth Yellow features ──────────────────────────────────────
  for (feat in YELLOW_FEATURES) {
    bins <- tryCatch(
      extract_bins_with_diagnostics(ebm_new, feat,
                                     strict_continuous_cuts = TRUE, verbose = FALSE),
      error = function(e) list(success = FALSE, error_msg = e$message)
    )
    if (!bins$success) { cat(sprintf("  [Skip] %s: %s\n", feat, bins$error_msg)); next }

    original_scores <- bins$dt$score
    x_vec  <- bins$dt$x
    bin_w  <- if ("bin_weight" %in% names(bins$dt)) bins$dt$bin_weight else rep(1, length(x_vec))

    # Domain restriction: use PHYS_RANGES if available, skip if index-based feature
    phys_range_use       <- NULL
    update_only_phys_use <- GOV_UPDATE_ONLY_PHYS_RANGE
    if (feat %in% names(PHYS_RANGES)) phys_range_use <- PHYS_RANGES[[feat]]
    if (isTRUE(bins$x_is_index)) {
      phys_range_use       <- NULL
      update_only_phys_use <- FALSE
    }

    feat_spar <- get_spar(feat)

    sm <- smooth_scores_spline_fixed(
      x = x_vec, y_old = original_scores, w = bin_w,
      spar = feat_spar,
      phys_range = phys_range_use,
      update_only_phys = update_only_phys_use,
      clip_to_old = GOV_CLIP_TO_ORIGINAL_RANGE,
      clip_margin = GOV_CLIP_MARGIN
    )
    y_new <- sm$y_new

    # Special-bin protection: preserve missing/unknown/non-finite x bins
    idx_no_x <- which(!is.finite(x_vec))
    if (length(idx_no_x) > 0) y_new[idx_no_x] <- original_scores[idx_no_x]
    if (!GOV_UPDATE_MISSING_BIN) {
      idx_m <- which(bins$dt$is_missing)
      if (length(idx_m) > 0) y_new[idx_m] <- original_scores[idx_m]
    }
    if (!GOV_UPDATE_UNKNOWN_BIN) {
      idx_u <- which(bins$dt$is_unknown)
      if (length(idx_u) > 0) y_new[idx_u] <- original_scores[idx_u]
    }

    mi <- if (is.finite(bins$missing_idx)) as.integer(bins$missing_idx) else as.integer(-1)
    ui <- if (is.finite(bins$unknown_idx)) as.integer(bins$unknown_idx) else NULL

    delta <- tryCatch({
      d <- py_update_term(ebm_new, as.integer(bins$term_index),
                          reticulate::r_to_py(y_new),
                          mean_preserve = TRUE,
                          preserve_missing = !GOV_UPDATE_MISSING_BIN,
                          missing_index = mi,
                          preserve_unknown = !GOV_UPDATE_UNKNOWN_BIN,
                          unknown_index = ui,
                          shift_bagged_scores = TRUE,
                          shift_bagged_intercept = TRUE)
      suppressWarnings(as.numeric(reticulate::py_to_r(d)))
    }, error = function(e) NA_real_)

    if (is.finite(delta)) {
      shift_total <- shift_total + delta
      n_smoothed  <- n_smoothed + 1L
      cat(sprintf("  [SMOOTH] %-25s spar=%.2f shift=%+.6f %s\n",
                  feat, feat_spar, delta,
                  ifelse(sm$fit_ok, "", "[fit_failed]")))
    }
  }

  # ── Phase 2: Zero Red features ──────────────────────────────────────────
  for (feat in RED_FEATURES) {
    bins <- tryCatch(
      extract_bins_with_diagnostics(ebm_new, feat, strict_continuous_cuts = TRUE),
      error = function(e) list(success = FALSE, error_msg = e$message)
    )
    if (!bins$success) { cat(sprintf("  [Skip] %s: %s\n", feat, bins$error_msg)); next }

    new_scores <- rep(0.0, bins$n_scores)
    mi <- if (is.finite(bins$missing_idx)) as.integer(bins$missing_idx) else -1L
    ui <- if (is.finite(bins$unknown_idx)) as.integer(bins$unknown_idx) else NULL

    delta <- tryCatch({
      d <- py_update_term(ebm_new, as.integer(bins$term_index),
                          reticulate::r_to_py(new_scores),
                          mean_preserve = TRUE, preserve_missing = FALSE, missing_index = mi,
                          preserve_unknown = FALSE, unknown_index = ui,
                          shift_bagged_scores = TRUE, shift_bagged_intercept = TRUE)
      suppressWarnings(as.numeric(reticulate::py_to_r(d)))
    }, error = function(e) NA_real_)

    if (!is.finite(delta)) { cat(sprintf("  [Error] %s: non-finite delta\n", feat)); next }
    shift_total <- shift_total + delta
    n_zeroed    <- n_zeroed + 1L
    cat(sprintf("  [ZEROED] %-25s shift=%+.6f\n", feat, delta))
  }

  if (n_smoothed == 0 && n_zeroed == 0) { cat("  [Skip] No features edited\n"); next }

  # ── Mean-preserve verification ──
  mp <- verify_mean_preserve_consistency(ebm_baseline, ebm_new, shift_total)
  cat(sprintf("  [MP] %s (disc=%.2e)\n", mp$status, mp$discrepancy))
  all_mp <- rbind(all_mp, data.table(iter_tag = iter_tag, status = mp$status,
                                      discrepancy = mp$discrepancy), fill = TRUE)

  # ── Save patched model ──
  model_prefix <- sub("_Final_EBM.*\\.joblib$", "", basename(model_path))
  joblib$dump(ebm_new, file.path(out_iter,
              paste0(model_prefix, "_Final_EBM_POSTHOC_v3.0_Governance.joblib")))

  # ── E-only external validation ──
  if (is.null(site_cache[[ext_site]]))
    site_cache[[ext_site]] <- load_site_eonly(ext_site, verbose = FALSE)
  dt_ext <- site_cache[[ext_site]]
  X_ext  <- get_predictor_df(dt_ext)
  y_ext  <- as.integer(dt_ext[[TARGET_COL]])

  p_pre  <- predict_proba_pos_r(ebm_baseline, X_ext)
  p_post <- predict_proba_pos_r(ebm_new, X_ext)

  auc_pre  <- auc_pair(y_ext, p_pre);  auc_post <- auc_pair(y_ext, p_post)
  b_pre    <- brier_score(y_ext, p_pre); b_post <- brier_score(y_ext, p_post)

  perf_row <- data.table(
    iter_tag = iter_tag, external = ext_site,
    n_smoothed = n_smoothed, n_zeroed = n_zeroed,
    AUPRC_pre = auc_pre["auprc"], AUPRC_post = auc_post["auprc"],
    d_AUPRC = auc_post["auprc"] - auc_pre["auprc"],
    AUROC_pre = auc_pre["auroc"], AUROC_post = auc_post["auroc"],
    d_AUROC = auc_post["auroc"] - auc_pre["auroc"],
    Brier_pre = b_pre, Brier_post = b_post, d_Brier = b_post - b_pre,
    mp_status = mp$status
  )
  all_perf <- rbind(all_perf, perf_row, fill = TRUE)

  cat(sprintf("  AUROC: %.4f → %.4f (Δ=%+.4f)\n", perf_row$AUROC_pre, perf_row$AUROC_post, perf_row$d_AUROC))
  cat(sprintf("  AUPRC: %.4f → %.4f (Δ=%+.4f)\n", perf_row$AUPRC_pre, perf_row$AUPRC_post, perf_row$d_AUPRC))
  cat(sprintf("  Brier: %.4f → %.4f (Δ=%+.4f)\n", perf_row$Brier_pre, perf_row$Brier_post, perf_row$d_Brier))

  # Save predictions for downstream NRI (10_reclassification.R)
  dt_nri <- data.table(patient_id = dt_ext[[ID_COL]], cluster_id = dt_ext[[CLUSTER_COL]],
                        outcome = y_ext, p_pre = p_pre, p_post = p_post,
                        iter_tag = iter_tag, held_out_center = ext_site)
  if (SESSION_DATE_COL %in% names(dt_ext)) {
    dt_nri[, Session_Date := as.Date(dt_ext[[SESSION_DATE_COL]])]
  }
  feat_cols_nri <- intersect(ALL_EDIT_FEATURES, names(dt_ext))
  if (length(feat_cols_nri) > 0) {
    for (fc in feat_cols_nri) dt_nri[, (fc) := dt_ext[[fc]]]
  }
  fst::write_fst(dt_nri, file.path(out_iter,
                 paste0(model_prefix, "_predictions_redtier.fst")))

  rm(ebm_baseline, ebm_new)
  if (FORCE_GC_EACH_ITER) gc(verbose = FALSE)
}

# =============================================================================
# Summary: Table 4 Panel B
# =============================================================================

fwrite(all_perf, file.path(OUT_ROOT, "Table_4B_Performance.csv"))
fwrite(all_mp,   file.path(OUT_ROOT, "Mean_Preserve_Verification.csv"))

cat("\n\n  ─── Table 4 Panel B: Post-hoc refinement verification ───\n")
if (nrow(all_perf) > 0) {
  cat(sprintf("  ΔAUPRC range: %+.6f to %+.6f (mean %+.4f ± %.4f)\n",
              min(all_perf$d_AUPRC), max(all_perf$d_AUPRC),
              mean(all_perf$d_AUPRC), sd(all_perf$d_AUPRC)))
  cat(sprintf("  ΔAUROC range: %+.6f to %+.6f (mean %+.4f ± %.4f)\n",
              min(all_perf$d_AUROC), max(all_perf$d_AUROC),
              mean(all_perf$d_AUROC), sd(all_perf$d_AUROC)))
  cat(sprintf("  ΔBrier range: %+.6f to %+.6f (mean %+.4f ± %.4f)\n",
              min(all_perf$d_Brier), max(all_perf$d_Brier),
              mean(all_perf$d_Brier), sd(all_perf$d_Brier)))
  cat(sprintf("  Mean-preserve: %d/%d pass\n",
              sum(all_mp$status == "PASS"), nrow(all_mp)))

  cat("\n  [Verify] Manuscript Table 4B expects:\n")
  cat("    ΔAUPRC: +0.0019 ± 0.0012 | ΔAUROC: +0.0010 ± 0.0017\n")
  cat("    Mean-preserve: 6/6 pass (max disc 5.99e-16)\n")
  cat("    NOTE: v3.0 uses E-only evaluation; expected values may shift\n")
  cat("          slightly from prior full-cohort results.\n")
}

cat(sprintf("\n  Output: %s\n\n", OUT_ROOT))
