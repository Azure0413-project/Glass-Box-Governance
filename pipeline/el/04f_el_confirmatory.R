# =============================================================================
# 04f_el_confirmatory.R — L-only Confirmatory Scoring
# IDH EBM Governance — Plan A Phase 4
#
# v3.0 update: Models and thresholds now read from 03_iecv_ebm.R v3.0
#   E-only output (RUN_DIR), not from 04d. L-partition derived on-the-fly
#   via apply_el_split(), removing dependency on 04c split table.
#
# Outputs:
#   Per-iteration: pre-edit & post-edit predictions on L
#   Confirmatory performance tables (per-center, per-policy)
#   Bootstrap CI (patient-level, 1000 resamples)
#   NRI / IDI on L
#
# Plan A Step 5 → "L-only Confirmatory Scoring"
# Plan A Step 6 → "Confirmatory Reporting"
#
# CRITICAL: This script uses LOCKED governance parameters.
# After Protocol Lock, no parameter may be changed based on L results.
#
# Gold standard smoothing: smooth_scores_spline_fixed() ported from
# Post-hoc v2.1, including 3-layer domain restriction.
#
# FIX LOG:
#   - Issue 3: uses shared apply_el_split() from 00_utils_r.R
#   - Issue 5: get_predictor_df() with column validation guard
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       00_config.R, 00_config_el.R, 00_utils_r.R, 00_utils_python.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(reticulate)
  library(precrec)
})

if (!exists("SITE_FILES"))        src("R/00_config.R")
if (!exists("read_fst_dt"))       src("R/00_utils_r.R")
if (!exists("init_python_full"))  src("R/00_utils_python.R")
if (!exists("EL_TEMPORAL_RATIO")) src("R/00_config_el.R")

# Verify shared helpers are available
if (!exists("apply_el_split", mode = "function"))
  stop("apply_el_split() not found. Ensure 00_utils_r.R Section 19 is loaded.")
if (!exists("smooth_scores_spline_fixed", mode = "function"))
  stop("smooth_scores_spline_fixed() not found. Ensure 00_utils_r.R Section 20 is intact.")
if (!exists("load_spar_map", mode = "function"))
  stop("load_spar_map() not found. Ensure 00_utils_r.R Section 20 is intact.")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 04f: L-only Confirmatory Scoring (Protocol Lock)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

ANALYSIS_TS <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
init_python_full()
joblib <- reticulate::import("joblib", convert = TRUE)

# =============================================================================
# SECTION A: Gold Standard Smoothing — now in 00_utils_r.R (Section 20)
# smooth_scores_spline_fixed() is shared by 04f and 07.
# =============================================================================

# =============================================================================
# SECTION B: Load locked parameters
# =============================================================================

cat("[Section B] Loading locked governance parameters ...\n\n")

# Load E-only SPAR map
spar_map_path <- file.path(EL_DIR_SHAPEQC, "EL_Eonly_SPAR_Map.csv")
if (!file.exists(spar_map_path)) {
  cat("  [v3.0] Deriving SPAR map from 05_shapeqc.R output\n")
  SPAR_MAP <- load_spar_map(RUN_DIR)
} else {
  SPAR_MAP <- fread(spar_map_path)
}

# ── Locked tier assignments: E-only (default) or full-data (sensitivity) ──
if (!exists("EL_TIER_SOURCE")) EL_TIER_SOURCE <- "eonly"

if (EL_TIER_SOURCE == "eonly") {
  # Option 1: Read E-only tiers from 04e output (legacy compatibility)
  qc_tiers_path <- file.path(EL_DIR_SHAPEQC, "EL_Eonly_QC_Tiers.csv")
  if (!file.exists(qc_tiers_path)) {
    # Option 1b: Auto-fallback to 05_shapeqc.R E-only output (v3.0 design)
    qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
    qc_dirs <- sort(qc_dirs[grepl("ShapeQC", basename(qc_dirs))], decreasing = TRUE)
    if (length(qc_dirs) > 0) {
      qc_csv <- file.path(qc_dirs[1], "Table_4A_Decision_Matrix.csv")
      if (file.exists(qc_csv)) {
        cat("  [v3.0] Reading tiers from 05_shapeqc.R output (E-only models)\n")
        EL_QC_TIERS <- fread(qc_csv)
        setnames(EL_QC_TIERS, "light", "tier", skip_absent = TRUE)
      } else {
        stop("No tier source found. Run 05_shapeqc.R on E-only models first, OR set EL_TIER_SOURCE='full' to use 00_config.R tiers.")
      }
    } else {
      stop("No ShapeQC output found. Run 05_shapeqc.R first, OR set EL_TIER_SOURCE='full'.")
    }
  } else {
    EL_QC_TIERS <- fread(qc_tiers_path)
  }

  LOCKED_GREEN  <- EL_QC_TIERS[tier == "Green",  feature]
  LOCKED_YELLOW <- EL_QC_TIERS[tier == "Yellow", feature]
  LOCKED_RED    <- EL_QC_TIERS[tier == "Red",    feature]
  LOCKED_GRAY   <- EL_QC_TIERS[tier == "Gray",   feature]

  cat("  Tier source: E-ONLY\n")

  # Update SPAR map to use E-only tiers for spar lookup
  # (SPAR values in the map are already computed from E-only J values)

} else {
  # Legacy: use full-data tiers from 00_config.R
  LOCKED_GREEN  <- GREEN_FEATURES
  LOCKED_YELLOW <- YELLOW_FEATURES
  LOCKED_RED    <- RED_FEATURES
  LOCKED_GRAY   <- GRAY_FEATURES

  cat("  Tier source: FULL-DATA (from 00_config.R — sensitivity only)\n")
}

cat(sprintf("  Green  (%d): %s\n", length(LOCKED_GREEN),  paste(LOCKED_GREEN, collapse = ", ")))
cat(sprintf("  Yellow (%d): %s\n", length(LOCKED_YELLOW), paste(LOCKED_YELLOW, collapse = ", ")))
cat(sprintf("  Red    (%d): %s\n", length(LOCKED_RED),    paste(LOCKED_RED, collapse = ", ")))
cat(sprintf("  Gray   (%d): %s\n", length(LOCKED_GRAY),   paste(LOCKED_GRAY, collapse = ", ")))

# ── Tier-tagged output directories (avoid overwrite between runs) ──
EL_TIER_TAG  <- if (EL_TIER_SOURCE == "eonly") "eonly" else "full"
EL_DIR_CONFIRM <- file.path(EL_OUT_ROOT, paste0("05_L_Only_Confirmatory_", EL_TIER_TAG))
EL_DIR_REPORT  <- file.path(EL_OUT_ROOT, paste0("06_Report_", EL_TIER_TAG))
dir.create(EL_DIR_CONFIRM, recursive = TRUE, showWarnings = FALSE)
dir.create(EL_DIR_REPORT,  recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Output tag : %s → %s\n\n", EL_TIER_TAG, EL_DIR_CONFIRM))

# Locked SPAR values
get_locked_spar <- function(feat) {
  val <- SPAR_MAP[feature == feat, spar]
  if (length(val) == 0 || is.na(val)) {
    warning(sprintf("[LOCK] Feature '%s' not in SPAR map; using fallback 0.70", feat))
    return(0.70)
  }
  val
}

# Under pure E-only design (v3.0), L partition is derived on-the-fly
# via apply_el_split() using SESSION_DATE_COL / E_SPLIT_QUANTILE from 00_config.R.
# EL_SPLIT_TABLE from 04c is no longer required.

# =============================================================================
# SECTION C: Apply governance + score on L (per iteration)
# =============================================================================

cat("\n[Section C] Running confirmatory scoring ...\n\n")

site_cache <- list()
all_confirm <- data.table()
all_mp      <- data.table()
all_nri     <- data.table()

for (it in ITERATIONS) {
  iter_tag <- sprintf("Iter%d_External_%s_Seed%d", it$iter, it$external, it$seed_id)
  ext_site <- it$external

  cat(sprintf("  === %s [External L: %s] ===\n", iter_tag, ext_site))

  iter_dir <- file.path(EL_EONLY_IECV_DIR, iter_tag)
  out_iter <- file.path(EL_DIR_CONFIRM, iter_tag)
  dir.create(out_iter, recursive = TRUE, showWarnings = FALSE)

  # --- Load E-only model (from 03 v3.0 E-only IECV) ---
  model_files <- list.files(iter_dir, pattern = "_Final_EBM\\.joblib$",
                             full.names = TRUE)
  if (length(model_files) == 0) { cat("  [Skip] No E-only model\n\n"); next }

  # --- Load E-only artifacts (threshold from 03 v3.0) ---
  art_files <- list.files(iter_dir, pattern = "_artifacts\\.rds$", full.names = TRUE)
  if (length(art_files) == 0) { cat("  [Skip] No E-only artifacts\n\n"); next }
  art <- readRDS(art_files[length(art_files)])
  thr_E    <- art$threshold
  thr_full <- art$threshold  # Under pure E-only design, no separate full-cohort threshold

  # --- Load external L data (apply_el_split from 00_utils_r.R) ---
  if (is.null(site_cache[[ext_site]])) {
    dt_full <- prep_site_dt(read_fst_dt(SITE_FILES[[ext_site]]), ext_site)
    dt_full <- apply_el_split(dt_full, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
    site_cache[[ext_site]] <- dt_full
  }
  dt_L <- site_cache[[ext_site]][get(COHORT_COL) == "L"]
  cat(sprintf("    L partition: %d sessions\n", nrow(dt_L)))

  # FIX Issue 5: validate predictor columns before extraction
  missing_cols <- setdiff(ALL_PREDICTORS, names(dt_L))
  if (length(missing_cols) > 0) {
    warning(sprintf("[%s] Missing predictor columns in L data: %s",
                    ext_site, paste(missing_cols, collapse = ", ")))
  }

  X_L  <- get_predictor_df(dt_L)
  y_L  <- as.integer(dt_L[[TARGET_COL]])
  cat(sprintf("    L data: %s sessions, IDH=%.2f%%\n",
              format(nrow(dt_L), big.mark = ","), mean(y_L) * 100))

  # ===========================================================
  # For each governance policy (A = primary, B = sensitivity)
  # ===========================================================

  for (pol in EL_POLICIES) {
    cat(sprintf("    --- Policy %s: %s ---\n", pol$id, pol$label))

    # Load fresh model copy
    ebm_pre  <- joblib$load(model_files[length(model_files)])
    ebm_post <- joblib$load(model_files[length(model_files)])

    # p_pre: before any governance
    p_pre <- predict_proba_pos_r(ebm_pre, X_L)

    # ── BUG 2 FIX: Cache original bins from unmodified model ──
    original_bins_cache <- list()
    for (feat_cache in c(LOCKED_YELLOW, LOCKED_RED)) {
      bins_orig <- tryCatch(
        extract_bins_with_diagnostics(ebm_pre, feat_cache,
                                       strict_continuous_cuts = TRUE, verbose = FALSE),
        error = function(e) list(success = FALSE)
      )
      if (bins_orig$success) original_bins_cache[[feat_cache]] <- bins_orig
    }

    shift_total <- 0.0
    n_smoothed  <- 0L
    n_zeroed    <- 0L

    # --- Apply Yellow action ---
    for (feat in LOCKED_YELLOW) {
      bins <- tryCatch(
        extract_bins_with_diagnostics(ebm_post, feat,
                                       strict_continuous_cuts = TRUE, verbose = FALSE),
        error = function(e) list(success = FALSE, error_msg = e$message)
      )
      if (!bins$success) { cat(sprintf("      [Skip] %s: %s\n", feat, bins$error_msg)); next }

      orig_bins <- original_bins_cache[[feat]]
      if (is.null(orig_bins) || !orig_bins$success) {
        cat(sprintf("      [Skip] %s: not in original cache\n", feat)); next
      }
      original_scores <- orig_bins$dt$score
      x_vec  <- orig_bins$dt$x
      bin_w  <- if ("bin_weight" %in% names(orig_bins$dt)) orig_bins$dt$bin_weight else rep(1, length(x_vec))

      phys_range <- NULL
      if (feat %in% names(PHYS_RANGES)) phys_range <- PHYS_RANGES[[feat]]
      phys_range_use <- phys_range
      update_only_phys_use <- isTRUE(EL_UPDATE_ONLY_PHYS_RANGE)
      if (isTRUE(bins$x_is_index)) {
        phys_range_use <- NULL
        update_only_phys_use <- FALSE
      }

      feat_spar <- get_locked_spar(feat)

      sm <- smooth_scores_spline_fixed(
        x = x_vec, y_old = original_scores, w = bin_w,
        spar = feat_spar,
        phys_range = phys_range_use,
        update_only_phys = update_only_phys_use,
        clip_to_old = isTRUE(EL_CLIP_TO_ORIGINAL_RANGE),
        clip_margin = EL_CLIP_MARGIN
      )
      y_new <- sm$y_new

      # Special-bin protection
      idx_no_x <- which(!is.finite(x_vec))
      if (length(idx_no_x) > 0) y_new[idx_no_x] <- original_scores[idx_no_x]
      if (!isTRUE(EL_UPDATE_MISSING_BIN)) {
        idx_m <- which(bins$dt$is_missing)
        if (length(idx_m) > 0) y_new[idx_m] <- original_scores[idx_m]
      }
      if (!isTRUE(EL_UPDATE_UNKNOWN_BIN)) {
        idx_u <- which(bins$dt$is_unknown)
        if (length(idx_u) > 0) y_new[idx_u] <- original_scores[idx_u]
      }

      mi <- if (is.finite(bins$missing_idx)) as.integer(bins$missing_idx) else as.integer(-1)
      ui <- if (is.finite(bins$unknown_idx)) as.integer(bins$unknown_idx) else NULL

      delta <- tryCatch({
        d <- py_update_term(ebm_post, as.integer(bins$term_index),
                            reticulate::r_to_py(y_new),
                            mean_preserve = TRUE,
                            preserve_missing = !isTRUE(EL_UPDATE_MISSING_BIN),
                            missing_index = mi,
                            preserve_unknown = !isTRUE(EL_UPDATE_UNKNOWN_BIN),
                            unknown_index = ui,
                            shift_bagged_scores = TRUE,
                            shift_bagged_intercept = TRUE)
        suppressWarnings(as.numeric(reticulate::py_to_r(d)))
      }, error = function(e) NA_real_)

      if (is.finite(delta)) {
        shift_total <- shift_total + delta
        n_smoothed <- n_smoothed + 1L
        cat(sprintf("      [SMOOTH] %-25s spar=%.2f shift=%+.6f %s\n",
                    feat, feat_spar, delta,
                    ifelse(sm$fit_ok, "", "[fit_failed]")))
      }
    }

    # --- Apply Red action ---
    red_action <- pol$red  # "zero" or "smooth"

    for (feat in LOCKED_RED) {
      bins <- tryCatch(
        extract_bins_with_diagnostics(ebm_post, feat,
                                       strict_continuous_cuts = TRUE, verbose = FALSE),
        error = function(e) list(success = FALSE, error_msg = e$message)
      )
      if (!bins$success) { cat(sprintf("      [Skip] %s: %s\n", feat, bins$error_msg)); next }

      if (red_action == "zero") {
        y_new <- rep(0.0, bins$n_scores)
        mi <- if (is.finite(bins$missing_idx)) as.integer(bins$missing_idx) else as.integer(-1)
        ui <- if (is.finite(bins$unknown_idx)) as.integer(bins$unknown_idx) else NULL

        delta <- tryCatch({
          d <- py_update_term(ebm_post, as.integer(bins$term_index),
                              reticulate::r_to_py(y_new),
                              mean_preserve = TRUE,
                              preserve_missing = FALSE, missing_index = mi,
                              preserve_unknown = FALSE, unknown_index = ui,
                              shift_bagged_scores = TRUE,
                              shift_bagged_intercept = TRUE)
          suppressWarnings(as.numeric(reticulate::py_to_r(d)))
        }, error = function(e) NA_real_)

        if (is.finite(delta)) {
          shift_total <- shift_total + delta
          n_zeroed <- n_zeroed + 1L
          cat(sprintf("      [ZEROED] %-25s shift=%+.6f\n", feat, delta))
        }

      } else if (red_action == "smooth") {
        orig_bins <- original_bins_cache[[feat]]
        if (is.null(orig_bins) || !orig_bins$success) {
          cat(sprintf("      [Skip] %s: not in original cache\n", feat)); next
        }
        original_scores <- orig_bins$dt$score
        x_vec  <- orig_bins$dt$x
        bin_w  <- if ("bin_weight" %in% names(orig_bins$dt)) orig_bins$dt$bin_weight else rep(1, length(x_vec))

        phys_range <- NULL
        if (feat %in% names(PHYS_RANGES)) phys_range <- PHYS_RANGES[[feat]]
        phys_range_use <- phys_range
        update_only_phys_use <- isTRUE(EL_UPDATE_ONLY_PHYS_RANGE)
        if (isTRUE(bins$x_is_index)) {
          phys_range_use <- NULL
          update_only_phys_use <- FALSE
        }

        feat_spar <- get_locked_spar(feat)
        sm <- smooth_scores_spline_fixed(
          x = x_vec, y_old = original_scores, w = bin_w,
          spar = feat_spar,
          phys_range = phys_range_use,
          update_only_phys = update_only_phys_use,
          clip_to_old = isTRUE(EL_CLIP_TO_ORIGINAL_RANGE),
          clip_margin = EL_CLIP_MARGIN
        )
        y_new_sm <- sm$y_new
        idx_no_x <- which(!is.finite(x_vec))
        if (length(idx_no_x) > 0) y_new_sm[idx_no_x] <- original_scores[idx_no_x]
        if (!isTRUE(EL_UPDATE_MISSING_BIN)) {
          idx_m <- which(bins$dt$is_missing)
          if (length(idx_m) > 0) y_new_sm[idx_m] <- original_scores[idx_m]
        }
        if (!isTRUE(EL_UPDATE_UNKNOWN_BIN)) {
          idx_u <- which(bins$dt$is_unknown)
          if (length(idx_u) > 0) y_new_sm[idx_u] <- original_scores[idx_u]
        }

        mi <- if (is.finite(bins$missing_idx)) as.integer(bins$missing_idx) else as.integer(-1)
        ui <- if (is.finite(bins$unknown_idx)) as.integer(bins$unknown_idx) else NULL

        delta <- tryCatch({
          d <- py_update_term(ebm_post, as.integer(bins$term_index),
                              reticulate::r_to_py(y_new_sm),
                              mean_preserve = TRUE,
                              preserve_missing = !isTRUE(EL_UPDATE_MISSING_BIN),
                              missing_index = mi,
                              preserve_unknown = !isTRUE(EL_UPDATE_UNKNOWN_BIN),
                              unknown_index = ui,
                              shift_bagged_scores = TRUE,
                              shift_bagged_intercept = TRUE)
          suppressWarnings(as.numeric(reticulate::py_to_r(d)))
        }, error = function(e) NA_real_)

        if (is.finite(delta)) {
          shift_total <- shift_total + delta
          n_smoothed <- n_smoothed + 1L
          cat(sprintf("      [SMOOTH] %-25s spar=%.2f shift=%+.6f\n",
                      feat, feat_spar, delta))
        }
      }
    }

    # --- p_post: after governance ---
    p_post <- predict_proba_pos_r(ebm_post, X_L)

    # --- Mean-preserve verification ---
    mp <- verify_mean_preserve_consistency(ebm_pre, ebm_post, shift_total)
    cat(sprintf("      [MP] %s (disc=%.2e)\n", mp$status, mp$discrepancy))

    # --- Performance metrics ---
    auc_pre  <- auc_pair(y_L, p_pre);  auc_post <- auc_pair(y_L, p_post)
    b_pre    <- brier_score(y_L, p_pre); b_post <- brier_score(y_L, p_post)

    perf_row <- data.table(
      iter_tag = iter_tag, external = ext_site, policy = pol$id,
      n_smoothed = n_smoothed, n_zeroed = n_zeroed,
      AUPRC_pre = auc_pre["auprc"], AUPRC_post = auc_post["auprc"],
      d_AUPRC = auc_post["auprc"] - auc_pre["auprc"],
      AUROC_pre = auc_pre["auroc"], AUROC_post = auc_post["auroc"],
      d_AUROC = auc_post["auroc"] - auc_pre["auroc"],
      Brier_pre = b_pre, Brier_post = b_post,
      d_Brier = b_post - b_pre,
      mp_status = mp$status, mp_disc = mp$discrepancy,
      threshold_E = thr_E, threshold_full = thr_full
    )

    # --- Threshold-based metrics (E-only locked threshold) ---
    bm_pre_E  <- bin_metrics(y_L, p_pre,  thr_E)
    bm_post_E <- bin_metrics(y_L, p_post, thr_E)
    perf_row[, `:=`(
      Sens_pre_E = bm_pre_E$Sensitivity, Sens_post_E = bm_post_E$Sensitivity,
      Spec_pre_E = bm_pre_E$Specificity, Spec_post_E = bm_post_E$Specificity,
      PPV_pre_E  = bm_pre_E$PPV,         PPV_post_E  = bm_post_E$PPV,
      NPV_pre_E  = bm_pre_E$NPV,         NPV_post_E  = bm_post_E$NPV,
      F1_pre_E   = bm_pre_E$F1,          F1_post_E   = bm_post_E$F1
    )]

    # --- Threshold-based metrics (full-IECV threshold, sensitivity) ---
    bm_pre_F  <- bin_metrics(y_L, p_pre,  thr_full)
    bm_post_F <- bin_metrics(y_L, p_post, thr_full)
    perf_row[, `:=`(
      Sens_pre_full = bm_pre_F$Sensitivity, Sens_post_full = bm_post_F$Sensitivity,
      Spec_pre_full = bm_pre_F$Specificity, Spec_post_full = bm_post_F$Specificity
    )]

    # --- Calibration (key secondary endpoint) ---
    safe_logit <- function(p) log(pmax(p, 1e-10) / pmax(1 - p, 1e-10))

    cal_post <- tryCatch({
      lp <- safe_logit(p_post)
      cal_fit <- glm(y_L ~ lp, family = binomial(link = "logit"))
      cal_slope <- unname(coef(cal_fit)[2])
      cal_int   <- unname(coef(cal_fit)[1])
      cal_OE    <- mean(y_L) / mean(p_post)
      list(slope = cal_slope, intercept = cal_int, OE = cal_OE)
    }, error = function(e) list(slope = NA_real_, intercept = NA_real_, OE = NA_real_))

    cal_pre <- tryCatch({
      lp <- safe_logit(p_pre)
      cal_fit <- glm(y_L ~ lp, family = binomial(link = "logit"))
      list(slope = unname(coef(cal_fit)[2]),
           intercept = unname(coef(cal_fit)[1]),
           OE = mean(y_L) / mean(p_pre))
    }, error = function(e) list(slope = NA_real_, intercept = NA_real_, OE = NA_real_))

    perf_row[, `:=`(
      cal_slope_pre = cal_pre$slope, cal_slope_post = cal_post$slope,
      cal_int_pre   = cal_pre$intercept, cal_int_post = cal_post$intercept,
      cal_OE_pre    = cal_pre$OE, cal_OE_post = cal_post$OE
    )]

    # Flag calibration NA
    if (is.na(cal_post$slope) || is.na(cal_pre$slope)) {
      cat(sprintf("      [WARN] Calibration returned NA for %s\n", iter_tag))
    }

    # --- DCA (secondary endpoint, 5-20%) ---
    dca_range <- seq(0.05, 0.20, by = 0.01)
    dca_post <- tryCatch(dca_curve(y_L, p_post, dca_range), error = function(e) NULL)
    dca_pre  <- tryCatch(dca_curve(y_L, p_pre,  dca_range), error = function(e) NULL)
    if (!is.null(dca_post) && !is.null(dca_pre) && pol$primary) {
      dca_dt <- data.table(
        iter_tag = iter_tag, external = ext_site,
        threshold = dca_range,
        NB_pre  = dca_pre$NB_model,
        NB_post = dca_post$NB_model,
        d_NB    = dca_post$NB_model - dca_pre$NB_model
      )
      fst::write_fst(dca_dt, file.path(out_iter,
                     sprintf("dca_L_%s.fst", pol$id)))
    }

    all_confirm <- rbind(all_confirm, perf_row, fill = TRUE)
    all_mp <- rbind(all_mp, data.table(iter_tag = iter_tag, policy = pol$id,
                                        status = mp$status, disc = mp$discrepancy))

    cat(sprintf("      AUROC: %.4f → %.4f (Δ=%+.4f)\n",
                perf_row$AUROC_pre, perf_row$AUROC_post, perf_row$d_AUROC))
    cat(sprintf("      AUPRC: %.4f → %.4f (Δ=%+.4f)\n",
                perf_row$AUPRC_pre, perf_row$AUPRC_post, perf_row$d_AUPRC))

    # --- NRI / IDI (for primary policy only) ---
    if (pol$primary) {
      classify <- function(p) {
        ifelse(p < NRI_THRESHOLDS[1], 1L,
               ifelse(p < NRI_THRESHOLDS[2], 2L, 3L))
      }
      cat_pre  <- classify(p_pre)
      cat_post <- classify(p_post)

      events    <- which(y_L == 1)
      nonevents <- which(y_L == 0)

      nri_up_event   <- mean(cat_post[events] > cat_pre[events]) -
                          mean(cat_post[events] < cat_pre[events])
      nri_up_nonevent <- mean(cat_post[nonevents] < cat_pre[nonevents]) -
                           mean(cat_post[nonevents] > cat_pre[nonevents])
      cat_nri <- nri_up_event + nri_up_nonevent

      reclass_rate <- mean(cat_pre != cat_post) * 100

      idi <- mean(p_post[events]) - mean(p_pre[events]) -
             (mean(p_post[nonevents]) - mean(p_pre[nonevents]))

      nri_row <- data.table(
        iter_tag = iter_tag, external = ext_site,
        cat_NRI = cat_nri, NRI_event = nri_up_event,
        NRI_nonevent = nri_up_nonevent,
        reclass_rate_pct = reclass_rate, IDI = idi
      )
      all_nri <- rbind(all_nri, nri_row, fill = TRUE)
      cat(sprintf("      NRI: %+.4f (event=%+.4f, nonevent=%+.4f) | reclass=%.2f%%\n",
                  cat_nri, nri_up_event, nri_up_nonevent, reclass_rate))
    }

    # Save predictions
    pred_dt <- data.table(
      patient_id = dt_L[[ID_COL]], cluster_id = dt_L[[CLUSTER_COL]],
      outcome = y_L, p_pre = p_pre, p_post = p_post,
      iter_tag = iter_tag, policy = pol$id
    )
    fst::write_fst(pred_dt, file.path(out_iter,
                   sprintf("predictions_L_%s.fst", pol$id)))

    rm(ebm_pre, ebm_post)
    if (FORCE_GC_EACH_ITER) gc(verbose = FALSE)
  }
  cat("\n")
}

# =============================================================================
# SECTION D: Bootstrap CI (patient-level, primary policy only)
# =============================================================================

cat("[Section D] Bootstrap CI on L ...\n\n")

boot_results <- list()

for (it in ITERATIONS) {
  iter_tag <- sprintf("Iter%d_External_%s_Seed%d", it$iter, it$external, it$seed_id)
  ext_site <- it$external
  out_iter <- file.path(EL_DIR_CONFIRM, iter_tag)

  pred_file <- list.files(out_iter, pattern = "predictions_L_A\\.fst$", full.names = TRUE)
  if (length(pred_file) == 0) next

  pred <- read_fst_dt(pred_file[1])
  set.seed(EL_BOOT_SEED + it$iter)

  boot_d_auprc <- boot_d_auroc <- rep(NA_real_, EL_BOOT_N)
  clusters <- unique(pred$cluster_id)

  for (b in seq_len(EL_BOOT_N)) {
    boot_clusters <- sample(clusters, length(clusters), replace = TRUE)

    boot_key <- data.table(cluster_id = boot_clusters)
    boot_rows <- pred[boot_key, on = "cluster_id", allow.cartesian = TRUE, nomatch = NULL]

    if (length(unique(boot_rows$outcome)) < 2) next

    auc_pre_b  <- tryCatch(auc_pair(boot_rows$outcome, boot_rows$p_pre),
                            error = function(e) c(auroc = NA, auprc = NA))
    auc_post_b <- tryCatch(auc_pair(boot_rows$outcome, boot_rows$p_post),
                            error = function(e) c(auroc = NA, auprc = NA))

    boot_d_auprc[b] <- auc_post_b["auprc"] - auc_pre_b["auprc"]
    boot_d_auroc[b] <- auc_post_b["auroc"] - auc_pre_b["auroc"]
  }

  valid_n <- sum(is.finite(boot_d_auprc))
  ci_auprc <- quantile(boot_d_auprc, c(0.025, 0.975), na.rm = TRUE)
  ci_auroc <- quantile(boot_d_auroc, c(0.025, 0.975), na.rm = TRUE)

  boot_results[[iter_tag]] <- data.table(
    iter_tag = iter_tag, external = ext_site,
    d_AUPRC_lo = ci_auprc[1], d_AUPRC_hi = ci_auprc[2],
    d_AUROC_lo = ci_auroc[1], d_AUROC_hi = ci_auroc[2],
    boot_valid = valid_n
  )

  cat(sprintf("  %s: ΔAUPRC 95%% CI [%+.4f, %+.4f] | ΔAUROC [%+.4f, %+.4f]\n",
              iter_tag, ci_auprc[1], ci_auprc[2], ci_auroc[1], ci_auroc[2]))
}

boot_dt <- rbindlist(boot_results, fill = TRUE)

# =============================================================================
# SECTION E: Save all results
# =============================================================================

fwrite(all_confirm, file.path(EL_DIR_CONFIRM, "EL_Confirmatory_Performance.csv"))
fwrite(all_mp,      file.path(EL_DIR_CONFIRM, "EL_Mean_Preserve_Verification.csv"))
fwrite(all_nri,     file.path(EL_DIR_CONFIRM, "EL_NRI_IDI.csv"))
fwrite(boot_dt,     file.path(EL_DIR_CONFIRM, "EL_Bootstrap_CI.csv"))

# =============================================================================
# SECTION F: Success Criterion Evaluation
# =============================================================================

cat("\n\n  ═══ CONFIRMATORY RESULTS (Policy A — Primary) ═══\n\n")

pol_A <- all_confirm[policy == "A"]

if (nrow(pol_A) > 0) {
  primary <- pol_A[grepl("Seed1", iter_tag)]
  cat("  Per-center (primary seed):\n")
  for (i in seq_len(nrow(primary))) {
    r <- primary[i]
    harm_flag <- ifelse(r$d_AUPRC < EL_HARM_MARGIN, " *** HARM ***", "")
    cat(sprintf("    %s: ΔAUPRC=%+.4f  ΔAUROC=%+.4f%s\n",
                r$external, r$d_AUPRC, r$d_AUROC, harm_flag))
  }

  mean_d_auprc <- mean(primary$d_AUPRC)
  mean_d_auroc <- mean(primary$d_AUROC)
  any_harm <- any(primary$d_AUPRC < EL_HARM_MARGIN)

  cat(sprintf("\n  Center-unweighted mean: ΔAUPRC=%+.4f  ΔAUROC=%+.4f\n",
              mean_d_auprc, mean_d_auroc))
  cat(sprintf("  Non-inferiority (ΔAUPRC ≥ 0): %s\n",
              ifelse(mean_d_auprc >= EL_NI_MARGIN, "PASS", "FAIL")))
  cat(sprintf("  No material harm (all > %.3f): %s\n",
              EL_HARM_MARGIN, ifelse(!any_harm, "PASS", "FAIL")))
  cat(sprintf("  Mean-preserve: %d/%d pass\n",
              sum(all_mp[policy == "A", status == "PASS"]),
              nrow(all_mp[policy == "A"])))
}

# Policy B comparison
pol_B <- all_confirm[policy == "B"]
if (nrow(pol_B) > 0) {
  cat("\n  ═══ SENSITIVITY ANALYSIS (Policy B — Smooth All) ═══\n\n")
  primary_B <- pol_B[grepl("Seed1", iter_tag)]
  for (i in seq_len(nrow(primary_B))) {
    r <- primary_B[i]
    cat(sprintf("    %s: ΔAUPRC=%+.4f  ΔAUROC=%+.4f\n",
                r$external, r$d_AUPRC, r$d_AUROC))
  }
  cat(sprintf("  Mean ΔAUPRC: A=%+.4f vs B=%+.4f → %s\n",
              mean(primary$d_AUPRC), mean(primary_B$d_AUPRC),
              ifelse(mean(primary$d_AUPRC) >= mean(primary_B$d_AUPRC),
                     "Policy A ≥ B (confirms primary)", "Policy B > A")))
}

cat(sprintf("\n  Output: %s\n", EL_DIR_CONFIRM))
cat("\n  Module 04f complete.\n\n")
