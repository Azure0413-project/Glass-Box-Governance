# =============================================================================
# 04b_temporal_holdout.R — Within-Center Temporal Hold-out Validation
# IDH EBM Governance Reproducibility Pipeline
#
# v3.0 update: Best HP are now auto-detected from 03_iecv_ebm.R v3.0
#   E-only IECV output. The within-center temporal split is applied
#   AFTER restricting to E-cohort, ensuring consistency with the E-only
#   pipeline — HP and training data both originate from the same epoch.
#
# Outputs:
#   Supp Table S6 — Within-center temporal hold-out validation performance
#
# Source: consolidated from original analysis scripts
# Design: each center is split by Session_Date at 75/25
#       Train EBM using the center's best HP from IECV and evaluate
#       No grid search re-run, no bootstrap
#
# Methods → "A supplementary within-center temporal hold-out analysis was
#  also conducted within each center's E-cohort (earlier 75% of overall
#  sessions): the earlier 75% of E-cohort sessions were used for training
#  and the later 25% for testing"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 03_iecv_ebm.R (needs best HP records)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(reticulate)
  library(precrec)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")
if (!exists("init_python_full")) src("R/00_utils_python.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 04b: Within-Center Temporal Hold-out\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- Initialize Python environment ---
init_python_full()

# --- Output directory ---
OUTPUT_DIR <- file.path(RUN_DIR, "temporal_validation")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Section 1: Best HP from IECV (per center, seed_id=1)
# Methods -> "Use IECV best HP to train EBM"
# These are the best hyperparameters selected when each center was held out.
# For temporal validation we use HP from the iteration where that center
# was held out (i.e., trained on the OTHER two centers), seed_id=1.
#
# In practice these are read from the IECV artifacts. Here we hardcode
# the values from the actual IECV run for exact reproducibility.
# =============================================================================

# NOTE: Update these values from your actual 03_iecv_ebm.R output
# (found in <RUN_DIR>/Iter*_artifacts.rds → best_hp)
TEMPORAL_HP <- list(
  TN = list(learning_rate = 0.005, interactions = 0L, max_bins = 128L, max_leaves = 2L,
            min_samples_leaf = 600L, smoothing_rounds = 150L, pos_multiplier = 6.0),
  D6 = list(learning_rate = 0.005, interactions = 0L, max_bins = 128L, max_leaves = 2L,
            min_samples_leaf = 600L, smoothing_rounds = 150L, pos_multiplier = 6.0),
  CY = list(learning_rate = 0.005, interactions = 0L, max_bins = 128L, max_leaves = 2L,
            min_samples_leaf = 600L, smoothing_rounds = 75L,  pos_multiplier = 6.0)
)

TEMPORAL_RATIO  <- 0.75   # 75% train / 25% test
TEMPORAL_SEED   <- 2024L  # seed for EBM training

cat(sprintf("  Temporal split: %.0f%% train / %.0f%% test\n",
            TEMPORAL_RATIO * 100, (1 - TEMPORAL_RATIO) * 100))
cat(sprintf("  Output: %s\n\n", OUTPUT_DIR))

# =============================================================================
# Section 2: Auto-detect IECV HP (if artifacts exist)
# =============================================================================

# Try to read actual HP from IECV artifacts
for (site in c("TN", "D6", "CY")) {
  # Look for seed_id=1 iteration where this site was held out
  iter_pattern <- sprintf("Iter.*External_%s_Seed1$", site)
  iter_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
  iter_dirs <- iter_dirs[grepl(iter_pattern, basename(iter_dirs))]

  if (length(iter_dirs) > 0) {
    art_files <- list.files(iter_dirs[1], pattern = "artifacts\\.rds$", full.names = TRUE)
    if (length(art_files) > 0) {
      art <- readRDS(art_files[length(art_files)])
      hp  <- art$best_hp
      if (!is.null(hp)) {
        # 03_iecv_ebm.R stores the full HP list in best_hp (including fixed +
        # tuned fields + random_state/outer_bags).  Extract the 7 core training
        # fields; use modifyList so any missing field keeps its default.
        core_fields <- c("learning_rate", "interactions", "max_bins",
                         "max_leaves", "min_samples_leaf", "smoothing_rounds",
                         "pos_multiplier")
        hp_patch <- hp[intersect(core_fields, names(hp))]
        # Ensure integer types for EBM Python API
        for (fn in intersect(c("interactions", "max_bins", "max_leaves",
                                "min_samples_leaf", "smoothing_rounds"), names(hp_patch)))
          hp_patch[[fn]] <- as.integer(hp_patch[[fn]])
        if (length(hp_patch) > 0) {
          TEMPORAL_HP[[site]] <- utils::modifyList(TEMPORAL_HP[[site]], hp_patch)
        }
        missing_hp <- setdiff(core_fields, names(hp))
        if (length(missing_hp) > 0)
          cat(sprintf("  %s: [WARN] Artifact best_hp missing fields: %s; keeping defaults\n",
                      site, paste(missing_hp, collapse = ", ")))
        cat(sprintf("  %s: HP loaded from IECV artifacts (leaf=%d, smooth=%d)\n",
                     site, TEMPORAL_HP[[site]]$min_samples_leaf, TEMPORAL_HP[[site]]$smoothing_rounds))
      }
    }
  }
}

# =============================================================================
# Section 3: Main Loop — Temporal validation per center
# =============================================================================

cat("\n[Step 1] Running temporal hold-out for each center ...\n\n")

results_list <- list()

for (site in c("TN", "D6", "CY")) {
  cat(sprintf("  === %s ===\n", site))

  # Load and preprocess
  dt <- prep_site_dt(read_fst_dt(SITE_FILES[[site]]), site)

  # Temporal split by configured session-date column
  date_col <- SESSION_DATE_COL
  if (!(date_col %in% names(dt))) {
    cat(sprintf("    [SKIP] %s: no %s column\n\n", site, date_col))
    next
  }
  dt[, (date_col) := as.Date(get(date_col))]

  # --- Fix: restrict to E-cohort BEFORE within-center temporal split ---
  # HP were tuned on E-only data (03_iecv_ebm.R v3.0). Using all sessions

  # (including L-period) for temporal train would mix data that the HP
  # selection never saw, creating an inconsistent evaluation context.
  # By filtering to E-only first, the within-center 75/25 temporal split
  # tests temporal generalization strictly within the derivation epoch.
  dt <- apply_el_split(dt, date_col, E_SPLIT_QUANTILE, COHORT_COL)
  n_before <- nrow(dt)
  dt <- dt[get(COHORT_COL) == "E"]
  if (n_before == 0L || nrow(dt) == 0L) {
    cat(sprintf("    [SKIP] %s: E-cohort is empty after temporal split\n\n", site))
    next
  }
  cat(sprintf("    E-cohort filter: %s → %s sessions (%.1f%%)\n",
              format(n_before, big.mark = ","),
              format(nrow(dt), big.mark = ","),
              nrow(dt) / n_before * 100))

  data.table::setorderv(dt, date_col)

  # Within-E temporal split (quantile-based)
  cutoff   <- as.Date(quantile(dt[[date_col]], probs = TEMPORAL_RATIO, type = 1, na.rm = TRUE))
  dt_train <- dt[get(date_col) <= cutoff]
  dt_test  <- dt[get(date_col) >  cutoff]
  if (nrow(dt_train) == 0L || nrow(dt_test) == 0L) {
    cat(sprintf("    [SKIP] %s: temporal split produced train=%d, test=%d (cutoff=%s)\n\n",
                site, nrow(dt_train), nrow(dt_test), as.character(cutoff)))
    next
  }

  cat(sprintf("    Train: %s sessions (%s to %s)\n",
              format(nrow(dt_train), big.mark = ","),
              min(dt_train[[date_col]]), max(dt_train[[date_col]])))
  cat(sprintf("    Test:  %s sessions (%s to %s)\n",
              format(nrow(dt_test), big.mark = ","),
              min(dt_test[[date_col]]), max(dt_test[[date_col]])))

  # Get HP for this center
  hp <- TEMPORAL_HP[[site]]
  hp$random_state <- TEMPORAL_SEED
  hp$outer_bags   <- OUTER_BAGS_FINAL

  # Train EBM
  X_train <- as.data.frame(dt_train[, ..ALL_PREDICTORS])
  y_train <- as.integer(dt_train[[TARGET_COL]])
  w_train <- make_sample_weight(y_train, hp$pos_multiplier)

  # OOF threshold tuning (5-fold grouped CV on training set)
  # Matches 03_iecv_ebm.R methodology — avoids resubstitution bias
  hp_oof <- hp
  hp_oof$outer_bags <- OUTER_BAGS_OOF
  fold_seed_temp <- TEMPORAL_SEED + 1000L
  folds_train <- make_group_folds(dt_train[[CLUSTER_COL]], k = K_INNER, seed = fold_seed_temp)
  p_oof_train <- rep(NA_real_, nrow(dt_train))

  for (f in seq_len(K_INNER)) {
    tr_idx <- which(folds_train != f); va_idx <- which(folds_train == f)
    X_tr_f <- as.data.frame(dt_train[tr_idx, ..ALL_PREDICTORS])
    y_tr_f <- y_train[tr_idx]
    w_tr_f <- make_sample_weight(y_tr_f, hp$pos_multiplier)
    ebm_f  <- build_ebm(hp_oof, EBM_N_JOBS, interpret_glassbox)
    ebm_f$fit(X_tr_f, y_tr_f, sample_weight = w_tr_f)
    X_va_f <- as.data.frame(dt_train[va_idx, ..ALL_PREDICTORS])
    p_oof_train[va_idx] <- predict_proba_pos_r(ebm_f, X_va_f)
  }
  thr_obj  <- tune_threshold_constraint(y_train, p_oof_train, THR_GRID, SENS_TARGET)
  best_thr <- thr_obj$threshold
  cat(sprintf("    OOF threshold: %.4f (met_sens80=%s)\n", best_thr, thr_obj$met_constraint))

  # Train final model on full training set
  ebm <- build_ebm(hp, EBM_N_JOBS, interpret_glassbox)
  ebm$fit(X_train, y_train, sample_weight = w_train)

  # Predict on temporal test
  X_test <- as.data.frame(dt_test[, ..ALL_PREDICTORS])
  y_test <- as.integer(dt_test[[TARGET_COL]])
  p_test <- predict_proba_pos_r(ebm, X_test)

  # Evaluate
  aucs <- auc_pair(y_test, p_test)
  bm   <- bin_metrics(y_test, p_test, best_thr)

  cat(sprintf("    AUROC=%.4f | AUPRC=%.4f | Sens=%.2f%% | Spec=%.2f%%\n",
              aucs["auroc"], aucs["auprc"],
              bm$Sensitivity * 100, bm$Specificity * 100))

  results_list[[site]] <- data.table(
    Center       = site,
    Train_n      = nrow(dt_train),
    Test_n       = nrow(dt_test),
    Train_prev   = mean(y_train),
    Test_prev    = mean(y_test),
    Train_period = sprintf("%s to %s", min(dt_train[[date_col]]), max(dt_train[[date_col]])),
    Test_period  = sprintf("%s to %s", min(dt_test[[date_col]]), max(dt_test[[date_col]])),
    Threshold    = best_thr,
    AUROC        = aucs["auroc"],
    AUPRC        = aucs["auprc"],
    Brier        = brier_score(y_test, p_test),
    Sensitivity  = bm$Sensitivity,
    Specificity  = bm$Specificity,
    PPV          = bm$PPV,
    NPV          = bm$NPV,
    F1           = bm$F1,
    HP_leaf      = hp$min_samples_leaf,
    HP_smooth    = hp$smoothing_rounds,
    HP_pos_mult  = hp$pos_multiplier
  )

  # Save predictions
  pred_dt <- data.table(
    Patient_ID = dt_test[[ID_COL]],
    Session_Date = dt_test[[date_col]],
    y_true = y_test, y_prob = p_test,
    y_hat = as.integer(p_test >= best_thr)
  )
  fst::write_fst(pred_dt, file.path(OUTPUT_DIR, sprintf("%s_temporal_predictions.fst", site)))

  cat("\n")
}

# =============================================================================
# Section 4: Output — Supp Table S6
# =============================================================================

supp_s6 <- rbindlist(results_list, fill = TRUE)
fwrite(supp_s6, file.path(OUTPUT_DIR, "Supp_S6_Temporal_Holdout.csv"))

cat("\n===== Supp Table S6: Within-Center Temporal Hold-out =====\n")
print(supp_s6[, .(Center, Train_n, Test_n,
                   AUROC = round(AUROC, 4), AUPRC = round(AUPRC, 4),
                   Sensitivity = round(Sensitivity, 4),
                   Specificity = round(Specificity, 4))])

cat("\n[Verify] Manuscript expects comparable AUROC across all three centers.\n")
cat(sprintf("  Got: TN=%.3f, D6=%.3f, CY=%.3f\n",
            supp_s6[Center == "TN", AUROC],
            supp_s6[Center == "D6", AUROC],
            supp_s6[Center == "CY", AUROC]))

cat(sprintf("\nOutput: %s\n", OUTPUT_DIR))
cat("\n  Module 04b complete.\n\n")
