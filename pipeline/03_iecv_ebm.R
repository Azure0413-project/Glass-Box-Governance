# =============================================================================
# 03_iecv_ebm.R — E-only IECV-EBM Training & External Validation
# IDH EBM Governance Reproducibility Pipeline
#
# Design change (v3.0): Pure E-only pipeline
#   All hyperparameter tuning, model training, and external validation use
#   exclusively E-cohort (earlier 75%) sessions. L-cohort (later 25%) is
#   completely untouched, ensuring a COMPLETE protocol lock for downstream
#   confirmation (script 12+). This eliminates the "partial protocol lock"
#   limitation of the previous full-cohort hyperparameter inheritance design.
#
# Outputs:
#   Table 2    — E-only external validation (AUPRC, AUROC, calibration, operating chars)
#   Supp S2    — Aggregated feature importance across iterations
#   Fig 4c     — Decision curve analysis
#   Per-iteration: trained EBM models (.joblib), predictions (.fst),
#                  shape functions, calibration, DCA, bootstrap CI
#   NEW: E/L split metadata (per-site session counts, split dates)
#
# Source: consolidated from original analysis scripts
# Design: 3 centers × 2 seeds = 6 parallel workers
#   Each worker: E-only data → grid search (60% subsample, 5-fold grouped CV) →
#                refine on full E-only dev → train final model → E-only external validation
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
# Hardware: recommended 6 cores + 30 threads (6 workers x 5 threads each)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(reticulate)
  library(precrec)
  library(ggplot2)
  library(gridExtra)
  library(parallel)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")
if (!exists("init_python_full")) src("R/00_utils_python.R")

cat("\n")
cat(paste(rep("#", 70), collapse = ""), "\n")
cat("# Module 03: E-only IECV-EBM v3.0\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n")

# --- Output directory ---
OUT_ROOT <- file.path(BASE_DIR, paste0("Run_", RUN_TS))
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# Keep downstream modules aligned to the current IECV run when this script is
# sourced inside run_all.R or followed immediately by Modules 04b–11.
RUN_DIR <- OUT_ROOT

# Additional config flags used in original code
DO_PHYS_CHECK <- TRUE
DO_SCALE      <- FALSE
B_BOOT        <- B_BOOTSTRAP
BOOT_SEED     <- BOOTSTRAP_SEED
GRID          <- as.data.table(EBM_GRID)
ENFORCE_SENS_CONSTRAINT <- TRUE
SITE_ORDER    <- c("TN", "D6", "CY")

cat(sprintf("  OUT_ROOT: %s\n", OUT_ROOT))
cat(sprintf("  Grid: %d combinations (3 min_samples_leaf × 4 smoothing_rounds × 3 pos_multiplier)\n", nrow(GRID)))
cat(sprintf("  Workers: %d × %d threads = %d total\n",
            N_WORKERS, N_THREADS_PER_WORKER, N_WORKERS * N_THREADS_PER_WORKER))

# =============================================================================
# Phase 6: Data preprocessing function (worker-compatible version)
# (accepts explicit arguments so it can be exported to parallel workers)
# =============================================================================

prep_site_dt_worker <- function(dt_raw, site_name, id_col, cluster_col, target_col,
                                binary_cols, cont_cols, compat_rename_map, phys_ranges,
                                do_phys_check) {
  dt <- data.table::copy(dt_raw)
  dt <- apply_compat_rename(dt, compat_rename_map)
  dt[, site := site_name]
  if (!(target_col %in% names(dt))) stop(sprintf("[%s] TARGET_COL not found: %s", site_name, target_col))
  dt[, (target_col) := to_target01(get(target_col))]
  dt <- dt[!is.na(get(target_col))]
  if (!(id_col %in% names(dt))) stop(sprintf("[%s] ID_COL not found: %s", site_name, id_col))
  dt[, (cluster_col) := paste0(site_name, "__", as.character(get(id_col)))]
  for (b in binary_cols) if (b %in% names(dt)) dt[, (b) := to_binary01(get(b))]
  for (cn in cont_cols) if (cn %in% names(dt)) dt[, (cn) := as.numeric(get(cn))]
  if (isTRUE(do_phys_check)) dt <- apply_physio_ranges(dt, phys_ranges)
  dt
}

# =============================================================================
# Phase 6.1: E/L Temporal Split
# Uses apply_el_split() from 00_utils_r.R (Section 19).
# Exported to parallel workers via clusterExport in Phase 8.
# =============================================================================

# =============================================================================
# Phase 7: Worker execution function
# Each worker handles one IECV iteration independently.
# =============================================================================

run_iteration_worker <- function(it) {

  suppressPackageStartupMessages({
    library(data.table); library(fst); library(reticulate)
    library(precrec); library(ggplot2); library(gridExtra)
  })

  Sys.setenv(
    OMP_NUM_THREADS      = as.character(N_THREADS_PER_WORKER),
    MKL_NUM_THREADS      = as.character(N_THREADS_PER_WORKER),
    OPENBLAS_NUM_THREADS  = as.character(N_THREADS_PER_WORKER)
  )
  reticulate::use_virtualenv(VENV_PATH, required = TRUE)
  glassbox <- reticulate::import("interpret.glassbox", convert = TRUE)
  joblib   <- reticulate::import("joblib",             convert = TRUE)

  iter_id     <- it$iter
  ext_site    <- it$external
  dev_sites   <- it$dev
  seed_id     <- it$seed_id
  master_seed <- it$master_seed

  iter_tag <- sprintf("Iter%d_External_%s_Seed%d", iter_id, ext_site, seed_id)
  out_dir  <- file.path(OUT_ROOT, iter_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  iter_start <- Sys.time()
  iter_ts    <- format(iter_start, "%Y_%m_%d_%H_%M")

  log_file <- file.path(out_dir, sprintf("%s_log.txt", iter_ts))
  log_con  <- file(log_file, open = "wt")
  log_msg  <- function(msg) { cat(msg, "\n", file = log_con); flush(log_con) }

  log_msg(sprintf("[%s] Started at %s", iter_tag, iter_ts))
  log_msg(sprintf("[%s] Master seed: %d", iter_tag, master_seed))

  tryCatch({

    # ---- Step 1: Data Loading & E/L Split ----
    log_msg("\n[Step 1] Data Loading & E/L Temporal Split")

    dt_dev1 <- prep_site_dt_worker(read_fst_dt(SITE_FILES[[dev_sites[1]]]), dev_sites[1],
                                    ID_COL, CLUSTER_COL, TARGET_COL, BINARY_COLS, CONT_COLS,
                                    COMPAT_RENAME_MAP, PHYS_RANGES, DO_PHYS_CHECK)
    dt_dev2 <- prep_site_dt_worker(read_fst_dt(SITE_FILES[[dev_sites[2]]]), dev_sites[2],
                                    ID_COL, CLUSTER_COL, TARGET_COL, BINARY_COLS, CONT_COLS,
                                    COMPAT_RENAME_MAP, PHYS_RANGES, DO_PHYS_CHECK)
    dt_ext  <- prep_site_dt_worker(read_fst_dt(SITE_FILES[[ext_site]]), ext_site,
                                    ID_COL, CLUSTER_COL, TARGET_COL, BINARY_COLS, CONT_COLS,
                                    COMPAT_RENAME_MAP, PHYS_RANGES, DO_PHYS_CHECK)

    # Apply E/L split to each site
    dt_dev1 <- apply_el_split(dt_dev1, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
    dt_dev2 <- apply_el_split(dt_dev2, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
    dt_ext  <- apply_el_split(dt_ext,  SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)

    # Log full counts before filtering
    n_dev1_all <- nrow(dt_dev1); n_dev1_E <- sum(dt_dev1[[COHORT_COL]] == "E")
    n_dev2_all <- nrow(dt_dev2); n_dev2_E <- sum(dt_dev2[[COHORT_COL]] == "E")
    n_ext_all  <- nrow(dt_ext);  n_ext_E  <- sum(dt_ext[[COHORT_COL]] == "E")

    log_msg(sprintf("  %s: %d total → %d E (%.1f%%)", dev_sites[1], n_dev1_all, n_dev1_E, n_dev1_E/n_dev1_all*100))
    log_msg(sprintf("  %s: %d total → %d E (%.1f%%)", dev_sites[2], n_dev2_all, n_dev2_E, n_dev2_E/n_dev2_all*100))
    log_msg(sprintf("  %s: %d total → %d E (%.1f%%)", ext_site,     n_ext_all,  n_ext_E,  n_ext_E/n_ext_all*100))

    # Save split metadata (for downstream L-cohort confirmation scripts)
    split_meta <- data.table(
      site    = c(dev_sites, ext_site),
      n_total = c(n_dev1_all, n_dev2_all, n_ext_all),
      n_E     = c(n_dev1_E, n_dev2_E, n_ext_E),
      n_L     = c(n_dev1_all - n_dev1_E, n_dev2_all - n_dev2_E, n_ext_all - n_ext_E),
      E_cutdate = c(
        as.character(dt_dev1[get(COHORT_COL) == "E", max(get(SESSION_DATE_COL), na.rm = TRUE)]),
        as.character(dt_dev2[get(COHORT_COL) == "E", max(get(SESSION_DATE_COL), na.rm = TRUE)]),
        as.character(dt_ext[get(COHORT_COL)  == "E", max(get(SESSION_DATE_COL), na.rm = TRUE)])
      )
    )
    fwrite(split_meta, file.path(out_dir, sprintf("%s_EL_split_meta.csv", iter_ts)))

    # ---- Filter to E-cohort ONLY ----
    dt_dev1 <- dt_dev1[get(COHORT_COL) == "E"]
    dt_dev2 <- dt_dev2[get(COHORT_COL) == "E"]
    dt_ext  <- dt_ext[get(COHORT_COL)  == "E"]

    dev_pool <- rbindlist(list(dt_dev1, dt_dev2), use.names = TRUE, fill = TRUE)
    check_required_columns(dev_pool, c(ID_COL, CLUSTER_COL, "site", TARGET_COL, ALL_PREDICTORS), "Dev_Pool_E")
    check_required_columns(dt_ext, c(ID_COL, CLUSTER_COL, "site", TARGET_COL, ALL_PREDICTORS), ext_site)

    log_msg(sprintf("  Dev_Pool (E-only): n=%d, pos=%.2f%%", nrow(dev_pool), mean(dev_pool[[TARGET_COL]]) * 100))
    log_msg(sprintf("  External (E-only): n=%d, pos=%.2f%%", nrow(dt_ext), mean(dt_ext[[TARGET_COL]]) * 100))

    # ---- Step 2: Hyperparameter Optimization (Grid Search on 60% subsample) ----
    log_msg("\n[Step 2] Hyperparameter Optimization")
    log_msg(sprintf("  Grid combinations: %d", nrow(GRID)))

    y_full <- as.integer(dev_pool[[TARGET_COL]])

    set.seed(master_seed)
    cluster_ids      <- unique(dev_pool[[CLUSTER_COL]])
    n_sample         <- ceiling(length(cluster_ids) * GRID_SUBSAMPLE_RATIO)
    sampled_clusters <- sample(cluster_ids, n_sample)
    dt_sub           <- dev_pool[get(CLUSTER_COL) %in% sampled_clusters]
    y_sub            <- as.integer(dt_sub[[TARGET_COL]])

    log_msg(sprintf("  Subsampling: %.0f%% (%d rows)", GRID_SUBSAMPLE_RATIO * 100, nrow(dt_sub)))

    res_list   <- vector("list", nrow(GRID))
    grid_start <- Sys.time()

    for (i in seq_len(nrow(GRID))) {
      hp <- c(EBM_FIXED_HP, as.list(GRID[i]))
      hp$outer_bags   <- OUTER_BAGS_OOF
      hp$random_state <- master_seed

      fold_seed <- master_seed + 1000L
      folds     <- make_group_folds(dt_sub[[CLUSTER_COL]], k = K_INNER, seed = fold_seed)
      p_oof     <- rep(NA_real_, nrow(dt_sub))

      for (f in seq_len(K_INNER)) {
        tr_idx <- which(folds != f); va_idx <- which(folds == f)
        X_tr <- as.data.frame(dt_sub[tr_idx, ..ALL_PREDICTORS])
        y_tr <- as.integer(dt_sub[tr_idx][[TARGET_COL]])
        w_tr <- make_sample_weight(y_tr, hp$pos_multiplier)
        mdl  <- build_ebm(hp, EBM_N_JOBS, glassbox)
        mdl$fit(X_tr, y_tr, sample_weight = w_tr)
        X_va <- as.data.frame(dt_sub[va_idx, ..ALL_PREDICTORS])
        p_oof[va_idx] <- predict_proba_pos(mdl, X_va)
      }

      aucs    <- auc_pair(y_sub, p_oof)
      thr_obj <- tune_threshold_constraint(y_sub, p_oof, THR_GRID, SENS_TARGET)
      bm      <- bin_metrics(y_sub, p_oof, thr_obj$threshold)

      res_list[[i]] <- data.table(
        GRID[i], random_state = master_seed,
        oof_auroc = aucs["auroc"], oof_auprc = aucs["auprc"],
        oof_brier = brier_score(y_sub, p_oof), best_thr = thr_obj$threshold,
        met_sens80 = thr_obj$met_constraint, oof_f1_at_thr = bm$F1,
        oof_sens_at_thr = bm$Sensitivity, oof_spec_at_thr = bm$Specificity
      )

      log_msg(sprintf("  Grid %d/%d (leaf=%d, smooth=%d, pos=%.1f)",
                       i, nrow(GRID), hp$min_samples_leaf, hp$smoothing_rounds, hp$pos_multiplier))
    }

    res_dt       <- rbindlist(res_list, fill = TRUE)
    grid_elapsed <- as.numeric(difftime(Sys.time(), grid_start, units = "mins"))
    fwrite(res_dt, file.path(out_dir, sprintf("%s_tuning_results.csv", iter_ts)))

    # Select best: among those meeting sensitivity constraint, max AUPRC
    res_sel <- if (any(res_dt$met_sens80, na.rm = TRUE)) res_dt[met_sens80 == TRUE] else res_dt
    setorder(res_sel, -oof_auprc, -oof_auroc)
    best <- res_sel[1]

    log_msg(sprintf("  Grid complete: %.1f min", grid_elapsed))
    log_msg(sprintf("  Best: AUPRC=%.4f, leaf=%d, smooth=%d",
                     best$oof_auprc, best$min_samples_leaf, best$smoothing_rounds))

    # ---- Step 2.1: Refine on full E-only dev pool ----
    log_msg("\n[Step 2.1] Refine on full E-only dev pool")
    refine_start <- Sys.time()

    best_hp <- c(EBM_FIXED_HP, as.list(best[, .(min_samples_leaf, smoothing_rounds, pos_multiplier)]))
    best_hp$random_state <- master_seed
    best_hp$outer_bags   <- OUTER_BAGS_OOF

    fold_seed  <- master_seed + 1000L
    folds_full <- make_group_folds(dev_pool[[CLUSTER_COL]], k = K_INNER, seed = fold_seed)
    p_oof_full <- rep(NA_real_, nrow(dev_pool))

    for (f in seq_len(K_INNER)) {
      tr_idx <- which(folds_full != f); va_idx <- which(folds_full == f)
      X_tr <- as.data.frame(dev_pool[tr_idx, ..ALL_PREDICTORS])
      y_tr <- as.integer(dev_pool[tr_idx][[TARGET_COL]])
      w_tr <- make_sample_weight(y_tr, best_hp$pos_multiplier)
      mdl  <- build_ebm(best_hp, EBM_N_JOBS, glassbox)
      mdl$fit(X_tr, y_tr, sample_weight = w_tr)
      X_va <- as.data.frame(dev_pool[va_idx, ..ALL_PREDICTORS])
      p_oof_full[va_idx] <- predict_proba_pos(mdl, X_va)
    }

    aucs_full <- auc_pair(y_full, p_oof_full)
    thr_full  <- tune_threshold_constraint(y_full, p_oof_full, THR_GRID, SENS_TARGET)
    best_thr  <- thr_full$threshold

    refine_elapsed <- as.numeric(difftime(Sys.time(), refine_start, units = "mins"))
    log_msg(sprintf("  Final OOF: AUPRC=%.4f, AUROC=%.4f, Thr=%.4f",
                     aucs_full["auprc"], aucs_full["auroc"], best_thr))

    # Save OOF predictions
    oof_dt <- dev_pool[, .(Patient_ID = get(ID_COL), Patient_Cluster_ID = get(CLUSTER_COL),
                            site, y_true = as.integer(get(TARGET_COL)))]
    oof_dt[, y_prob := p_oof_full]
    fst::write_fst(oof_dt, file.path(out_dir, sprintf("%s_dev_oof_predictions.fst", iter_ts)))

    # ---- Step 3: Train Final Model (E-only dev, outer_bags=64) ----
    log_msg("\n[Step 3] Training Final Model (E-only dev pool)")

    X_dev <- as.data.frame(dev_pool[, ..ALL_PREDICTORS])
    y_dev <- as.integer(dev_pool[[TARGET_COL]])
    w_dev <- make_sample_weight(y_dev, best$pos_multiplier)

    hp_final <- c(EBM_FIXED_HP, as.list(best[, .(min_samples_leaf, smoothing_rounds, pos_multiplier)]))
    hp_final$random_state <- master_seed
    hp_final$outer_bags   <- OUTER_BAGS_FINAL

    ebm_final <- build_ebm(hp_final, EBM_N_JOBS, glassbox)
    ebm_final$fit(X_dev, y_dev, sample_weight = w_dev)

    model_path <- file.path(out_dir, sprintf("%s_Final_EBM.joblib", iter_ts))
    joblib$dump(ebm_final, model_path)
    log_msg(sprintf("  Model saved: %s", basename(model_path)))

    # Save artifacts
    saveRDS(list(
      iter = iter_id, external = ext_site, dev = dev_sites,
      seed_id = seed_id, master_seed = master_seed, timestamp = iter_ts,
      best_hp = best_hp, threshold = best_thr,
      internal_oof = list(AUROC = aucs_full["auroc"], AUPRC = aucs_full["auprc"],
                          Brier = brier_score(y_full, p_oof_full))
    ), file.path(out_dir, sprintf("%s_artifacts.rds", iter_ts)))

    # ---- Step 4: External Validation (E-only) ----
    log_msg(sprintf("\n[Step 4] External Validation — E-only (%s)", ext_site))

    X_ext <- as.data.frame(dt_ext[, ..ALL_PREDICTORS])
    y_ext <- as.integer(dt_ext[[TARGET_COL]])
    p_ext <- predict_proba_pos(ebm_final, X_ext)

    ext_pred <- dt_ext[, .(Patient_ID = get(ID_COL), Patient_Cluster_ID = get(CLUSTER_COL),
                            site, y_true = as.integer(get(TARGET_COL)))]
    ext_pred[, y_prob := p_ext]
    ext_pred[, y_hat := as.integer(y_prob >= best_thr)]
    fst::write_fst(ext_pred, file.path(out_dir, sprintf("%s_external_predictions.fst", iter_ts)))

    pt <- evaluate_point(ext_pred$y_true, ext_pred$y_prob, best_thr)
    fwrite(pt, file.path(out_dir, sprintf("%s_external_point_metrics.csv", iter_ts)))

    log_msg(sprintf("  AUROC=%.4f | AUPRC=%.4f | Sens=%.2f%%", pt$AUROC, pt$AUPRC, pt$Sensitivity * 100))

    # Calibration
    cal_dt <- calibration_bins(ext_pred$y_true, ext_pred$y_prob, CAL_BINS)
    fwrite(cal_dt, file.path(out_dir, sprintf("%s_calibration.csv", iter_ts)))

    # DCA
    dca_dt <- dca_curve(ext_pred$y_true, ext_pred$y_prob, DCA_PTS)
    fwrite(dca_dt, file.path(out_dir, sprintf("%s_dca.csv", iter_ts)))

    # Bootstrap CI
    log_msg(sprintf("\n[Step 4.1] Bootstrap %d", B_BOOT))
    boot <- bootstrap_external(ext_pred, best_thr, B_BOOT, BOOT_SEED, CLUSTER_COL, DCA_PTS)
    fwrite(boot$ci, file.path(out_dir, sprintf("%s_boot_ci.csv", iter_ts)))
    log_msg(sprintf("  Valid: %d/%d", boot$valid_count, B_BOOT))

    # DCA with CI
    if (!is.null(boot$dca_mat)) {
      dca_ci <- copy(dca_dt)
      dca_ci[, NB_lo := apply(boot$dca_mat, 2, quantile, 0.025, na.rm = TRUE)]
      dca_ci[, NB_hi := apply(boot$dca_mat, 2, quantile, 0.975, na.rm = TRUE)]
      fwrite(dca_ci, file.path(out_dir, sprintf("%s_dca_ci.csv", iter_ts)))
    }

    # Feature importance
    imp <- get_term_importances(ebm_final)
    if (!is.null(imp)) {
      term_names <- reticulate::py_to_r(ebm_final$term_names_)
      imp_dt <- data.table(feature = term_names, importance = as.numeric(imp))
      setorder(imp_dt, -importance)
      imp_dt[, rank := .I]
      fwrite(imp_dt, file.path(out_dir, sprintf("%s_feature_importance.csv", iter_ts)))
    }

    iter_elapsed <- as.numeric(difftime(Sys.time(), iter_start, units = "mins"))
    log_msg(sprintf("\n[Complete] %.2f min", iter_elapsed))
    close(log_con)

    # Return summary row
    return(data.table(
      iter = iter_id, external = ext_site, dev = paste(dev_sites, collapse = "+"),
      seed_id = seed_id, master_seed = master_seed,
      hp_min_samples_leaf = best$min_samples_leaf,
      hp_smoothing_rounds = best$smoothing_rounds,
      hp_pos_multiplier = best$pos_multiplier,
      internal_oof_auroc = aucs_full["auroc"], internal_oof_auprc = aucs_full["auprc"],
      internal_oof_brier = brier_score(y_full, p_oof_full), threshold = best_thr,
      external_AUROC = pt$AUROC, external_AUPRC = pt$AUPRC, external_Brier = pt$Brier,
      external_Sens = pt$Sensitivity, external_Spec = pt$Specificity, external_F1 = pt$F1,
      boot_valid_n = boot$valid_count, elapsed_min = iter_elapsed,
      grid_elapsed_min = grid_elapsed, refine_elapsed_min = refine_elapsed,
      timestamp = iter_ts
    ))

  }, error = function(e) {
    log_msg(sprintf("\n[ERROR] %s", e$message))
    close(log_con)
    return(data.table(iter = iter_id, external = ext_site, seed_id = seed_id, error = e$message))
  })
}

# =============================================================================
# Phase 8: Parallel Execution (6 workers)
# =============================================================================

cat("\n[START] Parallel execution of 6 iterations (3 sites × 2 seeds)\n")
cat(sprintf("[TIME] %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

start_time <- Sys.time()

cl <- makeCluster(N_WORKERS, type = "PSOCK")

clusterExport(cl, c(
  # Config
  "OUT_ROOT", "VENV_PATH", "SITE_FILES", "GRID", "EBM_FIXED_HP",
  "ID_COL", "CLUSTER_COL", "TARGET_COL", "BINARY_COLS", "CONT_COLS", "ALL_PREDICTORS",
  "COMPAT_RENAME_MAP", "PHYS_RANGES",
  "DO_PHYS_CHECK", "DO_SCALE", "K_INNER", "CV_FOLD_SEED", "SENS_TARGET", "THR_GRID",
  "DCA_PTS", "CAL_BINS", "B_BOOT", "BOOT_SEED",
  "OUTER_BAGS_OOF", "OUTER_BAGS_FINAL", "GRID_SUBSAMPLE_RATIO", "EBM_N_JOBS",
  "ENFORCE_SENS_CONSTRAINT", "N_THREADS_PER_WORKER",
  # E/L split config
  "SESSION_DATE_COL", "E_SPLIT_QUANTILE", "COHORT_COL",
  # Functions from 00_utils_r.R
  "read_fst_dt", "apply_compat_rename", "to_target01", "to_binary01",
  "apply_physio_ranges", "check_required_columns",
  "make_group_folds", "prep_site_dt_worker", "apply_el_split",
  "make_sample_weight", "get_positive_index", "build_ebm", "predict_proba_pos",
  "auc_pair", "brier_score", "bin_metrics", "confusion_by_thresholds", "tune_threshold_constraint",
  "calibration_bins", "dca_curve", "evaluate_point", "bootstrap_external",
  "get_term_importances"
))

results <- parLapply(cl, ITERATIONS, run_iteration_worker)
stopCluster(cl)

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

# =============================================================================
# Phase 9: Merge Results → Table 2 & Supp S2
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Merging Results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

summ_dt <- rbindlist(results, fill = TRUE)

# Handle errors
if ("error" %in% names(summ_dt)) {
  errors <- summ_dt[!is.na(error)]
  if (nrow(errors) > 0) {
    cat("[WARNING] Some iterations failed:\n")
    print(errors[, .(iter, external, seed_id, error)])
  }
  summ_dt <- summ_dt[is.na(error)]
  if ("error" %in% names(summ_dt)) summ_dt[, error := NULL]
}

fwrite(summ_dt, file.path(OUT_ROOT, "IECV_iteration_summary.csv"))

# --- Table 2 metrics ---
if (nrow(summ_dt) > 0 && "external_AUROC" %in% names(summ_dt)) {
  stab <- summ_dt[, .(
    n = .N,
    ext_AUROC_mean = mean(external_AUROC, na.rm = TRUE),
    ext_AUROC_sd   = sd(external_AUROC, na.rm = TRUE),
    ext_AUPRC_mean = mean(external_AUPRC, na.rm = TRUE),
    ext_AUPRC_sd   = sd(external_AUPRC, na.rm = TRUE),
    ext_Brier_mean = mean(external_Brier, na.rm = TRUE),
    ext_Sens_mean  = mean(external_Sens, na.rm = TRUE),
    ext_Spec_mean  = mean(external_Spec, na.rm = TRUE)
  )]
  fwrite(stab, file.path(OUT_ROOT, "GLOBAL_stability.csv"))

  # Transportability gap
  gap <- summ_dt[, .(iter, external, seed_id,
                      gap_AUROC = internal_oof_auroc - external_AUROC,
                      gap_AUPRC = internal_oof_auprc - external_AUPRC)]
  fwrite(gap, file.path(OUT_ROOT, "GLOBAL_transportability_gap.csv"))

  # Seed stability
  seed_stab <- summ_dt[, .(
    AUROC_range = diff(range(external_AUROC, na.rm = TRUE)),
    AUPRC_range = diff(range(external_AUPRC, na.rm = TRUE))
  ), by = external]
  fwrite(seed_stab, file.path(OUT_ROOT, "GLOBAL_seed_stability.csv"))
}

# --- Shape function export (for 05_shapeqc.R) ---
cat("  Exporting shape functions for downstream ShapeQC ...\n")
iter_tags <- vapply(ITERATIONS, function(it)
  sprintf("Iter%d_External_%s_Seed%d", it$iter, it$external, it$seed_id), character(1))
tryCatch({
  if (!exists("init_python_full")) src("R/00_utils_python.R")
  init_python_full()
  joblib_main <- reticulate::import("joblib", convert = TRUE)

  for (tag in iter_tags) {
    tag_dir <- file.path(OUT_ROOT, tag)
    model_files <- list.files(tag_dir, pattern = "_Final_EBM\\.joblib$", full.names = TRUE)
    if (length(model_files) == 0) next

    ebm <- tryCatch(joblib_main$load(model_files[length(model_files)]),
                     error = function(e) NULL)
    if (is.null(ebm)) next

    shapes_dir <- file.path(tag_dir, "Shapes")
    dir.create(shapes_dir, recursive = TRUE, showWarnings = FALSE)

    shape_rows <- list()
    for (feat in CONT_COLS) {
      bins_info <- tryCatch(
        extract_bins_with_diagnostics(ebm, feat, strict_continuous_cuts = FALSE, verbose = FALSE),
        error = function(e) list(success = FALSE))
      if (!bins_info$success) next
      dt_b <- bins_info$dt
      shape_rows[[length(shape_rows) + 1L]] <- data.table::data.table(
        feature = feat, bin_index = dt_b$bin_index,
        x = dt_b$x, score = dt_b$score, bin_weight = dt_b$bin_weight,
        is_missing = dt_b$is_missing, is_unknown = dt_b$is_unknown
      )
    }
    if (length(shape_rows) > 0) {
      shape_dt <- data.table::rbindlist(shape_rows, fill = TRUE)
      data.table::fwrite(shape_dt, file.path(shapes_dir, "feature_data.csv"))
    }
    rm(ebm); gc(verbose = FALSE)
  }
  cat("  Shape export complete.\n")
}, error = function(e) {
  cat(sprintf("  [WARN] Shape export failed: %s\n", e$message))
  cat("  Downstream 05_shapeqc.R may need manual shape extraction.\n")
})

# --- Supp S2: Feature importance aggregation ---
# (iter_tags already computed above for shape export)

imp_all <- rbindlist(lapply(iter_tags, function(tag) {
  fs <- list.files(file.path(OUT_ROOT, tag), pattern = "feature_importance\\.csv$", full.names = TRUE)
  if (length(fs) == 0) return(data.table())
  dt <- fread(fs[length(fs)])
  dt[, model := tag]
  dt
}), fill = TRUE)

if (nrow(imp_all) > 0) {
  fwrite(imp_all, file.path(OUT_ROOT, "GLOBAL_feature_importance.csv"))
  cat("  Supp S2: Feature importance saved.\n")
}

# Save complete results
saveRDS(list(
  summary = summ_dt, total_time_min = total_time, run_timestamp = RUN_TS
), file.path(OUT_ROOT, "IECV_Complete_Results.rds"))

# =============================================================================
# Phase 10: Final Report
# =============================================================================

cat("\n")
cat(paste(rep("#", 70), collapse = ""), "\n")
cat("# IECV-EBM Complete\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n")

cat("===== Table 2: E-only External Validation Summary =====\n")
if (nrow(summ_dt) > 0 && "external_AUROC" %in% names(summ_dt)) {
  print(summ_dt[, .(iter, external, seed_id, hp_min_samples_leaf, hp_smoothing_rounds,
                     external_AUROC = round(external_AUROC, 4),
                     external_AUPRC = round(external_AUPRC, 4),
                     external_Sens = round(external_Sens, 4))])

  cat(sprintf("\nMean ± SD: AUROC %.4f ± %.4f | AUPRC %.4f ± %.4f\n",
              mean(summ_dt$external_AUROC), sd(summ_dt$external_AUROC),
              mean(summ_dt$external_AUPRC), sd(summ_dt$external_AUPRC)))

  cat("\n[Verify] Previous full-cohort values were:\n")
  cat("  Mean AUROC 0.895 ± 0.018 | Mean AUPRC 0.583 ± 0.082\n")
  cat("  E-only values may differ slightly due to reduced sample size.\n")
  cat("  If hyperparameters changed → compare below with full-cohort selections.\n")
}

cat(sprintf("\nWall-clock: %.1f min\n", total_time))
cat(sprintf("Output: %s\n", OUT_ROOT))

# --- Hyperparameter comparison: E-only vs full-cohort (from Supp Methods) ---
cat("\n===== Hyperparameter Comparison: E-only vs Full-cohort =====\n")
cat("Full-cohort reference (from prior run):\n")
cat("  TN: min_samples_leaf=300,  smoothing_rounds=600\n")
cat("  D6: min_samples_leaf=1000, smoothing_rounds=600\n")
cat("  CY: min_samples_leaf=600,  smoothing_rounds=75\n")
if (nrow(summ_dt) > 0 && "hp_min_samples_leaf" %in% names(summ_dt)) {
  cat("E-only results (this run):\n")
  for (site in c("TN", "D6", "CY")) {
    row <- summ_dt[external == site & seed_id == 1]
    if (nrow(row) > 0) {
      cat(sprintf("  %s: min_samples_leaf=%d, smoothing_rounds=%d\n",
                  site, row$hp_min_samples_leaf, row$hp_smoothing_rounds))
    }
  }
  cat("\nIf hyperparameters match → downstream scripts need text-only changes.\n")
  cat("If hyperparameters differ → full pipeline re-run required (scripts 05–12).\n")
}

cat("\n")
