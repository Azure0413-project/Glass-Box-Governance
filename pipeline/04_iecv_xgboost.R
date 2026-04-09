# =============================================================================
# 04_iecv_xgboost.R — E-only XGBoost Matched Benchmark (IECV)
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Supp Table S1 — Head-to-head E-only external validation: EBM vs. XGBoost
#
# Source: consolidated from original analysis scripts
# Design: uses identical E-only IECV splits / predictors / preprocessing as 03_iecv_ebm.R v3.0
#       Only the model is replaced with XGBoost. Serves as a discrimination-only benchmark.
#
# Design change (v3.3.2): Pure E-only pipeline
#   All data loading, tuning, training, and validation use exclusively
#   E-cohort (earlier 75%) sessions, matching 03_iecv_ebm.R v3.0.
#
# Methods → "A matched XGBoost classifier served as a secondary benchmark
#  for discrimination only. It used the same predictors, preprocessing rules,
#  E-only IECV splits, and patient-level separation as the EBM"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(xgboost)
  library(precrec)
  library(ggplot2)
  library(gridExtra)
  library(parallel)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")

cat("\n")
cat(paste(rep("#", 70), collapse = ""), "\n")
cat("# Module 04: E-only IECV-XGBoost v3.3.2 (Matched Benchmark)\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n")

# --- Output directory (parallel to EBM output) ---
XGB_BASE_DIR <- gsub("IECV_EBM", "IECV_XGBoost", BASE_DIR)
OUT_ROOT_XGB <- file.path(XGB_BASE_DIR, paste0("Run_", RUN_TS))
dir.create(OUT_ROOT_XGB, recursive = TRUE, showWarnings = FALSE)

# --- XGBoost Grid Search (optimized: 16 combinations) ---
XGB_GRID <- expand.grid(
  max_depth      = c(3L, 6L),
  min_child_weight = c(50L, 200L),
  subsample      = c(0.7),
  colsample_bytree = c(0.7),
  eta            = c(0.01, 0.05),
  gamma          = c(0, 1),
  stringsAsFactors = FALSE
)
XGB_NROUNDS_MAX <- 2000L
XGB_EARLY_STOP  <- 50L

cat(sprintf("  OUT_ROOT: %s\n", OUT_ROOT_XGB))
cat(sprintf("  Grid: %d combinations (optimized from 162)\n", nrow(XGB_GRID)))

# =============================================================================
# XGBoost Core Functions
# =============================================================================

#' v3.3.1 fix: serialize/deserialize to avoid ALTREP issues in parallel
fix_xgb_altrep <- function(mdl) {
  tmp <- tempfile(fileext = ".xgb")
  on.exit(unlink(tmp), add = TRUE)
  xgboost::xgb.save(mdl, tmp)
  xgboost::xgb.load(tmp)
}

#' Train XGBoost with early stopping
train_xgb <- function(X_train, y_train, X_val, y_val, params,
                       nrounds = XGB_NROUNDS_MAX, early_stop = XGB_EARLY_STOP,
                       scale_pos_weight = 1.0) {
  dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  dval   <- xgb.DMatrix(data = as.matrix(X_val),   label = y_val)

  full_params <- c(params, list(
    objective        = "binary:logistic",
    eval_metric      = "aucpr",
    scale_pos_weight = scale_pos_weight,
    nthread          = N_THREADS_PER_WORKER
  ))

  mdl <- xgb.train(
    params    = full_params,
    data      = dtrain,
    nrounds   = nrounds,
    watchlist  = list(val = dval),
    early_stopping_rounds = early_stop,
    verbose   = 0
  )

  fix_xgb_altrep(mdl)
}

#' Predict probability from XGBoost
predict_proba_xgb <- function(mdl, X) {
  dmat <- xgb.DMatrix(data = as.matrix(X))
  as.numeric(predict(mdl, dmat))
}

# =============================================================================
# Worker Function
# =============================================================================

run_xgb_iteration <- function(it) {

  suppressPackageStartupMessages({
    library(data.table); library(fst); library(xgboost)
    library(precrec); library(ggplot2)
  })

  iter_id     <- it$iter
  ext_site    <- it$external
  dev_sites   <- it$dev
  seed_id     <- it$seed_id
  master_seed <- it$master_seed

  iter_tag <- sprintf("Iter%d_External_%s_Seed%d", iter_id, ext_site, seed_id)
  out_dir  <- file.path(OUT_ROOT_XGB, iter_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  iter_start <- Sys.time()

  tryCatch({

    # Step 1: Data (identical to EBM pipeline — E-only)
    dt_dev1 <- prep_site_dt(read_fst_dt(SITE_FILES[[dev_sites[1]]]), dev_sites[1])
    dt_dev2 <- prep_site_dt(read_fst_dt(SITE_FILES[[dev_sites[2]]]), dev_sites[2])
    dt_ext  <- prep_site_dt(read_fst_dt(SITE_FILES[[ext_site]]),     ext_site)

    # Apply E/L split and filter to E-only (matching 03_iecv_ebm.R v3.0)
    dt_dev1 <- apply_el_split(dt_dev1, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
    dt_dev2 <- apply_el_split(dt_dev2, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
    dt_ext  <- apply_el_split(dt_ext,  SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
    dt_dev1 <- dt_dev1[get(COHORT_COL) == "E"]
    dt_dev2 <- dt_dev2[get(COHORT_COL) == "E"]
    dt_ext  <- dt_ext[get(COHORT_COL)  == "E"]

    dev_pool <- rbindlist(list(dt_dev1, dt_dev2), use.names = TRUE, fill = TRUE)

    y_full <- as.integer(dev_pool[[TARGET_COL]])
    n1 <- sum(y_full == 1L); n0 <- sum(y_full == 0L)
    scale_pos_weight <- n0 / n1

    # Step 2: Grid Search (subsampled, 5-fold)
    set.seed(master_seed)
    cluster_ids      <- unique(dev_pool[[CLUSTER_COL]])
    n_sample         <- ceiling(length(cluster_ids) * GRID_SUBSAMPLE_RATIO)
    sampled_clusters <- sample(cluster_ids, n_sample)
    dt_sub           <- dev_pool[get(CLUSTER_COL) %in% sampled_clusters]
    y_sub            <- as.integer(dt_sub[[TARGET_COL]])

    fold_seed <- master_seed + 1000L
    folds     <- make_group_folds(dt_sub[[CLUSTER_COL]], k = K_INNER, seed = fold_seed)

    best_auprc <- -Inf; best_params <- NULL

    for (gi in seq_len(nrow(XGB_GRID))) {
      params <- as.list(XGB_GRID[gi])
      p_oof  <- rep(NA_real_, nrow(dt_sub))

      for (f in seq_len(K_INNER)) {
        tr_idx <- which(folds != f); va_idx <- which(folds == f)
        X_tr <- as.data.frame(dt_sub[tr_idx, ..ALL_PREDICTORS])
        X_va <- as.data.frame(dt_sub[va_idx, ..ALL_PREDICTORS])
        y_tr <- y_sub[tr_idx]; y_va <- y_sub[va_idx]

        mdl <- train_xgb(X_tr, y_tr, X_va, y_va, params,
                          scale_pos_weight = scale_pos_weight)
        p_oof[va_idx] <- predict_proba_xgb(mdl, X_va)
      }

      aucs <- auc_pair(y_sub, p_oof)
      if (!is.na(aucs["auprc"]) && aucs["auprc"] > best_auprc) {
        best_auprc  <- aucs["auprc"]
        best_params <- params
      }
    }

    # Step 2.1: Refine on full dev
    folds_full <- make_group_folds(dev_pool[[CLUSTER_COL]], k = K_INNER, seed = fold_seed)
    p_oof_full <- rep(NA_real_, nrow(dev_pool))

    for (f in seq_len(K_INNER)) {
      tr_idx <- which(folds_full != f); va_idx <- which(folds_full == f)
      X_tr <- as.data.frame(dev_pool[tr_idx, ..ALL_PREDICTORS])
      X_va <- as.data.frame(dev_pool[va_idx, ..ALL_PREDICTORS])
      y_tr <- y_full[tr_idx]; y_va <- y_full[va_idx]

      mdl <- train_xgb(X_tr, y_tr, X_va, y_va, best_params,
                        scale_pos_weight = scale_pos_weight)
      p_oof_full[va_idx] <- predict_proba_xgb(mdl, X_va)
    }

    aucs_full <- auc_pair(y_full, p_oof_full)
    thr_obj   <- tune_threshold_constraint(y_full, p_oof_full, THR_GRID, SENS_TARGET)
    best_thr  <- thr_obj$threshold

    # Step 3: Train final model
    X_dev <- as.data.frame(dev_pool[, ..ALL_PREDICTORS])
    X_ext <- as.data.frame(dt_ext[, ..ALL_PREDICTORS])
    y_ext <- as.integer(dt_ext[[TARGET_COL]])

    # Use last fold as pseudo-validation for early stopping
    va_idx_last <- which(folds_full == K_INNER)
    xgb_final <- train_xgb(X_dev, y_full,
                             as.data.frame(dev_pool[va_idx_last, ..ALL_PREDICTORS]),
                             y_full[va_idx_last],
                             best_params, scale_pos_weight = scale_pos_weight)
    xgb.save(xgb_final, file.path(out_dir, "Final_XGBoost.xgb"))

    # Step 4: External validation
    p_ext <- predict_proba_xgb(xgb_final, X_ext)
    pt    <- evaluate_point(y_ext, p_ext, best_thr)
    fwrite(pt, file.path(out_dir, "external_point_metrics.csv"))

    # Bootstrap
    ext_pred <- data.table(
      Patient_Cluster_ID = dt_ext[[CLUSTER_COL]],
      y_true = y_ext, y_prob = p_ext
    )
    boot <- bootstrap_external(ext_pred, best_thr, B_BOOTSTRAP, BOOTSTRAP_SEED,
                                "Patient_Cluster_ID", DCA_PTS)
    fwrite(boot$ci, file.path(out_dir, "boot_ci.csv"))

    iter_elapsed <- as.numeric(difftime(Sys.time(), iter_start, units = "mins"))

    data.table(
      iter = iter_id, external = ext_site, seed_id = seed_id, master_seed = master_seed,
      hp_max_depth = best_params$max_depth, hp_eta = best_params$eta,
      internal_oof_auroc = aucs_full["auroc"], internal_oof_auprc = aucs_full["auprc"],
      threshold = best_thr,
      external_AUROC = pt$AUROC, external_AUPRC = pt$AUPRC, external_Brier = pt$Brier,
      external_Sens = pt$Sensitivity, external_Spec = pt$Specificity,
      external_PPV = pt$PPV, external_NPV = pt$NPV, external_F1 = pt$F1,
      elapsed_min = iter_elapsed
    )

  }, error = function(e) {
    data.table(iter = iter_id, external = ext_site, seed_id = seed_id, error = e$message)
  })
}

# =============================================================================
# Parallel Execution
# =============================================================================

cat("\n[START] Parallel XGBoost (6 iterations)\n")
start_time <- Sys.time()

cl <- makeCluster(N_WORKERS, type = "PSOCK")
clusterExport(cl, c(
  "OUT_ROOT_XGB", "SITE_FILES", "XGB_GRID", "XGB_NROUNDS_MAX", "XGB_EARLY_STOP",
  "ID_COL", "CLUSTER_COL", "TARGET_COL", "BINARY_COLS", "CONT_COLS", "ALL_PREDICTORS",
  "COMPAT_RENAME_MAP", "PHYS_RANGES",
  "SESSION_DATE_COL", "E_SPLIT_QUANTILE", "COHORT_COL",
  "K_INNER", "SENS_TARGET", "THR_GRID", "GRID_SUBSAMPLE_RATIO",
  "N_THREADS_PER_WORKER", "B_BOOTSTRAP", "BOOTSTRAP_SEED", "DCA_PTS", "CAL_BINS",
  "read_fst_dt", "apply_compat_rename", "to_target01", "to_binary01",
  "apply_physio_ranges", "check_required_columns", "prep_site_dt",
  "make_group_folds", "apply_el_split", "auc_pair", "brier_score", "bin_metrics",
  "confusion_by_thresholds", "tune_threshold_constraint",
  "calibration_bins", "dca_curve", "evaluate_point", "bootstrap_external",
  "fix_xgb_altrep", "train_xgb", "predict_proba_xgb"
))

results <- parLapply(cl, ITERATIONS, run_xgb_iteration)
stopCluster(cl)

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

# =============================================================================
# Results: Supp Table S1
# =============================================================================

summ_xgb <- rbindlist(results, fill = TRUE)
if ("error" %in% names(summ_xgb) && any(!is.na(summ_xgb$error))) {
  cat("[WARNING] Errors:\n")
  print(summ_xgb[!is.na(error), .(iter, external, error)])
  summ_xgb <- summ_xgb[is.na(error)]
}

fwrite(summ_xgb, file.path(OUT_ROOT_XGB, "IECV_XGB_iteration_summary.csv"))

cat("\n===== Supp Table S1: XGBoost Benchmark =====\n")
if (nrow(summ_xgb) > 0) {
  print(summ_xgb[, .(iter, external, seed_id,
                       external_AUROC = round(external_AUROC, 4),
                       external_AUPRC = round(external_AUPRC, 4),
                       external_Sens = round(external_Sens, 4))])

  cat(sprintf("\nXGB Mean: AUROC %.4f ± %.4f | AUPRC %.4f ± %.4f\n",
              mean(summ_xgb$external_AUROC, na.rm = TRUE), sd(summ_xgb$external_AUROC, na.rm = TRUE),
              mean(summ_xgb$external_AUPRC, na.rm = TRUE), sd(summ_xgb$external_AUPRC, na.rm = TRUE)))

  cat("\n[Verify] Previous full-cohort values were:\n")
  cat("  XGB: AUROC 0.898 ± 0.018, AUPRC 0.596 ± 0.086\n")
  cat("  Difference (XGB-EBM): AUROC +0.003, AUPRC +0.012\n")
  cat("  E-only values may differ slightly due to reduced sample size.\n")
}

cat(sprintf("\nWall-clock: %.1f min | Output: %s\n\n", total_time, OUT_ROOT_XGB))
