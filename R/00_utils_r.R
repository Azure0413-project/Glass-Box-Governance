# =============================================================================
# 00_utils_r.R — Shared R Utility Functions
# IDH EBM Governance Reproducibility Pipeline
#
# Shared R utility functions used by all pipeline modules.
# Function logic is preserved from the original IECV-EBM v2.5.10 / Post-hoc v2.0-v2.2 / NRI v1.1.
# Duplicate functions were deduplicated; no computational logic was changed.
#
# Usage: src("R/00_config.R"); source("00_utils_r.R")
# =============================================================================

cat("  Loading 00_utils_r.R ... ")

# =============================================================================
# Section 1: Data I/O & Preprocessing
# =============================================================================

#' Read a .fst file and return as data.table
read_fst_dt <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  fst::read_fst(path, as.data.table = TRUE)
}

#' Apply cross-site column name harmonization
apply_compat_rename <- function(dt, rename_map = COMPAT_RENAME_MAP) {
  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(dt) && !(new %in% names(dt))) {
      data.table::setnames(dt, old, new)
    }
  }
  invisible(dt)
}

#' Convert outcome column to 0/1 integer (handles various text formats)
to_target01 <- function(x) {
  if (is.null(x)) return(NA_integer_)
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(ifelse(is.na(x), NA, x)))
  x2 <- toupper(trimws(as.character(x)))
  x2[x2 %in% c("TRUE", "T", "YES", "Y", "POS", "POSITIVE", "1")] <- "1"
  x2[x2 %in% c("FALSE", "F", "NO", "N", "NEG", "NEGATIVE", "0")] <- "0"
  suppressWarnings(as.integer(x2))
}

#' Convert binary predictor (Sex) to 0/1 integer
to_binary01 <- function(x) {
  if (is.null(x)) return(NA_integer_)
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(ifelse(is.na(x), NA, x)))
  x2 <- toupper(trimws(as.character(x)))
  x2[x2 %in% c("M", "MALE", "MAN", "BOY", "1", "TRUE", "T", "\u7537")] <- "1"
  x2[x2 %in% c("F", "FEMALE", "WOMAN", "GIRL", "0", "FALSE", "\u5973")] <- "0"
  suppressWarnings(as.integer(x2))
}

#' Physiological range filter: replace out-of-range values with NA
apply_physio_ranges <- function(dt, range_map = PHYS_RANGES) {
  for (nm in names(range_map)) {
    if (nm %in% names(dt)) {
      lo <- range_map[[nm]][1]
      hi <- range_map[[nm]][2]
      dt[!is.na(get(nm)) & (get(nm) < lo | get(nm) > hi), (nm) := NA_real_]
    }
  }
  invisible(dt)
}

#' Validate that required columns exist in a data.table
check_required_columns <- function(dt, required_cols, site_name = "Unknown") {
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("[%s] Missing columns: %s", site_name, paste(missing_cols, collapse = ", ")))
  }
  invisible(TRUE)
}

#' Complete single-site preprocessing pipeline
#' (rename -> target encoding -> cluster ID -> binary encoding -> physio filter)
prep_site_dt <- function(dt_raw, site_name) {
  dt <- data.table::copy(dt_raw)
  dt <- apply_compat_rename(dt, COMPAT_RENAME_MAP)
  dt[, site := site_name]

  # Target
  if (!(TARGET_COL %in% names(dt)))
    stop(sprintf("[%s] TARGET_COL '%s' not found", site_name, TARGET_COL))
  dt[, (TARGET_COL) := to_target01(get(TARGET_COL))]
  dt <- dt[!is.na(get(TARGET_COL))]

  # Cluster ID (patient-level, site-prefixed)
  if (!(ID_COL %in% names(dt)))
    stop(sprintf("[%s] ID_COL '%s' not found", site_name, ID_COL))
  dt[, (CLUSTER_COL) := paste0(site_name, "__", as.character(get(ID_COL)))]

  # Binary encoding
  for (b in BINARY_COLS) {
    if (b %in% names(dt)) dt[, (b) := to_binary01(get(b))]
  }

  # Continuous type coercion
  for (cn in CONT_COLS) {
    if (cn %in% names(dt)) dt[, (cn) := as.numeric(get(cn))]
  }

  # Physiological range filtering
  dt <- apply_physio_ranges(dt, PHYS_RANGES)

  dt
}

#' Extract predictor data.frame (for EBM/XGBoost input)
get_predictor_df <- function(dt) {
  X <- data.table::copy(dt)
  for (b in BINARY_COLS) if (!(b %in% names(X))) X[, (b) := NA_integer_]
  for (cc in CONT_COLS)  if (!(cc %in% names(X))) X[, (cc) := NA_real_]
  as.data.frame(X[, ..ALL_PREDICTORS])
}

# =============================================================================
# Section 2: Cross-Validation Utilities
# =============================================================================

#' Patient-level grouped k-fold splitting
make_group_folds <- function(group_ids, k = 5L, seed = 2024L) {
  set.seed(seed)
  uid <- unique(as.character(group_ids))
  uid <- sample(uid, length(uid))
  fold_id <- rep(seq_len(k), length.out = length(uid))
  map <- setNames(fold_id, uid)
  as.integer(map[as.character(group_ids)])
}

#' Positive-class reweighted sample weights
#' pos_multiplier=1.0 (or NULL) → uniform weights (no reweighting).
#' For inverse-prevalence balanced weighting, pass pos_multiplier = n0/n1.
make_sample_weight <- function(y, pos_multiplier = NULL) {
  if (any(!y %in% c(0L, 1L))) stop("y must be 0/1")
  n  <- length(y)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) stop("Single-class in training set.")

  if (!is.null(pos_multiplier) && pos_multiplier != 1.0) {
    base_w <- c(1, as.numeric(pos_multiplier))
    scale_factor <- n / (n0 * base_w[1] + n1 * base_w[2])
    w0 <- scale_factor * base_w[1]
    w1 <- scale_factor * base_w[2]
    return(ifelse(y == 1L, w1, w0))
  }
  rep(1, n)
}

# =============================================================================
# Section 3: Discrimination & Classification Metrics
# =============================================================================

#' AUROC + AUPRC via precrec
auc_pair <- function(y_true, y_prob) {
  valid_idx <- is.finite(y_prob) & !is.na(y_true)
  if (sum(valid_idx) < 10 || length(unique(y_true[valid_idx])) < 2)
    return(c(auroc = NA_real_, auprc = NA_real_))
  tryCatch({
    mm <- precrec::evalmod(scores = y_prob[valid_idx], labels = y_true[valid_idx])
    aa <- precrec::auc(mm)
    auroc <- as.numeric(subset(aa, curvetypes == "ROC")$aucs[1])
    auprc <- as.numeric(subset(aa, curvetypes == "PRC")$aucs[1])
    c(auroc = auroc, auprc = auprc)
  }, error = function(e) c(auroc = NA_real_, auprc = NA_real_))
}

#' Brier score
brier_score <- function(y_true, y_prob) {
  valid_idx <- is.finite(y_prob) & !is.na(y_true)
  if (sum(valid_idx) == 0) return(NA_real_)
  mean((y_prob[valid_idx] - y_true[valid_idx])^2)
}

#' Binary metrics at a single threshold
bin_metrics <- function(y_true, y_prob, thr) {
  if (!is.finite(thr))
    return(list(TP=NA,FP=NA,TN=NA,FN=NA,Sensitivity=NA,Specificity=NA,PPV=NA,NPV=NA,F1=NA))
  valid_idx <- is.finite(y_prob) & !is.na(y_true)
  y_true <- y_true[valid_idx]; y_prob <- y_prob[valid_idx]
  y_hat <- ifelse(y_prob >= thr, 1L, 0L)
  TP <- sum(y_true == 1L & y_hat == 1L)
  FP <- sum(y_true == 0L & y_hat == 1L)
  TN <- sum(y_true == 0L & y_hat == 0L)
  FN <- sum(y_true == 1L & y_hat == 0L)
  Sens <- if ((TP+FN)>0) TP/(TP+FN) else NA_real_
  Spec <- if ((TN+FP)>0) TN/(TN+FP) else NA_real_
  PPV  <- if ((TP+FP)>0) TP/(TP+FP) else NA_real_
  NPV  <- if ((TN+FN)>0) TN/(TN+FN) else NA_real_
  F1   <- if (!is.na(PPV) && !is.na(Sens) && (PPV+Sens)>0)
             2*PPV*Sens/(PPV+Sens) else NA_real_
  list(TP=TP, FP=FP, TN=TN, FN=FN,
       Sensitivity=Sens, Specificity=Spec, PPV=PPV, NPV=NPV, F1=F1)
}

#' Confusion matrix across a grid of thresholds (vectorised)
confusion_by_thresholds <- function(y_true, y_prob, thr_grid) {
  ord   <- order(y_prob, decreasing = TRUE)
  p_dec <- y_prob[ord]
  y_dec <- y_true[ord]
  n1 <- sum(y_true == 1L)
  n0 <- sum(y_true == 0L)

  cumTP <- c(0L, cumsum(y_dec == 1L))
  cumFP <- c(0L, cumsum(y_dec == 0L))
  q     <- -p_dec

  m  <- findInterval(-thr_grid, q)
  TP <- cumTP[m + 1L]
  FP <- cumFP[m + 1L]
  FN <- n1 - TP
  TN <- n0 - FP

  Sens <- ifelse((TP+FN)>0, TP/(TP+FN), NA_real_)
  Spec <- ifelse((TN+FP)>0, TN/(TN+FP), NA_real_)
  PPV  <- ifelse((TP+FP)>0, TP/(TP+FP), NA_real_)
  F1   <- ifelse(!is.na(PPV) & !is.na(Sens) & (PPV+Sens)>0,
                 2*PPV*Sens/(PPV+Sens), NA_real_)

  data.table::data.table(
    thr = thr_grid, TP = TP, FP = FP, TN = TN, FN = FN,
    sensitivity = Sens, specificity = Spec, precision = PPV, f1 = F1
  )
}

#' Threshold tuning: sensitivity ≥ target, then max F1
tune_threshold_constraint <- function(y_true, y_prob, thr_grid, sens_target = 0.80) {
  dt <- confusion_by_thresholds(y_true, y_prob, thr_grid)
  cand <- dt[!is.na(sensitivity) & sensitivity >= sens_target & !is.na(f1)]

  if (nrow(cand) == 0) {
    best <- dt[which.max(sensitivity)]
    return(list(threshold = best$thr, met_constraint = FALSE, table = dt))
  }

  best <- cand[order(-f1, -sensitivity, -precision)][1]
  list(threshold = best$thr, met_constraint = TRUE, table = dt)
}

# =============================================================================
# Section 4: Calibration & Decision Curve Analysis
# =============================================================================

#' Quantile-based calibration bins
calibration_bins <- function(y_true, y_prob, n_bins = 10L) {
  dt <- data.table::data.table(y = y_true, p = y_prob)
  qs <- unique(stats::quantile(dt$p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(qs) < 3) {
    dt[, bin := "all"]
  } else {
    dt[, bin := cut(p, breaks = qs, include.lowest = TRUE)]
  }
  dt[, .(n = .N, mean_pred = mean(p), obs_rate = mean(y)), by = bin][order(mean_pred)]
}

#' Decision curve analysis (net benefit)
dca_curve <- function(y_true, y_prob, pts = seq(0.01, 0.99, by = 0.01)) {
  valid <- complete.cases(y_true, y_prob)
  y_true <- y_true[valid]
  y_prob <- y_prob[valid]
  n <- length(y_true)

  if (n == 0) {
    return(data.table::data.table(
      threshold = pts, NB_model = NA_real_, NB_all = NA_real_, NB_none = 0
    ))
  }

  prev    <- mean(y_true)
  nb_all  <- prev - (1 - prev) * (pts / (1 - pts))
  nb_none <- rep(0, length(pts))

  nb_model <- sapply(pts, function(pt) {
    y_hat <- ifelse(y_prob >= pt, 1L, 0L)
    TP <- sum(y_true == 1L & y_hat == 1L, na.rm = TRUE)
    FP <- sum(y_true == 0L & y_hat == 1L, na.rm = TRUE)
    TP / n - FP / n * (pt / (1 - pt))
  })

  data.table::data.table(threshold = pts, NB_model = nb_model, NB_all = nb_all, NB_none = nb_none)
}

#' Point estimate of all external validation metrics
evaluate_point <- function(y_true, y_prob, thr) {
  aucs <- auc_pair(y_true, y_prob)
  bm   <- bin_metrics(y_true, y_prob, thr)
  data.table::data.table(
    AUROC = aucs["auroc"], AUPRC = aucs["auprc"],
    Brier = brier_score(y_true, y_prob), Threshold = thr,
    Sensitivity = bm$Sensitivity, Specificity = bm$Specificity,
    PPV = bm$PPV, NPV = bm$NPV, F1 = bm$F1,
    TP = bm$TP, FP = bm$FP, TN = bm$TN, FN = bm$FN
  )
}

# =============================================================================
# Section 5: Patient-Level Bootstrap
# =============================================================================

#' Bootstrap external validation with cluster (patient-level) resampling
bootstrap_external <- function(dt_ext_pred, thr, B = B_BOOTSTRAP, seed = BOOTSTRAP_SEED,
                               cluster_col = CLUSTER_COL, dca_pts = DCA_PTS) {
  set.seed(seed)
  y <- as.integer(dt_ext_pred$y_true)
  p <- as.numeric(dt_ext_pred$y_prob)

  idx_by_id <- split(seq_len(nrow(dt_ext_pred)), as.character(dt_ext_pred[[cluster_col]]))
  uid   <- names(idx_by_id)
  n_uid <- length(uid)

  met_mat <- matrix(NA_real_, nrow = B, ncol = 8)
  colnames(met_mat) <- c("AUROC", "AUPRC", "Brier", "Sensitivity", "Specificity", "PPV", "NPV", "F1")
  dca_mat <- matrix(NA_real_, nrow = B, ncol = length(dca_pts))

  for (b in seq_len(B)) {
    boot_uid <- sample(uid, size = n_uid, replace = TRUE)
    idx <- unlist(idx_by_id[boot_uid], use.names = FALSE)
    yb <- y[idx]; pb <- p[idx]
    if (length(unique(yb)) < 2) next

    tryCatch({
      aucs <- auc_pair(yb, pb)
      bm   <- bin_metrics(yb, pb, thr)
      met_mat[b, ] <- c(aucs["auroc"], aucs["auprc"], brier_score(yb, pb),
                         bm$Sensitivity, bm$Specificity, bm$PPV, bm$NPV, bm$F1)
      dca_dt <- dca_curve(yb, pb, pts = dca_pts)
      dca_mat[b, ] <- dca_dt$NB_model
    }, error = function(e) NULL)
  }

  met_dt <- data.table::as.data.table(met_mat)
  ci_dt  <- data.table::rbindlist(lapply(names(met_dt), function(m) {
    v <- met_dt[[m]]; v <- v[!is.na(v)]
    if (length(v) < 100) return(data.table::data.table(metric = m, ci_low = NA_real_, ci_high = NA_real_))
    data.table::data.table(
      metric  = m,
      ci_low  = as.numeric(stats::quantile(v, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(stats::quantile(v, 0.975, na.rm = TRUE))
    )
  }))

  list(boot_metrics = met_dt, ci = ci_dt, dca_mat = dca_mat,
       valid_count = sum(!is.na(met_mat[, 1])))
}

# =============================================================================
# Section 6: NRI / IDI Reclassification Functions
# =============================================================================

#' Assign risk category based on NRI thresholds
assign_risk_cat <- function(p, thresholds = NRI_THRESHOLDS, labels = NRI_LABELS) {
  cut(p, breaks = c(-Inf, thresholds, Inf), labels = labels, right = FALSE)
}

#' Categorical NRI (3-class: Low / Medium / High)
compute_categorical_nri <- function(dt) {
  dt <- data.table::copy(dt)
  dt[, cat_pre  := assign_risk_cat(p_pre)]
  dt[, cat_post := assign_risk_cat(p_post)]
  dt[, moved      := as.integer(cat_pre != cat_post)]
  dt[, moved_up   := as.integer(as.numeric(cat_post) > as.numeric(cat_pre))]
  dt[, moved_down := as.integer(as.numeric(cat_post) < as.numeric(cat_pre))]

  n_ev <- dt[outcome == 1, .N]
  n_ne <- dt[outcome == 0, .N]

  ev_up <- dt[outcome == 1, sum(moved_up)]
  ev_dn <- dt[outcome == 1, sum(moved_down)]
  ne_up <- dt[outcome == 0, sum(moved_up)]
  ne_dn <- dt[outcome == 0, sum(moved_down)]

  nri_ev <- (ev_up - ev_dn) / n_ev
  nri_ne <- (ne_dn - ne_up) / n_ne

  list(NRI = nri_ev + nri_ne, NRI_events = nri_ev, NRI_nonevents = nri_ne,
       n_total = nrow(dt), n_events = n_ev,
       n_moved = dt[, sum(moved)], pct_moved = dt[, 100 * mean(moved)],
       events_up = ev_up, events_down = ev_dn,
       nonevents_up = ne_up, nonevents_down = ne_dn)
}

#' Continuous NRI (category-free)
compute_continuous_nri <- function(dt) {
  ev <- dt[outcome == 1]
  ne <- dt[outcome == 0]
  nri_ev <- mean(ev$p_post > ev$p_pre) - mean(ev$p_post < ev$p_pre)
  nri_ne <- mean(ne$p_post < ne$p_pre) - mean(ne$p_post > ne$p_pre)
  list(cNRI = nri_ev + nri_ne, cNRI_events = nri_ev, cNRI_nonevents = nri_ne)
}

#' Integrated Discrimination Improvement (IDI)
compute_idi <- function(dt) {
  is_post <- dt[outcome == 1, mean(p_post)] - dt[outcome == 0, mean(p_post)]
  is_pre  <- dt[outcome == 1, mean(p_pre)]  - dt[outcome == 0, mean(p_pre)]
  idi <- is_post - is_pre
  list(IDI = idi, IS_pre = is_pre, IS_post = is_post, relative_IDI = idi / is_pre,
       mean_p_change_events    = dt[outcome == 1, mean(p_post - p_pre)],
       mean_p_change_nonevents = dt[outcome == 0, mean(p_post - p_pre)])
}

#' Bootstrap CI for NRI/IDI (patient-level cluster bootstrap)
bootstrap_nri_idi <- function(dt, B = B_BOOTSTRAP, seed = BOOTSTRAP_SEED) {
  set.seed(seed)
  patients <- unique(dt$cluster_id)
  n_pat    <- length(patients)

  boot_dt <- data.table::data.table(
    b = seq_len(B),
    cat_NRI = NA_real_, cat_NRI_ev = NA_real_, cat_NRI_ne = NA_real_,
    cNRI = NA_real_, cNRI_ev = NA_real_, cNRI_ne = NA_real_,
    IDI = NA_real_, rel_IDI = NA_real_
  )

  for (i in seq_len(B)) {
    bp <- sample(patients, n_pat, replace = TRUE)
    bi <- data.table::data.table(cluster_id = bp, .boot_id = seq_along(bp))
    bd <- merge(bi, dt, by = "cluster_id", allow.cartesian = TRUE)

    cr <- tryCatch(compute_categorical_nri(bd), error = function(e) NULL)
    cn <- tryCatch(compute_continuous_nri(bd),   error = function(e) NULL)
    id <- tryCatch(compute_idi(bd),              error = function(e) NULL)

    if (!is.null(cr)) data.table::set(boot_dt, i, c("cat_NRI", "cat_NRI_ev", "cat_NRI_ne"),
                                       list(cr$NRI, cr$NRI_events, cr$NRI_nonevents))
    if (!is.null(cn)) data.table::set(boot_dt, i, c("cNRI", "cNRI_ev", "cNRI_ne"),
                                       list(cn$cNRI, cn$cNRI_events, cn$cNRI_nonevents))
    if (!is.null(id)) data.table::set(boot_dt, i, c("IDI", "rel_IDI"),
                                       list(id$IDI, id$relative_IDI))
    if (i %% 200 == 0) cat(sprintf("    bootstrap %d/%d\n", i, B))
  }

  ci <- function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  list(boot = boot_dt,
       CI = list(
         cat_NRI    = ci(boot_dt$cat_NRI),
         cat_NRI_ev = ci(boot_dt$cat_NRI_ev),
         cat_NRI_ne = ci(boot_dt$cat_NRI_ne),
         cNRI       = ci(boot_dt$cNRI),
         IDI        = ci(boot_dt$IDI),
         rel_IDI    = ci(boot_dt$rel_IDI)
       ))
}

# =============================================================================
# Section 7: Cross-site Heterogeneity Functions
# (used in 02_heterogeneity.R for Table 1B, Supp Tables 1–2)
# =============================================================================

#' Cliff's delta (rank-based, efficient for large n)
cliffs_delta_fast <- function(x, y, max_n = SAMPLE_CONFIG$cliff_max_n) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 5 || length(y) < 5) return(NA_real_)
  if (length(x) > max_n) x <- sample(x, max_n)
  if (length(y) > max_n) y <- sample(y, max_n)

  n1 <- as.double(length(x)); n2 <- as.double(length(y))  # as.double prevents integer overflow when n1*n2 > 2^31
  combined <- c(x, y)
  ranks <- rank(combined, ties.method = "average")
  R1 <- sum(ranks[1:n1])
  U1 <- R1 - n1 * (n1 + 1) / 2
  (2 * U1) / (n1 * n2) - 1
}

#' Cliff's delta magnitude classification (Romano et al. thresholds)
cliffs_delta_magnitude <- function(delta) {
  if (is.na(delta)) return(NA_character_)
  ad <- abs(delta)
  if (ad < 0.147) return("Negligible")
  if (ad < 0.33)  return("Small")
  if (ad < 0.474) return("Medium")
  return("Large")
}

# =============================================================================
# Section 8: File Discovery Helpers
# =============================================================================

#' Parse iteration directory name into structured info
parse_iter_info <- function(iter_dir) {
  tag  <- basename(iter_dir)
  ext  <- sub(".*External_([^_]+)_.*", "\\1", tag)
  seed <- sub(".*Seed([0-9]+)$", "\\1", tag)
  list(iter_tag = tag, external = ext, seed_id = seed)
}

#' Find the most recent file matching a pattern in a directory
find_latest_file <- function(dir_path, pattern) {
  fs <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (!length(fs)) return(NA_character_)
  fs[which.max(file.mtime(fs))]
}

#' Logging helper
log_msg <- function(..., level = "INFO") {
  cat(sprintf("[%s] %s %s\n", level, format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# =============================================================================
# Section 18: E/L Temporal Split — Shared temporal_filter()
# Used by: 04c, 04d, 04e, 04f (Plan A pipeline)
# =============================================================================

#' Temporal E/L partition filter
#'
#' Splits a site data.table into E (earlier) or L (later) based on split_table.
#' Requires dt to have columns: site and the configured session-date column.
#'
#' @param dt          data.table with a single site
#' @param partition   "E" or "L"
#' @param split_table data.table with columns: Center, Split_Date
#' @param verbose     if TRUE, print partition summary
#' @param date_col    session-date column name; defaults to SESSION_DATE_COL
#' @return filtered data.table
temporal_filter <- function(dt, partition, split_table, verbose = TRUE,
                            date_col = SESSION_DATE_COL) {
  site <- unique(dt$site)
  if (length(site) != 1) stop("temporal_filter(): dt must contain exactly one site")

  if (!(date_col %in% names(dt)))
    stop(sprintf("temporal_filter(): %s column not found in %s", date_col, site))

  split_date <- as.Date(split_table[Center == site, Split_Date])
  if (length(split_date) != 1 || is.na(split_date))
    stop(sprintf("temporal_filter(): no valid Split_Date for center %s", site))

  dt[, (date_col) := as.Date(get(date_col))]

  if (partition == "E") {
    out <- dt[get(date_col) <= split_date]
  } else if (partition == "L") {
    out <- dt[get(date_col) > split_date]
  } else {
    stop("temporal_filter(): partition must be 'E' or 'L'")
  }

  if (isTRUE(verbose)) {
    cat(sprintf("  [%s] %s: %d/%d rows (split @ %s)\n",
                site, partition, nrow(out), nrow(dt), as.character(split_date)))
  }
  out
}

# =============================================================================
# Section 19: E-only Data Loading Helpers (v3.0 Pure E-only Design)
# Used by: 03, 04, 08, 09 — any script that loads raw data and needs E-only
# =============================================================================

#' Apply per-site chronological E/L split
#'
#' Adds a cohort column ("E" or "L") to dt based on per-site date quantile.
#' Sessions on or before the quantile cutpoint → "E"; after → "L".
#'
#' @param dt             data.table with 'site' column
#' @param date_col       name of the date column (e.g., "HD_Date")
#' @param split_quantile numeric, e.g. 0.75 for 75% E / 25% L
#' @param cohort_col     name of the output cohort column
#' @return dt with cohort_col added (modified by reference)
apply_el_split <- function(dt, date_col, split_quantile, cohort_col) {
  if (!(date_col %in% names(dt)))
    stop(sprintf("apply_el_split(): column '%s' not found. Check SESSION_DATE_COL in 00_config.R.", date_col))

  if (!inherits(dt[[date_col]], c("Date", "POSIXct", "IDate")))
    dt[, (date_col) := as.Date(get(date_col))]

  dt[, (cohort_col) := {
    dates <- get(date_col)
    cutpoint <- quantile(dates, probs = split_quantile, type = 1, na.rm = TRUE)
    ifelse(dates <= cutpoint, "E", "L")
  }, by = site]

  dt
}

#' Load a site, preprocess, apply E/L split, and return E-only data
#'
#' Convenience wrapper for scripts that need E-only raw data.
#' Requires SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL from 00_config.R
#'
#' @param site_name  character, e.g. "TN"
#' @param verbose    if TRUE, print E/L counts
#' @return data.table (E-cohort only)
load_site_eonly <- function(site_name, verbose = TRUE) {
  dt <- prep_site_dt(read_fst_dt(SITE_FILES[[site_name]]), site_name)
  dt <- apply_el_split(dt, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
  n_all <- nrow(dt)
  dt <- dt[get(COHORT_COL) == "E"]
  if (verbose) {
    cat(sprintf("    %s: %d total → %d E-only (%.1f%%)\n",
                site_name, n_all, nrow(dt), nrow(dt) / n_all * 100))
  }
  dt
}

cat("done.\n")

# =============================================================================
# Section 20: Gold Standard Spline Smoothing (shared by 04f, 07)
# Ported from Post-hoc v2.1, including 3-layer domain restriction.
# =============================================================================

#' Gold-standard spline smoothing with domain restriction
#'
#' Three-layer domain restriction:
#'   L1: Fitting only within phys_range
#'   L2: Updating only within phys_range
#'   L3: Clipping to original score range ± margin
#'
#' @param x         Bin centroid x-values
#' @param y_old     Original bin scores
#' @param w         Bin weights (from ebm.bin_weights_)
#' @param spar      Spline smoothing parameter
#' @param phys_range Physiological range [lo, hi] for this feature
#' @param update_only_phys If TRUE, restrict fitting & update to phys_range
#' @param clip_to_old If TRUE, clip smoothed values to original range ± margin
#' @param clip_margin Fraction of original range to extend clip bounds
#' @return list with y_new, fit_ok, n_fit, msg
smooth_scores_spline_fixed <- function(x, y_old, w, spar = 0.6,
                                        phys_range = NULL,
                                        update_only_phys = TRUE,
                                        clip_to_old = TRUE,
                                        clip_margin = 0.25) {
  n <- length(y_old)
  if (length(x) != n || length(w) != n)
    stop("length mismatch in smooth_scores_spline_fixed()")

  y_new <- as.numeric(y_old)

  # L1: Determine fitting indices
  idx_fit <- which(is.finite(x) & is.finite(y_old) & is.finite(w) & w > 0)

  if (!is.null(phys_range) && isTRUE(update_only_phys)) {
    idx_fit <- idx_fit[x[idx_fit] >= phys_range[1] & x[idx_fit] <= phys_range[2]]
  }
  if (length(idx_fit) < 5) {
    return(list(y_new = y_new, fit_ok = FALSE, n_fit = length(idx_fit),
                msg = "Insufficient data points for spline fitting"))
  }

  # Weighted aggregation by unique x (for smooth.spline)
  dt <- data.table::data.table(x = x[idx_fit], y = y_old[idx_fit], w = w[idx_fit])
  dt_u <- dt[, .(y = weighted.mean(y, w, na.rm = TRUE),
                  w = sum(w, na.rm = TRUE)), by = x][order(x)]
  if (nrow(dt_u) < 5) {
    return(list(y_new = y_new, fit_ok = FALSE, n_fit = nrow(dt_u),
                msg = "Insufficient unique x values for spline fitting"))
  }

  # Fit spline
  fit <- tryCatch(
    stats::smooth.spline(dt_u$x, dt_u$y, w = dt_u$w, spar = spar),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(y_new = y_new, fit_ok = FALSE, n_fit = nrow(dt_u),
                msg = "smooth.spline failed"))
  }

  # L2: Determine update indices
  idx_update <- which(is.finite(x) & is.finite(y_old))
  if (!is.null(phys_range) && isTRUE(update_only_phys)) {
    idx_update <- idx_update[x[idx_update] >= phys_range[1] &
                              x[idx_update] <= phys_range[2]]
  }
  if (length(idx_update) > 0) {
    y_hat <- stats::predict(fit, x[idx_update])$y
    y_new[idx_update] <- as.numeric(y_hat)
  }

  # L3: Clip to original range ± margin
  if (isTRUE(clip_to_old)) {
    valid_for_clip <- is.finite(x) & is.finite(y_old)
    if (sum(valid_for_clip) >= 2) {
      yr <- range(y_old[valid_for_clip], na.rm = TRUE)
      if (all(is.finite(yr))) {
        span <- diff(yr)
        if (!is.finite(span) || span <= 0) span <- abs(yr[1]) * 0.1 + 0.1
        lo <- yr[1] - clip_margin * span
        hi <- yr[2] + clip_margin * span
        y_new <- pmin(pmax(y_new, lo), hi)
      }
    }
  }

  list(y_new = y_new, fit_ok = TRUE, n_fit = nrow(dt_u), msg = "OK")
}

#' Load SPAR map from 05_shapeqc.R output
#'
#' Reads Table_4A_Decision_Matrix.csv and computes spar via SPAR_FORMULA(J).
#' @param run_dir  RUN_DIR path
#' @param spar_fn  Function mapping J → spar (default: SPAR_FORMULA from config)
#' @return data.table with columns: feature, tier, median_J, spar
load_spar_map <- function(run_dir, spar_fn = SPAR_FORMULA) {
  qc_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  qc_dirs <- sort(qc_dirs[grepl("ShapeQC", basename(qc_dirs))], decreasing = TRUE)
  if (length(qc_dirs) == 0)
    stop("[load_spar_map] No ShapeQC output found in RUN_DIR. Run 05_shapeqc.R first.")
  qc_csv <- file.path(qc_dirs[1], "Table_4A_Decision_Matrix.csv")
  if (!file.exists(qc_csv))
    stop(sprintf("[load_spar_map] %s not found.", qc_csv))
  dt <- data.table::fread(qc_csv)
  if ("light" %in% names(dt)) data.table::setnames(dt, "light", "tier")
  if (!("median_J" %in% names(dt)))
    stop("[load_spar_map] median_J column not found in ShapeQC decision matrix.")
  dt[, spar := sapply(median_J, spar_fn)]
  dt[, .(feature, tier, median_J, spar)]
}
