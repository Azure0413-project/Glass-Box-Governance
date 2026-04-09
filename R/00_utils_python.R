# =============================================================================
# 00_utils_python.R — Python/Reticulate Helper Functions
# IDH EBM Governance Reproducibility Pipeline
#
# All functions for interfacing with InterpretML EBM models via reticulate.
# Includes: model prediction, shape function extraction, bin parsing,
#        post-hoc term update (smoothing & zeroing), and mean-preserve verification.
#
# Source: consolidated from original analysis scripts
# These functions were identical across the original post-hoc scripts; deduplicated here.
#
# Usage: src("R/00_config.R"); source("00_utils_r.R"); source("00_utils_python.R")
#
# Note: this script calls reticulate::use_virtualenv() to initialize the Python environment.
#        Modules that do not require EBM (01/02) can skip this file.
# =============================================================================

cat("  Loading 00_utils_python.R ... ")

# =============================================================================
# Section 1: Python Environment Setup
# =============================================================================

init_python_env <- function(venv_path = VENV_PATH) {
  if (!dir.exists(venv_path)) {
    stop(sprintf("[Error] Python venv not found: %s\n  Set R_INTERPRET_VENV env var or edit VENV_PATH in 00_config.R", venv_path))
  }
  reticulate::use_virtualenv(venv_path, required = TRUE)
  cat(sprintf("  Python env: %s\n", reticulate::py_config()$python))
  invisible(TRUE)
}

# =============================================================================
# Section 2: Import Python Modules
# =============================================================================

init_python_modules <- function() {
  np     <<- reticulate::import("numpy",   convert = FALSE)
  joblib <<- reticulate::import("joblib",  convert = FALSE)
  interpret_glassbox <<- reticulate::import("interpret.glassbox", convert = FALSE)
  cat("  Python modules: numpy, joblib, interpret.glassbox\n")
  invisible(TRUE)
}

# =============================================================================
# Section 3: Core Python Functions (embedded via py_run_string)
#
# These functions directly operate on the internal structure of InterpretML EBM objects:
#   - predict_proba_pos: extract positive-class predicted probability
#   - get_intercept: retrieve the model intercept
#   - update_term_scores_with_shift: update shape function scores with automatic mean preservation
#   - extract_raw_bins_info: extract bin boundaries, scores, and weights
# =============================================================================

init_python_functions <- function() {

reticulate::py_run_string("
import numpy as np

def _wmean(x, w):
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(m):
        return np.nan
    return float(np.sum(x[m] * w[m]) / np.sum(w[m]))

def get_positive_index(classes_py):
    candidates = [1, '1', 'positive', 'pos', True, 'True', 'true', 'yes', 'Y']
    for i, c in enumerate(classes_py):
        if c in candidates:
            return int(i)
    if len(classes_py) == 2:
        return 1
    raise ValueError('Cannot determine positive class index')

def predict_proba_pos(model, X_df):
    try:
        import pandas as pd
        if not isinstance(X_df, pd.DataFrame):
            X_df = pd.DataFrame(X_df)
    except Exception:
        pass
    proba = model.predict_proba(X_df)
    classes = list(model.classes_)
    pos_idx = get_positive_index(classes)
    return np.asarray(proba)[:, pos_idx].astype(float)

def get_intercept(ebm):
    try:
        return float(ebm.intercept_[0])
    except Exception:
        try:
            return float(np.asarray(ebm.intercept_, dtype=float).flat[0])
        except Exception:
            return float('nan')

def update_term_scores_with_shift(ebm, term_index, new_scores,
                                   mean_preserve=True,
                                   preserve_missing=True,
                                   missing_index=-1,
                                   preserve_unknown=True,
                                   unknown_index=None,
                                   shift_bagged_scores=True,
                                   shift_bagged_intercept=True):
    old = np.asarray(ebm.term_scores_[term_index], dtype=float)
    new = np.asarray(new_scores, dtype=float)
    if old.shape != new.shape:
        raise ValueError(f'shape mismatch: old {old.shape} vs new {new.shape}')

    # Preserve missing bin
    if preserve_missing:
        try: mi = int(missing_index)
        except: mi = -1
        if mi >= 0 and mi < old.size:
            new = new.copy(); new.flat[mi] = old.flat[mi]

    # Preserve unknown bin
    if preserve_unknown and unknown_index is not None:
        try: ui = int(unknown_index)
        except: ui = -1
        if ui >= 0 and ui < old.size:
            new = new.copy(); new.flat[ui] = old.flat[ui]

    # Mean-preserve: adjust intercept so population-average predicted risk unchanged
    delta = 0.0
    if mean_preserve and hasattr(ebm, 'bin_weights_'):
        w = np.asarray(ebm.bin_weights_[term_index], dtype=float)
        if w.shape == old.shape:
            mu_old = _wmean(old, w)
            mu_new = _wmean(new, w)
            if np.isfinite(mu_old) and np.isfinite(mu_new):
                delta = float(mu_old - mu_new)
                try: ebm.intercept_[0] = float(ebm.intercept_[0]) + delta
                except:
                    inter = np.asarray(ebm.intercept_, dtype=float).copy()
                    inter.flat[0] = float(inter.flat[0]) + delta
                    ebm.intercept_ = inter
                if shift_bagged_intercept and hasattr(ebm, 'bagged_intercept_'):
                    try:
                        bi = np.asarray(ebm.bagged_intercept_, dtype=float)
                        ebm.bagged_intercept_ = bi + delta
                    except: pass

    # Shift bagged scores in sync
    if shift_bagged_scores and hasattr(ebm, 'bagged_scores_'):
        try:
            bag = np.asarray(ebm.bagged_scores_[term_index], dtype=float)
            d = new - old
            if bag.ndim == 2 and old.ndim == 1 and bag.shape[1] == old.size:
                ebm.bagged_scores_[term_index] = bag + d.reshape(1, -1)
            elif bag.ndim == old.ndim and bag.shape == old.shape:
                ebm.bagged_scores_[term_index] = bag + d
            elif bag.ndim > old.ndim:
                ebm.bagged_scores_[term_index] = bag + d.reshape((1,) * (bag.ndim - old.ndim) + old.shape)
        except: pass

    ebm.term_scores_[term_index] = new
    return float(delta)

def _safe_float_array(a):
    try:
        arr = np.asarray(a, dtype=float).ravel()
        arr = arr[np.isfinite(arr)]
        return arr if arr.size > 0 else None
    except: return None

def _first_ndarray(obj):
    if isinstance(obj, np.ndarray): return obj
    if isinstance(obj, (list, tuple)):
        for it in obj:
            found = _first_ndarray(it)
            if found is not None: return found
    return None

def extract_raw_bins_info(ebm, feature_name):
    try:
        t_names = list(ebm.term_names_)
        t_idx = t_names.index(feature_name)
    except: return None

    # Find feature index
    f_idx = None
    if hasattr(ebm, 'term_features_') and ebm.term_features_ is not None:
        try:
            tf = ebm.term_features_[t_idx]
            if isinstance(tf, np.ndarray): tf = tf.tolist()
            if isinstance(tf, (list, tuple)) and len(tf) >= 1: f_idx = int(tf[0])
            elif tf is not None: f_idx = int(tf)
        except: f_idx = None
    if f_idx is None:
        for attr in ['feature_names_in_', 'feature_names', 'feature_names_']:
            if hasattr(ebm, attr):
                try:
                    f_names = list(getattr(ebm, attr))
                    if feature_name in f_names: f_idx = f_names.index(feature_name); break
                except: continue

    # Feature type
    feature_type = None
    if f_idx is not None:
        for attr in ['feature_types_in_', 'feature_types', 'feature_types_']:
            if hasattr(ebm, attr):
                try:
                    ft = getattr(ebm, attr)
                    if len(ft) > f_idx: feature_type = str(ft[f_idx]); break
                except: continue

    scores  = np.asarray(ebm.term_scores_[t_idx], dtype=float).flatten()
    weights = np.asarray(ebm.bin_weights_[t_idx], dtype=float).flatten()

    # Feature bounds
    bounds = []
    if f_idx is not None and hasattr(ebm, 'feature_bounds_') and ebm.feature_bounds_ is not None:
        try:
            b = np.asarray(ebm.feature_bounds_[f_idx], dtype=float).ravel()
            if b.size == 2: bounds = [float(b[0]), float(b[1])]
        except: pass

    # Extract cuts from bins_
    cuts = None; extraction_method = 'unknown'; bins_summary = ''
    if f_idx is not None and hasattr(ebm, 'bins_'):
        try:
            raw_bins = ebm.bins_[f_idx]
            try: bins_len = len(raw_bins)
            except: bins_len = None
            bins_summary = f'type={type(raw_bins).__name__}; len={bins_len}'
            arr = _first_ndarray(raw_bins)
            if arr is not None:
                arr2 = _safe_float_array(arr)
                if arr2 is not None: cuts = arr2; extraction_method = f'bins_[{f_idx}]'
        except: pass

    # Fallback: preprocessor_
    if cuts is None and f_idx is not None and hasattr(ebm, 'preprocessor_') and ebm.preprocessor_ is not None:
        try:
            if hasattr(ebm.preprocessor_, 'col_bin_edges_'):
                edges = _safe_float_array(ebm.preprocessor_.col_bin_edges_[f_idx])
                if edges is not None: cuts = edges; extraction_method = f'preprocessor_.col_bin_edges_[{f_idx}]'
        except: pass

    cuts_list = []
    if cuts is not None:
        c = np.unique(np.asarray(cuts, dtype=float).ravel())
        c = c[np.isfinite(c)]
        cuts_list = c.tolist()

    return {
        'term_index': int(t_idx),
        'feature_index': int(f_idx) if f_idx is not None else -1,
        'feature_type': feature_type if feature_type is not None else '',
        'feature_bounds': bounds,
        'n_scores': int(len(scores)),
        'scores': scores.tolist(),
        'weights': weights.tolist(),
        'cuts': cuts_list,
        'extraction_method': extraction_method,
        'bins_summary': bins_summary
    }
")

  # Export Python functions as R-accessible handles
  py_predict_proba_pos <<- reticulate::py$predict_proba_pos
  py_update_term       <<- reticulate::py$update_term_scores_with_shift
  py_extract_raw_bins  <<- reticulate::py$extract_raw_bins_info
  py_get_intercept     <<- reticulate::py$get_intercept

  cat("  Python EBM functions: loaded\n")
  invisible(TRUE)
}

# =============================================================================
# Section 4: R-side Wrappers for Python Functions
# =============================================================================

#' Predict positive-class probability via Python EBM
predict_proba_pos_r <- function(mdl, X_df) {
  as.numeric(reticulate::py_to_r(py_predict_proba_pos(mdl, X_df)))
}

#' Get EBM intercept
get_intercept_r <- function(mdl) {
  as.numeric(reticulate::py_to_r(py_get_intercept(mdl)))
}

# =============================================================================
# Section 5: R-side EBM Builder (for training, used in 03_iecv_ebm.R)
# =============================================================================

#' Construct an EBM classifier from hyperparameters
build_ebm <- function(hp, n_jobs = EBM_N_JOBS, glassbox_module = interpret_glassbox) {
  args <- list(
    interactions     = as.integer(hp$interactions),
    random_state     = as.integer(hp$random_state),
    learning_rate    = as.numeric(hp$learning_rate),
    max_bins         = as.integer(hp$max_bins),
    max_leaves       = as.integer(hp$max_leaves),
    min_samples_leaf = as.integer(hp$min_samples_leaf),
    smoothing_rounds = as.integer(hp$smoothing_rounds),
    outer_bags       = as.integer(hp$outer_bags),
    n_jobs           = as.integer(n_jobs)
  )
  tryCatch(
    do.call(glassbox_module$ExplainableBoostingClassifier, args),
    error = function(e) {
      args$n_jobs <- NULL
      do.call(glassbox_module$ExplainableBoostingClassifier, args)
    }
  )
}

#' Get term importances (multiple attribute name fallbacks)
get_term_importances <- function(ebm_model) {
  imp <- tryCatch({
    result <- ebm_model$term_importances()
    reticulate::py_to_r(result)
  }, error = function(e) NULL)
  if (!is.null(imp) && is.numeric(imp)) return(imp)

  imp <- tryCatch({
    result <- ebm_model$term_importances_
    reticulate::py_to_r(result)
  }, error = function(e) NULL)
  if (!is.null(imp) && is.numeric(imp)) return(imp)

  imp <- tryCatch({
    result <- ebm_model$feature_importances_
    reticulate::py_to_r(result)
  }, error = function(e) NULL)
  if (!is.null(imp) && is.numeric(imp)) return(imp)

  NULL
}

#' R-side positive class index detection (for parallel workers)
#' This version works without Python helpers — used inside parLapply workers
get_positive_index <- function(classes_py) {
  candidates <- c(1, "1", "positive", "pos", "True", "true", "yes", "Y")
  hit <- which(classes_py %in% candidates)
  if (length(hit)) return(hit[1])
  if (length(classes_py) == 2) return(2L)
  stop("Cannot determine positive class index")
}

#' R-side predict_proba for positive class (for PSOCK parallel workers)
#' Calls mdl$predict_proba() directly via reticulate — no Python helper needed.
#' NOTE: This function is distinct from predict_proba_pos_r() (Section 4) which
#' wraps the Python-side py_predict_proba_pos helper. Use this version inside
#' parLapply workers where init_python_functions() has not been called.
predict_proba_pos <- function(mdl, X) {
  proba   <- mdl$predict_proba(X)
  proba_r <- reticulate::py_to_r(proba)
  classes <- reticulate::py_to_r(mdl$classes_)
  pos_idx <- get_positive_index(classes)
  as.numeric(proba_r[, pos_idx])
}

# =============================================================================
# Section 6: Bins Extraction & Centroid Computation
# =============================================================================

#' Compute bin centroids from cut points
compute_centroids_from_cuts <- function(n_scores, cuts, bounds = NULL, verbose = FALSE) {
  x_vec <- rep(NA_real_, n_scores)
  align_method <- "unknown"
  missing_idx <- NA_integer_; unknown_idx <- NA_integer_; x_is_index <- FALSE

  cuts <- sort(unique(as.numeric(cuts[is.finite(as.numeric(cuts))])))
  n_cuts <- length(cuts)

  lo <- hi <- NA_real_
  if (!is.null(bounds) && length(bounds) == 2 && all(is.finite(bounds))) {
    lo <- bounds[1]; hi <- bounds[2]
  }

  if (n_cuts == 0) {
    x_is_index <- TRUE; align_method <- "no_cuts_fallback_index"
    if (n_scores >= 2) {
      missing_idx <- 0L
      x_vec[2:n_scores] <- seq_len(n_scores - 1L)
    }
    return(list(x_vec = x_vec, align_method = align_method, n_with_x = sum(is.finite(x_vec)),
                missing_idx = missing_idx, unknown_idx = unknown_idx, x_is_index = x_is_index))
  }

  mids <- if (n_cuts > 1) (cuts[-n_cuts] + cuts[-1]) / 2 else numeric(0)
  avg_width <- if (n_cuts > 1) mean(diff(cuts), na.rm = TRUE) else 1
  if (!is.finite(avg_width) || avg_width <= 0) avg_width <- 1

  first_center <- if (is.finite(lo)) (lo + cuts[1]) / 2 else cuts[1] - avg_width / 2
  last_center  <- if (is.finite(hi)) (cuts[n_cuts] + hi) / 2 else cuts[n_cuts] + avg_width / 2
  x_centers <- c(first_center, mids, last_center)
  n_centers <- length(x_centers)

  if (n_scores == n_centers + 2) {
    missing_idx <- 0L; unknown_idx <- n_scores - 1L
    x_vec[2:(n_scores - 1L)] <- x_centers
    align_method <- sprintf("two_special_bins (cuts=%d, centers=%d)", n_cuts, n_centers)
  } else if (n_scores == n_centers + 1) {
    missing_idx <- 0L
    x_vec[2:n_scores] <- x_centers
    align_method <- sprintf("standard_missing_idx0 (cuts=%d, centers=%d)", n_cuts, n_centers)
  } else if (n_scores == n_centers) {
    x_vec <- x_centers
    align_method <- sprintf("no_special_bin (cuts=%d, centers=%d)", n_cuts, n_centers)
  } else {
    align_method <- sprintf("mismatch_fallback (n_scores=%d, cuts=%d, centers=%d)",
                            n_scores, n_cuts, n_centers)
    missing_idx <- 0L
    n_fill <- min(n_centers, n_scores - 1L)
    if (n_fill > 0) x_vec[2:(1 + n_fill)] <- x_centers[1:n_fill]
  }

  list(x_vec = x_vec, align_method = align_method, n_with_x = sum(is.finite(x_vec)),
       missing_idx = missing_idx, unknown_idx = unknown_idx, x_is_index = x_is_index)
}

#' Full bins extraction with diagnostics
extract_bins_with_diagnostics <- function(ebm, feature_name, strict_continuous_cuts = TRUE, verbose = FALSE) {
  raw_info <- reticulate::py_to_r(py_extract_raw_bins(ebm, feature_name))
  if (is.null(raw_info))
    return(list(success = FALSE, error_msg = sprintf("Feature '%s' not found", feature_name)))

  n_scores     <- raw_info$n_scores
  scores       <- as.numeric(raw_info$scores)
  weights      <- as.numeric(raw_info$weights)
  cuts         <- as.numeric(raw_info$cuts)
  feature_type <- as.character(raw_info$feature_type)
  bounds       <- as.numeric(raw_info$feature_bounds)
  if (length(bounds) != 2 || !all(is.finite(bounds))) bounds <- NULL

  if (isTRUE(strict_continuous_cuts) && nzchar(feature_type) &&
      grepl("continuous", feature_type, ignore.case = TRUE) && length(cuts) == 0)
    return(list(success = FALSE, error_msg = sprintf("Continuous '%s' has 0 cuts", feature_name)))

  cr <- compute_centroids_from_cuts(n_scores, cuts, bounds, verbose)

  dt <- data.table::data.table(
    bin_index  = 0:(n_scores - 1L),
    x          = cr$x_vec,
    score      = scores,
    bin_weight = weights
  )
  dt[, is_missing := if (is.finite(cr$missing_idx)) (bin_index == cr$missing_idx) else FALSE]
  dt[, is_unknown := if (is.finite(cr$unknown_idx)) (bin_index == cr$unknown_idx) else FALSE]

  list(success          = TRUE,
       dt               = dt,
       term_index       = as.integer(raw_info$term_index),
       feature_index    = as.integer(raw_info$feature_index),
       feature_type     = feature_type,
       feature_bounds   = bounds,
       x_is_index       = cr$x_is_index,
       n_scores         = n_scores,
       n_cuts           = length(cuts),
       n_with_x         = cr$n_with_x,
       n_with_weight    = sum(is.finite(weights) & weights > 0),
       extraction_method = raw_info$extraction_method,
       align_method      = cr$align_method,
       missing_idx       = cr$missing_idx,
       unknown_idx       = cr$unknown_idx)
}

# =============================================================================
# Section 7: Mean-Preserve Verification
# =============================================================================

#' Verify intercept shift consistency after term update
verify_mean_preserve_consistency <- function(ebm_old, ebm_new, accumulated_shift,
                                              tolerance = MEAN_PRESERVE_TOLERANCE) {
  intercept_old <- get_intercept_r(ebm_old)
  intercept_new <- get_intercept_r(ebm_new)
  actual_shift  <- intercept_new - intercept_old
  discrepancy   <- actual_shift - accumulated_shift

  list(intercept_old     = intercept_old,
       intercept_new     = intercept_new,
       actual_shift      = actual_shift,
       accumulated_shift = accumulated_shift,
       discrepancy       = discrepancy,
       tolerance         = tolerance,
       is_consistent     = abs(discrepancy) < tolerance,
       status            = ifelse(abs(discrepancy) < tolerance, "PASS", "WARN"))
}

# =============================================================================
# Section 8: Safe Python Object Conversion Helpers
# (v2.5.8 — handles list/array/matrix returns from Python)
# =============================================================================

#' Safely convert Python object to numeric vector
safe_to_numeric <- function(x) {
  if (is.null(x)) return(NULL)
  tryCatch({
    if (is.list(x))   x <- unlist(x, recursive = TRUE, use.names = FALSE)
    if (is.array(x) || is.matrix(x)) x <- as.vector(x)
    result <- as.numeric(x)
    if (length(result) == 0 || all(is.na(result))) return(NULL)
    result
  }, error = function(e) NULL)
}

#' Safely extract Python attribute by trying multiple names
safe_extract_py_attr <- function(py_obj, attr_names) {
  for (attr_name in attr_names) {
    result <- tryCatch({
      raw <- py_obj[[attr_name]]
      if (is.null(raw)) next
      converted <- reticulate::py_to_r(raw)
      safe_to_numeric(converted)
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
  }
  NULL
}

# =============================================================================
# Section 9: One-Shot Initialization
# =============================================================================

#' Initialize full Python environment + modules + functions
#' Call this once at the start of any script that needs EBM access.
init_python_full <- function(venv_path = VENV_PATH) {
  init_python_env(venv_path)
  init_python_modules()
  init_python_functions()
  cat("  Python environment fully initialized.\n")
  invisible(TRUE)
}

cat("done.\n")
