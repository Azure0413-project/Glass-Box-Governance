# =============================================================================
# 05_shapeqc.R — Shape Function Post-Hoc QC v3.0
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Table 4 Panel A — Feature-level transportability: C_core, J, S_seed, tier
#   Supp Table S3   — Cross-model shape-function correlation summary
#   Supp Table S14  — C_core sensitivity across 5 trimming configurations
#   Fig 3           — Governance quadrant scatter + individual shape overlays
#
# Source: consolidated from original analysis scripts
# Core logic:
#   Step 1: Load shape data -> 150-pt common grid projection -> mean-centering
#   Step 2: C_core (trimmed cross-fold concordance), J (jaggedness), S_seed
#   Step 3: Bootstrap CI for C_core
#   Step 4: 4-tier decision matrix (Green/Yellow/Red/Gray)
#   Step 5: Visualizations
#   Step 6: Excel report + Supp_Sensitivity (T1-T5)
#
# Methods → "Cross-fold core concordance (C_core) ... pairwise Spearman rank
#  correlations restricted to the trimmed domain"
# Methods → "Shape jaggedness (J) ... normalized total-variation index"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 03_iecv_ebm.R (needs shape function exports)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(ggplot2)
  library(openxlsx)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")
HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 05: ShapeQC v3.0\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# Local config (mirrors original v3.0 exactly)
# =============================================================================

G_GRID          <- SHAPEQC_N_GRID   # 150
EPSILON         <- 1e-8
BINARY_FEATS    <- c("Sex")
TRIM_P_LO      <- TRIM_CONFIGS[[PRIMARY_TRIM]][1]  # 0.005
TRIM_P_HI      <- TRIM_CONFIGS[[PRIMARY_TRIM]][2]  # 0.995
C0_THRESHOLD    <- SHAPEQC_TAU      # 0.70
J_THRESHOLD     <- SHAPEQC_KAPPA    # 1.50
BOOT_N          <- SHAPEQC_BOOT_N   # 2000
BOOT_CI         <- 0.95
BOOT_ALPHA      <- (1 - BOOT_CI) / 2
MIN_EFFECTIVE_N <- SHAPEQC_EFF_N_MIN  # from 00_config.R (=10); tied-rank safeguard threshold
MIN_GRID_FOR_CORR <- 30L
DISCRETE_THR    <- 15L

MODEL_STRUCTURE <- data.table(
  model_tag     = c("Iter1_External_TN_Seed1", "Iter2_External_TN_Seed2",
                     "Iter3_External_D6_Seed1", "Iter4_External_D6_Seed2",
                     "Iter5_External_CY_Seed1", "Iter6_External_CY_Seed2"),
  fold = c(1L, 1L, 2L, 2L, 3L, 3L),
  seed = c(1L, 2L, 1L, 2L, 1L, 2L),
  external_site = c("TN", "TN", "D6", "D6", "CY", "CY")
)

FOLD_TRAIN_SITES <- list("1" = c("D6", "CY"), "2" = c("TN", "CY"), "3" = c("TN", "D6"))

FEATURE_RANK_ORDER <- c(
  "IDH_N_28D", "Pre_HD_SBP", "UF_BW_Perc", "Target_UF_Volume",
  "Blood_Flow_Rate", "IDH_N_7D", "Heart_Rate", "Dry_Weight",
  "Age", "Pre_HD_Weight", "Start_DBP", "Body_Temperature",
  "Sex", "Dialysate_Temperature", "Dialysate_Flow_Rate", "Respiratory_Rate"
)

QC_TIMESTAMP <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
QC_OUT <- file.path(RUN_DIR, paste0("ShapeQC_v3.0_", QC_TIMESTAMP))
dir.create(QC_OUT, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  Grid: %d pts | Trim: T3 [P%.1f, P%.1f]\n", G_GRID, TRIM_P_LO*100, TRIM_P_HI*100))
cat(sprintf("  Decision: C0 >= %.2f & J <= %.2f\n", C0_THRESHOLD, J_THRESHOLD))
cat(sprintf("  Bootstrap: %d iters | Output: %s\n\n", BOOT_N, basename(QC_OUT)))

# =============================================================================
# Helper Functions
# =============================================================================

#' Spearman with tied-rank protection
safe_spearman_v2 <- function(x, y, min_eff_n = MIN_EFFECTIVE_N) {
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]; y <- y[valid]
  if (length(x) < 5) return(list(rho = NA_real_, tied_rank_flag = TRUE))
  eff_n_x <- length(unique(rank(x)))
  eff_n_y <- length(unique(rank(y)))
  eff_n   <- min(eff_n_x, eff_n_y)
  if (eff_n < min_eff_n) {
    rho <- suppressWarnings(cor(x, y, method = "spearman"))
    return(list(rho = rho, tied_rank_flag = TRUE, eff_n = eff_n))
  }
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  list(rho = rho, tied_rank_flag = FALSE, eff_n = eff_n)
}

safe_spearman <- function(x, y) safe_spearman_v2(x, y)$rho

#' Total-variation-based jaggedness: J = TV(f) / (Range(f) + eps) - 1
calc_jaggedness <- function(f) {
  f <- f[is.finite(f)]
  if (length(f) < 3) return(NA_real_)
  tv    <- sum(abs(diff(f)))
  rng   <- diff(range(f))
  tv / (rng + EPSILON) - 1
}

#' Step function from shape data
make_shape_stepfun <- function(dt_model) {
  dt_model <- dt_model[is.finite(x) & is.finite(score)][order(x)]
  if (nrow(dt_model) < 2) return(NULL)
  approxfun(dt_model$x, dt_model$score, method = "constant", rule = 2, ties = mean)
}

# =============================================================================
# Step 1: Load shape data + common grid projection
# =============================================================================

cat("[Step 1] Loading shape data & projecting onto common grid ...\n")

shape_all <- rbindlist(lapply(MODEL_STRUCTURE$model_tag, function(tag) {
  # Try xlsx first, then csv
  fp <- file.path(RUN_DIR, tag, "Shapes", "feature_data.xlsx")
  if (!file.exists(fp)) fp <- file.path(RUN_DIR, tag, "Shapes", "feature_data.csv")
  if (!file.exists(fp)) { warning(sprintf("  Missing: %s", tag)); return(data.table()) }
  dt <- if (grepl("\\.xlsx$", fp)) as.data.table(openxlsx::read.xlsx(fp)) else fread(fp)
  dt[, model := tag]
  dt
}), fill = TRUE)

shape_all <- merge(shape_all, MODEL_STRUCTURE[, .(model_tag, fold, seed)],
                   by.x = "model", by.y = "model_tag", all.x = TRUE)

cont_features <- setdiff(unique(shape_all$feature), BINARY_FEATS)
cont_features <- cont_features[nzchar(cont_features)]

# Load site data for pooled quantile computation
# NOTE: Uses the full temporal cohort (E+L) intentionally for trim bound
# computation. Trim bounds define the clinically plausible feature range
# (population reference range), not a model-training quantity.
# This is an analytic-frame choice, not predictive leakage: L-cohort
# outcomes are never used, and trim quantiles (P0.5–P99.5) are insensitive
# to the additional 25% of sessions. If reviewers require strict E-only
# purity, replace the loop below with load_site_eonly() per site.
site_data_cache <- list()
for (s in names(SITE_FILES)) {
  fp <- SITE_FILES[[s]]
  if (!file.exists(fp)) next
  dt_raw <- read_fst_dt(fp)
  dt_raw <- apply_compat_rename(dt_raw, COMPAT_RENAME_MAP)
  site_data_cache[[s]] <- dt_raw
  cat(sprintf("  %s: %s rows\n", s, format(nrow(dt_raw), big.mark = ",")))
}

# Build pooled data cache per fold
pooled_cache <- list()
pooled_all   <- list()
for (feat in cont_features) {
  pooled_cache[[feat]] <- list()
  all_vals <- numeric(0)
  for (fold_id in names(FOLD_TRAIN_SITES)) {
    fold_vals <- numeric(0)
    for (s in FOLD_TRAIN_SITES[[fold_id]]) {
      if (!is.null(site_data_cache[[s]]) && feat %in% names(site_data_cache[[s]])) {
        v <- as.numeric(site_data_cache[[s]][[feat]])
        fold_vals <- c(fold_vals, v[!is.na(v)])
      }
    }
    pooled_cache[[feat]][[fold_id]] <- fold_vals
    all_vals <- c(all_vals, fold_vals)
  }
  pooled_all[[feat]] <- all_vals
}

# Project onto common grid + mean-center
projected_list <- list()
feature_meta   <- list()

for (feat in cont_features) {
  dt_feat <- shape_all[feature == feat & !is.na(x)]
  if (nrow(dt_feat) < 10) next
  models_present <- unique(dt_feat$model)
  if (length(models_present) < 4) next

  pv <- pooled_all[[feat]]
  n_unique <- length(unique(pv))
  is_discrete <- n_unique < DISCRETE_THR

  x_grid <- seq(min(dt_feat$x), max(dt_feat$x), length.out = G_GRID)

  proj_dt <- rbindlist(lapply(models_present, function(m) {
    dt_m <- dt_feat[model == m]
    sfun <- make_shape_stepfun(dt_m)
    if (is.null(sfun)) return(data.table())
    info <- MODEL_STRUCTURE[model_tag == m]
    data.table(feature = feat, model = m, fold = info$fold, seed = info$seed,
               g = seq_len(G_GRID), x_g = x_grid, f_raw = sfun(x_grid))
  }), fill = TRUE)

  if (nrow(proj_dt) == 0) next
  proj_dt[, f_mean := mean(f_raw, na.rm = TRUE), by = model]
  proj_dt[, f_aligned := f_raw - f_mean]

  projected_list[[feat]] <- proj_dt
  feature_meta[[feat]] <- list(n_unique = n_unique, is_discrete = is_discrete)
}

cat(sprintf("  Projected %d features onto G=%d grid\n\n", length(projected_list), G_GRID))

# =============================================================================
# Step 2: Core QC metrics — C_core, J, S_seed
# =============================================================================

cat("[Step 2] Computing core QC metrics ...\n")

#' Compute fold-wise trim bounds using pooled training data
compute_trim_bounds <- function(feat) {
  fold_bounds <- lapply(names(FOLD_TRAIN_SITES), function(fold_id) {
    pv <- pooled_cache[[feat]][[fold_id]]
    if (length(pv) < 100) return(c(-Inf, Inf))
    as.numeric(quantile(pv, probs = c(TRIM_P_LO, TRIM_P_HI), na.rm = TRUE))
  })
  # Fold-wise intersection
  lo <- max(sapply(fold_bounds, `[`, 1), na.rm = TRUE)
  hi <- min(sapply(fold_bounds, `[`, 2), na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(list(trim_lo = -Inf, trim_hi = Inf, n_trim_grid = G_GRID))
  list(trim_lo = lo, trim_hi = hi,
       n_trim_grid = sum(projected_list[[feat]]$x_g >= lo & projected_list[[feat]]$x_g <= hi) / 6)
}

qc_list <- list()
for (feat in names(projected_list)) {
  proj <- projected_list[[feat]]
  meta <- feature_meta[[feat]]
  tb   <- compute_trim_bounds(feat)

  # Jaggedness J (per model, take median)
  j_per_model <- proj[, .(J = calc_jaggedness(f_aligned)), by = model]
  median_J    <- median(j_per_model$J, na.rm = TRUE)

  # Seed stability S_seed (within-fold median Spearman)
  seed_corrs <- numeric(0)
  for (fld in unique(proj$fold)) {
    f_s1 <- proj[fold == fld & seed == 1, f_aligned]
    f_s2 <- proj[fold == fld & seed == 2, f_aligned]
    if (length(f_s1) == length(f_s2) && length(f_s1) > 0)
      seed_corrs <- c(seed_corrs, safe_spearman(f_s1, f_s2))
  }
  S_seed <- median(seed_corrs, na.rm = TRUE)

  # Cross-fold concordance C_core (fold-averaged, trimmed)
  # Average across seeds within each fold → one curve per fold
  fold_avg <- proj[, .(f_fold = mean(f_aligned, na.rm = TRUE)), by = .(fold, g, x_g)]
  fold_ids <- sort(unique(fold_avg$fold))
  fold_pairs <- combn(fold_ids, 2, simplify = FALSE)

  c0_trim_vals   <- numeric(0)
  tied_rank_flag <- FALSE

  for (fp in fold_pairs) {
    f_k1_dt <- fold_avg[fold == fp[1]]
    f_k2_dt <- fold_avg[fold == fp[2]]
    if (nrow(f_k1_dt) != nrow(f_k2_dt)) next

    mask <- f_k1_dt$x_g >= tb$trim_lo & f_k1_dt$x_g <= tb$trim_hi
    if (sum(mask) >= MIN_GRID_FOR_CORR) {
      res <- safe_spearman_v2(f_k1_dt$f_fold[mask], f_k2_dt$f_fold[mask])
      if (isTRUE(res$tied_rank_flag)) tied_rank_flag <- TRUE
    } else {
      res <- safe_spearman_v2(f_k1_dt$f_fold, f_k2_dt$f_fold)
      if (isTRUE(res$tied_rank_flag)) tied_rank_flag <- TRUE
    }
    c0_trim_vals <- c(c0_trim_vals, res$rho)
  }

  C0_trim <- median(c0_trim_vals, na.rm = TRUE)

  qc_list[[feat]] <- data.table(
    feature = feat,
    C0_trim = round(C0_trim, 4),
    C0_min  = round(min(c0_trim_vals, na.rm = TRUE), 4),
    C0_max  = round(max(c0_trim_vals, na.rm = TRUE), 4),
    median_J = round(median_J, 3),
    S_seed   = round(S_seed, 4),
    is_discrete    = meta$is_discrete,
    tied_rank_flag = tied_rank_flag,
    trim_lo = tb$trim_lo,
    trim_hi = tb$trim_hi
  )
}

qc_dt <- rbindlist(qc_list)
cat(sprintf("  Computed QC for %d features\n", nrow(qc_dt)))

# =============================================================================
# Step 3: Bootstrap CI for C_core
# =============================================================================

cat(sprintf("\n[Step 3] Bootstrap CI for C_core (%d iterations) ...\n", BOOT_N))

boot_ci_list <- list()
for (feat in names(projected_list)) {
  proj <- projected_list[[feat]]
  tb   <- compute_trim_bounds(feat)

  fold_avg <- proj[, .(f_fold = mean(f_aligned, na.rm = TRUE)), by = .(fold, g, x_g)]
  fold_ids <- sort(unique(fold_avg$fold))
  fold_pairs <- combn(fold_ids, 2, simplify = FALSE)
  n_pairs <- length(fold_pairs)

  # Compute all fold-pair C0 values
  pair_c0 <- numeric(n_pairs)
  for (p in seq_len(n_pairs)) {
    fp <- fold_pairs[[p]]
    f_k1 <- fold_avg[fold == fp[1]]$f_fold
    f_k2 <- fold_avg[fold == fp[2]]$f_fold
    mask <- fold_avg[fold == fp[1]]$x_g >= tb$trim_lo & fold_avg[fold == fp[1]]$x_g <= tb$trim_hi
    if (sum(mask) >= MIN_GRID_FOR_CORR) {
      pair_c0[p] <- safe_spearman(f_k1[mask], f_k2[mask])
    } else {
      pair_c0[p] <- safe_spearman(f_k1, f_k2)
    }
  }

  # Bootstrap: resample fold-pairs with replacement
  set.seed(BOOTSTRAP_SEED)
  boot_medians <- replicate(BOOT_N, {
    idx <- sample(seq_len(n_pairs), n_pairs, replace = TRUE)
    median(pair_c0[idx], na.rm = TRUE)
  })

  boot_lo <- as.numeric(quantile(boot_medians, BOOT_ALPHA, na.rm = TRUE))
  boot_hi <- as.numeric(quantile(boot_medians, 1 - BOOT_ALPHA, na.rm = TRUE))
  crosses <- boot_lo < C0_THRESHOLD & boot_hi > C0_THRESHOLD

  boot_ci_list[[feat]] <- data.table(
    feature = feat, boot_lo = round(boot_lo, 4), boot_hi = round(boot_hi, 4),
    boot_se = round(sd(boot_medians, na.rm = TRUE), 4),
    crosses_threshold = crosses
  )
}
boot_ci_dt <- rbindlist(boot_ci_list)

qc_final <- merge(qc_dt, boot_ci_dt, by = "feature", all.x = TRUE)
cat(sprintf("  Bootstrap complete: %d features\n", nrow(boot_ci_dt)))

# =============================================================================
# Step 4: Traffic-light decision matrix
# =============================================================================

cat("\n[Step 4] Applying traffic-light decision matrix ...\n")

qc_final[, light := fifelse(
  tied_rank_flag == TRUE, "Gray",                              # Gray: rank-based QC lacks resolution
  fifelse(C0_trim >= C0_THRESHOLD & median_J <= J_THRESHOLD, "Green",
  fifelse(C0_trim >= C0_THRESHOLD & median_J >  J_THRESHOLD, "Yellow",
  fifelse(C0_trim <  C0_THRESHOLD, "Red", "Gray")))
)]

# Borderline flags
qc_final[, borderline := FALSE]
qc_final[light == "Red" & C0_trim >= (C0_THRESHOLD - 0.15) & C0_trim < C0_THRESHOLD,
         borderline := TRUE]
qc_final[light == "Green" & crosses_threshold == TRUE, borderline := TRUE]

# Action descriptions
qc_final[, action := fifelse(
  light == "Green" & borderline == FALSE,
  "Trust: shape is cross-site concordant and smooth. Use directly.",
  fifelse(light == "Green" & borderline == TRUE,
  "Trust with caution: concordant but CI crosses threshold. Report CI in supplement.",
  fifelse(light == "Yellow" & borderline == FALSE,
  "Refine: trend is reproducible but jagged. Apply post-hoc smoothing or reduce bins.",
  fifelse(light == "Yellow" & borderline == TRUE,
  "Refine with caution: jagged and CI crosses threshold. Smooth conservatively.",
  fifelse(light == "Red" & borderline == FALSE,
  "Site-specific: cross-site shapes diverge. Report per-site interpretations only.",
  fifelse(light == "Red" & borderline == TRUE,
  "Site-specific (borderline): near threshold but divergent. Investigate case-mix.",
  "Manual review: discrete/tied-rank feature. Inspect shape plot directly.")))))
)]

# Rank by importance order
qc_final[, final_rank := match(feature, FEATURE_RANK_ORDER)]
setorder(qc_final, final_rank)

cat("\n  ─── Table 4A: Feature-level transportability audit ───\n")
print(qc_final[, .(feature, light, C0_trim, boot_lo, boot_hi, median_J, S_seed, borderline)])

cat(sprintf("\n  Green=%d, Yellow=%d, Red=%d, Gray=%d\n",
            sum(qc_final$light == "Green"), sum(qc_final$light == "Yellow"),
            sum(qc_final$light == "Red"), sum(qc_final$light == "Gray")))

# =============================================================================
# Step 5: Visualizations
# =============================================================================

cat("\n[Step 5] Generating visualizations ...\n")
viz_dir <- file.path(QC_OUT, "Figures")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# 5a. Quadrant scatter: C_core vs J
quad_dt <- qc_final[!is.na(C0_trim) & !is.na(median_J) & is_discrete == FALSE]

if (nrow(quad_dt) > 0) {
  p_quad <- ggplot(quad_dt, aes(x = C0_trim, y = median_J)) +
    annotate("rect", xmin = C0_THRESHOLD, xmax = Inf, ymin = -Inf, ymax = J_THRESHOLD,
             fill = "#27ae60", alpha = 0.12) +
    annotate("rect", xmin = C0_THRESHOLD, xmax = Inf, ymin = J_THRESHOLD, ymax = Inf,
             fill = "#f1c40f", alpha = 0.12) +
    annotate("rect", xmin = -Inf, xmax = C0_THRESHOLD, ymin = -Inf, ymax = Inf,
             fill = "#e74c3c", alpha = 0.12) +
    geom_vline(xintercept = C0_THRESHOLD, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = J_THRESHOLD,  linetype = "dashed", color = "gray40") +
    geom_point(aes(color = light), size = 3)

  if (HAS_GGREPEL) {
    p_quad <- p_quad + ggrepel::geom_text_repel(aes(label = feature), size = 3, max.overlaps = 20)
  } else {
    warning("Package 'ggrepel' not installed; ShapeQC quadrant labels drawn without repel.")
    p_quad <- p_quad + geom_text(aes(label = feature), size = 3, vjust = -0.6, check_overlap = TRUE)
  }

  p_quad <- p_quad +
    scale_color_manual(values = TIER_COLORS) +
    labs(title = "ShapeQC v3.0: Feature-Level Governance Decision",
         subtitle = sprintf("Thresholds: C_core >= %.2f, J <= %.2f | Primary trim: %s",
                            C0_THRESHOLD, J_THRESHOLD, PRIMARY_TRIM),
         x = expression(C[core] ~ "(trimmed cross-fold concordance)"),
         y = expression(J ~ "(jaggedness index)"),
         color = "Tier") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(viz_dir, "QC_00_Quadrant_Decision.png"), p_quad,
         width = 10, height = 8, dpi = 300)
  cat("  Quadrant plot saved.\n")
}

# =============================================================================
# Step 6: Excel report + Supp sensitivity (T1-T5)
# =============================================================================

cat("\n[Step 6] Writing Excel report ...\n")

wb <- createWorkbook()
add_sheet <- function(wb, name, dt) {
  addWorksheet(wb, name)
  writeData(wb, name, dt)
}

# Decision matrix
decision_cols <- c("feature", "light", "borderline", "action",
                   "C0_trim", "boot_lo", "boot_hi", "median_J", "S_seed",
                   "is_discrete", "tied_rank_flag", "trim_lo", "trim_hi")
decision_dt <- qc_final[, .SD, .SDcols = intersect(decision_cols, names(qc_final))]
add_sheet(wb, "Decision_Matrix", decision_dt)

# Bootstrap CI
add_sheet(wb, "Bootstrap_CI", boot_ci_dt)

# Projected shapes
proj_all <- rbindlist(projected_list, fill = TRUE)
add_sheet(wb, "Projected_Shapes", proj_all[, .(feature, model, fold, seed, x_g, f_aligned)])

# Supp sensitivity: C_core across T1-T5 trim configs
cat("  Computing Supp_Sensitivity (T1-T5) ...\n")
sens_list <- list()
for (config_name in names(TRIM_CONFIGS)) {
  tc <- TRIM_CONFIGS[[config_name]]
  for (feat in names(projected_list)) {
    proj <- projected_list[[feat]]
    fold_avg <- proj[, .(f_fold = mean(f_aligned, na.rm = TRUE)), by = .(fold, g, x_g)]
    fold_ids <- sort(unique(fold_avg$fold))
    fold_pairs <- combn(fold_ids, 2, simplify = FALSE)

    # Compute trim bounds for this config
    fold_bounds <- lapply(names(FOLD_TRAIN_SITES), function(fid) {
      pv <- pooled_cache[[feat]][[fid]]
      if (length(pv) < 100) return(c(-Inf, Inf))
      as.numeric(quantile(pv, probs = tc, na.rm = TRUE))
    })
    lo_s <- max(sapply(fold_bounds, `[`, 1), na.rm = TRUE)
    hi_s <- min(sapply(fold_bounds, `[`, 2), na.rm = TRUE)
    if (!is.finite(lo_s) || !is.finite(hi_s)) { lo_s <- -Inf; hi_s <- Inf }

    c0_vals_s <- numeric(length(fold_pairs))
    for (pp in seq_along(fold_pairs)) {
      fp <- fold_pairs[[pp]]
      f_k1 <- fold_avg[fold == fp[1]]$f_fold
      f_k2 <- fold_avg[fold == fp[2]]$f_fold
      mask_s <- fold_avg[fold == fp[1]]$x_g >= lo_s & fold_avg[fold == fp[1]]$x_g <= hi_s
      if (sum(mask_s) >= MIN_GRID_FOR_CORR) {
        c0_vals_s[pp] <- safe_spearman(f_k1[mask_s], f_k2[mask_s])
      } else {
        c0_vals_s[pp] <- safe_spearman(f_k1, f_k2)
      }
    }

    sens_list[[length(sens_list) + 1]] <- data.table(
      feature = feat, config = config_name,
      p_lo = tc[1], p_hi = tc[2],
      C0_trim = round(median(c0_vals_s, na.rm = TRUE), 4)
    )
  }
}
sens_dt <- rbindlist(sens_list)
add_sheet(wb, "Supp_Sensitivity", sens_dt)

sens_wide <- dcast(sens_dt, feature ~ config, value.var = "C0_trim")
add_sheet(wb, "Supp_Sensitivity_Wide", sens_wide)

xlsx_path <- file.path(QC_OUT, "ShapeQC_v3.0_Report.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("  Excel saved: %s\n", basename(xlsx_path)))

# Also save CSVs
fwrite(decision_dt, file.path(QC_OUT, "Table_4A_Decision_Matrix.csv"))
fwrite(sens_wide,   file.path(QC_OUT, "Supp_S14_Trim_Sensitivity.csv"))

# Supp S3: Cross-model correlations (all 15 pairs)
cat("  Computing Supp S3 (cross-model correlations) ...\n")
s3_list <- list()
for (feat in names(projected_list)) {
  proj <- projected_list[[feat]]
  models <- unique(proj$model)
  if (length(models) < 2) next
  pairs <- combn(models, 2, simplify = FALSE)
  tb <- compute_trim_bounds(feat)

  for (pr in pairs) {
    f1 <- proj[model == pr[1] & x_g >= tb$trim_lo & x_g <= tb$trim_hi, f_aligned]
    f2 <- proj[model == pr[2] & x_g >= tb$trim_lo & x_g <= tb$trim_hi, f_aligned]
    if (length(f1) == length(f2) && length(f1) > 5) {
      s3_list[[length(s3_list) + 1]] <- data.table(
        feature = feat, pair = paste(pr, collapse = " vs "),
        pearson_r = round(cor(f1, f2, method = "pearson"), 4),
        spearman_rho = round(cor(f1, f2, method = "spearman"), 4)
      )
    }
  }
}
s3_dt <- rbindlist(s3_list)
fwrite(s3_dt, file.path(QC_OUT, "Supp_S3_Cross_Model_Correlations.csv"))

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 05 complete.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("  Output: %s\n", QC_OUT))
cat("  Files: ShapeQC_v3.0_Report.xlsx, Table_4A, Supp_S3, Supp_S14\n")

cat("\n  [Verify] Manuscript Table 4A key values:\n")
cat("    IDH_N_28D: C0=1.00, J=1.00 → Green\n")
cat("    Dialysate_Flow_Rate: C0=0.15, J=1.84 → Red\n")
cat("    Blood_Flow_Rate: C0=0.89, J=2.30 → Yellow\n")
for (feat_check in c("IDH_N_28D", "Dialysate_Flow_Rate", "Blood_Flow_Rate")) {
  row <- qc_final[feature == feat_check]
  if (nrow(row) > 0) {
    cat(sprintf("    Got %s: C0=%.2f, J=%.2f → %s\n",
                feat_check, row$C0_trim, row$median_J, row$light))
  }
}
cat("\n")
