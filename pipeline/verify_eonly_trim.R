# =============================================================================
# verify_eonly_trim_alignment.R — v3.0 -> v3.0.2 Trim-Bound Sanity Check
# IDH EBM Governance Reproducibility Pipeline
#
# Purpose:
#   Before switching to E-only rerun of 05_shapeqc.R (v3.0.2), quantify the
#   impact of changing trim bounds from E+L to E-only on each predictor's
#   P0.5 / P99.5 quantiles, and predict tier-change risk.
#
# Output:
#   Trim_Bounds_Comparison.csv  — per feature x fold comparison of E+L vs E-only
#   Tier_Risk_Report.csv        — estimated tier-change risk per feature
#   Console summary
#
# Usage:
#   source("00_config.R")
#   source("00_utils_r.R")
#   source("verify_eonly_trim_alignment.R")
#
# Design rationale:
#   Trim bounds are computed from pooled training sites' P0.5 / P99.5
#   with fold-wise intersection ([max of low, min of high]).
#   Removing the 25% L-cohort samples should have minimal impact on
#   extreme quantiles, but per-feature confirmation is needed to check
#   whether any feature crosses a tier boundary.
#
# Prerequisites:
#   - 00_config.R (SITE_FILES, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
#   - 00_utils_r.R (read_fst_dt, prep_site_dt, apply_el_split, load_site_eonly)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
})

if (!exists("SITE_FILES"))      source("00_config.R")
if (!exists("load_site_eonly")) source("00_utils_r.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  verify_eonly_trim_alignment — E+L vs E-only trim bounds\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ---- Configuration (mirrors 05_shapeqc.R) ----
TRIM_P_LO <- 0.005
TRIM_P_HI <- 0.995
FOLD_TRAIN_SITES <- list(
  "1" = c("D6", "CY"),
  "2" = c("TN", "CY"),
  "3" = c("TN", "D6")
)
FEATURES_TO_CHECK <- c(
  "IDH_N_28D", "Pre_HD_SBP", "UF_BW_Perc", "Target_UF_Volume",
  "Blood_Flow_Rate", "IDH_N_7D", "Heart_Rate", "Dry_Weight",
  "Age", "Pre_HD_Weight", "Start_DBP", "Body_Temperature",
  "Dialysate_Temperature", "Dialysate_Flow_Rate", "Respiratory_Rate"
)

OUT_DIR <- file.path(RUN_DIR,
                     paste0("Verify_Eonly_Trim_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Output: %s\n\n", OUT_DIR))

# =============================================================================
# Step 1: Load both E+L and E-only site data
# =============================================================================

cat("[Step 1] Loading site data under both frames ...\n")

site_EL    <- list()   # full cohort (old frame)
site_Eonly <- list()   # E-cohort only (new frame)

for (s in names(SITE_FILES)) {
  fp <- SITE_FILES[[s]]
  if (!file.exists(fp)) {
    cat(sprintf("  [skip] %s: file not found\n", s)); next
  }
  dt_full <- prep_site_dt(read_fst_dt(fp), s)
  dt_full <- apply_el_split(dt_full, SESSION_DATE_COL, E_SPLIT_QUANTILE, COHORT_COL)
  site_EL[[s]]    <- dt_full
  site_Eonly[[s]] <- dt_full[get(COHORT_COL) == "E"]
  cat(sprintf("  %s: %s E+L rows | %s E-only rows (%.1f%%)\n",
              s,
              format(nrow(site_EL[[s]]),    big.mark = ","),
              format(nrow(site_Eonly[[s]]), big.mark = ","),
              100 * nrow(site_Eonly[[s]]) / nrow(site_EL[[s]])))
}

# =============================================================================
# Step 2: Compute trim bounds under both frames (feature x fold)
# =============================================================================

cat("\n[Step 2] Computing trim bounds (P0.5, P99.5) per feature x fold ...\n")

#' Pooled quantiles for a (feature, fold) given a site-data map.
pooled_q <- function(site_map, feat, fold_id) {
  vals <- numeric(0)
  for (s in FOLD_TRAIN_SITES[[fold_id]]) {
    dt <- site_map[[s]]
    if (is.null(dt) || !(feat %in% names(dt))) next
    v <- as.numeric(dt[[feat]])
    vals <- c(vals, v[!is.na(v)])
  }
  if (length(vals) < 100) return(list(lo = NA_real_, hi = NA_real_, n = length(vals)))
  list(lo = as.numeric(quantile(vals, TRIM_P_LO, na.rm = TRUE)),
       hi = as.numeric(quantile(vals, TRIM_P_HI, na.rm = TRUE)),
       n  = length(vals))
}

rows <- list()
for (feat in FEATURES_TO_CHECK) {
  for (fold_id in names(FOLD_TRAIN_SITES)) {
    q_el <- pooled_q(site_EL,    feat, fold_id)
    q_e  <- pooled_q(site_Eonly, feat, fold_id)
    rows[[length(rows) + 1]] <- data.table(
      feature   = feat,
      fold      = fold_id,
      train_sites = paste(FOLD_TRAIN_SITES[[fold_id]], collapse = "+"),
      n_EL      = q_el$n,
      n_Eonly   = q_e$n,
      P0_5_EL   = q_el$lo,
      P0_5_Eonly= q_e$lo,
      delta_lo  = q_e$lo - q_el$lo,
      P99_5_EL  = q_el$hi,
      P99_5_Eonly = q_e$hi,
      delta_hi  = q_e$hi - q_el$hi
    )
  }
}
bounds_dt <- rbindlist(rows)

# Fold-wise intersection (what 05_shapeqc.R actually uses for each feature)
intersect_dt <- bounds_dt[, .(
  trim_lo_EL    = max(P0_5_EL,    na.rm = TRUE),
  trim_hi_EL    = min(P99_5_EL,   na.rm = TRUE),
  trim_lo_Eonly = max(P0_5_Eonly, na.rm = TRUE),
  trim_hi_Eonly = min(P99_5_Eonly,na.rm = TRUE)
), by = feature]
intersect_dt[, delta_lo := trim_lo_Eonly - trim_lo_EL]
intersect_dt[, delta_hi := trim_hi_Eonly - trim_hi_EL]
intersect_dt[, range_EL    := trim_hi_EL    - trim_lo_EL]
intersect_dt[, range_Eonly := trim_hi_Eonly - trim_lo_Eonly]
intersect_dt[, pct_range_shrink := 100 * (range_EL - range_Eonly) / range_EL]

fwrite(bounds_dt,    file.path(OUT_DIR, "Trim_Bounds_Per_Fold.csv"))
fwrite(intersect_dt, file.path(OUT_DIR, "Trim_Bounds_Intersection.csv"))

cat("  Saved: Trim_Bounds_Per_Fold.csv (feature x fold detail)\n")
cat("  Saved: Trim_Bounds_Intersection.csv (fold-wise intersection)\n")

# =============================================================================
# Step 3: Flag features with non-trivial trim changes
# =============================================================================

cat("\n[Step 3] Flagging features with non-trivial changes ...\n")

# Threshold: flag if trim bound shifts by > 1% of feature range,
# or if the intersection range shrinks by > 2%.
THRESHOLD_BOUND_SHIFT_PCT <- 1.0   # % of range
THRESHOLD_RANGE_SHRINK_PCT <- 2.0  # %

intersect_dt[, lo_shift_pct := 100 * abs(delta_lo) / range_EL]
intersect_dt[, hi_shift_pct := 100 * abs(delta_hi) / range_EL]
intersect_dt[, any_flag := lo_shift_pct > THRESHOLD_BOUND_SHIFT_PCT |
                         hi_shift_pct > THRESHOLD_BOUND_SHIFT_PCT |
                         pct_range_shrink > THRESHOLD_RANGE_SHRINK_PCT]

flagged <- intersect_dt[any_flag == TRUE]

if (nrow(flagged) == 0) {
  cat("  No feature exceeds shift thresholds (P0.5/P99.5 move <1% of range,\n")
  cat("    intersection range shrinks <2%). Tier assignments should be stable.\n")
} else {
  cat(sprintf("  %d feature(s) flagged -- inspect below:\n\n", nrow(flagged)))
  print(flagged[, .(feature,
                    trim_lo_EL, trim_lo_Eonly, lo_shift_pct,
                    trim_hi_EL, trim_hi_Eonly, hi_shift_pct,
                    pct_range_shrink)])
}

# =============================================================================
# Step 4: Compare against existing ShapeQC output (if available)
# =============================================================================

cat("\n[Step 4] Comparing against existing ShapeQC output (if any) ...\n")

qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
qc_dirs <- qc_dirs[grepl("ShapeQC_v3\\.0(r|_Eonly)?_", basename(qc_dirs))]

if (length(qc_dirs) == 0) {
  cat("  No prior ShapeQC output found. Run 05_shapeqc.R at least once.\n")
} else {
  latest <- sort(qc_dirs, decreasing = TRUE)[1]
  dm_path <- file.path(latest, "Table_4A_Decision_Matrix.csv")
  if (file.exists(dm_path)) {
    dm <- fread(dm_path)
    if ("light" %in% names(dm) && !("tier" %in% names(dm))) setnames(dm, "light", "tier")
    cat(sprintf("  Latest ShapeQC: %s\n", basename(latest)))

    # Identify tier-boundary-sensitive features:
    #   C0_trim within +/-0.05 of 0.70, OR median_J within +/-0.10 of 1.50.
    at_risk <- dm[!is.na(C0_trim) & !is.na(median_J) &
                  (abs(C0_trim  - 0.70) < 0.05 |
                   abs(median_J - 1.50) < 0.10)]
    if (nrow(at_risk) > 0) {
      cat("\n  Boundary-sensitive features (existing assignment may flip under E-only):\n")
      print(at_risk[, .(feature, tier, C0_trim, median_J, trim_lo, trim_hi)])
      cat("\n  -> Prioritize re-inspection of these after 05_shapeqc.R v3.0.2 rerun.\n")
    } else {
      cat("  No feature is near any tier threshold (+/-0.05 on C0, +/-0.10 on J).\n")
      cat("    Tier assignments almost certainly stable under E-only.\n")
    }
  } else {
    cat(sprintf("  [Note] Decision matrix not found in %s\n", basename(latest)))
  }
}

# =============================================================================
# Step 5: Quick summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Summary\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat(sprintf("  Features evaluated : %d\n", length(FEATURES_TO_CHECK)))
cat(sprintf("  Flagged features   : %d\n", nrow(flagged)))
cat(sprintf("  Median |delta P0.5|    : %.4g (feature-scale units)\n",
            median(abs(intersect_dt$delta_lo), na.rm = TRUE)))
cat(sprintf("  Median |delta P99.5|   : %.4g\n",
            median(abs(intersect_dt$delta_hi), na.rm = TRUE)))
cat(sprintf("  Max range shrink   : %.2f%%\n",
            max(intersect_dt$pct_range_shrink, na.rm = TRUE)))

cat("\n  Interpretation guidance:\n")
cat("    - Expected: |delta quantile| dominated by sample-size noise; <1% of range.\n")
cat("    - If a boundary-sensitive feature shifts >1%, re-check its tier after rerun.\n")
cat("    - Gray-tier (tied-rank) features are unaffected (trim bounds irrelevant).\n")
cat(sprintf("\n  Output directory: %s\n\n", OUT_DIR))
