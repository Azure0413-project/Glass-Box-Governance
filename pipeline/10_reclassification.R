# =============================================================================
# 10_reclassification.R — NRI / IDI Reclassification Analysis
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Supp Table S13 Panel A — Per-iteration NRI/IDI with bootstrap CI
#   Supp Table S13 Panel B — Reclassification cross-tables (events/non-events)
#   Supp Table S13 Panel C — Subgroup reclassification by edited feature
#   Manuscript placeholder values for Results section
#
# Source: consolidated from original analysis scripts
#
# Design:
#   p_pre  = original EBM predictions (before any governance)
#   p_post = Yellow-smoothed + Red-zeroed predictions (after full governance)
#   Comparison at 5%/10% thresholds → 3 risk categories
#
# Methods → "Categorical NRI was computed using prespecified risk thresholds
#  of 5% and 10%, defining low-, intermediate-, and high-risk categories.
#  Continuous NRI and IDI were also calculated."
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 07_posthoc_governance.R (needs predictions_redtier.fst)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(fst)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt") || !exists("load_site_eonly", mode = "function")) src("R/00_utils_r.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 10: NRI / IDI Reclassification Analysis\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

ANALYSIS_TS <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

# --- Auto-detect governance output (v3.0 or legacy v2.2) ---
gov_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
gov_dirs <- gov_dirs[grepl("Posthoc_v3\\.0_Governance|Posthoc_v2\\.2_4Tier|Posthoc_v2\\.2_RedTier",
                           basename(gov_dirs))]
if (length(gov_dirs) == 0) stop("No governance output found. Run 07_posthoc_governance.R first.")
REDTIER_DIR <- sort(gov_dirs, decreasing = TRUE)[1]

OUT_ROOT <- file.path(RUN_DIR, paste0("NRI_Reclassification_v1.1_", ANALYSIS_TS))
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  Governance dir: %s\n", basename(REDTIER_DIR)))
cat(sprintf("  Output: %s\n", basename(OUT_ROOT)))
cat(sprintf("  Thresholds: %s\n", paste0(NRI_THRESHOLDS * 100, "%", collapse = " / ")))
cat("  p_pre = original EBM | p_post = smoothed + Red zeroed\n\n")

CURRENT_TIER_MAP <- FEATURE_TIER_MAP
tryCatch({
  qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
  qc_dirs <- sort(qc_dirs[grepl("ShapeQC", basename(qc_dirs))], decreasing = TRUE)
  if (length(qc_dirs) > 0) {
    dm_path <- file.path(qc_dirs[1], "Table_4A_Decision_Matrix.csv")
    if (file.exists(dm_path)) {
      dm <- fread(dm_path)
      if ("light" %in% names(dm) && !("tier" %in% names(dm))) setnames(dm, "light", "tier")
      if (all(c("feature", "tier") %in% names(dm))) {
        CURRENT_TIER_MAP <- setNames(dm$tier, dm$feature)
        cat(sprintf("  [INFO] Tier map loaded from %s\n", basename(qc_dirs[1])))
      }
    }
  }
}, error = function(e) {
  cat(sprintf("  [WARN] Tier map fallback to 00_config.R (%s)\n", e$message))
})

EDIT_FEATURES_CURRENT <- unique(c(
  ALL_EDIT_FEATURES,
  names(CURRENT_TIER_MAP)[CURRENT_TIER_MAP %in% c("Yellow", "Red")]
))

# =============================================================================
# Phase 1: Load predictions
# Strategy: 07's output (v3.0 or v2.2) contains both p_pre and p_post.
# =============================================================================

cat("[Phase 1] Loading predictions ...\n")

iter_dirs_v22 <- list.dirs(REDTIER_DIR, recursive = FALSE, full.names = TRUE)
iter_dirs_v22 <- iter_dirs_v22[grepl("^Iter[0-9]+_External_", basename(iter_dirs_v22))]

all_data <- list()

for (iter_dir_v22 in iter_dirs_v22) {
  tag  <- basename(iter_dir_v22)
  info <- parse_iter_info(iter_dir_v22)

  # Primary: load from v2.2 (07's output has both p_pre and p_post)
  fst_v22 <- list.files(iter_dir_v22, pattern = "_predictions_redtier\\.fst$|_4tier_predictions\\.fst$",
                         full.names = TRUE)
  if (length(fst_v22) == 0) { cat(sprintf("  [Skip] No predictions: %s\n", tag)); next }
  dt_v22 <- read_fst_dt(fst_v22[1])

  # Normalize column names
  if ("y_true" %in% names(dt_v22)) setnames(dt_v22, "y_true", "outcome")
  if ("p_old" %in% names(dt_v22)) setnames(dt_v22, "p_old", "p_pre")
  if ("p_post_redtier" %in% names(dt_v22)) setnames(dt_v22, "p_post_redtier", "p_post")
  if ("patient_id" %in% names(dt_v22) && !(ID_COL %in% names(dt_v22)))
    setnames(dt_v22, "patient_id", ID_COL)

  if (!("p_pre" %in% names(dt_v22)) || !("p_post" %in% names(dt_v22))) {
    cat(sprintf("  [Skip] Missing p_pre/p_post: %s\n", tag)); next
  }

  dt_pred <- dt_v22
  dt_pred[, iter_tag := tag]
  dt_pred[, held_out_center := info$external]
  dt_pred[, seed := info$seed_id]

  # Rename cluster column
  if (CLUSTER_COL %in% names(dt_pred)) setnames(dt_pred, CLUSTER_COL, "cluster_id")
  if ("Patient_Cluster_ID" %in% names(dt_pred) && !("cluster_id" %in% names(dt_pred)))
    setnames(dt_pred, "Patient_Cluster_ID", "cluster_id")
  if ("cluster_id" %in% names(dt_v22) && !("cluster_id" %in% names(dt_pred)))
    setnames(dt_pred, "cluster_id", "cluster_id")  # already correct

  # Attach feature values for subgroup analysis.
  # Preferred path: use features already saved by 07_posthoc_governance.R.
  # Fallback path: align to E-only site data by row order only when needed.
  ext_site <- info$external
  missing_feat_cols <- setdiff(EDIT_FEATURES_CURRENT, names(dt_pred))
  if (length(missing_feat_cols) > 0) {
    dt_site <- load_site_eonly(ext_site, verbose = FALSE)
    feat_cols <- intersect(missing_feat_cols, names(dt_site))
    if (length(feat_cols) > 0) {
      # Positional alignment: both dt_pred and dt_site originate from
      # load_site_eonly() with identical row order.  No session-level unique
      # key exists (Patient_ID / Patient_Cluster_ID are patient-level),
      # so keyed merge would produce many-to-many inflation.
      if (nrow(dt_pred) == nrow(dt_site)) {
        for (fc in feat_cols) dt_pred[, (fc) := dt_site[[fc]]]
        cat(sprintf("  [%s] subgroup features attached by positional alignment\n", tag))
      } else {
        cat(sprintf("  [WARN] %s: row mismatch (pred=%d, site=%d); subgroup features skipped\n",
                    tag, nrow(dt_pred), nrow(dt_site)))
      }
    }
  }

  all_data[[tag]] <- dt_pred
  cat(sprintf("  %s: %s sessions | p_pre/p_post aligned\n",
              tag, format(nrow(dt_pred), big.mark = ",")))
}

if (length(all_data) == 0) stop("No usable governance prediction files found for NRI/IDI analysis.")
dat <- rbindlist(all_data, fill = TRUE)
cat(sprintf("\n  Total: %s sessions across %d iterations\n\n",
            format(nrow(dat), big.mark = ","), length(all_data)))

# =============================================================================
# Phase 2: Compute NRI / IDI per iteration
# =============================================================================

cat("[Phase 2] Computing NRI / IDI ...\n")

iterations  <- unique(dat$iter_tag)
all_results <- list()

for (it in iterations) {
  dt_iter <- dat[iter_tag == it]
  center  <- unique(dt_iter$held_out_center)

  cat(sprintf("\n  --- %s (Held-out: %s) ---\n", it, center))
  cat(sprintf("  N = %s | Events = %d (%.1f%%)\n",
              format(nrow(dt_iter), big.mark = ","),
              dt_iter[, sum(outcome)], dt_iter[, 100 * mean(outcome)]))

  cat_nri  <- compute_categorical_nri(dt_iter)
  cont_nri <- compute_continuous_nri(dt_iter)
  idi      <- compute_idi(dt_iter)

  cat(sprintf("  Cat NRI: %+.4f (ev: %+.4f, ne: %+.4f)\n",
              cat_nri$NRI, cat_nri$NRI_events, cat_nri$NRI_nonevents))
  cat(sprintf("  Cont NRI: %+.4f | IDI: %+.6f\n", cont_nri$cNRI, idi$IDI))
  cat(sprintf("  Reclassified: %d / %s (%.2f%%)\n",
              cat_nri$n_moved, format(cat_nri$n_total, big.mark = ","), cat_nri$pct_moved))

  cat("  Running bootstrap ...\n")
  boot <- bootstrap_nri_idi(dt_iter, B = B_BOOTSTRAP)

  all_results[[it]] <- list(
    iter_tag = it, held_out_center = center, n = nrow(dt_iter),
    n_events = dt_iter[, sum(outcome)],
    cat_nri = cat_nri, cont_nri = cont_nri, idi = idi, bootstrap = boot
  )
}

# =============================================================================
# Phase 3: Subgroup reclassification by edited feature
# =============================================================================

cat("\n\n[Phase 3] Subgroup reclassification ...\n")

subgroup_results <- list()
for (feat in EDIT_FEATURES_CURRENT) {
  tier <- CURRENT_TIER_MAP[feat]
  if (length(tier) == 0 || is.na(tier)) tier <- FEATURE_TIER_MAP[feat]
  if (!(feat %in% names(dat)) || all(is.na(dat[[feat]]))) next
  cat(sprintf("  %s (%s) ... ", feat, tier))

  per_iter <- list()
  for (it in iterations) {
    dt_iter   <- dat[iter_tag == it]
    feat_vals <- dt_iter[[feat]]
    if (mean(is.na(feat_vals)) > 0.5) next

    valid_idx <- !is.na(feat_vals)
    q_breaks  <- unique(quantile(feat_vals[valid_idx], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
    if (length(q_breaks) < 3) next

    dt_valid <- dt_iter[valid_idx]
    dt_valid[, feat_q := cut(get(feat), breaks = q_breaks, include.lowest = TRUE, labels = FALSE)]
    dt_valid[, delta_p := p_post - p_pre]
    dt_valid[, cat_pre := assign_risk_cat(p_pre)]
    dt_valid[, cat_post := assign_risk_cat(p_post)]
    dt_valid[, moved := as.integer(cat_pre != cat_post)]

    q_summ <- dt_valid[, .(
      n = .N, n_events = sum(outcome),
      mean_delta_p = mean(delta_p), median_delta_p = median(delta_p),
      mean_abs_delta_p = mean(abs(delta_p)),
      pct_moved = 100 * mean(moved)
    ), by = feat_q][order(feat_q)]

    per_iter[[it]] <- list(
      iter_tag = it, held_out = unique(dt_iter$held_out_center),
      pct_moved = dt_valid[, 100 * mean(moved)],
      mean_abs_delta_p = dt_valid[, mean(abs(delta_p))],
      quartile_summary = q_summ
    )
  }

  if (length(per_iter) == 0) { cat("skipped\n"); next }
  pct_vals <- sapply(per_iter, `[[`, "pct_moved")
  abs_vals <- sapply(per_iter, `[[`, "mean_abs_delta_p")
  subgroup_results[[feat]] <- list(
    feature = feat, tier = tier, per_iteration = per_iter,
    summary = list(mean_pct_moved = mean(pct_vals), range_pct_moved = range(pct_vals),
                   mean_abs_delta_p = mean(abs_vals))
  )
  cat(sprintf("reclass=%.2f%% | |Δp|=%.5f\n", mean(pct_vals), mean(abs_vals)))
}

# =============================================================================
# Phase 4: Output — Supp Table S13
# =============================================================================

cat("\n\n[Phase 4] Generating Supp Table S13 ...\n")

fmt_ci <- function(point, ci) sprintf("%.4f [%.4f, %.4f]", point, ci[1], ci[2])

# Panel A: Per-iteration NRI/IDI
panel_a <- rbindlist(lapply(all_results, function(res) {
  bc <- res$bootstrap$CI
  data.table(
    Held_out = res$held_out_center, N = res$n, Events = res$n_events,
    Pct_Reclass = sprintf("%.2f", res$cat_nri$pct_moved),
    Cat_NRI = fmt_ci(res$cat_nri$NRI, bc$cat_NRI),
    NRI_plus = sprintf("%.4f", res$cat_nri$NRI_events),
    NRI_minus = sprintf("%.4f", res$cat_nri$NRI_nonevents),
    Cont_NRI = fmt_ci(res$cont_nri$cNRI, bc$cNRI),
    IDI = fmt_ci(res$idi$IDI, bc$IDI),
    Rel_IDI = sprintf("%.4f", res$idi$relative_IDI)
  )
}))
fwrite(panel_a, file.path(OUT_ROOT, "Table_S13_Panel_A.csv"))
cat("  Saved: Table_S13_Panel_A.csv\n")

# Panel B: Cross-tables
panel_b_list <- list()
for (it in iterations) {
  dt_iter <- dat[iter_tag == it]
  dt_iter[, cat_pre := assign_risk_cat(p_pre)]
  dt_iter[, cat_post := assign_risk_cat(p_post)]
  for (oc in c(0L, 1L)) {
    ct <- dt_iter[outcome == oc, .N, by = .(cat_pre, cat_post)]
    ct_wide <- dcast(ct, cat_pre ~ cat_post, value.var = "N", fill = 0L)
    ct_wide[, iter_tag := it]
    ct_wide[, outcome := oc]
    ct_wide[, held_out := unique(dt_iter$held_out_center)]
    panel_b_list <- c(panel_b_list, list(ct_wide))
  }
}
panel_b <- rbindlist(panel_b_list, fill = TRUE)
fwrite(panel_b, file.path(OUT_ROOT, "Table_S13_Panel_B.csv"))
cat("  Saved: Table_S13_Panel_B.csv\n")

# Panel C: Subgroup summary
if (length(subgroup_results) > 0) {
  panel_c <- rbindlist(lapply(names(subgroup_results), function(fn) {
    sr <- subgroup_results[[fn]]
    data.table(Feature = fn, Tier = sr$tier,
               Mean_AbsDp = sprintf("%.5f", sr$summary$mean_abs_delta_p),
               Pct_Reclass = sprintf("%.2f", sr$summary$mean_pct_moved),
               Range = sprintf("%.2f-%.2f", sr$summary$range_pct_moved[1], sr$summary$range_pct_moved[2]))
  }))
  fwrite(panel_c, file.path(OUT_ROOT, "Table_S13_Panel_C.csv"))
  cat("  Saved: Table_S13_Panel_C.csv\n")

  # Quartile detail
  q_list <- list()
  for (fn in names(subgroup_results)) {
    sr <- subgroup_results[[fn]]
    for (itn in names(sr$per_iteration)) {
      pi <- sr$per_iteration[[itn]]
      if (!is.null(pi$quartile_summary)) {
        qs <- copy(pi$quartile_summary)
        qs[, feature := fn]; qs[, tier := sr$tier]; qs[, iter_tag := itn]
        q_list <- c(q_list, list(qs))
      }
    }
  }
  if (length(q_list) > 0) {
    fwrite(rbindlist(q_list, fill = TRUE), file.path(OUT_ROOT, "Subgroup_quartile_detail.csv"))
    cat("  Saved: Subgroup_quartile_detail.csv\n")
  }
}

# Complete results
saveRDS(list(nri_idi = all_results, subgroups = subgroup_results,
             config = list(thresholds = NRI_THRESHOLDS, B = B_BOOTSTRAP)),
        file.path(OUT_ROOT, "NRI_Complete_Results.rds"))

# =============================================================================
# Phase 5: Manuscript placeholder values
# =============================================================================

cat("\n\n  ─── Manuscript Values ───\n")

pct_moved     <- sapply(all_results, function(x) x$cat_nri$pct_moved)
cat_nri_vals  <- sapply(all_results, function(x) x$cat_nri$NRI)
cat_nri_ne    <- sapply(all_results, function(x) x$cat_nri$NRI_nonevents)
cont_nri_vals <- sapply(all_results, function(x) x$cont_nri$cNRI)

cat(sprintf("  Reclassified: %.2f%% (range %.2f–%.2f%%)\n",
            mean(pct_moved), min(pct_moved), max(pct_moved)))
cat(sprintf("  Categorical NRI: %+.4f (mean across 6 iters)\n", mean(cat_nri_vals)))
cat(sprintf("  NRI- (nonevents): %+.4f\n", mean(cat_nri_ne)))
cat(sprintf("  Continuous NRI: %+.4f\n", mean(cont_nri_vals)))

cat("\n  [Verify] Supp Table S13 expects:\n")
cat("    7.77% crossed thresholds (range 5.88–11.85%)\n")
cat("    Cat NRI: +0.0083 | NRI-: +0.0177 | Cont NRI: -0.2173\n")
cat("    NOTE: If the manuscript cites different numbers (e.g. 2.87%%),\n")
cat("          check consistency between main text and S13 values.
")

cat(sprintf("\n  Output: %s\n\n", OUT_ROOT))
