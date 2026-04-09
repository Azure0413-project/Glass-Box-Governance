# =============================================================================
# 08_sensitivity_lofo.R — Smoothing Sensitivity + LOFO + Spar Robustness
# IDH EBM Governance Reproducibility Pipeline
#
# v3.0 update: External validation data now filtered to E-cohort only,
#   matching 03_iecv_ebm.R v3.0 pure E-only design.
#
# Outputs:
#   Supp Table S8  — Six-scenario sensitivity analysis
#   Supp Table S9  — Pairwise cross-scenario ΔAUROC comparisons
#   Supp Table S10 — Leave-one-feature-out marginal contributions
#   Supp Table S11 — Three-scenario spar robustness summary
#
# Source: consolidated from original analysis scripts
#
# Design:
#   Part A — 6 scenarios varying feature scope × spar formula:
#     S1: Yellow-only, J-adaptive
#     S2: Yellow+Red, J-adaptive (= v2.0 baseline; "Orange" in legacy 5-tier)
#     S3: Yellow+Red, J-adaptive (lower slope b=0.10)
#     S4: Yellow+Red, J-adaptive (higher slope b=0.20)
#     S5: Yellow+Red, Fixed spar=0.60
#     S6: Yellow+Red, Fixed spar=0.80
#   Part B — LOFO: leave-one-feature-out marginal contribution
#   Part C — 3-scenario spar (clinical framing: S2 vs S3 vs S4)
#
# Methods → "Alternative smoothing specifications, zeroing-versus-smoothing
#  comparisons, and leave-one-feature-out attribution all supported
#  the interpretation"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 03_iecv_ebm.R and 05_shapeqc.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(fst); library(reticulate); library(precrec)
  library(ggplot2); library(openxlsx)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")
if (!exists("init_python_full")) src("R/00_utils_python.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 08: Sensitivity + LOFO + Spar Robustness\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

ANALYSIS_TS <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
init_python_full()
joblib <- reticulate::import("joblib", convert = TRUE)

# --- Auto-detect ShapeQC and IECV directories ---
qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
qc_dirs <- qc_dirs[grepl("ShapeQC_v3\\.0_", basename(qc_dirs))]
if (length(qc_dirs) == 0) stop("ShapeQC not found. Run 05_shapeqc.R first.")
qc_dir <- sort(qc_dirs, decreasing = TRUE)[1]

xlsx_path <- file.path(qc_dir, "ShapeQC_v3.0_Report.xlsx")
csv_path  <- file.path(qc_dir, "Table_4A_Decision_Matrix.csv")
qc_decision <- if (file.exists(csv_path)) {
  fread(csv_path)
} else if (file.exists(xlsx_path)) {
  as.data.table(openxlsx::read.xlsx(xlsx_path, sheet = "Decision_Matrix"))
} else {
  NULL
}
if (!is.null(qc_decision) && "light" %in% names(qc_decision) && !("tier" %in% names(qc_decision))) {
  setnames(qc_decision, "light", "tier")
}

MASTER_OUT <- file.path(RUN_DIR, paste0("Sensitivity_LOFO_", ANALYSIS_TS))
dir.create(MASTER_OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Part A: 6-Scenario Definitions
# =============================================================================

# Features by original 5-tier (pre-simplified) for sensitivity scenario compatibility.
# In the final 4-tier manuscript, Orange was merged into Red and zeroed.
# These labels are used ONLY for the sensitivity smoothing scenarios (S2–S6)
# which test the "what if we had smoothed instead of zeroed?" counterfactual.

# --- Dynamic tier loading from 05_shapeqc.R (overrides hardcoded defaults) ---
YELLOW_ORIG    <- NULL
RED_SMOOTHABLE <- NULL
tryCatch({
  if (!is.null(qc_decision) && "tier" %in% names(qc_decision)) {
    YELLOW_ORIG    <- qc_decision[tier == "Yellow", feature]
    # Exclude Dialysate_Flow_Rate: C0 too low for meaningful smoothing comparison
    RED_SMOOTHABLE <- setdiff(qc_decision[tier == "Red", feature], "Dialysate_Flow_Rate")
    cat(sprintf("  [INFO] Sensitivity tiers loaded dynamically: Yellow=%d, Red_smoothable=%d\n",
                length(YELLOW_ORIG), length(RED_SMOOTHABLE)))
  }
}, error = function(e) {
  cat(sprintf("  [WARN] Dynamic tier load failed (%s); using hardcoded sensitivity tiers\n", e$message))
})

# Fallback to hardcoded if dynamic load failed
if (is.null(YELLOW_ORIG) || length(YELLOW_ORIG) == 0)
  YELLOW_ORIG <- c("Blood_Flow_Rate", "Heart_Rate", "Start_DBP")
if (is.null(RED_SMOOTHABLE) || length(RED_SMOOTHABLE) == 0)
  RED_SMOOTHABLE <- c("Target_UF_Volume", "Dry_Weight", "Age", "Pre_HD_Weight", "Dialysate_Temperature")
# Dialysate_Flow_Rate excluded: C0 too low for a meaningful smoothing counterfactual.
# Legacy alias for backward compatibility
ORANGE_ORIG <- RED_SMOOTHABLE

# ⚠ NOTE: S2–S4 include a +0.05 spar penalty for Red-smoothable ("orange") features
#   (spar_c = 0.05), inherited from the legacy 5-tier sensitivity design. This means
#   S2 does NOT exactly replicate the primary governance formula in 00_config.R
#   (SPAR_FORMULA: 0.30 + 0.15*J, no orange penalty) used by 07 and 04f.
#   S2 is a legacy sensitivity baseline; Table 4B results come from 07, not S2.
SCENARIOS <- list(
  list(id = "S1_YellowOnly", label = "S1: Yellow-only, J-adaptive",
       features = YELLOW_ORIG,
       spar_mode = "jadaptive", spar_a = 0.30, spar_b = 0.15, spar_c = 0,
       spar_lo = 0.50, spar_hi = 0.90, spar_fixed = NA),
  list(id = "S2_Full_v2.0", label = "S2: Yellow+Red, J-adaptive (baseline)",
       features = c(YELLOW_ORIG, ORANGE_ORIG),
       spar_mode = "jadaptive", spar_a = 0.30, spar_b = 0.15, spar_c = 0.05,
       spar_lo = 0.50, spar_hi = 0.90, spar_fixed = NA),
  list(id = "S3_LowSlope", label = "S3: Yellow+Red, J-adaptive (b=0.10)",
       features = c(YELLOW_ORIG, ORANGE_ORIG),
       spar_mode = "jadaptive", spar_a = 0.30, spar_b = 0.10, spar_c = 0.05,
       spar_lo = 0.50, spar_hi = 0.90, spar_fixed = NA),
  list(id = "S4_HighSlope", label = "S4: Yellow+Red, J-adaptive (b=0.20)",
       features = c(YELLOW_ORIG, ORANGE_ORIG),
       spar_mode = "jadaptive", spar_a = 0.30, spar_b = 0.20, spar_c = 0.05,
       spar_lo = 0.50, spar_hi = 0.90, spar_fixed = NA),
  list(id = "S5_Fixed060", label = "S5: Yellow+Red, Fixed spar=0.60",
       features = c(YELLOW_ORIG, ORANGE_ORIG),
       spar_mode = "fixed", spar_fixed = 0.60,
       spar_a = NA, spar_b = NA, spar_c = NA, spar_lo = NA, spar_hi = NA),
  list(id = "S6_Fixed080", label = "S6: Yellow+Red, Fixed spar=0.80",
       features = c(YELLOW_ORIG, ORANGE_ORIG),
       spar_mode = "fixed", spar_fixed = 0.80,
       spar_a = NA, spar_b = NA, spar_c = NA, spar_lo = NA, spar_hi = NA)
)

#' Compute spar for a given scenario + feature
compute_scenario_spar <- function(scenario, feature_name) {
  if (scenario$spar_mode == "fixed") return(scenario$spar_fixed)
  J <- NA_real_
  if (!is.null(qc_decision)) {
    rows <- qc_decision[feature == feature_name]
    if (nrow(rows) > 0) J <- rows$median_J[1]
  }
  if (is.na(J)) J <- 1.5  # default
  is_orange <- feature_name %in% ORANGE_ORIG
  base_spar <- scenario$spar_a + scenario$spar_b * J + scenario$spar_c * as.numeric(is_orange)
  min(scenario$spar_hi, max(scenario$spar_lo, base_spar))
}

cat(sprintf("  Scenarios: %d\n", length(SCENARIOS)))
cat(sprintf("  Output: %s\n\n", basename(MASTER_OUT)))

# =============================================================================
# Core smoothing function (shared by scenarios and LOFO)
# Applies spline smoothing to selected features, evaluates on external set.
# =============================================================================

run_smoothing_scenario <- function(sc_features, spar_map, iter_dir, scenario_label,
                                    collect_mp = TRUE) {
  info <- parse_iter_info(iter_dir)
  iter_tag <- info$iter_tag; ext_site <- info$external

  # Load original model
  orig_iter <- file.path(RUN_DIR, iter_tag)
  model_path <- find_latest_file(orig_iter, "_Final_EBM\\.joblib$")
  if (is.na(model_path)) return(NULL)

  ebm <- tryCatch(joblib$load(model_path), error = function(e) NULL)
  if (is.null(ebm)) return(NULL)

  thr <- NA_real_
  art_path <- find_latest_file(orig_iter, "_artifacts\\.rds$")
  if (!is.na(art_path)) {
    art <- tryCatch(readRDS(art_path), error = function(e) NULL)
    if (!is.null(art$threshold)) thr <- as.numeric(art$threshold)
  }

  # p_pre (original)
  if (is.null(site_cache[[ext_site]]))
    site_cache[[ext_site]] <<- load_site_eonly(ext_site, verbose = FALSE)
  dt_ext <- site_cache[[ext_site]]
  X_ext  <- get_predictor_df(dt_ext)
  y_ext  <- as.integer(dt_ext[[TARGET_COL]])
  p_pre  <- predict_proba_pos_r(ebm, X_ext)

  # Apply smoothing
  shift_total <- 0.0
  for (feat in sc_features) {
    spar_val <- spar_map[[feat]]
    if (is.na(spar_val)) next

    bins <- tryCatch(extract_bins_with_diagnostics(ebm, feat), error = function(e) list(success = FALSE))
    if (!bins$success) next

    dt_bins <- bins$dt[is.finite(x) & !is_missing & !is_unknown]
    if (nrow(dt_bins) < 4) next

    # Spline smoothing
    # NOTE: This simplified smoothing intentionally omits physio-range restriction
    # and bin_weight weighting (unlike the gold-standard smooth_scores_spline_fixed
    # in 04f). This is acceptable for cross-scenario sensitivity analysis (S1-S6)
    # where the goal is to compare spar formulas, not replicate exact governance.
    ss <- tryCatch(smooth.spline(dt_bins$x, dt_bins$score, spar = spar_val),
                   error = function(e) NULL)
    if (is.null(ss)) next

    new_scores <- as.numeric(bins$dt$score)  # preserve original scores for non-fitted bins
    smoothed_vals <- predict(ss, dt_bins$x)$y
    new_scores[dt_bins$bin_index + 1L] <- smoothed_vals

    # Preserve missing/unknown bins
    if (is.finite(bins$missing_idx)) new_scores[bins$missing_idx + 1L] <- bins$dt$score[bins$missing_idx + 1L]
    if (is.finite(bins$unknown_idx)) new_scores[bins$unknown_idx + 1L] <- bins$dt$score[bins$unknown_idx + 1L]

    delta <- tryCatch({
      d <- py_update_term(ebm, as.integer(bins$term_index),
                          reticulate::r_to_py(new_scores),
                          mean_preserve = TRUE, preserve_missing = TRUE,
                          missing_index = as.integer(bins$missing_idx),
                          shift_bagged_scores = TRUE, shift_bagged_intercept = TRUE)
      suppressWarnings(as.numeric(reticulate::py_to_r(d)))
    }, error = function(e) NA_real_)

    if (is.finite(delta)) shift_total <- shift_total + delta
  }

  # p_post
  p_post <- predict_proba_pos_r(ebm, X_ext)
  auc_pre  <- auc_pair(y_ext, p_pre)
  auc_post <- auc_pair(y_ext, p_post)

  data.table(
    iter_tag = iter_tag, external = ext_site, scenario_label = scenario_label,
    AUROC_pre = auc_pre["auroc"], AUROC_post = auc_post["auroc"],
    d_AUROC = auc_post["auroc"] - auc_pre["auroc"],
    AUPRC_pre = auc_pre["auprc"], AUPRC_post = auc_post["auprc"],
    d_AUPRC = auc_post["auprc"] - auc_pre["auprc"],
    Brier_pre = brier_score(y_ext, p_pre), Brier_post = brier_score(y_ext, p_post),
    d_Brier = brier_score(y_ext, p_post) - brier_score(y_ext, p_pre)
  )
}

# =============================================================================
# Execute scenarios
# =============================================================================

iter_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
iter_dirs <- iter_dirs[grepl("^Iter[0-9]+_External_", basename(iter_dirs))]
iter_dirs <- iter_dirs[!grepl("Posthoc|ShapeQC|Sensitivity|NRI|TV|temporal", basename(iter_dirs))]

# Lazy-loading site data cache (<<- used inside run_smoothing_scenario to
# avoid reloading per iteration; safe because scenarios run sequentially)
site_cache <- list()
MASTER_RESULTS <- data.table()

cat("[Part A] Running 6 scenarios ...\n")

for (sc_idx in seq_along(SCENARIOS)) {
  sc <- SCENARIOS[[sc_idx]]
  cat(sprintf("\n  --- Scenario %d/%d: %s ---\n", sc_idx, length(SCENARIOS), sc$label))

  spar_map <- setNames(
    sapply(sc$features, function(f) compute_scenario_spar(sc, f)),
    sc$features
  )

  for (iter_dir in iter_dirs) {
    row <- tryCatch(
      run_smoothing_scenario(sc$features, spar_map, iter_dir, sc$id),
      error = function(e) NULL
    )
    if (!is.null(row)) {
      row[, scenario_id := sc$id]
      MASTER_RESULTS <- rbind(MASTER_RESULTS, row, fill = TRUE)
    }
  }
}

# =============================================================================
# Supp S8: Cross-scenario comparison
# =============================================================================

cat("\n\n[Supp S8] Cross-scenario comparison ...\n")

cross_overall <- MASTER_RESULTS[, .(
  n_iter = .N,
  mean_d_AUROC = mean(d_AUROC, na.rm = TRUE), sd_d_AUROC = sd(d_AUROC, na.rm = TRUE),
  mean_d_AUPRC = mean(d_AUPRC, na.rm = TRUE), sd_d_AUPRC = sd(d_AUPRC, na.rm = TRUE),
  mean_d_Brier = mean(d_Brier, na.rm = TRUE)
), by = scenario_id]

fwrite(cross_overall, file.path(MASTER_OUT, "Supp_S8_Cross_Scenario.csv"))
fwrite(MASTER_RESULTS, file.path(MASTER_OUT, "Supp_S8_All_Iterations.csv"))

cat("\n  Supp Table S8:\n")
print(cross_overall)

# =============================================================================
# Supp S9: Pairwise scenario comparisons
# =============================================================================

cat("\n[Supp S9] Pairwise comparisons ...\n")

sc_ids <- unique(MASTER_RESULTS$scenario_id)
pairwise_list <- list()
for (i in seq_along(sc_ids)) {
  for (j in seq_along(sc_ids)) {
    if (i >= j) next
    d1 <- MASTER_RESULTS[scenario_id == sc_ids[i]]
    d2 <- MASTER_RESULTS[scenario_id == sc_ids[j]]
    merged <- merge(d1[, .(iter_tag, d_AUROC_1 = d_AUROC)],
                    d2[, .(iter_tag, d_AUROC_2 = d_AUROC)], by = "iter_tag")
    if (nrow(merged) >= 3) {
      tt <- tryCatch(t.test(merged$d_AUROC_1, merged$d_AUROC_2, paired = TRUE),
                     error = function(e) list(p.value = NA, estimate = NA))
      pairwise_list[[length(pairwise_list) + 1]] <- data.table(
        sc_A = sc_ids[i], sc_B = sc_ids[j],
        mean_diff = as.numeric(tt$estimate), p_value = tt$p.value, n = nrow(merged)
      )
    }
  }
}
pairwise_dt <- rbindlist(pairwise_list)
fwrite(pairwise_dt, file.path(MASTER_OUT, "Supp_S9_Pairwise.csv"))

# =============================================================================
# Part B: LOFO (Leave-One-Feature-Out)
# =============================================================================

cat("\n[Part B] Leave-One-Feature-Out ...\n")

TARGET_FEATURES_LOFO <- c(YELLOW_ORIG, ORANGE_ORIG)  # 8 features
LOFO_RESULTS <- data.table()

# Reference: S2 (full Yellow+Red)
s2_ref <- MASTER_RESULTS[scenario_id == "S2_Full_v2.0"]

for (fi in seq_along(TARGET_FEATURES_LOFO)) {
  leave_out <- TARGET_FEATURES_LOFO[fi]
  remaining <- setdiff(TARGET_FEATURES_LOFO, leave_out)

  cat(sprintf("  [LOFO %d/%d] Without: %-25s (keeping %d)\n",
              fi, length(TARGET_FEATURES_LOFO), leave_out, length(remaining)))

  # Use S2 spar settings for remaining features
  sc_s2 <- SCENARIOS[[2]]
  spar_map <- setNames(sapply(remaining, function(f) compute_scenario_spar(sc_s2, f)), remaining)

  for (iter_dir in iter_dirs) {
    row <- tryCatch(
      run_smoothing_scenario(remaining, spar_map, iter_dir, sprintf("LOFO_without_%s", leave_out)),
      error = function(e) NULL
    )
    if (!is.null(row)) {
      row[, left_out_feature := leave_out]
      LOFO_RESULTS <- rbind(LOFO_RESULTS, row, fill = TRUE)
    }
  }
}

fwrite(LOFO_RESULTS, file.path(MASTER_OUT, "LOFO_raw_results.csv"))

# Supp S10: Marginal contributions
cat("\n[Supp S10] Marginal contributions ...\n")

MARGINAL_TABLE <- data.table()
for (feat in TARGET_FEATURES_LOFO) {
  lofo_feat <- LOFO_RESULTS[left_out_feature == feat,
                              .(iter_tag, external, d_AUROC_lofo = d_AUROC, d_AUPRC_lofo = d_AUPRC)]
  merged <- merge(s2_ref[, .(iter_tag, external, d_AUROC_full = d_AUROC, d_AUPRC_full = d_AUPRC)],
                  lofo_feat, by = c("iter_tag", "external"))
  if (nrow(merged) == 0) next
  merged[, marginal_AUROC := d_AUROC_full - d_AUROC_lofo]
  merged[, marginal_AUPRC := d_AUPRC_full - d_AUPRC_lofo]

  MARGINAL_TABLE <- rbind(MARGINAL_TABLE, data.table(
    feature = feat,
    tier = if (feat %in% YELLOW_ORIG) "Yellow" else "Red",
    mean_marg_AUROC = mean(merged$marginal_AUROC, na.rm = TRUE),
    sd_marg_AUROC   = sd(merged$marginal_AUROC, na.rm = TRUE),
    mean_marg_AUPRC = mean(merged$marginal_AUPRC, na.rm = TRUE),
    n_iter = nrow(merged)
  ), fill = TRUE)
}
setorder(MARGINAL_TABLE, -mean_marg_AUROC)
fwrite(MARGINAL_TABLE, file.path(MASTER_OUT, "Supp_S10_LOFO_Marginal.csv"))

cat("\n  Supp Table S10:\n")
print(MARGINAL_TABLE)

# =============================================================================
# Part C: Supp S11 — 3-scenario spar robustness (S2 vs S3 vs S4)
# =============================================================================

cat("\n[Supp S11] 3-scenario spar robustness ...\n")

spar3 <- cross_overall[scenario_id %in% c("S2_Full_v2.0", "S3_LowSlope", "S4_HighSlope")]
fwrite(spar3, file.path(MASTER_OUT, "Supp_S11_Spar_Robustness.csv"))
cat("\n  Supp Table S11:\n")
print(spar3)

# =============================================================================
# Visualization
# =============================================================================

if (ENABLE_PLOTTING && nrow(MARGINAL_TABLE) > 0) {
  cat("\n  Generating LOFO forest plot ...\n")
  mt <- copy(MARGINAL_TABLE)
  mt[, feature := factor(feature, levels = mt[order(mean_marg_AUROC), feature])]

  p_lofo <- ggplot(mt, aes(x = mean_marg_AUROC, y = feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = mean_marg_AUROC - 1.96 * sd_marg_AUROC,
                       xmax = mean_marg_AUROC + 1.96 * sd_marg_AUROC),
                   height = 0.3, color = "gray40") +
    geom_point(aes(color = tier), size = 3) +
    scale_color_manual(values = c("Yellow" = "#f1c40f", "Red" = "#e74c3c")) +
    labs(title = "LOFO: Marginal ΔAUROC contribution per feature",
         x = "Marginal ΔAUROC (S2_full - S2_without)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(MASTER_OUT, "Fig_LOFO_Forest.png"), p_lofo, width = 10, height = 6, dpi = 300)
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 08 complete.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("  Output: %s\n", MASTER_OUT))
cat(sprintf("  Total scenarios: %d | LOFO features: %d\n",
            length(SCENARIOS), length(TARGET_FEATURES_LOFO)))
cat(sprintf("  Total result rows: %d (scenarios) + %d (LOFO)\n\n",
            nrow(MASTER_RESULTS), nrow(LOFO_RESULTS)))
