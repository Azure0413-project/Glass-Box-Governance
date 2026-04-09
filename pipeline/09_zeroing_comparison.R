# =============================================================================
# 09_zeroing_comparison.R — Zeroing-vs-Smoothing & τ Threshold Sensitivity
# IDH EBM Governance Reproducibility Pipeline
#
# v3.0 update: External validation data now filtered to E-cohort only,
#   matching 03_iecv_ebm.R v3.0 pure E-only design.
#
# Outputs:
#   Supp Table S12 — Tier-specific exclusion (zeroing out) vs smoothing
#   Supp Table S15 — Threshold sensitivity (τ=0.50 vs τ=0.70)
#
# Source: consolidated from original analysis scripts
#   S12: Post-hoc v2.3 (Red+Yellow Zeroing Out)
#   S15: Post-hoc v2.4 (Maximum Parsimony Zeroing Out)
#
# Design:
#   Scenario A (S12): Zero out all 9 non-Green features (6 Red + 3 Yellow)
#     → Compare with v2.2 (which zeros only 6 Red) to assess
#       "zeroing-versus-smoothing" for Yellow features
#   Scenario B (S15): Zero out all 10 non-Green features (6R + 3Y + 1 Gray)
#     → "Maximum parsimony" — what if we used τ=0.50 instead of τ=0.70?
#     → Methods: "|ΔAUROC difference| = 0.0007"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 03_iecv_ebm.R (original models)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(fst); library(reticulate); library(precrec)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")
if (!exists("init_python_full")) src("R/00_utils_python.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 09: Zeroing Comparison & τ Threshold Sensitivity\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

ANALYSIS_TS <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
init_python_full()
joblib <- reticulate::import("joblib", convert = TRUE)

# =============================================================================
# Feature lists for two scenarios
# =============================================================================

# S12: Zero Red + Yellow (9 features)
ZERO_RED_YELLOW <- c(
  # Red (6)
  "Dialysate_Flow_Rate", "Dry_Weight", "Pre_HD_Weight",
  "Age", "Target_UF_Volume", "Dialysate_Temperature",
  # Yellow (3)
  "Blood_Flow_Rate", "Heart_Rate", "Start_DBP"
)

# S15: Maximum parsimony — zero all non-Green (10 features)
ZERO_MAX_PARSIMONY <- c(
  ZERO_RED_YELLOW,
  "IDH_N_7D"  # Gray
)

SCENARIOS_ZERO <- list(
  list(id = "S12_RedYellow", label = "S12: Zero Red+Yellow (9 features)",
       features = ZERO_RED_YELLOW, supp_table = "S12"),
  list(id = "S15_MaxParsimony", label = "S15: Maximum Parsimony (10 features)",
       features = ZERO_MAX_PARSIMONY, supp_table = "S15")
)

# =============================================================================
# Core zeroing function (shared across scenarios)
# =============================================================================

run_zeroing_scenario <- function(zero_features, iter_dir, scenario_label) {
  info     <- parse_iter_info(iter_dir)
  iter_tag <- info$iter_tag
  ext_site <- info$external

  # Load ORIGINAL model (not v2.0 patched — zeroing from scratch)
  orig_iter  <- file.path(RUN_DIR, iter_tag)
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

  # External data
  if (is.null(site_cache[[ext_site]]))
    site_cache[[ext_site]] <<- load_site_eonly(ext_site, verbose = FALSE)
  dt_ext <- site_cache[[ext_site]]
  X_ext  <- get_predictor_df(dt_ext)
  y_ext  <- as.integer(dt_ext[[TARGET_COL]])

  # p_pre (original)
  p_pre <- predict_proba_pos_r(ebm, X_ext)

  # Zero out features
  shift_total <- 0.0
  n_zeroed    <- 0L
  for (feat in zero_features) {
    bins <- tryCatch(extract_bins_with_diagnostics(ebm, feat, strict_continuous_cuts = FALSE),
                     error = function(e) list(success = FALSE))
    if (!bins$success) next

    new_scores <- rep(0.0, bins$n_scores)
    mi <- if (is.finite(bins$missing_idx)) as.integer(bins$missing_idx) else -1L
    ui <- if (is.finite(bins$unknown_idx)) as.integer(bins$unknown_idx) else NULL

    delta <- tryCatch({
      d <- py_update_term(ebm, as.integer(bins$term_index),
                          reticulate::r_to_py(new_scores),
                          mean_preserve = TRUE, preserve_missing = FALSE, missing_index = mi,
                          preserve_unknown = FALSE, unknown_index = ui,
                          shift_bagged_scores = TRUE, shift_bagged_intercept = TRUE)
      suppressWarnings(as.numeric(reticulate::py_to_r(d)))
    }, error = function(e) NA_real_)

    if (is.finite(delta)) { shift_total <- shift_total + delta; n_zeroed <- n_zeroed + 1L }
  }

  # p_post
  p_post   <- predict_proba_pos_r(ebm, X_ext)
  auc_pre  <- auc_pair(y_ext, p_pre)
  auc_post <- auc_pair(y_ext, p_post)

  data.table(
    iter_tag = iter_tag, external = ext_site, scenario = scenario_label,
    n_zeroed = n_zeroed,
    AUROC_pre = auc_pre["auroc"], AUROC_post = auc_post["auroc"],
    d_AUROC = auc_post["auroc"] - auc_pre["auroc"],
    AUPRC_pre = auc_pre["auprc"], AUPRC_post = auc_post["auprc"],
    d_AUPRC = auc_post["auprc"] - auc_pre["auprc"],
    Brier_pre = brier_score(y_ext, p_pre), Brier_post = brier_score(y_ext, p_post),
    d_Brier = brier_score(y_ext, p_post) - brier_score(y_ext, p_pre)
  )
}

# =============================================================================
# Execute both scenarios
# =============================================================================

iter_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
iter_dirs <- iter_dirs[grepl("^Iter[0-9]+_External_", basename(iter_dirs))]
iter_dirs <- iter_dirs[!grepl("Posthoc|ShapeQC|Sensitivity|NRI|TV|temporal", basename(iter_dirs))]

# Lazy-loading site data cache (<<- used inside run_zeroing_scenario to
# avoid reloading per iteration; safe because scenarios run sequentially)
site_cache  <- list()
all_results <- data.table()

for (sc in SCENARIOS_ZERO) {
  cat(sprintf("\n  === %s ===\n", sc$label))

  for (iter_dir in iter_dirs) {
    row <- tryCatch(
      run_zeroing_scenario(sc$features, iter_dir, sc$id),
      error = function(e) NULL
    )
    if (!is.null(row)) all_results <- rbind(all_results, row, fill = TRUE)
  }
}

# =============================================================================
# Output: Supp S12 and S15
# =============================================================================

OUT_09 <- file.path(RUN_DIR, paste0("Zeroing_Comparison_", ANALYSIS_TS))
dir.create(OUT_09, recursive = TRUE, showWarnings = FALSE)

fwrite(all_results, file.path(OUT_09, "All_Zeroing_Results.csv"))

# Per-scenario summaries
for (sc in SCENARIOS_ZERO) {
  sc_dt <- all_results[scenario == sc$id]
  if (nrow(sc_dt) == 0) next

  summary_dt <- sc_dt[, .(
    n_iter = .N,
    mean_d_AUROC = mean(d_AUROC, na.rm = TRUE), sd_d_AUROC = sd(d_AUROC, na.rm = TRUE),
    mean_d_AUPRC = mean(d_AUPRC, na.rm = TRUE), sd_d_AUPRC = sd(d_AUPRC, na.rm = TRUE),
    mean_d_Brier = mean(d_Brier, na.rm = TRUE),
    mean_AUROC_post = mean(AUROC_post, na.rm = TRUE)
  )]

  fwrite(summary_dt, file.path(OUT_09, sprintf("Supp_%s_Summary.csv", sc$supp_table)))

  cat(sprintf("\n  --- Supp Table %s: %s ---\n", sc$supp_table, sc$label))
  cat(sprintf("  Mean ΔAUROC: %+.4f ± %.4f\n", summary_dt$mean_d_AUROC, summary_dt$sd_d_AUROC))
  cat(sprintf("  Mean ΔAUPRC: %+.4f ± %.4f\n", summary_dt$mean_d_AUPRC, summary_dt$sd_d_AUPRC))
  cat(sprintf("  Mean AUROC_post: %.4f\n", summary_dt$mean_AUROC_post))
}

# =============================================================================
# τ threshold sensitivity comparison (S15 key result)
# =============================================================================

cat("\n\n  ─── τ Threshold Sensitivity (Supp S15) ───\n")

# Compare Red-only zeroing (τ=0.70, from 07 output) with S15 (all non-Green, τ=0.50)
# Load governance results from 07_posthoc_governance output (v3.0 or legacy v2.2)
v22_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
v22_dirs <- v22_dirs[grepl("Posthoc_v3\\.0_Governance|Posthoc_v2\\.2_4Tier", basename(v22_dirs))]

if (length(v22_dirs) > 0) {
  v22_perf_file <- file.path(sort(v22_dirs, decreasing = TRUE)[1], "Table_4B_Performance.csv")
  if (file.exists(v22_perf_file)) {
    v22_perf <- fread(v22_perf_file)
    s15_dt   <- all_results[scenario == "S15_MaxParsimony"]

    if (nrow(v22_perf) > 0 && nrow(s15_dt) > 0) {
      auroc_v22 <- mean(v22_perf$AUROC_post, na.rm = TRUE)
      auroc_s15 <- mean(s15_dt$AUROC_post, na.rm = TRUE)
      diff_auroc <- abs(auroc_v22 - auroc_s15)

      cat(sprintf("  τ=0.70 (v2.2, 6 Red zeroed):     mean AUROC_post = %.4f\n", auroc_v22))
      cat(sprintf("  τ=0.50 (S15, 10 non-Green zeroed): mean AUROC_post = %.4f\n", auroc_s15))
      cat(sprintf("  |ΔAUROC difference| = %.4f\n", diff_auroc))

      cat("\n  [Verify] Manuscript states: |ΔAUROC difference| = 0.0007\n")
      cat(sprintf("  Got: %.4f\n", diff_auroc))

      tau_comparison <- data.table(
        threshold = c("tau_0.70 (standard)", "tau_0.50 (parsimony)"),
        n_features_zeroed = c(6, 10),
        mean_AUROC_post = c(auroc_v22, auroc_s15),
        abs_AUROC_diff = diff_auroc
      )
      fwrite(tau_comparison, file.path(OUT_09, "Supp_S15_Tau_Comparison.csv"))
    }
  }
} else {
  cat("  [Note] v2.2 results not found — run 07_posthoc_governance.R first for comparison.\n")
  cat("  S15 standalone results saved.\n")
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 09 complete.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("  Output: %s\n\n", OUT_09))
