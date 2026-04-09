# =============================================================================
# 04c_el_split_table.R — Build E/L Temporal Split Table
# IDH EBM Governance — Plan A Phase 1
#
# Outputs:
#   Split table CSV — per-center E/L session counts, patient counts,
#                     date ranges, IDH prevalence
#   Frozen Analysis Plan skeleton (TXT)
#
# Plan A Step 1: build the temporal split table for TN, D6, CY
# Uses the same center-specific 75/25 split as Supplementary Table S6
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       (temporal_filter() now lives in 00_utils_r.R)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
})

if (!exists("SITE_FILES"))   src("R/00_config.R")
if (!exists("read_fst_dt"))  src("R/00_utils_r.R")
if (!exists("EL_TEMPORAL_RATIO")) src("R/00_config_el.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 04c: E/L Temporal Split Table\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# Step 1: Build split table for each center
# =============================================================================

cat("[Step 1] Building temporal split table ...\n\n")

split_results <- list()

for (site in c("TN", "D6", "CY")) {
  cat(sprintf("  === %s ===\n", site))

  # Load and preprocess
  dt <- prep_site_dt(read_fst_dt(SITE_FILES[[site]]), site)

  date_col <- SESSION_DATE_COL
  if (!(date_col %in% names(dt))) {
    stop(sprintf("[%s] %s column not found", site, date_col))
  }
  dt[, (date_col) := as.Date(get(date_col))]
  data.table::setorderv(dt, date_col)

  # 75th percentile split (by date quantile, matching apply_el_split() in 00_utils_r.R)
  cutoff_date <- as.Date(quantile(dt[[date_col]], probs = EL_TEMPORAL_RATIO,
                                   type = 1, na.rm = TRUE))

  dt_E <- dt[get(date_col) <= cutoff_date]
  dt_L <- dt[get(date_col) >  cutoff_date]

  y_E <- as.integer(dt_E[[TARGET_COL]])
  y_L <- as.integer(dt_L[[TARGET_COL]])

  # Patient overlap check
  pids_E <- unique(dt_E[[ID_COL]])
  pids_L <- unique(dt_L[[ID_COL]])
  pids_overlap <- intersect(pids_E, pids_L)

  row <- data.table(
    Center           = site,
    Split_Date       = as.character(cutoff_date),
    E_sessions       = nrow(dt_E),
    E_patients       = length(pids_E),
    E_date_start     = as.character(min(dt_E[[date_col]])),
    E_date_end       = as.character(max(dt_E[[date_col]])),
    E_IDH_events     = sum(y_E, na.rm = TRUE),
    E_IDH_prevalence = round(mean(y_E, na.rm = TRUE) * 100, 2),
    L_sessions       = nrow(dt_L),
    L_patients       = length(pids_L),
    L_date_start     = as.character(min(dt_L[[date_col]])),
    L_date_end       = as.character(max(dt_L[[date_col]])),
    L_IDH_events     = sum(y_L, na.rm = TRUE),
    L_IDH_prevalence = round(mean(y_L, na.rm = TRUE) * 100, 2),
    Overlap_patients = length(pids_overlap),
    Overlap_pct      = round(length(pids_overlap) / length(pids_L) * 100, 1)
  )
  split_results[[site]] <- row

  cat(sprintf("    Split date : %s\n", cutoff_date))
  cat(sprintf("    E: %s sessions, %d patients (%s to %s), IDH=%.2f%%\n",
              format(nrow(dt_E), big.mark = ","), length(pids_E),
              min(dt_E[[date_col]]), max(dt_E[[date_col]]),
              mean(y_E, na.rm = TRUE) * 100))
  cat(sprintf("    L: %s sessions, %d patients (%s to %s), IDH=%.2f%%\n",
              format(nrow(dt_L), big.mark = ","), length(pids_L),
              min(dt_L[[date_col]]), max(dt_L[[date_col]]),
              mean(y_L, na.rm = TRUE) * 100))
  cat(sprintf("    Patient overlap: %d (%.1f%% of L patients)\n\n",
              length(pids_overlap),
              length(pids_overlap) / length(pids_L) * 100))
}

split_table <- rbindlist(split_results, fill = TRUE)
fwrite(split_table, file.path(EL_DIR_SPLIT, "EL_Split_Table.csv"))

# =============================================================================
# Step 2: Verify consistency with S6
# =============================================================================

cat("[Step 2] Verifying consistency with S6 ...\n\n")

s6_expected <- data.table(
  Center     = c("TN",     "D6",    "CY"),
  S6_train_n = c(111335L,  31856L,  293154L),
  S6_test_n  = c(37080L,   10565L,  97479L)
)

for (i in seq_len(nrow(s6_expected))) {
  site  <- s6_expected$Center[i]
  el    <- split_table[Center == site]
  s6_tn <- s6_expected$S6_train_n[i]
  s6_ts <- s6_expected$S6_test_n[i]

  match_E <- (el$E_sessions == s6_tn)
  match_L <- (el$L_sessions == s6_ts)

  cat(sprintf("  %s: E=%s vs S6_train=%s [%s] | L=%s vs S6_test=%s [%s]\n",
              site,
              format(el$E_sessions, big.mark = ","),
              format(s6_tn, big.mark = ","),
              ifelse(match_E, "MATCH", "DIFF"),
              format(el$L_sessions, big.mark = ","),
              format(s6_ts, big.mark = ","),
              ifelse(match_L, "MATCH", "DIFF")))
}

# =============================================================================
# Step 3: IECV E-only development pool summary
# For each IECV fold, show the E-only dev pool composition
# =============================================================================

cat("\n[Step 3] E-only IECV development pool sizes ...\n\n")

iecv_pool_summary <- list()
for (it in ITERATIONS) {
  dev_E_n <- sum(split_table[Center %in% it$dev, E_sessions])
  dev_E_events <- sum(split_table[Center %in% it$dev, E_IDH_events])
  ext_L_n <- split_table[Center == it$external, L_sessions]
  ext_L_events <- split_table[Center == it$external, L_IDH_events]

  row <- data.table(
    iter        = it$iter,
    external    = it$external,
    seed_id     = it$seed_id,
    dev_sites   = paste(it$dev, collapse = "+"),
    dev_E_n     = dev_E_n,
    dev_E_prev  = round(dev_E_events / dev_E_n * 100, 2),
    ext_L_n     = ext_L_n,
    ext_L_prev  = round(ext_L_events / ext_L_n * 100, 2)
  )
  iecv_pool_summary[[length(iecv_pool_summary) + 1]] <- row

  if (it$seed_id == 1L) {
    cat(sprintf("  Fold %s (ext=%s): dev_E = %s (%.2f%%) → ext_L = %s (%.2f%%)\n",
                it$external, it$external,
                format(dev_E_n, big.mark = ","), row$dev_E_prev,
                format(ext_L_n, big.mark = ","), row$ext_L_prev))
  }
}
iecv_pool_dt <- rbindlist(iecv_pool_summary, fill = TRUE)
fwrite(iecv_pool_dt, file.path(EL_DIR_SPLIT, "EL_IECV_Pool_Summary.csv"))

# =============================================================================
# Step 4: Generate Frozen Analysis Plan skeleton
# =============================================================================

cat("\n[Step 4] Writing Frozen Analysis Plan skeleton ...\n\n")

fap_text <- sprintf(
"=========================================================================
FROZEN ANALYSIS PLAN — E/L Temporal Split with Protocol Lock
Date: %s
=========================================================================

1. TEMPORAL SPLIT
   Ratio: %.0f%% Earlier (E) / %.0f%% Later (L), center-specific
   Split method: session-date percentile (matching Supplementary Table S6)
   Patient overlap across E/L: allowed (session-level prediction)

   Center   Split Date    E sessions   L sessions   E IDH%%   L IDH%%
   TN       %s       %s       %s       %s      %s
   D6       %s       %s       %s       %s      %s
   CY       %s       %s       %s       %s      %s

2. HYPERPARAMETERS
   Source: inherited from full-cohort IECV (no re-tuning on E-only)
   Justification: HP selection did not involve L data or governance outcomes

3. GOVERNANCE TIER ASSIGNMENT
   C_core threshold (τ): %.2f
   Jaggedness threshold (κ): %.2f
   Trim: %s (P%.1f–P%.1f)
   Tied-rank safeguard: min effective unique ranks = %d

4. CANDIDATE POLICIES (evaluated on E-only derivation)
   Policy A (primary): Smooth Yellow + Zero Red
   Policy B (sensitivity): Smooth All Eligible

5. SMOOTHING SPECIFICATION
   Function: smooth_scores_spline_fixed()
   Spar formula: clamp(0.30 + 0.15 × J, 0.50, 0.90)
   UPDATE_ONLY_PHYS_RANGE: TRUE
   CLIP_TO_ORIGINAL_RANGE: TRUE (margin = 0.25)
   UPDATE_MISSING_BIN: FALSE
   UPDATE_UNKNOWN_BIN: FALSE
   Min unique x for fitting: 5
   Weighting: ebm.bin_weights_[term_index]

6. THRESHOLD SELECTION
   Rule: development OOF sensitivity ≥ 0.80 + F1 maximization
   Source: E-only development pool OOF predictions
   Sensitivity analysis: full-IECV thresholds (Supp)

7. PRIMARY CONFIRMATORY ENDPOINTS
   Primary: ΔAUPRC (post-edit − pre-edit) on L cohorts
   Success: center-unweighted mean ΔAUPRC ≥ 0
   Harm: no single center with ΔAUPRC < %.3f

8. KEY SECONDARY ENDPOINTS
   ΔAUROC, post-edit calibration (slope/intercept/O:E),
   threshold-based operating characteristics (locked threshold)

9. EXPLORATORY ENDPOINTS
   Categorical NRI, continuous NRI, IDI, DCA (5-20%%)

10. BOOTSTRAP
    Method: patient-level, per center
    Resamples: %d
    Seed: %d

11. AGGREGATION
    Center-unweighted (not sample-size-weighted)
    Primary seed (seed_id=1) for confirmatory; seed_id=2 for sensitivity

=========================================================================
LOCK TIMESTAMP: [TO BE FILLED AFTER E-ONLY DERIVATION]
This document must be frozen before any L data is accessed.
=========================================================================
",
format(Sys.time(), "%Y-%m-%d %H:%M"),
EL_TEMPORAL_RATIO * 100, (1 - EL_TEMPORAL_RATIO) * 100,
split_table[1, Split_Date], format(split_table[1, E_sessions], big.mark = ","),
format(split_table[1, L_sessions], big.mark = ","),
split_table[1, E_IDH_prevalence], split_table[1, L_IDH_prevalence],
split_table[2, Split_Date], format(split_table[2, E_sessions], big.mark = ","),
format(split_table[2, L_sessions], big.mark = ","),
split_table[2, E_IDH_prevalence], split_table[2, L_IDH_prevalence],
split_table[3, Split_Date], format(split_table[3, E_sessions], big.mark = ","),
format(split_table[3, L_sessions], big.mark = ","),
split_table[3, E_IDH_prevalence], split_table[3, L_IDH_prevalence],
SHAPEQC_TAU, SHAPEQC_KAPPA, PRIMARY_TRIM,
TRIM_CONFIGS[[PRIMARY_TRIM]][1]*100, TRIM_CONFIGS[[PRIMARY_TRIM]][2]*100,
SHAPEQC_EFF_N_MIN,
EL_HARM_MARGIN,
EL_BOOT_N, EL_BOOT_SEED
)

writeLines(fap_text, file.path(EL_DIR_SPLIT, "Frozen_Analysis_Plan.txt"))

cat(sprintf("  Split table : %s\n", file.path(EL_DIR_SPLIT, "EL_Split_Table.csv")))
cat(sprintf("  IECV pools  : %s\n", file.path(EL_DIR_SPLIT, "EL_IECV_Pool_Summary.csv")))
cat(sprintf("  FAP skeleton: %s\n", file.path(EL_DIR_SPLIT, "Frozen_Analysis_Plan.txt")))
cat("\n  Module 04c complete.\n\n")
