# =============================================================================
# 04g_el_report.R — Confirmatory Report Generation
# IDH EBM Governance — Plan A Phase 5
#
# Outputs:
#   Confirmatory Summary Table (for manuscript Results section)
#   Protocol Lock Audit Trail
#   Manuscript language suggestions
#
# FIX LOG:
#   - Issue 6: Calibration NA values now flagged with warnings in report
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists("SITE_FILES"))        src("R/00_config.R")
if (!exists("EL_TEMPORAL_RATIO")) src("R/00_config_el.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 04g: Confirmatory Report\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# Load all results
# =============================================================================

# ── Tier-tagged directories (must match 04f output) ──
if (!exists("EL_TIER_SOURCE")) EL_TIER_SOURCE <- "eonly"
EL_TIER_TAG    <- if (EL_TIER_SOURCE == "eonly") "eonly" else "full"
EL_DIR_CONFIRM <- file.path(EL_OUT_ROOT, paste0("05_L_Only_Confirmatory_", EL_TIER_TAG))
EL_DIR_REPORT  <- file.path(EL_OUT_ROOT, paste0("06_Report_", EL_TIER_TAG))
dir.create(EL_DIR_REPORT, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Tier source: %s | Reading from: %s\n\n", EL_TIER_TAG, EL_DIR_CONFIRM))

perf   <- fread(file.path(EL_DIR_CONFIRM, "EL_Confirmatory_Performance.csv"))
mp     <- fread(file.path(EL_DIR_CONFIRM, "EL_Mean_Preserve_Verification.csv"))
nri    <- fread(file.path(EL_DIR_CONFIRM, "EL_NRI_IDI.csv"))
boot   <- fread(file.path(EL_DIR_CONFIRM, "EL_Bootstrap_CI.csv"))
split  <- fread(file.path(EL_DIR_SPLIT, "EL_Split_Table.csv"))

# v3.0: Agreement and decision files are optional (04e deprecated)
agree_path <- file.path(EL_DIR_SHAPEQC, "EL_Morphological_Agreement.csv")
decision_path <- file.path(EL_DIR_SHAPEQC, "EL_Expert_Review_Decision.csv")

if (file.exists(agree_path)) {
  agree <- fread(agree_path)
} else {
  cat("  [v3.0] Morphological agreement file not found (04e deprecated) — skipping.\n")
  agree <- NULL
}

if (file.exists(decision_path)) {
  decision <- fread(decision_path)
} else {
  cat("  [v3.0] Expert review decision not found (04e deprecated) — using defaults.\n")
  decision <- data.table(
    decision = "N/A (pure E-only design — 04e deprecated)",
    all_tiers_match = NA,
    min_curve_rho = NA_real_
  )
}

# Load locked tier assignments (match what 04f actually used)
if (EL_TIER_SOURCE == "eonly") {
  qc_tiers_path <- file.path(EL_DIR_SHAPEQC, "EL_Eonly_QC_Tiers.csv")
  if (file.exists(qc_tiers_path)) {
    el_qc <- fread(qc_tiers_path)
  } else {
    # v3.0 fallback: read from 05_shapeqc.R output
    qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
    qc_dirs <- qc_dirs[grepl("ShapeQC", basename(qc_dirs))]
    if (length(qc_dirs) > 0) {
      qc_csv <- file.path(qc_dirs[length(qc_dirs)], "Table_4A_Decision_Matrix.csv")
      if (file.exists(qc_csv)) {
        cat("  [v3.0] Reading tiers from 05_shapeqc.R output\n")
        el_qc <- fread(qc_csv)
        setnames(el_qc, "light", "tier", skip_absent = TRUE)
      } else {
        stop("No tier source found. Run 05_shapeqc.R first.")
      }
    } else {
      stop("No ShapeQC output found.")
    }
  }
  REPORT_YELLOW <- el_qc[tier == "Yellow", feature]
  REPORT_RED    <- el_qc[tier == "Red",    feature]
  REPORT_GREEN  <- el_qc[tier == "Green",  feature]
  REPORT_GRAY   <- el_qc[tier == "Gray",   feature]
  REPORT_TIER_SRC <- "E-only (05_shapeqc.R)"
} else {
  REPORT_YELLOW <- YELLOW_FEATURES
  REPORT_RED    <- RED_FEATURES
  REPORT_GREEN  <- GREEN_FEATURES
  REPORT_GRAY   <- GRAY_FEATURES
  REPORT_TIER_SRC <- "Full-data (00_config.R)"
}

# =============================================================================
# Table: Confirmatory Summary (Policy A, primary seed)
# =============================================================================

cat("[Report 1] Confirmatory Summary Table\n\n")

pol_A_s1 <- perf[policy == "A" & grepl("Seed1", iter_tag)]

report_table <- pol_A_s1[, .(
  Center  = external,
  N_L     = NA_integer_,  # will fill from split table
  AUPRC_pre  = round(AUPRC_pre, 4),
  AUPRC_post = round(AUPRC_post, 4),
  d_AUPRC    = sprintf("%+.4f", d_AUPRC),
  AUROC_pre  = round(AUROC_pre, 4),
  AUROC_post = round(AUROC_post, 4),
  d_AUROC    = sprintf("%+.4f", d_AUROC),
  d_Brier    = sprintf("%+.4f", d_Brier),
  Sens_post  = sprintf("%.1f%%", Sens_post_E * 100),
  Spec_post  = sprintf("%.1f%%", Spec_post_E * 100),
  Threshold_E = round(threshold_E, 4)
)]

for (i in seq_len(nrow(report_table))) {
  report_table$N_L[i] <- split[Center == report_table$Center[i], L_sessions]
}

# Add bootstrap CI
if (nrow(boot) > 0) {
  boot_s1 <- boot[grepl("Seed1", iter_tag)]
  for (i in seq_len(nrow(report_table))) {
    ctr <- report_table$Center[i]
    b_row <- boot_s1[external == ctr]
    if (nrow(b_row) > 0) {
      report_table[i, d_AUPRC_CI := sprintf("[%+.4f, %+.4f]",
                                             b_row$d_AUPRC_lo, b_row$d_AUPRC_hi)]
    }
  }
}

# Summary row
mean_row <- data.table(
  Center = "Mean (unweighted)",
  N_L = sum(report_table$N_L, na.rm = TRUE),
  AUPRC_pre = round(mean(pol_A_s1$AUPRC_pre), 4),
  AUPRC_post = round(mean(pol_A_s1$AUPRC_post), 4),
  d_AUPRC = sprintf("%+.4f", mean(pol_A_s1$d_AUPRC)),
  AUROC_pre = round(mean(pol_A_s1$AUROC_pre), 4),
  AUROC_post = round(mean(pol_A_s1$AUROC_post), 4),
  d_AUROC = sprintf("%+.4f", mean(pol_A_s1$d_AUROC)),
  d_Brier = sprintf("%+.4f", mean(pol_A_s1$d_Brier)),
  Sens_post = "",
  Spec_post = "",
  Threshold_E = NA_real_
)
report_full <- rbind(report_table, mean_row, fill = TRUE)

fwrite(report_full, file.path(EL_DIR_REPORT, "Confirmatory_Summary_Table.csv"))
cat("  Confirmatory Summary Table:\n")
print(report_full[, .(Center, N_L, AUPRC_pre, AUPRC_post, d_AUPRC, d_AUROC)])

# =============================================================================
# Table: Calibration Summary (FIX Issue 6: flag NA values)
# =============================================================================

cat("\n[Report 1b] Calibration Summary\n\n")

cal_cols <- c("cal_slope_pre", "cal_slope_post", "cal_int_pre", "cal_int_post",
              "cal_OE_pre", "cal_OE_post")
has_cal <- all(cal_cols %in% names(pol_A_s1))

if (has_cal) {
  cal_summary <- pol_A_s1[, .(
    Center = external,
    slope_pre  = round(cal_slope_pre, 3),
    slope_post = round(cal_slope_post, 3),
    int_pre    = round(cal_int_pre, 3),
    int_post   = round(cal_int_post, 3),
    OE_pre     = round(cal_OE_pre, 3),
    OE_post    = round(cal_OE_post, 3)
  )]

  # FIX Issue 6: Check for NA calibration values
  cal_na_count <- sum(is.na(unlist(cal_summary[, -1])))
  if (cal_na_count > 0) {
    cat(sprintf("  [WARN] %d NA values in calibration summary (glm may have failed)\n", cal_na_count))
  }

  fwrite(cal_summary, file.path(EL_DIR_REPORT, "Calibration_Summary.csv"))
  print(cal_summary)
} else {
  cat("  [WARN] Calibration columns not found in performance table\n")
}

# =============================================================================
# Table: NRI Summary
# =============================================================================

if (nrow(nri) > 0) {
  cat("\n[Report 2] NRI Summary\n\n")

  nri_s1 <- nri[grepl("Seed1", iter_tag)]
  nri_summary <- nri_s1[, .(
    Center = external,
    cat_NRI = sprintf("%+.4f", cat_NRI),
    NRI_event = sprintf("%+.4f", NRI_event),
    NRI_nonevent = sprintf("%+.4f", NRI_nonevent),
    reclass_pct = sprintf("%.2f%%", reclass_rate_pct),
    IDI = sprintf("%+.6f", IDI)
  )]

  nri_mean <- data.table(
    Center = "Mean",
    cat_NRI = sprintf("%+.4f", mean(nri_s1$cat_NRI)),
    NRI_event = sprintf("%+.4f", mean(nri_s1$NRI_event)),
    NRI_nonevent = sprintf("%+.4f", mean(nri_s1$NRI_nonevent)),
    reclass_pct = sprintf("%.2f%%", mean(nri_s1$reclass_rate_pct)),
    IDI = sprintf("%+.6f", mean(nri_s1$IDI))
  )
  nri_full <- rbind(nri_summary, nri_mean, fill = TRUE)
  fwrite(nri_full, file.path(EL_DIR_REPORT, "NRI_Summary.csv"))
  print(nri_full)
}

# =============================================================================
# Table: Seed sensitivity (primary vs. seed 2)
# =============================================================================

cat("\n[Report 3] Seed Sensitivity\n\n")

pol_A_all <- perf[policy == "A"]
seed_comp <- pol_A_all[, .(
  iter_tag, external,
  seed = ifelse(grepl("Seed1", iter_tag), "Primary", "Sensitivity"),
  d_AUPRC = round(d_AUPRC, 4),
  d_AUROC = round(d_AUROC, 4)
)]
fwrite(seed_comp, file.path(EL_DIR_REPORT, "Seed_Sensitivity.csv"))
print(seed_comp)

# =============================================================================
# Table: Policy comparison (A vs. B)
# =============================================================================

cat("\n[Report 4] Policy Comparison\n\n")

pol_comp <- perf[grepl("Seed1", iter_tag), .(
  external, policy,
  d_AUPRC = round(d_AUPRC, 4),
  d_AUROC = round(d_AUROC, 4)
)]
pol_comp_wide <- dcast(pol_comp, external ~ policy, value.var = c("d_AUPRC", "d_AUROC"))
fwrite(pol_comp_wide, file.path(EL_DIR_REPORT, "Policy_Comparison.csv"))
print(pol_comp_wide)

# =============================================================================
# Protocol Lock Audit Trail
# =============================================================================

cat("\n[Report 5] Protocol Lock Audit Trail\n\n")

audit <- sprintf(
"=========================================================================
PROTOCOL LOCK AUDIT TRAIL
Generated: %s
=========================================================================

TEMPORAL SPLIT
  Ratio: %.0f%% E / %.0f%% L
  TN: split %s | E=%s | L=%s
  D6: split %s | E=%s | L=%s
  CY: split %s | E=%s | L=%s

HYPERPARAMETERS: Tuned within E-only IECV (03 v3.0 — complete protocol lock)

E-ONLY SHAPEQC
  Expert review decision: %s
  Tiers all match: %s
  Min curve rho: %s

LOCKED GOVERNANCE (Policy A)
  Tier source: %s
  Green  (%d): %s
  Yellow (smoothed): %s
  Red (zeroed): %s
  Gray   (%d): %s
  Spar formula: clamp(0.30 + 0.15*J, 0.50, 0.90)
  Domain: UPDATE_ONLY_PHYS_RANGE=TRUE, CLIP=TRUE(0.25)

MEAN-PRESERVE VERIFICATION
  Policy A: %d/%d pass (max disc = %.2e)

SUCCESS CRITERIA EVALUATION
  Primary metric: ΔAUPRC (post - pre) on L
  Center-unweighted mean ΔAUPRC: %+.4f
  Non-inferiority (≥ 0): %s
  Material harm (< %.3f): %s
  MP verification: %s

=========================================================================
",
format(Sys.time(), "%Y-%m-%d %H:%M"),
EL_TEMPORAL_RATIO * 100, (1 - EL_TEMPORAL_RATIO) * 100,
split$Split_Date[1], format(split$E_sessions[1], big.mark=","), format(split$L_sessions[1], big.mark=","),
split$Split_Date[2], format(split$E_sessions[2], big.mark=","), format(split$L_sessions[2], big.mark=","),
split$Split_Date[3], format(split$E_sessions[3], big.mark=","), format(split$L_sessions[3], big.mark=","),
decision$decision,
ifelse(is.na(decision$all_tiers_match), "N/A", ifelse(decision$all_tiers_match, "YES", "NO")),
ifelse(is.na(decision$min_curve_rho), "N/A (04e deprecated)", sprintf("%.3f", decision$min_curve_rho)),
REPORT_TIER_SRC,
length(REPORT_GREEN), paste(REPORT_GREEN, collapse = ", "),
paste(REPORT_YELLOW, collapse = ", "),
paste(REPORT_RED, collapse = ", "),
length(REPORT_GRAY), paste(if (length(REPORT_GRAY) > 0) REPORT_GRAY else "(none)", collapse = ", "),
sum(mp[policy == "A", status == "PASS"]),
nrow(mp[policy == "A"]),
max(mp[policy == "A", disc]),
mean(pol_A_s1$d_AUPRC),
ifelse(mean(pol_A_s1$d_AUPRC) >= 0, "PASS", "FAIL"),
EL_HARM_MARGIN,
ifelse(all(pol_A_s1$d_AUPRC >= EL_HARM_MARGIN), "PASS (no center harmed)", "FAIL"),
ifelse(all(mp[policy == "A", status == "PASS"]), "ALL PASS", "SOME FAILED")
)

writeLines(audit, file.path(EL_DIR_REPORT, "Protocol_Lock_Audit_Trail.txt"))

# =============================================================================
# Manuscript Language Suggestions
# =============================================================================

cat("\n[Report 6] Manuscript language suggestions\n\n")

methods_text <- '
--- METHODS: Two-stage design (new paragraph) ---

"To address the potential circularity between governance-rule derivation
and verification, we implemented a two-stage temporal split design.
Within each center, sessions were divided at the 75th session-date
percentile into an earlier derivation set (E, 75%) and a later
confirmation set (L, 25%), matching the supplementary temporal hold-out
(Table S6). In Stage 1, E-only IECV was conducted: for each held-out
center, the development pool comprised E sessions from the remaining
two centers, with hyperparameters inherited from the full-cohort IECV.
Shape-function QC, tier assignment, and the governance protocol (Policy A:
smooth Yellow, zero Red) were finalized exclusively on E-only derivation
results and locked in a frozen analysis plan before any L data were
examined. In Stage 2, the locked governance protocol was applied to E-only
models and evaluated on the L partition of each held-out center. Primary
confirmatory endpoints were ΔAUPRC and ΔAUROC (post-edit minus pre-edit);
success required center-unweighted mean ΔAUPRC ≥ 0 with no single center
exhibiting ΔAUPRC < −0.01."

--- RESULTS: Confirmatory subsection ---

"The E-only shape functions showed high morphological agreement with
full-cohort shapes (per-feature Spearman ρ ≥ [min_rho]; all tier
assignments identical), and the original expert review was retained.
On independent temporal confirmation (L cohorts), the locked governance
protocol [preserved/improved] discrimination (mean ΔAUPRC [value];
mean ΔAUROC [value]) with [X]% of sessions reclassified across
actionable thresholds (categorical NRI [value])."

--- DISCUSSION: Language calibration ---

Replace: "prespecified" for tier boundaries and zeroing decisions
With:    "empirically calibrated during the derivation stage and locked
          before confirmatory evaluation"

Replace: "deployment readiness"
With:    "retrospective derivation with independent temporal confirmation"

S8 framing: "Supplementary Table S8 reports the derivation-stage
scenario analysis that informed governance decisions; these results
were generated exclusively from E-only data."
'

writeLines(methods_text, file.path(EL_DIR_REPORT, "Manuscript_Language_Suggestions.txt"))
cat("  Written to: Manuscript_Language_Suggestions.txt\n")

cat(sprintf("\n  All reports: %s\n", EL_DIR_REPORT))

# =============================================================================
# Cross-comparison: E-only vs Full-data tiers (if both exist)
# =============================================================================

dir_eonly <- file.path(EL_OUT_ROOT, "05_L_Only_Confirmatory_eonly")
dir_full  <- file.path(EL_OUT_ROOT, "05_L_Only_Confirmatory_full")
both_exist <- file.exists(file.path(dir_eonly, "EL_Confirmatory_Performance.csv")) &&
              file.exists(file.path(dir_full,  "EL_Confirmatory_Performance.csv"))

if (both_exist) {
  cat("\n[Report 7] Cross-comparison: E-only tiers vs Full-data tiers\n\n")

  perf_e <- fread(file.path(dir_eonly, "EL_Confirmatory_Performance.csv"))
  perf_f <- fread(file.path(dir_full,  "EL_Confirmatory_Performance.csv"))

  # Policy A, primary seed only
  pe <- perf_e[policy == "A" & grepl("Seed1", iter_tag),
               .(Center = external, d_AUPRC_eonly = round(d_AUPRC, 4),
                 d_AUROC_eonly = round(d_AUROC, 4),
                 n_smooth_eonly = n_smoothed, n_zero_eonly = n_zeroed)]
  pf <- perf_f[policy == "A" & grepl("Seed1", iter_tag),
               .(Center = external, d_AUPRC_full = round(d_AUPRC, 4),
                 d_AUROC_full = round(d_AUROC, 4),
                 n_smooth_full = n_smoothed, n_zero_full = n_zeroed)]

  cross <- merge(pe, pf, by = "Center", all = TRUE)
  cross_report_dir <- file.path(EL_OUT_ROOT, "07_Cross_Comparison")
  dir.create(cross_report_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(cross, file.path(cross_report_dir, "Tier_Source_Comparison.csv"))

  cat("  Cross-comparison (Policy A, primary seed):\n")
  print(cross[, .(Center, d_AUPRC_eonly, d_AUPRC_full, d_AUROC_eonly, d_AUROC_full)])

  # Summary row
  cat(sprintf("\n  E-only mean: ΔAUPRC=%+.4f  ΔAUROC=%+.4f  (smooth=%d, zero=%d)\n",
              mean(pe$d_AUPRC_eonly), mean(pe$d_AUROC_eonly),
              pe$n_smooth_eonly[1], pe$n_zero_eonly[1]))
  cat(sprintf("  Full   mean: ΔAUPRC=%+.4f  ΔAUROC=%+.4f  (smooth=%d, zero=%d)\n",
              mean(pf$d_AUPRC_full), mean(pf$d_AUROC_full),
              pf$n_smooth_full[1], pf$n_zero_full[1]))

  both_pass <- mean(pe$d_AUPRC_eonly) >= 0 && mean(pf$d_AUPRC_full) >= 0
  cat(sprintf("  Both pass non-inferiority: %s\n", ifelse(both_pass, "YES", "NO")))

  cat(sprintf("\n  Saved to: %s\n", cross_report_dir))
}

cat("\n  Module 04g complete.\n\n")
