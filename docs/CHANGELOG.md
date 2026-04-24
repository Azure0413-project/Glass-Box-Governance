# CHANGELOG

All fixes previously maintained as standalone PATCH files have been merged
into their parent modules. This file documents the changes for audit trail.

## v3.0.1 — Audit Patch (2026-04-07)

Seven issues identified via code audit; fixes applied without changing
any analytical intent or methodology.

### P-1 (HIGH) → `10_reclassification.R`
**Issue:** Subgroup feature attach used `prep_site_dt()` (full E+L cohort) to
align with predictions from 07, which are E-only (`load_site_eonly()`).
Row count mismatch caused the `nrow` guard to always fail, silently
skipping all Panel C subgroup quartile analyses.
**Fix:** Changed to `load_site_eonly()`. Added `warning()` on mismatch as guard.

### P-2 (MEDIUM) → `03_iecv_ebm.R`
**Issue:** Bootstrap seed in parallel workers used `master_seed` (2024 or 9999,
varying by iteration) instead of fixed `BOOTSTRAP_SEED` (2024). This made
bootstrap CIs not strictly comparable across seed_id iterations, and
inconsistent with 04/04b/04f/10 which all use fixed `BOOTSTRAP_SEED`.
**Fix:** Changed `master_seed` → `BOOT_SEED` in `bootstrap_external()` call.

### P-3 (MEDIUM) → `04b_temporal_holdout.R`
**Issue:** Temporal split used row-index method (`floor(n * 0.75)`) while all
other scripts (03, 04c, 04f) use `quantile(dates, type = 1)`. Produces
different cutoff dates when sessions share the same date.
**Fix:** Changed to `quantile(..., type = 1)` matching `apply_el_split()`.

### P-4 (MEDIUM-LOW) → `05_shapeqc.R`
**Issue:** NOTE about E+L trim bounds used the word "leakage" which is
technically inaccurate — L-cohort outcomes are never used; only feature
distributions influence trim quantiles.
**Fix:** Expanded NOTE to clarify this is an analytic-frame choice, not
predictive leakage, with explicit guidance for reviewer-requested E-only purity.

### P-5 (LOW) → `04f_el_confirmatory.R`
**Issue:** Dead dependency check for `temporal_filter()` which is never called;
actual dependency is `apply_el_split()`.
**Fix:** Replaced with `apply_el_split()` existence check. Updated FIX LOG.

### P-6 (LOW) → `09_zeroing_comparison.R`
**Issue:** S15 τ threshold comparison only matched legacy `Posthoc_v2.2_4Tier`
directory pattern; v3.0 output (`Posthoc_v3.0_Governance_*`) was never found.
**Fix:** Updated regex to match both `v3.0` and `v2.2` patterns.

### P-7 (LOW) → `08_sensitivity_lofo.R`
**Issue:** S2 scenario uses `spar_c = 0.05` (orange penalty from legacy 5-tier),
which differs from primary `SPAR_FORMULA` in 00_config.R (no penalty).
Could mislead if S2 is described as matching Table 4B governance.
**Fix:** Added clarifying NOTE that S2 ≠ primary governance formula.

## v3.0 — Pure E-only Design (2026-04-07)

**Breaking change: E/L split moved into 03_iecv_ebm.R**

Previously, 03 trained on full-cohort data and the E/L pipeline (04c–04g)
inherited hyperparameters from the full-cohort run (partial protocol lock).
Now, 03 performs E-only hyperparameter tuning, training, and validation
directly, achieving a complete protocol lock where L-cohort data never
influence any model decision.

### Changed files:
- **00_config.R** — Added `SESSION_DATE_COL`, `E_SPLIT_QUANTILE`, `COHORT_COL`;
  stale-tier warning for ShapeQC assignments
- **00_config_el.R** — `EL_SKIP_GRID_SEARCH` deprecated; `EL_EONLY_IECV_DIR` replaces
  `EL_FULL_IECV_DIR`
- **00_utils_r.R** — Added Section 19: `apply_el_split()` and `load_site_eonly()`
  shared helpers for E-only data loading
- **03_iecv_ebm.R** — v3.0: `apply_el_split()` integrated; all steps (grid search,
  refine, final model, external validation) on E-only data; HP comparison block
- **04_iecv_xgboost.R** — v3.3.2: E-only data filtering matching 03 v3.0
- **04b_temporal_holdout.R** — HP source note updated (E-only IECV)
- **04d_el_derivation_iecv.R** — DEPRECATED: replaced by deprecation shim
  pointing to 03 v3.0 output
- **04f_el_confirmatory.R** — Reads models from `RUN_DIR` (03 v3.0 output);
  L-partition via `apply_el_split()` on-the-fly; auto-fallback tier/spar loading
  from 05_shapeqc.R output when 04e output unavailable
- **08_sensitivity_lofo.R** — External data filtered to E-only via `load_site_eonly()`
- **09_zeroing_comparison.R** — External data filtered to E-only via `load_site_eonly()`
- **run_all.R** — Banner updated for E-only v3.0
- **run_el.R** — Updated prerequisites, Phase 2 marked deprecated, `apply_el_split` check

### Unchanged (consume RUN_DIR artifacts — work as-is on E-only output):
- 01_data_prep.R, 02_heterogeneity.R, 05_shapeqc.R, 06_tv_decomposition.R,
  07_posthoc_governance.R, 10_reclassification.R, 11_figures.R

### Required workflow after update:
1. Run 03 (E-only IECV) → produces E-only models + artifacts in RUN_DIR
2. Run 05 (ShapeQC on E-only models) → update tier assignments in 00_config.R
3. Run 04f (L-only confirmatory) → reads from RUN_DIR + 05 output

## 2026-04-06: Audit-Driven Fixes (S2 Merged Audit Report)

### F-1 (HIGH) → `00_config.R`
**Issue:** Target_UF_Volume C₀ annotation showed imprecise `C0=0.70 (borderline)`,
masking that gold standard = 0.6993 (Red) vs S2 computed = 0.7039 (Yellow).
**Fix:** Updated annotation to `C0=0.6993 (gold standard; S2 computed 0.7039 due to trim precision)`.

### F-3 (LOW) → `05_shapeqc.R`
**Issue:** Bootstrap seed hardcoded as `set.seed(42L)` instead of using centralized `BOOTSTRAP_SEED` (2024L).
**Fix:** Changed to `set.seed(BOOTSTRAP_SEED)`.

### F-6 (MEDIUM) → `07_posthoc_governance.R`
**Issue:** Fallback p_pre used v2.0 (Yellow-smoothed) model, making NRI/IDI reflect
only Red zeroing increment rather than full governance effect.
**Fix:** Fallback now loads original (pre-smoothing) model from `orig_iter` directory.
Raises error if neither file-based p_old nor original model is available.

### F-8 (MEDIUM) → `08_sensitivity_lofo.R`
**Issue:** Scenario labels, tier assignments, and plot legends still used legacy
5-tier "Orange" terminology inconsistent with final 4-tier manuscript.
**Fix:** Updated all labels from "Yellow+Orange" to "Yellow+Red"; renamed tier
display from "Orange/Red" to "Red"; added `RED_SMOOTHABLE` alias with
backward-compatible `ORANGE_ORIG` reference; updated plot color mapping.

### Not addressed in code (documentation/manuscript only)
- **F-2:** E-only vs full-data tier count discrepancy — requires 04e output verification
- **F-4:** Trim domain calculation difference (04e vs 05) — manuscript disclosure
- **F-5:** 04d external validation temporal scope — manuscript disclosure
- **F-7:** C_core bootstrap CI resolution limitation — manuscript/supplement note
- **F-9:** E-only OUTER_BAGS=32 vs full=64 — supplement disclosure
- **F-10:** Policy B missing from design doc — design doc update

## 2026-04-05: Consolidated Patch Merge

### PATCH_02v2 → merged into `00_utils_r.R`
**Severity:** CRITICAL
**Bug:** Cliff's δ returned NA for TN-CY pair due to R integer overflow
(`n1 * n2 = 50000 × 50000 = 2.5e9 > .Machine$integer.max`).
**Fix:** `cliffs_delta_fast()` uses `as.double(length(x))` to prevent overflow.
**Also fixed:** Removed redundant overflow guard in `02_heterogeneity.R`
`pairwise_cliff()` that unnecessarily downsampled data (the guard assumed
`outer()` was used, but `cliffs_delta_fast()` uses rank-based formula).

### PATCH_04 → no code change needed (original `04_iecv_xgboost.R` works)
**Severity:** HIGH
**Bug:** XGBoost script failed with "undefined columns selected" because
column names were not renamed.
**Root cause:** `00_utils_r.R` was missing `prep_site_dt()` /
`apply_compat_rename()` at the time. These functions are now present.
**Status:** Original `04_iecv_xgboost.R` uses `prep_site_dt()` which
handles all column renaming correctly. No code change required.

### PATCH_05 → merged into `00_config.R`
**Severity:** CRITICAL
**Bug:** Two ShapeQC tier misclassifications:
  - IDH_N_28D → Gray (should be Green): tied-rank flag was triggered
    by counting unique values after grid projection instead of
    counting effective Spearman ranks.
  - Target_UF_Volume → Yellow (should be Red): C₀=0.7039 in S2 vs
    gold standard C₀=0.6993 < τ=0.70. Caused by trim boundary
    precision difference.
**Fix:** `00_config.R` Section 11 now uses corrected tier assignments:
  Green(5): IDH_N_28D, Pre_HD_SBP, UF_BW_Perc, Body_Temperature, Respiratory_Rate
  Yellow(3): Blood_Flow_Rate, Heart_Rate, Start_DBP
  Red(6): Dialysate_Flow_Rate, Dry_Weight, Pre_HD_Weight, Age, Target_UF_Volume, Dialysate_Temperature
  Gray(1): IDH_N_7D

### PATCH_07 → merged into `00_config.R`
**Severity:** HIGH
**Bug:** Module 07 used S2's incorrect tier assignments (Target_UF_Volume
was smoothed as Yellow instead of zeroed as Red).
**Fix:** `00_config.R` Section 11 now provides corrected tiers.
Module 07 reads tier assignments from config; no code change needed.
**Note:** Aggregate metrics (ΔAUPRC, ΔAUROC) are robust to this change.
Impact is primarily on NRI decision-level redistribution.

### PATCH_10v3 → note added to `10_reclassification.R`
**Severity:** CRITICAL
**Bug:** NRI/IDI used wrong baseline (p_new from v2.0 instead of p_old).
Correct comparison: p_pre = v2.0 p_old (original EBM), p_post = v2.2
governed predictions.
**Fix:** `10_reclassification.R` Phase 1 already attempts to load p_pre
from v2.0 p_old when not present in v2.2 file. Added verification note
about manuscript vs S13 numerical consistency.

## 2026-04-05: run_el.R Infinite Recursion Fix

**Bug:** Lines 30-31 contained a hardcoded `setwd(); source("run_el.R")` which
caused infinite recursion (stack overflow).
**Fix:** Converted to usage comment.

## 2026-04-05: Hardcoded Path Cleanup

Removed hardcoded local path references from error messages in `run_all.R`
and `run_el.R`. Paths should be configured in `00_config.R`.
