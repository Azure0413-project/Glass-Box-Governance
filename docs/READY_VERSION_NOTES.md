# Ready-to-use integration patch

This version was built from `S2_0407_routeB_patched_v2.zip` and then patched to remove the remaining high-risk issues that blocked direct use.

## Main fixes

1. `03_iecv_ebm.R`
   - restored the original refinement fold seed (`master_seed + 1000L`)
   - updates `RUN_DIR <- OUT_ROOT` so downstream modules read the current run

2. `run_all.R`
   - refreshes `RUN_DIR` immediately after Module 03

3. `04b_temporal_holdout.R`
   - keeps the E-only temporal-holdout design
   - fixes hyperparameter inheritance so only tuned fields overwrite defaults
   - respects `SESSION_DATE_COL`
   - adds guards for empty E-cohort / empty temporal test split

4. `00_utils_r.R`, `04c_el_split_table.R`, `04f_el_confirmatory.R`
   - align E/L helpers with `SESSION_DATE_COL`
   - make ShapeQC/SPAR discovery more robust

5. `05_shapeqc.R` + `README.md`
   - add `ggrepel` to install instructions
   - fall back gracefully when `ggrepel` is unavailable

6. `07_posthoc_governance.R`
   - imports `joblib` explicitly
   - keeps dynamic tier loading from the latest ShapeQC output
   - saves edited-feature values (and session date when available) into downstream prediction files

7. `08_sensitivity_lofo.R`
   - imports `joblib` explicitly
   - fixes dynamic tier loading from ShapeQC (including `light` → `tier` rename)
   - removes the accidental hard overwrite of `RED_SMOOTHABLE`

8. `10_reclassification.R`
   - prefers feature columns already saved by Module 07
   - falls back to E-only row-order attachment only when necessary
   - loads the latest ShapeQC tier map for subgroup labeling
   - stops early if no usable governance prediction files are found

9. `09_zeroing_comparison.R`
   - imports `joblib` explicitly before loading original EBM models

10. `run_el.R`
   - protocol lock is now a real stop condition unless `PROCEED_PAST_LOCK <- TRUE`

## Remaining caveat

This patch was validated statically in the current environment. An end-to-end R execution test was not possible here because `R` / `Rscript` are not available in the container.
