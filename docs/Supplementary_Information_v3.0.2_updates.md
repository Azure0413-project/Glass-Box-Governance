# Supplementary Information — v3.0.2 E-only Alignment Patch

**Patch companion for** `05_shapeqc.R` v3.0.2 and `Shape_Function_Post-Hoc_QC_v3_0r_Eonly_self-contained.R` (v3.0r.E).

**Purpose.** The revised ShapeQC pipeline computes trim bounds (P0·5–P99·5) and discrete-feature checks on the E-cohort only, closing the last remaining L-cohort dependency in the tier-assignment pipeline. The Supplementary Information text below must be updated in four places to reflect this.

Text is supplied verbatim — copy each **REPLACE** block into the corresponding location in `Supplementary_Information.docx` (or the main Supplementary file used for submission).

---

## Patch 1 — Supp Methods Part A · "Shape-function projection and computation domain"

**Location:** First full paragraph after the heading *"Shape-function projection and computation domain"* (approximately line 61 of the extracted text).

**CURRENT:**
> For each predictor, the six learned EBM shape functions (three IECV folds × two independent random seeds) were projected onto a common 150-point grid and mean-centered. All cross-fold comparisons were conducted within a robust computation domain defined by the fold-wise intersection of the pooled 0·5th–99·5th percentiles, ensuring that concordance estimates were not distorted by extrapolation into data-sparse regions.

**REPLACE WITH:**
> For each predictor, the six learned EBM shape functions (three IECV folds × two independent random seeds) were projected onto a common 150-point grid and mean-centered. All cross-fold comparisons were conducted within a robust computation domain defined by the fold-wise intersection of the pooled 0·5th–99·5th percentiles **of the E-cohort training distribution**, ensuring that concordance estimates were not distorted by extrapolation into data-sparse regions **and preserving the pure E-only protocol lock: no L-cohort session contributed to any tier-assignment computation.**

---

## Patch 2 — Supp Table S3 caption · "Cross-model shape-function correlation summary"

**Location:** Table S3 caption paragraph beginning *"Each feature's six IECV-trained shape functions were projected onto a common 150-point grid..."* (approximately line 271).

**CURRENT:**
> Each feature's six IECV-trained shape functions were projected onto a common 150-point grid and mean-centered. Pearson *r* and Spearman *ρ* were computed for each of the C(6, 2) = 15 model pairs within the quantile-trimmed domain (T3: 0.5th–99.5th percentiles of pooled data; see Methods). Tier labels follow the full-data configuration (see Analysis overview in Results); primary tier assignments are in Table 3, Panel A. …

**REPLACE WITH:**
> Each feature's six IECV-trained shape functions were projected onto a common 150-point grid and mean-centered. Pearson *r* and Spearman *ρ* were computed for each of the C(6, 2) = 15 model pairs within the quantile-trimmed domain (T3: 0.5th–99.5th percentiles of the **pooled E-cohort training distribution**; see Methods). Tier labels follow the **E-only** configuration (see Analysis overview in Results); primary tier assignments are in Table 3, Panel A. …

*(Remainder of caption unchanged.)*

---

## Patch 3 — Supp Table S14 footnote · "C_core sensitivity across five trimming configurations"

**Location:** The sentence immediately following the Table S14 body (line ~658):

**CURRENT:**
> *C__core_* = trimmed cross-fold pairwise Spearman rank concordance, computed on P0·5–P99·5 of the pooled training distribution. *C__med_*, C*__min__*, C*__max__* = median, minimum, and maximum across all fold-pairs.

**REPLACE WITH:**
> *C__core_* = trimmed cross-fold pairwise Spearman rank concordance, computed on P0·5–P99·5 of the **pooled E-cohort training distribution**. *C__med_*, C*__min__*, C*__max__* = median, minimum, and maximum across all fold-pairs.

---

## Patch 4 — Supp Table S5 footnote (optional; for complete consistency)

**Location:** Footnote under Table S5 TV decomposition (line ~330):

**CURRENT:**
> Each cell reports the percentage of total shape-function variation retained within the central [*p*, *1 – p*] weighted-sample range (core TV retained, %), where *p* is the tail trim threshold (each side). Values are mean across six IECV models; ranges in parentheses denote minimum–maximum across models.

**REPLACE WITH:**
> Each cell reports the percentage of total shape-function variation retained within the central [*p*, *1 – p*] weighted-sample range **of the E-cohort training distribution** (core TV retained, %), where *p* is the tail trim threshold (each side). Values are mean across six IECV models; ranges in parentheses denote minimum–maximum across models.

*Rationale: The TV decomposition (S5) uses the same pooled weights as the ShapeQC trim domain (S3/S14). Once the ShapeQC frame is tightened to E-only, S5 must be described the same way for internal consistency even if the numerical values change only marginally.*

---

## Patch 5 — Add to CHANGELOG.md (code repository, not the manuscript)

```markdown
## v3.0.2 — E-only Trim-Bound Alignment (2026-04-17)

Closes the last remaining L-cohort dependency in the ShapeQC tier-assignment
pipeline. Previously, trim bounds (P0·5–P99·5) for C_core and discrete-feature
checks were computed on the full E+L cohort — an analytic-frame choice that
did not involve outcomes but was structurally inconsistent with the
E-only protocol lock already enforced in modules 03 / 04 / 08 / 09 / 10.

### Changed files:
- **05_shapeqc.R** — `site_data_cache` loading switched from raw E+L read
  to `load_site_eonly()`; trim bounds and discrete-feature detection now
  use E-cohort data only. Output directory renamed
  `ShapeQC_v3.0_Eonly_<ts>` (keeps the "v3.0_" prefix so existing regex in
  07_posthoc_governance.R and 08_sensitivity_lofo.R continues to match).
- **Shape_Function_Post-Hoc_QC_v3_0r_Eonly_self-contained.R** — mirror change
  with embedded `apply_el_split_inline()` / `load_site_eonly_inline()`
  helpers (required because the self-contained script does not source
  00_utils_r.R).
- **verify_eonly_trim_alignment.R** — new diagnostic script that quantifies
  E+L vs E-only trim-bound differences per feature × fold and flags
  tier-boundary-sensitive features requiring re-inspection.
- **Supplementary_Information.docx** — four textual patches
  (Part A paragraph, Table S3 caption, Table S5 footnote, Table S14
  footnote) specifying E-cohort training distribution. See
  `Supplementary_Information_v3.0.2_updates.md`.

### Unchanged (already on E-only data):
- 08_sensitivity_lofo.R (S10)
- 09_zeroing_comparison.R / 09_zeroing_comparison_based_on_E.R (S12)
- 10_reclassification.R (S13 Panel A / B / C)

### Verification workflow:
1. `source("verify_eonly_trim_alignment.R")` to preview trim-bound shifts
   and flag any tier-boundary-sensitive features.
2. Rerun `05_shapeqc.R` (v3.0.2) → produces
   `ShapeQC_v3.0_Eonly_<ts>/Table_4A_Decision_Matrix.csv`.
3. Diff against previous `Table_4A_Decision_Matrix.csv`; confirm no tier
   flip (expected — tier assignments are robust to the 25% sample
   reduction because P0·5 / P99·5 are insensitive to moderate sample
   change).
4. If any tier flips, re-run downstream (07, 08, 09, 10) and update
   Table 4A / Figure 3 / Table 4B / S10 / S12 / S13.
5. Apply the four manuscript patches in
   `Supplementary_Information_v3.0.2_updates.md`.
```

---

## Summary of Manuscript Impact

| Location | Change | Why |
|----------|--------|-----|
| Supp Methods Part A | +"E-cohort training distribution" + pure E-only lock statement | S3 / S14 / Table 4A trim bounds now E-only |
| Supp S3 caption | "pooled data" → "pooled E-cohort training distribution"; "full-data" → "E-only" | S3 numerics derive from E-only |
| Supp S14 footnote | "pooled training distribution" → "pooled E-cohort training distribution" | S14 is C_core under 5 trim configs — same frame |
| Supp S5 footnote | Add "E-cohort" qualifier to weighted-sample range | S5 TV decomposition uses the same weights |
| Main text (Methods) | Verify no conflicting "full cohort" or "all sessions" language around ShapeQC trim bounds | Should already be consistent; check v25_clean.docx Methods §ShapeQC |

**Note on numerical impact.** Because P0·5 / P99·5 are quantile statistics at the extreme tails, removing the latest 25% of sessions (L-cohort) will perturb trim bounds by well under 1% of each feature's range in almost all cases. Tier assignments are therefore expected to be numerically near-identical to the previous run. The `verify_eonly_trim_alignment.R` script produces a per-feature quantitative confirmation before the full rerun.
