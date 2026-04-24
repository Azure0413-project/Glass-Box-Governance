# IDH EBM Governance — Reproducibility Pipeline

Code repository for:

> **"Audit and Edit: Feature-Level Governance for Transportable
> Intradialytic Hypotension Prediction in Multicenter Hemodialysis"**

---

## Repository Structure

```
EBM_Paper/
|
|-- _init.R                           # Project bootstrap (sets PROJ_ROOT, provides src())
|-- run_all.R                         # Master script: runs modules 01-11 in order
|-- run_el.R                          # Master script: E/L confirmatory pipeline
|-- run_03_eonly.R                    # Standalone launcher for E-only IECV tuning
|-- run_Eonly_update.R                # v3.0.2 E-only alignment update orchestrator
|-- extract_and_plot.py               # EBM shape function extraction & visualization
|-- README.md
|-- .gitignore
|
|-- R/                                # Shared infrastructure & utilities
|   |-- 00_config.R                   #   Global settings (paths, seeds, features, tiers)
|   |-- 00_config_el.R                #   E/L temporal split settings
|   |-- 00_utils_r.R                  #   Shared R functions (metrics, bootstrap, NRI/IDI)
|   |-- 00_utils_python.R             #   Python/reticulate EBM interface functions
|   |-- extract_ebm_importance.R      #   Extract feature importance from trained EBMs
|   |-- compute_kendall_w.R           #   Kendall's W rank stability analysis (Supp S2)
|
|-- pipeline/                         # Analysis modules (execute in order)
|   |-- 01_data_prep.R                #   Data loading & preprocessing      -> Table 1A, Supp S7
|   |-- 02_heterogeneity.R            #   Cross-site heterogeneity           -> Table 1B, Fig 2
|   |-- 03_iecv_ebm.R                 #   E-only IECV-EBM training           -> Table 2, Supp S2
|   |-- 04_iecv_xgboost.R             #   XGBoost benchmark                  -> Supp S1
|   |-- 04b_temporal_holdout.R         #   Temporal hold-out validation       -> Supp S6
|   |-- 05_shapeqc.R                  #   Shape function QC (4-tier)         -> Table 4A
|   |-- 06_tv_decomposition.R         #   Total variation decomposition      -> Supp S5
|   |-- 07_posthoc_governance.R        #   Yellow smoothing + Red zeroing     -> Table 4B
|   |-- 08_sensitivity_lofo.R          #   Sensitivity & LOFO analysis        -> Supp S8-S11
|   |-- 09_zeroing_comparison.R        #   Zero vs smooth comparison          -> Supp S12, S15
|   |-- 10_reclassification.R          #   NRI / IDI reclassification         -> Supp S13
|   |-- 11_figures.R                   #   Publication figures                -> Fig 1-5
|   |-- verify_eonly_trim.R            #   E+L vs E-only trim-bound sanity check
|   |
|   |-- el/                            #   E/L confirmatory sub-pipeline
|       |-- 04c_el_split_table.R       #     Phase 1: E/L split table
|       |-- 04d_el_derivation_iecv.R   #     Phase 2: DEPRECATED (-> 03 v3.0)
|       |-- 04e_el_shapeqc_agreement.R #     Phase 3: DEPRECATED (-> 05)
|       |-- 04f_el_confirmatory.R      #     Phase 4: L-only confirmatory scoring
|       |-- 04g_el_report.R            #     Phase 5: Report generation
|
|-- docs/                             # Development documentation
    |-- CHANGELOG.md                  #   Bug fix and patch history
    |-- READY_VERSION_NOTES.md        #   Integration notes
    |-- Supplementary_Information_v3.0.2_updates.md  # v3.0.2 manuscript patches
```

---

## Quick Start

```r
# From the project root directory:
source("run_all.R")

# Or run individual modules:
source("_init.R")              # must be sourced first
src("R/00_config.R")
src("R/00_utils_r.R")
src("R/00_utils_python.R")     # only needed for modules 03+
src("pipeline/03_iecv_ebm.R")
```

The `_init.R` bootstrap script sets the project root and provides the `src()`
function, which resolves file paths relative to the project root. All runner
scripts (`run_all.R`, `run_el.R`, `run_03_eonly.R`) automatically source
`_init.R` on startup.

---

## Pipeline Overview

### Main Pipeline (`run_all.R`)

| Module | Script | Paper Output | Runtime |
|--------|--------|-------------|---------|
| 01 | `pipeline/01_data_prep.R` | Table 1A, Supp S7 | < 1 min |
| 02 | `pipeline/02_heterogeneity.R` | Table 1B, Fig 2, Supp 1-2 | ~5 min |
| 03 | `pipeline/03_iecv_ebm.R` | Table 2, Supp S2 | **2-6 hours** |
| 04 | `pipeline/04_iecv_xgboost.R` | Supp S1 | **1-3 hours** |
| 04b | `pipeline/04b_temporal_holdout.R` | Supp S6 | ~30 min |
| 05 | `pipeline/05_shapeqc.R` | Table 4A, Supp S3, S14 | ~10 min |
| 06 | `pipeline/06_tv_decomposition.R` | Supp S5 | < 5 min |
| 07 | `pipeline/07_posthoc_governance.R` | Table 4B | ~10 min |
| 08 | `pipeline/08_sensitivity_lofo.R` | Supp S8-S11 | ~30 min |
| 09 | `pipeline/09_zeroing_comparison.R` | Supp S12, S15 | ~10 min |
| 10 | `pipeline/10_reclassification.R` | Supp S13 | ~5 min |
| 11 | `pipeline/11_figures.R` | Fig 1-5 | < 5 min |

### E/L Confirmatory Pipeline (`run_el.R`)

Separates governance-rule derivation (E = earlier 75%) from verification
(L = later 25%) to address potential circularity. A **Protocol Lock**
checkpoint ensures no L-cohort data influences governance decisions.

```
Phase 1: Build E/L split table           (pipeline/el/04c)
Phase 2: [DEPRECATED] E-only derivation  (pipeline/el/04d)
Phase 3: [DEPRECATED] ShapeQC agreement  (pipeline/el/04e)
         ======== PROTOCOL LOCK ========
Phase 4: L-only confirmatory scoring     (pipeline/el/04f)
Phase 5: Report generation               (pipeline/el/04g)
```

---

## Key Parameters

| Paper Description | Value | Config Variable |
|-------------------|-------|-----------------|
| Random seeds | 2024, 9999 | `RANDOM_SEEDS` |
| Inner CV folds | 5 | `K_INNER` |
| Sensitivity floor | >= 0.80 | `SENS_TARGET` |
| Bootstrap resamples | 1,000 | `B_BOOTSTRAP` |
| ShapeQC concordance tau | 0.70 | `SHAPEQC_TAU` |
| ShapeQC jaggedness kappa | 1.50 | `SHAPEQC_KAPPA` |
| Primary trim bounds | P0.5-P99.5 | `TRIM_CONFIGS[["T3"]]` |
| NRI thresholds | 5%, 10% | `NRI_THRESHOLDS` |
| Shape grid resolution | 150 pts | `SHAPEQC_N_GRID` |
| E/L temporal split | 75% / 25% | `EL_TEMPORAL_RATIO` |
| Harm margin | dAUPRC < -0.01 | `EL_HARM_MARGIN` |

---

## Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| R | 4.5.2 | Main analysis environment |
| Python | 3.12.12 | EBM training via reticulate |
| reticulate | 1.44.1 | R-Python bridge |
| InterpretML | 0.7.4 | Explainable Boosting Machine |
| scikit-learn | 1.8.0 | Supporting ML utilities |
| data.table | 1.18.0 | High-performance data manipulation |

### R Packages

```r
install.packages(c(
  "data.table", "fst", "reticulate", "precrec",
  "ggplot2", "gridExtra", "parallel",
  "dplyr", "tidyr", "ggridges", "openxlsx", "ggrepel"
))
```

### Python Virtual Environment

```bash
python -m venv ~/.virtualenvs/r-interpret
source ~/.virtualenvs/r-interpret/bin/activate  # Linux/Mac
# or: .virtualenvs\r-interpret\Scripts\activate  # Windows
pip install interpret==0.7.4 scikit-learn==1.8.0 numpy pandas joblib
```

---

## Data

De-identified session-level `.fst` files from three hemodialysis centers
(TN, D6, CY). File paths are configured in `R/00_config.R` via `SITE_FILES`.

---

## License

[To be specified upon publication]
