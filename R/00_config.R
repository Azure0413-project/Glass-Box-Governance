# =============================================================================
# 00_config.R — Global Configuration
# IDH EBM Governance Reproducibility Pipeline
#
# All shared settings used by downstream modules are defined here.
# Each script loads this file via src("R/00_config.R").
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  IDH EBM Governance — 00_config.R loaded\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# Section 1: Environment
# =============================================================================

set.seed(1107)   # NOTE: downstream scripts set their own seeds; this is a fallback
options(stringsAsFactors = FALSE)

RUN_TS <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

# =============================================================================
# Section 2: Paths (edit to match your local environment)
# =============================================================================

# Dynamic output root directory (created per run)
# Edit OUT_ROOT to your preferred output location
OUT_ROOT <- Sys.getenv("EBM_OUT_ROOT", unset = file.path(getwd(), "output"))
OUT_DIR  <- file.path(OUT_ROOT, paste0("Exp_", format(Sys.time(), "%Y_%m_%d_%H_%M")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DIR_TABLE <- file.path(OUT_DIR, "tables")
DIR_FIG   <- file.path(OUT_DIR, "figures")
DIR_SUPP  <- file.path(OUT_DIR, "supplement")
DIR_MODEL <- file.path(OUT_DIR, "models")
lapply(c(DIR_TABLE, DIR_FIG, DIR_SUPP, DIR_MODEL),
       dir.create, recursive = TRUE, showWarnings = FALSE)

# Data directory: set EBM_DATA_DIR env variable or edit this path
DATA_DIR <- Sys.getenv("EBM_DATA_DIR", unset = file.path(getwd(), "data"))

# IECV-EBM base directory for parallel runs
BASE_DIR <- Sys.getenv("EBM_BASE_DIR", unset = file.path(OUT_ROOT, "Run_Parallel"))

# IECV-EBM output directory (produced by 03_iecv_ebm.R; consumed by downstream modules)
# After running 03_iecv_ebm.R, update this path to the latest Run_* directory.
RUN_DIR <- Sys.getenv("EBM_RUN_DIR", unset = file.path(BASE_DIR, "Run_YYYY_MM_DD_HH_MM_SS"))

# Python virtual environment
VENV_PATH <- Sys.getenv(
  "R_INTERPRET_VENV",
  unset = file.path(Sys.getenv("HOME", unset = Sys.getenv("USERPROFILE", unset = "~")),
                    ".virtualenvs", "r-interpret")
)

# Raw data files (.fst) from three hemodialysis centers
# File names use anonymized center codes (TN, D6, CY)
SITE_FILES <- list(
  TN = file.path(DATA_DIR, "center_TN.fst"),
  D6 = file.path(DATA_DIR, "center_D6.fst"),
  CY = file.path(DATA_DIR, "center_CY.fst")
)

# =============================================================================
# Section 3: Column Definitions
# =============================================================================

ID_COL      <- "Patient_ID"
CLUSTER_COL <- "Patient_Cluster_ID"
TARGET_COL  <- "Nadir90/100"

BINARY_COLS <- c("Sex")

CONT_COLS <- c(
  "Age", "IDH_N_7D", "IDH_N_28D",
  "Pre_HD_SBP", "Start_DBP", "Heart_Rate", "Respiratory_Rate", "Body_Temperature",
  "Pre_HD_Weight", "Dry_Weight", "Target_UF_Volume", "UF_BW_Perc",
  "Blood_Flow_Rate", "Dialysate_Flow_Rate", "Dialysate_Temperature"
)

ALL_PREDICTORS <- c(BINARY_COLS, CONT_COLS)  # 16 predictors total

# Cross-site column name harmonization map
COMPAT_RENAME_MAP <- c(
  "blood-speed"         = "Blood_Flow_Rate",
  "Dialysis-blood-rate" = "Blood_Flow_Rate",
  "Dialysis-blood-temp" = "Dialysate_Temperature",
  "start-weight"        = "Pre_HD_Weight",
  "\u9ad4\u6eab_New"    = "Body_Temperature",
  "\u9810\u4f30\u8131\u6c34\u91cf" = "Target_UF_Volume",
  "HR"                  = "Heart_Rate",
  "RR"                  = "Respiratory_Rate",
  "UF_BW_perc"          = "UF_BW_Perc",
  "DP-start"            = "Start_DBP",
  "ID"                  = "Patient_ID"
)

# =============================================================================
# Section 4: Physiological Range Filters
# Methods -> "Values outside prespecified physiological ranges were recoded
#  as missing"
# =============================================================================

PHYS_RANGES <- list(
  Pre_HD_SBP            = c(60, 250),
  Start_DBP             = c(25, 150),
  Heart_Rate            = c(30, 150),
  Respiratory_Rate      = c(5, 40),
  Body_Temperature      = c(34, 42),
  Target_UF_Volume      = c(0, 7),
  UF_BW_Perc            = c(-0.05, 0.095),
  Blood_Flow_Rate       = c(60, 400),
  Dialysate_Flow_Rate   = c(80, 1000),
  Dialysate_Temperature = c(34, 39),
  Dry_Weight            = c(25, 200),
  Pre_HD_Weight         = c(25, 200)
)

# =============================================================================
# Section 5: IECV Iterations (3 centers x 2 seeds = 6 models)
# =============================================================================

RANDOM_SEEDS <- c(2024L, 9999L)

ITERATIONS <- list(
  list(iter = 1L, external = "TN", dev = c("D6", "CY"), seed_id = 1L, master_seed = 2024L),
  list(iter = 2L, external = "TN", dev = c("D6", "CY"), seed_id = 2L, master_seed = 9999L),
  list(iter = 3L, external = "D6", dev = c("TN", "CY"), seed_id = 1L, master_seed = 2024L),
  list(iter = 4L, external = "D6", dev = c("TN", "CY"), seed_id = 2L, master_seed = 9999L),
  list(iter = 5L, external = "CY", dev = c("TN", "D6"), seed_id = 1L, master_seed = 2024L),
  list(iter = 6L, external = "CY", dev = c("TN", "D6"), seed_id = 2L, master_seed = 9999L)
)

# =============================================================================
# Section 5.1: E/L Temporal Split
# Methods -> "Each center was split chronologically at the 75th percentile
#  of sessions ordered by session date."
# Under the pure E-only design, hyperparameter tuning, model training, and
# external validation all use E-cohort (earlier 75%) exclusively.
# L-cohort (later 25%) is reserved for downstream confirmation (04f).
# =============================================================================

SESSION_DATE_COL  <- "Session_Date"  # adjust to match actual column name
E_SPLIT_QUANTILE  <- 0.75            # E = earlier 75%, L = later 25%
COHORT_COL        <- "EL_cohort"     # added column: "E" or "L"

# =============================================================================
# Section 6: EBM Hyperparameter Grid (v2.5.10)
# Methods -> "Hyperparameters governing leaf size, smoothing rounds,
#  and positive-class weighting were tuned by five-fold grouped
#  patient-level cross-validation within the E-cohort development pool."
# NOTE: Under the pure E-only design, hyperparameter tuning uses only
#  E-cohort sessions. This ensures a complete protocol lock -- L-cohort
#  data do not influence any model decision, including hyperparameters.
# =============================================================================

EBM_FIXED_HP <- list(
  learning_rate = 0.005,
  interactions  = 0L,
  max_bins      = 128L,
  max_leaves    = 2L
)

EBM_GRID <- expand.grid(
  min_samples_leaf = c(300L, 600L, 1000L),
  smoothing_rounds = c(75L, 150L, 600L, 1200L),
  pos_multiplier   = c(1.5, 6.0, 9.0),
  stringsAsFactors = FALSE
)

# =============================================================================
# Section 7: Cross-Validation & Threshold Tuning
# Methods -> "Thresholds were selected on development-set out-of-fold
#  predictions using a prespecified sensitivity floor of 0.80
#  and F1 maximisation among eligible thresholds."
# =============================================================================

K_INNER              <- 5L        # 5-fold grouped CV within dev pool
CV_FOLD_SEED         <- 2024L
SENS_TARGET          <- 0.80
THR_GRID             <- seq(0.001, 0.999, by = 0.001)
OUTER_BAGS_OOF       <- 10L       # outer_bags for grid search
OUTER_BAGS_FINAL     <- 64L       # outer_bags for final model
GRID_SUBSAMPLE_RATIO <- 0.6       # grid search subsampling
EBM_N_JOBS           <- 5L        # threads per EBM

# =============================================================================
# Section 8: Bootstrap & Calibration
# Methods -> "Uncertainty was estimated by 1,000 patient-level
#  bootstrap resamples."
# =============================================================================

B_BOOTSTRAP    <- 1000L
BOOTSTRAP_SEED <- 2024L
BOOT_CONF      <- c(0.025, 0.975)
CAL_BINS       <- 10L
DCA_PTS        <- seq(0.01, 0.99, by = 0.01)

# =============================================================================
# Section 9: Parallel Processing
# =============================================================================

N_WORKERS             <- 6L
N_THREADS_PER_WORKER  <- 5L

# =============================================================================
# Section 10: ShapeQC Governance Thresholds (v3.0)
# Methods -> "Green: Median C_core >= 0.70 and J <= 1.50;
#  Yellow: Median C_core >= 0.70 and J > 1.50;
#  Red: Median C_core < 0.70;
#  Gray: manual review (tied-rank safeguard)"
# =============================================================================

SHAPEQC_TAU       <- 0.70    # C_core concordance threshold
SHAPEQC_KAPPA     <- 1.50    # Jaggedness J threshold
SHAPEQC_N_GRID    <- 150L    # common grid points
SHAPEQC_BOOT_N    <- 2000L   # bootstrap iterations for CI
SHAPEQC_EFF_N_MIN <- 10L     # minimum effective unique ranks (tied-rank safeguard)

# Trim configurations for sensitivity analysis (Supp Table S14)
# Methods -> "Sensitivity of tier assignment to trim boundaries was assessed
#  across five prespecified configurations"
TRIM_CONFIGS <- list(
  T1 = c(0.025, 0.975),    # P2.5-P97.5
  T2 = c(0.010, 0.990),    # P1.0-P99.0
  T3 = c(0.005, 0.995),    # P0.5-P99.5 (PRIMARY)
  T4 = c(0.001, 0.999),    # P0.1-P99.9
  T5 = c(0.000, 1.000)     # Untrimmed
)
PRIMARY_TRIM <- "T3"

# =============================================================================
# Section 11: Four-Tier Feature Governance (final manuscript classification)
# Table 4A -> assignment based on ShapeQC v3.0 Decision_Matrix
# Note: Sex (binary predictor) is excluded from ShapeQC -- shape functions
#       are not defined for binary features.
#
# IMPORTANT: The tier assignments below were derived from the FULL-COHORT
#   IECV configuration. After switching to E-only hyperparameter tuning,
#   re-run the ShapeQC pipeline (scripts 05-07) on E-only models and
#   UPDATE these assignments with the new E-only tier results.
#   Expected E-only tiers (from prior analysis): 8 Green, 2 Yellow, 5 Red, 0 Gray.
#   Until updated, downstream governance scripts will use stale tiers.
# =============================================================================

GREEN_FEATURES <- c(
  "IDH_N_28D",               # C0=1.00, J=1.00
  "Pre_HD_SBP",              # C0=0.97, J=1.05
  "UF_BW_Perc",              # C0=0.96, J=1.05
  "Body_Temperature",        # C0=0.91, J=1.43
  "Respiratory_Rate"         # C0=0.85, J=1.49
)

YELLOW_FEATURES <- c(
  "Blood_Flow_Rate",         # C0=0.89, J=2.30
  "Heart_Rate",              # C0=0.81, J=1.53
  "Start_DBP"               # C0=0.97, J=2.09
)

RED_FEATURES <- c(
  "Dialysate_Flow_Rate",     # C0=0.15
  "Dry_Weight",              # C0=-0.55
  "Pre_HD_Weight",           # C0=-0.12
  "Age",                     # C0=0.61
  "Target_UF_Volume",        # C0=0.6993 (gold standard; S2 computed 0.7039 due to trim precision)
  "Dialysate_Temperature"    # C0=0.58
)

GRAY_FEATURES <- c(
  "IDH_N_7D"                 # tied-rank safeguard; retained after manual review
)

ALL_EDIT_FEATURES <- c(YELLOW_FEATURES, RED_FEATURES)

FEATURE_TIER_MAP <- setNames(
  c(rep("Green",  length(GREEN_FEATURES)),
    rep("Yellow", length(YELLOW_FEATURES)),
    rep("Red",    length(RED_FEATURES)),
    rep("Gray",   length(GRAY_FEATURES))),
  c(GREEN_FEATURES, YELLOW_FEATURES, RED_FEATURES, GRAY_FEATURES)
)

# =============================================================================
# Section 12: Post-hoc Smoothing Parameters
# Methods -> "Smoothing used cubic splines with a prespecified J-adaptive
#  mapping from jaggedness to spar"
# =============================================================================

SPAR_FORMULA <- function(J) {
  # Applied to Yellow-tier features only (Red features are zeroed, not smoothed)
  base_spar <- 0.30 + 0.15 * J
  max(0.50, min(base_spar, 0.90))
}

MEAN_PRESERVE_TOLERANCE <- 1e-6

# =============================================================================
# Section 13: NRI / IDI Reclassification
# Methods -> "Categorical NRI was computed using prespecified risk thresholds
#  of 5% and 10%, defining low-, intermediate-, and high-risk categories."
# =============================================================================

NRI_THRESHOLDS <- c(0.05, 0.10)
NRI_LABELS     <- c("Low (<5%)", "Medium (5-10%)", "High (>10%)")

# =============================================================================
# Section 14: Cross-site Heterogeneity Sampling
# =============================================================================

SAMPLE_CONFIG <- list(
  kruskal_max_n = 100000L,
  ks_max_n      = 100000L,
  cliff_max_n   = 50000L,
  plot_max_n    = 50000L
)

# =============================================================================
# Section 15: Visualization Palette
# =============================================================================

SITE_COLORS <- c("TN" = "#3498db", "D6" = "#e74c3c", "CY" = "#2ecc71")

# 6-model IECV iteration colors (for shape function overlays)
ITER_COLORS <- c(
  "Iter1_TN_s1" = "#1f77b4",
  "Iter2_TN_s2" = "#aec7e8",
  "Iter3_D6_s1" = "#ff7f0e",
  "Iter4_D6_s2" = "#ffbb78",
  "Iter5_CY_s1" = "#2ca02c",
  "Iter6_CY_s2" = "#98df8a"
)

# Governance tier colors (for ShapeQC figures)
TIER_COLORS <- c(
  "Green"  = "#27ae60",
  "Yellow" = "#f1c40f",
  "Red"    = "#e74c3c",
  "Gray"   = "#95a5a6"
)

# =============================================================================
# Section 16: Miscellaneous
# =============================================================================

ENABLE_PLOTTING    <- TRUE
FORCE_GC_EACH_ITER <- TRUE

# =============================================================================
# Section 17: Software Versions (for reproducibility statement)
# Methods -> "All analyses were conducted using R 4.5.2 interfacing with
#  Python 3.12.12 via reticulate 1.44.1. EBMs were trained with
#  InterpretML 0.7.4 and supporting analyses used scikit-learn 1.8.0
#  and data.table 1.18.0."
# =============================================================================

EXPECTED_VERSIONS <- list(
  R              = "4.5.2",
  python         = "3.12.12",
  reticulate     = "1.44.1",
  interpret      = "0.7.4",
  scikit_learn   = "1.8.0",
  data.table     = "1.18.0"
)

# Optional runtime version check
check_versions <- function(verbose = TRUE) {
  ok <- TRUE
  r_ver <- paste0(R.version$major, ".", R.version$minor)
  if (r_ver != EXPECTED_VERSIONS$R) {
    if (verbose) cat(sprintf("  [WARN] R version: got %s, expected %s\n", r_ver, EXPECTED_VERSIONS$R))
    ok <- FALSE
  }
  dt_ver <- tryCatch(as.character(packageVersion("data.table")), error = function(e) "?")
  if (dt_ver != EXPECTED_VERSIONS$data.table) {
    if (verbose) cat(sprintf("  [WARN] data.table: got %s, expected %s\n", dt_ver, EXPECTED_VERSIONS$data.table))
    ok <- FALSE
  }
  if (ok && verbose) cat("  Version check: OK\n")
  invisible(ok)
}

# =============================================================================
# Print summary
# =============================================================================

cat(sprintf("  OUT_DIR    : %s\n", OUT_DIR))
cat(sprintf("  DATA_DIR   : %s\n", DATA_DIR))
cat(sprintf("  RUN_DIR    : %s\n", RUN_DIR))
cat(sprintf("  Seeds      : %s\n", paste(RANDOM_SEEDS, collapse = ", ")))
cat(sprintf("  E/L split  : %.0f%% / %.0f%% (date col: %s)\n",
            E_SPLIT_QUANTILE * 100, (1 - E_SPLIT_QUANTILE) * 100, SESSION_DATE_COL))
cat(sprintf("  Predictors : %d (%d binary + %d continuous)\n",
            length(ALL_PREDICTORS), length(BINARY_COLS), length(CONT_COLS)))
cat(sprintf("  Tiers      : Green=%d, Yellow=%d, Red=%d, Gray=%d\n",
            length(GREEN_FEATURES), length(YELLOW_FEATURES),
            length(RED_FEATURES), length(GRAY_FEATURES)))
cat(sprintf("  ShapeQC    : tau=%.2f, kappa=%.2f, trim=%s (P%.1f-P%.1f)\n",
            SHAPEQC_TAU, SHAPEQC_KAPPA, PRIMARY_TRIM,
            TRIM_CONFIGS[[PRIMARY_TRIM]][1]*100, TRIM_CONFIGS[[PRIMARY_TRIM]][2]*100))
cat(sprintf("  NRI        : thresholds %s\n",
            paste0(NRI_THRESHOLDS * 100, "%", collapse = " / ")))
cat("\n")