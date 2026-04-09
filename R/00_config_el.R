# =============================================================================
# 00_config_el.R — E/L Temporal Split Configuration Supplement
# IDH EBM Governance — Plan A: E/L Temporal Split with Protocol Lock
#
# Loaded after 00_config.R to add E/L temporal split and confirmatory analysis settings.
# Usage: src("R/00_config.R"); source("00_config_el.R")
# =============================================================================

# --- Auto-load prerequisite dependencies ---
if (!exists("RUN_DIR"))    src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  E/L Temporal Split Configuration Loaded\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# Section EL-1: Temporal Split Parameters
# Plan A Step 1 — center-specific 75/25 temporal split (matching S6)
# =============================================================================

EL_TEMPORAL_RATIO <- 0.75   # 75% Earlier / 25% Later (matches S6)

# =============================================================================
# Section EL-2: Hyperparameter Source
# Plan A Step 2 — Pure E-only design (v3.0)
# Methods → "Hyperparameters were tuned within the E-cohort IECV pipeline
#  (03_iecv_ebm.R v3.0), ensuring complete protocol lock — L-cohort data
#  do not influence any model decision, including hyperparameters."
#
# Under this design, 04d_el_derivation_iecv.R is REDUNDANT because 03
# already performs E-only tuning, training, and validation.
# The confirmatory pipeline (04f) reads models directly from RUN_DIR.
# =============================================================================

EL_SKIP_GRID_SEARCH <- FALSE  # Legacy flag — no longer used.
                               # 03 v3.0 tunes on E-only by design.

# E-only IECV artifacts source = same as main RUN_DIR
EL_EONLY_IECV_DIR <- RUN_DIR

# =============================================================================
# Section EL-2b: Memory Management for E-only Training
# Original pipeline (03_iecv_ebm.R) used parallel mclapply across iterations,
# so each R worker had its own memory space. The EL pipeline runs sequentially,
# accumulating memory. Reduce EBM internal parallelism to avoid OOM kills.
# =============================================================================

EL_EBM_N_JOBS       <- 3L     # Original used 5 with max_bins=64;
                               # we use max_bins=128, so 3 is conservative.
                               # Safe with PSOCK process isolation.
EL_OUTER_BAGS_OOF   <- 10L    # Same as full-IECV (PSOCK = clean memory)
EL_OUTER_BAGS_FINAL <- 64L    # Matched to full-IECV for consistency.

# =============================================================================
# Section EL-3: Governance Policy Candidates
# Plan A Step 3 — candidate governance policies (limited to 2)
# Policy A = primary; Policy B = sensitivity analysis
# =============================================================================

EL_POLICIES <- list(
  A = list(
    id       = "A",
    label    = "Smooth Yellow + Zero Red (Primary)",
    yellow   = "smooth",
    red      = "zero",
    primary  = TRUE
  ),
  B = list(
    id       = "B",
    label    = "Smooth All Eligible (Sensitivity)",
    yellow   = "smooth",
    red      = "smooth",
    primary  = FALSE
  )
)

# =============================================================================
# Section EL-4: Gold Standard Smoothing Parameters
# Ported from Post-hoc v2.1 gold standard
# =============================================================================

EL_UPDATE_ONLY_PHYS_RANGE <- TRUE
EL_UPDATE_MISSING_BIN     <- FALSE
EL_UPDATE_UNKNOWN_BIN     <- FALSE
EL_CLIP_TO_ORIGINAL_RANGE <- TRUE
EL_CLIP_MARGIN            <- 0.25
EL_SPLINE_MIN_UNIQUE_X    <- 5L

# =============================================================================
# Section EL-5: Protocol Lock — Success Criteria
# Plan A Step 4 — "Primary confirmatory endpoints vs. secondary"
# =============================================================================

EL_PRIMARY_SEED <- 1L   # seed_id=1 (master_seed=2024) is primary
                         # seed_id=2 (master_seed=9999) is sensitivity

# Material harm definition: ΔAUPRC < -0.01 at any single center
EL_HARM_MARGIN <- -0.01

# Non-inferiority: center-unweighted mean ΔAUPRC ≥ 0
EL_NI_MARGIN   <- 0.0

# Bootstrap
EL_BOOT_N    <- 1000L
EL_BOOT_SEED <- 2024L

# =============================================================================
# Section EL-6: Morphological Agreement Threshold
# Plan A Step 3 — if agreement >= 0.95 and tiers unchanged, retain expert review
# =============================================================================

EL_AGREEMENT_THRESHOLD <- 0.95  # min Spearman rho per feature

# =============================================================================
# Section EL-6b: Tier Source for Protocol Lock
# "eonly" = E-only tiers from 04e (clean derivation–confirmation separation)
# "full"  = full-data tiers from 00_config.R (legacy, not recommended)
#
# Rationale for "eonly" (default):
#   Full-data tiers were informed by S8 scenario analysis on external data
#   that includes L-period sessions. Using them in confirmatory scoring
#   does not fully eliminate the governance-derivation–verification circularity
#   that Plan A was designed to address.
# =============================================================================

EL_TIER_SOURCE <- "eonly"  # "eonly" or "full"

# =============================================================================
# Section EL-7: Output Directories
# =============================================================================

EL_OUT_ROOT <- file.path(OUT_DIR, "EL_Confirmatory")
EL_DIR_SPLIT     <- file.path(EL_OUT_ROOT, "01_Split_Table")
EL_DIR_DERIVATION <- file.path(EL_OUT_ROOT, "02_E_Only_IECV")
EL_DIR_SHAPEQC   <- file.path(EL_OUT_ROOT, "03_E_Only_ShapeQC")
EL_DIR_LOCK      <- file.path(EL_OUT_ROOT, "04_Protocol_Lock")
EL_DIR_CONFIRM   <- file.path(EL_OUT_ROOT, "05_L_Only_Confirmatory")
EL_DIR_REPORT    <- file.path(EL_OUT_ROOT, "06_Report")

lapply(c(EL_DIR_SPLIT, EL_DIR_DERIVATION, EL_DIR_SHAPEQC,
         EL_DIR_LOCK, EL_DIR_CONFIRM, EL_DIR_REPORT),
       dir.create, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Print summary
# =============================================================================

cat(sprintf("  Temporal split : %.0f%% E / %.0f%% L\n",
            EL_TEMPORAL_RATIO * 100, (1 - EL_TEMPORAL_RATIO) * 100))
cat(sprintf("  Skip grid search: %s\n", EL_SKIP_GRID_SEARCH))
cat(sprintf("  EBM training   : n_jobs=%d, oof_bags=%d, final_bags=%d\n",
            EL_EBM_N_JOBS, EL_OUTER_BAGS_OOF, EL_OUTER_BAGS_FINAL))
cat(sprintf("  Policies       : %d (A=primary, B=sensitivity)\n", length(EL_POLICIES)))
cat(sprintf("  Harm margin    : ΔAUPRC < %.3f\n", EL_HARM_MARGIN))
cat(sprintf("  Agreement thr  : Spearman ρ ≥ %.2f\n", EL_AGREEMENT_THRESHOLD))
cat(sprintf("  Output root    : %s\n", EL_OUT_ROOT))
cat("\n")
