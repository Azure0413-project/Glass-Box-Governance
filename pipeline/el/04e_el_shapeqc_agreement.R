# =============================================================================
# 04e_el_shapeqc_agreement.R — DEPRECATED (v3.0 Pure E-only Design)
# IDH EBM Governance — Plan A Phase 3
#
# ⚠ THIS SCRIPT IS NO LONGER NEEDED IN ITS ORIGINAL FORM.
#
# Under the pure E-only design (v3.0):
#   - 03_iecv_ebm.R produces E-only models directly
#   - 05_shapeqc.R runs on these E-only models → produces E-only tier
#     assignments (Table_4A_Decision_Matrix.csv)
#   - There is no separate "full-data" run, so the E-only vs full-data
#     morphological agreement comparison is no longer applicable.
#
# Outputs previously produced by this script:
#   - EL_Eonly_QC_Tiers.csv  → Now use 05_shapeqc.R Table_4A_Decision_Matrix.csv
#   - EL_Eonly_SPAR_Map.csv  → Derived on-the-fly by 04f from 05's output + SPAR_FORMULA
#
# If you need to compare E-only shapes against a HISTORICAL full-cohort run
# (e.g., from a pre-v3.0 pipeline), retain the original 04e script and point
# EL_DIR_DERIVATION to the historical run directory.
#
# History:
#   v2.x — E-only vs full-data shape agreement + tier derivation
#   v3.0 — Superseded by 05_shapeqc.R on E-only models
# =============================================================================

if (!exists("SITE_FILES")) src("R/00_config.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 04e: DEPRECATED — Superseded by 05_shapeqc.R on E-only models\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("  Under the pure E-only design (v3.0), 05_shapeqc.R now produces\n")
cat("  E-only tier assignments directly from 03's E-only models.\n")
cat("  The E-only vs full-data agreement comparison is no longer applicable.\n\n")

# Check if 05's output exists
qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
qc_dirs <- qc_dirs[grepl("ShapeQC", basename(qc_dirs))]

if (length(qc_dirs) > 0) {
  qc_csv <- file.path(qc_dirs[length(qc_dirs)], "Table_4A_Decision_Matrix.csv")
  if (file.exists(qc_csv)) {
    qc_dt <- data.table::fread(qc_csv)
    tier_col <- if ("light" %in% names(qc_dt)) "light" else "tier"
    cat(sprintf("  05_shapeqc.R output found: %s\n", basename(qc_dirs[length(qc_dirs)])))
    cat(sprintf("    Green=%d, Yellow=%d, Red=%d, Gray=%d\n",
                sum(qc_dt[[tier_col]] == "Green"), sum(qc_dt[[tier_col]] == "Yellow"),
                sum(qc_dt[[tier_col]] == "Red"),   sum(qc_dt[[tier_col]] == "Gray")))
    cat("\n  → 04f can read tiers directly from 05's output.\n")
  } else {
    cat("  [WARNING] Table_4A not found in ShapeQC output. Run 05_shapeqc.R.\n")
  }
} else {
  cat("  [WARNING] No ShapeQC output found. Run 05_shapeqc.R on E-only models.\n")
}

cat("\n  Module 04e (deprecated shim) complete.\n\n")
