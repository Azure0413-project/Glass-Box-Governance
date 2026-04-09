# =============================================================================
# 01_data_prep.R — Data Loading, Preprocessing & Cohort Description
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Table 1 Panel A — Study population and outcome prevalence
#   Supp Table S7  — Missing data rates per predictor per center
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 01: Data Preparation\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- Output directory ---
OUT_01 <- file.path(RUN_DIR, paste0("01_Data_Prep_", RUN_TS))
dir.create(OUT_01, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Output: %s\n\n", OUT_01))

# =============================================================================
# Step 1: Load raw data from three centers
# =============================================================================

cat("[Step 1] Loading raw data ...\n")

dt_raw <- list()
for (site in names(SITE_FILES)) {
  fp <- SITE_FILES[[site]]
  if (!file.exists(fp)) stop(sprintf("File not found: %s", fp))
  dt_raw[[site]] <- read_fst_dt(fp)
  cat(sprintf("  %s: %s rows, %d cols\n", site, format(nrow(dt_raw[[site]]), big.mark = ","),
              ncol(dt_raw[[site]])))
}

# =============================================================================
# Step 2: Preprocessing — rename, target encoding, physio range filtering
# Methods → "Values outside prespecified physiological ranges were recoded
#  as missing, units were harmonized across sites, and no imputation or
#  scaling was applied"
# =============================================================================

cat("\n[Step 2] Preprocessing ...\n")

dt_sites <- list()
for (site in names(dt_raw)) {
  dt_sites[[site]] <- prep_site_dt(dt_raw[[site]], site)
  cat(sprintf("  %s: %s sessions after filtering (removed %s with NA target)\n",
              site,
              format(nrow(dt_sites[[site]]), big.mark = ","),
              format(nrow(dt_raw[[site]]) - nrow(dt_sites[[site]]), big.mark = ",")))
}

# Combine all sites
dt_all <- rbindlist(dt_sites, use.names = TRUE, fill = TRUE)
cat(sprintf("\n  Combined: %s sessions, %s patients\n",
            format(nrow(dt_all), big.mark = ","),
            format(uniqueN(dt_all[[ID_COL]]), big.mark = ",")))

# =============================================================================
# Step 3: Table 1 Panel A — Study population and outcome prevalence
# =============================================================================

cat("\n[Step 3] Generating Table 1 Panel A ...\n")

table1a <- dt_all[, .(
  Patients    = uniqueN(get(ID_COL)),
  Sessions    = .N,
  IDH_events  = sum(get(TARGET_COL) == 1L, na.rm = TRUE),
  IDH_rate    = 100 * mean(get(TARGET_COL) == 1L, na.rm = TRUE)
), by = site]

# Add overall row
overall <- data.table(
  site       = "Overall",
  Patients   = uniqueN(dt_all[[ID_COL]]),
  Sessions   = nrow(dt_all),
  IDH_events = sum(dt_all[[TARGET_COL]] == 1L, na.rm = TRUE),
  IDH_rate   = 100 * mean(dt_all[[TARGET_COL]] == 1L, na.rm = TRUE)
)
table1a <- rbindlist(list(table1a, overall), use.names = TRUE)

# Add share column
table1a[, Share_pct := round(100 * Sessions / sum(Sessions[site != "Overall"]), 1)]
table1a[site == "Overall", Share_pct := 100.0]

# Format
table1a[, IDH_rate := round(IDH_rate, 2)]

cat("\n  ──────────────────────────────────────────────────────────────\n")
cat("  Table 1 Panel A: Study population and outcome prevalence\n")
cat("  ──────────────────────────────────────────────────────────────\n")
print(table1a)
cat("\n")

# Verify against manuscript values
cat("  [Verify] Manuscript expects:\n")
cat("    Overall: 1,695 patients, 581,469 sessions, 41,786 IDH (7.19%)\n")
cat(sprintf("    Got:     %s patients, %s sessions, %s IDH (%.2f%%)\n",
            format(overall$Patients, big.mark = ","),
            format(overall$Sessions, big.mark = ","),
            format(overall$IDH_events, big.mark = ","),
            overall$IDH_rate))

# Save
fwrite(table1a, file.path(OUT_01, "Table_1A_Cohort.csv"))

# =============================================================================
# Step 4: Supplementary Table S7 — Missing data rates
# Methods → "Missing values include both originally absent data and values
#  recoded as missing by prespecified physiological plausibility filtering"
# =============================================================================

cat("\n[Step 4] Generating Supp Table S7 (missing data rates) ...\n")

# Predictor display order (matches manuscript S7 grouping)
predictor_order <- c(
  "Age", "Sex",
  "Pre_HD_SBP", "Start_DBP", "Heart_Rate", "Respiratory_Rate",
  "IDH_N_7D", "IDH_N_28D",
  "Dry_Weight", "Pre_HD_Weight", "Body_Temperature",
  "Blood_Flow_Rate", "Dialysate_Flow_Rate", "Dialysate_Temperature",
  "Target_UF_Volume", "UF_BW_Perc"
)

# Predictor display names (for publication)
predictor_labels <- c(
  Age = "Age", Sex = "Sex",
  Pre_HD_SBP = "Pre-HD SBP", Start_DBP = "Start DBP",
  Heart_Rate = "Heart rate", Respiratory_Rate = "Respiratory rate",
  IDH_N_7D = "IDH count (7-day)", IDH_N_28D = "IDH count (28-day)",
  Dry_Weight = "Dry weight", Pre_HD_Weight = "Pre-HD weight",
  Body_Temperature = "Body temperature",
  Blood_Flow_Rate = "Blood flow rate", Dialysate_Flow_Rate = "Dialysate flow rate",
  Dialysate_Temperature = "Dialysate temperature",
  Target_UF_Volume = "Target UF volume", UF_BW_Perc = "UF/BW ratio"
)

# Compute missing rates per site
missing_list <- list()
for (site in c("TN", "D6", "CY")) {
  dt_site <- dt_sites[[site]]
  n_site  <- nrow(dt_site)
  for (pred in predictor_order) {
    if (!(pred %in% names(dt_site))) {
      n_miss <- n_site  # column entirely absent
    } else {
      n_miss <- sum(is.na(dt_site[[pred]]))
    }
    pct <- 100 * n_miss / n_site
    missing_list[[length(missing_list) + 1]] <- data.table(
      Predictor  = pred,
      Site       = site,
      n_sessions = n_site,
      n_missing  = n_miss,
      pct_missing = pct
    )
  }
}
missing_dt <- rbindlist(missing_list)

# Wide format
missing_wide <- dcast(missing_dt, Predictor ~ Site,
                      value.var = c("n_missing", "pct_missing"))

# Overall missing
n_total <- nrow(dt_all)
missing_wide[, n_missing_Overall := n_missing_TN + n_missing_D6 + n_missing_CY]
missing_wide[, pct_missing_Overall := round(100 * n_missing_Overall / n_total, 1)]

# Round percentages
for (col in grep("^pct_", names(missing_wide), value = TRUE)) {
  missing_wide[, (col) := round(get(col), 1)]
}

# Add display labels and order
missing_wide[, Label := predictor_labels[Predictor]]
missing_wide[, Predictor := factor(Predictor, levels = predictor_order)]
setorder(missing_wide, Predictor)

cat("\n  Supp Table S7: Missing data rates\n")
print(missing_wide[, .(Label, n_missing_TN, pct_missing_TN,
                        n_missing_D6, pct_missing_D6,
                        n_missing_CY, pct_missing_CY,
                        pct_missing_Overall)])
cat("\n")

# Save
fwrite(missing_wide, file.path(OUT_01, "Supp_S7_Missing_Data.csv"))

# =============================================================================
# Step 5: Descriptive statistics per site (Mean ± SD, Median [IQR])
# =============================================================================

cat("[Step 5] Descriptive statistics per site ...\n")

desc_list <- list()
for (site in c("TN", "D6", "CY")) {
  dt_site <- dt_sites[[site]]
  for (pred in CONT_COLS) {
    if (!(pred %in% names(dt_site))) next
    x <- dt_site[[pred]]
    x <- x[!is.na(x)]
    desc_list[[length(desc_list) + 1]] <- data.table(
      Site     = site,
      Variable = pred,
      N        = length(x),
      N_NA     = sum(is.na(dt_site[[pred]])),
      Mean     = mean(x),
      SD       = sd(x),
      Median   = median(x),
      Q1       = quantile(x, 0.25),
      Q3       = quantile(x, 0.75),
      Min      = min(x),
      Max      = max(x)
    )
  }
  # Sex (binary)
  if ("Sex" %in% names(dt_site)) {
    sex_vals <- dt_site[["Sex"]]
    sex_vals <- sex_vals[!is.na(sex_vals)]
    desc_list[[length(desc_list) + 1]] <- data.table(
      Site = site, Variable = "Sex",
      N = length(sex_vals), N_NA = sum(is.na(dt_site[["Sex"]])),
      Mean = mean(sex_vals), SD = sd(sex_vals),
      Median = median(sex_vals), Q1 = NA_real_, Q3 = NA_real_,
      Min = min(sex_vals), Max = max(sex_vals)
    )
  }
}
desc_stats <- rbindlist(desc_list)

# Save
fwrite(desc_stats, file.path(OUT_01, "Descriptive_Statistics.csv"))
cat(sprintf("  Saved: Descriptive_Statistics.csv (%d rows)\n", nrow(desc_stats)))

# =============================================================================
# Step 6: Save preprocessed site data for downstream modules
# =============================================================================

cat("\n[Step 6] Saving preprocessed data ...\n")

for (site in names(dt_sites)) {
  out_path <- file.path(OUT_01, sprintf("%s_preprocessed.fst", site))
  fst::write_fst(dt_sites[[site]], out_path, compress = 50)
  cat(sprintf("  %s: %s → %s\n", site, format(nrow(dt_sites[[site]]), big.mark = ","),
              basename(out_path)))
}

# Also save combined
combined_path <- file.path(OUT_01, "All_sites_preprocessed.fst")
fst::write_fst(dt_all, combined_path, compress = 50)
cat(sprintf("  Combined: %s → %s\n", format(nrow(dt_all), big.mark = ","),
            basename(combined_path)))

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 01 complete.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

output_files <- list.files(OUT_01, full.names = FALSE)
cat(sprintf("  Output directory: %s\n", OUT_01))
for (f in output_files) {
  fsize <- file.size(file.path(OUT_01, f))
  cat(sprintf("    %s (%.1f KB)\n", f, fsize / 1024))
}
cat("\n")
