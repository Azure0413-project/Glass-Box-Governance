# =============================================================================
# run_Eonly_update.R — Master Index for v3.0.2 E-only Alignment Update
# IDH EBM Governance Reproducibility Pipeline
#
# Purpose:
#   With minimal rerun cost, produce all supplementary table CSVs required
#   for the v3.0.2 E-only analysis, collect them into a timestamped folder
#   <UPDATE_ROOT_BASE>/Update_YYYYMMDD_HHMM, and create a zip bundle.
#   The folder structure is designed for direct upload to update the paper's
#   Supp Tables (S3 / S10 / S12 / S13 / S14).
#
# Key clarification:
#   - 05_shapeqc.R is POST-HOC ONLY and does NOT retrain any EBM. It reads
#     RUN_DIR/Iter*/Shapes/feature_data.{xlsx,csv} to compute C_core/J/S_seed
#     and tier assignments. Pre-trained E-only models can be safely reused.
#   - This script provides an EXISTING_RUN_DIR override to point directly
#     at an existing trained run, avoiding any retraining concern.
#
# Pre-run setup (one-time):
#   1. Place the updated R files in the pipeline directory (where 00_config.R
#      lives):
#        - 05_shapeqc.R                    (v3.0.2 — replaces old version)
#        - verify_eonly_trim_alignment.R   (new)
#        - run_Eonly_update.R              (this file)
#   2. (Optional) Rename old `ShapeQC_v3.0_<ts>` directories to
#      `ShapeQC_v3.0_legacy_full_<ts>` to prevent 07/08 regex from
#      matching stale results.
#
# Execution:
#   setwd("/path/to/pipeline/")
#   source("run_Eonly_update.R")
#
# Output (example):
#   <UPDATE_ROOT_BASE>/Update_20260417_1530/      <- collected folder
#   <UPDATE_ROOT_BASE>/Update_20260417_1530.zip   <- uploadable zip
#
# Revision history:
#   2026-04-17 — Fixed file.path double-slash issue on Windows;
#                added join_path() helper and EXISTING_RUN_DIR override.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
})

# =============================================================================
# Section 0 — Configuration
# =============================================================================

# Master switches. Turn these off if you just want to re-bundle existing outputs.
RUN_VERIFY     <- TRUE   # verify_eonly_trim_alignment.R
RUN_SHAPEQC    <- TRUE   # 05_shapeqc.R v3.0.2 — POST-HOC ONLY, no retraining.
                         # Reads shape functions from RUN_DIR/Iter*/Shapes/ ;
                         # does NOT touch EBM training. Takes ~10-30 min
                         # (dominated by C_core bootstrap), not hours.
BUNDLE_DOWNSTREAM <- TRUE  # copy latest S10 / S12 / S13 from their existing dirs

# ---- Output base path for the update bundle ----
# Edit this to your preferred output location. NO trailing slash.
UPDATE_ROOT_BASE <- Sys.getenv("EBM_UPDATE_ROOT", unset = file.path(getwd(), "output"))

# ---- Existing trained E-only run (skip retraining) ----
# If set to a non-empty path, this OVERRIDES the RUN_DIR produced by 00_config.R,
# forcing the pipeline to use models you already trained.
#   - Leave ""  to use RUN_DIR from 00_config.R.
#   - Set to your E-only trained directory to reuse existing shape exports.
# Requirements of the target directory:
#   - Contains Iter{1..6}_External_*/Shapes/feature_data.xlsx (or .csv)
#   - Contains Iter{1..6}_External_*/*_Final_EBM.joblib (needed by 08/09 if rerun)
EXISTING_RUN_DIR <- Sys.getenv("EBM_RUN_DIR", unset = "")

# Path helper: join a base path and a child, collapsing any double slashes
# (fixes the "//" problem that breaks dir.create on some R/Windows combos).
join_path <- function(base, child) {
  base <- sub("[/\\\\]+$", "", base)   # strip trailing / or \
  child <- sub("^[/\\\\]+", "", child) # strip leading  / or \
  out <- paste0(base, "/", child)
  # collapse any residual //
  while (grepl("//", out, fixed = TRUE)) out <- gsub("//", "/", out, fixed = TRUE)
  out
}

# =============================================================================
# Section 1 — Setup & sanity checks
# =============================================================================

# ---- Auto-detect the directory of THIS script (where the new files live) ----
# Same pattern used by run_all.R. Works when the user calls
# `source("/path/to/run_Eonly_update.R")` from any CWD.
SCRIPT_DIR <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) dirname(normalizePath(ofile, winslash = "/"))
  else getwd()
}, error = function(e) getwd())
SCRIPT_DIR <- gsub("\\\\", "/", SCRIPT_DIR)

# Optional manual override: set this BEFORE source("run_Eonly_update.R") if you
# want to point to a different directory containing 05_shapeqc.R (v3.0.2) and
# verify_eonly_trim_alignment.R. e.g.
#   NEW_FILES_DIR <- "/path/to/new_scripts"
#   source("/path/to/new_scripts/run_Eonly_update.R")
if (exists("NEW_FILES_DIR", inherits = TRUE) && nzchar(NEW_FILES_DIR)) {
  SCRIPT_DIR <- gsub("\\\\", "/", NEW_FILES_DIR)
}

if (!exists("SITE_FILES"))      source("00_config.R")
if (!exists("load_site_eonly")) source("00_utils_r.R")

# ---- Honour EXISTING_RUN_DIR override (no retraining) ----
if (nzchar(EXISTING_RUN_DIR)) {
  EXISTING_RUN_DIR <- gsub("\\\\", "/", EXISTING_RUN_DIR)
  if (!dir.exists(EXISTING_RUN_DIR))
    stop(sprintf("EXISTING_RUN_DIR does not exist: %s", EXISTING_RUN_DIR))

  # Verify it's a trained run: must have Iter*_External_* subdirs with Shapes/
  iter_dirs <- list.dirs(EXISTING_RUN_DIR, recursive = FALSE, full.names = TRUE)
  iter_dirs <- iter_dirs[grepl("^Iter[0-9]+_External_", basename(iter_dirs))]
  if (length(iter_dirs) < 6)
    warning(sprintf("EXISTING_RUN_DIR has %d Iter*_External_ subdirs (expected 6).",
                    length(iter_dirs)))
  shape_ok <- vapply(iter_dirs, function(d)
    file.exists(file.path(d, "Shapes", "feature_data.xlsx")) ||
    file.exists(file.path(d, "Shapes", "feature_data.csv")),
    logical(1))
  if (!all(shape_ok))
    stop(sprintf("Shape exports missing in some iter dirs: %s",
                 paste(basename(iter_dirs)[!shape_ok], collapse = ", ")))

  RUN_DIR <- EXISTING_RUN_DIR
  cat(sprintf("  [OVERRIDE] RUN_DIR redirected to existing trained run:\n    %s\n",
              RUN_DIR))
}

TS           <- format(Sys.time(), "%Y%m%d_%H%M")
UPDATE_ROOT  <- join_path(UPDATE_ROOT_BASE, paste0("Update_", TS))

dir_ok <- dir.create(UPDATE_ROOT, recursive = TRUE, showWarnings = FALSE)
if (!dir_ok && !dir.exists(UPDATE_ROOT))
  stop(sprintf("Failed to create update directory: %s\nCheck that %s exists and is writable.",
               UPDATE_ROOT, UPDATE_ROOT_BASE))

# Canonicalise after creation (gives clean path; mustWork=TRUE now safe)
UPDATE_ROOT  <- normalizePath(UPDATE_ROOT, winslash = "/", mustWork = TRUE)

LOG_FILE <- join_path(UPDATE_ROOT, "00_RUN_LOG.txt")
# Touch the log file immediately to fail fast if unwritable
tryCatch(cat("", file = LOG_FILE),
         error = function(e) stop(sprintf("Cannot write LOG_FILE at %s: %s",
                                          LOG_FILE, conditionMessage(e))))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG_FILE, append = TRUE)
}

cat("\n", paste(rep("#", 70), collapse = ""), "\n", sep = "")
cat("# v3.0.2 E-only Alignment — Master Update Script\n")
cat(paste(rep("#", 70), collapse = ""), "\n\n", sep = "")

log_msg("Script dir  : ", SCRIPT_DIR)
log_msg("Update root : ", UPDATE_ROOT)
log_msg("Pipeline dir: ", getwd())
log_msg("RUN_DIR     : ", RUN_DIR)

# ---- Required-files check: look in SCRIPT_DIR first, then CWD ----
required_files <- c("05_shapeqc.R", "verify_eonly_trim_alignment.R")

resolve_file <- function(nm) {
  # Try SCRIPT_DIR first, then CWD
  p1 <- file.path(SCRIPT_DIR, nm)
  p2 <- file.path(getwd(), nm)
  if (file.exists(p1)) return(normalizePath(p1, winslash = "/"))
  if (file.exists(p2)) return(normalizePath(p2, winslash = "/"))
  return(NA_character_)
}
resolved <- setNames(vapply(required_files, resolve_file, character(1)),
                     required_files)
missing  <- required_files[is.na(resolved)]
if (length(missing) > 0) {
  stop(sprintf(
    "Missing required files: %s\n  Expected in: %s\n            or: %s\nPlease place the new files in one of these directories.",
    paste(missing, collapse = ", "), SCRIPT_DIR, getwd()))
}
log_msg("Resolved 05_shapeqc.R            : ", resolved["05_shapeqc.R"])
log_msg("Resolved verify_eonly_trim_*.R   : ", resolved["verify_eonly_trim_alignment.R"])

# Sanity: is 05_shapeqc.R actually the v3.0.2 patched version?
hdr <- readLines(resolved["05_shapeqc.R"], n = 5L)
if (!any(grepl("v3\\.0\\.2", hdr))) {
  warning("05_shapeqc.R header does not contain 'v3.0.2'. Did you replace it with the revised version?")
  log_msg("[WARN] 05_shapeqc.R may still be the legacy version (header check failed).")
}

# =============================================================================
# Section 2 — Run verify (cheap)
# =============================================================================

safe_source <- function(path, label) {
  log_msg("Running ", label, " ...")
  res <- tryCatch({
    t0 <- Sys.time()
    source(path, local = FALSE)
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    log_msg("  OK (", dt, " min)")
    TRUE
  }, error = function(e) {
    log_msg("  FAILED: ", conditionMessage(e))
    FALSE
  })
  res
}

ok_verify <- if (RUN_VERIFY) safe_source(resolved["verify_eonly_trim_alignment.R"],
                                         "verify") else {
  log_msg("[Skip] RUN_VERIFY = FALSE")
  TRUE
}

# =============================================================================
# Section 3 — Run ShapeQC v3.0.2 (heavy)
# =============================================================================

if (RUN_SHAPEQC) {
  # Note: 05 internally creates a fresh `ShapeQC_v3.0_Eonly_<ts>` dir under RUN_DIR.
  ok_shapeqc <- safe_source(resolved["05_shapeqc.R"],
                            "05_shapeqc.R (v3.0.2 E-only)")
} else {
  log_msg("[Skip] RUN_SHAPEQC = FALSE")
  ok_shapeqc <- TRUE
}

# =============================================================================
# Section 4 — Locate latest output directories
# =============================================================================

log_msg("Locating latest output directories ...")

latest_dir_like <- function(root, pattern) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(pattern, basename(dirs))]
  if (length(dirs) == 0) return(NA_character_)
  sort(dirs, decreasing = TRUE)[1]
}

dir_verify  <- latest_dir_like(RUN_DIR, "^Verify_Eonly_Trim_")
dir_shapeqc <- latest_dir_like(RUN_DIR, "^ShapeQC_v3\\.0_Eonly_")
dir_lofo    <- latest_dir_like(RUN_DIR, "^Sensitivity_LOFO_")
dir_zero    <- latest_dir_like(RUN_DIR, "^Zeroing_Comparison_")
dir_nri     <- latest_dir_like(RUN_DIR, "^NRI_Reclassification_")

# For the tier-flip diff, find the MOST RECENT pre-v3.0.2 ShapeQC output.
dir_shapeqc_prev <- {
  all <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
  cand <- all[grepl("^ShapeQC_v3\\.0(r|_)?[0-9_]+$", basename(all))]
  # exclude the new Eonly dir
  cand <- cand[!grepl("_Eonly_", basename(cand))]
  if (length(cand) == 0) NA_character_ else sort(cand, decreasing = TRUE)[1]
}

for (nm in c("dir_verify", "dir_shapeqc", "dir_shapeqc_prev",
             "dir_lofo", "dir_zero", "dir_nri")) {
  val <- get(nm)
  log_msg(sprintf("  %-20s : %s", nm,
                  if (is.na(val)) "-- (none found)" else basename(val)))
}

# =============================================================================
# Section 5 — Harvest files into UPDATE_ROOT
# =============================================================================

harvest <- function(src_dir, dest_subfolder, file_patterns = NULL, source_label = NULL) {
  if (is.na(src_dir) || !dir.exists(src_dir)) {
    log_msg("  [skip] ", dest_subfolder, " -- source directory missing")
    return(FALSE)
  }
  dest <- file.path(UPDATE_ROOT, dest_subfolder)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  # Write a SOURCE.txt for traceability
  writeLines(c(
    sprintf("Source directory: %s", src_dir),
    sprintf("Harvested at    : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    if (!is.null(source_label)) sprintf("Produced by     : %s", source_label) else ""
  ), file.path(dest, "SOURCE.txt"))

  files <- if (is.null(file_patterns)) {
    list.files(src_dir, full.names = TRUE, recursive = FALSE)
  } else {
    unlist(lapply(file_patterns, function(p)
      list.files(src_dir, pattern = p, full.names = TRUE, recursive = FALSE)))
  }
  files <- files[!dir.exists(files)]   # files only, not subdirs
  if (length(files) == 0) {
    log_msg("  [warn] ", dest_subfolder, " -- no files matched")
    return(FALSE)
  }
  file.copy(files, dest, overwrite = TRUE)
  log_msg(sprintf("  %-30s <- %d file(s)", dest_subfolder, length(files)))
  TRUE
}

log_msg("Harvesting ...")

harvest(dir_verify,  "01_Verify",
        file_patterns = "\\.csv$",
        source_label = "verify_eonly_trim_alignment.R")

harvest(dir_shapeqc, "02_Supp_S3_S14_ShapeQC_Eonly",
        file_patterns = c("Table_4A.*\\.csv$",
                          "Supp_S3.*\\.csv$",
                          "Supp_S14.*\\.csv$",
                          "ShapeQC.*Report\\.xlsx$"),
        source_label = "05_shapeqc.R v3.0.2 (E-only)")

if (BUNDLE_DOWNSTREAM) {
  harvest(dir_lofo, "03_Supp_S10_LOFO",
          file_patterns = c("Supp_S10.*\\.csv$",
                            "LOFO.*\\.csv$",
                            "Supp_S11.*\\.csv$",
                            "Supp_S8.*\\.csv$",
                            "Supp_S9.*\\.csv$"),
          source_label = "08_sensitivity_lofo.R (E-only filter already applied)")

  harvest(dir_zero,  "04_Supp_S12_Zeroing",
          file_patterns = c("Supp_S12.*\\.csv$",
                            "Supp_S15.*\\.csv$",
                            "All_Zeroing.*\\.csv$"),
          source_label = "09_zeroing_comparison.R (E-only filter already applied)")

  harvest(dir_nri,   "05_Supp_S13_NRI",
          file_patterns = "Table_S13.*\\.csv$",
          source_label = "10_reclassification.R (E-only filter already applied)")
}

# =============================================================================
# Section 6 — Tier-flip detection (diff new Table_4A vs previous)
# =============================================================================

log_msg("Checking for tier flips ...")

tier_flip_report <- data.table()
if (!is.na(dir_shapeqc) && !is.na(dir_shapeqc_prev)) {
  f_new  <- file.path(dir_shapeqc,      "Table_4A_Decision_Matrix.csv")
  f_prev <- file.path(dir_shapeqc_prev, "Table_4A_Decision_Matrix.csv")
  if (file.exists(f_new) && file.exists(f_prev)) {
    dn <- fread(f_new);  if ("light" %in% names(dn)) setnames(dn, "light", "tier")
    dp <- fread(f_prev); if ("light" %in% names(dp)) setnames(dp, "light", "tier")

    key_cols <- intersect(c("feature", "tier", "C0_trim", "median_J"), names(dn))
    key_cols_p <- intersect(c("feature", "tier", "C0_trim", "median_J"), names(dp))
    m <- merge(dn[, ..key_cols], dp[, ..key_cols_p],
               by = "feature", suffixes = c("_new", "_prev"), all = TRUE)
    m[, tier_flip := tier_new != tier_prev]
    m[, delta_C0  := round(C0_trim_new  - C0_trim_prev,  4)]
    m[, delta_J   := round(median_J_new - median_J_prev, 3)]
    tier_flip_report <- m
    fwrite(m, file.path(UPDATE_ROOT, "06_TierFlip_Diff.csv"))

    n_flips <- sum(m$tier_flip, na.rm = TRUE)
    if (n_flips == 0) {
      log_msg("  [OK] No tier flips -- S10 / S12 / S13 remain valid without rerun.")
    } else {
      log_msg("  [ACTION] ", n_flips, " feature(s) flipped tier -- rerun 07/08/09/10 recommended.")
      print(m[tier_flip == TRUE])
    }
  } else {
    log_msg("  [skip] Cannot diff: one of the Table_4A files is missing.")
  }
} else {
  log_msg("  [skip] No prior ShapeQC output found -- diff not possible.")
}

# =============================================================================
# Section 7 — Write README.md / MANIFEST.csv
# =============================================================================

log_msg("Writing README.md and MANIFEST.csv ...")

manifest <- data.table(
  file = character(), supp_table = character(),
  produced_by = character(), frame = character(), notes = character()
)
add_mf <- function(file, supp, by, frame, notes = "") {
  manifest <<- rbind(manifest, data.table(
    file = file, supp_table = supp, produced_by = by, frame = frame, notes = notes
  ))
}

# Build manifest programmatically from what was actually harvested
for (sub in list.dirs(UPDATE_ROOT, recursive = FALSE)) {
  for (f in list.files(sub, pattern = "\\.csv$|\\.xlsx$", full.names = FALSE)) {
    if (f == "SOURCE.txt") next
    supp <- {
      # Heuristic mapping from filename to Supp table
      if      (grepl("S3_", f))  "Supp_S3"
      else if (grepl("S10_", f)) "Supp_S10"
      else if (grepl("S12_", f)) "Supp_S12"
      else if (grepl("S13_", f)) "Supp_S13"
      else if (grepl("S14_", f)) "Supp_S14"
      else if (grepl("S11_", f)) "Supp_S11"
      else if (grepl("S15_", f)) "Supp_S15"
      else if (grepl("S8_", f))  "Supp_S8"
      else if (grepl("S9_", f))  "Supp_S9"
      else if (grepl("Table_4A", f)) "Table_4A_main"
      else if (grepl("Trim_Bounds", f)) "diagnostic"
      else if (grepl("TierFlip", f))    "diagnostic"
      else "--"
    }
    add_mf(file = file.path(basename(sub), f),
           supp = supp,
           by   = basename(sub),
           frame = "E-only (v3.0.2)",
           notes = "")
  }
}

fwrite(manifest, file.path(UPDATE_ROOT, "00_MANIFEST.csv"))

readme_lines <- c(
  "# v3.0.2 E-only Alignment — Update Bundle",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Pipeline : %s", getwd()),
  sprintf("Source RUN_DIR: %s", RUN_DIR),
  "",
  "## Purpose",
  "",
  "This bundle contains the CSV / XLSX files needed to update the paper's",
  "Supplementary Tables (S3 / S10 / S12 / S13 / S14) after the v3.0.2",
  "E-only analytic-frame alignment.",
  "",
  "## Bundle Layout",
  "",
  "| Folder | Contains | Supp Table |",
  "|--------|----------|------------|",
  "| 01_Verify/ | Trim-bound E+L vs E-only comparison (sanity check) | -- |",
  "| 02_Supp_S3_S14_ShapeQC_Eonly/ | New ShapeQC v3.0.2 outputs | **S3, S14, Table 4A** |",
  "| 03_Supp_S10_LOFO/ | LOFO + smoothing sensitivity (already E-only) | **S10** (+S8, S9, S11) |",
  "| 04_Supp_S12_Zeroing/ | Zeroing vs smoothing (already E-only) | **S12** (+S15) |",
  "| 05_Supp_S13_NRI/ | NRI / IDI panels (already E-only) | **S13 A/B/C** |",
  "| 06_TierFlip_Diff.csv | Diff of new vs previous Table 4A | diagnostic |",
  "| 00_MANIFEST.csv | Machine-readable file to Supp map | -- |",
  "| 00_RUN_LOG.txt | Execution log | -- |",
  "",
  "## Tier-Flip Status",
  "",
  if (nrow(tier_flip_report) == 0)
    "No previous ShapeQC output was available for diff (first v3.0.2 run)."
  else {
    nf <- sum(tier_flip_report$tier_flip, na.rm = TRUE)
    if (nf == 0)
      "**No tier flips.** Existing S10 / S12 / S13 outputs remain numerically valid."
    else
      sprintf("**%d tier flip(s) detected.** Rerun 07 -> 08 -> 09 -> 10 before finalising Supp tables. See 06_TierFlip_Diff.csv.", nf)
  },
  ""
)
writeLines(readme_lines, file.path(UPDATE_ROOT, "00_README.md"))

# =============================================================================
# Section 8 — Zip the bundle
# =============================================================================

log_msg("Zipping ...")

zip_path <- paste0(UPDATE_ROOT, ".zip")
zip_ok <- FALSE

# Preferred: zip::zip (cross-platform, no external dependency)
if (requireNamespace("zip", quietly = TRUE)) {
  tryCatch({
    zip::zip(zipfile = zip_path,
             files   = list.files(UPDATE_ROOT, recursive = TRUE, full.names = FALSE),
             root    = UPDATE_ROOT)
    zip_ok <- TRUE
    log_msg("  [zip::zip] Created: ", zip_path)
  }, error = function(e) {
    log_msg("  [zip::zip] failed: ", conditionMessage(e))
  })
}

# Fallback: utils::zip (requires system zip on Windows)
if (!zip_ok) {
  tryCatch({
    old_wd <- getwd()
    setwd(dirname(UPDATE_ROOT))
    on.exit(setwd(old_wd), add = TRUE)
    utils::zip(zipfile = basename(zip_path),
               files   = basename(UPDATE_ROOT),
               flags   = "-r9Xq")
    zip_ok <- file.exists(zip_path)
    if (zip_ok) log_msg("  [utils::zip] Created: ", zip_path)
  }, error = function(e) {
    log_msg("  [utils::zip] failed: ", conditionMessage(e),
            " -- folder is still available at ", UPDATE_ROOT)
  })
}

# =============================================================================
# Section 9 — Final summary
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("  Update bundle complete\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

cat(sprintf("  Folder : %s\n", UPDATE_ROOT))
if (zip_ok) cat(sprintf("  Zip    : %s\n", zip_path))
cat(sprintf("  Files  : %d\n",
            length(list.files(UPDATE_ROOT, recursive = TRUE))))

if (nrow(tier_flip_report) > 0) {
  nf <- sum(tier_flip_report$tier_flip, na.rm = TRUE)
  if (nf == 0) {
    cat("  Status : No tier flips -- upload zip for Supp update.\n")
  } else {
    cat(sprintf("  Status : %d tier flip(s) -- rerun downstream before finalising.\n", nf))
  }
} else {
  cat("  Status : First v3.0.2 run (no diff baseline).\n")
}

cat("\n")
