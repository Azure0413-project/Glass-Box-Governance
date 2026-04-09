# =============================================================================
# 06_tv_decomposition.R — Total-Variation Tail vs Core Decomposition
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Supp Table S5 — TV decomposition: tail- vs core-localized irregularity
#   Fig (Supp)    — TV retention curve, heatmap, stacked bar
#
# Source: consolidated from original analysis scripts
# Design: For each feature × IECV iteration, decompose total variation (TV)
#       into tail-localized vs core-localized contributions using weighted
#       quantile trimming at 5 prespecified thresholds.
#
# Methods → "Total-variation decomposition revealed that Yellow-tier
#  irregularity was predominantly tail-localized—the outer 5% tails
#  contained up to 63% of total variation"
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 05_shapeqc.R (needs ShapeQC projected shapes)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(ggplot2)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 06: TV Tail vs Core Decomposition\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# Configuration
# =============================================================================

# Trim thresholds (each side percentile)
TRIM_PCTS <- c(0.01, 0.05, 0.10, 0.20, 0.40)

# Output directory
FIG_OUT_DIR <- file.path(RUN_DIR, paste0("TV_TailAnalysis_", RUN_TS))
dir.create(FIG_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  Trim thresholds: %s\n", paste0(TRIM_PCTS * 100, "%", collapse = ", ")))
cat(sprintf("  Output: %s\n\n", basename(FIG_OUT_DIR)))

# =============================================================================
# Load shape data from ShapeQC
# =============================================================================

cat("[Step 1] Loading shape data ...\n")

# Auto-detect ShapeQC output
qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
qc_dirs <- qc_dirs[grepl("ShapeQC_v3\\.0_", basename(qc_dirs))]
if (length(qc_dirs) == 0) stop("ShapeQC output not found. Run 05_shapeqc.R first.")
qc_dir <- sort(qc_dirs, decreasing = TRUE)[1]

# Load projected shapes from Excel or try individual feature_data files
xlsx_path <- file.path(qc_dir, "ShapeQC_v3.0_Report.xlsx")
if (file.exists(xlsx_path)) {
  proj_all <- as.data.table(openxlsx::read.xlsx(xlsx_path, sheet = "Projected_Shapes"))
  cat(sprintf("  Loaded from ShapeQC report: %d rows\n", nrow(proj_all)))
} else {
  # Fallback: load from iteration Shapes directories
  iter_tags <- c("Iter1_External_TN_Seed1", "Iter2_External_TN_Seed2",
                  "Iter3_External_D6_Seed1", "Iter4_External_D6_Seed2",
                  "Iter5_External_CY_Seed1", "Iter6_External_CY_Seed2")
  proj_all <- rbindlist(lapply(iter_tags, function(tag) {
    fp <- file.path(RUN_DIR, tag, "Shapes", "feature_data.csv")
    if (!file.exists(fp)) fp <- file.path(RUN_DIR, tag, "Shapes", "feature_data.xlsx")
    if (!file.exists(fp)) return(data.table())
    dt <- if (grepl("\\.xlsx$", fp)) as.data.table(openxlsx::read.xlsx(fp)) else fread(fp)
    dt[, model := tag]
    dt
  }), fill = TRUE)
  cat(sprintf("  Loaded from individual Shapes: %d rows\n", nrow(proj_all)))
}

# Ensure columns
# proj_all should have: feature, model, x (or x_g), score (or f_aligned)
if ("x_g" %in% names(proj_all) && !("x" %in% names(proj_all))) setnames(proj_all, "x_g", "x")
if ("f_aligned" %in% names(proj_all) && !("score" %in% names(proj_all))) setnames(proj_all, "f_aligned", "score")

# Load site data for weighted quantiles
site_data_cache <- list()
for (s in names(SITE_FILES)) {
  fp <- SITE_FILES[[s]]
  if (!file.exists(fp)) next
  dt_raw <- read_fst_dt(fp)
  dt_raw <- apply_compat_rename(dt_raw, COMPAT_RENAME_MAP)
  site_data_cache[[s]] <- dt_raw
}

# Build weights from training data
# Each iteration (model) was trained on 2 sites; use those sites' data as weights
model_train_sites <- list(
  Iter1_External_TN_Seed1 = c("D6", "CY"), Iter2_External_TN_Seed2 = c("D6", "CY"),
  Iter3_External_D6_Seed1 = c("TN", "CY"), Iter4_External_D6_Seed2 = c("TN", "CY"),
  Iter5_External_CY_Seed1 = c("TN", "D6"), Iter6_External_CY_Seed2 = c("TN", "D6")
)

# =============================================================================
# Weighted Quantile Function
# =============================================================================

weighted_quantile <- function(x, w, probs) {
  valid <- is.finite(x) & is.finite(w) & w > 0
  x <- x[valid]; w <- w[valid]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  sapply(probs, function(p) {
    idx <- which(cum_w >= p)[1]
    if (is.na(idx)) return(x[length(x)])
    x[idx]
  })
}

# =============================================================================
# Main Analysis: TV decomposition by trim threshold
# =============================================================================

cat("\n[Step 2] Computing TV decomposition ...\n")

features <- setdiff(unique(proj_all$feature), c("Sex", ""))
models   <- unique(proj_all$model)

res_list <- list()
for (feat in features) {
  for (mdl in models) {
    sub <- proj_all[feature == feat & model == mdl & is.finite(x) & is.finite(score)]
    if (nrow(sub) < 10) next

    # Determine x and score columns
    x_val <- sub$x
    s_val <- sub$score

    # Get weights from training data
    train_sites <- model_train_sites[[mdl]]
    w_vals <- numeric(0)
    x_train <- numeric(0)
    for (ts in train_sites) {
      if (!is.null(site_data_cache[[ts]]) && feat %in% names(site_data_cache[[ts]])) {
        v <- as.numeric(site_data_cache[[ts]][[feat]])
        v <- v[!is.na(v)]
        x_train <- c(x_train, v)
      }
    }

    # Create weights by binning training data onto shape grid
    # Simple approach: use uniform weights if can't match
    w_val <- rep(1, nrow(sub))
    if (length(x_train) > 0) {
      h <- hist(x_train, breaks = c(-Inf, head(sub$x, -1) + diff(sub$x)/2, Inf), plot = FALSE)
      if (length(h$counts) == nrow(sub)) w_val <- h$counts + 1  # +1 smoothing
    }

    # Total TV
    TV_total <- sum(abs(diff(s_val)))
    if (TV_total < 1e-12) next

    # For each trim threshold, compute core/tail TV
    for (p in TRIM_PCTS) {
      q_lo <- weighted_quantile(x_val, w_val, p)
      q_hi <- weighted_quantile(x_val, w_val, 1 - p)

      core_mask <- x_val >= q_lo & x_val <= q_hi
      if (sum(core_mask) < 3) next

      # Core TV: TV within [q_lo, q_hi]
      s_core <- s_val[core_mask]
      TV_core <- sum(abs(diff(s_core)))

      # Tail TV = Total - Core
      TV_tail <- TV_total - TV_core

      res_list[[length(res_list) + 1]] <- data.table(
        feature   = feat,
        model     = mdl,
        trim_pct  = p,
        TV_total  = TV_total,
        TV_core   = TV_core,
        TV_tail   = TV_tail,
        core_frac = TV_core / TV_total,
        tail_frac = TV_tail / TV_total
      )
    }
  }
}

res_dt <- rbindlist(res_list)
cat(sprintf("  Done: %d rows (%d features × %d iterations × %d thresholds)\n",
            nrow(res_dt), length(features), length(models), length(TRIM_PCTS)))

# =============================================================================
# Summary statistics
# =============================================================================

res_summary <- res_dt[, .(
  core_frac_mean = mean(core_frac, na.rm = TRUE),
  core_frac_lo   = quantile(core_frac, 0.025, na.rm = TRUE),
  core_frac_hi   = quantile(core_frac, 0.975, na.rm = TRUE),
  tail_frac_mean = mean(tail_frac, na.rm = TRUE),
  tail_frac_lo   = quantile(tail_frac, 0.025, na.rm = TRUE),
  tail_frac_hi   = quantile(tail_frac, 0.975, na.rm = TRUE),
  TV_total_mean  = mean(TV_total, na.rm = TRUE)
), by = .(feature, trim_pct)]

fwrite(res_dt,      file.path(FIG_OUT_DIR, "TV_tail_analysis_raw.csv"))
fwrite(res_summary, file.path(FIG_OUT_DIR, "Supp_S5_TV_Summary.csv"))

# =============================================================================
# Figure 1: TV Retention Curve
# =============================================================================

cat("\n[Step 3] Generating figures ...\n")

# Add anchor at trim=0 (100% retained)
anchor_dt <- res_summary[, .(trim_pct = 0,
                              core_frac_mean = 1, core_frac_lo = 1, core_frac_hi = 1,
                              tail_frac_mean = 0, tail_frac_lo = 0, tail_frac_hi = 0,
                              TV_total_mean = TV_total_mean[1]),
                          by = feature]
ret_data <- rbind(anchor_dt, res_summary, fill = TRUE)[order(feature, trim_pct)]

# Classify archetypes for color-coding
ret_data[, archetype := fifelse(
  feature %in% c("Blood_Flow_Rate", "Heart_Rate", "Start_DBP"), "Yellow (jagged, tail-dominated)",
  fifelse(feature %in% RED_FEATURES, "Red (non-transportable)",
  fifelse(feature %in% GREEN_FEATURES, "Green (transportable)", "Other"))
)]

p1 <- ggplot(ret_data, aes(x = trim_pct * 100, y = core_frac_mean * 100, color = archetype)) +
  geom_line(aes(group = feature), linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_color_manual(values = c("Green (transportable)" = "#27ae60",
                                 "Yellow (jagged, tail-dominated)" = "#f1c40f",
                                 "Red (non-transportable)" = "#e74c3c",
                                 "Other" = "gray60")) +
  labs(title = "TV Retention: How much TV remains after trimming data-sparse tails?",
       x = "Tail trimmed (each side, % of weighted sample)",
       y = "Core TV retained (%)",
       color = "Governance tier") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

ggsave(file.path(FIG_OUT_DIR, "Fig_TV_Retention_Curve.png"), p1, width = 10, height = 7, dpi = 300)

# =============================================================================
# Figure 2: Heatmap — features × trim thresholds (tail fraction)
# =============================================================================

heat_data <- res_summary[, .(feature, trim_pct, tail_frac_mean)]
heat_data[, trim_label := factor(paste0(trim_pct * 100, "%"),
                                  levels = paste0(sort(TRIM_PCTS) * 100, "%"))]

# Order features by TV_total
feat_order <- res_summary[trim_pct == TRIM_PCTS[1], .(TV_total_mean = mean(TV_total_mean)),
                           by = feature][order(-TV_total_mean)]$feature
heat_data[, feature := factor(feature, levels = rev(feat_order))]

p2 <- ggplot(heat_data, aes(x = trim_label, y = feature, fill = tail_frac_mean * 100)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", tail_frac_mean * 100)), size = 3) +
  scale_fill_gradient2(low = "#2ecc71", mid = "#f1c40f", high = "#e74c3c",
                       midpoint = 30, name = "% TV in tails",
                       limits = c(0, 100)) +
  labs(title = "TV concentration in data-sparse tails",
       subtitle = "Cell = % of total TV outside the central [p, 1-p] weighted sample range",
       x = "Tail trim threshold (each side)",
       y = "Feature") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank())

ggsave(file.path(FIG_OUT_DIR, "Fig_TV_Tail_Heatmap.png"), p2, width = 9, height = 8, dpi = 300)
cat("  Figures saved.\n")

# =============================================================================
# Console summary
# =============================================================================

cat("\n  ─── Supp S5: TV Tail vs Core Summary ───\n")
for (feat in c("Blood_Flow_Rate", "Heart_Rate", "Start_DBP",
                "Dialysate_Flow_Rate", "Dry_Weight", "Pre_HD_SBP")) {
  sub5 <- res_summary[feature == feat & trim_pct == 0.05]
  if (nrow(sub5) > 0) {
    cat(sprintf("  %s: %.0f%% tail TV at 5%% trim (tier: %s)\n",
                feat, sub5$tail_frac_mean * 100,
                ifelse(feat %in% names(FEATURE_TIER_MAP), FEATURE_TIER_MAP[feat], "?")))
  }
}

cat("\n  [Verify] Manuscript states Yellow-tier: outer 5% tails contain up to 63% TV\n")

cat(sprintf("\n  Output: %s\n", FIG_OUT_DIR))
cat("\n  Module 06 complete.\n\n")
