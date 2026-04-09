# =============================================================================
# 02_heterogeneity.R — Cross-site Distributional Heterogeneity
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Table 1 Panel B — Cliff's δ effect sizes for all pairwise comparisons
#   Supp Table 1   — Kruskal-Wallis tests for between-center differences
#   Supp Table 2   — KS tests for pairwise between-center differences
#   Fig 2a         — Cliff's δ heatmap
#   Fig 2b         — Density overlays for key predictors
#
# Source: consolidated from original analysis scripts
# All sampling parameters are retained at original values (KW 100k, KS 100k, Cliff 50k)
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
#       Requires 01_data_prep.R (or raw .fst data files)
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
cat("  Module 02: Cross-site Heterogeneity\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- Output directory ---
OUT_02 <- file.path(RUN_DIR, paste0("02_Heterogeneity_", RUN_TS))
dir.create(OUT_02, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Output: %s\n", OUT_02))
cat(sprintf("  Sampling: KW=%dk, KS=%dk, Cliff=%dk, Plot=%dk\n\n",
            SAMPLE_CONFIG$kruskal_max_n / 1000, SAMPLE_CONFIG$ks_max_n / 1000,
            SAMPLE_CONFIG$cliff_max_n / 1000, SAMPLE_CONFIG$plot_max_n / 1000))

# =============================================================================
# Step 1: Load preprocessed data (or raw if needed)
# =============================================================================

cat("[Step 1] Loading data ...\n")

# Try preprocessed from 01_data_prep.R first
prep_dir <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
prep_dir <- prep_dir[grepl("01_Data_Prep", basename(prep_dir))]

dt_list <- list()
if (length(prep_dir) > 0) {
  prep_dir <- sort(prep_dir, decreasing = TRUE)[1]
  cat(sprintf("  Using preprocessed data from: %s\n", basename(prep_dir)))
  for (site in c("TN", "D6", "CY")) {
    fp <- file.path(prep_dir, sprintf("%s_preprocessed.fst", site))
    if (file.exists(fp)) {
      dt_list[[site]] <- read_fst_dt(fp)
    } else {
      dt_list[[site]] <- prep_site_dt(read_fst_dt(SITE_FILES[[site]]), site)
    }
    cat(sprintf("  %s: %s sessions\n", site, format(nrow(dt_list[[site]]), big.mark = ",")))
  }
} else {
  cat("  Preprocessing from raw .fst files ...\n")
  for (site in c("TN", "D6", "CY")) {
    dt_list[[site]] <- prep_site_dt(read_fst_dt(SITE_FILES[[site]]), site)
    cat(sprintf("  %s: %s sessions\n", site, format(nrow(dt_list[[site]]), big.mark = ",")))
  }
}

# =============================================================================
# Step 2: Kruskal-Wallis Test (3-group comparison)
# Methods → "Formal Kruskal–Wallis ... tests ... are reported in
#  the Supplementary Information"
# Supp Table 1
# =============================================================================

cat("\n[Step 2] Kruskal-Wallis tests ...\n")

kruskal_test <- function(dt_list, var, max_n = SAMPLE_CONFIG$kruskal_max_n) {
  vals <- lapply(names(dt_list), function(nm) {
    x <- dt_list[[nm]][[var]]
    x <- x[!is.na(x)]
    if (length(x) > max_n) x <- sample(x, max_n)
    x
  })
  df <- rbindlist(lapply(seq_along(vals), function(i)
    data.table(Site = names(dt_list)[i], Value = vals[[i]])))
  if (nrow(df) < 10) return(list(statistic = NA, p.value = NA))
  test <- kruskal.test(Value ~ Site, data = df)
  list(statistic = test$statistic, p.value = test$p.value)
}

kw_results <- rbindlist(lapply(ALL_PREDICTORS, function(var) {
  res <- kruskal_test(dt_list, var)
  data.table(Variable = var, KW_statistic = res$statistic, KW_pvalue = res$p.value)
}))
kw_results[, KW_sig := fifelse(KW_pvalue < 0.001, "***",
                        fifelse(KW_pvalue < 0.01, "**",
                         fifelse(KW_pvalue < 0.05, "*", "")))]

cat(sprintf("  %d variables tested. Significant (p<0.001): %d\n",
            nrow(kw_results), sum(kw_results$KW_pvalue < 0.001, na.rm = TRUE)))

# Save Supp Table 1
fwrite(kw_results, file.path(OUT_02, "Supp_Table_1_Kruskal_Wallis.csv"))

# =============================================================================
# Step 3: Pairwise KS Tests
# Supp Table 2
# =============================================================================

cat("\n[Step 3] Pairwise KS tests ...\n")

pairwise_ks <- function(dt_list, var, max_n = SAMPLE_CONFIG$ks_max_n) {
  sites <- names(dt_list)
  pairs <- combn(sites, 2, simplify = FALSE)
  rbindlist(lapply(pairs, function(pr) {
    x1 <- dt_list[[pr[1]]][[var]]; x1 <- x1[!is.na(x1)]
    x2 <- dt_list[[pr[2]]][[var]]; x2 <- x2[!is.na(x2)]
    if (length(x1) > max_n) x1 <- sample(x1, max_n)
    if (length(x2) > max_n) x2 <- sample(x2, max_n)
    if (length(x1) < 5 || length(x2) < 5)
      return(data.table(Pair = paste(pr, collapse = " vs "), D = NA, p.value = NA))
    ks <- suppressWarnings(ks.test(x1, x2))
    data.table(Pair = paste(pr, collapse = " vs "), D = ks$statistic, p.value = ks$p.value)
  }))
}

ks_results <- rbindlist(lapply(ALL_PREDICTORS, function(var) {
  res <- pairwise_ks(dt_list, var)
  res[, Variable := var]
  res
}))

# Wide format for Supp Table 2
ks_wide <- dcast(ks_results, Variable ~ Pair, value.var = "D")
fwrite(ks_wide, file.path(OUT_02, "Supp_Table_2_KS_Tests.csv"))
cat(sprintf("  %d variable × %d pairs = %d tests\n",
            length(ALL_PREDICTORS), 3, nrow(ks_results)))

# =============================================================================
# Step 4: Cliff's Delta (pairwise effect sizes)
# Table 1 Panel B
# Methods → "Cliff's δ ... interpreted using standard negligible, small,
#  medium, and large effect-size thresholds"
# =============================================================================

cat("\n[Step 4] Cliff's delta effect sizes ...\n")

pairwise_cliff <- function(dt_list, var, max_n = SAMPLE_CONFIG$cliff_max_n) {
  sites <- names(dt_list)
  pairs <- combn(sites, 2, simplify = FALSE)
  rbindlist(lapply(pairs, function(pr) {
    x1 <- dt_list[[pr[1]]][[var]]; x1 <- x1[!is.na(x1)]
    x2 <- dt_list[[pr[2]]][[var]]; x2 <- x2[!is.na(x2)]
    if (length(x1) > max_n) x1 <- sample(x1, max_n)
    if (length(x2) > max_n) x2 <- sample(x2, max_n)
    # Note: integer overflow (n1*n2 > 2^31) is handled inside cliffs_delta_fast()
    #       via as.double() — no additional downsampling needed here.
    if (length(x1) < 5 || length(x2) < 5)
      return(data.table(Pair = paste(pr, collapse = " vs "),
                        Cliffs_Delta = NA_real_,
                        Magnitude    = NA_character_))
    delta <- cliffs_delta_fast(x1, x2)
    data.table(
      Pair         = paste(pr, collapse = " vs "),
      Cliffs_Delta = round(delta, 3),
      Magnitude    = cliffs_delta_magnitude(delta)
    )
  }))
}

pb <- txtProgressBar(min = 0, max = length(ALL_PREDICTORS), style = 3)
effect_results <- rbindlist(lapply(seq_along(ALL_PREDICTORS), function(i) {
  var <- ALL_PREDICTORS[i]
  res <- pairwise_cliff(dt_list, var)
  res[, Variable := var]
  setTxtProgressBar(pb, i)
  res
}))
close(pb)

# --- Build Table 1 Panel B (ordered by max |δ|) ---
delta_wide <- dcast(effect_results, Variable ~ Pair, value.var = "Cliffs_Delta")
mag_wide   <- dcast(effect_results, Variable ~ Pair, value.var = "Magnitude")

# Rename columns
delta_cols <- setdiff(names(delta_wide), "Variable")
for (col in delta_cols) {
  new_name <- paste0("Delta_", gsub(" vs ", "_", col))
  setnames(delta_wide, col, new_name)
}
mag_cols <- setdiff(names(mag_wide), "Variable")
for (col in mag_cols) {
  new_name <- paste0("Mag_", gsub(" vs ", "_", col))
  setnames(mag_wide, col, new_name)
}

table1b <- merge(delta_wide, mag_wide, by = "Variable")

# Compute max |δ| and the pair with max
delta_col_names <- grep("^Delta_", names(table1b), value = TRUE)
table1b[, Max_abs_delta := apply(.SD, 1, function(r) max(abs(r), na.rm = TRUE)),
        .SDcols = delta_col_names]
table1b[, Max_pair := apply(.SD, 1, function(r) {
  idx <- which.max(abs(r))
  gsub("Delta_", "", delta_col_names[idx])
}), .SDcols = delta_col_names]
table1b[, Max_magnitude := vapply(Max_abs_delta, cliffs_delta_magnitude, character(1))]

# Order by descending max |δ| (matching manuscript Table 1B)
setorder(table1b, -Max_abs_delta)

cat("\n\n  ──────────────────────────────────────────────────────────────\n")
cat("  Table 1 Panel B: Cliff's δ effect sizes\n")
cat("  ──────────────────────────────────────────────────────────────\n")
print(table1b[, .(Variable, Delta_TN_D6, Delta_TN_CY, Delta_D6_CY,
                   Max_abs_delta = round(Max_abs_delta, 2), Max_magnitude, Max_pair)])
cat("\n")

# Verify against manuscript
cat("  [Verify] Top 3 by max |δ| should be:\n")
cat("    1. Dialysate_Temperature: 0.81 (Large, TN vs D6)\n")
cat("    2. Age: 0.45 (Medium, TN vs CY)\n")
cat("    3. Respiratory_Rate: 0.33 (Medium, TN vs CY)\n")
cat(sprintf("    Got: 1. %s: %.2f (%s)\n", table1b$Variable[1],
            table1b$Max_abs_delta[1], table1b$Max_magnitude[1]))
cat(sprintf("         2. %s: %.2f (%s)\n", table1b$Variable[2],
            table1b$Max_abs_delta[2], table1b$Max_magnitude[2]))
cat(sprintf("         3. %s: %.2f (%s)\n", table1b$Variable[3],
            table1b$Max_abs_delta[3], table1b$Max_magnitude[3]))

# Save Table 1B
fwrite(table1b, file.path(OUT_02, "Table_1B_Cliffs_Delta.csv"))

# =============================================================================
# Step 5: Merge all heterogeneity stats
# =============================================================================

cat("\n[Step 5] Merging heterogeneity summary ...\n")

ks_wide2 <- dcast(ks_results, Variable ~ Pair, value.var = "D")
ks_cols2 <- setdiff(names(ks_wide2), "Variable")
for (col in ks_cols2) setnames(ks_wide2, col, paste0("KS_", gsub(" vs ", "_", col)))

het_summary <- merge(kw_results, ks_wide2, by = "Variable", all.x = TRUE)
het_summary <- merge(het_summary, delta_wide, by = "Variable", all.x = TRUE)
het_summary <- merge(het_summary, mag_wide, by = "Variable", all.x = TRUE)

ks_cn <- grep("^KS_", names(het_summary), value = TRUE)
delta_cn <- grep("^Delta_", names(het_summary), value = TRUE)
het_summary[, Avg_KS    := rowMeans(.SD, na.rm = TRUE), .SDcols = ks_cn]
het_summary[, Avg_Delta := rowMeans(abs(.SD), na.rm = TRUE), .SDcols = delta_cn]
setorder(het_summary, -Avg_KS)

fwrite(het_summary, file.path(OUT_02, "Heterogeneity_Summary.csv"))

# =============================================================================
# Step 6: Visualizations
# =============================================================================

if (ENABLE_PLOTTING) {
  cat("\n[Step 6] Generating visualizations ...\n")

  # --- 6.1 Sampling for plots ---
  dt_combined <- rbindlist(lapply(names(dt_list), function(nm) {
    dt <- dt_list[[nm]][, c(..CONT_COLS, "site"), with = FALSE]
    n_sample <- min(nrow(dt), SAMPLE_CONFIG$plot_max_n)
    if (nrow(dt) > n_sample) dt <- dt[sample(.N, n_sample)]
    dt
  }))
  dt_melt <- melt(dt_combined, id.vars = "site",
                  variable.name = "Variable", value.name = "Value")

  # --- 6.2 Cliff's δ Heatmap (Fig 2a) ---
  cat("  Fig 2a: Cliff's δ heatmap ...\n")

  effect_long <- melt(het_summary[, c("Variable", delta_cn), with = FALSE],
                      id.vars = "Variable", variable.name = "Pair", value.name = "Delta")
  effect_long[, Pair := gsub("Delta_", "", Pair)]
  effect_long[, Pair := gsub("_", " vs ", Pair)]
  effect_long[, Variable := factor(Variable, levels = rev(het_summary$Variable))]

  p_delta <- ggplot(effect_long, aes(x = Pair, y = Variable, fill = Delta)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Delta)), size = 3) +
    scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c",
                         midpoint = 0, name = expression("Cliff's " * delta),
                         limits = c(-1, 1)) +
    labs(title = expression("Cross-site Cliff's " * delta * " effect sizes"),
         subtitle = expression("|" * delta * "|: <0.147 Negligible, <0.33 Small, <0.474 Medium, " >= "0.474 Large"),
         x = "Center pair", y = "Predictor") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = 8),
          panel.grid = element_blank())

  ggsave(file.path(OUT_02, "Fig_2a_Cliffs_Delta_Heatmap.png"), p_delta,
         width = 8, height = max(6, nrow(het_summary) * 0.4), dpi = 300)

  # --- 6.3 Density overlays (Fig 2b) ---
  cat("  Fig 2b: Density overlays ...\n")

  # Top variables by distributional shift (medium/large δ first, then small)
  top_vars <- het_summary$Variable[1:min(12, nrow(het_summary))]
  dt_top <- dt_melt[Variable %in% top_vars]
  n_rows_facet <- ceiling(length(top_vars) / 4)

  p_density <- ggplot(dt_top, aes(x = Value, fill = site, color = site)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values = SITE_COLORS) +
    scale_color_manual(values = SITE_COLORS) +
    facet_wrap(~ Variable, scales = "free", ncol = 4) +
    labs(title = "Cross-site distributional comparison (density overlays)",
         subtitle = sprintf("Top %d predictors by KS distance (sampled n = %dk per center)",
                            length(top_vars), SAMPLE_CONFIG$plot_max_n / 1000),
         x = "Value", y = "Density", fill = "Center", color = "Center") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "top")

  ggsave(file.path(OUT_02, "Fig_2b_Density_Overlays.png"), p_density,
         width = 14, height = 3.5 * n_rows_facet, dpi = 300)

  # --- 6.4 KS Heatmap ---
  cat("  KS heatmap ...\n")

  ks_long <- melt(het_summary[, c("Variable", ks_cn), with = FALSE],
                  id.vars = "Variable", variable.name = "Pair", value.name = "KS_D")
  ks_long[, Pair := gsub("KS_", "", Pair)]
  ks_long[, Pair := gsub("_", " vs ", Pair)]
  ks_long[, Variable := factor(Variable, levels = rev(het_summary$Variable))]

  p_ks <- ggplot(ks_long, aes(x = Pair, y = Variable, fill = KS_D)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", KS_D)), size = 3) +
    scale_fill_gradient2(low = "#2ecc71", mid = "#f1c40f", high = "#e74c3c",
                         midpoint = 0.1, name = "KS D") +
    labs(title = "KS statistic heatmap", x = "Center pair", y = "Predictor") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid = element_blank())

  ggsave(file.path(OUT_02, "KS_Heatmap.png"), p_ks,
         width = 8, height = max(6, nrow(het_summary) * 0.4), dpi = 300)

  # --- 6.5 ECDF plots ---
  cat("  ECDF plots ...\n")

  p_ecdf <- ggplot(dt_top, aes(x = Value, color = site)) +
    stat_ecdf(linewidth = 0.8) +
    scale_color_manual(values = SITE_COLORS) +
    facet_wrap(~ Variable, scales = "free_x", ncol = 4) +
    labs(title = "Empirical CDF comparison",
         x = "Value", y = "Cumulative probability", color = "Center") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "top")

  ggsave(file.path(OUT_02, "ECDF_Overlays.png"), p_ecdf,
         width = 14, height = 3.5 * n_rows_facet, dpi = 300)

  cat("  Plots saved.\n")
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 02 complete.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

output_files <- list.files(OUT_02, full.names = FALSE)
cat(sprintf("  Output directory: %s\n", OUT_02))
for (f in output_files) {
  fsize <- file.size(file.path(OUT_02, f))
  cat(sprintf("    %s (%.1f KB)\n", f, fsize / 1024))
}
cat("\n")
