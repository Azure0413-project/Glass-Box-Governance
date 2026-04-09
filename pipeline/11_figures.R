# =============================================================================
# 11_figures.R — Publication Figures
# IDH EBM Governance Reproducibility Pipeline
#
# Outputs:
#   Fig 1  — Study design / governance framework (conceptual diagram)
#   Fig 2  — Cross-site heterogeneity (Cliff's δ heatmap + density overlays)
#   Fig 3  — ShapeQC quadrant + individual shape overlays
#   Fig 4  — External validation (ROC/PRC curves, calibration, DCA)
#   Fig 5  — Reclassification waterfall / probability shift
#
# This script aggregates data from all modules and generates publication-quality figures.
# Most figures are already produced in individual modules; this script creates composites.
#
# Prerequisites: src("R/00_config.R"); source("00_utils_r.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(gridExtra)
})

if (!exists("SITE_FILES")) src("R/00_config.R")
if (!exists("read_fst_dt")) src("R/00_utils_r.R")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 11: Publication Figures\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

FIG_DIR <- file.path(RUN_DIR, paste0("Figures_Publication_", RUN_TS))
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  Output: %s\n\n", basename(FIG_DIR)))

# Common theme for all figures
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
    plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "gray40"),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
}

# =============================================================================
# Fig 2: Cross-site Heterogeneity
# (Composite: Panel A = Cliff's δ heatmap, Panel B = density overlays)
# =============================================================================

cat("[Fig 2] Cross-site Heterogeneity ...\n")

het_dir <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
het_dir <- het_dir[grepl("02_Heterogeneity", basename(het_dir))]

if (length(het_dir) > 0) {
  het_dir <- sort(het_dir, decreasing = TRUE)[1]
  delta_file <- file.path(het_dir, "Table_1B_Cliffs_Delta.csv")

  if (file.exists(delta_file)) {
    table1b <- fread(delta_file)
    delta_cols <- grep("^Delta_", names(table1b), value = TRUE)

    # Panel A: Cliff's δ heatmap
    delta_long <- melt(table1b[, c("Variable", delta_cols), with = FALSE],
                       id.vars = "Variable", variable.name = "Pair", value.name = "Delta")
    delta_long[, Pair := gsub("Delta_", "", Pair)]
    delta_long[, Pair := gsub("_", " vs ", Pair)]
    delta_long[, Variable := factor(Variable, levels = rev(table1b$Variable))]

    p2a <- ggplot(delta_long, aes(x = Pair, y = Variable, fill = Delta)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", Delta)), size = 2.8) +
      scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c",
                           midpoint = 0, name = expression("Cliff's " * delta),
                           limits = c(-1, 1)) +
      labs(title = "(a) Cross-site Cliff's \u03b4 effect sizes",
           x = NULL, y = NULL) +
      theme_pub() +
      theme(panel.grid = element_blank(), legend.position = "right",
            axis.text.y = element_text(size = 9))

    ggsave(file.path(FIG_DIR, "Fig2a_Cliffs_Delta.png"), p2a,
           width = 7, height = max(5, nrow(table1b) * 0.35), dpi = 300)
    cat("  Fig 2a saved.\n")
  }
} else {
  cat("  [Skip] Module 02 output not found.\n")
}

# =============================================================================
# Fig 3: ShapeQC Governance Quadrant
# =============================================================================

cat("[Fig 3] ShapeQC Governance ...\n")

qc_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
qc_dirs <- qc_dirs[grepl("ShapeQC_v3\\.0_", basename(qc_dirs))]

if (length(qc_dirs) > 0) {
  qc_dir <- sort(qc_dirs, decreasing = TRUE)[1]
  qc_csv <- file.path(qc_dir, "Table_4A_Decision_Matrix.csv")

  if (file.exists(qc_csv)) {
    qc_dt <- fread(qc_csv)
    quad_dt <- qc_dt[!is.na(C0_trim) & !is.na(median_J) & is_discrete == FALSE]

    if (nrow(quad_dt) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
      p3 <- ggplot(quad_dt, aes(x = C0_trim, y = median_J)) +
        annotate("rect", xmin = SHAPEQC_TAU, xmax = Inf, ymin = -Inf, ymax = SHAPEQC_KAPPA,
                 fill = "#27ae60", alpha = 0.10) +
        annotate("rect", xmin = SHAPEQC_TAU, xmax = Inf, ymin = SHAPEQC_KAPPA, ymax = Inf,
                 fill = "#f1c40f", alpha = 0.10) +
        annotate("rect", xmin = -Inf, xmax = SHAPEQC_TAU, ymin = -Inf, ymax = Inf,
                 fill = "#e74c3c", alpha = 0.10) +
        geom_vline(xintercept = SHAPEQC_TAU, linetype = "dashed", color = "gray40") +
        geom_hline(yintercept = SHAPEQC_KAPPA, linetype = "dashed", color = "gray40") +
        geom_point(aes(color = light), size = 3.5) +
        ggrepel::geom_text_repel(aes(label = feature), size = 3.2, max.overlaps = 20,
                                  seed = 42) +
        scale_color_manual(values = TIER_COLORS) +
        annotate("text", x = 0.95, y = 0.5, label = "Green\n(Trust)", color = "#27ae60",
                 fontface = "bold", size = 3.5, alpha = 0.6) +
        annotate("text", x = 0.95, y = 3.5, label = "Yellow\n(Refine)", color = "#c5a000",
                 fontface = "bold", size = 3.5, alpha = 0.6) +
        annotate("text", x = 0.20, y = 2.0, label = "Red\n(Non-\ntransportable)",
                 color = "#c0392b", fontface = "bold", size = 3.5, alpha = 0.6) +
        labs(title = "Feature-Level Governance Decision (ShapeQC v3.0)",
             subtitle = sprintf("\u03c4 = %.2f (concordance) | \u03ba = %.2f (jaggedness) | Trim: %s",
                                SHAPEQC_TAU, SHAPEQC_KAPPA, PRIMARY_TRIM),
             x = expression(C[core] ~ "(trimmed cross-fold concordance)"),
             y = expression(J ~ "(jaggedness index)"),
             color = "Tier") +
        theme_pub() +
        theme(legend.position = "right")

      ggsave(file.path(FIG_DIR, "Fig3_ShapeQC_Quadrant.png"), p3,
             width = 9, height = 7, dpi = 300)
      cat("  Fig 3 saved.\n")
    }
  }
} else {
  cat("  [Skip] Module 05 output not found.\n")
}

# =============================================================================
# Fig 4: External Validation Composite
# (Panel A: AUROC forest, Panel B: AUPRC forest, Panel C: DCA overlay)
# =============================================================================

cat("[Fig 4] External Validation ...\n")

summ_file <- file.path(RUN_DIR, "IECV_iteration_summary.csv")
if (!file.exists(summ_file)) {
  # Try within Run_* directories
  run_dirs <- list.dirs(BASE_DIR, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("^Run_", basename(run_dirs))]
  for (rd in sort(run_dirs, decreasing = TRUE)) {
    sf <- file.path(rd, "IECV_iteration_summary.csv")
    if (file.exists(sf)) { summ_file <- sf; break }
  }
}

if (file.exists(summ_file)) {
  summ <- fread(summ_file)

  if ("external_AUROC" %in% names(summ) && nrow(summ) > 0) {
    summ[, label := sprintf("%s (Seed %s)", external, seed_id)]

    # Panel A: AUROC by iteration
    p4a <- ggplot(summ, aes(x = external_AUROC, y = reorder(label, external_AUROC))) +
      geom_point(aes(color = external), size = 3) +
      geom_vline(xintercept = mean(summ$external_AUROC), linetype = "dashed", color = "gray50") +
      scale_color_manual(values = SITE_COLORS) +
      labs(title = "(a) External AUROC", x = "AUROC", y = NULL) +
      theme_pub() + theme(legend.position = "none")

    # Panel B: AUPRC by iteration
    p4b <- ggplot(summ, aes(x = external_AUPRC, y = reorder(label, external_AUPRC))) +
      geom_point(aes(color = external), size = 3) +
      geom_vline(xintercept = mean(summ$external_AUPRC), linetype = "dashed", color = "gray50") +
      scale_color_manual(values = SITE_COLORS) +
      labs(title = "(b) External AUPRC", x = "AUPRC", y = NULL) +
      theme_pub() + theme(legend.position = "none")

    p4_combined <- gridExtra::grid.arrange(p4a, p4b, ncol = 2)
    ggsave(file.path(FIG_DIR, "Fig4_External_Validation.png"),
           p4_combined, width = 14, height = 5, dpi = 300)
    cat("  Fig 4 saved.\n")
  }
} else {
  cat("  [Skip] IECV summary not found.\n")
}

# =============================================================================
# Fig 5: Reclassification Probability Shift
# =============================================================================

cat("[Fig 5] Reclassification ...\n")

nri_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
nri_dirs <- nri_dirs[grepl("NRI_Reclassification", basename(nri_dirs))]

if (length(nri_dirs) > 0) {
  nri_dir <- sort(nri_dirs, decreasing = TRUE)[1]
  rds_file <- file.path(nri_dir, "NRI_Complete_Results.rds")

  if (file.exists(rds_file)) {
    nri_res <- readRDS(rds_file)

    # Summary bar: NRI components across iterations
    nri_summ <- rbindlist(lapply(nri_res$nri_idi, function(x) {
      data.table(
        iter = x$iter_tag, center = x$held_out_center,
        NRI_events = x$cat_nri$NRI_events,
        NRI_nonevents = x$cat_nri$NRI_nonevents,
        cat_NRI = x$cat_nri$NRI,
        pct_moved = x$cat_nri$pct_moved
      )
    }))

    nri_long <- melt(nri_summ[, .(iter, center, NRI_events, NRI_nonevents)],
                     id.vars = c("iter", "center"),
                     variable.name = "Component", value.name = "Value")
    nri_long[, Component := fifelse(Component == "NRI_events", "NRI+ (events)",
                                     "NRI\u2212 (non-events)")]

    p5 <- ggplot(nri_long, aes(x = iter, y = Value, fill = Component)) +
      geom_col(position = "dodge", width = 0.6) +
      geom_hline(yintercept = 0, color = "gray50") +
      scale_fill_manual(values = c("NRI+ (events)" = "#e74c3c",
                                    "NRI\u2212 (non-events)" = "#3498db")) +
      labs(title = "Categorical NRI decomposition across IECV iterations",
           subtitle = sprintf("Thresholds: %s | Positive NRI = improvement",
                              paste0(NRI_THRESHOLDS * 100, "%", collapse = " / ")),
           x = "IECV Iteration", y = "NRI Component", fill = NULL) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    ggsave(file.path(FIG_DIR, "Fig5_NRI_Decomposition.png"), p5,
           width = 10, height = 6, dpi = 300)
    cat("  Fig 5 saved.\n")
  }
} else {
  cat("  [Skip] Module 10 output not found.\n")
}

# =============================================================================
# Supp Figure: Trim Sensitivity Trend (from ShapeQC)
# =============================================================================

cat("[Supp] Trim sensitivity ...\n")

if (length(qc_dirs) > 0) {
  qc_dir <- sort(qc_dirs, decreasing = TRUE)[1]
  sens_csv <- file.path(qc_dir, "Supp_S14_Trim_Sensitivity.csv")

  if (file.exists(sens_csv)) {
    sens_wide <- fread(sens_csv)
    sens_long <- melt(sens_wide, id.vars = "feature",
                      variable.name = "config", value.name = "C0_trim")
    sens_long[, config := factor(config, levels = c("T1", "T2", "T3", "T4", "T5"))]

    # Color by tier
    tier_map <- FEATURE_TIER_MAP[sens_long$feature]
    sens_long[, tier := tier_map]

    p_sens <- ggplot(sens_long, aes(x = config, y = C0_trim, group = feature, color = tier)) +
      geom_line(linewidth = 0.8, alpha = 0.7) +
      geom_point(size = 2) +
      geom_hline(yintercept = SHAPEQC_TAU, linetype = "dashed", color = "gray40") +
      scale_color_manual(values = TIER_COLORS, na.value = "gray60") +
      labs(title = "C_core sensitivity across trim configurations",
           subtitle = "T3 (P0.5\u2013P99.5) = primary configuration",
           x = "Trim configuration", y = expression(C[core]),
           color = "Governance tier") +
      theme_pub()

    ggsave(file.path(FIG_DIR, "Supp_Trim_Sensitivity.png"), p_sens,
           width = 9, height = 6, dpi = 300)
    cat("  Supp trim sensitivity saved.\n")
  }
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Module 11 complete.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

output_files <- list.files(FIG_DIR, full.names = FALSE)
cat(sprintf("  Output: %s\n", FIG_DIR))
for (f in output_files) {
  fsize <- file.size(file.path(FIG_DIR, f))
  cat(sprintf("    %s (%.1f KB)\n", f, fsize / 1024))
}
cat("\n  Pipeline complete.\n\n")
