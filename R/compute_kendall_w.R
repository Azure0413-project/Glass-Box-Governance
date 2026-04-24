# =============================================================================
# compute_kendall_w.R
# ---------------------------------------------------------------------------
# Cross-IECV iteration feature-importance rank stability analysis
#
# Verifies the manuscript claim:
#   "importance rankings were highly stable across IECV iterations
#    (Kendall's W = 0.88; Supplementary Table S2)"
#
# Computes two complementary indicators:
#   [Option A] Kendall's W -- 16 features x 6 iterations overall rank concordance
#   [Option B] Top-3 rank invariance -- whether the top 3 features are identical
#
# Prerequisites:
#   - RUN_DIR must contain GLOBAL_feature_importance_raw.csv
#     (produced by extract_ebm_importance.R; columns: feature, model, importance)
#   - Falls back to GLOBAL_feature_importance.csv (03_iecv_ebm.R output)
#
# Usage:
#   source("00_config.R")
#   source("compute_kendall_w.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

stopifnot("RUN_DIR must be set" = exists("RUN_DIR"))

# --- locate the per-iteration importance file ---
imp_path <- file.path(RUN_DIR, "GLOBAL_feature_importance_raw.csv")
if (!file.exists(imp_path)) {
  imp_path <- file.path(RUN_DIR, "GLOBAL_feature_importance.csv")
}
if (!file.exists(imp_path)) stop(sprintf("Feature importance CSV not found under %s", RUN_DIR))

cat(sprintf("\n[Input] %s\n", imp_path))
imp <- fread(imp_path)

# Defensive: some versions store the rank as 'rank', but we recompute here
# so the results are reproducible regardless of upstream ties-handling.
imp[, rank_in_model := frank(-importance, ties.method = "average"), by = model]

# --- build features x models rank matrix ---
rank_mat   <- dcast(imp, feature ~ model, value.var = "rank_in_model")
feat_names <- rank_mat$feature
M          <- as.matrix(rank_mat[, -1, with = FALSE])
rownames(M) <- feat_names

n <- nrow(M)   # subjects  = features (expected 16)
m <- ncol(M)   # raters    = iterations (expected 6)

cat(sprintf("[Shape] %d features x %d iterations\n", n, m))
stopifnot("Need at least 2 iterations and 3 features" = (m >= 2 && n >= 3))

# =============================================================================
# [A] Kendall's W -- manual computation with ties correction
# =============================================================================
# W = 12 S / (m^2 (n^3 - n) - m * sum(T_j))
# where:
#   R_i = sum of ranks for item i across m raters
#   S   = sum_i (R_i - mean(R_i))^2
#   T_j = sum over tied groups g in rater j of (t_g^3 - t_g)
# =============================================================================

R_i   <- rowSums(M)
R_bar <- mean(R_i)
S     <- sum((R_i - R_bar)^2)

T_j <- apply(M, 2, function(col) {
  tab <- table(col)
  sum(tab^3 - tab)
})
T_sum <- sum(T_j)

W_ties <- (12 * S) / (m^2 * (n^3 - n) - m * T_sum)
W_unc  <- (12 * S) / (m^2 * (n^3 - n))

chi2 <- m * (n - 1) * W_ties
pval <- pchisq(chi2, df = n - 1, lower.tail = FALSE)

cat("\n=============================================================\n")
cat("[Option A] Kendall's W -- 16 features x 6 iterations\n")
cat("=============================================================\n")
cat(sprintf("  W (ties-corrected) : %.4f\n", W_ties))
cat(sprintf("  W (uncorrected)    : %.4f\n", W_unc))
cat(sprintf("  chi-square(df=%d) : %.2f\n", n - 1, chi2))
cat(sprintf("  p-value            : %.3e\n", pval))

# Sanity cross-check with irr if present
if (requireNamespace("irr", quietly = TRUE)) {
  res_irr <- irr::kendall(M, correct = TRUE)
  cat(sprintf("  irr::kendall W     : %.4f   (independent check)\n",
              res_irr$value))
} else {
  cat("  (install.packages('irr') for an independent cross-check)\n")
}

# =============================================================================
# [B] Top-3 invariance -- same 3 features? same order?
# =============================================================================

top3_per_iter <- apply(M, 2, function(col) rownames(M)[order(col)[1:3]])
# columns = iterations; rows = position 1/2/3

top3_sets    <- lapply(seq_len(m), function(j) sort(top3_per_iter[, j]))
top3_ordered <- lapply(seq_len(m), function(j) top3_per_iter[, j])
same_set   <- length(unique(top3_sets))    == 1L
same_order <- length(unique(top3_ordered)) == 1L

cat("\n=============================================================\n")
cat("[Option B] Top-3 rank structure per iteration\n")
cat("=============================================================\n")
print(top3_per_iter)

cat(sprintf("\n  Top-3 feature set identical across all %d iters : %s\n",
            m, ifelse(same_set, "YES", "NO")))
cat(sprintf("  Top-3 order identical (positions 1/2/3 match)  : %s\n",
            ifelse(same_order, "YES", "NO")))

# Rank range of each top-3 member (how far off the podium did they ever drop?)
top3_consensus <- names(sort(rowMeans(M))[1:3])
cat("\n  Per-iteration rank range of top-3 consensus members:\n")
for (f in top3_consensus) {
  rks <- M[f, ]
  cat(sprintf("    %-25s  rank range = [%d, %d]  (mean %.2f)\n",
              f, as.integer(min(rks)), as.integer(max(rks)), mean(rks)))
}

# =============================================================================
# [Supplementary] Per-iteration Spearman rho vs consensus rank
# =============================================================================
# Useful as a backup sentence: "each iteration's ranking agreed with the
# consensus at rho >= X" -- complements Kendall's W.
# =============================================================================

consensus_rank <- rank(R_i)   # 1 = most important (smallest rank sum)
spearmans <- apply(M, 2, function(col) cor(col, consensus_rank, method = "spearman"))

cat("\n=============================================================\n")
cat("[Supplementary] Per-iteration Spearman rho vs consensus rank\n")
cat("=============================================================\n")
for (j in seq_len(m)) {
  cat(sprintf("  %-32s rho = %.3f\n", colnames(M)[j], spearmans[j]))
}
cat(sprintf("  Mean rho = %.3f   Min rho = %.3f\n",
            mean(spearmans), min(spearmans)))

# =============================================================================
# [Save]
# =============================================================================

out_csv <- file.path(RUN_DIR, "Supp_S2_rank_stability.csv")

rank_long <- melt(rank_mat, id.vars = "feature",
                  variable.name = "iteration", value.name = "rank")
setorder(rank_long, feature, iteration)
fwrite(rank_long, out_csv)

summary_out <- data.table(
  metric = c("kendall_W_ties_corrected", "kendall_W_uncorrected",
             "chi_square", "df", "p_value",
             "top3_set_invariant", "top3_order_invariant",
             "mean_spearman_vs_consensus", "min_spearman_vs_consensus",
             "n_features", "n_iterations"),
  value  = c(sprintf("%.4f", W_ties), sprintf("%.4f", W_unc),
             sprintf("%.3f", chi2), as.character(n - 1),
             sprintf("%.3e", pval),
             as.character(same_set), as.character(same_order),
             sprintf("%.3f", mean(spearmans)),
             sprintf("%.3f", min(spearmans)),
             as.character(n), as.character(m))
)
fwrite(summary_out, file.path(RUN_DIR, "Supp_S2_rank_stability_summary.csv"))

cat(sprintf("\n[Saved]\n  %s\n  %s\n",
            out_csv,
            file.path(RUN_DIR, "Supp_S2_rank_stability_summary.csv")))

# =============================================================================
# [Decision guide]
# =============================================================================

cat("\n=============================================================\n")
cat("[Decision guide -- which claim to make in the manuscript]\n")
cat("=============================================================\n")
if (same_order) {
  cat("  -> Top-3 are RANK-INVARIANT across all iterations.\n")
  cat("     Strongest claim available. Prefer Option B wording:\n")
  cat("       'The three most influential predictors retained\n")
  cat("        ranks 1, 2, and 3 in every IECV iteration.'\n")
  cat(sprintf("     Keep Kendall's W = %.2f as aggregate backup in S2.\n", W_ties))
} else if (same_set) {
  cat("  -> Top-3 feature SET is invariant but internal order shifts.\n")
  cat("     Safer wording:\n")
  cat("       'The same three predictors occupied the top three ranks\n")
  cat("        in every iteration (internal order varied).'\n")
  cat(sprintf("     Pair with Kendall's W = %.2f for aggregate stability.\n", W_ties))
} else {
  cat("  -> Top-3 set NOT invariant across all iterations.\n")
  cat("     Stick with Option A:\n")
  cat(sprintf("       'importance rankings were highly stable across IECV\n"))
  cat(sprintf("        iterations (Kendall's W = %.2f)'\n", W_ties))
  cat("     Do NOT claim 'invariant top-3'.\n")
}
cat("\n")
