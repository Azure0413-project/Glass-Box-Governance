#!/usr/bin/env python3
# =============================================================================
# extract_and_plot.py — EBM Feature Importance Extraction & Visualization
# IDH EBM Governance Reproducibility Pipeline
#
# Dynamically discovers iteration directories and model files.
# Extracts feature importance using the official InterpretML weighted method.
# Produces:
#   1. Feature importance bar chart (grouped by governance tier, colored by tier)
#   2. Shape function overlay plots (colored by IECV iteration, grouped by site)
#   3. CSV exports of raw + summary importance tables
#
# Usage (inside Docker):
#   cd /path/to/EBM_Paper
#   python extract_and_plot.py                          # auto-detect RUN_DIR
#   python extract_and_plot.py --run-dir /path/to/dir   # explicit path
#
# No hardcoded paths, version strings, or iteration names.
# =============================================================================

import os
import sys
import re
import glob
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import joblib

# Suppress non-critical warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Matplotlib backend for headless rendering
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec


# =============================================================================
# Section 1: Configuration — all derived dynamically, no hardcoding
# =============================================================================

# Governance tier assignments — read from the EBM model's feature list
# These tiers are determined by ShapeQC (C_core concordance × Jaggedness)
# and reflect the CURRENT pipeline run, not a stale past version.
#
# Tier logic (from 00_config.R Section 10):
#   Green:  C_core >= tau AND J <= kappa  → Trust
#   Yellow: C_core >= tau AND J >  kappa  → Refine (smooth)
#   Red:    C_core <  tau                → Non-transportable (zero)
#   Gray:   Manual review (tied-rank safeguard)

TIER_DEFINITIONS = {
    "Green":  {"color": "#27ae60", "label": "Green (Trust)",              "order": 0},
    "Yellow": {"color": "#f1c40f", "label": "Yellow (Refine)",            "order": 1},
    "Red":    {"color": "#e74c3c", "label": "Red (Non-transportable)",    "order": 2},
    "Gray":   {"color": "#95a5a6", "label": "Gray (Manual review)",       "order": 3},
}

# Site color families for IECV iteration overlays (2 shades per site)
def build_site_palette(sites):
    """
    Dynamically build a color palette for sites.
    Each site gets a base color + lighter shade for seed variants.
    """
    base_colors = [
        ("#1f77b4", "#aec7e8"),  # Blue family
        ("#ff7f0e", "#ffbb78"),  # Orange family
        ("#2ca02c", "#98df8a"),  # Green family
        ("#d62728", "#ff9896"),  # Red family
        ("#9467bd", "#c5b0d5"),  # Purple family
        ("#8c564b", "#c49c94"),  # Brown family
    ]
    palette = {}
    for i, site in enumerate(sorted(set(sites))):
        idx = i % len(base_colors)
        palette[site] = base_colors[idx]
    return palette


# =============================================================================
# Section 2: Auto-discovery of iteration directories and model files
# =============================================================================

def discover_run_dir(base_path):
    """
    Auto-discover the RUN_DIR containing IECV iteration directories.
    Searches for directories matching the pattern Iter*_External_*.
    """
    # First check if base_path itself contains iteration dirs
    iter_dirs = glob.glob(os.path.join(base_path, "Iter*_External_*"))
    if iter_dirs:
        return base_path

    # Search subdirectories (e.g., checking_joblib, output/Run_*)
    for subdir in ["checking_joblib", "output"]:
        candidate = os.path.join(base_path, subdir)
        if os.path.isdir(candidate):
            iter_dirs = glob.glob(os.path.join(candidate, "Iter*_External_*"))
            if iter_dirs:
                return candidate

    # Search for Run_* directories
    run_dirs = sorted(glob.glob(os.path.join(base_path, "output", "Run_*")))
    for rd in reversed(run_dirs):
        iter_dirs = glob.glob(os.path.join(rd, "Iter*_External_*"))
        if iter_dirs:
            return rd

    return None


def parse_iter_tag(dirname):
    """
    Parse iteration metadata from directory name.
    Pattern: Iter{N}_External_{SITE}_Seed{S}
    Returns dict with iter, site, seed, or None if unparseable.
    """
    m = re.match(r"Iter(\d+)_External_(\w+)_Seed(\d+)", dirname)
    if m:
        return {
            "iter": int(m.group(1)),
            "site": m.group(2),
            "seed": int(m.group(3)),
            "tag": dirname,
        }
    # Fallback: try more flexible pattern
    m = re.match(r"Iter(\d+)_External_(\w+)", dirname)
    if m:
        return {
            "iter": int(m.group(1)),
            "site": m.group(2),
            "seed": 0,
            "tag": dirname,
        }
    return None


def find_models(run_dir):
    """
    Discover all IECV iteration directories and their EBM model files.
    Returns list of dicts with keys: iter_dir, model_path, meta.
    """
    iter_dirs = sorted(glob.glob(os.path.join(run_dir, "Iter*_External_*")))
    if not iter_dirs:
        raise FileNotFoundError(
            f"No iteration directories (Iter*_External_*) found in {run_dir}"
        )

    models = []
    for d in iter_dirs:
        if not os.path.isdir(d):
            continue
        meta = parse_iter_tag(os.path.basename(d))
        if meta is None:
            print(f"  [WARN] Skipping unparseable dir: {os.path.basename(d)}")
            continue

        # Find model file — prefer POSTHOC governance models, fall back to Final_EBM
        model_files = glob.glob(os.path.join(d, "*_Final_EBM*.joblib"))
        # Exclude any LOFO/sensitivity variants
        model_files = [f for f in model_files if "LOFO" not in f and "Sensitivity" not in f]

        if not model_files:
            print(f"  [WARN] No EBM model found in {os.path.basename(d)}")
            continue

        # Prefer POSTHOC governance model if available
        posthoc = [f for f in model_files if "POSTHOC" in f or "Governance" in f]
        model_path = posthoc[0] if posthoc else model_files[0]

        models.append({
            "iter_dir": d,
            "model_path": model_path,
            "meta": meta,
        })

    return models


# =============================================================================
# Section 3: Feature importance extraction (official InterpretML method)
# =============================================================================

def extract_importance(model_path, verbose=True):
    """
    Extract feature importance from a trained EBM model.
    Official method: weighted mean of |score| per feature.
    """
    if verbose:
        print(f"    Loading: {os.path.basename(model_path)} ...", end="", flush=True)

    ebm = joblib.load(model_path)
    feature_names = list(ebm.term_names_)
    n_terms = len(feature_names)

    rows = []
    for i in range(n_terms):
        fname = feature_names[i]
        # Skip interaction terms (contain ' x ' or ' & ')
        if " x " in fname or " & " in fname:
            continue

        scores = np.array(ebm.term_scores_[i], dtype=np.float64).flatten()
        weights = np.array(ebm.bin_weights_[i], dtype=np.float64).flatten()
        abs_scores = np.abs(scores)

        total_w = weights.sum()
        if total_w > 0:
            importance = np.sum(abs_scores * weights) / total_w
        else:
            importance = np.mean(abs_scores)

        rows.append({
            "feature": fname,
            "importance": importance,
            "max_abs_score": abs_scores.max(),
            "score_range": scores.max() - scores.min(),
            "n_bins": len(scores),
            "n_bins_nonzero_wt": int((weights > 0).sum()),
            "total_samples": total_w,
        })

    if verbose:
        print(f" {len(rows)} features")

    return pd.DataFrame(rows)


def infer_tiers_from_importance(summary_df):
    """
    Infer governance tiers from model behavior rather than hardcoding.
    Features with near-zero importance (score_range ≈ 0) were likely zeroed (Red).
    Features with moderate smoothing remain Yellow.
    Stable, high-concordance features are Green.
    Binary features (Sex) are excluded from shape-based tier assignment.

    This is a heuristic for visualization only; actual tier assignments come
    from ShapeQC pipeline output.
    """
    tiers = {}

    # Check if we have a SPAR map or decision matrix CSV alongside
    # (This makes it truly pipeline-driven rather than heuristic)
    return tiers  # Return empty — will try CSV-based lookup first


def load_tiers_from_pipeline(run_dir):
    """
    Try to load tier assignments from ShapeQC pipeline output.
    Searches for Table_4A_Decision_Matrix.csv in ShapeQC output dirs.
    Falls back to SPAR_Map.csv or config-derived tiers.
    """
    # Search for ShapeQC output
    qc_dirs = sorted(glob.glob(os.path.join(run_dir, "*ShapeQC*")))
    for qd in reversed(qc_dirs):
        csv_path = os.path.join(qd, "Table_4A_Decision_Matrix.csv")
        if os.path.exists(csv_path):
            dt = pd.read_csv(csv_path)
            if "feature" in dt.columns and "light" in dt.columns:
                tier_map = dict(zip(dt["feature"], dt["light"]))
                print(f"  [Tiers] Loaded from {os.path.basename(csv_path)}")
                return tier_map

    # Search for SPAR map
    for qd in reversed(qc_dirs):
        spar_files = glob.glob(os.path.join(qd, "*SPAR_Map*.csv"))
        for sf in spar_files:
            dt = pd.read_csv(sf)
            if "feature" in dt.columns and "tier" in dt.columns:
                tier_map = dict(zip(dt["feature"], dt["tier"]))
                print(f"  [Tiers] Loaded from {os.path.basename(sf)}")
                return tier_map

    return None


def assign_tiers_from_config():
    """
    Fallback: derive tiers from the R config file definitions.
    Reads 00_config.R to parse GREEN/YELLOW/RED/GRAY feature lists.
    This ensures consistency with the R pipeline without hardcoding tier values.
    """
    config_paths = [
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "R", "00_config.R"),
    ]

    for config_path in config_paths:
        if not os.path.exists(config_path):
            continue

        tier_map = {}
        current_tier = None
        with open(config_path, "r") as f:
            for line in f:
                line = line.strip()
                if "GREEN_FEATURES" in line and "<-" in line:
                    current_tier = "Green"
                elif "YELLOW_FEATURES" in line and "<-" in line:
                    current_tier = "Yellow"
                elif "RED_FEATURES" in line and "<-" in line:
                    current_tier = "Red"
                elif "GRAY_FEATURES" in line and "<-" in line:
                    current_tier = "Gray"
                elif "ALL_EDIT_FEATURES" in line or "FEATURE_TIER_MAP" in line:
                    current_tier = None

                if current_tier and line.startswith('"'):
                    # Extract feature name from quoted string
                    m = re.match(r'"([^"]+)"', line)
                    if m:
                        tier_map[m.group(1)] = current_tier

        if tier_map:
            print(f"  [Tiers] Loaded from {os.path.basename(config_path)}")
            return tier_map

    return None


# =============================================================================
# Section 4: Visualization
# =============================================================================

def plot_importance_by_tier(summary_df, tier_map, output_path, title_suffix=""):
    """
    Bar chart of feature importance grouped and colored by governance tier.
    Features are ordered: Green → Yellow → Red → Gray, then by importance within tier.
    """
    df = summary_df.copy()

    # Assign tiers
    df["tier"] = df["feature"].map(tier_map).fillna("Unknown")

    # Sort: by tier order, then by importance within tier
    tier_order_map = {t: d["order"] for t, d in TIER_DEFINITIONS.items()}
    tier_order_map["Unknown"] = 99
    df["tier_order"] = df["tier"].map(tier_order_map)
    df = df.sort_values(["tier_order", "mean_importance"], ascending=[True, False])
    df = df.reset_index(drop=True)

    # Build color list
    colors = []
    for _, row in df.iterrows():
        tier = row["tier"]
        if tier in TIER_DEFINITIONS:
            colors.append(TIER_DEFINITIONS[tier]["color"])
        else:
            colors.append("#bdc3c7")

    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.4)))

    bars = ax.barh(
        range(len(df)), df["mean_importance"],
        xerr=df.get("sd_importance", None),
        color=colors, edgecolor="white", linewidth=0.5,
        capsize=3, error_kw={"elinewidth": 1, "capthick": 1}
    )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["feature"], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Mean Weighted Importance (|score| × sample weight)", fontsize=11)
    ax.set_title(
        f"EBM Feature Importance by Governance Tier{title_suffix}",
        fontsize=13, fontweight="bold", pad=15
    )

    # Add tier group separators
    prev_tier = None
    for i, (_, row) in enumerate(df.iterrows()):
        if prev_tier is not None and row["tier"] != prev_tier:
            ax.axhline(y=i - 0.5, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
        prev_tier = row["tier"]

    # Legend
    legend_handles = []
    tiers_present = df["tier"].unique()
    for tier_name in ["Green", "Yellow", "Red", "Gray"]:
        if tier_name in tiers_present:
            td = TIER_DEFINITIONS[tier_name]
            legend_handles.append(
                mpatches.Patch(color=td["color"], label=td["label"])
            )

    ax.legend(handles=legend_handles, loc="lower right", fontsize=9,
              framealpha=0.9, edgecolor="gray")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(output_path)}")


def plot_importance_per_iteration(all_imp_df, tier_map, output_path):
    """
    Grouped bar chart showing per-iteration importance, colored by site.
    Iterations are grouped by site with visual separators.
    """
    df = all_imp_df.copy()
    df["tier"] = df["feature"].map(tier_map).fillna("Unknown")
    tier_order_map = {t: d["order"] for t, d in TIER_DEFINITIONS.items()}
    tier_order_map["Unknown"] = 99
    df["tier_order"] = df["tier"].map(tier_order_map)

    # Get feature order (by mean importance within tier groups)
    feat_order = (
        df.groupby("feature")
        .agg(tier_order=("tier_order", "first"), mean_imp=("importance", "mean"))
        .sort_values(["tier_order", "mean_imp"], ascending=[True, False])
        .index.tolist()
    )

    # Build site palette
    sites = sorted(df["site"].unique())
    site_palette = build_site_palette(sites)

    # Build iteration color map
    iter_colors = {}
    for _, row in df[["model", "site", "seed"]].drop_duplicates().iterrows():
        base, light = site_palette[row["site"]]
        color = base if row["seed"] == 1 else light
        iter_colors[row["model"]] = color

    models = sorted(df["model"].unique())
    n_models = len(models)
    n_features = len(feat_order)

    fig, ax = plt.subplots(figsize=(12, max(6, n_features * 0.5)))

    bar_height = 0.8 / n_models
    for j, model in enumerate(models):
        model_data = df[df["model"] == model].set_index("feature")
        vals = [model_data.loc[f, "importance"] if f in model_data.index else 0
                for f in feat_order]
        y_positions = [i + (j - n_models / 2 + 0.5) * bar_height for i in range(n_features)]
        ax.barh(y_positions, vals, height=bar_height,
                color=iter_colors[model], edgecolor="white", linewidth=0.3,
                label=model)

    ax.set_yticks(range(n_features))
    ax.set_yticklabels(feat_order, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Feature Importance", fontsize=11)
    ax.set_title("Per-Iteration Feature Importance (colored by validation site)",
                 fontsize=13, fontweight="bold", pad=15)

    # Site legend
    legend_handles = []
    for site in sites:
        base, light = site_palette[site]
        legend_handles.append(mpatches.Patch(color=base, label=f"{site} (Seed 1)"))
        legend_handles.append(mpatches.Patch(color=light, label=f"{site} (Seed 2)"))

    ax.legend(handles=legend_handles, loc="lower right", fontsize=8,
              ncol=len(sites), framealpha=0.9, edgecolor="gray")

    # Tier separators
    prev_tier = None
    for i, feat in enumerate(feat_order):
        tier = tier_map.get(feat, "Unknown")
        if prev_tier is not None and tier != prev_tier:
            ax.axhline(y=i - 0.5, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
        prev_tier = tier

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(output_path)}")


def plot_shape_functions(models_info, tier_map, output_dir):
    """
    Overlay shape function plots for all iterations, one subplot per feature.
    Lines colored by site (iteration grouping).
    """
    # Load all models and extract shape data
    print("  Loading models for shape function extraction ...")
    shape_data = {}  # feature -> list of (x_values, scores, meta)

    sites = set()
    for mi in models_info:
        meta = mi["meta"]
        sites.add(meta["site"])
        print(f"    Loading {meta['tag']} ...", end="", flush=True)
        ebm = joblib.load(mi["model_path"])

        for i, fname in enumerate(ebm.term_names_):
            if " x " in fname or " & " in fname:
                continue

            scores = np.array(ebm.term_scores_[i], dtype=np.float64).flatten()

            # Get bin edges for x-axis
            if hasattr(ebm, "bins_") and len(ebm.bins_) > i:
                bins = np.array(ebm.bins_[i], dtype=np.float64).flatten()
                # bins has len = scores - 1 for continuous features
                if len(bins) == len(scores) - 1:
                    x_vals = np.concatenate([[bins[0] - 1], bins])
                elif len(bins) == len(scores):
                    x_vals = bins
                else:
                    x_vals = np.arange(len(scores))
            else:
                x_vals = np.arange(len(scores))

            if fname not in shape_data:
                shape_data[fname] = []
            shape_data[fname].append((x_vals, scores, meta))

        del ebm
        print(" done")

    if not shape_data:
        print("  [WARN] No shape data extracted")
        return

    site_palette = build_site_palette(list(sites))
    features = sorted(shape_data.keys())

    # Determine grid layout
    n_feat = len(features)
    n_cols = 4
    n_rows = (n_feat + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1) if n_cols > 1 else np.array([[axes]])

    for idx, feat in enumerate(features):
        r, c = divmod(idx, n_cols)
        ax = axes[r, c]

        tier = tier_map.get(feat, "Unknown")
        tier_color = TIER_DEFINITIONS.get(tier, {"color": "#bdc3c7"})["color"]

        for x_vals, scores, meta in shape_data[feat]:
            base, light = site_palette[meta["site"]]
            color = base if meta["seed"] == 1 else light
            label = f"{meta['site']} s{meta['seed']}"
            ax.step(x_vals, scores, color=color, alpha=0.7,
                    linewidth=1.5, where="mid", label=label)

        ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="-")
        ax.set_title(feat, fontsize=10, fontweight="bold",
                     color=tier_color if tier != "Unknown" else "black")
        ax.tick_params(labelsize=8)
        ax.set_ylabel("Score", fontsize=8)

        # Add tier badge
        ax.text(0.02, 0.98, tier, transform=ax.transAxes, fontsize=7,
                verticalalignment="top", fontweight="bold",
                color="white",
                bbox=dict(boxstyle="round,pad=0.2", facecolor=tier_color, alpha=0.8))

    # Hide empty subplots
    for idx in range(n_feat, n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r, c].set_visible(False)

    # Build shared legend
    legend_handles = []
    for site in sorted(sites):
        base, light = site_palette[site]
        legend_handles.append(plt.Line2D([0], [0], color=base, lw=2, label=f"{site} Seed 1"))
        legend_handles.append(plt.Line2D([0], [0], color=light, lw=2, label=f"{site} Seed 2"))

    fig.legend(handles=legend_handles, loc="lower center", ncol=len(sites) * 2,
               fontsize=9, framealpha=0.9, edgecolor="gray")

    fig.suptitle("Shape Functions Across IECV Iterations (colored by validation site)",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout(rect=[0, 0.03, 1, 0.99])

    out_path = os.path.join(output_dir, "Shape_Functions_Overlay.png")
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: Shape_Functions_Overlay.png")


# =============================================================================
# Section 5: Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extract EBM feature importance and generate governance-tier visualizations."
    )
    parser.add_argument(
        "--run-dir", type=str, default=None,
        help="Path to the IECV run directory containing Iter*_External_* subdirs. "
             "Auto-discovered if not specified."
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Output directory for figures and CSVs. Defaults to <run_dir>/figures_python/."
    )
    parser.add_argument(
        "--skip-shapes", action="store_true",
        help="Skip shape function overlay plots (faster)."
    )
    args = parser.parse_args()

    # --- Discover RUN_DIR ---
    project_root = os.path.dirname(os.path.abspath(__file__))
    if args.run_dir:
        run_dir = args.run_dir
    else:
        run_dir = discover_run_dir(project_root)
        if run_dir is None:
            # Try environment variable
            env_dir = os.environ.get("EBM_RUN_DIR")
            if env_dir and os.path.isdir(env_dir):
                run_dir = env_dir
            else:
                print("ERROR: Could not auto-discover RUN_DIR.")
                print("  Use --run-dir to specify explicitly,")
                print("  or set EBM_RUN_DIR environment variable.")
                sys.exit(1)

    print("=" * 70)
    print("  EBM Feature Importance Extraction & Visualization")
    print("=" * 70)
    print(f"  RUN_DIR: {run_dir}")

    # --- Discover models ---
    models_info = find_models(run_dir)
    print(f"  Found {len(models_info)} iterations:")
    for mi in models_info:
        m = mi["meta"]
        print(f"    Iter{m['iter']}: External={m['site']}, Seed={m['seed']}")
        print(f"      Model: {os.path.basename(mi['model_path'])}")

    # --- Load tier assignments ---
    tier_map = load_tiers_from_pipeline(run_dir)
    if tier_map is None:
        tier_map = assign_tiers_from_config()
    if tier_map is None:
        print("  [WARN] No tier assignments found; all features will be 'Unknown'")
        tier_map = {}
    else:
        tier_counts = {}
        for t in tier_map.values():
            tier_counts[t] = tier_counts.get(t, 0) + 1
        print(f"  Tiers: {tier_counts}")

    # --- Output directory ---
    output_dir = args.output_dir or os.path.join(run_dir, "figures_python")
    os.makedirs(output_dir, exist_ok=True)
    print(f"  Output: {output_dir}\n")

    # --- Extract importance from all models ---
    print("  Extracting feature importance ...")
    all_imp = []
    for mi in models_info:
        meta = mi["meta"]
        df = extract_importance(mi["model_path"])
        df["model"] = meta["tag"]
        df["site"] = meta["site"]
        df["seed"] = meta["seed"]
        all_imp.append(df)

    all_imp_df = pd.concat(all_imp, ignore_index=True)

    # Normalized share within each model
    for model in all_imp_df["model"].unique():
        mask = all_imp_df["model"] == model
        total = all_imp_df.loc[mask, "importance"].sum()
        all_imp_df.loc[mask, "norm_share"] = all_imp_df.loc[mask, "importance"] / total

    # --- Summary table ---
    summary = all_imp_df.groupby("feature").agg(
        mean_importance=("importance", "mean"),
        sd_importance=("importance", "std"),
        min_importance=("importance", "min"),
        max_importance=("importance", "max"),
        mean_max_abs=("max_abs_score", "mean"),
        mean_range=("score_range", "mean"),
        mean_norm_share=("norm_share", "mean"),
    ).reset_index()

    summary["cv"] = summary["sd_importance"] / summary["mean_importance"]
    summary["consensus_rank"] = summary["mean_importance"].rank(ascending=False, method="min").astype(int)
    summary = summary.sort_values("consensus_rank")

    # Print summary table
    print(f"\n  Total: {len(all_imp_df)} rows ({len(models_info)} models × {all_imp_df['feature'].nunique()} features)")
    print("\n  --- Feature Importance Summary (weighted, official) ---\n")
    print(f"  {'Feature':25s} {'Mean':>8s} {'SD':>8s} {'CV':>6s} {'Share%':>8s} {'Rank':>5s}")
    print(f"  {'-' * 68}")
    for _, r in summary.iterrows():
        print(f"  {r['feature']:25s} {r['mean_importance']:8.4f} {r['sd_importance']:8.4f} "
              f"{r['cv']:6.2f} {r['mean_norm_share']*100:7.1f}% {r['consensus_rank']:5d}")

    top3_share = summary[summary["consensus_rank"] <= 3]["mean_norm_share"].sum() * 100
    print(f"\n  Top-3 normalized share: {top3_share:.1f}%")

    # --- Save CSVs ---
    raw_path = os.path.join(output_dir, "feature_importance_raw.csv")
    summ_path = os.path.join(output_dir, "feature_importance_summary.csv")
    all_imp_df.to_csv(raw_path, index=False)
    summary.to_csv(summ_path, index=False)
    print(f"\n  Saved: {os.path.basename(raw_path)}")
    print(f"  Saved: {os.path.basename(summ_path)}")

    # --- Plot 1: Importance by tier ---
    print("\n  Generating plots ...")
    plot_importance_by_tier(
        summary, tier_map,
        os.path.join(output_dir, "Feature_Importance_by_Tier.png")
    )

    # --- Plot 2: Per-iteration importance ---
    plot_importance_per_iteration(
        all_imp_df, tier_map,
        os.path.join(output_dir, "Feature_Importance_per_Iteration.png")
    )

    # --- Plot 3: Shape function overlays ---
    if not args.skip_shapes:
        plot_shape_functions(models_info, tier_map, output_dir)
    else:
        print("  [Skip] Shape function plots (--skip-shapes)")

    # --- Done ---
    print(f"\n{'=' * 70}")
    print(f"  Done. All outputs saved to: {output_dir}")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
