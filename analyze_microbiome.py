from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


def safe_log2(x: float) -> float:
    return float(np.log2(x))


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini–Hochberg FDR correction.
    Returns q-values (adjusted p-values).
    """
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    if n == 0:
        return pvals

    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(1, n + 1))

    # enforce monotonicity
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    out = np.empty_like(q)
    out[order] = q
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Explore associations between microbiome features and chemotherapy toxicity severity."
    )

    parser.add_argument(
        "--abundance",
        default="data/example_microbiome.csv",
        help="Path to microbiome abundance table (CSV/TSV). Default: data/example_microbiome.csv",
    )
    parser.add_argument(
        "--metadata",
        default="data/example_metadata.csv",
        help="Path to clinical metadata table (CSV/TSV). Default: data/example_metadata.csv",
    )
    parser.add_argument(
        "--sep",
        default=",",
        help="Delimiter for input files. Use ',' for CSV or '\\t' for TSV. Default: ','",
    )
    parser.add_argument(
        "--group_col",
        default="Severity",
        help="Column name in metadata that defines the groups. Default: Severity",
    )
    parser.add_argument(
        "--group1",
        default="Mild",
        help="First group label (baseline). Default: Mild",
    )
    parser.add_argument(
        "--group2",
        default="Severe",
        help="Second group label (compared against group1). Default: Severe",
    )
    parser.add_argument(
        "--out",
        default="results",
        help="Output directory. Default: results",
    )

    # Volcano labeling mode
    parser.add_argument(
        "--label_mode",
        choices=["threshold", "top", "none"],
        default="threshold",
        help="How to label points on the volcano plot: threshold/top/none. Default: threshold",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=10,
        help="If --label_mode top, label the top N features by q_value (FDR). Default: 10",
    )

    # Thresholds for labeling + threshold lines
    parser.add_argument(
        "--use_q",
        action="store_true",
        help="If set, thresholding uses q_value (FDR) instead of raw p_value.",
    )
    parser.add_argument(
        "--p_thresh",
        type=float,
        default=0.05,
        help="Threshold for p_value (used if --use_q is NOT set). Default: 0.05",
    )
    parser.add_argument(
        "--q_thresh",
        type=float,
        default=0.10,
        help="Threshold for q_value/FDR (used if --use_q IS set). Default: 0.10",
    )
    parser.add_argument(
        "--fc_thresh",
        type=float,
        default=1.0,
        help="Label features with |log2FC| > fc_thresh and draw vertical lines at ±fc_thresh. Default: 1.0",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    abundance_df = pd.read_csv(args.abundance, sep=args.sep)
    metadata_df = pd.read_csv(args.metadata, sep=args.sep)

    # --- Validate required columns ---
    if "SampleID" not in abundance_df.columns:
        raise ValueError("Abundance table must contain a 'SampleID' column.")
    if "SampleID" not in metadata_df.columns:
        raise ValueError("Metadata table must contain a 'SampleID' column.")
    if args.group_col not in metadata_df.columns:
        raise ValueError(f"Metadata table must contain '{args.group_col}' column.")

    # --- Check SampleID overlap ---
    abundance_ids = set(abundance_df["SampleID"])
    metadata_ids = set(metadata_df["SampleID"])

    missing_in_abundance = sorted(metadata_ids - abundance_ids)
    missing_in_metadata = sorted(abundance_ids - metadata_ids)

    if missing_in_abundance or missing_in_metadata:
        print("WARNING: SampleID mismatch detected!\n")
        if missing_in_abundance:
            print("Samples present in metadata but missing in abundance:")
            print(missing_in_abundance)
        if missing_in_metadata:
            print("\nSamples present in abundance but missing in metadata:")
            print(missing_in_metadata)
        raise SystemExit("\nFix SampleID mismatch before continuing.")

    merged = pd.merge(metadata_df, abundance_df, on="SampleID", how="inner")

    # --- Define groups ---
    group_col = args.group_col
    group1 = args.group1
    group2 = args.group2

    g1 = merged[merged[group_col] == group1]
    g2 = merged[merged[group_col] == group2]

    if g1.empty or g2.empty:
        raise ValueError(
            f"One of the groups is empty. Found sizes: {group1}={len(g1)}, {group2}={len(g2)}"
        )

    # --- Feature columns (abundance columns only) ---
    feature_cols = [c for c in abundance_df.columns if c != "SampleID"]

    # --- Compute stats per feature ---
    eps = 1e-9
    rows = []

    for feat in feature_cols:
        x1 = g1[feat].astype(float).values
        x2 = g2[feat].astype(float).values

        mean1 = float(np.mean(x1))
        mean2 = float(np.mean(x2))

        fold_change = (mean2 + eps) / (mean1 + eps)
        log2fc = safe_log2(fold_change)

        _, pval = mannwhitneyu(x1, x2, alternative="two-sided")

        rows.append(
            {
                "feature": feat,
                f"mean_{group1}": mean1,
                f"mean_{group2}": mean2,
                f"log2FC_({group2}_vs_{group1})": log2fc,
                "p_value": float(pval),
            }
        )

    results_df = pd.DataFrame(rows)

    # --- Multiple testing correction (BH / FDR) ---
    results_df["q_value"] = benjamini_hochberg(results_df["p_value"].values)

    # Sort after q-values computed
    results_df = results_df.sort_values(["q_value", "p_value"], ascending=True).reset_index(drop=True)

    # --- Save results table ---
    out_csv = out_dir / "results_table.csv"
    results_df.to_csv(out_csv, index=False)

    # --- Volcano plot ---
    log2fc_col = f"log2FC_({group2}_vs_{group1})"
    results_df["neg_log10_p"] = -np.log10(results_df["p_value"] + 1e-12)

    plt.figure(figsize=(7, 5))
    plt.scatter(
        results_df[log2fc_col],
        results_df["neg_log10_p"],
        alpha=0.7,
    )

    # Threshold lines: vertical lines at ±fc_thresh
    plt.axvline(args.fc_thresh, linestyle="--", color="grey")
    plt.axvline(-args.fc_thresh, linestyle="--", color="grey")
    plt.axvline(0, linestyle=":", color="grey")

    # Horizontal line:
    # If using q-value thresholds, we approximate by converting q_thresh to a p-threshold line is not correct.
    # So we draw the horizontal line according to p_thresh (since y-axis is -log10(p)).
    # We still *label* by q if --use_q is set (below).
    plt.axhline(-np.log10(args.p_thresh + 1e-12), linestyle="--", color="grey")

    plt.xlabel(f"log2 Fold Change ({group2} vs {group1})")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano Plot: Microbiome Associations")

    # --- Label points on volcano ---
    if args.label_mode == "top":
        to_label = results_df.nsmallest(args.top_n, "q_value")
    elif args.label_mode == "threshold":
        if args.use_q:
            to_label = results_df[
                (results_df["q_value"] < args.q_thresh)
                & (results_df[log2fc_col].abs() > args.fc_thresh)
            ]
        else:
            to_label = results_df[
                (results_df["p_value"] < args.p_thresh)
                & (results_df[log2fc_col].abs() > args.fc_thresh)
            ]
    else:
        to_label = results_df.iloc[0:0]

    for _, row in to_label.iterrows():
        x = row[log2fc_col]
        y = row["neg_log10_p"]
        name = row["feature"]
        plt.text(x + 0.02, y + 0.02, name, fontsize=8)

    volcano_path = out_dir / "volcano_plot.png"
    plt.tight_layout()
    plt.savefig(volcano_path)
    plt.close()

    # --- Boxplots for top 3 features (by q-value) ---
    top_features = results_df["feature"].head(3).tolist()
    for feat in top_features:
        plt.figure(figsize=(4, 4))
        merged.boxplot(column=feat, by=group_col)
        plt.title(f"{feat} abundance by {group_col}")
        plt.suptitle("")
        plt.ylabel("Abundance")

        boxplot_path = out_dir / f"boxplot_{feat}.png"
        plt.tight_layout()
        plt.savefig(boxplot_path)
        plt.close()

    # --- Summary prints ---
    print("Analysis complete ✅")
    print(f"Samples: {group1}={len(g1)}, {group2}={len(g2)}")
    print(f"Features analyzed: {len(feature_cols)}")
    print(f"Saved results table: {out_csv}")
    print(f"Saved volcano plot: {volcano_path}")
    print(f"Labeled {len(to_label)} features on volcano plot (mode: {args.label_mode}).")
    print("Saved boxplots for:", ", ".join(top_features))
    print("\nTop results:")
    print(results_df.head(10))


if __name__ == "__main__":
    main()