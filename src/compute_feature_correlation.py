#!/usr/bin/env python3
import argparse
import sys, os
import numpy as np
import numpy.typing as npt
from scipy.fftpack import sc_diff
from hotspotter.io import load_bed_files, load_recombination_maps
from hotspotter.transform import make_windows, sort_windows, concatenate_windows, preprocess_vc_stats
from hotspotter.recombination_rate_analyzer import compute_feature_correlation
from pathlib import Path
import pandas as pd
from pybedtools import BedTool
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="compute_feature_correlation.py", 
        usage="Compute and Plot the Pearson correlation between the recombination rate and several population genomics statistics"
    )
    parser.add_argument('--rho', type=str, nargs="+", required=True)
    parser.add_argument('--feature-a', type=str, nargs="+", required=True)
    parser.add_argument('--feature-b', type=str, nargs="+", required=True)
    parser.add_argument('--feature-c', type=str, nargs="+", required=True)
    parser.add_argument('--rho-name', type=str, required=True)
    parser.add_argument('--feature-a-name', type=str, required=True)
    parser.add_argument('--feature-b-name', type=str, required=True)
    parser.add_argument('--feature-c-name', type=str, required=True)
    parser.add_argument('--window-size', '-ws', type=int, required=True)
    parser.add_argument('--output-directory', '-o', type=str, required=False)
    return parser.parse_args()


def compute_intersections(beds: dict[str, BedTool]) -> dict[str, BedTool]:
    return {
        f"{n1}_x_{n2}": b1.intersect(b2, wa=True, wb=True)
        for (n1, b1), (n2, b2) in combinations(beds.items(), 2)
    }


def fill_correlation_matrices(results: dict[str, tuple[float, float]], labels: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    idx = {name: i for i, name in enumerate(labels)}
    R = np.full((len(labels), len(labels)), np.nan, dtype=np.float64)
    P = np.full((len(labels), len(labels)), np.nan, dtype=np.float64)
    np.fill_diagonal(R, 1.0)
    np.fill_diagonal(P, 0.0000)

    for key, r in results.items():
        n1, n2 = key.split("_x_")
        i, j = idx[n1], idx[n2]
        if i > j:
            R[i, j], P[i, j] = r[0], r[1]
        else:
            R[j, i], P[j, i] = r[0], r[1]

    return pd.DataFrame(R, index=labels, columns=labels), pd.DataFrame(P, index=labels, columns=labels)


def plot_correlation_matrix(R_df: pd.DataFrame, P_df: pd.DataFrame, labels: list[str], output_dir: str) -> None:
    # Format p-values in scientific notation
    f = np.vectorize(lambda p: f"p:{p:.1e}" if not np.isnan(p) else "")
    annot_df = pd.DataFrame(f(P_df.values), index=P_df.index, columns=P_df.columns)
    
    plt.figure(figsize=(8, 6), dpi=200)
    ax = sns.heatmap(
        R_df, vmin=0, vmax=1,
        cmap='coolwarm',
        xticklabels=labels,
        yticklabels=labels,
        square=True,
        annot=annot_df, fmt="",
        annot_kws={"size": 12},
        cbar_kws={"label": "Pearson Correlation"}  # ✅ only valid param here
    )

    # ✅ Customize colorbar label font
    cbar = ax.collections[0].colorbar
    cbar.set_label("Pearson Correlation", size=12, weight='bold', labelpad=25)
    plt.xticks(fontsize=14, rotation=45, fontweight='bold')
    plt.yticks(fontsize=14, rotation=0, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/feature_correlation_plot.png", format="png")


def main() -> None:
    args = parse_arguments()
    bin_size = args.window_size

    if args.output_directory:
        os.makedirs(args.output_directory, exist_ok=True)

    try:
        rec_maps = load_recombination_maps(args.rho)
        loaded_feature_a = load_bed_files(args.feature_a)
        loaded_feature_b = load_bed_files(args.feature_b)
        loaded_feature_c = load_bed_files(args.feature_c)
    except FileNotFoundError as e:
        print(f"File error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error {e}: unknown error occurred")
        sys.exit(1)

    sorted_rec_map_dataframes = sort_windows(rec_maps)
    rec_maps_windowed = [make_windows(df, args.window_size) for df in sorted_rec_map_dataframes]

    concatenated_rec_maps = concatenate_windows(rec_maps_windowed)
    concatenated_feature_a = preprocess_vc_stats(loaded_feature_a)
    concatenated_feature_b = preprocess_vc_stats(loaded_feature_b)
    concatenated_feature_c = preprocess_vc_stats(loaded_feature_c)

    rec_maps_bed = BedTool.from_dataframe(concatenated_rec_maps, header=True)
    feature_a_bed = BedTool.from_dataframe(concatenated_feature_a, header=True)
    feature_b_bed = BedTool.from_dataframe(concatenated_feature_b, header=True)
    feature_c_bed = BedTool.from_dataframe(concatenated_feature_c, header=True)

    beds = {
        args.rho_name: rec_maps_bed,
        args.feature_a_name: feature_a_bed,
        args.feature_b_name: feature_b_bed,
        args.feature_c_name: feature_c_bed,
    }

    pairwise_intersections = compute_intersections(beds)
    
    results = {
        name: compute_feature_correlation(bt, score_A_idx=3, score_B_idx=9)
        for name, bt in pairwise_intersections.items()
    }

    labels = list(beds.keys())
    R_df, P_df = fill_correlation_matrices(results, labels)
    plot_correlation_matrix(R_df, P_df, labels, args.output_directory)


if __name__ == "__main__":
    main()

