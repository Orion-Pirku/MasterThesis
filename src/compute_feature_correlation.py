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
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="compute_feature_correlation.py", 
        usage="Compute and Plot the Pearson correlation between the recombination rate and several population genomics statistics"
        )
    parser.add_argument('--rho',type=str, nargs="+", help='recombination rate file', required=True)
    parser.add_argument('--feature-a', type=str, nargs="+", help='first feature files', required=True)
    parser.add_argument('--feature-b', type=str, nargs="+", help='second feature files', required=True)
    parser.add_argument('--feature-c', type=str, nargs="+", help='third feature file', required=True)
    parser.add_argument('--feature-a-name', type=str, help='first feature name', required=True)
    parser.add_argument('--feature-b-name', type=str, help='second feature name', required=True)
    parser.add_argument('--feature-c-name', type=str, help='third feature name', required=True)
    parser.add_argument('--window-size', '-ws', type=int, required=True, help="Window Size by which to interpolate the recombianation rate")
    parser.add_argument('--output-directory', '-o', type=str, required=False)
    return parser.parse_args() 

def main() -> None:
    args: argparse.Namespace = parse_arguments()
    
    bin_size: int = args.window_size    
    pairwise_intersections: dict[str, BedTool] = {}
    if args.output_directory:
        os.makedirs(args.output_directory, exist_ok=True)
    try:
        rec_maps = load_recombination_maps(args.rho)
        loaded_feature_a: list[pd.DataFrame] = load_bed_files(args.feature_a)
        loaded_feature_b = load_bed_files(args.feature_b)
        loaded_feature_c = load_bed_files(args.feature_c)
    except FileNotFoundError as e:  # Usually file not found not file exists
        print(f"File error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error {e}: unknown error occurred")
        sys.exit(1)
        
    sorted_rec_map_dataframes: list[pd.DataFrame] = sort_windows(rec_maps)
    rec_maps_windowed: list[pd.DataFrame] = [make_windows(df, args.window_size) for df in sorted_rec_map_dataframes]
    
    concatenated_rec_maps: pd.DataFrame = concatenate_windows(rec_maps_windowed)
    concatenated_feature_a: pd.DataFrame = preprocess_vc_stats(loaded_feature_a)
    concatenated_feature_b: pd.DataFrame = preprocess_vc_stats(loaded_feature_b)
    concatenated_feature_c: pd.DataFrame = preprocess_vc_stats(loaded_feature_c)
     
    rec_maps_bed: BedTool = BedTool.from_dataframe(concatenated_rec_maps, header=True)
    feature_a_bed: BedTool = BedTool.from_dataframe(concatenated_feature_a, header=True)
    feature_b_bed: BedTool = BedTool.from_dataframe(concatenated_feature_b, header=True)
    feature_c_bed: BedTool = BedTool.from_dataframe(concatenated_feature_c, header=True)     
    
    beds = {
        "rho": rec_maps_bed,
        args.feature_a_name: feature_a_bed,
        args.feature_b_name: feature_b_bed,
        args.feature_c_name: feature_c_bed,
    }
    
    pairwise_intersections = {
        f"{n1}_x_{n2}": b1.intersect(b2, wa=True, wb=True)
        for (n1, b1), (n2, b2) in combinations(beds.items(), 2)
    }
    
    results = {
        name: compute_feature_correlation(bt, score_A_idx=3, score_B_idx=9)
        for name, bt in pairwise_intersections.items()  # or intersections_vs_rho
    }

    labels = ["rho", args.feature_a_name, args.feature_b_name, args.feature_c_name]
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
    
    R_df = pd.DataFrame(R, index=labels, columns=labels)
    P_df = pd.DataFrame(P, index=labels, columns=labels)
    annot_df = P_df.applymap(lambda p: f"p={p:.3f}" if not np.isnan(p) else "")
    
    plt.figure(figsize=(12,8))
    sns.heatmap(
        R_df, vmin=0, vmax=1, 
        cmap='coolwarm', 
        xticklabels=labels, 
        yticklabels=labels, 
        square=True,
        annot=annot_df, fmt="",
        annot_kws={"size": 12, "color": "black"}
        )
    
    plt.xticks(fontsize=18, rotation=45)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
    
