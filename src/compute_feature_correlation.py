#!/usr/bin/env python3
import argparse
import sys
import pyranges as pr

# your package imports
from hotspotter.io import load_and_prepare_feature, load_recombination_maps
from hotspotter.plotting import plot_correlation_matrix
from hotspotter.recombination_rate_analyzer import (
    compute_feature_correlation,
    compute_intersections,
    fill_correlation_matrices,
)
from hotspotter.transform import (
    concatenate_windows,
    make_windows,
    sort_windows
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="compute_feature_correlation.py",
        description="Compute & plot Pearson correlation between recombination"
        + "rate and any number of population-genomics features.",
    )
    parser.add_argument(
        "--rho",
        type=str,
        nargs="+",
        required=True,
        help="Recombination map BED(s). Globs allowed.",
    )
    parser.add_argument(
        "--rho-name",
        type=str,
        required=True,
        help='Display name for the recombination rate layer (e.g., "rho").',
    )
    # Repeated --feature flags: first token = label, rest = BED files
    parser.add_argument(
        "--feature",
        action="append",
        nargs="+",
        metavar=("LABEL", "BED..."),
        help="Repeatable. Example: --feature tajima-d tajima_chr1.bed tajima_chr2.bed",
    )
    parser.add_argument(
        "--window-size",
        "-ws",
        type=int,
        required=True,
        help="Window size in bp for recombination map windowing.",
    )
    parser.add_argument(
        "--output-directory",
        "-o",
        type=str,
        default=".",
        help="Directory to write outputs (default: current dir).",
    )
    parser.add_argument(
        "--genome-sizes", type=str, help="File containing chromosome sizes"
    )
    args = parser.parse_args()

    # Require at least one feature
    if not args.feature or len(args.feature) == 0:
        parser.error("Provide at least one --feature <label> <beds...>")
        sys.exit(1)
    return args


def main() -> None:
    args = parse_arguments()
    outdir = args.output_directory
    ws = args.window_size

    try:
        # Recombination maps (windowed)
        rec_maps = load_recombination_maps(args.rho)
        sorted_rec_map_dataframes = sort_windows(rec_maps)
        rec_maps_windowed = [make_windows(df, ws) for df in sorted_rec_map_dataframes]
        concatenated_rec_maps = concatenate_windows(rec_maps_windowed)
        concatenated_rec_maps = concatenated_rec_maps.rename(
            columns={"CHROM": "Chromosome", "START": "Start", "END": "End"}
        )
        rec_maps_bed = pr.PyRanges(concatenated_rec_maps)

        # Generic features
        beds: dict[str, pr.PyRanges] = {args.rho_name: rec_maps_bed}
        for feat_spec in args.feature:
            label, *files = feat_spec
            if len(files) == 0:
                raise ValueError(f"--feature {label} has no files")
            beds[label] = load_and_prepare_feature(
                label,
                files,
                window_size=args.window_size,
                genome_sizes_file=args.genome_sizes,
            )

    except FileNotFoundError as e:
        print(f"[ERROR] file not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] {type(e).__name__}: {e}")
        sys.exit(1)

    pairwise_intersections = compute_intersections(beds)
    results = {
        name: compute_feature_correlation(bt)
        for name, bt in pairwise_intersections.items()
    }

    labels = list(beds.keys())
    R_df, P_df = fill_correlation_matrices(results, labels)
    plot_correlation_matrix(R_df, P_df, labels, outdir, window_size=args.window_size)


if __name__ == "__main__":
    main()
