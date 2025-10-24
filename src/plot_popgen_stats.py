#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import seaborn
from glob import glob
import pandas as pd
import argparse
import sys
from hotspotter.io import load_and_prepare_feature, load_recombination_maps
from hotspotter.plotting import plot_pop_gen_stats
from hotspotter.transform import (
    preprocess_popgen_stats,
    concatenate_windows,
    sort_windows,
    compute_gene_density,
    make_windows
)
import pyranges as pr

def _check_file_type(input_files: list[str]) -> None:
    allowed_suffixes: tuple[str, ...] = (
        ".pi",
        ".windowpi",
        ".D",
        ".bed",
        ".tajimaD",
        ".snpden",
        ".snpdens"
    )
    if not any(file.endswith(allowed_suffixes) for file in input_files):
        print(f"Error: Only files ending with {allowed_suffixes} are allowed.")
        sys.exit(1)


def parse_arguments():
    plot_stats = argparse.ArgumentParser(
        prog="Plot Pop-Gen Stats",
        usage='%(prog)s [options]', 
        description="Population Genomics Statistics Plotter"
    )

    plot_stats.add_argument(
            '--feature',
            action='append',
            nargs='+',
            metavar=('LABEL', 'BED...'),
            help='Repeatable. Example: --feature tajima-d tajima_chr1.bed tajima_chr2.bed'
            )
    plot_stats.add_argument(
        '-o',
        '--output-file', 
        required=False, 
        type=str,
        default="plot.png",
        help='Name of output file (plot), default: plot.png'
    )
    plot_stats.add_argument(
        '-g',
        '--genome-sizes',
        type=str,
        help='Tab delimited file containing chromosomes and their sizes'
    )
    plot_stats.add_argument(
        '-w',
        '--window-size',
        type=int,
        required=True,
        default=100_000,
    )
    plot_stats.add_argument(
        '--output-file-format', 
        required=False, 
        type=str, 
        choices=["png", "svg", "jpeg", "pdf"],
        default="png",
        help="Format of output figure, default: png"
    )
    plot_stats.add_argument(
        '-y',
        '--y-axis-title',
        type=str,
        required=False,
        default="Value",
        help='Title of the Y-axis'
    )
    plot_stats.add_argument(
        '-t',
        '--figure-title',
        type=str,
        required=False,
        default="",
        help='Title of the Whole figure'
    )
    plot_stats.add_argument(
        '-c',
        '--plot-color',
        required=False,
        type=str,
        help='Color of the plot. default: black',
        default="black"
    )
    if len(sys.argv) == 1:
        plot_stats.print_help()
        sys.exit(1)
    return plot_stats.parse_args()


def main():
    try:
        args = parse_arguments()
    except Exception as e:
        print(f"Could not parse args, error: {e}")
        sys.exit(1)

    try:
        features: dict[str, pr.PyRanges | pd.DataFrame] = {}
        for feature_label in args.feature:
            label, *files = feature_label
            if len(files) == 0:
                raise ValueError(f"--feature {label} has no files")
            if label in any("rho", "cM/Mb", "rec-rate"):
                rec_rates = load_recombination_maps(files)
                rec_rates = [make_windows(df, args.window_size) for df in rec_rates]
                sorted_rec_rates = sort_windows(rec_rates)
                concat_rec_rates = concatenate_windows(sorted_rec_rates)
                features[label] = concat_rec_rates
            features[label] = load_and_prepare_feature(
                label, 
                files, 
                window_size=args.window_size,
                genome_sizes_file=args.genome_sizes
            ).df

    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
        sys.exit(1)

    except IsADirectoryError as e:
        print(f"Expected a file but found a directory: {e.filename}")
        sys.exit(1)

    except PermissionError as e:
        print(f"Permission denied: {e.filename}")
        sys.exit(1)

    except ValueError as e:
        print(f"Invalid input: {e}")
        sys.exit(1)

    except TypeError as e:
        print(f"Type error in input: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
    try:
        for label, df in features.items():
            plot_pop_gen_stats(
                data_frame=df,
                y_axis_title=label, 
                plotLineColor=args.plot_color,
                outFileName=f"{label}_{args.output_file}",
                outFileFormat=args.output_file_format,
                title=args.figure_title
                )
    except Exception as e:
        print(f"Could not generate plot! Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()