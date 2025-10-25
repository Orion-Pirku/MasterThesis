#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
from hotspotter.io import load_recombination_maps
from hotspotter.transform import (
    make_windows,
    sort_windows,
    concatenate_windows
)
from hotspotter.plotting import (
    plot_genome_wide_rho,
    plot_score_distribution,
    plot_score_strength_per_chrom
)

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="plot_rec_rate.py",
        description=(
            "CLI tool to plot the distribution of the recombination rate "
            "over the genome."
        ),
        usage="%(prog)s -i <input_files> -w <window_size> -o <output_dir>"

    )
    parser.add_argument(
        "-i",
        "--input-files",
        type=str,
        nargs="+",
        help="input raw recombination rate files inferred from pyrho"
        )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        help="Window size over which to normalize recombination rate",
        required=True
    )
    parser.add_argument(
        '-f',
        '--feature',
        type=str,
        required=True
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        type=str,
        required=True,
        help="Directory where plots will be saved."
    )


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    try:
        loaded_bed_files: list[pd.DataFrame] = load_recombination_maps(
            args.input_files
        )
    except FileNotFoundError as e:
        print (f"Error {e}: Files {args.input_files} not found!")
        sys.exit(1)
    except IsADirectoryError as e:
        print (f"Error {e}: Files are a directory!")
        sys.exit(1)
    except Exception as e:
        print (f"Unknown error {e} occured!")
        sys.exit(1)

    sorted_rec_rate: list[pd.DataFrame] = sort_windows(loaded_bed_files)
    rec_rate_windows: list[pd.DataFrame] = [make_windows(df,args.window_size)for df in sorted_rec_rate]
    conc_rec_rate_windows: pd.DataFrame = concatenate_windows(rec_rate_windows)

    plot_genome_wide_rho(rec_rate_windows, f"{args.output_directory}/genome_wide_{args.feature}.png")
    plot_score_distribution(rec_rate_windows, f"{args.output_directory}/{args.feature}_distribution.png")
    plot_score_strength_per_chrom(conc_rec_rate_windows, f"{args.output_directory}/{args.feature}strength_per_chrom.png")


if __name__ == "__main__":
    main()
