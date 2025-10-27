#!/usr/bin/env python3
import pandas as pd
import sys
import traceback
import argparse
from hotspotter.io import load_recombination_maps
from hotspotter.transform import (
        make_windows,
        sort_windows,
        concatenate_windows
        )

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            prog="clean_raw_rho.py",
            usage="%(prog)s [options]",
            description="CLI Tool to clean and transform raw pyrho output"
            )
    parser.add_argument(
            "-i",
            "--input-files",
            type=str,
            nargs="+",
            help="Input files: raw pyrho output",
            required=True
            )
    
    parser.add_argument(
            "-w",
            "--window-size",
            type=int,
            help="Window size over which the recombination rate is averaged",
            default=10000,
            required=False
            )
    parser.add_argument(
            "-e",
            "--effective-pop-size",
            type=float,
            help="Effective population size to normalize recombination rate." +
            "Pyrho Already includes this in the recombination rate inference" +
            "Needed only if the recombination rate was inferred by LDhelment",
            required=False,
            default=0.0
            ),
    parser.add_argument(
            "-c",
            "--concatenate",
            type=bool,
            required=False,
            default=False
            )
    parser.add_argument(
            "-o",
            "--output-file",
            type=str,
            required=True
            )
    if len(sys.argv) == 1:
        parser.print_help()
    return parser.parse_args()

def main() -> None:
    args = parse_arguments()
    
    raw_rec_data: list[pd.DataFrame] = load_recombination_maps(args.input_files)
    sorted_raw_rec_data: list[pd.DataFrame] = sort_windows(raw_rec_data)
    recombination_windows_list: list[pd.DataFrame] = [make_windows(
        df,
        args.window_size,
        args.effective_pop_size
        ) for df in sorted_raw_rec_data]
    if args.concatenate is True:
        concatenated_recombination_windows: pd.DataFrame = concatenate_windows(recombination_windows_list)
        concatenated_recombination_windows.to_csv(
            args.output_file,
            sep="\t",
            header=True,
            index=False
            )
    else:
        for df in recombination_windows_list:
            chr = df.iloc[0,0]
            df.to_csv(
                    f"{chr}_{args.output_file}",
                    sep='\t',
                    header=True,
                    index=False
                    )
if __name__ == "__main__":
    main()

