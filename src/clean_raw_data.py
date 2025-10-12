from json import load
from shlex import join
from hotspotter.io import load_bed_files, load_recombination_maps
from hotspotter.transform import sort_windows, concatenate_windows
import argparse
from pathlib import Path
import pandas as pd
import sys

def preprocess_rho(args):

    rec_maps: list[pd.DataFrame] = load_recombination_maps(args.input_files)
    sorted_rec_maps: list[pd.DataFrame] = sort_windows(rec_maps)
    joined_rec_maps: pd.DataFrame = concatenate_windows(sorted_rec_maps)
    joined_rec_maps.to_csv(args.output, sep="\t", header=True, index=False, chunksize=1000)

def preprocess_pop_gen_stats(args):
    pop_gen_files: list[pd.DataFrame] = load_bed_files(args.input_files)
    sorted_pop_gen: list[pd.DataFrame] = sort_windows(pop_gen_files)
    
    for df in sorted_pop_gen:
        col_names: pd.Index = df.columns
        # Case-insensitive check for 'end' column
        if "end" not in df.columns.str.lower().tolist() and "bin_end" not in df.columns.str.lower().tolist():
            # Calculate window sizes (difference between next start and current start)
            window_size: pd.Series = df[col_names[1]].shift(-1) - df[col_names[1]]            
            # Fill last window size with the previous window size (to avoid NaN)
            window_size.iloc[-1] = window_size.iloc[-2]
            window_size = window_size.astype(int)
            # Create new BIN_END column (end = start + window_size - 1)
            bin_end: pd.Series = df[col_names[1]].astype(int) + window_size - 1
            # Insert BIN_END at position 2 (3rd column)
            df.insert(loc=2, column="BIN_END", value=bin_end)
    joined_pop_gen: pd.DataFrame = concatenate_windows(sorted_pop_gen)
    joined_pop_gen.to_csv(args.output, sep="\t", header=True, index=False)
        
def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clean_raw_data", 
        usage="A CLI tools written in python to pre-process the vcfstats and pyrho outputs"
        )
    parser.add_argument(
        "--verbose", 
        "-v", 
        action="count", 
        default=0, 
        required=False, 
        help='verbosity of logging program for debug purposes'
    )
    
    sub = parser.add_subparsers(dest="command")
    rho_parser = sub.add_parser("rho", help="pre-processing pyrho output files")
    rho_parser.add_argument("--input-files", '-i', nargs="+", type=Path, help="Path to input files", required=True)
    rho_parser.add_argument("--output", "-o", type=str, help="Name of the output file", required=True)
    rho_parser.add_argument("--processes", "-p", default=1, type=int, help="Number of processes for concurrent file parsing. Default is 1")
    rho_parser.set_defaults(func=preprocess_rho)
    
    tajimaD_parser = sub.add_parser("pop_stats", help=r"pre-processing vcf_stats outputs such as Window $\pi$, Tajima D, snp density")
    tajimaD_parser.add_argument("--input-files", '-i', nargs="+", type=Path, help="Path to input files", required=True)
    tajimaD_parser.add_argument("--output", "-o", type=str, help="Name of the output file", required=True)
    tajimaD_parser.add_argument("--processes", "-p", default=1, type=int, help="Number of processes for concurrent file parsing. Default is 1")
    tajimaD_parser.set_defaults(func=preprocess_pop_gen_stats)
    return parser


def main() -> None:
    parser = argument_parser()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)
    
if __name__ == "__main__":
    main()