#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
from hotspotter.io import load_hotspots_bed_files
from hotspotter.transform import concatenate_windows, sort_windows


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="concatenate_bed_files",
        usage="%(prog)s [options]",
        description="CLI Tool to concantenate individual bed files",
    )
    parser.add_argument(
        "-i",
        "--input-files",
        type=str,
        nargs="+",
        help="Files to be concatenated",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        help="Name of the concatenated file",
        default="output.bed",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def main():
    args = parse_arguments()
    bed_files = load_hotspots_bed_files(args.input_files)
    sorted_bed_files = sort_windows(bed_files)
    concantenated_bed_file = concatenate_windows(sorted_bed_files)
    concantenated_bed_file.to_csv(args.output_file, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
