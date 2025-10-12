#!/usr/bin/env python3

import argparse
from argparse import ArgumentParser, Namespace
import numpy as np
import glob
from typing import  List, Optional

def compute_table_average(input_files: List[str], output_name: str = "output_table.txt") -> Optional[np.ndarray]:
    arrays = []

    # Read the header from the first file
    header = None
    try:
        with open(input_files[0], 'r') as f:
            header = f.readline().strip()  # Read the first line (header)
    except Exception as e:
        print(f"Failed to read header from {input_files[0]}: {e}")
        return None

    # Load each file (excluding the header)
    for file in input_files:
        try:
            arr = np.loadtxt(file, delimiter="\t", skiprows=1)  # Skip the header row
            arrays.append(arr)
        except Exception as e:
            print(f"Failed to load {file}: {e}")

    # If no valid files are loaded
    if not arrays:
        print("No valid input files loaded.")
        return None

    # Compute the average of the arrays
    try:
        average_table = np.mean(arrays, axis=0)

        # Save the header and the computed average table to the output file
        with open(output_name, 'w') as f:
            if header:  # Save the header
                f.write(header + "\n")
            # Save the averaged data
            np.savetxt(f, average_table, delimiter="\t", fmt="%.6f")

        return average_table
    except Exception as e:
        print(f"Failed to compute or save average: {e}")
        return None

def parse_args() -> Namespace:
    
    parser: ArgumentParser= argparse.ArgumentParser(
        prog="merge_msmc2_tables.py",
        description="Joins multiple msmc2 population size history tables\ninto one by calculating the element-wise average of those tables")    
    parser.add_argument(
        "-i", "--input_directory",
        type=str,
        required=True,
        help="input the name of the directory with wildcard * for file names e.g., files/*.txt. Dont forget quotation marks!")
    
    parser.add_argument(
        "-o", "--output_name",
        type=str,
        required=False,
        default="output_msmc2_table.txt",
        help="input name of output file. Default name is: output_msmc2_table.txt. Dont forget quotation marks around the name"
    )


    return parser.parse_args()

def main():
    try:
        args = parse_args()
        input_files = glob.glob(args.input_directory)
    except Exception as e:
        print(f"Error parsing arguments or finding files: {e}")
        return

    if not input_files:
        print("No files matched the input pattern.")
        return

    compute_table_average(input_files, args.output_name)

if __name__ == "__main__":
    main()


