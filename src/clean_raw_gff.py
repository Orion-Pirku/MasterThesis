#!/usr/bin/env python3    
import argparse
from pathlib import Path
import pandas as pd
import sys
import re
import os
from hotspotter.transform import clean_raw_gff
from hotspotter.io import load_accession_mapping, save_per_chromosome


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='clean_gff_data',
        usage='%(prog)s [options]',
        description='Tool to clean raw GFF.'
    )

    parser.add_argument('-gff', required=True, help='Path to the GFF file')
    parser.add_argument('--genome-size', required=True, help='Path to genome sizes file. Contains (Chrom - Chrom Size)')
    parser.add_argument('--accession-map', required=True, help='Path to genome accession file')
    parser.add_argument('-feature', default='gene', help='Feature to extract (e.g., gene, exon)')
    parser.add_argument('-output', required=True, default="genes.gff", help='Name of the output file')
    
    if len(sys.argv) == 1:
        parser.print_help()
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    os.makedirs(args.output, exist_ok=True)

    for file_arg in [args.gff, args.genome_size, args.accession_map]:
        if not Path(file_arg).is_file():
            print(f"[ERROR] File does not exist: {file_arg}")
            sys.exit(1)

    accession_map = load_accession_mapping(args.accession_map)
    df = clean_raw_gff(args.gff, args.feature, accession_map)

    try:
        df.to_csv(f"{args.output}/bSylAtri_{args.feature}.bed", sep="\t", index=False, header=True)
        save_per_chromosome(df, args.output, args.feature)
        print(f"[INFO] Cleaned data saved to {Path(args.output).resolve()}")
    except Exception as e:
        print(f"[ERROR] Failed to write output: {e}")
        sys.exit(1)
        
if __name__ == "__main__":
    main()
