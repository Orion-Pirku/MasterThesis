#!/usr/bin/env python3    
import argparse
from pathlib import Path
import pandas as pd
import sys
import re
import os

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='clean_gff_data',
        usage='%(prog)s [options]',
        description='Tool to clean raw GFF.'
    )

    parser.add_argument('--gff', required=True, help='Path to the GFF file')
    parser.add_argument('--genome-size', required=True, help='Path to genome sizes file. Contains (Chrom - Chrom Size)')
    parser.add_argument('--accession-map', required=True, help='Path to genome accession file')
    parser.add_argument('--feature', default='gene', help='Feature to extract (e.g., gene, exon)')
    parser.add_argument('--output', required=True, default="genes.gff", help='Name of the output file')
    
    if len(sys.argv) == 1:
        parser.print_help()
    return parser.parse_args()

def _extract_feature_name(attribute) -> str | None:
    match = re.search(r"Name=([^;]+)", str(attribute))
    return match.group(1) if match else None

def load_accession_mapping(file: str) -> dict[str, int]:
    try:
        dataframe = pd.read_csv(file, sep='\t', header=None, names=["CHROM_NUM", "ACCESSION"])
        dataframe["CHROM_NUM"] = "chr"+dataframe["CHROM_NUM"].astype(str)
        return dict(zip(dataframe["ACCESSION"].astype(str), dataframe["CHROM_NUM"].astype(str)))
    except Exception as e:
        print(f"[ERROR] Failed to load accession mapping: {e}")
        sys.exit(1)
        
def clean_raw_gff(
    file_path: str, 
    genomic_feature: str,
    genome_accession: dict[str, str]) -> pd.DataFrame:
    gff: pd.DataFrame = pd.read_csv(
            file_path, 
            sep="\t", 
            comment="#", 
            header=None,
            names=["ACCESSION", "ASSEMBLY", "TYPE", "START", "END", "SCORE", "STRAND", "DOT", "SEQ_ID"]
            )
    feature = gff[gff["TYPE"].str.lower() == genomic_feature.lower()].copy()
    feature["NAME"] = feature["SEQ_ID"].apply(_extract_feature_name)
    feature["CHROM"] = feature["ACCESSION"].astype(str).map(genome_accession)
    feature = feature.dropna(subset=["CHROM"])
    feature['CHROM'] = feature["CHROM"].astype(str)
    feature["START"] = feature["START"].astype(int)
    feature["END"] = feature["END"].astype(int)
    return feature[["CHROM", "START", "END", "STRAND", "TYPE", "NAME"]]

def save_per_chromosome(dataframe: pd.DataFrame, output_path: str, feature_name: str) -> None: 
    for chrom, sub in dataframe.groupby("CHROM"):
        chrom_str = "NA" if pd.isna(chrom) else str(chrom)
        chrom_number = re.sub(r'[a-zA-Z_-]+', "", chrom_str)
        file_path = os.path.join(output_path, f"chr_{chrom_number}_{feature_name}.bed")
        sub.to_csv(file_path, sep='\t', index=False)

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
        df.to_csv(f"{args.output}/bSylAtri_genes.bed", sep="\t", index=False, header=True)
        save_per_chromosome(df, args.output, args.feature)
        print(f"[INFO] Cleaned data saved to {Path(args.output).resolve()}")
    except Exception as e:
        print(f"[ERROR] Failed to write output: {e}")
        sys.exit(1)
        
if __name__ == "__main__":
    main()
