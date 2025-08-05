# io.py

import pandas as pd
from pathlib import Path
from glob import glob
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Union, Literal
from pybedtools import BedTool
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_rmap_file(rmap_file: str) -> pd.DataFrame:
    """
    Parse a single recombination map file into a standardized DataFrame.

    The chromosome is derived from the filename, assuming it follows the pattern:
    something_chr10_1.txt → chromosome = "chr101"
    """
    rmap_file_path = Path(rmap_file)
    parts = rmap_file_path.name.split("_")
    chromosome = "".join(parts[1:3])  # e.g., ['chr10', '1'] → 'chr101'

    df = pd.read_csv(rmap_file_path, sep="\t", header=None)
    df.insert(0, "CHROM", chromosome)
    df.rename(columns={0: "START", 1: "END", 2: "RHO"}, inplace=True)
    return df

def load_recombination_maps(rmap_files: str) -> List[pd.DataFrame]:
    """
    Load and parse multiple recombination map files in parallel.

    Args:
        rmap_files (str): Glob pattern matching recombination map files.

    Returns:
        List[pd.DataFrame]: List of DataFrames with parsed map data.
    """
    rmap_file_names = glob(rmap_files, recursive=True)
    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dataframes = pool.map(parse_rmap_file, rmap_file_names)
    return dataframes

def save_transformed_data(datasets: List[pd.DataFrame], prefix="transformed") -> None:
    for data in datasets:
        chrom_name = data["chrom"].iat[0]
        window_size = data["end"].iat[0] - data["start"].iat[0]
        filename = f"{prefix}_{chrom_name}_w{window_size}.tsv"
        data.to_csv(filename, sep="\t", header=True, index=False)

def parse_bed_file(file: str) -> pd.DataFrame:
    """
    Parse a BED file into a DataFrame with standardized columns.
    Assumes the first line is a header to skip.
    """
    return pd.read_csv(
        file,
        sep='\t',
        header=None,
        names=['Chromosome', 'Start', 'End', 'Score', 'Mindpoint'],
        skiprows=1  # Adjust if your files don't have headers
    )

def load_bed_files(
    input_files: str,
    return_type: Literal["bed", "dataframe"] = "dataframe") -> Union[BedTool, pd.DataFrame]:
    """
    Load multiple BED files in parallel and return either a merged DataFrame or BedTool object.
    
    Args:
        input_files (str): Glob pattern for input BED files.
        return_type (str): One of 'bed' or 'dataframe'.

    Returns:
        Union[BedTool, pd.DataFrame]
    """
    file_paths = glob(input_files, recursive=True)
    
    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dfs = pool.map(parse_bed_file, file_paths)

    merged_df = (
        pd.concat(dfs, ignore_index=True)
        .sort_values(by=["Chromosome", "Start"])
        .reset_index(drop=True)
    )

    if return_type == "bed":
        return BedTool.from_dataframe(merged_df)
    elif return_type == "dataframe":
        return merged_df
    else:
        raise ValueError(f"Unknown return_type: {return_type}")

def load_fasta_file(fasta_file: str) -> List[SeqRecord]:
    fasta_file_path = Path(fasta_file)
    return list(SeqIO.parse(fasta_file_path, "fasta"))

def parse_fasta_file(
    fasta_file: str,
    mapping: Dict[str, str]) -> Union[int, Exception]:
    
    loaded_fasta = load_fasta_file(fasta_file)
    modified_fasta: List[SeqRecord] = []
    try:
        for record in loaded_fasta:
            if record.id in mapping:
                new_id = mapping[record.id]
                record.id = new_id
            modified_fasta.append(record)
        output_path = Path("./blackcap.fasta")
        return SeqIO.write(modified_fasta, output_path, "fasta")
    except Exception as e:
        print(f"Error renaming fasta file {e}")
        raise
