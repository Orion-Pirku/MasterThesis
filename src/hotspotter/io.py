# io.py
from tkinter.tix import InputOnly
import pandas as pd
from pathlib import Path
from glob import glob
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Union, Literal
from pybedtools import BedTool
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy.typing as npt
import numpy as np
import pandas as pd
from regex import D
import pyranges as pr
from .transform import (
    concatenate_windows, 
    compute_gene_density,
    preprocess_popgen_stats
)
import sys

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


def _classify_feature(label: str) -> str | None:
    label_lower = label.lower().replace("_", " ")
    if "gene density" in label_lower or "genes" in label_lower:
        return "gene density"
    if "gc content" in label_lower:
        return "gc content"
    return next((k for k in ("snp density", "tajima d", "window pi") if k in label_lower), None)


def _normalize_paths(
    paths: list[str] | str
) -> list[str]:
    if isinstance(paths, str):
        return [paths]
    return paths


def load_and_prepare_feature(
    label: str,
    paths: list[str] | str,
    genome_sizes_file: str,
    window_size: int,
) -> pr.PyRanges:

    feature_kind = _classify_feature(label)
    if feature_kind is None:
        raise ValueError(f"Unrecognized feature label: {label!r}")

    file_paths = _normalize_paths(paths)
    if not file_paths:
        raise ValueError("No input paths provided.")

    # Load all files into a list
    loaded_dfs: list[pd.DataFrame]
    if len(file_paths) > 1:
        # prefer your project’s loader if available
        try:
            loaded_dfs = load_bed_files(file_paths)  # type: ignore[name-defined]
        except NameError:
            loaded_dfs = [pd.read_csv(p, sep="\t", header=0) for p in file_paths]
    else:
        loaded_dfs = [pd.read_csv(file_paths[0], sep="\t", header=0)]

    multiple = len(loaded_dfs) > 1

    # ---- Dispatch per feature, concatenating ONLY if multiple files ----
    if feature_kind == "gene density":
        df_in = concatenate_windows(loaded_dfs) if multiple else loaded_dfs[0]  # type: ignore[name-defined]
        processed = compute_gene_density(  # type: ignore[name-defined]
            df_in,
            genome_sizes=genome_sizes_file,
            window_size=window_size
        )

    elif feature_kind == "gc content":
        processed = concatenate_windows(loaded_dfs) if multiple else loaded_dfs[0]  # type: ignore[name-defined]

    elif feature_kind in ("snp density", "tajima d", "window pi"):
        preprocessed_list = preprocess_popgen_stats(  # type: ignore[name-defined]
            feature_dataframes=loaded_dfs
        )
        processed = (
            concatenate_windows(preprocessed_list)  # type: ignore[name-defined]
            if len(preprocessed_list) > 1
            else preprocessed_list[0]
        )

    else:
        # should be unreachable
        raise RuntimeError(f"Unhandled feature kind: {feature_kind!r}")

    # Normalize expected column names (no-op if already correct)
    processed = processed.rename(
        columns={"CHROM": "Chromosome", "START": "Start", "END": "End"},
        errors="ignore",
    )

    # Validate required columns exist before building PyRanges
    required = {"Chromosome", "Start", "End"}
    missing = required.difference(processed.columns)
    if missing:
        raise ValueError(
            f"Processed feature missing columns: {sorted(missing)}"
        )

    return pr.PyRanges(df=processed)


def save_hotspots_as_bed(
    chrom_name: str, 
    genomic_midpoint: npt.NDArray[np.int64],
    smoothed_signal: npt.NDArray[np.float64],
    output_file_name: str,
    hotspot_index: npt.NDArray[np.int64]) -> None:
    start: npt.NDArray[np.int64] = genomic_midpoint - 50
    end: npt.NDArray[np.int64] = genomic_midpoint + 50
    bed_object: pd.DataFrame = pd.DataFrame(
        {
            'Chrom': [chrom_name] * len(hotspot_index), 
            'Start': start[hotspot_index],
            'End': end[hotspot_index], 
            'cM/Mb': smoothed_signal[hotspot_index]
        }).astype(
            {
                'Chrom': str,
                'Start': int,
                'End': int,
                'cM/Mb': float
                
            })
    return bed_object.to_csv(output_file_name, sep='\t', index=False)


def _parse_rmap_file(rmap_file: str) -> pd.DataFrame:
    """
    Parse a single recombination map file into a standardized DataFrame.

    The chromosome is derived from the filename, assuming it follows the pattern:
    something_chr10_1.txt → chromosome = "chr101"
    """
    rmap_file_path = Path(rmap_file)
    parts = rmap_file_path.name.split(".")
    chromosome = "".join(parts[1:3])  # e.g., ['chr10', '1'] → 'chr101'

    df = pd.read_csv(rmap_file_path, sep="\t", header=None)
    df.insert(0, "CHROM", chromosome)
    df.rename(columns={0: "START", 1: "END", 2: "RHO"}, inplace=True)
    return df

def load_recombination_maps(rmap_files: str | list[str]) -> list[pd.DataFrame]:
    """
    Load and parse multiple recombination map files in parallel.

    Args:
        rmap_files (str or list[str]): Glob pattern matching recombination map files.

    Returns:
        List[pd.DataFrame]: List of DataFrames with parsed map data.
    """
    
    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dataframes = pool.map(_parse_rmap_file, rmap_files)
    return dataframes

def save_transformed_data(datasets: List[pd.DataFrame], prefix="transformed") -> None:
    for data in datasets:
        chrom_name = data["chrom"].iat[0]
        window_size = data["end"].iat[0] - data["start"].iat[0]
        filename = f"{prefix}_{chrom_name}_w{window_size}.tsv"
        data.to_csv(filename, sep="\t", header=True, index=False)

def _parse_bed_file(file: str) -> pd.DataFrame:
    """
    Parse a BED file into a DataFrame with standardized columns.
    Assumes the first line is a header to skip.
    """
    file_path = Path(file)
    data_frame = pd.read_csv(
            file_path,
            sep='\t',
            header=0
        )
    if data_frame.iloc[:,0].str.contains("_").any():
        data_frame.iloc[:,0] = data_frame.iloc[:,0].str.replace("_","", regex=False)
        # Add 'END' column if missing

    return data_frame

def load_bed_files(input_files: str | list[str]) -> list[pd.DataFrame]:
    
    if isinstance(input_files, list):  # already expanded by shell
        matched_files = [str(p) for p in input_files if Path(p).exists()]
    else:
        raise TypeError("input_files must be a str or list of str objects.")

    if not matched_files:
        raise FileNotFoundError(f"No BED files found for pattern(s): {input_files}")

    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dfs = pool.map(_parse_bed_file, matched_files)
    return dfs

def load_fasta_file(fasta_file: str) -> List[SeqRecord]:
    fasta_file_path = Path(fasta_file)
    return list(SeqIO.parse(fasta_file_path, "fasta"))

def parse_and_rename_fasta_file(
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
