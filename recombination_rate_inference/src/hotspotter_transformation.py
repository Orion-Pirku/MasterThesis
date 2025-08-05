import numpy as np
import pandas as pd
import re
import pyranges as pr
from pybedtools import BedTool
from typing import Literal, Union, List, Dict, overload, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def rename_fasta_header(
    fasta: List[SeqRecord],
    mapping: Dict[str, str]) -> List[SeqRecord]:

    modified_fasta: List[SeqRecord] = []
    try:
        for record in fasta:
            if record.id in mapping:
                new_id = mapping[record.id]
                record.id = new_id
            modified_fasta.append(record)
        return modified_fasta
    except Exception as e:
        print(f"Error renaming fasta file {e}")
        raise

def transform_feature_table(
    gff_object: Union[pd.DataFrame, BedTool, pr.PyRanges],
    return_type: Literal["bed", "dataframe", "pyranges"] = "bed",
    chr_pattern: str = r"chr([1-9]|[12][0-9]|3[0-3])$") -> Any:
    if isinstance(gff_object, BedTool):
        gff = gff_object.to_dataframe()
    elif isinstance(gff_object, pd.DataFrame):
        gff = gff_object
    elif isinstance(gff_object, pr.PyRanges):
        gff = gff_object.df 
    else:
        raise TypeError("Input must be a pandas DataFrame or a BedTool object")

    required_columns: List[str] = ["chromosome", "start", "end", "# feature", "name", "symbol"]
    missing = set(required_columns) - set(gff.columns)
    
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df_transformed = (
        gff
        .filter(items=required_columns)
        .assign(
            start=lambda df: df["start"].astype("int64"),
            end=lambda df: df["end"].astype("int64"),
            chromosome=lambda df: "chr" + df["chromosome"].astype(str)
        )
        .loc[lambda df: df["chromosome"].str.match(chr_pattern)]
    )
  
    # Rename 'chromosome' to 'chrom' because BedTool expects 'chrom' as the first column
    df_transformed = df_transformed.rename(
        columns={"chromosome": "Chromosome",
                 "start": "Start",
                 "end": "End",
                 "# feature": "Feature",
                 "symbol": "Symbol",
                 "name": "Name"})

    # Reorder columns so that chrom, start, end are first 3 columns (required for BedTool)
    bedtool_columns = ["Chromosome", "Start", "End", "Feature", "Name","Symbol"]
    df_transformed = df_transformed[bedtool_columns]
    df_sorted = df_transformed.sort_values(by=["Chromosome", "Start"]).reset_index(drop=True)

    match return_type:
        case "bed":
            return BedTool.from_dataframe(df_sorted)
        case "dataframe":
            return df_sorted
        case "pyranges":
            return pr.PyRanges(df_sorted)
        case _:
            raise ValueError(f"Unknown return_type: {return_type}")


def split_bedtool_by_feature(
    input_bed: pr.PyRanges | pd.DataFrame | BedTool,
    return_type: Literal["bed", "dataframe", "pyranges"]) -> dict[str, Union[pr.PyRanges, BedTool, pd.DataFrame]]:
    
    if isinstance(input_bed, pr.PyRanges):
        df = input_bed.df
    elif isinstance(input_bed, pd.DataFrame):
        df = input_bed.copy()
    elif isinstance(input_bed, BedTool):
        df = input_bed.to_dataframe()
    else:
        raise TypeError("Input must be a PyRanges, pandas DataFrame, or BedTool object")
    
    if df.shape[1] < 4:
        raise ValueError("Expected at least 4 columns in the BedTool Data Format (chrom, start, end, features)")

    feature_dict: Dict[str, Union[BedTool, pd.DataFrame, pr.PyRanges]] = {}

    unique_features = df.iloc[:, 3].unique()
    for feature in unique_features:
        df_feature = df[df.iloc[:, 3] == feature]
        if return_type == "bed":
            feature_dict[feature] = BedTool.from_dataframe(df_feature)
        elif return_type == "dataframe":
            feature_dict[feature] = df_feature
        elif return_type == "pyranges":
            feature_dict[feature] = pr.PyRanges(df_feature)
        else:
            raise ValueError(f"Unknown return_type: {return_type}")

    return feature_dict

def make_windows(df: pd.DataFrame, window_size: int) -> pd.DataFrame:
    raw = df.to_numpy()
    midpoints = (raw[:, 1].astype(int) + raw[:, 2].astype(int)) // 2
    bin_ids = midpoints // window_size

    unique_bins, inv_indices = np.unique(bin_ids, return_inverse=True)
    sums = np.zeros(len(unique_bins))
    counts = np.zeros(len(unique_bins))

    np.add.at(sums, inv_indices, raw[:, 3].astype(float))
    np.add.at(counts, inv_indices, 1)

    means = (sums / counts) * 100 * 1e6  # scaling factor

    starts = unique_bins * window_size
    ends = starts + window_size
    midpoints = (starts + ends) // 2
    chrom = raw[0, 0]

    result = np.column_stack([
        np.full(len(unique_bins), chrom),
        starts,
        ends,
        means,
        midpoints
    ])

    return pd.DataFrame(result, columns=["Chrom", "Start", "End", "Score(cM/Mb)", "Midpoint"]).astype({
        "Chrom": str,
        "Start": int,
        "End": int,
        "cM/Mb": float,
        "Midpoint": int
    })


def chr_sort_key(df: pd.DataFrame) -> int | float:
    chrom = str(df.iloc[0]["Chrom"])
    match = re.match(r"chr(\d+)", chrom)
    if match:
        return int(match.group(1))
    else:
        return float("inf")  # push non-numeric chroms (e.g., chrX) to end


def sort_windows(windows: List[pd.DataFrame]) -> List[pd.DataFrame]:
    windows.sort(key=chr_sort_key)
    return windows

def concatenate_windows(windows: List[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(windows, axis=0).reset_index(drop=True)

def create_chromosome_mapping(genome_sizes: pd.DataFrame) -> Dict[str, str]:
    genome_sizes["chrom"] = [f"chr{i+1}" for i in range(len(genome_sizes))]
    return dict(zip(genome_sizes.iloc[:, 0], genome_sizes.iloc[:, 2]))    