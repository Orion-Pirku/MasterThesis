from networkx import maximal_independent_set
import numpy as np
import pandas as pd
import re
import pyranges as pr
from pybedtools import BedTool
from typing import Literal, Union, List, Dict, overload, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from hotspotter.transform_utils import (
    _convert_to_dataframe,
    _validate_required_columns,
    _transform_dataframe,
    _convert_to_output,
    _split_by_feature,
    _compute_midpoint
)
import numpy.typing as npt


def _extract_feature_name(attribute: str) -> str | None:
    match = re.search(r"Name=([^;]+)", attribute)
    return match.group(1) if match else None 

def clean_raw_gff(file_path: str, genomic_feature: str) -> pd.DataFrame:
    gff: pd.DataFrame = pd.read_csv(
            file_path,
            sep="\t",
            comment="#",
            header=None,
            names=["SEQ_ID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRIBUTES"]
            )
    feature = gff[gff["TYPE"] == genomic_feature.lower()].copy()
    feature["NAME"] = feature["ATTRIBUTES"].astype(str).apply(_extract_feature_name)
    return feature[["SEQ_ID", "START", "END", "SCORE", "TYPE", "NAME"]]


def compute_gene_density(
        gene_dataframe: pd.DataFrame, 
        window_size: int, 
        genome_sizes: str | pd.DataFrame
        ) -> pd.DataFrame:
    if isinstance(genome_sizes, str):
        genome_df: pd.DataFrame = pd.read_csv(
            genome_sizes, sep=r"\s+", header=None, names=["CHROM", "SIZE"]
        )
    else:
        genome_df = genome_sizes.copy()
    
    valid_chroms = set(gene_dataframe["CHROM"].unique())
    genome_df = genome_df[genome_df["CHROM"].isin(valid_chroms)]

    if genome_df.empty:
        raise ValueError(
            "No overlapping chromosomes between gene_dataframe and genome_sizes!"
        )
    # --- Build fixed windows per chromosome ---
    window_dfs: list[pd.DataFrame] = []
    for _, row in genome_df.iterrows():
        chrom: str = str(row["CHROM"])
        length: int = int(row["SIZE"])
        starts = list(range(0, length, window_size))
        ends = [min(s + window_size, length) for s in starts]
        window_dfs.append(
            pd.DataFrame({
                "Chromosome": chrom,
                "Start": starts,
                "End": ends
            })
        )

    windows_df: pd.DataFrame = pd.concat(window_dfs, ignore_index=True)
    windows: pr.PyRanges = pr.PyRanges(windows_df)
    if not {"CHROM", "START", "END"}.issubset(gene_dataframe.columns):
        raise ValueError(
            "gene_dataframe must contain columns: CHROM, START, END"
            )

    gene_coords = gene_dataframe[["CHROM", "START", "END"]].copy()
    gene_coords = gene_coords.rename(columns={
        "CHROM": "Chromosome",
        "START": "Start",
        "END": "End"
    })
    genes: pr.PyRanges = pr.PyRanges(df=gene_coords)     

    overlap_counts: pr.PyRanges = windows.count_overlaps(genes)
    overlap_df: pd.DataFrame = overlap_counts.df.rename(
        columns={"NumberOverlaps": "GENE_COUNT"}
        )

    # --- Compute density ---
    overlap_df["WINDOW_SIZE"] = overlap_df["End"] - overlap_df["Start"]
    overlap_df["GENE_DENSITY"] = overlap_df["GENE_COUNT"] / overlap_df["WINDOW_SIZE"]

    # --- Standardize columns ---
    return overlap_df[["Chromosome", "Start", "End", "GENE_COUNT", "GENE_DENSITY"]]


def preprocess_popgen_stats(feature_dataframes: list[pd.DataFrame]) -> list[pd.DataFrame]:
    sorted_bed_object: list[pd.DataFrame] = sort_windows(feature_dataframes)
    
    for df in sorted_bed_object:
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
            df.insert(2, "BIN_END", bin_end)

        df = df.rename(columns={"CHROM": "Chromosome","BIN_START": "Start", "BIN_END": "End"}, inplace=True)
            
    return sorted_bed_object



def interpolate_recombination_rate(
        sorted_recombination_maps: list[pd.DataFrame] | pd.DataFrame
        ) -> tuple[list[npt.NDArray[np.float64]], list[npt.NDArray[np.int64]]]:

    if isinstance(sorted_recombination_maps, pd.DataFrame):
        sorted_recombination_maps = [sorted_recombination_maps]
    
    recombination_maps: list[pd.DataFrame] = [_compute_midpoint(df) for df in sorted_recombination_maps]
    recombination_rates: list[npt.NDArray[np.float64]] = [df['Rho'].to_numpy(dtype=float) for df in recombination_maps]
    midpoints: list[npt.NDArray[np.int64]] = [df['Midpoint'].to_numpy(dtype=int) for df in recombination_maps]
    genomic_bins: list[npt.NDArray[np.int64]] = [np.arange(mid.min(), mid.max() + 50, 50) for mid in midpoints]
    binned_rec_rates = [
        np.interp(bins, mid, rho)
        for bins, mid, rho in zip(genomic_bins, midpoints, recombination_rates)
    ] 
    binned_rec_rates = [rec_rate * 1e8 for rec_rate in binned_rec_rates]
    return binned_rec_rates, genomic_bins

def transform_gff_table(
    gff_object: Union[pd.DataFrame, BedTool, pr.PyRanges],
    chr_pattern: str = r"chr([1-9]|[12][0-9]|3[0-3])$",
    return_type: Literal["bed", "dataframe", "pyranges"] = "dataframe" 
    ) -> Union[BedTool, pd.DataFrame, pr.PyRanges]:
    
    df = _convert_to_dataframe(gff_object)
    _validate_required_columns(df)
    df = _transform_dataframe(df, chr_pattern)
    return _convert_to_output(df, return_type)

def split_bedtool_by_feature(
    input_bed: pr.PyRanges | pd.DataFrame | BedTool,
    return_type: Literal["bed", "dataframe", "pyranges"]
    ) -> dict[str, Union[pr.PyRanges, BedTool, pd.DataFrame]]:
    
    df = _convert_to_dataframe(input_bed) 
    if df.shape[1] < 4:
        raise ValueError("Expected at least 4 columns in the BedTool Data Format (chrom, start, end, features)")

    feature_dict: Dict[str, Union[BedTool, pd.DataFrame, pr.PyRanges]] = {}

    return _split_by_feature(df, return_type) 

def make_windows(df: pd.DataFrame, window_size: int) -> pd.DataFrame:
    raw = df.to_numpy()
    start_raw = raw[:,1].astype(int)
    end_raw = raw[:,2].astype(int)
    
    midpoints = (start_raw + end_raw) // 2
    
    max_end = end_raw.max()
    num_bins = (max_end//window_size) + 1
    
    bins = np.arange(0, num_bins *  window_size + 1, window_size)
    
    bin_indeces = np.digitize(midpoints, bins) - 1
    
    
    unique_bins = np.arange(len(bins)-1)
    sums = np.zeros(len(unique_bins))
    counts = np.zeros(len(unique_bins))

    np.add.at(sums, bin_indeces, raw[:, 3].astype(float))
    np.add.at(counts, bin_indeces, 1)

    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts!=0)
    means = means * 100 * 1e6

    starts = bins[:-1]
    ends = bins[1:] - 1
    midpoints = (starts + ends) //2
    chrom = raw[0, 0]

    result = np.column_stack([
        np.full(len(unique_bins), chrom),
        starts,
        ends,
        means,
        midpoints
    ])

    return pd.DataFrame(result, columns=["CHROM", "START", "END", "cM/Mb", "MIDPOINT"]).astype({
        "CHROM": str,
        "START": int,
        "END": int,
        "cM/Mb": float,
        "MIDPOINT": int
    })


def chr_sort_key(df: pd.DataFrame) -> int | float:
    chrom = str(df.iloc[0]["CHROM"])
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
