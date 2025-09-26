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
