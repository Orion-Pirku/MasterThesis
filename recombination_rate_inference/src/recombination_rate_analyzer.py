from pybedtools import BedTool
import pandas as pd
import numpy as np
from typing import Union, Mapping, List, Tuple
from multiprocessing import Pool, cpu_count
import pyranges as pr
import numpy.typing as npt
from scipy.stats import pearsonr
import scipy.signal as ss

def get_score_strength_per_chrom(
    input_windows: Union[BedTool, pd.DataFrame]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits the recombination rate per chromosome on three quantiles based on its strength

    Parameters
    -------------
    input_windows: 
    Genomic windows of BedTool or DataFrame format whose 4th column contains recombination rate or score

    Returns
    -------------
    Tuple[pd.DataFrame, pd.DataFrame]:
    sorted Score strength per chromosome and min max values per quantile

    """

    if not isinstance(input_windows, pd.DataFrame):
        genomic_windows = input_windows.to_dataframe(disable_auto_names=True, 
                                        names=["chrom", "start", "end", "score(cM/Mb)", "midpoint"])
    else:
        genomic_windows = input_windows.iloc[:, :5].copy()
        genomic_windows.columns = ["chrom", "start", "end", "score(cM/Mb)", "midpoint"]
        
    genomic_windows["score_strength"] = pd.qcut(genomic_windows["score(cM/Mb)"], q=3, labels=["low", "medium", "high"])
    min_max_per_label = genomic_windows.groupby("score_strength")["score(cM/Mb)"].agg(['min', 'max']).reset_index()
    min_max_per_label.to_csv("Recombination_Rate_Strength.csv", sep = "\t") 
    counts = genomic_windows.groupby(["chrom", "score_strength"]).size().reset_index(name="count")
    total_per_chrom = genomic_windows.groupby("chrom").size().reset_index(name="total")
    merged = pd.merge(counts, total_per_chrom, on="chrom")
    merged['percent'] = (merged['count'] / merged['total']) * 100
    chrom_order = [f"chr{i}" for i in range(1, 34)]
    merged["chrom"] = pd.Categorical(merged["chrom"], categories=chrom_order, ordered=True)
    return merged.sort_values(['chrom', 'score_strength']), min_max_per_label


def compute_genomic_feature_correlation(
    bed_feature: BedTool,
    score_one_idx: int,
    score_two_idx: int
    ) -> Tuple[float, float]:
    """
    Computes the Pearson correlation between two feature score columns 
    from a BedTool intersect result.

    Parameters
    ----------
    bed_feature : BedTool
        A BedTool object, resulting from BedTool.intersect(a, b, wa=True, wb=True).
    score_one_idx : int
        Column index for the score found in bed_feature coming from the first BED input.
    score_two_idx : int
        Column index for the score found in bed_feature coming from the second BED input.

    Returns
    -------
    Tuple[float, float]
        Pearson correlation coefficient and two-tailed p-value.
    """

    df = bed_feature.to_dataframe(names=None)

    if df.empty:
        raise ValueError("The provided BedTool object resulted in an empty DataFrame.")

    if score_one_idx >= df.shape[1] or score_two_idx >= df.shape[1]:
        raise IndexError("Score index out of bounds")

    feature_one_array = df.iloc[:, score_one_idx].astype(float).to_numpy()
    feature_two_array = df.iloc[:, score_two_idx].astype(float).to_numpy()

    return pearsonr(
        x=feature_one_array,
        y=feature_two_array,
        alternative="two-sided"
    )

def call_recombination_hotspots(
    input_bed: Union[BedTool, pd.DataFrame], 
    score_column_idx: int
    ) -> npt.NDArray[np.int32]:
    
    
