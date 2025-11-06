from pandas.core import window
from pybedtools import BedTool
import pandas as pd
import numpy as np
from typing import Union, Mapping, List, Tuple
from multiprocessing import Pool, cpu_count
import pyranges as pr
import numpy.typing as npt
from scipy.stats import pearsonr
import scipy.signal as ss
import scipy.ndimage as sn
import pyranges as pr
from itertools import combinations
import pywt


def compute_intersections(beds: dict[str, pr.PyRanges]) -> dict[str, pr.PyRanges]:
    return {
        f"{n1}_x_{n2}": b1.join(b2)
        for (n1, b1), (n2, b2) in combinations(beds.items(), 2)
    }


def fill_correlation_matrices(
    results: dict[str, tuple[float, float]], labels: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    idx = {name: i for i, name in enumerate(labels)}
    R = np.full((len(labels), len(labels)), np.nan, dtype=np.float64)
    P = np.full((len(labels), len(labels)), np.nan, dtype=np.float64)
    np.fill_diagonal(R, 1.0)
    np.fill_diagonal(P, 0.0)

    for key, (rval, pval) in results.items():
        n1, n2 = key.split("_x_")
        i, j = idx[n1], idx[n2]
        if i > j:
            R[i, j], P[i, j] = rval, pval
        else:
            R[j, i], P[j, i] = rval, pval
    return (
        pd.DataFrame(R, index=labels, columns=labels),
        pd.DataFrame(P, index=labels, columns=labels),
    )


def get_score_strength_per_chrom(
    input_windows: Union[BedTool, pd.DataFrame],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not isinstance(input_windows, pd.DataFrame):
        genomic_windows = input_windows.to_dataframe(
            disable_auto_names=True,
            names=["chrom", "start", "end", "score(cM/Mb)", "midpoint"],
        )
    else:
        genomic_windows = input_windows.iloc[:, :5].copy()
        genomic_windows.columns = ["chrom", "start", "end", "score(cM/Mb)", "midpoint"]

    genomic_windows["score_strength"] = pd.qcut(
        genomic_windows["score(cM/Mb)"], q=3, labels=["low", "medium", "high"]
    )
    min_max_per_label = (
        genomic_windows.groupby("score_strength")["score(cM/Mb)"]
        .agg(["min", "max"])
        .reset_index()
    )
    min_max_per_label.to_csv("Recombination_Rate_Strength.csv", sep="\t")
    counts = (
        genomic_windows.groupby(["chrom", "score_strength"])
        .size()
        .reset_index(name="count")
    )
    total_per_chrom = genomic_windows.groupby("chrom").size().reset_index(name="total")
    merged = pd.merge(counts, total_per_chrom, on="chrom")
    merged["percent"] = (merged["count"] / merged["total"]) * 100
    chrom_order = [f"chr{i}" for i in range(1, 34)]
    merged["chrom"] = pd.Categorical(
        merged["chrom"], categories=chrom_order, ordered=True
    )
    return merged.sort_values(["chrom", "score_strength"]), min_max_per_label


def compute_feature_correlation(
    intersection_object: pr.PyRanges | pd.DataFrame,
) -> tuple[float, float]:
    if isinstance(intersection_object, pr.PyRanges):
        intersection_df = intersection_object.df
    elif isinstance(intersection_object, pd.DataFrame):
        intersection_df = intersection_object
    else:
        raise TypeError("intersection_object must be a BedTool or a DataFrame")

    if intersection_df.empty:
        raise ValueError("The provided BedTool object resulted in an empty DataFrame.")

    float_cols = intersection_df.select_dtypes(include=["float64"]).columns[:2]
    print(float_cols)
    if len(float_cols) < 2:
        raise ValueError(f"Need ≥2 numeric columns, found {list(float_cols)}")
    a = pd.to_numeric(intersection_df[float_cols[0]], errors="coerce")
    b = pd.to_numeric(intersection_df[float_cols[1]], errors="coerce")

    mask = a.notna() & b.notna()  # keep only rows where both are numeric

    feature_one_array = a[mask].to_numpy(dtype=float)
    feature_two_array = b[mask].to_numpy(dtype=float)

    return pearsonr(x=feature_one_array, y=feature_two_array, alternative="two-sided")


def call_hotspots_cwt(
    raw_signal: npt.NDArray[np.float64], prominence_cutoff: float
) -> tuple[
    npt.NDArray[np.int64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]:
    scales = np.unique(np.round(np.logspace(np.log10(50), np.log10(600), 200))).astype(
        int
    )

    smoothed_signal = sn.gaussian_filter1d(raw_signal, sigma=5)

    coeff_matrix, _ = pywt.cwt(smoothed_signal, scales, "mexh", method="fft")
    positive_coeff = np.maximum(coeff_matrix, 0.0)
    normalized_coeff = (positive_coeff - positive_coeff.mean(axis=1, keepdims=True)) / (
        positive_coeff.std(axis=1, keepdims=True) + 1e-9
    )
    strong_mask = normalized_coeff > 2.0
    persistence_ratio = strong_mask.mean(axis=0)

    persistence_gate = (persistence_ratio >= 0.3).astype(float)
    score = positive_coeff.max(axis=0) * persistence_gate

    positive_scores = score[score > 0]
    prominence_threshold = (
        np.percentile(positive_scores, prominence_cutoff)
        if positive_scores.size
        else 0.0
    )

    peaks, _ = ss.find_peaks(
        score,
        prominence=prominence_threshold,
    )

    return peaks.astype(np.int64), smoothed_signal, coeff_matrix, scales
