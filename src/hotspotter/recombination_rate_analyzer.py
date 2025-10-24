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


def compute_intersections(
        beds: dict[str, pr.PyRanges]
        ) -> dict[str, pr.PyRanges]:
    return {
        f"{n1}_x_{n2}": b1.join(b2)
        for (n1, b1), (n2, b2) in combinations(beds.items(), 2)
    }


def fill_correlation_matrices(
        results: dict[str, tuple[float, float]], 
        labels: list[str]
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
            pd.DataFrame(
                R,
                index=labels,
                columns=labels),
            pd.DataFrame(
                P,
                index=labels,
                columns=labels
            )
        )



def get_score_strength_per_chrom(
    input_windows: Union[BedTool, pd.DataFrame]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:

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


def compute_feature_correlation(
    intersection_object: pr.PyRanges | pd.DataFrame) -> tuple[float, float]:

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


    return pearsonr(
        x=feature_one_array,
        y=feature_two_array,
        alternative="two-sided"
    )

def _compute_optimal_threshold(
    data: npt.NDArray[np.float64],
    k_min: float,
    k_max: float,
    sigma: float) -> float:
    # Build k grid (ensure at least two points)
    k = np.arange(k_min, k_max, 1.0, dtype=float)
    if k.size < 2:
        # fall back to a single candidate
        k = np.array([k_min], dtype=float)

    n_peaks = np.empty(k.shape, dtype=np.int64)

    # Count peaks for each candidate
    for idx, ki in enumerate(k):
        # Tune wlen/distance as you like; ensure they are <= len(data) where applicable
        peaks, _ = ss.find_peaks(data, prominence=ki * sigma, wlen=50_001, distance=10_000)
        n_peaks[idx] = peaks.size

    # Normalize x in [0,1]; guard zero range
    denom_x = (k.max() - k.min())
    x = (k - k.min()) / denom_x if denom_x != 0 else np.zeros_like(k, dtype=float)

    # Normalize y in [0,1]; guard zero range
    denom_y = (n_peaks.max() - n_peaks.min())
    y = (n_peaks - n_peaks.min()) / denom_y if denom_y != 0 else np.zeros_like(n_peaks, dtype=float)

    # Line from first to last point
    x0, y0 = x[0], y[0]
    x1, y1 = x[-1], y[-1]

    # Perpendicular distance of each (x,y) to the first-last line (a standard "kneedle"/L-method style)
    num = np.abs((y1 - y0) * x - (x1 - x0) * y + x1 * y0 - y1 * x0)
    den = np.hypot(y1 - y0, x1 - x0)

    # If all points are identical, just pick the: first k
    if den == 0:
        return float(k[0])

    dist = num / den
    idx_star = int(np.argmax(dist))
    return float(k[idx_star])

def _calculate_robust_sigma(
    smoothed_signal: npt.NDArray[np.float64], 
    window_size: int) -> float:
    
    if window_size % 2 == 0:
        raise ValueError("Error: window_size must be an odd number")
    baseline: npt.NDArray[np.float64] = sn.median_filter(smoothed_signal, size=window_size, mode='reflect')
    contrast: npt.NDArray[np.float64] = smoothed_signal - baseline
    mad = np.median(np.abs(contrast - np.median(contrast))) + 1e-12
    sigma = 1.4826 * mad
    return float(sigma)

def call_hotspots(raw_signal: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64], dict[str, npt.NDArray[np.float64]]]:
    smoothed_signal: npt.NDArray[np.float64] = ss.savgol_filter(raw_signal, window_length=2001, polyorder=1)
    robust_sigma = float(_calculate_robust_sigma(smoothed_signal, window_size=50_001))
    optimal_threshold = float(_compute_optimal_threshold(smoothed_signal, k_min=1, k_max=30, sigma=robust_sigma))
    peaks, properties = ss.find_peaks(
        smoothed_signal, 
        prominence=optimal_threshold*robust_sigma, 
        wlen=50_001
        )
    peaks = peaks.astype(np.int64, copy=False)
    properties = {k: np.asarray(v, dtype=np.float64) for k, v in properties.items()}
    return smoothed_signal, peaks, properties
