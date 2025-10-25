import numpy as np
import pandas as pd
import re
import pyranges as pr
from pybedtools import BedTool
from typing import Literal, Union, List, Dict, overload, Any, Sequence, Optional

def _convert_to_dataframe(gff_object: Union[pd.DataFrame, BedTool, pr.PyRanges]) -> pd.DataFrame:
    if isinstance(gff_object, BedTool):
        return gff_object.to_dataframe()
    elif isinstance(gff_object, pd.DataFrame):
        return gff_object.copy()
    elif isinstance(gff_object, pr.PyRanges):
        return gff_object.df
    else:
        raise TypeError("Input must be a pandas DataFrame, BedTool, or PyRanges object")


def _validate_required_columns(df: pd.DataFrame):
    required = {"chromosome", "start", "end", "# feature", "name", "symbol"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _transform_dataframe(df: pd.DataFrame, chr_pattern: str) -> pd.DataFrame:
    df = (
        df[["chromosome", "start", "end", "# feature", "name", "symbol"]]
        .assign(
            start=lambda d: d["start"].astype("int64"),
            end=lambda d: d["end"].astype("int64"),
            chromosome=lambda d: "chr" + d["chromosome"].astype(str),
        )
        .loc[lambda d: d["chromosome"].str.match(chr_pattern)]
        .rename(columns={
            "chromosome": "Chromosome",
            "start": "Start",
            "end": "End",
            "# feature": "Feature",
            "name": "Name",
            "symbol": "Symbol"
        })
        .sort_values(by=["Chromosome", "Start"])
        .reset_index(drop=True)
    )
    return df


def _convert_to_output(df: pd.DataFrame, return_type: str) -> Union[BedTool, pd.DataFrame, pr.PyRanges]:
    columns_order = ["Chromosome", "Start", "End", "Feature", "Name", "Symbol"]
    df = df[columns_order]

    match return_type:
        case "bed":
            return BedTool.from_dataframe(df)
        case "dataframe":
            return df
        case "pyranges":
            return pr.PyRanges(df)
        case _:
            raise ValueError(f"Unsupported return_type: {return_type}")

def _split_by_feature(
    df: pd.DataFrame,
    return_type: str
) -> Dict[str, Union[BedTool, pd.DataFrame, pr.PyRanges]]:
    feature_dict: Dict[str, Union[BedTool, pr.PyRanges, pd.DataFrame]] = {}
    feature_col = df.columns[3]  # 4th column: feature

    for feature in df[feature_col].unique():
        df_subset = df[df[feature_col] == feature]

        match return_type:
            case "bed":
                feature_dict[feature] = BedTool.from_dataframe(df_subset)
            case "dataframe":
                feature_dict[feature] = df_subset
            case "pyranges":
                feature_dict[feature] = pr.PyRanges(df_subset)
            case _:
                raise ValueError(f"Unknown return_type: {return_type}")
    
    return feature_dict


def _compute_midpoint(DataFrame: pd.DataFrame) -> pd.DataFrame:
    DataFrame['Midpoint'] = (DataFrame['Start'] + DataFrame['End']) // 2
    return DataFrame

def _to_dataframe_any(input_windows: pd.DataFrame | pr.PyRanges) -> pd.DataFrame:
    if isinstance(input_windows, pd.DataFrame):
        df = input_windows.copy()
    elif isinstance(input_windows, pr.PyRanges):
        df = input_windows.df.copy()
    else:
        raise TypeError("input must be a pandas DataFrame or a PyRanges")
    if list(df.columns) == list(range(len(df.columns))):
        df.columns = [f"col{i}" for i in range(len(df.columns))]
    return df

# ----------------------------- Detection -----------------------------

def _detect_group_col(df: pd.DataFrame, group_col: Optional[str] = None) -> str:
    """
    Heuristically pick a 'group' column (chromosome-like) without using names:
    1) prefer non-numeric (object/string/category) with many repeats (categorical-like)
    2) else pick numeric/categorical column with lowest unique ratio (distinct/rows)
    """
    if group_col is not None:
        return group_col

    n = len(df)
    if n == 0:
        return df.columns[0]

    def unique_ratio(s: pd.Series) -> float:
        try:
            return s.nunique(dropna=True) / max(1, len(s))
        except Exception:
            return 1.0

    # Candidate sets
    non_numeric = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        # Choose the one with smallest unique ratio (most repeated labels)
        col = min(non_numeric, key=lambda c: unique_ratio(df[c]))
        return col

    # Fallback: among numeric, choose the most categorical-like (smallest unique ratio)
    col = min(df.columns, key=lambda c: unique_ratio(df[c]))
    return col


def _detect_score_col(df: pd.DataFrame, score_col: Optional[str] = None) -> str:
    """
    Pick a numeric 'score' column purely from data:
    - prefer float dtype
    - prefer higher variability (std > 0) and higher unique ratio
    - if none numeric, try coercing object columns to numeric and re-evaluate
    """
    if score_col is not None:
        return score_col

    n = len(df)
    if n == 0:
        return df.columns[0]

    def score_metric(s: pd.Series) -> tuple:
        # Higher is better: (is_float, unique_ratio, std_not_zero)
        is_float = int(pd.api.types.is_float_dtype(s))
        try:
            s_num = pd.to_numeric(s, errors="coerce")
        except Exception:
            s_num = pd.Series([np.nan] * len(s))
        ur = s_num.nunique(dropna=True) / max(1, len(s_num))
        std_ok = int(np.nanstd(s_num) > 0)
        return (is_float, ur, std_ok)

    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if numeric_cols:
        # Pick the column with the best metric
        best = max(numeric_cols, key=lambda c: score_metric(df[c]))
        if np.nanstd(pd.to_numeric(df[best], errors="coerce")) > 0:
            return best

    # Try coercing non-numeric to numeric
    coerce_candidates = []
    for c in df.columns:
        s_num = pd.to_numeric(df[c], errors="coerce")
        if s_num.notna().any():
            coerce_candidates.append((c, s_num))

    if coerce_candidates:
        # Evaluate metric on coerced series
        def metric_on_series(s_num: pd.Series) -> tuple:
            is_float = 1  # coerced numeric → treat as float-like
            ur = s_num.nunique(dropna=True) / max(1, len(s_num))
            std_ok = int(np.nanstd(s_num) > 0)
            return (is_float, ur, std_ok)

        best_c, best_s = max(coerce_candidates, key=lambda x: metric_on_series(x[1]))
        if np.nanstd(best_s) > 0:
            return best_c

    # Last resort
    return df.columns[0]


# ----------------------------- Binning -----------------------------

def _bin_scores_quantile(
    s: pd.Series,
    n_bins: int = 3,
    labels: Optional[Sequence[str]] = None,
) -> tuple[pd.Categorical, pd.DataFrame]:
    """
    Quantile-bins a numeric series robustly (handles ties/constant data).
    Returns (categorical bins, ranges_df).
    """
    x = pd.to_numeric(s, errors="coerce")
    x = x.dropna()
    if x.empty:
        cat = pd.Categorical(np.repeat(np.nan, len(s)), categories=[], ordered=True)
        ranges = pd.DataFrame(columns=["score_strength", "min", "max"])
        return cat, ranges

    # Try qcut with duplicates dropped; if too few unique values, fall back to cut
    try:
        binned, edges = pd.qcut(x, q=n_bins, retbins=True, duplicates="drop")
    except ValueError:
        # Fall back to equally spaced bins
        binned, edges = pd.cut(x, bins=n_bins, retbins=True, include_lowest=True)

    # Re-label
    n_real_bins = len(edges) - 1
    if labels is None:
        labels = ["low", "medium", "high"][:n_real_bins] if n_bins == 3 else [f"bin{i+1}" for i in range(n_real_bins)]
    else:
        labels = list(labels)[:n_real_bins]

    binned = pd.Categorical(pd.cut(x, bins=edges, include_lowest=True, labels=labels, ordered=True))

    ranges = (
        pd.DataFrame({"score_strength": binned, "val": x})
        .groupby("score_strength", observed=True)["val"]
        .agg(min="min", max="max")
        .reset_index()
    )
    return binned, ranges.rename(columns={"val": ""})


# ----------------------------- Sorting -----------------------------

def _natural_key(v) -> tuple:
    """Numbers first (by integer if present), then lexicographic."""
    s = str(v)
    # find trailing integer
    i = len(s) - 1
    while i >= 0 and s[i].isdigit():
        i -= 1
    num = None
    if i < len(s) - 1:
        try:
            num = int(s[i+1:])
        except Exception:
            num = None
    return (0, num) if num is not None else (1, s)


# ----------------------------- Public API -----------------------------

def get_score_strength_per_chrom(
    input_windows: pd.DataFrame | pr.PyRanges,
    group_col: Optional[str] = None,
    score_col: Optional[str] = None,
    n_bins: int = 3,
    labels: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute per-group percentages of score-strength bins and bin ranges.

    Returns:
      pct:   columns = ['CHROM','score_strength','percent'] (CHROM is categorical & ordered)
      ranges: columns = ['score_strength','min','max']
    """
    df = _to_dataframe_any(input_windows)

    gcol = _detect_group_col(df, group_col)
    scol = _detect_score_col(df, score_col)

    # Build working frame
    work = pd.DataFrame({ "CHROM": df[gcol], "score": pd.to_numeric(df[scol], errors="coerce") })
    work = work.dropna(subset=["score"]).reset_index(drop=True)
    if work.empty:
        # empty outputs with expected schema
        pct = pd.DataFrame(columns=["CHROM", "score_strength", "percent"])
        ranges = pd.DataFrame(columns=["score_strength", "min", "max"])
        return pct, ranges

    # Bin
    binned, ranges = _bin_scores_quantile(work["score"], n_bins=n_bins, labels=labels)
    work["score_strength"] = binned

    # Percentages per group
    pct = (
        pd.crosstab(work["CHROM"], work["score_strength"], normalize="index")
          .mul(100)
          .stack(dropna=False)
          .rename("percent")
          .reset_index()
    )
    pct.columns = ["CHROM", "score_strength", "percent"]

    # Order groups naturally; keep only groups that appear
    chroms = sorted(pct["CHROM"].dropna().unique().tolist(), key=_natural_key)
    pct["CHROM"] = pd.Categorical(pct["CHROM"], categories=chroms, ordered=True)
    pct = pct.sort_values(["CHROM", "score_strength"])

    return pct, ranges
