from pybedtools import BedTool
import pandas as pd
from multiprocessing import Pool, cpu_count
import numpy as np
import pyranges as pr
from typing import (
    Union,
    Mapping,
    Tuple,
    List
    )
import numpy.typing as npt
from jaccard_utils import (
    _load_chromsizes, 
    _to_bedtool, 
    _shuffle_jaccard_once, 
    _shuffle_intervals, 
    _empirical_pvalue_two_tailed
    )

def compute_jaccard_index(
    bed_input1: Union[pr.PyRanges, BedTool, pd.DataFrame],
    bed_input2: Union[pr.PyRanges, BedTool, pd.DataFrame, Mapping[str, Union[pr.PyRanges, BedTool, pd.DataFrame]]],
    chromsizes_file: str,
    iterations: int = 1000) -> pd.DataFrame:

    chromsizes = _load_chromsizes(chromsizes_file)
    bed1 = _to_bedtool(bed_input1)
    bed1_df = bed1.to_dataframe()

    if isinstance(bed_input2, Mapping):
        bed2_dict = {k: _to_bedtool(v) for k, v in bed_input2.items()}
    else:
        bed2 = _to_bedtool(bed_input2)
        name = (
            bed_input2.df.iloc[0, 3]
            if isinstance(bed_input2, pr.PyRanges) and bed_input2.df.shape[1] > 3
            else bed_input2.iloc[0, 3]
            if isinstance(bed_input2, pd.DataFrame) and bed_input2.shape[1] > 3
            else "feature"
        )
        bed2_dict = {name: bed2}

    records = []
    for feature, bed2 in bed2_dict.items():
        bed2_df = bed2.to_dataframe()
        observed = bed1.jaccard(bed2)["jaccard"]
        null_dist = _shuffle_intervals(bed1_df, bed2_df, chromsizes, iterations)

        records.append({
            "feature": feature,
            "observed_jaccard": observed,
            "shuffled_mean": null_dist.mean(),
            "shuffled_median": np.median(null_dist),
            "ci_lower_95%": np.percentile(null_dist, 2.5),
            "ci_upper_95%": np.percentile(null_dist, 97.5),
            "p_value_two_tailed": _empirical_pvalue_two_tailed(observed, null_dist)
        })

    return pd.DataFrame(records)
