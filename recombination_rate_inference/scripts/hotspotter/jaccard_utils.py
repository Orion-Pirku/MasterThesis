from pybedtools import BedTool
import pandas as pd
from multiprocessing import Pool, cpu_count
import numpy as np
import pyranges as pr
from typing import Union, Mapping, Tuple, List
import numpy.typing as npt

def _load_chromsizes(chromsizes_file: str) -> dict:
    df = pd.read_csv(chromsizes_file, sep="\t", header=None, names=["chrom", "size"])
    return {row["chrom"]: (0, row["size"]) for _, row in df.iterrows()}


def _to_bedtool(obj: Union[pr.PyRanges, BedTool, pd.DataFrame]) -> BedTool:
    if isinstance(obj, BedTool):
        return obj.sort()
    elif isinstance(obj, pr.PyRanges):
        return BedTool.from_dataframe(obj.df).sort()
    elif isinstance(obj, pd.DataFrame):
        return BedTool.from_dataframe(obj).sort()
    else:
        raise TypeError("Input must be PyRanges, BedTool, or pandas DataFrame")


def _shuffle_jaccard_once(args) -> float:
    bed1_df, bed2_df, chromsizes = args
    bed1 = BedTool.from_dataframe(bed1_df)
    bed2 = BedTool.from_dataframe(bed2_df)
    shuffled_bed1 = bed1.shuffle(genome=chromsizes, chrom=True).sort()
    return shuffled_bed1.jaccard(bed2)["jaccard"]


def _shuffle_intervals(bed1_df: pd.DataFrame, bed2_df: pd.DataFrame, chromsizes: dict, iterations: int) -> npt.NDArray[np.float64]:
    args = [(bed1_df, bed2_df, chromsizes)] * iterations
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(_shuffle_jaccard_once, args)
    return np.array(results, dtype=np.float64)


def _empirical_pvalue_two_tailed(observed: float, null: Union[List[float], npt.NDArray[np.float64]]) -> float:
    null = np.asarray(null, dtype=np.float64)
    diff_obs = abs(observed - null.mean())
    return (np.sum(np.abs(null - null.mean()) >= diff_obs) + 1) / (len(null) + 1)

