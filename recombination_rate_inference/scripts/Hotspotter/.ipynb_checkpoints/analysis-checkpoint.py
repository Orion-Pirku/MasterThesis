import pandas as pd
from typing import List
from pybedtools import BedTool
from typing import Union, List, Dict, Mapping
import numpy as np
import pyranges as pr
import numpy.typing as npt
from multiprocessing import Pool, cpu_count
from scipy.stats import pearsonr
from pathlib import Path

class GenomeStats:
      __slots__ = ("chromsizes")
    
      def __init__(self, chromsizes: str):
        self.chromsizes = self._load_chromsizes(chromsizes)
      @staticmethod
      def _load_chromsizes(chromsizes: str) -> dict:
          df = pd.read_csv(chromsizes, sep="\t", header=None, names=["chrom", "size"])
          return {row["chrom"]: (0, row["size"]) for _, row in df.iterrows()}
    
      @staticmethod
      def to_bedtool(obj: Union[pr.PyRanges, BedTool, pd.DataFrame]) -> BedTool:
          if isinstance(obj, BedTool):
              return obj.sort()
          elif isinstance(obj, pr.PyRanges):
              return BedTool.from_dataframe(obj.df).sort()
          elif isinstance(obj, pd.DataFrame):
              return BedTool.from_dataframe(obj).sort()
          else:
              raise TypeError("Input must be PyRanges, BedTool, or pandas DataFrame")
    
      @staticmethod
      def _single_iter(args):
          bed1_df, bed2_df, chromsizes = args
          bed1 = BedTool.from_dataframe(bed1_df)
          bed2 = BedTool.from_dataframe(bed2_df)
          shuffled_bed1 = bed1.shuffle(genome=chromsizes, chrom=True).sort()
          return shuffled_bed1.jaccard(bed2)["jaccard"]
    
      def shuffle_intervals(
          self, bed1_df: pd.DataFrame, 
          bed2_df: pd.DataFrame,
          iterations: int) -> npt.NDArray[np.float64]:
          
          args = [(bed1_df, bed2_df, self.chromsizes)] * iterations
          with Pool(processes=cpu_count()) as pool:
              results = pool.map(self._single_iter, args)
          return np.array(results, dtype=np.float64)
    
      @staticmethod
      def empirical_pvalue_two_tailed(
          observed: float, 
          null: List[float] | npt.NDArray[np.float64]) -> float:
          if not isinstance(null, np.ndarray):
              null = np.asarray(null, dtype=np.float64)
          null_mean = null.mean()
          diff_obs = abs(observed - null_mean)
          p_val = (np.sum(np.abs(null - null_mean) >= diff_obs) + 1) / (len(null) + 1)
          return p_val
    
      def compute_jaccard(
          self,
          bed_input1: Union[pr.PyRanges, BedTool, pd.DataFrame],
          bed_input2: Union[pr.PyRanges, BedTool, pd.DataFrame, Mapping[str, Union[pr.PyRanges, BedTool, pd.DataFrame]]],
          iterations: int = 1000 ) -> pd.DataFrame:
          bed1 = self.to_bedtool(bed_input1)
          bed1_df = bed1.to_dataframe()
    
          if isinstance(bed_input2, Mapping):
              bed2_dict = {k: self.to_bedtool(v) for k, v in bed_input2.items()}
          else:
              bed2 = self.to_bedtool(bed_input2)
              feature_name = (
                  bed_input2.df.iloc[0, 3]
                  if isinstance(bed_input2, pr.PyRanges) and bed_input2.df.shape[1] > 3
                  else bed_input2.iloc[0, 3]
                  if isinstance(bed_input2, pd.DataFrame) and bed_input2.shape[1] > 3
                  else "feature"
              )
              bed2_dict = {feature_name: bed2}
    
          records = []
          for feature, bed2 in bed2_dict.items():
              bed2_df: pd.DataFrame = bed2.to_dataframe()
              observed_jaccard: float = bed1.jaccard(bed2)["jaccard"]
              shuffled_jaccards: npt.NDArray[np.float64] = self.shuffle_intervals(bed1_df, bed2_df, iterations)
    
              ci_lower: np.float64 = np.percentile(shuffled_jaccards, 2.5)
              ci_upper: np.float64 = np.percentile(shuffled_jaccards, 97.5)
              mean: np.float64 = shuffled_jaccards.mean()
              median: np.float64 = np.median(shuffled_jaccards)
              pval: float = self.empirical_pvalue_two_tailed(observed_jaccard, shuffled_jaccards)
    
              records.append({
                  "feature": feature,
                  "observed_jaccard": observed_jaccard,
                  "shuffled_jaccard_mean": mean,
                  "shuffled_jaccard_median": median,
                  "ci_lower_95%": ci_lower,
                  "ci_upper_95%": ci_upper,
                  "p_value_two_tailed": pval
              })
    
          return pd.DataFrame(records)
    

def compute_gc_content(
    bed_object: BedTool | pd.DataFrame,
    genome_fna_path: str) -> BedTool:
    genome_fna_file = Path(genome_fna_path).expanduser()
    bed = (
        BedTool.from_dataframe(bed_object)
        if isinstance(bed_object, pd.DataFrame)
        else bed_object
    )
    return bed.nucleotide_content(fi=genome_fna_file)


def compute_genome_wide_correlation(
    bed_1: BedTool | np.ndarray,
    bed_2: BedTool ) -> np.ndarray:

    bed_1_np = bed1.to_dataframe().iloc[:, 2].to_numpy()
    bed_2_np = bed2.to_dataframe().iloc[:, 2].to_numpy()

def summarize_recombination_windows(
    windows: List[pd.DataFrame],
    to_latex: bool = False) -> pd.DataFrame:
    
    summary_table = pd.DataFrame([
        {
            "chrom": df["Chrom"].iloc[0],
            **pd.to_numeric(df["cM/Mb"], errors="coerce").describe()
        }
        for df in windows
    ])

    columns = ["chrom", "count", "mean", "std", "min", "25%", "50%", "75%", "max"]
    summary_table = summary_table[columns]
    summary_table[columns[2:]] = summary_table[columns[2:]].round(3)

    summary_table = summary_table.sort_values(
        "chrom", 
        key=lambda x: x.str.extract(r'chr(\d+|W|Z)')[0].map(lambda s: int(s) if s and s.isdigit() else float('inf'))
    )
    
    if to_latex:
        summary_table.to_latex("recombination_map_summary.tex", index=False, float_format="%.3f")
        
    return summary_table
