import pandas as pd
from typing import List
from pybedtools import BedTool
from typing import Union, List, Dict, Mapping
import numpy as np
import pyranges as pr
import numpy.typing as npt
from multiprocessing import Pool, cpu_count
from pathlib import Path
 

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
