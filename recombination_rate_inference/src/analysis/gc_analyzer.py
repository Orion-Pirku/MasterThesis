from pybedtools import BedTool
import pandas as pd
import numpy as np
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr
from typing import Tuple, Union

def compute_bed_object_gc_content(
        genome_file: Union[str, Path],
        bed_object: Union[BedTool, pd.DataFrame]) -> BedTool:
  
    genome_fna_file = Path(genome_file).expanduser()
    
    if isinstance(bed_object, pd.DataFrame):
        bed = BedTool.from_dataframe(bed_object)
    else:
        bed = bed_object
    return bed.nucleotide_content(fi=genome_fna_file) #type: ignore


def compute_genomic_feature_correlation(
  bed_feature: BedTool,
  score_one_idx: int, 
  score_two_idx: int) -> Tuple[float, float]:
    
    # bed_feature must be a intersection between two individual 
    # bed objects using wa and wb as arguments.
    # score_one/two_idx are the indeces of the columns in the 
    # intersected bed_feature that contain the scores of interest   
    
    
    df = bed_feature.to_dataframe(names=None)
    
    if score_one_idx >= df.shape[1] or score_two_idx >= df.shape[1]:
        raise IndexError("Score index out of bounds")
    
    feature_one_array = df.iloc[:, score_one_idx].astype(float).to_numpy()
    feature_two_array = df.iloc[:, score_two_idx].astype(float).to_numpy()
   
    return pearsonr(
        x=feature_one_array, 
        y=feature_two_array, 
        alternative="two-sided"
        )


