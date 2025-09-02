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

