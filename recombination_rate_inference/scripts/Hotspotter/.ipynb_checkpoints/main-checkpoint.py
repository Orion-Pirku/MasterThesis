from transformation import transform_bed, transform_rmap
from plotting import *
import io
from analysis import *
from pybedtools import BedTool
from typing import cast
from pathlib import Path

if __name__=="__main__":
    
    rec_rate_bed = io_bed.load_bed_files(
        "/home/orion/MasterThesis-main/recombination_rate_inference/output/processed/10kb_windows/*", 
        return_type="bed"
        )
    
    blackcap_features = pd.read_csv(
        "/home/orion/MasterThesis-main/recombination_rate_inference/input/GCF_009819655.1_bSylAtr1.pri_feature_table.txt", 
        sep = '\t', 
        header = 0, 
        dtype=str)
    
     
    blackcap_gff: BedTool = transform_bed.transform_feature_table(
            blackcap_features,
            return_type = "bed",
        )
     
    blackcap_gff_genes = BedTool.filter(
        blackcap_gff,
        lambda bed: bed.name in ("gene","mRNA","CDS")
    ).sort()
    
    blackcap_genome = (
        pd.read_csv("/home/orion/MasterThesis-main/recombination_rate_inference/input/SylAtri_genome.txt.tsv", sep="\t", usecols=['Chromosome name', 'Seq length'])
        .assign(**{'Chromosome name': lambda df: 'chr' + df['Chromosome name'].astype(str)})
        .head(33)
    )
    gc = compute_gc_content(rec_rate_bed, "~/MasterThesis-main/recombination_rate_inference/input/GCA_009819655.1_bSylAtr1.pri_genomic.fasta")
    print(gc.head(10))