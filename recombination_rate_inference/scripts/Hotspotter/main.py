from hotspotter_transformation import *
from hotspotter_plotting import *
import hotspotter_io
from hotspotter_analysis import *
from pybedtools import BedTool
from typing import cast
from pathlib import Path

# def compute_jaccard_index()

if __name__=="__main__":
    
    rec_rate_bed = hotspotter_io.load_bed_files(
        "/home/orion/MasterThesis-main/recombination_rate_inference/output/processed/10kb_windows/*", 
        return_type="bed"
        )
    
    blackcap_features = pd.read_csv(
        "/home/orion/MasterThesis-main/recombination_rate_inference/input/GCF_009819655.1_bSylAtr1.pri_feature_table.txt", 
        sep = '\t', 
        header = 0, 
        dtype=str)
    
     
    blackcap_gff: BedTool = transform_feature_table(
            blackcap_features,
            return_type = "bed",
        )
     
    blackcap_gff_genes = BedTool.filter(
        blackcap_gff,
        lambda bed: bed.name in ("gene","mRNA","CDS")
    ).sort()
    
    # blackcap_genome = pd.read_csv(
        # "../../input/original_genome_sizes.txt", 
        # sep="\t", 
        # header=None, 
        # names=['chromosome', 'size']
        # )
    
    #genome_mapping = create_chromosome_mapping(blackcap_genome)
    #hotspotter_io.parse_fasta_file("../../input/GCA_009819655.1_bSylAtr1.pri_genomic.fasta", genome_mapping) 
    gc = compute_gc_content(rec_rate_bed, genome_fna_path="./blackcap.fasta")
    print(gc.sort().head(1))
