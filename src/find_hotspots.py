from pathlib import Path
from glob import glob
import argparse
from hotspotter.io import load_recombination_maps
from hotspotter.transform import interpolate_recombination_rate, chr_sort_key
from hotspotter.recombination_rate_analyzer import call_hotspots
from hotspotter.plotting import plot_recombination_hotspots

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog='find_hotspots', 
                                     usage='python script for calling meiotic recombination hotspots')
    parser.add_argument('--input-file-directory', '-i', type=str, required=True, description='directory containing the raw recombinaiton map files')
    parser.add_argument()
    return parser.parse_args()