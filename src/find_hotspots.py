#!/usr/bin/env python3
from pathlib import Path
from glob import glob
import argparse
from hotspotter.io import load_recombination_maps, save_hotspots_as_bed
from hotspotter.transform import interpolate_recombination_rate, chr_sort_key
from hotspotter.recombination_rate_analyzer import call_hotspots
from hotspotter.plotting import plot_recombination_hotspots
import os, sys

def parse_arguments() -> argparse.Namespace:
    
    parser = argparse.ArgumentParser(
        prog='find_hotspots', 
        usage='python script for calling meiotic recombination hotspots'
        )
    parser.add_argument(
        '--input-files', 
        '-i', nargs='+',
        type=str, 
        required=True, 
        help='directory containing the raw recombination map files'
        )
    parser.add_argument(
        '--output-directory', 
        '-o', 
        type=str, 
        required=True, 
        help='directory where the output will be stored'
        )
    
    return parser.parse_args()

def main() -> None:
    arguments = parse_arguments()
    os.makedirs(f"{arguments.output_directory}/plots", exist_ok=True)
    os.makedirs(f"{arguments.output_directory}/bed", exist_ok=True)
    
    try:
        input_files = arguments.input_files 
        recombination_maps = load_recombination_maps(input_files)
    except Exception as e:
        print(f"Failed to load recombination maps: {e}", file=sys.stderr)
        sys.exit(1)    
    
    sorted_recombination_maps = sorted(recombination_maps, key=chr_sort_key)
    interp_rec_rate, genomic_bins = interpolate_recombination_rate(sorted_recombination_maps)
    smoothed_rec_rate, peaks, props = zip(*[call_hotspots(rec_rate) for rec_rate in interp_rec_rate])
    
    for i, (signal, peak_idxs, bins) in enumerate(zip(smoothed_rec_rate, peaks, genomic_bins)):
        output_path = f"{arguments.output_directory}/plots/recombination_hotspots_chr{i+1}.png"
        plot_recombination_hotspots(
            smoothed_signal=signal,
            midpoint_bins=bins,
            peak_indeces=peak_idxs,
            chromosome_number=i+1,
            output_name=output_path)
        
        output_file = f"{arguments.output_directory}/bed/recombination_hotspots_chr{i+1}.bed"
        save_hotspots_as_bed(f'chr{i+1}', bins, signal, output_file, peak_idxs)
    
if __name__ == "__main__":
    main()