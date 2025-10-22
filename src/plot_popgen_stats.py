#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import seaborn
from glob import glob
import pandas as pd
import argparse
import sys
from hotspotter.io import load_bed_files
from hotspotter.transform import preprocess_popgen_stats, concatenate_windows, sort_windows
import seaborn as sns

def _check_file_type(input_files: list[str]) -> None:
    allowed_suffixes: tuple[str, ...] = (".pi", ".windowpi", ".D", ".bed", ".tajimaD", ".snpden", ".snpdens")
    if not any(file.endswith(allowed_suffixes) for file in input_files):
        print(f"Error: Only files ending with {allowed_suffixes} are allowed.")
        sys.exit(1)
         
def plot_pop_gen_stats(
    data_frame: pd.DataFrame, 
    plotLineColor: str, 
    outFileName: str, 
    outFileFormat: str,
    y_axis_title: str,
    title: str
    ):
    # Extract sortable chromosome order and clean display label
    seaborn.set_style("dark")
    
    data_frame_colnames: pd.Index = data_frame.columns
    if "midpoint" not in data_frame_colnames.str.lower():
        midpoint: pd.Series = (data_frame[data_frame_colnames[1]] + data_frame[data_frame_colnames[2]]) // 2 
        data_frame.insert(loc=3, column="MIDPOINT", value=midpoint)
    
    # Sort chromosomes by the sortable value
    chromosomes = data_frame["CHROM"].unique()
    
    fig, axes = plt.subplots(
        11, 3,
        figsize=(10, 12),
        squeeze=False
    )
    
    axes = axes.flatten()
    print(data_frame.head())
    
    for i, chrom in enumerate(chromosomes):
        chromosome_data = data_frame[data_frame["CHROM"] == chrom]

        sns.lineplot(
            data=chromosome_data,
            x=chromosome_data.iloc[:, 3],
            y=chromosome_data.iloc[:, 5],
            color=plotLineColor,
            ax=axes[i],
            legend=False
        )
        axes[i].set_xlabel(f"Chromosome {chrom.replace("chr", "")}", fontsize=10, weight="bold")
        axes[i].set_ylabel(y_axis_title, fontsize=10, weight="bold")
    
    fig.suptitle(f"{title}", fontsize=16, weight='bold')
    plt.tight_layout()
    
    if outFileName and outFileFormat:
        fig.savefig(outFileName, format=outFileFormat, dpi=600, bbox_inches="tight")
    
    plt.show()
    
def parse_arguments():
    plot_stats = argparse.ArgumentParser(
        prog="Plot Pop-Gen Stats",
        usage='%(prog)s [options]', 
        description="Population Genomics Statistics Plotter"
    )
    
    plot_stats.add_argument('-i', '--input-files', 
                        required=True, 
                        nargs="+",
                        type=str, 
                        help='Input files of types snpdens, tajimaD, or window.pi')
    plot_stats.add_argument('-o', '--output-file', 
                        required=False, 
                        type=str,
                        default="plot.png",
                        help='Name of output file (plot), default: plot.png')
    plot_stats.add_argument('-f', '--output-file-format', 
                        required=False, 
                        type=str, 
                        choices=["png", "svg", "jpeg", "pdf"],
                        default="png",
                        help="Format of output figure, default: png")
    plot_stats.add_argument('-y', '--y-axis-title',
                        type=str,
                        required=False,
                        default="Value",
                        help='Title of the Y-axis')
    plot_stats.add_argument('-t', '--figure-title',
                        type=str,
                        required=False,
                        default="Value",
                        help='Title of the Whole figure')
    plot_stats.add_argument('-c', '--plot-color',
                        required=True,
                        type=str,
                        help='Color of the plot. default: black',
                        default="black")
    if len(sys.argv) == 1:
        plot_stats.print_help()
        sys.exit(1)
     
    return plot_stats.parse_args()


if __name__ == "__main__":
    
    try:
        args = parse_arguments()
    except Exception as e:
        print(f"Could not parse args, error: {e}")
        sys.exit(1)

    try:
        _check_file_type(args.input_files)
        loaded_files = load_bed_files(args.input_files)

    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
        sys.exit(1)

    except IsADirectoryError as e:
        print(f"Expected a file but found a directory: {e.filename}")
        sys.exit(1)

    except PermissionError as e:
        print(f"Permission denied: {e.filename}")
        sys.exit(1)

    except ValueError as e:
        print(f"Invalid input: {e}")
        sys.exit(1)

    except TypeError as e:
        print(f"Type error in input: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

    if not loaded_files:
        print("No files were loaded. Exiting.")
        sys.exit(1)

    try:
        preprocessed_data: list[pd.DataFrame] = preprocess_popgen_stats(loaded_files)
        sorted_data: list[pd.DataFrame] = sort_windows(preprocessed_data)
        concatenated_data: pd.DataFrame = concatenate_windows(sorted_data)
    except Exception as e:
        print(f"Could not merge data frames from files: {e}")
        sys.exit(1)

    try:
        plot_pop_gen_stats(
            data_frame=concatenated_data,
            y_axis_title=args.y_axis_title,
            plotLineColor=args.plot_color,
            outFileName=args.output_file,
            outFileFormat=args.output_file_format,
            title=args.figure_title
            )
    except Exception as e:
        print(f"Could not generate plot! Error: {e}")
        sys.exit(1)
