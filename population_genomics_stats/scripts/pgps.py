import os
import polars as pl
import matplotlib.pyplot as plt
import seaborn
from glob import glob
import pandas as pd
from typing import Dict, List, Union, Tuple, Optional, Pattern
import pyarrow
import argparse
from os import sep
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

def load_files(directory_path: str, file_pattern: str) -> Dict[str, pl.DataFrame]:
    """
    Optimized: Loads files matching a specific pattern from a directory using Polars
    and stores them with their chromosome label as key.
    """
    loaded_files_dict: Dict[str, pl.DataFrame] = {}
    files: List[str] = glob(os.path.join(directory_path, file_pattern))

    if not files:
        print(f"No files matching the pattern '{file_pattern}' were found in '{directory_path}'.")
        return loaded_files_dict

    # Precompile regex
    chr_pattern: Pattern[str] = re.compile(r"chr_(\d+|[A-Z])")

    def read_and_parse(file: str) -> Optional[Tuple[str, pl.DataFrame]]:
        if not os.path.isfile(file):
            return None

        basename = os.path.basename(file)
        file_basename = os.path.splitext(basename)[0]
        match = chr_pattern.search(file_basename)

        if not match:
            print(f"Skipping file without valid chr label: {basename}")
            return None

        chr_label = match.group()
        try:
            df = pl.read_csv(file, separator='\t')
            return chr_label, df
        except Exception as e:
            print(f"Error reading {file}: {e}")
            return None

    # Use threads for parallel I/O (safe for Polars)
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as executor:
        futures = {executor.submit(read_and_parse, f): f for f in files}
        for future in as_completed(futures):
            result = future.result()
            if result:
                chr_label, df = result
                loaded_files_dict[chr_label] = df

    return loaded_files_dict

def inspect_column_types(dataFrames: dict[str, pl.DataFrame], column: str = "CHROM"):
    print(f"Inspecting column types for '{column}':")
    for key, df in dataFrames.items():
        if column in df.columns:
            dtype = df.schema[column]
            print(f" - {key}: {column} -> {dtype}")
        else:
            print(f" - {key}: {column} column NOT FOUND")


def merge_dataframes_from_dict(dataFrames: dict[str, pl.DataFrame], first_file_keyName: str) -> pl.DataFrame:
    """
    Vertically stacks DataFrames stored in a dictionary.
    Assumes all DataFrames have the same schema.
    """

    if first_file_keyName not in dataFrames:
        raise ValueError(f"Key '{first_file_keyName}' not found in dataFrames.")

    # Optional: Cast all CHROM columns to string to prevent type mismatch
    #for key in dataFrames:
     #   if "CHROM" in dataFrames[key].columns:
      #      dataFrames[key] = dataFrames[key].with_columns([pl.col("CHROM").cast(pl.Utf8), pl.col("TajimaD").cast(pl.Float64)])
       # else:
        #    raise KeyError(f"'CHROM' column missing in DataFrame with key '{key}'.")

    merged_df = dataFrames[first_file_keyName]

    for key, df in dataFrames.items():
        if key != first_file_keyName:
            merged_df = merged_df.vstack(df)

    return merged_df

def sort_by_chromosome(dataFrame: pl.DataFrame) -> pl.DataFrame:
    
    dataFrame = dataFrame.with_columns([
        pl.col("CHROM")
        .str.extract(r"(\d+)$")
        .cast(pl.Int32)
        .alias("CHROM_NUM")
    ])
    
    dataFrame = dataFrame.sort("CHROM_NUM")
    
    return dataFrame.drop("CHROM_NUM")


def chrom_to_sortable(chrom: str) -> Union[int, float]:
    chrom = chrom.lower().replace("chr_", "").replace("chr", "")
    if chrom.isdigit():
        return int(chrom)
    elif chrom == "w":
        return 34
    elif chrom == "z":
        return 35
    else:
        return float("inf")
    

def chrom_numeric_label(chrom: str) -> Union[str, int]:
    """Map 'chr_X' -> 'X', 'chr_1' -> 1, etc. for clean display."""
    chrom = chrom.lower().replace("chr_", "").replace("chr", "")
    return chrom.upper() if not chrom.isdigit() else int(chrom)
         
def plot_popGen(dataFrame: pl.DataFrame, 
                            xValues: str, 
                            yValues: str,
                            yAxisTitle: str, 
                            plotLineColor: str, 
                            outFileName: str, 
                            outFileFormat: str):
    
    pandas_dataFrame: pd.DataFrame = dataFrame.to_pandas()

    # Extract sortable chromosome order and clean display label
    pandas_dataFrame["CHROM_SORT"] = pandas_dataFrame["CHROM"].apply(chrom_to_sortable)
    pandas_dataFrame["CHROM_LABEL"] = pandas_dataFrame["CHROM"].apply(chrom_numeric_label)

    seaborn.set_style("darkgrid")
    
    # Sort chromosomes by the sortable value
    chromosomes = pandas_dataFrame.drop_duplicates(subset="CHROM").sort_values("CHROM_SORT")

    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(14, min(6 * len(chromosomes), 60)))

    if len(chromosomes) == 1:
        axes = [axes]

    for i, (_, row) in enumerate(chromosomes.iterrows()):
        chrom_label = row["CHROM_LABEL"]
        chrom_sort_val = row["CHROM_SORT"]
        chrom_orig = row["CHROM"]
        chromosome_data = pandas_dataFrame[pandas_dataFrame["CHROM"] == chrom_orig]

        seaborn.lineplot(
            data=chromosome_data,
            x=xValues,
            y=yValues,
            color=plotLineColor,
            ax=axes[i] # type: ignore
        )
        axes[i].set_xlabel(f"Chromosome {chrom_label}", fontsize=12)
        axes[i].set_ylabel(yAxisTitle, fontsize=12)

    plt.tight_layout()
    
    return plt.savefig(outFileName, format=outFileFormat, dpi=600, bbox_inches="tight")
    
def CommandLineArguments():
    plot_stats = argparse.ArgumentParser(prog="pgsp",
                                     usage='%(prog)s [options]', 
                                     description="Population Genomics Statistics Plotter")
    plot_stats.add_argument('-i', '--input_directory', 
                        required=True, 
                        type=str, 
                        help='input directory of files of type snpdens, tajimaD or window.pi' )
    plot_stats.add_argument('-p', '--file_pattern', 
                        type=str, 
                        required=True,
                        help='input the pattern of the files whose stats you want to plot e.g., *_SNP_Density_10kb*')
    plot_stats.add_argument('-o', '--output_file', 
                        required=True, 
                        type=str,
                        help='name of output file (plot)')
    plot_stats.add_argument('-f', '--output_file_format', 
                        required=True, 
                        type=str, 
                        help="format of output figure e.g., svg, png, jpeg.")
    plot_stats.add_argument('-X', '--X_axis_values', 
                        type=str,
                        required=True,
                        help='Name of Column from input files to be plotted on the X-axis')
    plot_stats.add_argument('-Y', '--Y_axis_values',
                        type=str,
                        required=True,
                        help='Name of Column from input files to be plotted on the Y-axis')
    plot_stats.add_argument('-y', '--y_axis_title',
                        type=str,
                        required=True,
                        help='Name of y-axis')
    plot_stats.add_argument('-c', '--plot_line_color',
                        required=True,
                        type=str,
                        help='color of the lineplot')
     
    cl_arguments = plot_stats.parse_args()
    return cl_arguments

if __name__ == "__main__":
    
    try:
        parsedArguments = CommandLineArguments()
    except Exception as e:
        print(f"Could not parse args, error: {e}")
        sys.exit(1)

    try:
        loadedFile = load_files(
            parsedArguments.input_directory,
            parsedArguments.file_pattern
        )
    except Exception as e:
        print(f"Could not load files: {e}")
        sys.exit(1)

    if not loadedFile:
        print("No files were loaded. Exiting.")
        sys.exit(1)

    try:
        first_file_key = next(iter(loadedFile))
        mergedDataFrames = merge_dataframes_from_dict(loadedFile, first_file_key)
        sorted_merged_dataFrame = sort_by_chromosome(mergedDataFrames)
    except Exception as e:
        print(f"Could not merge data frames from files: {e}")
        sys.exit(1)

    try:
        plot_popGen(
            dataFrame=sorted_merged_dataFrame,
            xValues=parsedArguments.X_axis_values,
            yValues=parsedArguments.Y_axis_values,
            yAxisTitle=parsedArguments.y_axis_title,
            plotLineColor=parsedArguments.plot_line_color,
            outFileName=parsedArguments.output_file,
            outFileFormat=parsedArguments.output_file_format
        )
    except Exception as e:
        print(f"Could not generate plot! Error: {e}")
        sys.exit(1)
