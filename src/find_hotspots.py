#!/usr/bin/env python3
from pathlib import Path
import argparse
from tqdm import tqdm
from hotspotter.io import load_recombination_maps, load_bed_files, save_hotspots_as_bed
from hotspotter.transform import chr_sort_key
from hotspotter.recombination_rate_analyzer import call_hotspots_cwt
from hotspotter.plotting import plot_recombination_hotspots, plot_scalogram
import sys
import numpy as np
import numpy.typing as npt
import pandas as pd


def _get_rho(dataframe: pd.DataFrame) -> npt.NDArray[np.float64]:
    rho_column = dataframe.iloc[:, 3].to_numpy(dtype=float)
    return rho_column


def _get_midpoint(dataframe: pd.DataFrame) -> npt.NDArray[np.int64]:
    dataframe["MIDPOINT"] = (
        dataframe["START"].astype(int) + dataframe["END"].astype(int)
    ) // 2
    midpoint_array = dataframe["MIDPOINT"].to_numpy(dtype=int)
    return midpoint_array


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="find_hotspots",
        usage="%(prog)s [options]",
        description="CLI Tool to call meiotic recombination hotspots using wavelet transformation",
    )
    parser.add_argument(
        "--input-files",
        "-i",
        nargs="+",
        type=str,
        required=True,
        help="directory containing the raw recombination map files",
    )
    parser.add_argument(
        "-e",
        "--effective-population-size",
        type=float,
        default=0.0,
        required=False,
        help="effective population size for computation of cM/Mb. Default is 0.0",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=True,
        help="Threshold value over which to compute the peak prominence used to call the hotspot peaks",
    )
    parser.add_argument(
        "--output-directory",
        "-o",
        type=str,
        required=True,
        help="directory where the output will be stored",
    )
    if len(sys.argv) == 1:
        parser.print_help()

    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()

    plots_dir = Path(arguments.output_directory) / "plots"
    bed_dir = Path(arguments.output_directory) / "bed"
    plots_dir.mkdir(parents=True, exist_ok=True)
    bed_dir.mkdir(parents=True, exist_ok=True)

    try:
        input_files = arguments.input_files
        if not input_files:
            print("No input files provided.", file=sys.stderr)
            sys.exit(1)

        p0 = Path(input_files[0])
        suffixes = [s.lower() for s in p0.suffixes]
        is_bed = ".bed" in suffixes

        if is_bed:
            recombination_maps = load_bed_files(input_files)
        else:
            recombination_maps = load_recombination_maps(input_files)

    except Exception as e:
        print(f"Failed to load recombination maps: {e}", file=sys.stderr)
        sys.exit(1)

    sorted_recombination_maps = sorted(recombination_maps, key=chr_sort_key)

    outer = tqdm(
        total=len(sorted_recombination_maps),
        desc="Processing maps",
        unit="chr",
        dynamic_ncols=True,
    )

    for idx, rec_map in enumerate(sorted_recombination_maps, start=1):
        inner = tqdm(total=6, desc=f"chr{idx} steps", leave=False, dynamic_ncols=True)
        try:
            rho_signal = _get_rho(rec_map)
            inner.update(1)
            midpoint_bins = _get_midpoint(rec_map)
            inner.update(1)
            peak_idxs, smoothed_rho, coeff_matrix, scales = call_hotspots_cwt(
                rho_signal, prominence_cutoff=arguments.threshold
            )
            inner.update(1)

            plot_path = (
                plots_dir
                / f"recombination_hotspots_chr{idx}_t{arguments.threshold}.png"
            )
            plot_recombination_hotspots(
                smoothed_signal=smoothed_rho,
                midpoint_bins=midpoint_bins,
                peak_indeces=peak_idxs,
                chromosome_number=idx,
                output_name=str(plot_path),
            )
            inner.update(1)
            scalogram_path = (
                plots_dir / f"scalogram_chr{idx}_t{arguments.threshold}.png"
            )
            plot_scalogram(
                coeff_matrix,
                scales,
                xlabel=f"Chromosome{idx}",
                output_name=scalogram_path,
            )
            inner.update(1)

            bed_path = (
                bed_dir / f"recombination_hotspots_chr{idx}_t{arguments.threshold}.bed"
            )
            save_hotspots_as_bed(
                chrom_name=f"chr{idx}",
                genomic_midpoint=midpoint_bins,
                smoothed_signal=rho_signal,
                output_file_name=str(bed_path),
                hotspot_index=peak_idxs,  # type: ignore
            )
            inner.update(1)

        except Exception as e:
            tqdm.write(f"[ERROR] chr{idx}: {e}", file=sys.stderr)
        finally:
            inner.close()
            outer.update(1)

    outer.close()


if __name__ == "__main__":
    main()
