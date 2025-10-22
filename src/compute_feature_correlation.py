#!/usr/bin/env python3
import argparse
import sys, os
from pathlib import Path
import re
import pandas as pd
import pyranges as pr

# your package imports
from hotspotter.io import (
        load_bed_files, 
        load_recombination_maps,
        _classify_feature,
        _normalize_paths)


from hotspotter.plotting import (
    plot_correlation_matrix
)
from hotspotter.transform import (
    make_windows, 
    sort_windows, 
    concatenate_windows, 
    preprocess_popgen_stats, 
    compute_gene_density
    )
from hotspotter.recombination_rate_analyzer import (
    compute_feature_correlation,
    compute_intersections,
    fill_correlation_matrices
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="compute_feature_correlation.py",
        description="Compute & plot Pearson correlation between recombination rate and any number of population-genomics features."
    )
    parser.add_argument(
            '--rho',
            type=str,
            nargs='+',
            required=True,
            help='Recombination map BED(s). Globs allowed.'
            )
    parser.add_argument(
            '--rho-name',
            type=str,
            required=True,
            help='Display name for the recombination rate layer (e.g., "rho").')
    # Repeated --feature flags: first token = label, rest = BED files
    parser.add_argument(
            '--feature',
            action='append',
            nargs='+',
            metavar=('LABEL', 'BED...'),
            help='Repeatable. Example: --feature tajima-d tajima_chr1.bed tajima_chr2.bed'
            )
    parser.add_argument(
            '--window-size', 
            '-ws', 
            type=int, 
            required=True,
            help='Window size in bp for recombination map windowing.'
            )
    parser.add_argument(
            '--output-directory', 
            '-o', 
            type=str, 
            default='.',
            help='Directory to write outputs (default: current dir).'
            )
    parser.add_argument(
            '--genome-sizes', 
            type=str, 
            help='File containing chromosome sizes'
            )
    args = parser.parse_args()

    # Require at least one feature
    if not args.feature or len(args.feature) == 0:
        parser.error("Provide at least one --feature <label> <beds...>")

    return args


def load_and_prepare_feature(
    label: str,
    paths: list[str] | str,
    genome_sizes_file: str,
    window_size: int,
) -> pr.PyRanges:

    feature_kind = _classify_feature(label)
    if feature_kind is None:
        raise ValueError(f"Unrecognized feature label: {label!r}")

    file_paths = _normalize_paths(paths)
    if not file_paths:
        raise ValueError("No input paths provided.")

    # Load all files into a list
    loaded_dfs: list[pd.DataFrame]
    if len(file_paths) > 1:
        # prefer your project’s loader if available
        try:
            loaded_dfs = load_bed_files(file_paths)  # type: ignore[name-defined]
        except NameError:
            loaded_dfs = [pd.read_csv(p, sep="\t", header=0) for p in file_paths]
    else:
        loaded_dfs = [pd.read_csv(file_paths[0], sep="\t", header=0)]

    multiple = len(loaded_dfs) > 1

    # ---- Dispatch per feature, concatenating ONLY if multiple files ----
    if feature_kind == "gene_density":
        df_in = concatenate_windows(loaded_dfs) if multiple else loaded_dfs[0]  # type: ignore[name-defined]
        processed = compute_gene_density(  # type: ignore[name-defined]
            df_in, genome_sizes=genome_sizes_file, window_size=window_size
        )

    elif feature_kind == "gc_content":
        processed = concatenate_windows(loaded_dfs) if multiple else loaded_dfs[0]  # type: ignore[name-defined]

    elif feature_kind == "popgen":
        preprocessed_list = preprocess_popgen_stats(  # type: ignore[name-defined]
            feature_dataframes=loaded_dfs
        )
        processed = (
            concatenate_windows(preprocessed_list)  # type: ignore[name-defined]
            if len(preprocessed_list) > 1
            else preprocessed_list[0]
        )

    else:
        # should be unreachable
        raise RuntimeError(f"Unhandled feature kind: {feature_kind!r}")

    # Normalize expected column names (no-op if already correct)
    processed = processed.rename(
        columns={"CHROM": "Chromosome", "START": "Start", "END": "End"},
        errors="ignore",
    )

    # Validate required columns exist before building PyRanges
    required = {"Chromosome", "Start", "End"}
    missing = required.difference(processed.columns)
    if missing:
        raise ValueError(f"Processed feature missing columns: {sorted(missing)}")

    return pr.PyRanges(df=processed)


def main() -> None:
    args = parse_arguments()
    outdir = args.output_directory
    ws = args.window_size

    try:
        # Recombination maps (windowed)
        rec_maps = load_recombination_maps(args.rho)
        sorted_rec_map_dataframes = sort_windows(rec_maps)
        rec_maps_windowed = [make_windows(df, ws) for df in sorted_rec_map_dataframes]
        concatenated_rec_maps = concatenate_windows(rec_maps_windowed)
        concatenated_rec_maps = concatenated_rec_maps.rename(
                columns={
                    "CHROM": "Chromosome",
                    "START": "Start",
                    "END": "End"
                    }
                )
        rec_maps_bed = pr.PyRanges(concatenated_rec_maps)

        # Generic features
        beds: dict[str, pr.PyRanges] = {args.rho_name: rec_maps_bed}
        for feat_spec in args.feature:
            label, *files = feat_spec
            if len(files) == 0:
                raise ValueError(f"--feature {label} has no files")
            beds[label] = load_and_prepare_feature(
                label, 
                files, 
                window_size=args.window_size, 
                genome_sizes_file=args.genome_sizes
            )

    except FileNotFoundError as e:
        print(f"[ERROR] file not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] {type(e).__name__}: {e}")
        sys.exit(1)

    # Pairwise intersects and correlations
    pairwise_intersections = compute_intersections(beds)
    results = {
        name: compute_feature_correlation(
            bt, 
            score_A_idx=3, 
            score_B_idx=8
        )
        for name, bt in pairwise_intersections.items()
    }

    labels = list(beds.keys())
    R_df, P_df = fill_correlation_matrices(results, labels)
    plot_correlation_matrix(R_df, P_df, labels, outdir)


if __name__ == "__main__":
    main()

