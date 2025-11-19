# io.py
from __future__ import annotations
import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Union
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy.typing as npt
import numpy as np
import pandas as pd
import pyranges as pr
from .transform import (
    concatenate_windows,
    compute_gene_density,
    preprocess_popgen_stats,
    make_windows,
    sort_windows,
)
import sys
import os
import re
from collections.abc import Iterable, Sequence
from typing import Literal


FeatureKind = Literal[
    "gene density", "gc content", "rec rate", "snp density", "tajima d", "window pi"
]


def _check_file_type(input_files: Sequence[Path]) -> None:
    allowed_suffixes: tuple[str, ...] = (
        ".pi",
        ".windowpi",
        ".D",
        ".bed",
        ".tajimaD",
        ".snpden",
        ".snpdens",
        ".rmap",
    )
    if not any(f.name.endswith(allowed_suffixes) for f in input_files):
        suffix_str = ", ".join(allowed_suffixes)
        raise ValueError(f"Only files ending with one of ({suffix_str}) are allowed.")


def _classify_feature(label: str) -> FeatureKind | None:
    label_lower = label.lower().replace("_", " ")
    if "gene density" in label_lower or "genes" in label_lower:
        return "gene density"
    if "gc content" in label_lower:
        return "gc content"
    if "rec rate" in label_lower:
        return "rec rate"
    for key in ("snp density", "tajima d", "window pi"):
        if key in label_lower:
            return key
    return None


def _normalize_paths(paths: str | Path | Iterable[str | Path]) -> list[Path]:
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(p) for p in paths]


def _ensure_non_empty_dataframe(
    df: pd.DataFrame, feature_kind: FeatureKind
) -> pd.DataFrame:
    if df.empty:
        raise ValueError(
            f"Processed feature '{feature_kind}' produced an empty DataFrame"
        )
    return df


def _rename_genomic_columns(df: pd.DataFrame) -> pd.DataFrame:
    renamed = df.rename(
        columns={"CHROM": "Chromosome", "START": "Start", "END": "End"},
        errors="ignore",
    )
    required = {"Chromosome", "Start", "End"}
    missing = required.difference(renamed.columns)
    if missing:
        raise ValueError(f"Processed feature missing columns: {sorted(missing)}")
    return renamed


def _process_gene_density(
    file_paths: Sequence[Path],
    genome_sizes_file: str,
    window_size: int,
) -> pd.DataFrame:
    loaded_dfs = load_bed_files([str(p) for p in file_paths])
    df_in = concatenate_windows(loaded_dfs) if len(loaded_dfs) > 1 else loaded_dfs[0]
    processed = compute_gene_density(
        df_in,
        genome_sizes=genome_sizes_file,
        window_size=window_size,
    )
    processed.to_csv(
        "blackcap_gene_density.bed",
        sep="\t",
        index=False,
        header=True,
    )
    return processed


def _process_rec_rate(
    file_paths: Sequence[Path],
    window_size: int,
) -> pd.DataFrame:
    rec_rates = load_recombination_maps([str(p) for p in file_paths])
    rec_rates_windows = [make_windows(df, window_size) for df in rec_rates]
    sorted_rec_rates = sort_windows(rec_rates_windows)
    return concatenate_windows(sorted_rec_rates)


def _process_gc_content(file_paths: Sequence[Path]) -> pd.DataFrame:
    loaded_dfs = load_bed_files([str(p) for p in file_paths])
    return concatenate_windows(loaded_dfs) if len(loaded_dfs) > 1 else loaded_dfs[0]


def _process_popgen_stat(
    file_paths: Sequence[Path],
    feature_kind: Literal["snp density", "tajima d", "window pi"],
) -> pd.DataFrame:
    loaded_dfs = load_bed_files([str(p) for p in file_paths])
    preprocessed_list = preprocess_popgen_stats(feature_dataframes=loaded_dfs)
    processed = (
        concatenate_windows(preprocessed_list)
        if len(preprocessed_list) > 1
        else preprocessed_list[0]
    )
    out = feature_kind.replace(" ", "_")
    processed.to_csv(
        f"blackcap_{out}.bed",
        sep="\t",
        index=False,
    )
    return processed


def load_and_prepare_feature(
    label: str,
    paths: str | Path | Iterable[str | Path],
    genome_sizes_file: str,
    window_size: int,
) -> pr.PyRanges:
    feature_kind = _classify_feature(label)
    if feature_kind is None:
        raise ValueError(f"Unrecognized feature label: {label!r}")

    file_paths = _normalize_paths(paths)
    if not file_paths:
        raise ValueError("No input paths provided.")

    _check_file_type(file_paths)

    match feature_kind:
        case "gene density":
            processed = _process_gene_density(
                file_paths=file_paths,
                genome_sizes_file=genome_sizes_file,
                window_size=window_size,
            )
        case "rec rate":
            processed = _process_rec_rate(
                file_paths=file_paths,
                window_size=window_size,
            )
        case "gc content":
            processed = _process_gc_content(file_paths)
        case "snp density" | "tajima d" | "window pi":
            processed = _process_popgen_stat(
                file_paths=file_paths,
                feature_kind=feature_kind,
            )

    processed = _ensure_non_empty_dataframe(processed, feature_kind)
    processed = _rename_genomic_columns(processed)

    return pr.PyRanges(df=processed)


def save_hotspots_as_bed(
    chrom_name: str,
    genomic_midpoint: npt.NDArray[np.int64],
    smoothed_signal: npt.NDArray[np.float64],
    output_file_name: str,
    hotspot_index: npt.NDArray[np.int64],
) -> None:
    start: npt.NDArray[np.int64] = genomic_midpoint - 50
    end: npt.NDArray[np.int64] = genomic_midpoint + 50
    bed_object: pd.DataFrame = pd.DataFrame(
        {
            "Chrom": [chrom_name] * len(hotspot_index),
            "Start": start[hotspot_index],
            "End": end[hotspot_index],
            "cM/Mb": smoothed_signal[hotspot_index],
        }
    ).astype({"Chrom": str, "Start": int, "End": int, "cM/Mb": float})
    return bed_object.to_csv(output_file_name, sep="\t", index=False, header=False)


def _parse_rmap_file(rmap_file: str) -> pd.DataFrame:
    rmap_file_path = Path(rmap_file)
    if match := re.search(r"(chr\d{1,2}A?)", str(rmap_file_path)):
        chromosome = match.group(1)
    else:
        raise ValueError(f"Could not find chromosome in filename:{rmap_file_path}")
    df = pd.read_csv(rmap_file_path, sep="\t", header=None, comment="#")
    df.insert(0, "CHROM", chromosome)
    df.rename(columns={0: "START", 1: "END", 2: "RHO"}, inplace=True)
    return df


def load_recombination_maps(rmap_files: str | list[str]) -> list[pd.DataFrame]:
    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dataframes = pool.map(_parse_rmap_file, rmap_files)
    return dataframes


def save_transformed_data(datasets: List[pd.DataFrame], prefix="transformed") -> None:
    for data in datasets:
        chrom_name = data["chrom"].iat[0]
        window_size = data["end"].iat[0] - data["start"].iat[0]
        filename = f"{prefix}_{chrom_name}_w{window_size}.tsv"
        data.to_csv(filename, sep="\t", header=True, index=False)


def _parse_bed_file(file: str) -> pd.DataFrame:
    """
    Parse a BED file into a DataFrame with standardized columns.
    Assumes the first line is a header to skip.
    """
    file_path = Path(file)
    data_frame = pd.read_csv(file_path, sep="\t", header=0)
    if data_frame.iloc[:, 0].str.contains("_").any():
        data_frame.iloc[:, 0] = data_frame.iloc[:, 0].str.replace("_", "", regex=False)
    column_names = [col.lower() for col in data_frame.columns]
    if not any("end" in column for column in column_names):
        start = data_frame.iloc[:, 1]
        bin_size = start.iloc[1] - start.iloc[0]

        end = start.shift(-1) - 1
        end.iloc[-1] = start.iloc[-1] + bin_size - 1  # fix last NaN

        data_frame.insert(2, "End", end.astype(int))
    return data_frame


def load_bed_files(input_files: str | list[str]) -> list[pd.DataFrame]:
    if isinstance(input_files, list):  # already expanded by shell
        matched_files = [str(p) for p in input_files if Path(p).exists()]
    else:
        raise TypeError("input_files must be a str or list of str objects.")

    if not matched_files:
        raise FileNotFoundError(f"No BED files found for pattern(s): {input_files}")

    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dfs = pool.map(_parse_bed_file, matched_files)
    return dfs


def _parse_hotspots_bed_file(file: str) -> pd.DataFrame:
    file_path = Path(file)
    data_frame = pd.read_csv(file_path, sep="\t", header=None)
    return data_frame


def load_hotspots_bed_files(input_files: str | list[str]) -> list[pd.DataFrame]:
    if isinstance(input_files, list):
        matched_files = [str(p) for p in input_files if Path(p).exists()]
    else:
        raise TypeError("input_files must be a str or list of str objects.")

    if not matched_files:
        raise FileNotFoundError(f"No BED files found for pattern(s): {input_files}")

    with Pool(processes=max(1, cpu_count() // 2)) as pool:
        dfs = pool.map(_parse_hotspots_bed_file, matched_files)
    return dfs


def _load_fasta_file(fasta_file: str) -> List[SeqRecord]:
    fasta_file_path = Path(fasta_file)
    return list(SeqIO.parse(fasta_file_path, "fasta"))


def parse_and_rename_fasta_file(
    fasta_file: str, mapping: Dict[str, str]
) -> Union[int, Exception]:
    loaded_fasta = _load_fasta_file(fasta_file)
    modified_fasta: List[SeqRecord] = []
    try:
        for record in loaded_fasta:
            if record.id in mapping:
                new_id = mapping[record.id]
                record.id = new_id
            modified_fasta.append(record)
        output_path = Path("./blackcap.fasta")
        return SeqIO.write(modified_fasta, output_path, "fasta")
    except Exception as e:
        print(f"Error renaming fasta file {e}")
        raise


def load_accession_mapping(file: str) -> dict[str, str]:
    try:
        dataframe = pd.read_csv(
            file, sep="\t", header=None, names=["CHROM_NUM", "ACCESSION"]
        )
        dataframe["CHROM_NUM"] = "chr" + dataframe["CHROM_NUM"].astype(str)
        return dict(
            zip(dataframe["ACCESSION"].astype(str), dataframe["CHROM_NUM"].astype(str))
        )
    except Exception as e:
        print(f"[ERROR] Failed to load accession mapping: {e}")
        sys.exit(1)


def save_per_chromosome(
    dataframe: pd.DataFrame, output_path: str, feature_name: str
) -> None:
    for chrom, sub in dataframe.groupby("CHROM"):
        chrom_str = "NA" if pd.isna(chrom) else str(chrom)
        chrom_number = re.sub(r"[a-zA-Z_-]+", "", chrom_str)
        file_path = os.path.join(output_path, f"chr_{chrom_number}_{feature_name}.bed")
        sub.to_csv(file_path, sep="\t", index=False)
