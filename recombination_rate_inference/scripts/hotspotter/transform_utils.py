import numpy as np
import pandas as pd
import re
import pyranges as pr
from pybedtools import BedTool
from typing import Literal, Union, List, Dict, overload, Any

def _convert_to_dataframe(gff_object: Union[pd.DataFrame, BedTool, pr.PyRanges]) -> pd.DataFrame:
    if isinstance(gff_object, BedTool):
        return gff_object.to_dataframe()
    elif isinstance(gff_object, pd.DataFrame):
        return gff_object.copy()
    elif isinstance(gff_object, pr.PyRanges):
        return gff_object.df
    else:
        raise TypeError("Input must be a pandas DataFrame, BedTool, or PyRanges object")


def _validate_required_columns(df: pd.DataFrame):
    required = {"chromosome", "start", "end", "# feature", "name", "symbol"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _transform_dataframe(df: pd.DataFrame, chr_pattern: str) -> pd.DataFrame:
    df = (
        df[["chromosome", "start", "end", "# feature", "name", "symbol"]]
        .assign(
            start=lambda d: d["start"].astype("int64"),
            end=lambda d: d["end"].astype("int64"),
            chromosome=lambda d: "chr" + d["chromosome"].astype(str),
        )
        .loc[lambda d: d["chromosome"].str.match(chr_pattern)]
        .rename(columns={
            "chromosome": "Chromosome",
            "start": "Start",
            "end": "End",
            "# feature": "Feature",
            "name": "Name",
            "symbol": "Symbol"
        })
        .sort_values(by=["Chromosome", "Start"])
        .reset_index(drop=True)
    )
    return df


def _convert_to_output(df: pd.DataFrame, return_type: str) -> Union[BedTool, pd.DataFrame, pr.PyRanges]:
    columns_order = ["Chromosome", "Start", "End", "Feature", "Name", "Symbol"]
    df = df[columns_order]

    match return_type:
        case "bed":
            return BedTool.from_dataframe(df)
        case "dataframe":
            return df
        case "pyranges":
            return pr.PyRanges(df)
        case _:
            raise ValueError(f"Unsupported return_type: {return_type}")

def _split_by_feature(
    df: pd.DataFrame,
    return_type: str
) -> Dict[str, Union[BedTool, pd.DataFrame, pr.PyRanges]]:
    feature_dict: Dict[str, Union[BedTool, pr.PyRanges, pd.DataFrame]] = {}
    feature_col = df.columns[3]  # 4th column: feature

    for feature in df[feature_col].unique():
        df_subset = df[df[feature_col] == feature]

        match return_type:
            case "bed":
                feature_dict[feature] = BedTool.from_dataframe(df_subset)
            case "dataframe":
                feature_dict[feature] = df_subset
            case "pyranges":
                feature_dict[feature] = pr.PyRanges(df_subset)
            case _:
                raise ValueError(f"Unknown return_type: {return_type}")
    
    return feature_dict


def _compute_midpoint(DataFrame: pd.DataFrame) -> pd.DataFrame:
    DataFrame['Midpoint'] = (DataFrame['Start'] + DataFrame['End']) // 2
    return DataFrame