import pandas as pd
import math
import re
from typing import List, Dict, Union, Tuple
from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy.typing as npt
import seaborn as sns
import pyranges as pr
from pathlib import Path
import numpy as np
from matplotlib.gridspec import GridSpec
from .transform_utils import *
from .transform import compute_midpoint


def plot_correlation_matrix(
    R_df: pd.DataFrame,
    P_df: pd.DataFrame,
    labels: list[str],
    output_dir: str,
    window_size: int,
) -> None:
    # annotate with stars
    formatted_labels = [label.replace("_", " ") for label in labels]
    f = np.vectorize(lambda p: f"{pval_to_stars(p)}" if not np.isnan(p) else "")
    annot_df = pd.DataFrame(f(P_df.values), index=P_df.index, columns=P_df.columns)

    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        R_df,
        vmin=0,
        vmax=1,
        cmap="coolwarm",
        xticklabels=formatted_labels,
        yticklabels=formatted_labels,
        square=True,
        annot=annot_df,
        fmt="",
        annot_kws={"size": 114, "fontweight": "bold"},
        cbar_kws={"label": "Pearson Correlation"},
    )
    cbar = ax.collections[0].colorbar
    cbar.set_label(  # type: ignore
        "Pearson Correlation", size=12, weight="bold", labelpad=25
    )

    plt.xticks(fontsize=14, rotation=45, ha="right", fontweight="bold")

    plt.yticks(fontsize=14, rotation=0, fontweight="bold")

    plt.tight_layout()

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    out_png = Path(output_dir) / f"feature_correlation_plot_{window_size}.png"

    plt.savefig(out_png, dpi=200, format="png")

    print(f"[INFO] wrote {out_png}")


def plot_pop_gen_stats(
    data_frame: pd.DataFrame | pr.PyRanges,
    plotLineColor: str = "black",
    outFileName: str = "plot",
    outFileFormat: str = "png",
    y_axis_title: str = "",
    title: str = "",
):
    if isinstance(data_frame, pr.PyRanges):
        data_frame = data_frame.df

    cols_lower = data_frame.columns.str.lower()
    if "midpoint" not in cols_lower:
        data_frame = compute_midpoint(data_frame)

    chromosomes = data_frame.iloc[:, 0].unique()
    n_chromosomes = len(chromosomes)
    cols = list(data_frame.columns)
    chrom_col = cols[0]

    if data_frame.shape[1] > 5:

        def y_series(df):
            return df.iloc[:, 5]
    else:
        num_cols = data_frame.select_dtypes(include=["float64"]).columns
        ycol = num_cols[-1]

        def y_series(df):
            return df[ycol]

    plt.style.use("seaborn-v0_8-white")
    chrom_lengths = [
        df["MIDPOINT"].max() - df["MIDPOINT"].min() if not df.empty else 0
        for _, df in data_frame.groupby(chrom_col)
    ]

    max_len = max(chrom_lengths) if chrom_lengths else 0
    tick_step = 10_000_000
    tickvals = list(range(0, int(max_len) + tick_step, tick_step))
    ticklabels = [str(tv // 1_000_000) for tv in tickvals]

    if n_chromosomes == 1:
        chrom = chromosomes[0]
        chrom_df = data_frame[data_frame.iloc[:, 0] == chrom]
        x = chrom_df["MIDPOINT"]
        y = y_series(chrom_df)
        fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
        ax.plot(x, y, color=plotLineColor, lw=1)
        chrom_label = re.sub(r"[-_A-Za-z]+", "", str(chrom))
        ax.set_xlabel(f"Chromosome {chrom_label}", fontsize=14, weight="bold")
        ax.set_ylabel(re.sub(r"[-_]", " ", y_axis_title), fontsize=14, weight="bold")
        ax.grid(True, alpha=0.3)
        fig.suptitle(title, fontsize=14, weight="bold")

    else:
        fig, axes = plt.subplots(
            nrows=n_chromosomes,
            ncols=1,
            figsize=(12, 2 * n_chromosomes),
            squeeze=False,
            constrained_layout=True,
            sharey=True,
        )
        axes = axes.flatten()
        for i, chrom in enumerate(chromosomes):
            chrom_df = data_frame[data_frame.iloc[:, 0] == chrom]
            x = chrom_df["MIDPOINT"]
            y = y_series(chrom_df)
            ax = axes[i]
            ax.plot(x, y, color=plotLineColor, lw=1)
            ax.set_xlabel(
                f"Chromosome {str(chrom).replace('chr', '')}",
                fontsize=14,
                fontweight="bold",
            )
            ax.set_ylabel(
                re.sub(r"[-_]", " ", y_axis_title), fontsize=14, fontweight="bold"
            )
            ax.set_xlim(-tick_step, max_len + tick_step)
            ax.set_xticks(tickvals)
            ax.set_xticklabels(ticklabels)
            ax.grid(True, alpha=0.3)
        for j in range(len(chromosomes), len(axes)):
            fig.delaxes(axes[j])
        fig.suptitle(title, fontsize=14, weight="bold")

    if outFileName and outFileFormat:
        fig.savefig(outFileName, format=outFileFormat, dpi=300, bbox_inches="tight")


def pval_to_stars(p: float):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return ""


def plot_jaccard_test_results(jaccard_results: pd.DataFrame) -> None:
    jaccard_results["error_lower"] = (
        jaccard_results["shuffled_jaccard_mean"] - jaccard_results["ci_lower_95%"]
    )
    jaccard_results["error_upper"] = (
        jaccard_results["ci_upper_95%"] - jaccard_results["shuffled_jaccard_mean"]
    )

    features = jaccard_results["feature"]
    x = np.arange(len(features))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.bar(
        x - width / 2,
        jaccard_results["shuffled_jaccard_mean"],
        width,
        yerr=[jaccard_results["error_lower"], jaccard_results["error_upper"]],
        label="Shuffled Recombination Rate Windows",
        capsize=5,
        color="lightblue",
        edgecolor="black",
    )

    ax.bar(
        x + width / 2,
        jaccard_results["observed_jaccard"],
        width,
        label="True Recombination Rate Windows",
        color="steelblue",
        edgecolor="black",
    )

    for idx, (_, row) in enumerate(jaccard_results.iterrows()):
        stars = pval_to_stars(row["p_value_two_tailed"])
        if stars:
            ax.text(
                x[idx],
                0.55,
                stars,
                ha="center",
                va="bottom",
                fontsize=16,
                color="black",
            )
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Jaccard Index")
    ax.set_xlabel("Genomic Feature")
    ax.set_title(
        "Jaccard Index Comparison of True "
        + "vs. Shuffled\nRecombination Rate Windows with Genomic Features"
    )
    ax.legend()

    plt.tight_layout(pad=3.0)
    plt.savefig("True_vs_Shuffled_Jaccard.png", dpi=300)


def plot_score_strength_per_chrom(
    input_data: pd.DataFrame | pr.PyRanges,
    output_name: str,
    group_col: Optional[str] = None,
    score_col: Optional[str] = None,
    n_bins: int = 3,
    labels: Optional[Sequence[str]] = None,
) -> None:
    pct, ranges = get_score_strength_per_chrom(
        input_data,
        group_col=group_col,
        score_col=score_col,
        n_bins=n_bins,
        labels=labels,
    )

    if pct.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No numeric scores to plot.",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_name, dpi=300, bbox_inches="tight")
        return

    # Only present groups
    chroms = [c for c in pct["CHROM"].cat.categories if (pct["CHROM"] == c).any()]
    n = len(chroms)
    cols = 4
    rows = math.ceil(n / cols) + 1
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))
    axes = np.array(axes).flatten()

    # Bin labels/colors
    if labels is None:
        unique_bins = pct["score_strength"].dropna().unique().tolist()
        # Preserve categorical order if present
        if isinstance(unique_bins, pd.Categorical):
            bins_order = list(unique_bins.categories)
        else:
            bins_order = sorted(unique_bins)
    else:
        bins_order = list(labels)

    # fallbacks for color palette length
    default_colors = ["#8EC07C", "#FABD2F", "#FB4934", "#83A598", "#D3869B", "#B8BB26"]
    colors = default_colors[: len(bins_order)]

    for i, chrom in enumerate(chroms):
        ax = axes[i]
        vals = (
            pct.query("CHROM == @chrom")
            .set_index("score_strength")["percent"]
            .reindex(bins_order, fill_value=0)
            .tolist()
        )
        ax.pie(
            vals,
            labels=bins_order,
            autopct=lambda p: f"{p:.1f}%" if p > 0 else "",
            startangle=90,
            wedgeprops=dict(linewidth=0.5, edgecolor="white"),
            colors=colors,
            textprops={"fontsize": 14, "weight": "bold"},
        )
        ax.axis("equal")
        ax.set_title(f"Chromosome {i + 1}", fontsize=16, fontweight="bold")

    # Range table
    tbl_ax = axes[-cols]
    tbl_ax.axis("off")
    tbl_ax.set_title("Score Strength Range", fontsize=14, pad=10)

    r = ranges
    if not r.empty:
        # Ensure expected columns & rounding
        r = r.rename(columns={"min": "min", "max": "max"})
        r[["min", "max"]] = r[["min", "max"]].round(2)
        t = tbl_ax.table(
            cellText=r.values.tolist(),
            colLabels=list(r.columns),
            cellLoc="left",
            loc="center",
        )
        t.auto_set_font_size(False)
        t.set_fontsize(14)
        t.scale(1.0, 1.2)

    # Turn off any leftover axes
    for j in range(len(chroms), len(axes)):
        if axes[j] is not tbl_ax:
            axes[j].axis("off")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(output_name, dpi=300, bbox_inches="tight")


def plot_score_distribution(
    data: List[pd.DataFrame], plot_name: str | None = None
) -> None:
    n = len(data)
    if n == 0:
        return

    cols = 4
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(
        nrows=rows, ncols=cols, figsize=(6 * cols, 3.5 * rows), sharey=True
    )

    axes = axes.flatten()

    for i, df in enumerate(data):
        ax = axes[i]
        x_label = df.select_dtypes(include=["float64"]).columns[0]
        x = df.select_dtypes(include=["float64"])
        ax.hist(x, bins="auto", density=True, color="firebrick")
        title = str(df.iloc[0, 0]) if not df.empty else f"Chrom {i + 1}"
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(x_label)
        if i % cols == 0:
            ax.set_ylabel("Frequency")

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.tight_layout(rect=(0, 0, 1, 0.97))

    if plot_name:
        fig.savefig(plot_name, dpi=300, bbox_inches="tight")


def plot_genome_wide_rho(
    windows: List[pd.DataFrame], plot_name: str | None = None
) -> None:
    n = len(windows)
    if n == 0:
        return
    chrom_lengths = [
        (df["MIDPOINT"].max() - df["MIDPOINT"].min()) if not df.empty else 0
        for df in windows
    ]
    max_len = max(chrom_lengths)
    tick_step = 10_000_000
    tickvals = list(range(0, int(max_len) + tick_step, tick_step))
    ticklabels = [str(tv // 1_000_000) for tv in tickvals]
    fig, axes = plt.subplots(nrows=n, ncols=1, sharey=True, figsize=(19, 2 * n))
    if n == 1:
        axes = [axes]

    for i, (ax, df) in enumerate(zip(axes, windows), start=1):
        x = df["MIDPOINT"].astype(float).values
        y = df["cM/Mb"].astype(float).values
        ax.plot(x, y, linewidth=1, color="firebrick")
        ax.set_xlim(-tick_step, max_len + tick_step)
        ax.set_xticks(tickvals)
        ax.set_xticklabels(ticklabels)
        ax.set_title(f"Chromosome {i}", loc="right", fontsize=12)
        if i == n:
            ax.set_xlabel("Genomic Position (Mb)", fontsize=18)
        ax.tick_params(direction="out", length=4)

    fig.text(0.02, 0.5, "Recombination Rate cM/Mb", rotation="vertical", fontsize=18)
    fig.subplots_adjust(top=0.95, left=0.06, right=0.98, hspace=0.5)
    if plot_name:
        fig.savefig(plot_name, dpi=300, bbox_inches="tight")


def plot_recombination_hotspots(
    smoothed_signal: npt.NDArray[np.float64],
    midpoint_bins: npt.NDArray[np.int64],
    peak_indeces: npt.NDArray[np.int64],
    chromosome_number: int | None = None,
    output_name: str = "recombination_hotspots.png",
):
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.figure(figsize=(12, 6), dpi=200)

    plt.plot(midpoint_bins, smoothed_signal, color="black", lw=1, alpha=0.4)
    plt.scatter(
        midpoint_bins[peak_indeces],
        smoothed_signal[peak_indeces],
        marker="o",
        s=10,
        color="red",
    )
    plt.xlabel(f"Chromosome {chromosome_number}", fontsize=14, fontweight="bold")
    plt.ylabel(r"Recombination Rate  $\rho$", fontsize=14, fontweight="bold")
    plt.title(f"Chromosome {chromosome_number}")
    plt.tight_layout()
    plt.savefig(output_name)


def plot_scalogram(
    coeff_matrix,
    scales,
    output_name,
    fs=1.0,
    xlabel="Sample",
    ylabel="Scale (samples)",
    max_columns=5000,
):
    S = np.abs(coeff_matrix)
    if S.shape[1] > max_columns:
        stride = S.shape[1] // max_columns
        trimmed_length = (S.shape[1] // stride) * stride
        S = (
            S[:, :trimmed_length]
            .reshape(S.shape[0], trimmed_length // stride, stride)
            .max(axis=2)
        )
    scaleogram_extent = [0, S.shape[1] / fs, scales[0], scales[-1]]
    plt.figure(figsize=(12, 4))
    plt.imshow(
        S, aspect="auto", origin="lower", extent=scaleogram_extent, cmap="viridis"
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title("Scalogram")
    plt.colorbar(label="|CWT|")
    plt.tight_layout()
    plt.savefig(output_name, dpi=200)
