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

def plot_correlation_matrix(
        R_df: pd.DataFrame, 
        P_df: pd.DataFrame, 
        labels: list[str], 
        output_dir: str,
        window_size: int
        ) -> None:
    # annotate with stars
    formatted_labels = [label.replace("_", " ") for label in labels]
    f = np.vectorize(
        lambda p: f"{pval_to_stars(p)}" if not np.isnan(p) else ""
        )
    annot_df = pd.DataFrame(
            f(P_df.values), 
            index=P_df.index, 
            columns=P_df.columns
            )

    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        R_df, vmin=0, vmax=1, cmap='coolwarm',
        xticklabels=formatted_labels,
        yticklabels=formatted_labels,
        square=True,
        annot=annot_df, fmt="",
        annot_kws={"size": 12, "fontweight": "bold"},
        cbar_kws={"label": "Pearson Correlation"}
    )
    cbar = ax.collections[0].colorbar
    cbar.set_label(                     # type: ignore
            "Pearson Correlation", 
            size=12, 
            weight='bold', 
            labelpad=25
            )  
    
    plt.xticks(
            fontsize=12, 
            rotation=45, 
            ha='right', 
            fontweight='bold'
            )

    plt.yticks(
            fontsize=12, 
            rotation=0, 
            fontweight='bold'
            )
    
    plt.tight_layout()
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    out_png = Path(output_dir)/f"feature_correlation_plot_{window_size}.png"
    
    plt.savefig(out_png, dpi=200, format="png")
    
    print(f"[INFO] wrote {out_png}")

def plot_pop_gen_stats(
    data_frame: pd.DataFrame | pr.PyRanges,
    plotLineColor: str = "black",
    outFileName: str = "plot",
    outFileFormat: str = "png",
    y_axis_title: str = "",
    title: str = ""
):
    # Convert PyRanges → DataFrame
    if isinstance(data_frame, pr.PyRanges):
        data_frame = data_frame.df

    sns.set_style("dark")

    # Ensure a MIDPOINT column exists
    cols = data_frame.columns
    if "midpoint" not in cols.str.lower():
        midpoint = (data_frame.iloc[:, 1] + data_frame.iloc[:, 2]) // 2
        data_frame.insert(loc=3, column="MIDPOINT", value=midpoint)

    chromosomes = data_frame.iloc[:, 0].unique()
    n_chromosomes = len(chromosomes)

    # Pick Y column
    if data_frame.shape[1] > 5:
        def y_series(df): return df.iloc[:, 5]
    else:
        num_cols = data_frame.select_dtypes(include="number").columns
        ycol = num_cols[-1]
        def y_series(df): return df[ycol]

    # === CASE 1: only one chromosome ===
    if n_chromosomes == 1:
        chrom = chromosomes[0]
        fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
        chrom_df = data_frame[data_frame.iloc[:, 0] == chrom]

        sns.lineplot(
            data=chrom_df,
            x="MIDPOINT",
            y=y_series(chrom_df),
            color=plotLineColor,
            ax=ax,
            legend=False
        )

        ax.set_xlabel(f"Chromosome {re.sub(r"[-_A-Za-z]+", "", chrom)}",
                      fontsize=9, weight="bold")
        ax.set_ylabel(re.sub(r"[-_]", " ", y_axis_title),
                      fontsize=9, weight="bold")
        ax.grid(True, alpha=0.3)
        fig.suptitle(title, fontsize=14, weight='bold')

    # === CASE 2: multiple chromosomes ===
    else:
        n_cols = 3
        n_rows = math.ceil(n_chromosomes / n_cols)
        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(10, 3.5 * n_rows),
            squeeze=False,
            constrained_layout=True,
            sharey=True
        )
        axes = axes.flatten()

        for i, chrom in enumerate(chromosomes):
            chrom_df = data_frame[data_frame.iloc[:, 0] == chrom]
            sns.lineplot(
                data=chrom_df,
                x="MIDPOINT",
                y=y_series(chrom_df),
                color=plotLineColor,
                ax=axes[i],
                legend=False
            )
            axes[i].set_xlabel(f"Chromosome {str(chrom).replace('chr', '')}",
                               fontsize=8, weight="bold")
            axes[i].set_ylabel(re.sub(r"[-_]", " ", y_axis_title),
                               fontsize=8, weight="bold")
            axes[i].grid(True, alpha=0.3)

        # Remove unused panels
        for j in range(len(chromosomes), len(axes)):
            fig.delaxes(axes[j])

        fig.suptitle(title, fontsize=14, weight='bold')

    # === Save & show ===
    if outFileName and outFileFormat:
        fig.savefig(outFileName, format=outFileFormat, dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig)

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
        
def plot_jaccard_test_results_matplotlib(
    jaccard_results: pd.DataFrame
    ) -> None:
    jaccard_results['error_lower'] = (
        jaccard_results['shuffled_jaccard_mean'] - jaccard_results['ci_lower_95%']
        )
    jaccard_results['error_upper'] = (
        jaccard_results['ci_upper_95%'] - jaccard_results['shuffled_jaccard_mean']
    )
    
    features = jaccard_results['feature']
    x = np.arange(len(features))  
    width = 0.35  

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.bar(
        x - width/2,
        jaccard_results['shuffled_jaccard_mean'],
        width,
        yerr=[jaccard_results['error_lower'], jaccard_results['error_upper']],
        label='Shuffled Recombination Rate Windows',
        capsize=5,
        color='lightblue',
        edgecolor='black'
    )

    ax.bar(
        x + width/2,
        jaccard_results['observed_jaccard'],
        width,
        label='True Recombination Rate Windows',
        color='steelblue',
        edgecolor='black'
    )

    for idx, (_, row) in enumerate(jaccard_results.iterrows()):
        stars = pval_to_stars(row['p_value_two_tailed'])
        if stars:
            ax.text(
                x[idx],
                0.55,  
                stars,
                ha='center',
                va='bottom',
                fontsize=16,
                color='black'
            )
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha='right')
    ax.set_ylim(0, 1)
    ax.set_ylabel("Jaccard Index")
    ax.set_xlabel("Genomic Feature")
    ax.set_title(
        "Jaccard Index Comparison of True "
        + 
        "vs. Shuffled\nRecombination Rate Windows with Genomic Features"
        )
    ax.legend()

    plt.tight_layout(pad=3.0)
    plt.savefig("True_vs_Shuffled_Jaccard.png", dpi=300)


def get_score_strength_per_chrom(
    input_windows: BedTool | pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not isinstance(input_windows, pd.DataFrame):
        genomic_windows = input_windows.to_dataframe(
            disable_auto_names=True, 
            names=["chrom", "start", "end", "score(cM/Mb)", "midpoint"]
            )
    else:
        genomic_windows = input_windows.iloc[:, :5].copy()
        genomic_windows.columns = ["chrom", "start", "end", "score(cM/Mb)", "midpoint"]
        
    genomic_windows["score_strength"] = pd.qcut(
        genomic_windows["score(cM/Mb)"],
        q=3,
        labels=["low", "medium", "high"]
        )
    min_max_per_label = (
        genomic_windows.groupby("score_strength")["score(cM/Mb)"]
        .agg(['min', 'max'])
        .reset_index()
        )
    min_max_per_label.to_csv("Recombination_Rate_Strength.csv", sep = "\t") 
    counts = (
        genomic_windows
        .groupby(["chrom", "score_strength"])
        .size()
        .reset_index(name="count")
    )
    total_per_chrom = genomic_windows.groupby("chrom").size().reset_index(name="total")
    merged = pd.merge(counts, total_per_chrom, on="chrom")
    merged['percent'] = (merged['count'] / merged['total']) * 100
    chrom_order = [f"chr{i}" for i in range(1, 34)]
    merged["chrom"] = pd.Categorical(merged["chrom"], categories=chrom_order, ordered=True)
    return merged.sort_values(['chrom', 'score_strength']), min_max_per_label



def plot_score_strength_per_chrom(
    input_data: pd.DataFrame, 
    save_figure: bool = False) -> None:
    _, score_ranges = get_score_strength_per_chrom(input_data)
    unique_chromosomes = list(input_data.iloc[:, 0].unique())
    chromosome_number = len(unique_chromosomes)
    column_number = 4
    pie_number = math.ceil(chromosome_number / column_number)
    row_number = pie_number + 1

    fig = plt.figure(figsize=(6 * column_number, 6 * row_number))
    gs = GridSpec(row_number, column_number, figure=fig, hspace=0.35, wspace=0.25)

    for i, chrom in enumerate(unique_chromosomes):
        r = i // column_number
        c = i % column_number
        ax = fig.add_subplot(gs[r, c])
        chromosome_data = input_data[input_data.iloc[:, 0] == chrom]
        labels = chromosome_data.iloc[:, 1].astype(str).tolist()
        values = chromosome_data.iloc[:, 4].astype(float).tolist()
        ax.pie(
            values,
            labels=labels,
            autopct=lambda p: f"{p:.1f}%" if p > 0 else "",
            startangle=90,
            wedgeprops=dict(linewidth=0.5, edgecolor="white")
        )
        ax.axis('equal')
        ax.set_title(f"Chromosome {i+1}", fontsize=12)

    table_ax = fig.add_subplot(gs[row_number - 1, :2])
    table_ax.axis("off")
    table_ax.set_title("Recombination Rate Strength Range", fontsize=14, pad=10)
    cols = list(score_ranges.columns)
    table_df = score_ranges.copy()
    for col in cols[1:]:
        if pd.api.types.is_numeric_dtype(table_df[col]):
            table_df[col] = table_df[col].round(2)
    table = table_ax.table(
        cellText=table_df.values.tolist(),
        colLabels=cols,
        cellLoc="left",
        loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.2)

    if column_number > 2:
        for c in range(2, column_number):
            ax_empty = fig.add_subplot(gs[row_number - 1, c])
            ax_empty.axis("off")

    fig.suptitle("Score Strength Percentages per Chromosome", fontsize=16, y=0.995)
    if save_figure:
        fig.savefig("Rec_Rate_Strength_per_Chromosome.png", dpi=300, bbox_inches="tight")
        fig.savefig("Rec_Rate_Strength_per_Chromosome.pdf", bbox_inches="tight")
    plt.show()


def plot_genome_wide_rho(
    windows: List[pd.DataFrame],
    plot_name: str | None = None) -> None:
    n = len(windows)
    if n == 0:
        return
    chrom_lengths = [(df["midpoint"].max() - df["midpoint"].min()) if not df.empty else 0 for df in windows]
    max_len = max(chrom_lengths)
    tick_step = 10_000_000
    tickvals = list(range(0, int(max_len) + tick_step, tick_step))
    ticklabels = [str(tv // 1_000_000) for tv in tickvals]

    fig, axes = plt.subplots(nrows=n, ncols=1, sharey=True, figsize=(19, 3 * n))
    if n == 1:
        axes = [axes]

    for i, (ax, df) in enumerate(zip(axes, windows), start=1):
        x = df["midpoint"].astype(float).values
        y = df["cM/Mb"].astype(float).values
        ax.plot(x, y, linewidth=1, color="#722f37")
        ax.set_xlim(-tick_step, max_len + tick_step)
        ax.set_xticks(tickvals)
        ax.set_xticklabels(ticklabels)
        ax.set_ylabel("cM/Mb", fontsize=10)
        ax.set_title(f"Chromosome {i}", loc="left", fontsize=11)
        if i == n:
            ax.set_xlabel("Genomic Position (Mb)", fontsize=11)
        ax.tick_params(direction="out", length=4)
        ax.set_facecolor("#f5f5f5")

    fig.suptitle("Recombination Rate over Genomic Positions", fontsize=18, y=0.995)
    fig.subplots_adjust(top=0.95, left=0.06, right=0.98, hspace=0.25)
    if plot_name:
        fig.savefig(plot_name, dpi=300, bbox_inches="tight")
    plt.show()


def plot_rho_distribution(
    data: List[pd.DataFrame],
    plot_name: str | None = None) -> None:
    n = len(data)
    if n == 0:
        return
    cols = 4
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(6 * cols, 3.5 * rows), sharey=True)
    axes = np.atleast_2d(axes)
    for i, df in enumerate(data):
        r = i // cols
        c = i % cols
        ax = axes[r, c]
        x = df["cM/Mb"].astype(float).values
        ax.hist(x, bins="auto", density=True)
        title = str(df["chrom"].iat[0]) if not df.empty else f"Chrom {i+1}"
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("cM/Mb")
        if c == 0:
            ax.set_ylabel("Frequency")

    total_slots = rows * cols
    for j in range(n, total_slots):
        r = j // cols
        c = j % cols
        axes[r, c].axis("off")

    fig.suptitle("Recombination Rate (cM/Mb) Distribution by Chromosome", fontsize=16, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    if plot_name:
        fig.savefig(plot_name, dpi=300, bbox_inches="tight")
    plt.show()

def plot_recombination_hotspots(
    smoothed_signal: npt.NDArray[np.float64],
    midpoint_bins: npt.NDArray[np.int64], 
    peak_indeces: list[int], 
    chromosome_number: int | None = None,
    output_name: str = 'recombination_hotspots.png'
    ):
    
    plt.style.use('seaborn-v0_8-darkgrid')
    plt.figure(figsize=(12,6), dpi=200)
    
    plt.plot(midpoint_bins, smoothed_signal, color='black', lw=1, alpha=0.4)
    plt.scatter(midpoint_bins[peak_indeces], smoothed_signal[peak_indeces], marker='X', s=40, color='red', edgecolors='black')
    plt.xlabel("Genomic Position")
    plt.ylabel(r'Recombination Rate  $\frac{\mathrm{cM}}{\mathrm{Mb}}$')
    plt.title(f'Chromosome {chromosome_number}')
    plt.tight_layout()
    plt.savefig(output_name)