import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import math
from plotly.subplots import make_subplots
from typing import List, Dict, Union, Tuple
from polars import Boolean
from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy.typing as npt
import numpy as np

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
        
def plot_jaccard_test_results_matplotlib(jaccard_results: pd.DataFrame) -> None:
    # Calculate error bars
    jaccard_results['error_lower'] = jaccard_results['shuffled_jaccard_mean'] - jaccard_results['ci_lower_95%']
    jaccard_results['error_upper'] = jaccard_results['ci_upper_95%'] - jaccard_results['shuffled_jaccard_mean']
    
    features = jaccard_results['feature']
    x = np.arange(len(features))  # X locations for groups
    width = 0.35  # Width of bars

    fig, ax = plt.subplots(figsize=(12, 6))

    # Bar plot for shuffled values with error bars
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

    # Bar plot for observed values
    ax.bar(
        x + width/2,
        jaccard_results['observed_jaccard'],
        width,
        label='True Recombination Rate Windows',
        color='steelblue',
        edgecolor='black'
    )

    # Add p-value stars
    for idx, (_, row) in enumerate(jaccard_results.iterrows()):
        stars = pval_to_stars(row['p_value_two_tailed'])
        if stars:
            ax.text(
                x[idx],
                0.55,  # You can dynamically place this above the bars if needed
                stars,
                ha='center',
                va='bottom',
                fontsize=16,
                color='black'
            )

    # Formatting
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha='right')
    ax.set_ylim(0, 1)
    ax.set_ylabel("Jaccard Index")
    ax.set_xlabel("Genomic Feature")
    ax.set_title("Jaccard Index Comparison of True vs. Shuffled\nRecombination Rate Windows with Genomic Features")
    ax.legend()

    plt.tight_layout()
    plt.savefig("True_vs_Shuffled_Jaccard.png", dpi=300)


def get_score_strength_per_chrom(input_windows: BedTool | pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not isinstance(input_windows, pd.DataFrame):
        genomic_windows = input_windows.to_dataframe(disable_auto_names=True, 
                                     names=["chrom", "start", "end", "score(cM/Mb)", "midpoint"])
    else:
        genomic_windows = input_windows.iloc[:, :5].copy()
        genomic_windows.columns = ["chrom", "start", "end", "score(cM/Mb)", "midpoint"]
        
    genomic_windows["score_strength"] = pd.qcut(genomic_windows["score(cM/Mb)"], q=3, labels=["low", "medium", "high"])
    min_max_per_label = genomic_windows.groupby("score_strength")["score(cM/Mb)"].agg(['min', 'max']).reset_index()
    min_max_per_label.to_csv("Recombination_Rate_Strength.csv", sep = "\t") 
    counts = genomic_windows.groupby(["chrom", "score_strength"]).size().reset_index(name="count")
    total_per_chrom = genomic_windows.groupby("chrom").size().reset_index(name="total")
    merged = pd.merge(counts, total_per_chrom, on="chrom")
    merged['percent'] = (merged['count'] / merged['total']) * 100
    chrom_order = [f"chr{i}" for i in range(1, 34)]
    merged["chrom"] = pd.Categorical(merged["chrom"], categories=chrom_order, ordered=True)
    return merged.sort_values(['chrom', 'score_strength']), min_max_per_label


def plot_score_strength_per_chrom(
    input_data: pd.DataFrame, 
    save_figure: bool = False) -> None:
    
    score_strength, score_ranges = get_score_strength_per_chrom(input_data)
    from plotly.subplots import make_subplots
    unique_chromosomes = input_data.iloc[:, 0].unique()
    chromosome_number = len(unique_chromosomes)
    pie_number = math.ceil(chromosome_number / 4)
    row_number = pie_number + 1
    column_number = 4
    specs = [[{"type": "pie"}] * column_number for _ in range(pie_number)]
    specs.append([{"type": "table"}, {"type": "table"}, {}, {}])
    subplot_titles = [f"Chromosome {i+1}" for i in range(chromosome_number)] + ["Recombination Rate Strength Range"]
    
    figure = make_subplots(
        rows = row_number, 
        cols = column_number,
        specs=specs,
        subplot_titles=subplot_titles
        )
    
    for i, chrom in enumerate(unique_chromosomes):
        row = (i // column_number) + 1
        col = (i % column_number) + 1
        
        chromosome_data = input_data[input_data.iloc[:, 0] == chrom]
        score_strength = chromosome_data.iloc[:, 1].tolist()
        percent_per_chrom = chromosome_data.iloc[:, 4].tolist()
        
        figure.add_trace(
            go.Pie(
                labels=score_strength,
                values=percent_per_chrom,
                name=f"Chromosome {i+1}",
                automargin=True,
                textposition="inside"
            ),
            row = row,
            col = col
            )
    
    
    figure.add_trace(
        go.Table(
            header = dict(
                values=["Recombination Rate(cM/Mb) Strength", "Lower Bound", "Upper Bound"],
                align ="left",
                fill_color="lightblue",
                font=dict(size=13, color="black")),
            cells=dict(
                values=[score_ranges[col].round(2).tolist() for col in score_ranges.columns],
                align="left",
                font=dict(size=13, color="black")
            ),
        ),
        )
    figure.update_layout(
        height=300 * pie_number,
        width=300 * column_number,
        title_text="Score Strength Percentages per Chromosome",
        showlegend=True  # Turn off if you want individual legends per pie
    )
    if save_figure == True:
        figure.write_html("Rec_Rate_Strength_per_Chromosome.html")
        figure.write_image("Rec_Rate_Strength_per_Chromosome.png", height=1600, width=1200, scale=3)
    
    figure.show()
    
def plot_genome_wide_rho(windows: List[pd.DataFrame], plot_name: str | None = None) -> None:
    n = len(windows)
    cols = 1
    rows = n
    
    chrom_lengths = [df["midpoint"].max() - df["midpoint"].min() for df in windows]
    max_len = max(chrom_lengths)
    
    
    figure = make_subplots(rows=rows, 
                      cols=cols, 
                      subplot_titles=[f"Chromosome {i}" for i in range(1, 34)], 
                      shared_yaxes=True,
                      shared_xaxes=False)
    
    for i, df in enumerate(windows):
        tick_step = 10_000_000
        tickvals = list(range(0, max_len + tick_step, tick_step))
        ticktext = list(str(i // 10_000_000) for i in tickvals)
        row = i + 1
        
        figure.add_trace(
            go.Scatter(x = df["midpoint"].astype(int),
                       y = df["cM/Mb"].astype(float),
                       name = f"Chromosome {i}",
                       mode='lines',
                       line = dict(width=1, color="#722f37")),
            row,
            col = 1
        )
        
        figure.update_xaxes(
            title_text = "Genomic Position (Mb)",
            title_font = dict(size = 14, color = "black"),
            range=[0 - tick_step, max_len + tick_step],
            tickvals=tickvals,
            ticktext=ticktext,
            tickfont = dict(size=12, color='black'),
            ticks='outside',
            showline=False,
            linecolor='black',
            linewidth=2,
            row=row,
            col=1,
            gridcolor='white',
            showgrid=False
        )
        
    figure.update_layout(
        height=300 * rows,
        width=1900 * cols,
        title = {
            'text':'Recombination Rate over Genomic Positions',
            'x': 0.5,
            'y': 0.999,
            'yanchor': 'top',
            'xanchor': 'center',
            'font': dict(size=26, color = "black")
            },
        showlegend=False,
        margin=dict(t=80, l=20, r=20, b=20),
        plot_bgcolor='#f5f5f5'
    )
       
    figure.update_yaxes(
        title_text="cM/Mb",
        ticks='outside',
        tickfont=dict(size=12, color='black'),
        title_font = dict(size = 14, color = "black"),
        showline=False,
        linecolor='black',
        linewidth=2,
        gridcolor='white',
        showgrid=False)
    
    if plot_name:
        figure.write_html(plot_name, auto_open=True)
    figure.show()
    
   
def plot_rho_distribution(data: List[pd.DataFrame], plot_name: str | None = None):
    n = len(data)
    cols = 4
    rows = math.ceil(n/cols)
    figure = make_subplots(rows=rows,
                           cols=cols,
                           subplot_titles=[df["chrom"].iat[0] for df in data],
                           shared_yaxes=True)
    
    for i, df in enumerate(data):
        row = i // cols + 1
        col = i % cols + 1
        
        figure.add_trace(
            go.Histogram(x=df["cM/Mb"].astype(float), name=f"{df['chrom'].iat[0]}", showlegend=False, histnorm="probability"),
            row,
            col
        )
        
    figure.update_layout(
        height=300 * rows,
        width=300 * cols,
        title_text="Recombination Rate (cM/Mb) Distribution by Chromosome",
        bargap=0.2
    )
    figure.update_xaxes(title_text="cM/Mb")
    figure.update_yaxes(title_text="Frequency")
    
    if plot_name:
        try:
            figure.write_image(plot_name)
        except Exception as e:
            raise RuntimeError(f"Failed to save plot to {plot_name}: {e}")
    
    figure.show()

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