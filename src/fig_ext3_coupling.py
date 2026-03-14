#!/usr/bin/env python3
"""
Publication figure: Extension #3 — Enhancer-Driven Coupling Network.

Four panels:
  (A) Bar chart of top/bottom cell types by mean Jaccard coupling
  (B) Three-way Mantel comparison: eQTL-caQTL, eQTL-Jaccard, caQTL-Jaccard
  (C) Regulon breadth vs TOPPLE RI scatter
  (D) Stabilizer vs destabilizer comparison (n_unique_genes, mean_ct_breadth)

Input:  results/ext3_regulatory_coupling.csv, results/ext3_regulon_breadth.csv
Output: figures/fig_ext3_coupling.pdf + .png
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


# ---------- style ------------------------------------------------------------

def setup_style():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 7,
        'axes.linewidth': 0.6,
        'axes.labelsize': 7.5,
        'axes.titlesize': 9,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'lines.linewidth': 0.8,
        'savefig.dpi': 300,
    })


BLUE = '#4878CF'
RED = '#D65F5F'
TEAL = '#5F9EA0'
GREY = '#555555'
STAB_COLOR = '#4393C3'
DEST_COLOR = '#D6604D'


# ---------- figure -----------------------------------------------------------

def make_figure(coupling_df, breadth_df, outpath):
    setup_style()
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 155 / 25.4))

    # ---- Panel A: Top/bottom cell types by mean Jaccard coupling ----
    ax = axes[0, 0]
    coupling_sorted = coupling_df.sort_values('mean_jaccard', ascending=False)
    n_show = 10
    top = coupling_sorted.head(n_show)
    bottom = coupling_sorted.tail(n_show)
    show_df = pd.concat([top, bottom])

    # Shorten cell type names for display
    labels = show_df['cell_type'].values
    short_labels = [l.replace('_', ' ') for l in labels]

    y_pos = np.arange(len(show_df))
    colors_bar = [TEAL] * n_show + [GREY] * n_show
    ax.barh(y_pos, show_df['mean_jaccard'].values, height=0.7,
            color=colors_bar, edgecolor='none', alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(short_labels, fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel('Mean Jaccard coupling')
    ax.axvline(coupling_df['mean_jaccard'].mean(), color='black', ls='--',
               lw=0.6, alpha=0.5, label=f"Overall mean = {coupling_df['mean_jaccard'].mean():.3f}")
    ax.legend(fontsize=5.5, frameon=False, loc='lower right')

    # Separator line between top and bottom
    ax.axhline(n_show - 0.5, color='black', ls=':', lw=0.4, alpha=0.5)
    ax.text(0.50, n_show - 0.8, 'Top 10', fontsize=5.5, color=TEAL,
            fontweight='bold', transform=ax.get_yaxis_transform())
    ax.text(0.50, n_show + 0.2, 'Bottom 10', fontsize=5.5, color=GREY,
            fontweight='bold', transform=ax.get_yaxis_transform())

    ax.set_title('A', fontweight='bold', loc='left', pad=4)

    # ---- Panel B: Three-way Mantel comparison ----
    ax = axes[0, 1]
    mantel_labels = ['eQTL\nvs\ncaQTL', 'eQTL\nvs\nJaccard', 'caQTL\nvs\nJaccard']
    mantel_r = [0.482, 0.388, 0.566]
    bar_colors = [BLUE, '#8B7DBE', RED]  # blue, purple blend, red

    bars = ax.bar(mantel_labels, mantel_r, width=0.55, color=bar_colors,
                  alpha=0.8, edgecolor='none')

    # Add value labels and significance
    for bar, r_val in zip(bars, mantel_r):
        y = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, y + 0.015,
                f'r = {r_val:.3f}\nP < 0.0001',
                ha='center', va='bottom', fontsize=6)

    ax.set_ylabel('Mantel correlation (r)')
    ax.set_ylim(0, 0.75)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))

    # Reference line at zero
    ax.axhline(0, color='black', lw=0.4)

    ax.set_title('B', fontweight='bold', loc='left', pad=4)

    # ---- Panel C: Regulon breadth vs TOPPLE RI scatter ----
    ax = axes[1, 0]

    stab_mask = breadth_df['stability_class'] == 'stabilizer'
    dest_mask = breadth_df['stability_class'] == 'destabilizer'
    inter_mask = breadth_df['stability_class'] == 'intermediate'

    ax.scatter(breadth_df.loc[inter_mask, 'mean_ct_breadth'],
               breadth_df.loc[inter_mask, 'mean_RI'],
               s=10, alpha=0.4, color=GREY, linewidths=0.3, edgecolors='white',
               label='Intermediate', zorder=2)
    ax.scatter(breadth_df.loc[stab_mask, 'mean_ct_breadth'],
               breadth_df.loc[stab_mask, 'mean_RI'],
               s=18, alpha=0.7, color=STAB_COLOR, linewidths=0.3, edgecolors='white',
               label='Stabilizer', zorder=3)
    ax.scatter(breadth_df.loc[dest_mask, 'mean_ct_breadth'],
               breadth_df.loc[dest_mask, 'mean_RI'],
               s=18, alpha=0.7, color=DEST_COLOR, linewidths=0.3, edgecolors='white',
               label='Destabilizer', zorder=3)

    # Fit line (all points)
    x_all = breadth_df['mean_ct_breadth'].values
    y_all = breadth_df['mean_RI'].values
    slope, intercept, _, _, _ = stats.linregress(x_all, y_all)
    x_line = np.linspace(x_all.min(), x_all.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, 'k--', lw=0.7, alpha=0.5)

    # Spearman stats
    rho, pval = stats.spearmanr(x_all, y_all)
    ax.text(0.03, 0.97, f'Spearman $\\rho$ = {rho:.3f}\nP = {pval:.1e}',
            transform=ax.transAxes, ha='left', va='top', fontsize=6,
            bbox=dict(facecolor='white', edgecolor='#cccccc', alpha=0.9, lw=0.4))

    ax.set_xlabel('Mean cell-type breadth (enhancer accessibility)')
    ax.set_ylabel('TOPPLE redistribution index (RI)')
    ax.legend(frameon=False, fontsize=5.5, loc='upper right')

    # Label a few notable regulons
    notable = ['HSF1_+', 'JUNB_+', 'EGR1_+', 'GATA1_+', 'ERG_+']
    for reg in notable:
        row = breadth_df[breadth_df['regulon_clean'] == reg]
        if len(row) == 1:
            x_pt = row['mean_ct_breadth'].values[0]
            y_pt = row['mean_RI'].values[0]
            label_text = reg.replace('_+', '')
            ax.annotate(label_text, xy=(x_pt, y_pt),
                        xytext=(5, 3), textcoords='offset points',
                        fontsize=4.5, color='black', fontstyle='italic')

    ax.set_title('C', fontweight='bold', loc='left', pad=4)

    # ---- Panel D: Stabilizer vs destabilizer comparison (grouped bar) ----
    ax = axes[1, 1]
    from matplotlib.patches import Patch

    stab_data = breadth_df[breadth_df['stability_class'] == 'stabilizer']
    dest_data = breadth_df[breadth_df['stability_class'] == 'destabilizer']

    # Two metrics as grouped bars, with values normalized for display
    metric_labels = ['Unique target\ngenes', 'Mean cell-type\nbreadth']
    stab_means = [stab_data['n_unique_genes'].mean(),
                  stab_data['mean_ct_breadth'].mean()]
    dest_means = [dest_data['n_unique_genes'].mean(),
                  dest_data['mean_ct_breadth'].mean()]
    stab_sems = [stab_data['n_unique_genes'].sem(),
                 stab_data['mean_ct_breadth'].sem()]
    dest_sems = [dest_data['n_unique_genes'].sem(),
                 dest_data['mean_ct_breadth'].sem()]

    x = np.arange(len(metric_labels))
    w = 0.32

    bars_s = ax.bar(x - w/2, stab_means, w, yerr=stab_sems, capsize=3,
                    color=STAB_COLOR, edgecolor='black', linewidth=0.5, alpha=0.75,
                    label=f'Stabilizer (n = {len(stab_data)})')
    bars_d = ax.bar(x + w/2, dest_means, w, yerr=dest_sems, capsize=3,
                    color=DEST_COLOR, edgecolor='black', linewidth=0.5, alpha=0.75,
                    label=f'Destabilizer (n = {len(dest_data)})')

    # Value labels on top of bars
    for bar, val in zip(bars_s, stab_means):
        y = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, y + stab_sems[0] * 1.2 if val > 100 else y + stab_sems[1] * 1.2,
                f'{val:.0f}' if val > 100 else f'{val:.1f}',
                ha='center', va='bottom', fontsize=5.5, color=STAB_COLOR, fontweight='bold')
    for bar, val in zip(bars_d, dest_means):
        y = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, y + dest_sems[0] * 1.2 if val > 100 else y + dest_sems[1] * 1.2,
                f'{val:.0f}' if val > 100 else f'{val:.1f}',
                ha='center', va='bottom', fontsize=5.5, color=DEST_COLOR, fontweight='bold')

    # Significance brackets
    for i, metric in enumerate(['n_unique_genes', 'mean_ct_breadth']):
        u_stat, u_pval = stats.mannwhitneyu(
            stab_data[metric].values, dest_data[metric].values, alternative='greater')
        y_top = max(stab_means[i] + stab_sems[i], dest_means[i] + dest_sems[i]) * 1.18
        ax.plot([x[i] - w/2, x[i] + w/2], [y_top, y_top], 'k-', lw=0.6)
        ax.text(x[i], y_top * 1.01, f'P = {u_pval:.1e}',
                ha='center', va='bottom', fontsize=5)

    ax.set_xticks(x)
    ax.set_xticklabels(metric_labels, fontsize=6.5)
    ax.set_ylabel('Value')
    ax.legend(frameon=False, fontsize=5.5, loc='upper left')

    # Use log scale since gene counts (~400) vs breadth (~25) differ by order of magnitude
    # Actually, keep linear but note the scale difference is visible

    ax.set_title('D', fontweight='bold', loc='left', pad=4)

    # ---- Save ----
    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    png_path = str(outpath).replace('.pdf', '.png')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {outpath}")
    print(f"Saved: {png_path}")
    plt.close(fig)


# ---------- main -------------------------------------------------------------

def main():
    coupling_path = 'results/ext3_regulatory_coupling.csv'
    breadth_path = 'results/ext3_regulon_breadth.csv'
    outpath = 'figures/fig_ext3_coupling.pdf'

    print("Loading data ...")
    coupling_df = pd.read_csv(coupling_path)
    print(f"  Regulatory coupling: {len(coupling_df)} cell types")
    print(f"  Columns: {list(coupling_df.columns)}")

    breadth_df = pd.read_csv(breadth_path)
    print(f"  Regulon breadth: {len(breadth_df)} regulons")
    print(f"  Columns: {list(breadth_df.columns)}")

    # Quick stats
    stab = breadth_df[breadth_df['stability_class'] == 'stabilizer']
    dest = breadth_df[breadth_df['stability_class'] == 'destabilizer']
    print(f"\n  Stabilizers: n={len(stab)}, mean genes={stab['n_unique_genes'].mean():.0f}, "
          f"mean breadth={stab['mean_ct_breadth'].mean():.1f}")
    print(f"  Destabilizers: n={len(dest)}, mean genes={dest['n_unique_genes'].mean():.0f}, "
          f"mean breadth={dest['mean_ct_breadth'].mean():.1f}")

    rho, pval = stats.spearmanr(breadth_df['mean_ct_breadth'], breadth_df['mean_RI'])
    print(f"  Breadth-RI Spearman: rho={rho:.3f}, P={pval:.1e}")

    print("\nGenerating figure ...")
    make_figure(coupling_df, breadth_df, outpath)
    print("Done.")


if __name__ == '__main__':
    main()
