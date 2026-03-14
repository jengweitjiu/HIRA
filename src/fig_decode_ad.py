#!/usr/bin/env python3
"""
Publication figure: DECODE-AD — Dissecting AD GWAS Loci via CIMA QTL Overlap.

Four panels:
  (A) Heatmap: 25 AD loci x top 20 cell types (eQTL gene count)
  (B) Bar chart: cell type enrichment eQTL vs caQTL (side by side, top 15)
  (C) Manhattan-style dot plot of AD loci
  (D) Cell-type effect profile for the top AD locus (rs558269137)

Input:  results/decode_ad_*.csv
Output: figures/fig_decode_ad.pdf + .png
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as mcolors


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
GREY = '#888888'
GOLD = '#E8A838'


# ---------- helpers ----------------------------------------------------------

def shorten_ct(name, max_len=18):
    """Shorten cell type names for display."""
    s = name.replace('_', ' ')
    if len(s) > max_len:
        s = s[:max_len - 1] + '.'
    return s


# ---------- Panel A: Heatmap -------------------------------------------------

def panel_a(ax, eqtl_mat, summary_df):
    """Heatmap: 25 AD loci x top 20 cell types (eQTL gene count)."""
    # Select top 20 cell types by total gene count
    ct_cols = [c for c in eqtl_mat.columns if c != 'locus_id']
    ct_totals = eqtl_mat[ct_cols].sum().sort_values(ascending=False)
    top20_ct = ct_totals.head(20).index.tolist()

    # Merge with summary to get chromosomal position for sorting
    merged = eqtl_mat.merge(summary_df[['locus_id', 'chromosome', 'base_pair_location', 'variant_id']],
                            on='locus_id', how='left')
    merged = merged.sort_values(['chromosome', 'base_pair_location'])

    # Build matrix
    mat = merged.set_index('locus_id')[top20_ct].values.astype(float)
    row_labels = []
    for _, row in merged.iterrows():
        row_labels.append(f"{row['variant_id']}")
    col_labels = [shorten_ct(c, 16) for c in top20_ct]

    # Use sqrt scale for better visibility
    mat_plot = np.sqrt(mat)
    vmax = np.sqrt(max(mat.max().max(), 1))

    im = ax.imshow(mat_plot, aspect='auto', cmap='YlOrRd',
                   vmin=0, vmax=vmax, interpolation='nearest')

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=60, ha='right', fontsize=5)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=4.5)
    ax.set_xlabel('Cell type')
    ax.set_ylabel('AD GWAS locus (lead SNP)')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02, shrink=0.8)
    # Reconstruct ticks in original scale
    tick_vals = [0, 2, 5, 10, 18]
    tick_vals = [v for v in tick_vals if v <= mat.max().max()]
    cbar.set_ticks([np.sqrt(v) for v in tick_vals])
    cbar.set_ticklabels([str(v) for v in tick_vals])
    cbar.set_label('Gene count', fontsize=6)
    cbar.ax.tick_params(labelsize=5)


# ---------- Panel B: Side-by-side bar chart ----------------------------------

def panel_b(ax, eqtl_mat, caqtl_mat):
    """Bar chart: cell type enrichment eQTL vs caQTL (top 15 by caQTL)."""
    # Compute totals per cell type
    eqtl_cols = [c for c in eqtl_mat.columns if c != 'locus_id']
    caqtl_cols = [c for c in caqtl_mat.columns if c != 'locus_id']

    eqtl_totals = eqtl_mat[eqtl_cols].sum()
    caqtl_totals = caqtl_mat[caqtl_cols].sum()

    # Union of cell types, align
    all_ct = sorted(set(eqtl_totals.index) | set(caqtl_totals.index))
    eqtl_ser = eqtl_totals.reindex(all_ct, fill_value=0)
    caqtl_ser = caqtl_totals.reindex(all_ct, fill_value=0)

    # Combined score for ranking: use caQTL totals
    combined = caqtl_ser + eqtl_ser
    top15 = combined.sort_values(ascending=False).head(15).index.tolist()

    e_vals = eqtl_ser[top15].values
    c_vals = caqtl_ser[top15].values

    y = np.arange(len(top15))
    bar_h = 0.35

    ax.barh(y - bar_h / 2, e_vals, height=bar_h, color=BLUE, alpha=0.85,
            label='eQTL genes', edgecolor='none')
    ax.barh(y + bar_h / 2, c_vals, height=bar_h, color=RED, alpha=0.85,
            label='caQTL peaks', edgecolor='none')

    labels = [shorten_ct(c, 20) for c in top15]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.5)
    ax.invert_yaxis()
    ax.set_xlabel('Total associations across 25 AD loci')
    ax.legend(loc='lower right', frameon=False, fontsize=6)

    # Highlight cMono_CD14 if present
    for i, ct in enumerate(top15):
        if 'cMono_CD14' in ct:
            ax.get_yticklabels()[i].set_fontweight('bold')
            ax.get_yticklabels()[i].set_color(RED)


# ---------- Panel C: Manhattan-style dot plot --------------------------------

def panel_c(ax, summary_df, eqtl_mat, caqtl_mat):
    """Manhattan-style dot plot of AD loci."""
    df = summary_df.copy()
    df = df.sort_values(['chromosome', 'base_pair_location'])

    # Count eQTL genes per locus
    eqtl_cols = [c for c in eqtl_mat.columns if c != 'locus_id']
    eqtl_gene_counts = eqtl_mat.set_index('locus_id')[eqtl_cols].sum(axis=1)

    df['n_eqtl'] = df['locus_id'].map(eqtl_gene_counts).fillna(0)
    df['has_eqtl'] = df['has_eqtl'].astype(bool)
    df['has_caqtl'] = df['has_caqtl'].astype(bool)

    # Build x positions: cumulative across chromosomes
    chroms = sorted(df['chromosome'].unique())
    chrom_offsets = {}
    cum = 0
    chrom_mids = {}
    for ch in chroms:
        sub = df[df['chromosome'] == ch]
        chrom_offsets[ch] = cum - sub['base_pair_location'].min()
        mn = cum
        mx = cum + (sub['base_pair_location'].max() - sub['base_pair_location'].min())
        chrom_mids[ch] = (mn + mx) / 2
        cum = mx + 50_000_000  # gap between chromosomes

    df['x'] = df.apply(lambda r: r['base_pair_location'] + chrom_offsets[r['chromosome']], axis=1)
    df['y'] = -np.log10(df['p_value'].astype(float))

    # Color: both = gold, eQTL only = blue, caQTL only = red, neither = grey
    colors = []
    for _, r in df.iterrows():
        if r['has_eqtl'] and r['has_caqtl']:
            colors.append(GOLD)
        elif r['has_eqtl']:
            colors.append(BLUE)
        elif r['has_caqtl']:
            colors.append(RED)
        else:
            colors.append(GREY)

    # Dot size proportional to eQTL gene count
    sizes = 15 + df['n_eqtl'].values * 1.5

    ax.scatter(df['x'].values, df['y'].values, s=sizes, c=colors,
               edgecolors='black', linewidths=0.3, zorder=5, alpha=0.85)

    # Genome-wide significance line
    ax.axhline(-np.log10(5e-8), color='grey', ls='--', lw=0.5, zorder=1)

    # Chromosome labels
    for ch in chroms:
        ax.text(chrom_mids[ch], -1.8, str(ch), ha='center', fontsize=5, color='#555')

    # Alternating background shading
    for i, ch in enumerate(chroms):
        sub = df[df['chromosome'] == ch]
        xmin = sub['x'].min() - 20_000_000
        xmax = sub['x'].max() + 20_000_000
        if i % 2 == 0:
            ax.axvspan(xmin, xmax, color='#f0f0f0', alpha=0.5, zorder=0)

    ax.set_ylabel('$-\\log_{10}(P_{GWAS})$')
    ax.set_xlabel('Chromosome')
    ax.set_xlim(df['x'].min() - 30_000_000, df['x'].max() + 30_000_000)
    ax.set_xticks([])

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=GOLD,
               markersize=5, markeredgecolor='black', markeredgewidth=0.3, label='Both eQTL+caQTL'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=BLUE,
               markersize=5, markeredgecolor='black', markeredgewidth=0.3, label='eQTL only'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=RED,
               markersize=5, markeredgecolor='black', markeredgewidth=0.3, label='caQTL only'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=GREY,
               markersize=5, markeredgecolor='black', markeredgewidth=0.3, label='No overlap'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False,
              fontsize=5, handletextpad=0.3, borderpad=0.3)

    # Annotate top locus
    top_row = df.loc[df['y'].idxmax()]
    ax.annotate(top_row['variant_id'], xy=(top_row['x'], top_row['y']),
                xytext=(5, 5), textcoords='offset points', fontsize=4.5,
                fontstyle='italic', color='black')


# ---------- Panel D: Cell-type effect profile for top locus ------------------

def panel_d(ax, eqtl_df):
    """Cell-type effect profile for top AD locus (rs558269137)."""
    locus_data = eqtl_df[eqtl_df['locus_id'] == 'AD_locus_1'].copy()

    if locus_data.empty:
        ax.text(0.5, 0.5, 'No data for AD_locus_1', transform=ax.transAxes,
                ha='center')
        return

    # Mean |slope| per cell type
    ct_effects = locus_data.groupby('celltype').agg(
        mean_abs_slope=('slope', lambda x: np.abs(x).mean()),
        n_genes=('phenotype_id', 'nunique'),
        min_pval=('pval_nominal', 'min')
    ).reset_index()

    ct_effects = ct_effects.sort_values('mean_abs_slope', ascending=False)

    # Show all cell types that have data
    n_show = min(20, len(ct_effects))
    show = ct_effects.head(n_show).copy()

    y = np.arange(n_show)
    # Color by significance: highlight if min_pval < 1e-10
    colors_bar = [BLUE if p < 1e-10 else '#aec7e8' for p in show['min_pval'].values]

    ax.barh(y, show['mean_abs_slope'].values, height=0.7,
            color=colors_bar, edgecolor='none', alpha=0.85)

    # Add gene count annotation
    for i, (_, row) in enumerate(show.iterrows()):
        ax.text(row['mean_abs_slope'] + 0.005, i, f"n={int(row['n_genes'])}",
                va='center', fontsize=4.5, color='#555')

    labels = [shorten_ct(c, 22) for c in show['celltype'].values]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel('Mean |eQTL slope|')
    ax.set_title('Top AD locus: rs558269137 (chr1)', fontsize=7, fontweight='bold')

    # Legend for significance coloring
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=BLUE, label='$P < 10^{-10}$'),
        Patch(facecolor='#aec7e8', label='$P \\geq 10^{-10}$'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', frameon=False, fontsize=5.5)


# ---------- main figure ------------------------------------------------------

def make_figure(eqtl_mat, caqtl_mat, summary_df, eqtl_df, outpath):
    setup_style()
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 170 / 25.4))

    panel_a(axes[0, 0], eqtl_mat, summary_df)
    panel_b(axes[0, 1], eqtl_mat, caqtl_mat)
    panel_c(axes[1, 0], summary_df, eqtl_mat, caqtl_mat)
    panel_d(axes[1, 1], eqtl_df)

    # Panel labels
    for label, ax in zip(['A', 'B', 'C', 'D'], axes.flat):
        ax.text(-0.12, 1.08, label, transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top')

    plt.tight_layout(w_pad=2.5, h_pad=2.5)

    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(outpath.with_suffix('.png'), bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved: {outpath.with_suffix('.pdf')}")
    print(f"Saved: {outpath.with_suffix('.png')}")


# ---------- CLI ---------------------------------------------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='DECODE-AD publication figure')
    parser.add_argument('--outpath', default='figures/fig_decode_ad.pdf')
    args = parser.parse_args()

    base = Path(__file__).resolve().parent.parent

    # Load data
    print("Loading data...")
    eqtl_mat = pd.read_csv(base / 'results' / 'decode_ad_locus_celltype_eqtl.csv')
    print(f"  eQTL matrix: {eqtl_mat.shape}")
    print(f"  Columns: {list(eqtl_mat.columns[:5])}... ({len(eqtl_mat.columns)} total)")

    caqtl_mat = pd.read_csv(base / 'results' / 'decode_ad_locus_celltype_caqtl.csv')
    print(f"  caQTL matrix: {caqtl_mat.shape}")

    summary_df = pd.read_csv(base / 'results' / 'decode_ad_summary.csv')
    print(f"  Summary: {summary_df.shape}")
    print(f"  Summary columns: {list(summary_df.columns)}")

    eqtl_df = pd.read_csv(base / 'results' / 'decode_ad_eqtl.csv')
    print(f"  eQTL detail: {eqtl_df.shape}")
    print(f"  eQTL columns: {list(eqtl_df.columns)}")

    # Validate locus counts
    n_loci = summary_df['locus_id'].nunique()
    print(f"\n  {n_loci} AD lead loci")
    n_eqtl_loci = summary_df[summary_df['has_eqtl'] == True].shape[0]
    n_caqtl_loci = summary_df[summary_df['has_caqtl'] == True].shape[0]
    print(f"  {n_eqtl_loci} with eQTL overlap, {n_caqtl_loci} with caQTL overlap")

    outpath = base / args.outpath
    make_figure(eqtl_mat, caqtl_mat, summary_df, eqtl_df, outpath)
    print("\nDone.")
