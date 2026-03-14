#!/usr/bin/env python3
"""
DECODE-AD Extended Figure — 4-panel publication figure.

(A) Disease comparison: cell-type SMR association counts (AD vs As vs RA)
(B) Drug target cell-type resolution (IFNG match)
(C) SMR mediation gene breadth (eQTL genes × cell types)
(D) caQTL→eQTL mediation paths (top peaks → genes)

Output: figures/fig_decode_ad_extended.pdf + .png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from pathlib import Path

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
})

def load_data():
    """Load all results CSVs."""
    res = Path("results")
    disease = pd.read_csv(res / "decode_ad_disease_comparison.csv")
    drugs = pd.read_csv(res / "decode_ad_drug_targets.csv")
    mediation = pd.read_csv(res / "decode_ad_smr_mediation.csv")
    return disease, drugs, mediation


def panel_a(ax, disease_df):
    """Disease comparison: top 15 cell types by AD association count."""
    # Pivot: rows=celltype, cols=trait, vals=n_associations
    pivot = disease_df.pivot_table(
        index='celltype', columns='trait', values='n_associations', fill_value=0
    )

    # Sort by AD count, take top 15
    if 'AD' in pivot.columns:
        pivot = pivot.sort_values('AD', ascending=False).head(15)
    else:
        pivot = pivot.head(15)

    diseases = [d for d in ['AD', 'As', 'RA'] if d in pivot.columns]
    colors = {'AD': '#E74C3C', 'As': '#3498DB', 'RA': '#2ECC71'}

    x = np.arange(len(pivot))
    width = 0.25

    for i, dis in enumerate(diseases):
        offset = (i - 1) * width
        ax.bar(x + offset, pivot[dis], width, label=dis, color=colors.get(dis, '#999'),
               edgecolor='white', linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels([ct.split('_')[0] + '_' + '_'.join(ct.split('_')[1:3])
                        if len(ct.split('_')) > 2 else ct
                        for ct in pivot.index],
                       rotation=45, ha='right', fontsize=6)
    ax.set_ylabel('SMR associations')
    ax.set_title('Disease cell-type architecture', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_b(ax, drugs_df):
    """Drug target analysis — matched vs unmatched."""
    # drugs_df has only matched genes; build full target list
    all_targets = ['IL4R', 'IL13', 'IL31RA', 'JAK1', 'JAK2', 'JAK3', 'TYK2',
                   'IL5', 'IL5RA', 'TSLP', 'IL33', 'IL1RL1', 'PDE4', 'PDE4A',
                   'PDE4B', 'PDE4D', 'IL22', 'IL17A', 'IL12B', 'IFNG', 'IL2']
    matched_genes = set(drugs_df['gene'].values)

    colors = ['#E74C3C' if t in matched_genes else '#D5D8DC' for t in all_targets]

    y = np.arange(len(all_targets))
    ax.barh(y, [1]*len(all_targets), color=colors, edgecolor='white', linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(all_targets, fontsize=5.5)
    ax.set_xlim(0, 1.5)
    ax.set_xticks([])
    ax.set_title('AD drug target overlap', fontsize=9, fontweight='bold')

    # Annotate matched
    for i, t in enumerate(all_targets):
        if t in matched_genes:
            row = drugs_df[drugs_df['gene'] == t].iloc[0]
            ax.text(1.05, i, f"eQTL: {row['all_celltypes']}",
                    va='center', fontsize=6, color='#E74C3C', fontweight='bold')

    n_matched = len(matched_genes)
    legend_elements = [mpatches.Patch(facecolor='#E74C3C', label=f'Matched ({n_matched})'),
                       mpatches.Patch(facecolor='#D5D8DC', label=f'Not found ({len(all_targets)-n_matched})')]
    ax.legend(handles=legend_elements, fontsize=6, frameon=False, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()


def panel_c(ax, mediation_df):
    """SMR mediation gene breadth — horizontal bars showing CT count per gene."""
    # Filter for eQTL AD genes
    eqtl_rows = mediation_df[
        (mediation_df['QTL'] == 'eQTL') &
        (mediation_df['trait'] == 'AD')
    ].copy()

    if len(eqtl_rows) == 0:
        ax.text(0.5, 0.5, 'No eQTL SMR data', transform=ax.transAxes, ha='center')
        return

    # Count cell types per gene
    gene_cts = eqtl_rows.groupby('gene').agg(
        n_cts=('celltype', 'nunique'),
        min_p=('p_SMR', 'min')
    ).sort_values('n_cts', ascending=True)

    y = np.arange(len(gene_cts))
    colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(gene_cts)))

    ax.barh(y, gene_cts['n_cts'], color=colors, edgecolor='white', linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(gene_cts.index, fontsize=7)
    ax.set_xlabel('Number of cell types')
    ax.set_title('AD SMR mediating genes', fontsize=9, fontweight='bold')

    # Annotate P-values
    for i, (gene, row) in enumerate(gene_cts.iterrows()):
        ax.text(row['n_cts'] + 0.2, i, f"P={row['min_p']:.1e}", va='center', fontsize=5.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_d(ax, mediation_df):
    """caQTL→eQTL mediation paths — top cell types by caQTL peak count."""
    caqtl_rows = mediation_df[
        (mediation_df['QTL'] == 'caQTL') &
        (mediation_df['trait'] == 'AD')
    ].copy()

    if len(caqtl_rows) == 0:
        ax.text(0.5, 0.5, 'No caQTL SMR data', transform=ax.transAxes, ha='center')
        return

    ct_counts = caqtl_rows.groupby('celltype').size().sort_values(ascending=True).tail(10)

    y = np.arange(len(ct_counts))
    colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(ct_counts)))

    ax.barh(y, ct_counts.values, color=colors, edgecolor='white', linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([ct.replace('_', '\n', 1) if len(ct) > 15 else ct
                        for ct in ct_counts.index], fontsize=6)
    ax.set_xlabel('caQTL peaks')
    ax.set_title('AD caQTL SMR by cell type', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def main():
    disease, drugs, mediation = load_data()

    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 160/25.4))
    fig.suptitle('DECODE-AD Extended Analysis', fontsize=11, fontweight='bold', y=0.98)

    panel_a(axes[0, 0], disease)
    panel_b(axes[0, 1], drugs)
    panel_c(axes[1, 0], mediation)
    panel_d(axes[1, 1], mediation)

    # Panel labels
    for i, (ax, label) in enumerate(zip(axes.flat, ['A', 'B', 'C', 'D'])):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    fig_dir = Path("figures")
    fig_dir.mkdir(exist_ok=True)
    fig.savefig(fig_dir / "fig_decode_ad_extended.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(fig_dir / "fig_decode_ad_extended.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_decode_ad_extended.pdf + .png")


if __name__ == "__main__":
    main()
