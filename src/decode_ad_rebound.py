#!/usr/bin/env python3
"""
DECODE-AD Rebound Analysis.

For each AD eGene and drug target, compute:
  (a) stabilizer_load: fraction of regulating TFs that are TOPPLE stabilizers
  (b) coupling_disruption: mean disruption score across eQTL cell types

Compare drug targets vs AD eGenes. Generate scatter figure.

Output: results/decode_ad_rebound_analysis.csv
        figures/fig_decode_ad_rebound.pdf + .png
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def main():
    print("=" * 60)
    print("DECODE-AD Rebound Analysis")
    print("=" * 60)

    # -----------------------------------------------------------
    # Load data
    # -----------------------------------------------------------
    s4 = pd.read_csv('data/raw/CIMA_Table_S4.csv')
    print(f"  S4 GRN: {s4.shape}")

    topple = pd.read_csv('results/topple.csv')
    tf_class = topple.set_index(
        topple['regulon'].str.replace(r'_\+$', '', regex=True)
    )['stability_class'].to_dict()
    s4_tfs = set(s4['TF'].unique())
    print(f"  S4 TFs: {len(s4_tfs)}")
    print(f"  TOPPLE regulons: {len(topple)}")

    ad_eqtl = pd.read_csv('results/decode_ad_eqtl.csv')
    ad_egenes = sorted(ad_eqtl['phenotype_id'].unique())
    print(f"  AD eGenes: {len(ad_egenes)}")

    disruption = pd.read_csv('results/decode_ad_disruption_network.csv')
    print(f"  Disruption edges: {len(disruption)}")

    # Build per-cell-type mean disruption score
    ct_disruption = {}
    all_cts = set(disruption['cell_type_1']) | set(disruption['cell_type_2'])
    for ct in all_cts:
        edges = disruption[
            (disruption['cell_type_1'] == ct) | (disruption['cell_type_2'] == ct)
        ]
        ct_disruption[ct] = edges['disruption_score'].mean()

    # Build gene -> regulating TFs mapping from S4
    gene_tfs = s4.groupby('Gene')['TF'].apply(lambda x: set(x)).to_dict()

    # Drug targets
    drug_targets = ['IL4R', 'IL13', 'IL13RA1', 'JAK1', 'JAK2', 'JAK3', 'TYK2',
                    'IL31RA', 'TSLP', 'IL33', 'IL1RL1', 'PDE4A', 'PDE4B', 'PDE4D']

    # Load S6 for cell-type eQTL info on drug targets
    s6 = pd.read_csv('data/raw/CIMA_Table_S6.csv')
    s6_eqtl = s6[s6['analysis'] == 'cis-eQTL']

    # -----------------------------------------------------------
    # Step 1 & 2: Compute metrics for all AD eGenes + drug targets
    # -----------------------------------------------------------
    all_genes = sorted(set(ad_egenes) | set(drug_targets))
    results = []

    for gene in all_genes:
        is_ad = gene in ad_egenes
        is_drug = gene in drug_targets

        # (a) Stabilizer load
        tfs = gene_tfs.get(gene, set())
        if tfs:
            n_stab = sum(1 for tf in tfs if tf_class.get(tf, '') == 'stabilizer')
            n_destab = sum(1 for tf in tfs if tf_class.get(tf, '') == 'destabilizer')
            n_mapped = sum(1 for tf in tfs if tf in tf_class)
            stabilizer_load = n_stab / n_mapped if n_mapped > 0 else np.nan
        else:
            n_stab = n_destab = n_mapped = 0
            stabilizer_load = np.nan

        # (b) Coupling disruption: mean disruption across cell types with eQTL
        # For AD eGenes, use ad_eqtl; for drug targets, use S6
        if is_ad:
            gene_cts = ad_eqtl[ad_eqtl['phenotype_id'] == gene]['celltype'].unique()
        else:
            gene_cts = s6_eqtl[s6_eqtl['phenotype_id'] == gene]['celltype'].unique()

        if len(gene_cts) > 0:
            disrupt_vals = [ct_disruption[ct] for ct in gene_cts if ct in ct_disruption]
            coupling_disruption = np.mean(disrupt_vals) if disrupt_vals else np.nan
        else:
            coupling_disruption = np.nan

        results.append({
            'gene': gene,
            'is_ad_egene': is_ad,
            'is_drug_target': is_drug,
            'n_regulating_TFs': len(tfs),
            'n_stabilizer_TFs': n_stab,
            'n_destabilizer_TFs': n_destab,
            'n_mapped_TFs': n_mapped,
            'stabilizer_load': stabilizer_load,
            'n_eqtl_celltypes': len(gene_cts),
            'coupling_disruption': coupling_disruption,
        })

    result_df = pd.DataFrame(results)

    # -----------------------------------------------------------
    # Step 3: Compare drug targets vs AD eGenes
    # -----------------------------------------------------------
    ad_only = result_df[result_df['is_ad_egene'] & ~result_df['is_drug_target']].dropna(
        subset=['stabilizer_load', 'coupling_disruption'])
    drugs_only = result_df[result_df['is_drug_target']].dropna(
        subset=['coupling_disruption'])

    print(f"\n  AD eGenes with both metrics: {len(ad_only)}")
    print(f"  Drug targets with disruption data: {len(drugs_only)}")

    print(f"\n  Drug target details:")
    for _, row in drugs_only.sort_values('coupling_disruption').iterrows():
        print(f"    {row['gene']:12s}  stab_load={row['stabilizer_load']:.2f}"
              if not np.isnan(row['stabilizer_load'])
              else f"    {row['gene']:12s}  stab_load=N/A",
              end='')
        print(f"  disruption={row['coupling_disruption']:.3f}  "
              f"CTs={row['n_eqtl_celltypes']:.0f}  TFs={row['n_regulating_TFs']:.0f}")

    # Mann-Whitney: coupling_disruption
    if len(drugs_only) >= 2 and len(ad_only) >= 2:
        u1, p1 = stats.mannwhitneyu(
            drugs_only['coupling_disruption'].dropna(),
            ad_only['coupling_disruption'].dropna(),
            alternative='less'  # hypothesis: drugs have LOWER disruption
        )
        print(f"\n  Mann-Whitney (coupling_disruption):")
        print(f"    Drug targets mean: {drugs_only['coupling_disruption'].mean():.4f}")
        print(f"    AD eGenes mean:    {ad_only['coupling_disruption'].mean():.4f}")
        print(f"    U={u1:.0f}, P={p1:.2e} (one-tailed: drugs < AD)")

    # Mann-Whitney: stabilizer_load
    drugs_stab = drugs_only['stabilizer_load'].dropna()
    ad_stab = ad_only['stabilizer_load'].dropna()
    if len(drugs_stab) >= 2 and len(ad_stab) >= 2:
        u2, p2 = stats.mannwhitneyu(drugs_stab, ad_stab, alternative='two-sided')
        print(f"\n  Mann-Whitney (stabilizer_load):")
        print(f"    Drug targets mean: {drugs_stab.mean():.4f}")
        print(f"    AD eGenes mean:    {ad_stab.mean():.4f}")
        print(f"    U={u2:.0f}, P={p2:.2e} (two-tailed)")

    # Correlation: stabilizer_load vs coupling_disruption
    valid = result_df.dropna(subset=['stabilizer_load', 'coupling_disruption'])
    if len(valid) > 10:
        rho, p = stats.spearmanr(valid['stabilizer_load'], valid['coupling_disruption'])
        print(f"\n  Stabilizer load vs disruption correlation (n={len(valid)}):")
        print(f"    Spearman rho={rho:.3f}, P={p:.2e}")

    result_df.to_csv('results/decode_ad_rebound_analysis.csv', index=False)
    print(f"\nSaved to results/decode_ad_rebound_analysis.csv")

    # -----------------------------------------------------------
    # Step 4: Figure
    # -----------------------------------------------------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8,
    })

    fig, axes = plt.subplots(1, 2, figsize=(180/25.4, 85/25.4))

    # Panel A: Scatter — stabilizer_load vs coupling_disruption
    ax = axes[0]
    # AD eGenes (gray)
    ad_plot = ad_only.dropna(subset=['stabilizer_load', 'coupling_disruption'])
    ax.scatter(ad_plot['stabilizer_load'], ad_plot['coupling_disruption'],
               s=15, c='#BDC3C7', alpha=0.5, edgecolors='none', label='AD eGenes', zorder=2)

    # Drug targets (colored, labeled)
    drug_colors = {
        'IL4R': '#E74C3C', 'JAK1': '#3498DB', 'JAK2': '#3498DB', 'JAK3': '#3498DB',
        'TYK2': '#2ECC71', 'PDE4A': '#F39C12', 'PDE4B': '#F39C12', 'PDE4D': '#F39C12',
        'IL13': '#9B59B6', 'IL13RA1': '#9B59B6', 'IL31RA': '#E67E22',
        'TSLP': '#1ABC9C', 'IL33': '#1ABC9C', 'IL1RL1': '#1ABC9C',
    }
    drug_plot = drugs_only.dropna(subset=['stabilizer_load', 'coupling_disruption'])
    for _, row in drug_plot.iterrows():
        color = drug_colors.get(row['gene'], '#E74C3C')
        ax.scatter(row['stabilizer_load'], row['coupling_disruption'],
                   s=60, c=color, edgecolors='black', linewidth=0.8,
                   zorder=5, marker='*')
        ax.annotate(row['gene'], (row['stabilizer_load'], row['coupling_disruption']),
                    fontsize=5.5, xytext=(4, 4), textcoords='offset points',
                    fontweight='bold', color=color)

    # Drug targets without stabilizer_load (plot at x=0 margin)
    drug_no_stab = drugs_only[drugs_only['stabilizer_load'].isna() &
                               drugs_only['coupling_disruption'].notna()]
    for _, row in drug_no_stab.iterrows():
        color = drug_colors.get(row['gene'], '#E74C3C')
        ax.scatter(-0.03, row['coupling_disruption'],
                   s=60, c=color, edgecolors='black', linewidth=0.8,
                   zorder=5, marker='*')
        ax.annotate(row['gene'], (-0.03, row['coupling_disruption']),
                    fontsize=5.5, xytext=(4, 4), textcoords='offset points',
                    fontweight='bold', color=color)

    ax.set_xlabel('Stabilizer load (fraction)')
    ax.set_ylabel('Coupling disruption')
    ax.set_title('AD eGenes vs Drug Targets', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add quadrant annotations
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    mid_x = (xlim[0] + xlim[1]) / 2
    mid_y = (ylim[0] + ylim[1]) / 2
    ax.axhline(mid_y, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
    ax.axvline(mid_x, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)

    # Panel B: Box plot comparison
    ax = axes[1]
    data_box = []
    labels_box = []

    if len(ad_only) > 0:
        data_box.append(ad_only['coupling_disruption'].dropna().values)
        labels_box.append(f'AD eGenes\n(n={len(ad_only)})')
    if len(drugs_only) > 0:
        data_box.append(drugs_only['coupling_disruption'].dropna().values)
        labels_box.append(f'Drug targets\n(n={len(drugs_only)})')

    if data_box:
        bp = ax.boxplot(data_box, labels=labels_box, widths=0.5,
                        patch_artist=True, showfliers=True,
                        flierprops={'markersize': 3})
        colors_box = ['#BDC3C7', '#E74C3C']
        for patch, color in zip(bp['boxes'], colors_box[:len(data_box)]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Add individual drug target points
        if len(drugs_only) > 0:
            drug_vals = drugs_only['coupling_disruption'].dropna()
            jitter = np.random.RandomState(42).uniform(-0.1, 0.1, len(drug_vals))
            ax.scatter(np.full(len(drug_vals), 2) + jitter, drug_vals,
                       s=30, c='#E74C3C', edgecolors='black', linewidth=0.5, zorder=5)
            for val, name in zip(drug_vals, drugs_only.dropna(subset=['coupling_disruption'])['gene']):
                ax.annotate(name, (2.15, val), fontsize=5, color='#E74C3C')

    ax.set_ylabel('Coupling disruption')
    ax.set_title('Drug target positioning', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add P-value annotation
    if len(drugs_only) >= 2 and len(ad_only) >= 2:
        ax.text(0.5, 0.95, f'P={p1:.2e}', transform=ax.transAxes,
                ha='center', fontsize=7, style='italic')

    # Panel labels
    for ax_i, label in zip(axes, ['A', 'B']):
        ax_i.text(-0.15, 1.1, label, transform=ax_i.transAxes, fontsize=12,
                  fontweight='bold', va='top')

    plt.tight_layout()
    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_decode_ad_rebound.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_decode_ad_rebound.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_decode_ad_rebound.pdf + .png")


if __name__ == "__main__":
    main()
