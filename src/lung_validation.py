#!/usr/bin/env python3
"""
Lung eQTL Validation — HIRA's 4th Independent Tissue.

Apply DGSA and SICAI to lung single-cell eQTL data from Natri et al.
(Nature Genetics 2024) — 114 donors, 38 cell types, mashR effect estimates.

This validates whether HIRA's geometric decomposition principles
(non-additivity, coupling topology) hold in a completely different organ.

Input: data/lung/41588_2024_1702_MOESM4_ESM.xlsx (Tables S5, S8)
Output: results/lung_dgsa.csv, results/lung_sicai.csv, results/lung_validation_summary.csv
        figures/fig_lung_validation.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


# Lineage assignments for lung cell types
LUNG_LINEAGE = {
    'AT2': 'Epithelial', 'AT1': 'Epithelial', 'Transitional AT2': 'Epithelial',
    'Ciliated': 'Epithelial', 'Basal': 'Epithelial', 'KRT5-/KRT17+': 'Epithelial',
    'Differentiating ciliated': 'Epithelial', 'Proliferating - Epi': 'Epithelial',
    'Secretory - SCGB1A1+/SCGB3A2+': 'Epithelial', 'Secretory - SCGB3A2+': 'Epithelial',
    'Secretory - SCGB1A1+/MUC5B+': 'Epithelial',
    'Alveolar macrophage': 'Immune', 'Monocyte-derived macrophage': 'Immune',
    'Monocyte': 'Immune', 'Interstitial macrophage': 'Immune',
    'cDC1': 'Immune', 'Mast': 'Immune', 'Plasma': 'Immune',
    'Proliferating - Imm': 'Immune', 'Inflammatory monocyte': 'Immune',
    'CD4': 'Immune', 'CD8/NKT': 'Immune', 'pDC': 'Immune',
    'moDC': 'Immune', 'NK': 'Immune', 'B cells': 'Immune', 'cDC2': 'Immune',
    'Venule': 'Endothelial', 'Systemic venous': 'Endothelial',
    'Lymphatic': 'Endothelial', 'aCap': 'Endothelial',
    'gCap': 'Endothelial', 'Arteriole': 'Endothelial',
    'Alveolar FB': 'Mesenchymal', 'Pericyte': 'Mesenchymal',
    'SMC': 'Mesenchymal', 'Adventitial FB': 'Mesenchymal',
    'Mesothelial': 'Mesenchymal',
}


def load_lung_eqtl():
    """Load lung eQTL mashR effect estimates from Table S5."""
    s5_raw = pd.read_excel(
        'data/lung/41588_2024_1702_MOESM4_ESM.xlsx',
        sheet_name='Table S5', header=1
    )
    # Row 0 has actual column names
    actual_cols = s5_raw.iloc[0].values
    s5 = s5_raw.iloc[1:].copy()
    s5.columns = actual_cols

    # Clean up
    s5 = s5.rename(columns={
        'Gene symbol': 'gene',
        'Ensembl ID': 'ensembl_id',
        'rsID': 'rsid',
        'Chr': 'chr',
        'Gene start': 'gene_start',
        'Gene end': 'gene_end',
        'Varian location': 'variant_location',
        'Cluster': 'cluster',
    })

    # Cell type columns
    cell_types = [c for c in s5.columns if c not in
                  ['gene', 'ensembl_id', 'rsid', 'chr', 'gene_start', 'gene_end',
                   'variant_location', 'cluster']]

    # Convert to numeric
    for ct in cell_types:
        s5[ct] = pd.to_numeric(s5[ct], errors='coerce')

    print(f"  Loaded lung eQTL: {len(s5)} SNP-gene pairs, {len(cell_types)} cell types")
    print(f"  Unique genes: {s5['gene'].nunique()}")
    print(f"  Cell types: {cell_types}")

    return s5, cell_types


def load_lung_coloc():
    """Load GWAS colocalization results from Table S8."""
    s8_raw = pd.read_excel(
        'data/lung/41588_2024_1702_MOESM4_ESM.xlsx',
        sheet_name='Table S8', header=1
    )
    actual_cols = s8_raw.iloc[0].values
    s8 = s8_raw.iloc[1:].copy()
    s8.columns = actual_cols

    for col in ['PP H4', 'PP H3', 'H4/H3 ratio', '# SNPs tested']:
        s8[col] = pd.to_numeric(s8[col], errors='coerce')

    print(f"  Loaded lung coloc: {len(s8)} significant colocalizations")
    print(f"  GWAS datasets: {s8['GWAS dataset'].unique()}")
    print(f"  Lineages: {s8['Lineage'].value_counts().to_dict()}")

    return s8


def step1_dgsa(eqtl_df, cell_types):
    """
    DGSA — eQTL effect-size geometry decomposition for lung.

    For each gene, decompose the vector of mashR posterior effects across
    cell types into: non-additivity, magnitude, uniformity, Gini.
    """
    print("\n" + "=" * 60)
    print("STEP 1: DGSA — eQTL Geometry in Lung")
    print("=" * 60)

    # For each gene, take the lead eQTL (highest absolute effect)
    # Group by gene, pick the SNP with highest max |effect|
    eqtl_df['max_abs_effect'] = eqtl_df[cell_types].abs().max(axis=1)
    lead_eqtl = eqtl_df.sort_values('max_abs_effect', ascending=False).groupby('gene').first().reset_index()
    print(f"  Lead eQTLs: {len(lead_eqtl)} genes")

    # Filter genes with effects in >= 3 cell types
    min_cts = 3
    effect_threshold = 0.01  # mashR posterior effect threshold
    n_active = (lead_eqtl[cell_types].abs() > effect_threshold).sum(axis=1)
    lead_active = lead_eqtl[n_active >= min_cts].copy()
    print(f"  Genes with effects in >= {min_cts} cell types: {len(lead_active)}")

    results = []
    for _, row in lead_active.iterrows():
        effects = row[cell_types].values.astype(float)
        effects = np.nan_to_num(effects, 0)
        abs_effects = np.abs(effects)

        if abs_effects.sum() == 0:
            continue

        # Magnitude: L2 norm
        magnitude = np.sqrt(np.sum(effects ** 2))

        # Uniformity: cosine similarity to uniform vector
        n = len(effects)
        uniform = np.ones(n) / np.sqrt(n)
        norm_effects = abs_effects / np.sqrt(np.sum(abs_effects ** 2)) if np.sum(abs_effects ** 2) > 0 else abs_effects
        uniformity = np.dot(norm_effects, uniform)

        # Non-additivity: 1 - uniformity (deviation from equal effects)
        non_additivity = 1 - uniformity

        # Gini coefficient
        sorted_effects = np.sort(abs_effects)
        n_eff = len(sorted_effects)
        index = np.arange(1, n_eff + 1)
        gini = (2 * np.sum(index * sorted_effects) / (n_eff * np.sum(sorted_effects))) - (n_eff + 1) / n_eff if np.sum(sorted_effects) > 0 else 0

        # Sparsity: fraction of near-zero effects
        sparsity = np.mean(abs_effects < effect_threshold)

        # Number of active cell types
        n_active_ct = int(np.sum(abs_effects > effect_threshold))

        # Direction heterogeneity: fraction of effects with opposite sign to mean
        mean_sign = np.sign(np.mean(effects))
        if mean_sign != 0:
            direction_hetero = np.mean(np.sign(effects) != mean_sign)
        else:
            direction_hetero = 0.5

        results.append({
            'gene': row['gene'],
            'ensembl_id': row['ensembl_id'],
            'rsid': row['rsid'],
            'non_additivity': non_additivity,
            'magnitude': magnitude,
            'uniformity': uniformity,
            'gini': gini,
            'sparsity': sparsity,
            'direction_hetero': direction_hetero,
            'n_active_cell_types': n_active_ct,
            'max_effect': float(abs_effects.max()),
            'mean_effect': float(abs_effects.mean()),
        })

    dgsa_df = pd.DataFrame(results)
    print(f"  DGSA scored: {len(dgsa_df)} genes")
    print(f"\n  Summary statistics:")
    print(f"    Mean non-additivity: {dgsa_df['non_additivity'].mean():.4f}")
    print(f"    Mean Gini:           {dgsa_df['gini'].mean():.4f}")
    print(f"    Mean uniformity:     {dgsa_df['uniformity'].mean():.4f}")
    print(f"    Mean sparsity:       {dgsa_df['sparsity'].mean():.4f}")
    print(f"    Mean n_active CTs:   {dgsa_df['n_active_cell_types'].mean():.1f}")

    # Compare to CIMA DGSA
    print(f"\n  Comparison with CIMA blood eQTL DGSA:")
    print(f"    CIMA mean non-additivity: 0.8673 (5,253 genes, 69 cell types)")
    print(f"    Lung mean non-additivity: {dgsa_df['non_additivity'].mean():.4f} ({len(dgsa_df)} genes, {len(cell_types)} cell types)")

    return dgsa_df


def step2_sicai(eqtl_df, cell_types):
    """
    SICAI — Inter-cellular coupling from eQTL effect sharing.

    Compute pairwise Spearman correlation of mashR effects across genes
    as a proxy for r_b (genetic effect sharing).
    """
    print("\n" + "=" * 60)
    print("STEP 2: SICAI — Cell-Type Coupling in Lung")
    print("=" * 60)

    # Build cell-type effect matrix: for each cell type, the vector of gene effects
    # Use lead eQTL per gene
    eqtl_df['max_abs_effect'] = eqtl_df[cell_types].abs().max(axis=1)
    lead = eqtl_df.sort_values('max_abs_effect', ascending=False).groupby('gene').first().reset_index()

    # Filter to genes with sufficient effects
    n_active = (lead[cell_types].abs() > 0.01).sum(axis=1)
    lead_multi = lead[n_active >= 5].copy()
    print(f"  Genes with effects in >= 5 cell types: {len(lead_multi)}")

    # Effect matrix: genes × cell types
    effect_mat = lead_multi[cell_types].values.astype(float)
    effect_mat = np.nan_to_num(effect_mat, 0)

    # Pairwise Spearman correlation (proxy for r_b)
    n = len(cell_types)
    rb_matrix = pd.DataFrame(np.ones((n, n)), index=cell_types, columns=cell_types)

    for i in range(n):
        for j in range(i + 1, n):
            rho, _ = stats.spearmanr(effect_mat[:, i], effect_mat[:, j])
            rb_matrix.iloc[i, j] = rho
            rb_matrix.iloc[j, i] = rho

    # Summary
    upper_tri = rb_matrix.values[np.triu_indices(n, k=1)]
    print(f"\n  Pairwise coupling (Spearman rho proxy for r_b):")
    print(f"    Mean: {upper_tri.mean():.3f}")
    print(f"    Median: {np.median(upper_tri):.3f}")
    print(f"    Range: [{upper_tri.min():.3f}, {upper_tri.max():.3f}]")

    # Per-cell-type coupling
    ct_coupling = {}
    for ct in cell_types:
        others = [c for c in cell_types if c != ct]
        ct_coupling[ct] = rb_matrix.loc[ct, others].mean()

    ct_df = pd.DataFrame([
        {'cell_type': ct, 'mean_coupling': v, 'lineage': LUNG_LINEAGE.get(ct, 'Unknown')}
        for ct, v in ct_coupling.items()
    ]).sort_values('mean_coupling', ascending=False)

    print(f"\n  Per-cell-type coupling (top 10):")
    for _, row in ct_df.head(10).iterrows():
        print(f"    {row['cell_type']:35s} r_b={row['mean_coupling']:.3f}  [{row['lineage']}]")

    # Within-lineage vs between-lineage coupling
    within = []
    between = []
    for i in range(n):
        for j in range(i + 1, n):
            ct1, ct2 = cell_types[i], cell_types[j]
            lin1 = LUNG_LINEAGE.get(ct1, '')
            lin2 = LUNG_LINEAGE.get(ct2, '')
            if lin1 == lin2:
                within.append(rb_matrix.iloc[i, j])
            else:
                between.append(rb_matrix.iloc[i, j])

    if within and between:
        stat, p = stats.mannwhitneyu(within, between, alternative='greater')
        print(f"\n  Within-lineage coupling:  {np.mean(within):.3f} (n={len(within)})")
        print(f"  Between-lineage coupling: {np.mean(between):.3f} (n={len(between)})")
        print(f"  Mann-Whitney P (within > between): {p:.2e}")

    # Compare to CIMA
    print(f"\n  Comparison with CIMA:")
    print(f"    CIMA mean r_b: 0.818 (69 cell types, eQTL sharing)")
    print(f"    Lung mean rho: {upper_tri.mean():.3f} ({n} cell types, mashR effect correlation)")

    return rb_matrix, ct_df


def step3_disease_correlation(dgsa_df, coloc_df):
    """
    Correlate DGSA metrics with disease relevance (coloc PP H4).
    """
    print("\n" + "=" * 60)
    print("STEP 3: Disease Correlation (DGSA vs Colocalization)")
    print("=" * 60)

    # Count number of significant colocalizations per gene
    gene_coloc = coloc_df.groupby('Gene symbol').agg(
        n_coloc=('PP H4', 'count'),
        max_pp_h4=('PP H4', 'max'),
        n_cell_types=('Cell type/eQTL dataset', 'nunique'),
        gwas_traits=('GWAS dataset', lambda x: ', '.join(x.unique())),
    ).reset_index().rename(columns={'Gene symbol': 'gene'})

    print(f"  Genes with colocalizations: {len(gene_coloc)}")
    print(f"  GWAS traits: {coloc_df['GWAS dataset'].unique()}")

    # Merge with DGSA
    merged = dgsa_df.merge(gene_coloc, on='gene', how='left')
    merged['has_coloc'] = merged['n_coloc'].notna()
    merged['n_coloc'] = merged['n_coloc'].fillna(0)
    merged['n_cell_types_coloc'] = merged['n_cell_types'].fillna(0)

    n_with = merged['has_coloc'].sum()
    n_without = (~merged['has_coloc']).sum()
    print(f"  DGSA genes with colocalization: {n_with}/{len(merged)}")

    if n_with >= 5:
        # Compare non-additivity for genes with vs without colocalization
        na_with = merged[merged['has_coloc']]['non_additivity'].values
        na_without = merged[~merged['has_coloc']]['non_additivity'].values
        stat, p = stats.mannwhitneyu(na_with, na_without, alternative='two-sided')
        print(f"\n  Non-additivity: coloc genes vs non-coloc genes:")
        print(f"    With coloc:    {np.mean(na_with):.4f} (n={len(na_with)})")
        print(f"    Without coloc: {np.mean(na_without):.4f} (n={len(na_without)})")
        print(f"    Mann-Whitney P: {p:.2e}")

        # Gini
        g_with = merged[merged['has_coloc']]['gini'].values
        g_without = merged[~merged['has_coloc']]['gini'].values
        stat2, p2 = stats.mannwhitneyu(g_with, g_without, alternative='two-sided')
        print(f"\n  Gini: coloc genes vs non-coloc genes:")
        print(f"    With coloc:    {np.mean(g_with):.4f}")
        print(f"    Without coloc: {np.mean(g_without):.4f}")
        print(f"    Mann-Whitney P: {p2:.2e}")

        # Correlation: non-additivity vs number of colocalized cell types
        coloc_genes = merged[merged['has_coloc']].copy()
        if len(coloc_genes) >= 5:
            rho, p_rho = stats.spearmanr(coloc_genes['non_additivity'],
                                          coloc_genes['n_cell_types_coloc'])
            print(f"\n  Non-additivity vs n_coloc_cell_types: rho={rho:.3f}, P={p_rho:.2e}")

    return merged


def step4_figure(dgsa_df, rb_matrix, cell_types, merged, coloc_df):
    """Generate lung validation figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})
    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 160/25.4))

    lineage_colors = {
        'Epithelial': '#E74C3C',
        'Immune': '#3498DB',
        'Endothelial': '#2ECC71',
        'Mesenchymal': '#F39C12',
    }

    # Panel A: DGSA non-additivity distribution
    ax = axes[0, 0]
    ax.hist(dgsa_df['non_additivity'], bins=50, color='#2C3E50', alpha=0.7, edgecolor='white')
    ax.axvline(dgsa_df['non_additivity'].mean(), color='#E74C3C', linestyle='--', linewidth=1.5,
               label=f"Mean={dgsa_df['non_additivity'].mean():.3f}")
    ax.axvline(0.867, color='#3498DB', linestyle=':', linewidth=1.5,
               label='CIMA blood=0.867')
    ax.set_xlabel('Non-additivity')
    ax.set_ylabel('Gene count')
    ax.set_title('DGSA: Effect-Size Geometry', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6)

    # Panel B: Coupling heatmap
    ax = axes[0, 1]
    # Order by lineage
    ct_order = []
    for lin in ['Epithelial', 'Immune', 'Endothelial', 'Mesenchymal']:
        ct_order.extend([ct for ct in cell_types if LUNG_LINEAGE.get(ct) == lin])
    rb_ordered = rb_matrix.loc[ct_order, ct_order]
    im = ax.imshow(rb_ordered.values, cmap='RdYlBu_r', vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(ct_order)))
    ax.set_xticklabels([c[:12] for c in ct_order], rotation=90, fontsize=3)
    ax.set_yticks(range(len(ct_order)))
    ax.set_yticklabels([c[:12] for c in ct_order], fontsize=3)
    plt.colorbar(im, ax=ax, shrink=0.7, label='Coupling (rho)')
    ax.set_title('SICAI: Cell-Type Coupling', fontsize=9, fontweight='bold')
    # Add lineage separators
    cumsum = 0
    for lin in ['Epithelial', 'Immune', 'Endothelial']:
        cumsum += sum(1 for ct in ct_order if LUNG_LINEAGE.get(ct) == lin)
        ax.axhline(cumsum - 0.5, color='black', linewidth=0.5)
        ax.axvline(cumsum - 0.5, color='black', linewidth=0.5)

    # Panel C: Within vs between lineage coupling
    ax = axes[1, 0]
    n = len(cell_types)
    within_vals = []
    between_vals = []
    for i in range(n):
        for j in range(i + 1, n):
            ct1, ct2 = cell_types[i], cell_types[j]
            lin1 = LUNG_LINEAGE.get(ct1, '')
            lin2 = LUNG_LINEAGE.get(ct2, '')
            val = rb_matrix.iloc[i, j]
            if lin1 == lin2:
                within_vals.append(val)
            else:
                between_vals.append(val)
    bp = ax.boxplot([within_vals, between_vals],
                    labels=['Within-lineage', 'Between-lineage'],
                    patch_artist=True, widths=0.5)
    bp['boxes'][0].set_facecolor('#E74C3C')
    bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor('#3498DB')
    bp['boxes'][1].set_alpha(0.6)
    p = stats.mannwhitneyu(within_vals, between_vals, alternative='greater')[1]
    ax.text(0.5, 0.95, f'P={p:.2e}', transform=ax.transAxes, ha='center', va='top',
            fontsize=8, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_ylabel('Coupling (rho)')
    ax.set_title('Lineage Compartmentalization', fontsize=9, fontweight='bold')

    # Panel D: Disease colocalization vs non-additivity
    ax = axes[1, 1]
    if merged is not None and 'has_coloc' in merged.columns:
        coloc_g = merged[merged['has_coloc']]['non_additivity']
        non_coloc_g = merged[~merged['has_coloc']]['non_additivity']
        bp2 = ax.boxplot([non_coloc_g.values, coloc_g.values],
                         labels=['No coloc', 'GWAS coloc'],
                         patch_artist=True, widths=0.5)
        bp2['boxes'][0].set_facecolor('#95A5A6')
        bp2['boxes'][0].set_alpha(0.6)
        bp2['boxes'][1].set_facecolor('#E74C3C')
        bp2['boxes'][1].set_alpha(0.6)
        p2 = stats.mannwhitneyu(coloc_g, non_coloc_g, alternative='two-sided')[1]
        ax.text(0.5, 0.95, f'P={p2:.2e}', transform=ax.transAxes, ha='center', va='top',
                fontsize=8, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_ylabel('Non-additivity')
    ax.set_title('Disease Genes vs Geometry', fontsize=9, fontweight='bold')

    for ax, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('HIRA Validation in Human Lung (Natri et al. 2024, 38 Cell Types)',
                 fontsize=10, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_lung_validation.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_lung_validation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("\nSaved figures/fig_lung_validation.pdf + .png")


def main():
    print("=" * 60)
    print("Lung eQTL Validation — HIRA 4th Independent Tissue")
    print("Natri et al. (Nature Genetics 2024), 114 donors, 38 cell types")
    print("=" * 60)

    # Load data
    eqtl_df, cell_types = load_lung_eqtl()
    coloc_df = load_lung_coloc()

    # Step 1: DGSA
    dgsa_df = step1_dgsa(eqtl_df, cell_types)
    dgsa_df.to_csv('results/lung_dgsa.csv', index=False)
    print(f"  Saved results/lung_dgsa.csv")

    # Step 2: SICAI
    rb_matrix, ct_df = step2_sicai(eqtl_df, cell_types)
    rb_matrix.to_csv('results/lung_sicai_coupling.csv')
    ct_df.to_csv('results/lung_sicai_per_ct.csv', index=False)
    print(f"  Saved results/lung_sicai_*.csv")

    # Step 3: Disease correlation
    merged = step3_disease_correlation(dgsa_df, coloc_df)

    # Step 4: Figure
    step4_figure(dgsa_df, rb_matrix, cell_types, merged, coloc_df)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY — Lung vs CIMA (Blood)")
    print("=" * 60)

    upper_tri = rb_matrix.values[np.triu_indices(len(cell_types), k=1)]
    summary = {
        'metric': [
            'DGSA mean non-additivity',
            'DGSA mean Gini',
            'SICAI mean coupling',
            'Within-lineage coupling',
            'Between-lineage coupling',
            'N cell types',
            'N genes (DGSA)',
        ],
        'lung': [
            f"{dgsa_df['non_additivity'].mean():.4f}",
            f"{dgsa_df['gini'].mean():.4f}",
            f"{upper_tri.mean():.3f}",
            f"{np.mean([rb_matrix.iloc[i,j] for i in range(len(cell_types)) for j in range(i+1, len(cell_types)) if LUNG_LINEAGE.get(cell_types[i]) == LUNG_LINEAGE.get(cell_types[j])]):.3f}",
            f"{np.mean([rb_matrix.iloc[i,j] for i in range(len(cell_types)) for j in range(i+1, len(cell_types)) if LUNG_LINEAGE.get(cell_types[i]) != LUNG_LINEAGE.get(cell_types[j])]):.3f}",
            str(len(cell_types)),
            str(len(dgsa_df)),
        ],
        'cima_blood': [
            '0.8673',
            '0.1151',
            '0.818',
            '~0.85',
            '~0.79',
            '69',
            '5253',
        ],
    }
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv('results/lung_validation_summary.csv', index=False)
    print(summary_df.to_string(index=False))
    print(f"\nSaved results/lung_validation_summary.csv")

    print("\n" + "=" * 60)
    print("Lung validation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
