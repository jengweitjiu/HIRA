#!/usr/bin/env python3
"""
Layer 6: Cross-Atlas Spatial Integration (STRATA-Atlas)
=======================================================
Bridge CIMA blood atlas architecture → GSE202011 skin tissue.
Do cell types that are architecturally important in CIMA blood
also redistribute preferentially in psoriatic skin?

Mapping: 69 CIMA L4 cell types → 12 GSE202011 immune categories
         (many-to-one aggregation of CIMA subtypes)

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import anndata
import warnings
warnings.filterwarnings('ignore')

from utils import (setup_stdout, setup_figure_style, save_figure, scatter_panel,
                   format_sig, check_prerequisites, COLORS, RESULTS_DIR, FIG_DIR)
setup_stdout()


# ══════════════════════════════════════════════════════════════════════
# 1. Cell type mapping: CIMA L4 → GSE202011 immune categories
# ══════════════════════════════════════════════════════════════════════
# Many-to-one: multiple CIMA subtypes map to one GSE202011 category.
# Mapping based on lineage and known biology.

CIMA_TO_GSE = {
    # GSE202011 category: [list of CIMA L4 cell types]
    'CD8_TRM': [
        'CD8_CTL_GZMB', 'CD8_CTL_GZMK', 'CD8_CTL_IFI44L',
        'CD8_Tcm_IFI44L', 'CD8_Tem_CCR7neg',
    ],
    'CD4_Tconv': [
        'CD4_Tcm_CXCR5', 'CD4_Tcm_IFI44L', 'CD4_Tem_CCR5',
        'CD4_Tem_CCR7neg', 'CD4_Tfh-like_CXCR5',
        'CD4_Th1-like_GZMK', 'CD4_Th17-like_RORC',
        'CD4_Th22-like_CCR10', 'CD4_Th_CCR4', 'CD4_Th_LMNA',
        'CD4_Th_TNFRSF11A', 'CD4_Tn_CCR7', 'CD4_CTL_GZMH',
    ],
    'Treg': [
        'CD4_Treg_FOXP3', 'CD4_Treg_FCRL3', 'CD4_Tr1_IL10',
    ],
    'CD8_Teff': [
        'CD8_Tn_CCR7',  # naive, closest available
    ],
    'NK': [
        'Mature_NK_dim_FCGR3A', 'NK_bright_XCL1',
        'Terminal_NK_dim_CD160neg', 'Transitional_NK_GZMK',
        'Proliferative_NK_MKI67', 'NKT_NCR1',
    ],
    'Mast': [],  # No mast cells in CIMA blood atlas (tissue-resident)
    'Macrophage_M1': [
        'cMono_CD14', 'cMono_CXCL10', 'cMono_IFI44L', 'cMono_IL1B',
    ],
    'Macrophage_M2': [
        'ncMono_FCGR3A', 'ncMono_C1QA', 'ncMono_IFIT1',
        'intMono_GFRA2',
    ],
    'DC_myeloid': [
        'DC1_CLEC9A', 'DC2_CD1C', 'DC_CSF2RA', 'AS_DC',
    ],
    'pDC': [
        'pDC_IRF4',
    ],
    'Neutrophil': [],  # No neutrophils in CIMA PBMC atlas
    'ILC2': [
        'ILC2_IL2RA',
    ],
}


def load_cima_metrics():
    """Load CIMA per-cell-type metrics from all layers."""
    ri_matrix = pd.read_csv(RESULTS_DIR / "topple_ri_matrix.csv", index_col=0)
    sicai = pd.read_csv(RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv", index_col=0)
    dgsa = pd.read_csv(RESULTS_DIR / "dgsa_celltype_summary.csv", index_col=0)
    rb_matrix = pd.read_csv(RESULTS_DIR / "sicai_rb_matrix.csv", index_col=0)

    # Per-cell-type TOPPLE RI
    topple_ri = ri_matrix.mean(axis=0).rename('TOPPLE_RI')
    # SICAI mean r_b
    sicai_rb = sicai['mean_rb'].rename('SICAI_mean_rb')
    # DGSA non-additivity
    dgsa_na = dgsa['mean_non_additivity'].rename('DGSA_NonAdd')
    # Disease count
    disease_n = sicai['n_disease_associations'].rename('Disease_N')

    metrics = pd.concat([topple_ri, sicai_rb, dgsa_na, disease_n], axis=1)
    return metrics, rb_matrix


def aggregate_cima_to_gse(cima_metrics):
    """Aggregate CIMA metrics to GSE202011 immune categories."""
    records = []
    for gse_cat, cima_cts in CIMA_TO_GSE.items():
        if not cima_cts:
            continue
        # Find available CIMA cell types
        avail = [ct for ct in cima_cts if ct in cima_metrics.index]
        if not avail:
            continue
        sub = cima_metrics.loc[avail]
        rec = {
            'gse_category': gse_cat,
            'n_cima_subtypes': len(avail),
            'cima_subtypes': ', '.join(avail),
        }
        for col in cima_metrics.columns:
            rec[f'mean_{col}'] = sub[col].mean()
            rec[f'max_{col}'] = sub[col].max()
        records.append(rec)

    mapped = pd.DataFrame(records).set_index('gse_category')
    return mapped


def load_gse_immune_scores():
    """Load per-sample mean immune scores from STRATA results."""
    # Use the spot-level data from STRATA
    strata = pd.read_csv(RESULTS_DIR / "strata_sample_metrics.csv")

    # We need spot-level data — reload from h5ad files
    # But we can also compute redistribution from the STRATA script's data
    # More efficient: load the scored h5ad files directly

    # Find Google Drive path
    entries = os.listdir('G:/')
    drive_dirs = [d for d in entries
                  if d not in ('$RECYCLE.BIN', '.shortcut-targets-by-id',
                               '.Encrypted', 'System Volume Information')]
    base = os.path.join('G:/', drive_dirs[0], 'Fibroblast_Atlas')
    scored_dir = os.path.join(base, 'results_colab', 'immune_scored')
    meta = pd.read_csv(os.path.join(base, 'data', 'GSE202011_metadata.csv'))

    immune_cols = ['score_CD8_TRM', 'score_CD4_Tconv', 'score_Treg',
                   'score_CD8_Teff', 'score_NK', 'score_Mast',
                   'score_Macrophage_M1', 'score_Macrophage_M2',
                   'score_DC_myeloid', 'score_pDC', 'score_Neutrophil',
                   'score_ILC2']

    # Compute per-sample mean for each immune type
    sample_means = []
    for _, row in meta.iterrows():
        gsm = row['gsm']
        condition = row['condition']
        pasi = row['pasi']
        scored_files = [f for f in os.listdir(scored_dir)
                        if f.startswith(gsm) and f.endswith('.h5ad')]
        if not scored_files:
            continue
        adata = anndata.read_h5ad(os.path.join(scored_dir, scored_files[0]))
        avail = [c for c in immune_cols if c in adata.obs.columns]
        means = adata.obs[avail].mean().to_dict()
        means['gsm'] = gsm
        means['condition'] = condition
        means['pasi'] = pasi
        sample_means.append(means)

    df = pd.DataFrame(sample_means)
    # Clean column names
    df.columns = [c.replace('score_', '') if c.startswith('score_') else c
                  for c in df.columns]
    return df


def compute_spatial_redistribution(sample_means):
    """Compute fold change of immune scores: lesional vs healthy."""
    immune_cols = [c for c in sample_means.columns
                   if c not in ('gsm', 'condition', 'pasi')]

    healthy = sample_means[sample_means['condition'] == 'Healthy']
    lesional = sample_means[sample_means['condition'] == 'Lesional']

    # Mean across samples per condition
    h_mean = healthy[immune_cols].mean()
    l_mean = lesional[immune_cols].mean()

    # Fold change (handle near-zero and negative scores)
    # Shift so minimum is positive, then compute ratio
    shift = abs(min(h_mean.min(), l_mean.min())) + 0.01
    fc = (l_mean + shift) / (h_mean + shift)
    log2fc = np.log2(fc)

    # Also compute Mann-Whitney P for each cell type
    p_values = {}
    for col in immune_cols:
        if healthy[col].std() == 0 and lesional[col].std() == 0:
            p_values[col] = 1.0
            continue
        _, p = mannwhitneyu(healthy[col], lesional[col], alternative='two-sided')
        p_values[col] = p

    redist = pd.DataFrame({
        'healthy_mean': h_mean,
        'lesional_mean': l_mean,
        'fold_change': fc,
        'log2_FC': log2fc,
        'pvalue': pd.Series(p_values),
    })
    return redist


def compute_coupling_matrices(sample_means, rb_matrix, mapped):
    """Compare CIMA genetic coupling vs GSE202011 spatial co-occurrence."""
    immune_cols = [c for c in sample_means.columns
                   if c not in ('gsm', 'condition', 'pasi')]

    # GSE202011 spatial coupling: correlation of immune scores across spots
    # Use healthy samples for clean comparison
    healthy = sample_means[sample_means['condition'] == 'Healthy']
    gse_corr = healthy[immune_cols].corr(method='spearman')

    # CIMA genetic coupling: aggregate r_b to GSE categories
    gse_cats_with_cima = [cat for cat in immune_cols
                          if cat in mapped.index]

    # Build CIMA r_b at GSE category level
    n = len(gse_cats_with_cima)
    cima_coupling = pd.DataFrame(np.nan, index=gse_cats_with_cima,
                                  columns=gse_cats_with_cima)

    for i, cat_i in enumerate(gse_cats_with_cima):
        cima_i = [ct for ct in CIMA_TO_GSE[cat_i] if ct in rb_matrix.index]
        for j, cat_j in enumerate(gse_cats_with_cima):
            if i == j:
                cima_coupling.loc[cat_i, cat_j] = 1.0
                continue
            cima_j = [ct for ct in CIMA_TO_GSE[cat_j] if ct in rb_matrix.index]
            if not cima_i or not cima_j:
                continue
            # Mean r_b between all pairs of CIMA subtypes
            vals = []
            for ci in cima_i:
                for cj in cima_j:
                    if ci in rb_matrix.index and cj in rb_matrix.columns:
                        v = rb_matrix.loc[ci, cj]
                        if not np.isnan(v):
                            vals.append(v)
            if vals:
                cima_coupling.loc[cat_i, cat_j] = np.mean(vals)

    # GSE coupling for the same categories
    gse_coupling = gse_corr.loc[gse_cats_with_cima, gse_cats_with_cima]

    return cima_coupling, gse_coupling


def mantel_test(mat1, mat2, n_perm=9999):
    """Mantel test: correlation between two distance/similarity matrices."""
    # Extract upper triangle (excluding diagonal)
    n = mat1.shape[0]
    idx = np.triu_indices(n, k=1)
    v1 = mat1.values[idx]
    v2 = mat2.values[idx]

    # Remove NaN pairs
    mask = ~(np.isnan(v1) | np.isnan(v2))
    v1 = v1[mask]
    v2 = v2[mask]

    if len(v1) < 3:
        return np.nan, 1.0, 0

    r_obs, _ = pearsonr(v1, v2)

    # Permutation test
    count = 0
    for _ in range(n_perm):
        perm = np.random.permutation(n)
        mat2_perm = mat2.values[np.ix_(perm, perm)]
        v2_perm = mat2_perm[idx][mask]
        r_perm, _ = pearsonr(v1, v2_perm)
        if r_perm >= r_obs:
            count += 1
    p_val = (count + 1) / (n_perm + 1)

    return r_obs, p_val, len(v1)


def plot_figures(mapped, redist, cima_coupling, gse_coupling,
                 mantel_r, mantel_p, mantel_n):
    """Generate 3-panel cross-atlas figure."""
    setup_figure_style()

    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, 3, wspace=0.35)

    # ==================================================================
    # Panel A: Mapping table visualization
    # ==================================================================
    ax_a = fig.add_subplot(gs[0, 0])

    # Bar chart: number of CIMA subtypes per GSE category + mean TOPPLE RI
    plot_data = mapped.sort_values('mean_TOPPLE_RI', ascending=True).copy()
    colors = plt.cm.RdYlGn_r(
        (plot_data['mean_TOPPLE_RI'] - plot_data['mean_TOPPLE_RI'].min()) /
        (plot_data['mean_TOPPLE_RI'].max() - plot_data['mean_TOPPLE_RI'].min() + 1e-12)
    )

    bars = ax_a.barh(range(len(plot_data)), plot_data['mean_TOPPLE_RI'],
                     color=colors, edgecolor='white', height=0.6)
    ax_a.set_yticks(range(len(plot_data)))
    ax_a.set_yticklabels([f"{idx} (n={int(row['n_cima_subtypes'])})"
                          for idx, row in plot_data.iterrows()],
                         fontsize=8)
    ax_a.set_xlabel('Mean TOPPLE RI (from CIMA)')
    ax_a.set_title('A  CIMA Architecture per Tissue Cell Type',
                    fontweight='bold', loc='left', fontsize=10)

    # Annotate with disease count
    for i, (idx, row) in enumerate(plot_data.iterrows()):
        if not np.isnan(row.get('mean_Disease_N', np.nan)):
            ax_a.text(row['mean_TOPPLE_RI'] + 0.00005, i,
                      f"  D={row['mean_Disease_N']:.0f}",
                      va='center', fontsize=7, color='#7f8c8d')

    # ==================================================================
    # Panel B: CIMA architecture vs spatial redistribution
    # ==================================================================
    ax_b = fig.add_subplot(gs[0, 1])

    # Join mapped CIMA metrics with redistribution scores
    joint = mapped.join(redist[['log2_FC', 'pvalue']], how='inner')

    if len(joint) >= 4:
        # Use DGSA non-additivity as x (strongest disease predictor)
        x_col = 'mean_DGSA_NonAdd'
        if x_col in joint.columns and joint[x_col].notna().sum() >= 4:
            scatter_panel(
                ax_b,
                joint[x_col], joint['log2_FC'],
                'Mean CIMA Non-Additivity',
                'Spatial Redistribution (log$_2$ FC, L vs H)',
                'B  CIMA Architecture Predicts Tissue Redistribution',
                label_top=3, label_bot=2)
        else:
            # Fall back to TOPPLE RI
            scatter_panel(
                ax_b,
                joint['mean_TOPPLE_RI'], joint['log2_FC'],
                'Mean CIMA TOPPLE RI',
                'Spatial Redistribution (log$_2$ FC, L vs H)',
                'B  CIMA Architecture Predicts Tissue Redistribution',
                label_top=3, label_bot=2)
    else:
        ax_b.text(0.5, 0.5, f'n = {len(joint)} (insufficient)',
                  transform=ax_b.transAxes, ha='center')
        ax_b.set_title('B  CIMA → Tissue Redistribution',
                        fontweight='bold', loc='left')

    # ==================================================================
    # Panel C: Coupling matrix comparison (CIMA r_b vs GSE co-occurrence)
    # ==================================================================
    ax_c = fig.add_subplot(gs[0, 2])

    # Side-by-side: use upper triangle for CIMA, lower for GSE
    shared_cats = cima_coupling.index.intersection(gse_coupling.index)
    if len(shared_cats) >= 3:
        n = len(shared_cats)
        combined = np.zeros((n, n))
        cima_vals = cima_coupling.loc[shared_cats, shared_cats].values.astype(float)
        gse_vals = gse_coupling.loc[shared_cats, shared_cats].values.astype(float)

        # Upper triangle = CIMA r_b, lower triangle = GSE spatial corr
        for i in range(n):
            for j in range(n):
                if i == j:
                    combined[i, j] = 1.0
                elif i < j:
                    combined[i, j] = cima_vals[i, j] if not np.isnan(cima_vals[i, j]) else 0
                else:
                    combined[i, j] = gse_vals[i, j] if not np.isnan(gse_vals[i, j]) else 0

        short_names = [c.replace('Macrophage_', 'Mac_') for c in shared_cats]
        comb_df = pd.DataFrame(combined, index=short_names, columns=short_names)

        sns.heatmap(comb_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={'label': 'Coupling', 'shrink': 0.6},
                    linewidths=0.5, linecolor='white', ax=ax_c)

        ax_c.set_title('C  CIMA Genetic (upper) vs Tissue Spatial (lower)',
                        fontweight='bold', loc='left', fontsize=9)
        ax_c.tick_params(labelsize=7)

        # Add Mantel test result
        ax_c.text(0.5, -0.12,
                  f'Mantel r = {mantel_r:.3f}, P = {mantel_p:.3f} '
                  f'(n = {mantel_n} pairs)',
                  transform=ax_c.transAxes, fontsize=8, ha='center',
                  color='#2c3e50', fontweight='bold')

        # Add triangle labels
        ax_c.text(0.8, 0.2, 'CIMA\n$r_b$', transform=ax_c.transAxes,
                  fontsize=8, ha='center', va='center', color='#7f8c8d',
                  fontstyle='italic')
        ax_c.text(0.2, 0.8, 'GSE\nspatial', transform=ax_c.transAxes,
                  fontsize=8, ha='center', va='center', color='#7f8c8d',
                  fontstyle='italic')
    else:
        ax_c.text(0.5, 0.5, 'Insufficient shared categories',
                  transform=ax_c.transAxes, ha='center')

    fig.suptitle(
        'Cross-Atlas Integration: CIMA Blood Architecture → GSE202011 Tissue Organization',
        fontsize=13, fontweight='bold', y=1.02)

    save_figure(fig, FIG_DIR, 'fig6_strata_atlas')
    print("Figure saved: fig6_strata_atlas.png/pdf")


def main():
    print("=" * 70)
    print("Layer 6: Cross-Atlas Spatial Integration (STRATA-Atlas)")
    print("CIMA blood architecture → GSE202011 psoriasis tissue")
    print("=" * 70)

    required = [
        "topple_ri_matrix.csv",
        "sicai_alt_metrics_disease_merged.csv",
        "sicai_rb_matrix.csv",
        "dgsa_celltype_summary.csv",
        "strata_sample_metrics.csv",
    ]
    check_prerequisites(RESULTS_DIR, required)

    # Step 1: Load CIMA metrics and create mapping
    print("\n── Step 1: Cell Type Mapping ──")
    cima_metrics, rb_matrix = load_cima_metrics()
    mapped = aggregate_cima_to_gse(cima_metrics)
    print(f"Mapped {len(mapped)} GSE202011 categories to CIMA:")
    for cat, row in mapped.iterrows():
        print(f"  {cat}: {int(row['n_cima_subtypes'])} CIMA subtypes, "
              f"mean RI={row['mean_TOPPLE_RI']:.4f}, "
              f"mean NonAdd={row.get('mean_DGSA_NonAdd', float('nan')):.3f}")

    # Unmapped categories
    unmapped = [cat for cat, cts in CIMA_TO_GSE.items() if not cts]
    if unmapped:
        print(f"\n  Unmapped (no CIMA equivalent): {', '.join(unmapped)}")

    # Step 2: Load GSE202011 immune scores
    print("\n── Step 2: GSE202011 Immune Scores ──")
    sample_means = load_gse_immune_scores()
    print(f"  Loaded {len(sample_means)} samples")

    # Step 3: Compute spatial redistribution
    print("\n── Step 3: Spatial Redistribution (Lesional vs Healthy) ──")
    redist = compute_spatial_redistribution(sample_means)
    print("\n  Immune cell type redistribution in psoriasis:")
    print(f"  {'Cell Type':<20} {'Healthy':>10} {'Lesional':>10} "
          f"{'log2FC':>8} {'P':>10}")
    print(f"  {'-'*62}")
    for ct, row in redist.sort_values('log2_FC', ascending=False).iterrows():
        sig = format_sig(row['pvalue'])
        print(f"  {ct:<20} {row['healthy_mean']:>10.4f} "
              f"{row['lesional_mean']:>10.4f} "
              f"{row['log2_FC']:>8.3f} {row['pvalue']:>10.2e} {sig}")

    # Step 4: Core test — CIMA architecture vs spatial redistribution
    print("\n── Step 4: CIMA Architecture vs Spatial Redistribution ──")
    joint = mapped.join(redist[['log2_FC', 'pvalue']], how='inner')
    print(f"  Overlapping cell types: {len(joint)}")

    for metric_col, label in [
        ('mean_TOPPLE_RI', 'TOPPLE RI'),
        ('mean_SICAI_mean_rb', 'SICAI mean r_b'),
        ('mean_DGSA_NonAdd', 'DGSA Non-Additivity'),
        ('mean_Disease_N', 'Disease Count'),
    ]:
        if metric_col in joint.columns:
            valid = joint[[metric_col, 'log2_FC']].dropna()
            if len(valid) >= 4:
                rho, p = spearmanr(valid[metric_col], valid['log2_FC'])
                print(f"  {label} vs redistribution: "
                      f"rho={rho:.3f}, P={p:.3f}, n={len(valid)} {format_sig(p)}")

    # Step 5: Coupling matrix comparison
    print("\n── Step 5: Coupling Architecture Transfer (Mantel Test) ──")
    cima_coupling, gse_coupling = compute_coupling_matrices(
        sample_means, rb_matrix, mapped)

    shared = cima_coupling.index.intersection(gse_coupling.index)
    print(f"  Shared categories for coupling comparison: {len(shared)}")

    mantel_r, mantel_p, mantel_n = mantel_test(
        cima_coupling.loc[shared, shared],
        gse_coupling.loc[shared, shared],
        n_perm=9999)
    print(f"  Mantel test: r = {mantel_r:.3f}, P = {mantel_p:.4f}, "
          f"n = {mantel_n} unique pairs")

    # Save results
    mapped.to_csv(RESULTS_DIR / "strata_atlas_mapping.csv")
    redist.to_csv(RESULTS_DIR / "strata_atlas_redistribution.csv")
    joint_out = joint.copy()
    joint_out.to_csv(RESULTS_DIR / "strata_atlas_joint.csv")
    print(f"\n  Results saved to {RESULTS_DIR}")

    # Step 6: Figures
    print("\n── Step 6: Figures ──")
    plot_figures(mapped, redist, cima_coupling, gse_coupling,
                mantel_r, mantel_p, mantel_n)

    print(f"\n{'=' * 70}")
    print("CROSS-ATLAS INTEGRATION COMPLETE")
    print(f"{'=' * 70}")
    print(f"  Mapped cell types: {len(mapped)}")
    print(f"  Mantel test (genetic vs spatial coupling): "
          f"r={mantel_r:.3f}, P={mantel_p:.4f}")


if __name__ == "__main__":
    main()
