#!/usr/bin/env python3
"""
STRATA — Compositional Architecture Analysis of Psoriasis Visium Data
=====================================================================
Applies HIRA framework concepts (TOPPLE stabilizers, coupling architecture)
to GSE202011 Visium spatial transcriptomics data (psoriasis + healthy skin).

Data: 30 Visium samples (PSO lesional/non-lesional, PSA lesional/non-lesional,
      healthy), with cell2location fibroblast deconvolution (13 subtypes)
      and immune cell scoring (12 types) per spot.

No spatial XY coordinates — analysis uses compositional (spot-level) architecture.

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import h5py
import anndata
import warnings
warnings.filterwarnings('ignore')

from utils import (setup_stdout, setup_figure_style, save_figure, scatter_panel,
                   format_sig, ensure_dirs, COLORS, RESULTS_DIR, FIG_DIR)
setup_stdout()
ensure_dirs()


# ── Google Drive path resolution ──
def find_gdrive_path():
    """Find the Google Drive base path (handles Chinese folder name)."""
    if not os.path.exists('G:/'):
        raise FileNotFoundError("G: drive not found")
    entries = os.listdir('G:/')
    drive_dirs = [d for d in entries
                  if d not in ('$RECYCLE.BIN', '.shortcut-targets-by-id',
                               '.Encrypted', 'System Volume Information')]
    if not drive_dirs:
        raise FileNotFoundError("No Google Drive folder found on G:")
    return os.path.join('G:/', drive_dirs[0], 'Fibroblast_Atlas')


BASE_DIR = find_gdrive_path()
C2L_DIR = os.path.join(BASE_DIR, 'results_colab', 'cell2location')
SCORED_DIR = os.path.join(BASE_DIR, 'results_colab', 'immune_scored')
VISIUM_DIR = os.path.join(BASE_DIR, 'data', 'GSE202011_visium')
META_FILE = os.path.join(BASE_DIR, 'data', 'GSE202011_metadata.csv')

# TOPPLE stabilizer TFs (top 18 from CIMA analysis)
STABILIZER_TFS = [
    'HSF1', 'EGR1', 'KLF9', 'JUNB', 'JUN', 'PATZ1', 'ZEB1', 'PHF1',
    'ZNF101', 'FOSB', 'CEBPD', 'FOSL2', 'CREM', 'KLF2', 'ATF3',
    'FOS', 'STAT5B', 'IRF1',
]

# Destabilizer TFs (bottom of TOPPLE ranking)
DESTABILIZER_TFS = [
    'HOXA9', 'ZBTB8A', 'GATA1', 'GATA2', 'ERG', 'PBX1', 'MYCN',
    'TAL1', 'TCF7L1', 'HOXA6',
]


def load_metadata():
    """Load sample metadata with condition and PASI scores."""
    meta = pd.read_csv(META_FILE)
    # Create a mapping from GSM to metadata
    meta['gsm'] = meta['gsm'].astype(str)
    return meta


def load_all_samples(meta):
    """Load fibroblast deconvolution + immune scores + stabilizer expression per sample."""
    print("\nLoading 30 Visium samples...")
    all_samples = []

    for _, row in meta.iterrows():
        gsm = row['gsm']
        condition = row['condition']
        pasi = row['pasi']
        patient = row['patient']
        scored_id = row['scored_id']

        # Find c2l file
        c2l_candidates = [f for f in os.listdir(C2L_DIR)
                          if f.startswith(gsm) and f.endswith('.h5ad')]
        if not c2l_candidates:
            print(f"  WARNING: No c2l file for {gsm}, skipping")
            continue

        # Find scored file
        scored_candidates = [f for f in os.listdir(SCORED_DIR)
                             if f.startswith(gsm) and f.endswith('.h5ad')]
        if not scored_candidates:
            print(f"  WARNING: No scored file for {gsm}, skipping")
            continue

        # Find raw h5 file
        h5_candidates = [f for f in os.listdir(VISIUM_DIR)
                         if f.startswith(gsm) and f.endswith('.h5')]
        if not h5_candidates:
            print(f"  WARNING: No raw h5 for {gsm}, skipping")
            continue

        # Load c2l (fibroblast deconvolution from obsm)
        c2l_path = os.path.join(C2L_DIR, c2l_candidates[0])
        adata_c2l = anndata.read_h5ad(c2l_path)
        fib_df = adata_c2l.obsm['means_cell_abundance_w_sf'].copy()
        # Clean column names
        fib_df.columns = [c.replace('meanscell_abundance_w_sf_', '')
                          for c in fib_df.columns]
        fib_df.index = adata_c2l.obs_names

        # Load immune scores
        scored_path = os.path.join(SCORED_DIR, scored_candidates[0])
        adata_scored = anndata.read_h5ad(scored_path)
        score_cols = [c for c in adata_scored.obs.columns if c.startswith('score_')]
        immune_df = adata_scored.obs[score_cols].copy()
        immune_df.columns = [c.replace('score_', '') for c in immune_df.columns]
        immune_df.index = adata_scored.obs_names

        # Load raw expression for stabilizer/destabilizer TFs
        h5_path = os.path.join(VISIUM_DIR, h5_candidates[0])
        tf_expr = _load_tf_expression(h5_path, adata_c2l.obs_names)

        # Combine on shared barcodes
        shared = fib_df.index.intersection(immune_df.index).intersection(tf_expr.index)
        combined = pd.concat([fib_df.loc[shared], immune_df.loc[shared],
                              tf_expr.loc[shared]], axis=1)
        combined['gsm'] = gsm
        combined['condition'] = condition
        combined['pasi'] = pasi
        combined['patient'] = patient

        all_samples.append(combined)
        n_spots = len(combined)
        status = f"L" if condition == 'Lesional' else ("NL" if condition == 'Non-Lesional' else "H")
        print(f"  {gsm} [{status}] {n_spots} spots")

    df = pd.concat(all_samples, axis=0, ignore_index=False)
    print(f"\nTotal: {len(df):,} spots across {len(all_samples)} samples")
    return df


def _load_tf_expression(h5_path, barcodes):
    """Load stabilizer + destabilizer TF expression from raw h5."""
    all_tfs = STABILIZER_TFS + DESTABILIZER_TFS
    with h5py.File(h5_path, 'r') as h:
        gene_names = [g.decode() for g in h['matrix/features/name'][:]]
        bc_list = [b.decode() for b in h['matrix/barcodes'][:]]
        data = h['matrix/data'][:]
        indices = h['matrix/indices'][:]
        indptr = h['matrix/indptr'][:]
        shape = h['matrix/shape'][:]

    from scipy.sparse import csc_matrix
    mat = csc_matrix((data, indices, indptr), shape=shape)

    # Find TF indices
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    tf_data = {}
    for tf in all_tfs:
        if tf in gene_to_idx:
            idx = gene_to_idx[tf]
            tf_data[tf] = np.asarray(mat[idx, :].todense()).flatten()
        else:
            tf_data[tf] = np.zeros(mat.shape[1])

    tf_df = pd.DataFrame(tf_data, index=bc_list)
    # Match to requested barcodes
    shared = tf_df.index.intersection(barcodes)
    return tf_df.loc[shared]


def compute_sample_metrics(df):
    """Compute per-sample summary metrics."""
    print("\nComputing per-sample metrics...")

    fib_cols = ['F1: Superficial', 'F2/3: Stroma_PPARG+', 'F2: Universal',
                'F3: FRC-like', 'F4: DP_HHIP+', 'F4: DS_DPEP1+', 'F4: TNN+COCH+',
                'F5: NGFR+', 'F5: RAMP1+', 'F6: Inflammatory myofibroblast',
                'F6: Myofibroblast', 'F7: Fascia-like myofibroblast', 'F_Fascia']

    immune_cols = ['CD8_TRM', 'CD4_Tconv', 'Treg', 'CD8_Teff', 'NK', 'Mast',
                   'Macrophage_M1', 'Macrophage_M2', 'DC_myeloid', 'pDC',
                   'Neutrophil', 'ILC2']

    # Identify available columns
    avail_fib = [c for c in fib_cols if c in df.columns]
    avail_imm = [c for c in immune_cols if c in df.columns]
    avail_stab = [tf for tf in STABILIZER_TFS if tf in df.columns]
    avail_dest = [tf for tf in DESTABILIZER_TFS if tf in df.columns]
    all_ct = avail_fib + avail_imm

    print(f"  Fibroblast types: {len(avail_fib)}")
    print(f"  Immune types: {len(avail_imm)}")
    print(f"  Stabilizer TFs: {len(avail_stab)}")
    print(f"  Destabilizer TFs: {len(avail_dest)}")

    records = []
    for gsm, grp in df.groupby('gsm'):
        condition = grp['condition'].iloc[0]
        pasi = grp['pasi'].iloc[0]
        patient = grp['patient'].iloc[0]
        n_spots = len(grp)

        # Mean stabilizer expression (log1p-normalized)
        stab_expr = grp[avail_stab].mean(axis=1)  # mean across TFs per spot
        mean_stab = stab_expr.mean()               # mean across spots
        # Mean destabilizer expression
        dest_expr = grp[avail_dest].mean(axis=1)
        mean_dest = dest_expr.mean()
        # Ratio
        stab_dest_ratio = mean_stab / (mean_dest + 1e-6)

        # Immune diversity: Shannon entropy of mean immune scores per sample
        # Shift scores to positive, then normalize
        imm_means = grp[avail_imm].mean(axis=0)
        imm_shifted = imm_means - imm_means.min() + 0.01  # ensure positive
        imm_prop = imm_shifted / imm_shifted.sum()
        immune_entropy = -np.sum(imm_prop * np.log2(imm_prop + 1e-12))

        # Fibroblast diversity: Shannon entropy of mean fibroblast abundances
        fib_means = grp[avail_fib].mean(axis=0)
        fib_means_pos = np.maximum(fib_means, 0.01)
        fib_prop = fib_means_pos / fib_means_pos.sum()
        fib_entropy = -np.sum(fib_prop * np.log2(fib_prop + 1e-12))

        # Compositional coupling: mean abs correlation among all 25 cell types
        ct_data = grp[all_ct].dropna(axis=1, how='all')
        if ct_data.shape[1] > 2 and ct_data.shape[0] > 10:
            corr_mat = ct_data.corr(method='spearman').values
            # Mean off-diagonal absolute correlation
            mask = ~np.eye(corr_mat.shape[0], dtype=bool)
            coupling_strength = np.nanmean(np.abs(corr_mat[mask]))
            # Coupling entropy
            abs_corr = np.abs(corr_mat[mask])
            abs_corr_pos = abs_corr + 0.01
            abs_prop = abs_corr_pos / abs_corr_pos.sum()
            coupling_entropy = -np.sum(abs_prop * np.log2(abs_prop + 1e-12))
        else:
            coupling_strength = np.nan
            coupling_entropy = np.nan

        records.append({
            'gsm': gsm,
            'condition': condition,
            'pasi': pasi,
            'patient': patient,
            'n_spots': n_spots,
            'mean_stabilizer_expr': mean_stab,
            'mean_destabilizer_expr': mean_dest,
            'stab_dest_ratio': stab_dest_ratio,
            'immune_entropy': immune_entropy,
            'fib_entropy': fib_entropy,
            'coupling_strength': coupling_strength,
            'coupling_entropy': coupling_entropy,
        })

    metrics = pd.DataFrame(records)

    # Print summary by condition
    print("\n  Per-condition summary:")
    for cond in ['Healthy', 'Non-Lesional', 'Lesional']:
        sub = metrics[metrics['condition'] == cond]
        if len(sub) == 0:
            continue
        print(f"    {cond} (n={len(sub)}): "
              f"stab_expr={sub['mean_stabilizer_expr'].mean():.3f}, "
              f"immune_H={sub['immune_entropy'].mean():.3f}, "
              f"coupling={sub['coupling_strength'].mean():.3f}")

    return metrics


def run_statistical_tests(metrics):
    """Run STRATA core statistical tests."""
    print("\n" + "=" * 70)
    print("STRATA CORE TESTS")
    print("=" * 70)

    healthy = metrics[metrics['condition'] == 'Healthy']
    lesional = metrics[metrics['condition'] == 'Lesional']
    nonlesional = metrics[metrics['condition'] == 'Non-Lesional']

    results = {}

    # Test 1: Stabilizer expression — Healthy vs Lesional
    if len(healthy) > 0 and len(lesional) > 0:
        u, p = mannwhitneyu(healthy['mean_stabilizer_expr'],
                            lesional['mean_stabilizer_expr'],
                            alternative='two-sided')
        print(f"\n  1. Stabilizer expression: Healthy vs Lesional")
        print(f"     Healthy:  median = {healthy['mean_stabilizer_expr'].median():.4f} (n={len(healthy)})")
        print(f"     Lesional: median = {lesional['mean_stabilizer_expr'].median():.4f} (n={len(lesional)})")
        print(f"     Mann-Whitney U = {u:.0f}, P = {p:.2e} {format_sig(p)}")
        results['stab_healthy_vs_lesional'] = (u, p)

    # Test 2: Stabilizer expression vs PASI (diseased samples only)
    diseased = metrics[metrics['pasi'] > 0].copy()
    if len(diseased) >= 5:
        rho, p = spearmanr(diseased['mean_stabilizer_expr'], diseased['pasi'])
        print(f"\n  2. Stabilizer expression vs PASI (diseased only)")
        print(f"     Spearman rho = {rho:.3f}, P = {p:.2e}, n = {len(diseased)}")
        results['stab_vs_pasi'] = (rho, p)

    # Test 3: Immune diversity — Healthy vs Lesional
    if len(healthy) > 0 and len(lesional) > 0:
        u, p = mannwhitneyu(healthy['immune_entropy'],
                            lesional['immune_entropy'],
                            alternative='two-sided')
        print(f"\n  3. Immune diversity: Healthy vs Lesional")
        print(f"     Healthy:  median = {healthy['immune_entropy'].median():.4f}")
        print(f"     Lesional: median = {lesional['immune_entropy'].median():.4f}")
        print(f"     Mann-Whitney U = {u:.0f}, P = {p:.2e} {format_sig(p)}")
        results['immune_div_healthy_vs_lesional'] = (u, p)

    # Test 4: Coupling strength — Healthy vs Lesional
    if len(healthy) > 0 and len(lesional) > 0:
        u, p = mannwhitneyu(healthy['coupling_strength'],
                            lesional['coupling_strength'],
                            alternative='two-sided')
        print(f"\n  4. Coupling strength: Healthy vs Lesional")
        print(f"     Healthy:  median = {healthy['coupling_strength'].median():.4f}")
        print(f"     Lesional: median = {lesional['coupling_strength'].median():.4f}")
        print(f"     Mann-Whitney U = {u:.0f}, P = {p:.2e} {format_sig(p)}")
        results['coupling_healthy_vs_lesional'] = (u, p)

    # Test 5: Coupling strength vs PASI
    if len(diseased) >= 5:
        rho, p = spearmanr(diseased['coupling_strength'], diseased['pasi'])
        print(f"\n  5. Coupling strength vs PASI (diseased only)")
        print(f"     Spearman rho = {rho:.3f}, P = {p:.2e}, n = {len(diseased)}")
        results['coupling_vs_pasi'] = (rho, p)

    # Test 6: Kruskal-Wallis across all 3 conditions
    groups = [healthy['mean_stabilizer_expr'].values,
              nonlesional['mean_stabilizer_expr'].values,
              lesional['mean_stabilizer_expr'].values]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        h_stat, p = kruskal(*groups)
        print(f"\n  6. Kruskal-Wallis (Healthy vs NL vs L): H = {h_stat:.2f}, P = {p:.2e}")
        results['kruskal_3group'] = (h_stat, p)

    # Test 7: Stabilizer/destabilizer ratio vs PASI
    if len(diseased) >= 5:
        rho, p = spearmanr(diseased['stab_dest_ratio'], diseased['pasi'])
        print(f"\n  7. Stabilizer/Destabilizer ratio vs PASI")
        print(f"     Spearman rho = {rho:.3f}, P = {p:.2e}, n = {len(diseased)}")
        results['ratio_vs_pasi'] = (rho, p)

    return results


def plot_figures(metrics, df):
    """Generate 4-panel STRATA figure."""
    setup_figure_style()

    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.30)

    cond_palette = {'Healthy': '#27ae60', 'Non-Lesional': '#3498db', 'Lesional': '#e74c3c'}

    # ==================================================================
    # Panel A: Stabilizer expression by condition (box + swarm)
    # ==================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    order = ['Healthy', 'Non-Lesional', 'Lesional']
    avail_order = [c for c in order if c in metrics['condition'].values]

    sns.boxplot(data=metrics, x='condition', y='mean_stabilizer_expr',
                order=avail_order, palette=cond_palette,
                width=0.5, fliersize=0, ax=ax_a)
    sns.stripplot(data=metrics, x='condition', y='mean_stabilizer_expr',
                  order=avail_order, palette=cond_palette,
                  size=8, alpha=0.7, jitter=0.15, ax=ax_a)

    ax_a.set_xlabel('')
    ax_a.set_ylabel('Mean Stabilizer TF Expression')
    ax_a.set_title('A  TOPPLE Stabilizer Expression by Disease Status',
                    fontweight='bold', loc='left')

    # Add significance brackets
    healthy_vals = metrics[metrics['condition'] == 'Healthy']['mean_stabilizer_expr']
    lesional_vals = metrics[metrics['condition'] == 'Lesional']['mean_stabilizer_expr']
    if len(healthy_vals) > 0 and len(lesional_vals) > 0:
        _, p_hl = mannwhitneyu(healthy_vals, lesional_vals, alternative='two-sided')
        ymax = metrics['mean_stabilizer_expr'].max()
        h_idx = avail_order.index('Healthy') if 'Healthy' in avail_order else None
        l_idx = avail_order.index('Lesional') if 'Lesional' in avail_order else None
        if h_idx is not None and l_idx is not None:
            ax_a.plot([h_idx, h_idx, l_idx, l_idx],
                      [ymax * 1.05, ymax * 1.08, ymax * 1.08, ymax * 1.05],
                      color='black', linewidth=1)
            ax_a.text((h_idx + l_idx) / 2, ymax * 1.09,
                      f'P = {p_hl:.2e}', ha='center', fontsize=8)

    # ==================================================================
    # Panel B: Immune diversity by condition
    # ==================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    sns.boxplot(data=metrics, x='condition', y='immune_entropy',
                order=avail_order, palette=cond_palette,
                width=0.5, fliersize=0, ax=ax_b)
    sns.stripplot(data=metrics, x='condition', y='immune_entropy',
                  order=avail_order, palette=cond_palette,
                  size=8, alpha=0.7, jitter=0.15, ax=ax_b)

    ax_b.set_xlabel('')
    ax_b.set_ylabel('Immune Diversity (Shannon Entropy)')
    ax_b.set_title('B  Immune Compositional Diversity',
                    fontweight='bold', loc='left')

    # Significance
    if len(healthy_vals) > 0 and len(lesional_vals) > 0:
        h_ent = metrics[metrics['condition'] == 'Healthy']['immune_entropy']
        l_ent = metrics[metrics['condition'] == 'Lesional']['immune_entropy']
        _, p_ent = mannwhitneyu(h_ent, l_ent, alternative='two-sided')
        ymax_b = metrics['immune_entropy'].max()
        if h_idx is not None and l_idx is not None:
            ax_b.plot([h_idx, h_idx, l_idx, l_idx],
                      [ymax_b * 1.02, ymax_b * 1.04, ymax_b * 1.04, ymax_b * 1.02],
                      color='black', linewidth=1)
            ax_b.text((h_idx + l_idx) / 2, ymax_b * 1.045,
                      f'P = {p_ent:.2e}', ha='center', fontsize=8)

    # ==================================================================
    # Panel C: Stabilizer expression vs PASI
    # ==================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    diseased = metrics[metrics['pasi'] > 0].copy()
    if len(diseased) >= 5:
        scatter_panel(
            ax_c,
            diseased['mean_stabilizer_expr'], diseased['pasi'],
            'Mean Stabilizer TF Expression',
            'PASI Score',
            'C  Stabilizer Expression vs Disease Severity',
            label_top=0, label_bot=0)
    else:
        ax_c.text(0.5, 0.5, 'Insufficient data', transform=ax_c.transAxes,
                  ha='center', va='center')
        ax_c.set_title('C  Stabilizer Expression vs Disease Severity',
                        fontweight='bold', loc='left')

    # Color points by condition
    for cond, color in cond_palette.items():
        sub = diseased[diseased['condition'] == cond]
        if len(sub) > 0:
            ax_c.scatter(sub['mean_stabilizer_expr'], sub['pasi'],
                        c=color, s=60, zorder=5, edgecolors='white',
                        linewidth=0.5, label=cond)
    if len(diseased) > 0:
        ax_c.legend(fontsize=7, loc='best')

    # ==================================================================
    # Panel D: Compositional coupling heatmap (healthy samples)
    # ==================================================================
    ax_d = fig.add_subplot(gs[1, 1])

    fib_cols = [c for c in df.columns if c.startswith('F') and ':' in c]
    immune_cols = ['CD8_TRM', 'CD4_Tconv', 'Treg', 'CD8_Teff', 'NK', 'Mast',
                   'Macrophage_M1', 'Macrophage_M2', 'DC_myeloid', 'pDC',
                   'Neutrophil', 'ILC2']
    avail_imm = [c for c in immune_cols if c in df.columns]
    all_ct = fib_cols + avail_imm

    healthy_spots = df[df['condition'] == 'Healthy']
    if len(healthy_spots) > 20 and len(all_ct) > 2:
        ct_data = healthy_spots[all_ct].dropna(axis=1, how='all')
        # Shorten column names for display
        short_names = []
        for c in ct_data.columns:
            if ':' in c:
                short_names.append(c.split(': ')[1][:15] if ': ' in c else c[:15])
            else:
                short_names.append(c[:15])
        corr = ct_data.corr(method='spearman')
        corr.index = short_names
        corr.columns = short_names

        sns.heatmap(corr, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={'label': 'Spearman $\\rho$', 'shrink': 0.6},
                    ax=ax_d)
        ax_d.set_title('D  Compositional Coupling (Healthy Skin)',
                        fontweight='bold', loc='left', fontsize=9)
        ax_d.tick_params(labelsize=6)
        # Add divider between fibroblast and immune blocks
        n_fib = len([c for c in ct_data.columns if c.startswith('F')])
        ax_d.axhline(y=n_fib, color='black', linewidth=1.5)
        ax_d.axvline(x=n_fib, color='black', linewidth=1.5)
    else:
        ax_d.text(0.5, 0.5, 'Insufficient healthy data',
                  transform=ax_d.transAxes, ha='center', va='center')
        ax_d.set_title('D  Compositional Coupling (Healthy Skin)',
                        fontweight='bold', loc='left')

    fig.suptitle(
        'STRATA: Compositional Architecture of Psoriasis Tissue\n'
        'Stabilizer Expression + Cell-Type Coupling in GSE202011 Visium',
        fontsize=13, fontweight='bold', y=1.01)

    save_figure(fig, FIG_DIR, 'fig5_strata_spatial')
    print("Figure saved: fig5_strata_spatial.png/pdf")


def main():
    print("=" * 70)
    print("STRATA — Compositional Architecture Analysis")
    print("GSE202011 Visium: Psoriasis + Healthy Skin (30 samples)")
    print("=" * 70)

    meta = load_metadata()
    print(f"\nMetadata: {len(meta)} samples")
    print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")

    df = load_all_samples(meta)
    metrics = compute_sample_metrics(df)
    results = run_statistical_tests(metrics)

    # Save results
    metrics.to_csv(RESULTS_DIR / "strata_sample_metrics.csv", index=False)
    print(f"\nSample metrics saved: {RESULTS_DIR / 'strata_sample_metrics.csv'}")

    plot_figures(metrics, df)

    print(f"\n{'=' * 70}")
    print("STRATA ANALYSIS COMPLETE")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
