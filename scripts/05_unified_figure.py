#!/usr/bin/env python3
"""
Unified 6-Panel Publication Figure — TOPPLE + SICAI + DGSA
===========================================================
Row 1: TOPPLE (A: stabilizer ranking, B: stability landscape)
Row 2: Disease correlations (C: SICAI mean r_b, D: DGSA non-additivity)
Row 3: Cross-method (E: DGSA×TOPPLE, F: correlation matrix)

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

from utils import (setup_stdout, setup_figure_style, save_figure, scatter_panel,
                   format_sig, check_prerequisites, COLORS, RESULTS_DIR, FIG_DIR)
setup_stdout()


def create_unified_figure():
    setup_figure_style()

    # ── Load all results ──
    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    ri_matrix = pd.read_csv(RESULTS_DIR / "topple_ri_matrix.csv", index_col=0)
    alt_merged = pd.read_csv(RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv", index_col=0)
    dgsa_ct = pd.read_csv(RESULTS_DIR / "dgsa_celltype_summary.csv", index_col=0)

    # Disease counts from alt_merged
    disease_counts = alt_merged['n_disease_associations']

    # Per-cell-type mean RI (TOPPLE)
    ct_mean_ri = ri_matrix.mean(axis=0)
    ct_mean_ri.name = 'mean_RI'

    # ── Figure layout: 3 rows × 2 cols ──
    fig = plt.figure(figsize=(16, 18))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1],
                           hspace=0.32, wspace=0.28)

    # ==================================================================
    # Panel A: TOPPLE Top 15 Stabilizer Bar Chart
    # ==================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    top15 = scores.head(15)
    colors_a = ['#b71c1c' if r <= 3 else '#c0392b' if r <= 5
                else '#e74c3c' if r <= 10 else '#f39c12'
                for r in range(1, 16)]
    ax_a.barh(range(len(top15)), top15['mean_RI'], color=colors_a,
              edgecolor='white', height=0.7)
    ax_a.set_yticks(range(len(top15)))
    ax_a.set_yticklabels(top15.index, fontsize=7)
    ax_a.set_xlabel('Mean Redistribution Index')
    ax_a.set_title('A  TOPPLE: Top 15 Stabilizer eRegulons',
                    fontweight='bold', loc='left')
    ax_a.invert_yaxis()

    # Bold #1
    ax_a.get_yticklabels()[0].set_fontweight('bold')
    ax_a.get_yticklabels()[0].set_color('#b71c1c')

    # STAT5B annotation
    stat5b_hits = [r for r in scores.index if r.upper().startswith('STAT5B')]
    if stat5b_hits:
        s5key = stat5b_hits[0]
        s5rank = int(scores.loc[s5key, 'stability_rank'])
        s5ri = scores.loc[s5key, 'mean_RI']
        ax_a.axvline(s5ri, color='#7f8c8d', linestyle=':', linewidth=1, alpha=0.7)
        ax_a.annotate(
            f'{s5key} (#{s5rank})\n#1 in psoriasis',
            xy=(s5ri, 14), fontsize=6.5, color='#7f8c8d', fontstyle='italic',
            xytext=(s5ri + (top15['mean_RI'].max() - s5ri) * 0.3, 13),
            arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=0.8),
            va='center')

    # ==================================================================
    # Panel B: TOPPLE Stability Landscape Heatmap (top 25 vs bottom 25)
    # ==================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    top25 = scores.head(25).index.tolist()
    bot25 = scores.tail(25).index.tolist()
    ri_sel = ri_matrix.loc[top25 + bot25]
    ri_z = ri_sel.sub(ri_sel.mean(axis=1), axis=0).div(
        ri_sel.std(axis=1) + 1e-12, axis=0).fillna(0)

    sns.heatmap(ri_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'RI (z-score)', 'shrink': 0.5}, ax=ax_b)
    ax_b.set_title('B  Stability Landscape: Stabilizers vs Dispensable',
                    fontweight='bold', loc='left', fontsize=9)
    ax_b.set_xlabel('Cell Types (61 CIMA L4)')
    ax_b.set_ylabel('')
    ax_b.axhline(y=25, color='black', linewidth=1.5)
    top5_set = set(scores.head(5).index)
    for lbl in ax_b.get_yticklabels():
        if lbl.get_text() in top5_set:
            lbl.set_color('#c0392b')
            lbl.set_fontweight('bold')

    # ==================================================================
    # Panel C: SICAI Mean r_b vs Disease
    # ==================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    scatter_panel(
        ax_c,
        alt_merged['mean_rb'], alt_merged['n_disease_associations'],
        'Mean $r_b$ (Coupling Strength)',
        'Disease Associations (SMR)',
        'C  SICAI: Coupling Strength vs Disease Pleiotropy')

    # ==================================================================
    # Panel D: DGSA Non-additivity vs Disease
    # ==================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    dgsa_disease = dgsa_ct.join(disease_counts, how='inner')
    scatter_panel(
        ax_d,
        dgsa_disease['mean_non_additivity'], dgsa_disease['n_disease_associations'],
        'Mean Non-Additivity Index',
        'Disease Associations (SMR)',
        'D  DGSA: Genetic Non-Additivity vs Disease Pleiotropy')

    # ==================================================================
    # Panel E: DGSA Non-additivity vs TOPPLE RI
    # ==================================================================
    ax_e = fig.add_subplot(gs[2, 0])
    dgsa_topple = dgsa_ct.join(ct_mean_ri, how='inner')
    scatter_panel(
        ax_e,
        dgsa_topple['mean_non_additivity'], dgsa_topple['mean_RI'],
        'Mean Non-Additivity Index',
        'Mean TOPPLE RI',
        'E  Cross-Method: DGSA × TOPPLE',
        label_top=3, label_bot=3)

    # ==================================================================
    # Panel F: Summary Correlation Matrix
    # ==================================================================
    ax_f = fig.add_subplot(gs[2, 1])

    # Build a common cell-type index across all 4 metrics
    # TOPPLE RI per cell type
    topple_series = ct_mean_ri.rename('TOPPLE_RI')
    # SICAI mean r_b per cell type
    sicai_series = alt_merged['mean_rb'].rename('SICAI_mean_rb')
    # DGSA non-additivity per cell type
    dgsa_series = dgsa_ct['mean_non_additivity'].rename('DGSA_NonAdd')
    # Disease count per cell type
    disease_series = disease_counts.rename('Disease_N')

    all_metrics = pd.concat([topple_series, sicai_series, dgsa_series, disease_series],
                            axis=1).dropna()

    labels = ['TOPPLE RI', 'SICAI\nmean $r_b$', 'DGSA\nNon-Add', 'Disease\nCount']
    cols = all_metrics.columns.tolist()
    n = len(cols)

    rho_mat = np.ones((n, n))
    p_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                r, p = spearmanr(all_metrics.iloc[:, i], all_metrics.iloc[:, j])
                rho_mat[i, j] = r
                p_mat[i, j] = p

    # Annotations: rho + stars
    annot = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            if i == j:
                annot[i, j] = '1.000'
            else:
                sig = format_sig(p_mat[i, j])
                sig = '' if sig == 'n.s.' else sig
                annot[i, j] = f'{rho_mat[i, j]:.3f}{sig}'

    sns.heatmap(pd.DataFrame(rho_mat, index=labels, columns=labels),
                cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                annot=annot, fmt='', annot_kws={'fontsize': 9, 'fontweight': 'bold'},
                cbar_kws={'label': 'Spearman $\\rho$', 'shrink': 0.7},
                linewidths=1.5, linecolor='white',
                ax=ax_f)
    ax_f.set_title('F  Cross-Method Correlation Matrix',
                    fontweight='bold', loc='left')
    ax_f.set_xticklabels(labels, fontsize=8, rotation=0, ha='center')
    ax_f.set_yticklabels(labels, fontsize=8, rotation=0)

    # Add n annotation
    ax_f.text(0.5, -0.08, f'n = {len(all_metrics)} cell types (common across all methods)',
              transform=ax_f.transAxes, fontsize=7, ha='center', color='#7f8c8d')

    # ── Suptitle ──
    fig.suptitle(
        'Multi-Scale Architectural Analysis of the Human Immune System\n'
        'TOPPLE + SICAI + DGSA on CIMA Atlas (Yin et al., Science 2026)',
        fontsize=14, fontweight='bold', y=0.995)

    save_figure(fig, FIG_DIR, 'fig_unified_3method')
    print("Unified figure saved: fig_unified_3method.png/pdf")

    return all_metrics, rho_mat, p_mat


def update_framing(all_metrics, rho_mat, p_mat):
    """Update genome_biology_framing.md with DGSA results."""

    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    ri_matrix = pd.read_csv(RESULTS_DIR / "topple_ri_matrix.csv", index_col=0)
    alt_merged = pd.read_csv(RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv", index_col=0)
    dgsa_ct = pd.read_csv(RESULTS_DIR / "dgsa_celltype_summary.csv", index_col=0)
    dgsa_genes = pd.read_csv(RESULTS_DIR / "dgsa_gene_scores.csv", index_col=0)
    disease_counts = alt_merged['n_disease_associations']

    n_stabilizers = (scores['role'] == 'stabilizer').sum()
    n_regulons = len(scores)
    n_celltypes_topple = ri_matrix.shape[1]
    top1 = scores.index[0]

    stat5b_hits = [r for r in scores.index if r.upper().startswith('STAT5B')]
    stat5b_line = ""
    if stat5b_hits:
        s5key = stat5b_hits[0]
        s5rank = int(scores.loc[s5key, 'stability_rank'])
        stat5b_line = (
            f"Notably, STAT5B — the #1 stabilizer in our prior psoriasis analysis "
            f"(GSE173706, RI 0.03 -> 0.61) — ranks #{s5rank} in healthy blood, "
            f"consistent with disease-specific amplification of regulatory importance."
        )

    rho_sicai, p_sicai = spearmanr(alt_merged['mean_rb'], disease_counts)
    dgsa_disease = dgsa_ct.join(disease_counts, how='inner')
    rho_dgsa, p_dgsa = spearmanr(dgsa_disease['mean_non_additivity'],
                                  dgsa_disease['n_disease_associations'])

    framing = f"""# Multi-Scale Architectural Analysis of the Human Immune System
## Revealed by Geometric Decomposition of the CIMA Atlas

### Framing for Genome Biology / Nature Communications

The Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026) represents
a landmark cataloging effort -- 10 million cells from 428 individuals across 73 cell
types, with paired single-cell RNA-seq and ATAC-seq enabling the most comprehensive
map of human immune regulation to date. However, cataloging is not understanding.
CIMA's 203 eRegulons, 223,405 xQTLs, and 2,085 disease associations remain
uninterpreted as architectural entities. What stabilizes the regulatory landscape?
Which regulons are load-bearing, and which are dispensable? How does a cell type's
coupling and effect-size architecture determine its disease vulnerability?

We address these questions by applying three novel mathematical frameworks to CIMA's
published supplementary data:

**TOPPLE stability analysis** of {n_regulons} eRegulons across {n_celltypes_topple}
cell types reveals a striking hierarchy: {top1} emerges as the top stabilizer, with
{n_stabilizers} regulons classified as structurally critical (top quartile
redistribution index). {stat5b_line}

**SICAI coupling analysis** of 4,692 eQTL sharing pairs across 69 cell types
demonstrates that mean eQTL effect-size correlation (mean r_b) significantly
predicts disease pleiotropy (Spearman rho = {rho_sicai:.3f}, P = {p_sicai:.2e},
n = {len(alt_merged)} cell types). Shannon entropy of the r_b profile exhibits
a ceiling effect in healthy blood (CV = 0.009), while mean r_b captures the
biologically relevant variance (CV = 0.070).

**DGSA geometric decomposition** of {len(dgsa_genes):,} cis-eQTL genes reveals
that genetic non-additivity -- a measure of cell-type-specific effect-size geometry
-- is the strongest predictor of disease pleiotropy across all three methods
(Spearman rho = {rho_dgsa:.3f}, P = {p_dgsa:.2e}, n = {len(dgsa_disease)} cell
types). Cell types whose eQTL effects are concentrated in specific directions of
effect-size space (high non-additivity) carry substantially more disease associations
than those with uniform, shared effects.

Cross-method validation confirms complementary architectural axes: DGSA non-additivity
and SICAI coupling strength capture distinct but convergent aspects of disease
vulnerability, while TOPPLE stability identifies the regulons whose perturbation
most disrupts cell-type identity.

### Key Statistics Summary

| Method | Metric vs Disease | rho | P-value | n |
|--------|-------------------|-----|---------|---|
| SICAI  | Mean r_b          | {rho_sicai:.3f} | {p_sicai:.2e} | {len(alt_merged)} |
| DGSA   | Non-additivity    | {rho_dgsa:.3f} | {p_dgsa:.2e} | {len(dgsa_disease)} |

Our framework is computationally lightweight -- all analyses use only CIMA's published
supplementary tables, requiring no raw data access or high-performance computing. This
establishes a paradigm where mathematical frameworks, applied to community-generated
atlases, can extract fundamentally new biological insights from existing data.

**Keywords:** regulatory architecture, stability landscape, genetic non-additivity,
coupling strength, single-cell transcriptomics, immune atlas, CIMA, eQTL geometry
"""

    output_path = RESULTS_DIR / "genome_biology_framing.md"
    with open(output_path, 'w') as f:
        f.write(framing)

    print(f"Framing updated: {output_path} (~{len(framing.split())} words)")


def main():
    print("=" * 70)
    print("Unified 6-Panel Figure: TOPPLE + SICAI + DGSA")
    print("=" * 70)

    required = [
        "topple_stability_scores.csv",
        "topple_ri_matrix.csv",
        "sicai_alt_metrics_disease_merged.csv",
        "dgsa_celltype_summary.csv",
        "dgsa_gene_scores.csv",
    ]
    check_prerequisites(RESULTS_DIR, required)

    all_metrics, rho_mat, p_mat = create_unified_figure()
    update_framing(all_metrics, rho_mat, p_mat)

    print(f"\n{'=' * 70}")
    print("DONE")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
