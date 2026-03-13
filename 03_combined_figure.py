#!/usr/bin/env python3
"""
Day 3: Combined Proof-of-Concept Figure + Genome Biology Framing
================================================================
Combines TOPPLE (Day 1) + SICAI (Day 2) results into a 4-panel figure
suitable for Genome Biology / Nature Communications submission.

Panels:
  A: TOPPLE top 15 stabilizer bar chart (HSF1 #1, STAT5B #18 annotated)
  B: TOPPLE stability landscape heatmap (top 25 vs bottom 25 regulons)
  C: r_b heatmap (69×69, hierarchically clustered)
  D: Mean r_b vs disease associations scatter (rho=0.353, P=0.003, n=68)

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

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
                   check_prerequisites, COLORS, RESULTS_DIR, FIG_DIR)
setup_stdout()


def create_combined_figure():
    """Create 4-panel proof-of-concept figure for Genome Biology."""

    setup_figure_style()

    # Load results
    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    ri_matrix = pd.read_csv(RESULTS_DIR / "topple_ri_matrix.csv", index_col=0)
    rb_matrix = pd.read_csv(RESULTS_DIR / "sicai_rb_matrix.csv", index_col=0)
    alt_merged = pd.read_csv(RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv", index_col=0)

    # ── Figure layout: 2 rows × 2 columns ──
    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.1], hspace=0.35, wspace=0.3)

    # ══════════════════════════════════════════════════════════════════════
    # Panel A: TOPPLE Top 15 Stabilizer Bar Chart
    # ══════════════════════════════════════════════════════════════════════
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

    # Bold HSF1 (#1)
    ax_a.get_yticklabels()[0].set_fontweight('bold')
    ax_a.get_yticklabels()[0].set_color('#b71c1c')

    # Annotate STAT5B rank (#18, not in top 15 — add text annotation)
    stat5b_matches = [r for r in scores.index if r.upper().startswith('STAT5B')]
    if stat5b_matches:
        stat5b_key = stat5b_matches[0]
        stat5b_rank = int(scores.loc[stat5b_key, 'stability_rank'])
        stat5b_ri = scores.loc[stat5b_key, 'mean_RI']
        # Add a dashed reference line and annotation
        ax_a.axvline(stat5b_ri, color='#7f8c8d', linestyle=':', linewidth=1, alpha=0.7)
        ax_a.annotate(
            f'{stat5b_key} (#{stat5b_rank})\n#1 in psoriasis',
            xy=(stat5b_ri, 14),
            xytext=(stat5b_ri + (top15['mean_RI'].max() - stat5b_ri) * 0.3, 13),
            fontsize=6.5, color='#7f8c8d', fontstyle='italic',
            arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=0.8),
            va='center'
        )

    # ══════════════════════════════════════════════════════════════════════
    # Panel B: TOPPLE Stability Landscape Heatmap (top 25 vs bottom 25)
    # ══════════════════════════════════════════════════════════════════════
    ax_b = fig.add_subplot(gs[0, 1])
    top25 = scores.head(25).index.tolist()
    bot25 = scores.tail(25).index.tolist()
    selected_regs = top25 + bot25

    ri_subset = ri_matrix.loc[selected_regs]
    # Row-wise z-score for visualization
    ri_z = ri_subset.sub(ri_subset.mean(axis=1), axis=0).div(
        ri_subset.std(axis=1) + 1e-12, axis=0)
    ri_z = ri_z.fillna(0)

    sns.heatmap(
        ri_z,
        cmap='RdBu_r',
        center=0,
        vmin=-2, vmax=2,
        xticklabels=True, yticklabels=True,
        cbar_kws={'label': 'RI (z-score)', 'shrink': 0.5},
        ax=ax_b
    )
    ax_b.set_title('B  Stability Landscape: Stabilizers (top) vs Dispensable (bottom)',
                    fontweight='bold', loc='left', fontsize=9)
    ax_b.set_xlabel('Cell Types (61 CIMA L4)')
    ax_b.set_ylabel('')

    # Highlight top 5 y-labels
    top5_names = set(scores.head(5).index)
    for label in ax_b.get_yticklabels():
        if label.get_text() in top5_names:
            label.set_color(COLORS['stabilizer'])
            label.set_fontweight('bold')

    # Add divider line between top 25 and bottom 25
    ax_b.axhline(y=25, color='black', linewidth=1.5, linestyle='-')

    # ══════════════════════════════════════════════════════════════════════
    # Panel C: r_b Heatmap (69×69, hierarchically clustered)
    # ══════════════════════════════════════════════════════════════════════
    ax_c = fig.add_subplot(gs[1, 0])

    # Hierarchical clustering
    # Ensure symmetric, fill any NaN
    rb_vals = rb_matrix.values.copy()
    np.fill_diagonal(rb_vals, 1.0)
    rb_clean = pd.DataFrame(rb_vals, index=rb_matrix.index, columns=rb_matrix.columns)

    # Distance = 1 - r_b; clip to avoid negative distances
    dist = np.clip(1 - rb_clean.values, 0, 2)
    np.fill_diagonal(dist, 0)
    # Use condensed form
    dist_condensed = squareform(dist, checks=False)
    Z = linkage(dist_condensed, method='ward')
    order = leaves_list(Z)
    rb_ordered = rb_clean.iloc[order, order]

    sns.heatmap(
        rb_ordered,
        cmap='RdYlBu_r',
        vmin=0, vmax=1,
        xticklabels=True, yticklabels=True,
        cbar_kws={'label': '$r_b$ (effect-size correlation)', 'shrink': 0.5},
        ax=ax_c
    )
    ax_c.set_title('C  eQTL Effect-Size Correlation ($r_b$) — 69 Cell Types',
                    fontweight='bold', loc='left', fontsize=9)
    ax_c.set_xlabel('')
    ax_c.set_ylabel('')

    # ══════════════════════════════════════════════════════════════════════
    # Panel D: Mean r_b vs Disease Associations (the money plot)
    # ══════════════════════════════════════════════════════════════════════
    ax_d = fig.add_subplot(gs[1, 1])

    x = alt_merged['mean_rb']
    y = alt_merged['n_disease_associations']
    scatter_panel(ax_d, x, y,
                  xlabel='Mean $r_b$ (Coupling Strength)',
                  ylabel='Number of Disease Associations (SMR)',
                  title='D  Coupling Strength Predicts Disease Pleiotropy',
                  label_top=5, label_bot=3)

    # ── Suptitle ──
    fig.suptitle(
        'Regulatory Stability Architecture of the Human Immune System\n'
        'TOPPLE + SICAI Analysis of CIMA Atlas (Yin et al., Science 2026)',
        fontsize=13, fontweight='bold', y=1.02
    )

    save_figure(fig, FIG_DIR, 'fig_combined_proof_of_concept')
    print("✓ Combined figure saved: fig_combined_proof_of_concept.png/pdf")


def write_genome_biology_framing():
    """Generate framing text for Genome Biology cover letter / intro."""

    # Load results for stats
    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    ri_matrix = pd.read_csv(RESULTS_DIR / "topple_ri_matrix.csv", index_col=0)
    alt_merged = pd.read_csv(RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv", index_col=0)

    n_stabilizers = (scores['role'] == 'stabilizer').sum()
    n_regulons = len(scores)
    n_celltypes_topple = ri_matrix.shape[1]
    top1 = scores.index[0]
    n_celltypes_sicai = len(alt_merged)

    # STAT5B info
    stat5b_matches = [r for r in scores.index if r.upper().startswith('STAT5B')]
    stat5b_line = ""
    if stat5b_matches:
        stat5b_key = stat5b_matches[0]
        stat5b_rank = int(scores.loc[stat5b_key, 'stability_rank'])
        stat5b_line = (
            f"Notably, STAT5B — the #1 stabilizer in our prior psoriasis analysis "
            f"(GSE173706, RI 0.03→0.61) — ranks #{stat5b_rank} in healthy blood, "
            f"consistent with disease-specific amplification of regulatory importance."
        )

    rho, p_val = spearmanr(alt_merged['mean_rb'], alt_merged['n_disease_associations'])

    framing = f"""# Regulatory Stability Architecture of the Human Immune System
## Revealed by Multi-Scale Geometric Analysis of the CIMA Atlas

### Framing for Genome Biology

The Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026) represents
a landmark cataloging effort — 10 million cells from 428 individuals across 73 cell
types, with paired single-cell RNA-seq and ATAC-seq enabling the most comprehensive
map of human immune regulation to date. However, cataloging is not understanding.
CIMA's 203 eRegulons, 223,405 xQTLs, and 2,085 disease associations remain
uninterpreted as architectural entities. What stabilizes the regulatory landscape?
Which regulons are load-bearing, and which are dispensable? How does a cell type's
coupling architecture determine its disease vulnerability?

We address these questions by applying two novel mathematical frameworks — TOPPLE
(Transcriptomic Perturbation-based Phenotype Landscape Explorer) and SICAI
(Statistical Inference of Coupling Architecture Index) — to CIMA's published
supplementary data. These frameworks treat regulatory architecture itself as the
biomarker, shifting focus from individual regulon activities to the geometric
properties of the regulatory landscape.

**TOPPLE stability analysis** of {n_regulons} eRegulons across {n_celltypes_topple}
cell types reveals a striking hierarchy: {top1} emerges as the top stabilizer,
with {n_stabilizers} regulons classified as structurally critical (top quartile
redistribution index). {stat5b_line}

**SICAI coupling analysis** of 4,692 eQTL sharing pairs across 69 cell types
reveals that mean eQTL effect-size correlation (mean r_b) — a measure of coupling
strength — significantly predicts disease pleiotropy (Spearman rho = {rho:.3f},
P = {p_val:.2e}, n = {n_celltypes_sicai} cell types). Importantly, Shannon entropy
of the r_b profile shows a ceiling effect in healthy blood (CV = 0.009), while
mean r_b captures meaningful variance (CV = 0.070). This establishes that cell
types whose eQTL effects are broadly shared across the immune system carry greater
disease burden, connecting genetic architecture to clinical phenotype.

These results demonstrate that CIMA's catalog, when analyzed through the lens of
regulatory architecture, contains rich structural information invisible to standard
differential expression and enrichment approaches. The architecture IS the biomarker:
stability landscapes identify therapeutic targets (stabilizers as intervention nodes),
and coupling strength predicts disease susceptibility without requiring GWAS-level
sample sizes per cell type.

Our framework is computationally lightweight — all analyses use only CIMA's published
supplementary tables, requiring no raw data access or high-performance computing.
This establishes a paradigm where mathematical frameworks, applied to community-
generated atlases, can extract fundamentally new biological insights from existing
data.

**Keywords:** regulatory architecture, stability landscape, coupling strength,
single-cell transcriptomics, immune atlas, CIMA, perturbation analysis
"""

    output_path = RESULTS_DIR / "genome_biology_framing.md"
    with open(output_path, 'w') as f:
        f.write(framing)

    print(f"✓ Framing saved: {output_path}")
    print(f"  Word count: ~{len(framing.split())}")


def main():
    print("=" * 70)
    print("Day 3: Combined Figure + Genome Biology Framing")
    print("=" * 70)

    # Check prerequisites
    required = [
        "topple_stability_scores.csv",
        "topple_ri_matrix.csv",
        "sicai_rb_matrix.csv",
        "sicai_alt_metrics_disease_merged.csv",
    ]
    check_prerequisites(RESULTS_DIR, required)

    create_combined_figure()
    write_genome_biology_framing()

    print(f"\n{'=' * 70}")
    print("PROOF-OF-CONCEPT FIGURE COMPLETE")
    print(f"{'=' * 70}")
    print(f"  Figures: {list(FIG_DIR.glob('fig_combined*'))}")
    print(f"  Framing: {RESULTS_DIR / 'genome_biology_framing.md'}")
    print(f"\n  Next steps:")
    print(f"  1. Review figures -- iterate styling if needed")
    print(f"  2. Expand framing into full manuscript introduction")
    print(f"  3. Add DGSA analysis (Table S6) as third validation layer")
    print(f"  4. Submit to Genome Biology")


if __name__ == "__main__":
    main()
