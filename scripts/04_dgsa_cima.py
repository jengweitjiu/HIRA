#!/usr/bin/env python3
"""
DGSA Geometric Decomposition on CIMA Table S6
==============================================
Applies DGSA (Directional Geometric Decomposition of Genetic Signal Architecture)
to 71,530 cis-eQTLs across 69 cell types from Yin et al., Science 2026.

Core algorithm:
  1. Build gene × cell_type effect-size (beta) matrix from cis-eQTLs
  2. Per gene: decompose beta vector into magnitude, direction, uniformity
  3. Non-additivity index = 1 - uniformity (cell-type-specific regulation)
  4. Aggregate per cell type → correlate with disease & TOPPLE stability

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import io
import warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Configuration ──────────────────────────────────────────────────────
DATA_DIR = Path("data")
FIG_DIR = Path("figures")
RESULTS_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

COLORS = {
    'primary': '#2c3e50',
    'accent': '#e74c3c',
    'secondary': '#3498db',
    'tertiary': '#27ae60',
}


def setup_figure_style():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })


# ── Step 1: Load and filter eQTLs ─────────────────────────────────────
def load_eqtls(data_dir: Path) -> pd.DataFrame:
    """Load Table S6, filter to cis-eQTLs, keep only needed columns."""
    fpath = data_dir / "CIMA_Table_S6.csv"
    print(f"Loading {fpath} (filtering to cis-eQTL only)...")

    # Read only needed columns to save memory
    usecols = ['phenotype_id', 'celltype', 'slope', 'slope_se',
               'pval_nominal', 'variant_id', 'analysis']
    df = pd.read_csv(fpath, usecols=usecols)
    print(f"  Total rows: {len(df):,}")

    # Filter to cis-eQTL
    eqtl = df[df['analysis'] == 'cis-eQTL'].copy()
    print(f"  cis-eQTL rows: {len(eqtl):,}")
    print(f"  Unique genes: {eqtl['phenotype_id'].nunique():,}")
    print(f"  Unique cell types: {eqtl['celltype'].nunique()}")

    del df  # free memory
    return eqtl


# ── Step 2: Build effect-size matrix ──────────────────────────────────
def build_beta_matrix(eqtl: pd.DataFrame, min_celltypes: int = 3) -> pd.DataFrame:
    """Build gene × cell_type matrix of effect sizes (slope/beta).

    Only keep genes mapped in >= min_celltypes cell types.
    """
    print(f"\nBuilding beta matrix (min {min_celltypes} cell types per gene)...")

    # Pivot: gene × celltype, values = slope
    beta_matrix = eqtl.pivot_table(
        index='phenotype_id',
        columns='celltype',
        values='slope',
        aggfunc='first'  # one lead QTL per gene per cell type
    )

    n_genes_raw = len(beta_matrix)
    # Count non-NaN cell types per gene
    ct_counts = beta_matrix.notna().sum(axis=1)
    beta_matrix = beta_matrix[ct_counts >= min_celltypes]

    print(f"  Raw genes: {n_genes_raw:,}")
    print(f"  Genes with >= {min_celltypes} cell types: {len(beta_matrix):,}")
    print(f"  Matrix shape: {beta_matrix.shape}")
    print(f"  Sparsity: {beta_matrix.isna().sum().sum() / beta_matrix.size:.1%} NaN")

    return beta_matrix


# ── Step 3: Geometric decomposition per gene ──────────────────────────
def gini_coefficient(x):
    """Gini coefficient of array x (0 = perfect equality, 1 = max inequality)."""
    x = np.abs(x)
    x = np.sort(x)
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * x) / (n * np.sum(x))) - (n + 1) / n


def geometric_decomposition(beta_matrix: pd.DataFrame) -> pd.DataFrame:
    """Decompose each gene's beta vector into geometric components.

    For each gene:
      - magnitude: L2 norm of beta vector
      - direction: unit vector (implicit, not stored)
      - uniformity: cosine similarity to uniform vector
      - sparsity: fraction of cell types with non-zero (non-NaN) beta
      - gini: Gini coefficient of absolute betas
      - non_additivity: 1 - uniformity
    """
    n_genes = len(beta_matrix)
    n_ct = beta_matrix.shape[1]
    print(f"\nGeometric decomposition of {n_genes:,} genes across {n_ct} cell types...")

    # Uniform reference vector
    uniform = np.ones(n_ct) / np.sqrt(n_ct)

    results = []
    for gene in beta_matrix.index:
        betas = beta_matrix.loc[gene].values.copy()
        # Replace NaN with 0 for geometric computation
        betas_filled = np.nan_to_num(betas, nan=0.0)

        n_nonzero = np.count_nonzero(~np.isnan(beta_matrix.loc[gene].values))

        # Magnitude (L2 norm)
        magnitude = np.linalg.norm(betas_filled)

        if magnitude < 1e-15:
            results.append({
                'gene': gene,
                'magnitude': 0.0,
                'uniformity': 0.0,
                'sparsity': n_nonzero / n_ct,
                'gini': 0.0,
                'non_additivity': 1.0,
                'n_celltypes': n_nonzero,
                'mean_abs_beta': 0.0,
                'max_abs_beta': 0.0,
            })
            continue

        # Direction (unit vector)
        direction = betas_filled / magnitude

        # Uniformity: cosine similarity to uniform vector
        uniformity = np.dot(direction, uniform)
        # Clamp to [0, 1] — negative means net opposing directions
        uniformity = max(0.0, uniformity)

        # Sparsity
        sparsity = n_nonzero / n_ct

        # Gini of absolute betas (only non-NaN entries)
        nonnan_betas = betas_filled[~np.isnan(beta_matrix.loc[gene].values)]
        gini = gini_coefficient(nonnan_betas)

        # Non-additivity
        non_additivity = 1.0 - uniformity

        results.append({
            'gene': gene,
            'magnitude': magnitude,
            'uniformity': uniformity,
            'sparsity': sparsity,
            'gini': gini,
            'non_additivity': non_additivity,
            'n_celltypes': n_nonzero,
            'mean_abs_beta': np.abs(nonnan_betas).mean(),
            'max_abs_beta': np.abs(nonnan_betas).max(),
        })

    gene_scores = pd.DataFrame(results).set_index('gene')
    gene_scores = gene_scores.sort_values('non_additivity', ascending=False)

    print(f"  Non-additivity range: [{gene_scores['non_additivity'].min():.3f}, "
          f"{gene_scores['non_additivity'].max():.3f}]")
    print(f"  Mean non-additivity: {gene_scores['non_additivity'].mean():.3f}")
    print(f"  Median: {gene_scores['non_additivity'].median():.3f}")

    return gene_scores


# ── Step 4: Aggregate per cell type ───────────────────────────────────
def aggregate_per_celltype(beta_matrix: pd.DataFrame,
                           gene_scores: pd.DataFrame) -> pd.DataFrame:
    """Compute mean non-additivity for each cell type.

    For cell type j: average non_additivity of all genes with a QTL in j.
    """
    print("\nAggregating per cell type...")
    celltypes = beta_matrix.columns.tolist()
    results = []

    for ct in celltypes:
        # Genes with non-NaN beta in this cell type
        genes_in_ct = beta_matrix[ct].dropna().index
        genes_in_scores = genes_in_ct.intersection(gene_scores.index)

        if len(genes_in_scores) == 0:
            continue

        subset = gene_scores.loc[genes_in_scores]
        results.append({
            'cell_type': ct,
            'n_genes': len(genes_in_scores),
            'mean_non_additivity': subset['non_additivity'].mean(),
            'median_non_additivity': subset['non_additivity'].median(),
            'mean_magnitude': subset['magnitude'].mean(),
            'mean_gini': subset['gini'].mean(),
            'mean_uniformity': subset['uniformity'].mean(),
        })

    ct_summary = pd.DataFrame(results).set_index('cell_type')
    ct_summary = ct_summary.sort_values('mean_non_additivity', ascending=False)

    print(f"  Cell types: {len(ct_summary)}")
    print(f"  Mean non-additivity range: [{ct_summary['mean_non_additivity'].min():.4f}, "
          f"{ct_summary['mean_non_additivity'].max():.4f}]")

    return ct_summary


# ── Step 5: Cross-method correlations ─────────────────────────────────
def cross_method_correlations(ct_summary: pd.DataFrame,
                              results_dir: Path) -> dict:
    """Correlate DGSA cell-type metrics with TOPPLE stability and disease counts."""
    correlations = {}

    # Load disease associations from SICAI Day 2
    alt_merged_path = results_dir / "sicai_alt_metrics_disease_merged.csv"
    if alt_merged_path.exists():
        alt_merged = pd.read_csv(alt_merged_path, index_col=0)
        disease_counts = alt_merged['n_disease_associations']

        merged_d = ct_summary.join(disease_counts, how='inner')
        if len(merged_d) >= 5:
            rho, p = spearmanr(merged_d['mean_non_additivity'],
                               merged_d['n_disease_associations'])
            correlations['disease'] = {
                'rho': rho, 'p': p, 'n': len(merged_d), 'merged': merged_d
            }
            print(f"\n  DGSA non-additivity vs disease associations:")
            print(f"    Spearman rho = {rho:.3f}, P = {p:.2e}, n = {len(merged_d)}")

    # Load TOPPLE stability scores
    topple_path = results_dir / "topple_stability_scores.csv"
    if topple_path.exists():
        topple = pd.read_csv(topple_path, index_col=0)
        # TOPPLE is per-regulon, ct_summary is per-cell-type
        # We need per-cell-type mean RI from the RI matrix
        ri_path = results_dir / "topple_ri_matrix.csv"
        if ri_path.exists():
            ri_matrix = pd.read_csv(ri_path, index_col=0)
            # Mean RI per cell type (across all regulons)
            ct_mean_ri = ri_matrix.mean(axis=0)
            ct_mean_ri.name = 'mean_RI'

            merged_t = ct_summary.join(ct_mean_ri, how='inner')
            if len(merged_t) >= 5:
                rho, p = spearmanr(merged_t['mean_non_additivity'],
                                   merged_t['mean_RI'])
                correlations['topple'] = {
                    'rho': rho, 'p': p, 'n': len(merged_t), 'merged': merged_t
                }
                print(f"\n  DGSA non-additivity vs TOPPLE mean RI:")
                print(f"    Spearman rho = {rho:.3f}, P = {p:.2e}, n = {len(merged_t)}")

    return correlations


# ── Step 6: Figures ───────────────────────────────────────────────────
def find_top_gene_celltypes(beta_matrix, gene, top_n=3):
    """Find the top cell types by absolute beta for a gene."""
    row = beta_matrix.loc[gene].dropna().abs().sort_values(ascending=False)
    return list(row.head(top_n).index)


def plot_figures(gene_scores, ct_summary, beta_matrix, correlations, fig_dir):
    """Generate all DGSA figures."""
    setup_figure_style()

    # ── Fig 3a: Distribution of non-additivity ──
    fig, ax = plt.subplots(figsize=(7, 5))
    vals = gene_scores['non_additivity']
    ax.hist(vals, bins=60, color=COLORS['secondary'], edgecolor='white', alpha=0.85)
    ax.axvline(vals.median(), color=COLORS['accent'], linestyle='--', linewidth=1.5,
               label=f'Median = {vals.median():.3f}')
    ax.axvline(vals.mean(), color=COLORS['tertiary'], linestyle=':', linewidth=1.5,
               label=f'Mean = {vals.mean():.3f}')
    ax.set_xlabel('Non-Additivity Index (1 - uniformity)')
    ax.set_ylabel('Number of Genes')
    ax.set_title('DGSA: Distribution of Genetic Non-Additivity\n'
                 f'{len(gene_scores):,} cis-eQTL genes across 69 cell types',
                 fontweight='bold')
    ax.legend()
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig3a_non_additivity_distribution.{ext}')
    plt.close()
    print("  fig3a saved")

    # ── Fig 3b: Non-additivity vs disease associations ──
    if 'disease' in correlations:
        info = correlations['disease']
        merged = info['merged']
        fig, ax = plt.subplots(figsize=(8, 6))
        x = merged['mean_non_additivity']
        y = merged['n_disease_associations']

        ax.scatter(x, y, s=45, alpha=0.7, c=COLORS['primary'],
                   edgecolors='white', linewidth=0.5, zorder=3)

        # Regression line
        z = np.polyfit(x, y, 1)
        p_fit = np.poly1d(z)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, p_fit(x_line), '--', color=COLORS['accent'], linewidth=1.5)

        # Label top 5
        top5 = merged.nlargest(5, 'n_disease_associations')
        for idx, row in top5.iterrows():
            ax.annotate(idx, (row['mean_non_additivity'], row['n_disease_associations']),
                        fontsize=5.5, ha='left', va='bottom',
                        xytext=(4, 4), textcoords='offset points')

        stats = (f'Spearman $\\rho$ = {info["rho"]:.3f}\n'
                 f'P = {info["p"]:.2e}\nn = {info["n"]} cell types')
        ax.text(0.05, 0.95, stats, transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.8))

        ax.set_xlabel('Mean Non-Additivity Index per Cell Type')
        ax.set_ylabel('Number of Disease Associations (SMR)')
        ax.set_title('DGSA: Genetic Non-Additivity Predicts Disease Pleiotropy',
                      fontweight='bold')
        ax.grid(True, alpha=0.2)
        plt.tight_layout()
        for ext in ['png', 'pdf']:
            fig.savefig(fig_dir / f'fig3b_nonadditivity_vs_disease.{ext}')
        plt.close()
        print("  fig3b saved")

    # ── Fig 3c: Non-additivity vs TOPPLE stability ──
    if 'topple' in correlations:
        info = correlations['topple']
        merged = info['merged']
        fig, ax = plt.subplots(figsize=(8, 6))
        x = merged['mean_non_additivity']
        y = merged['mean_RI']

        ax.scatter(x, y, s=45, alpha=0.7, c=COLORS['primary'],
                   edgecolors='white', linewidth=0.5, zorder=3)

        z = np.polyfit(x, y, 1)
        p_fit = np.poly1d(z)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, p_fit(x_line), '--', color=COLORS['accent'], linewidth=1.5)

        # Label extremes
        for subset in [merged.nlargest(3, 'mean_RI'), merged.nsmallest(3, 'mean_RI')]:
            for idx, row in subset.iterrows():
                ax.annotate(idx, (row['mean_non_additivity'], row['mean_RI']),
                            fontsize=5.5, ha='left', va='bottom',
                            xytext=(4, 4), textcoords='offset points')

        stats = (f'Spearman $\\rho$ = {info["rho"]:.3f}\n'
                 f'P = {info["p"]:.2e}\nn = {info["n"]} cell types')
        ax.text(0.05, 0.95, stats, transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.8))

        ax.set_xlabel('Mean Non-Additivity Index per Cell Type')
        ax.set_ylabel('Mean TOPPLE Redistribution Index (RI)')
        ax.set_title('DGSA × TOPPLE: Non-Additivity vs Regulatory Stability',
                      fontweight='bold')
        ax.grid(True, alpha=0.2)
        plt.tight_layout()
        for ext in ['png', 'pdf']:
            fig.savefig(fig_dir / f'fig3c_nonadditivity_vs_topple.{ext}')
        plt.close()
        print("  fig3c saved")

    # ── Fig 3d: Top 20 most non-additive genes ──
    fig, ax = plt.subplots(figsize=(9, 6))
    top20 = gene_scores.head(20)
    colors_bar = [COLORS['accent']] * 20

    ax.barh(range(len(top20)), top20['non_additivity'], color=colors_bar,
            edgecolor='white', height=0.7)
    ax.set_yticks(range(len(top20)))

    # Annotate with top cell types
    ylabels = []
    for gene in top20.index:
        if gene in beta_matrix.index:
            top_cts = find_top_gene_celltypes(beta_matrix, gene, top_n=2)
            label = f"{gene}  ({', '.join(top_cts)})"
        else:
            label = gene
        ylabels.append(label)

    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel('Non-Additivity Index')
    ax.set_title('DGSA: Top 20 Most Cell-Type-Specific eQTL Genes\n'
                 '(top 2 cell types annotated)', fontweight='bold')
    ax.invert_yaxis()
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig3d_top20_nonadditivity.{ext}')
    plt.close()
    print("  fig3d saved")


# ── Main ───────────────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("DGSA Geometric Decomposition on CIMA Atlas (Science 2026)")
    print("71,530 cis-eQTLs across 69 cell types — Table S6")
    print("=" * 70)

    # 1. Load eQTLs
    eqtl = load_eqtls(DATA_DIR)

    # 2. Build beta matrix
    beta_matrix = build_beta_matrix(eqtl, min_celltypes=3)
    del eqtl  # free memory

    # 3. Geometric decomposition
    gene_scores = geometric_decomposition(beta_matrix)

    # 4. Print top results
    print(f"\n{'=' * 70}")
    print("TOP 10 MOST NON-ADDITIVE GENES (cell-type-specific eQTLs)")
    print("=" * 70)
    print(gene_scores.head(10)[['non_additivity', 'uniformity', 'magnitude',
                                 'gini', 'n_celltypes']].to_string())

    print(f"\n{'=' * 70}")
    print("TOP 10 MOST ADDITIVE/UNIFORM GENES (shared eQTLs)")
    print("=" * 70)
    print(gene_scores.tail(10)[['non_additivity', 'uniformity', 'magnitude',
                                 'gini', 'n_celltypes']].to_string())

    # 5. Aggregate per cell type
    ct_summary = aggregate_per_celltype(beta_matrix, gene_scores)

    print(f"\n  Top 5 cell types by mean non-additivity:")
    print(ct_summary.head()[['mean_non_additivity', 'n_genes', 'mean_gini']].to_string())
    print(f"\n  Bottom 5:")
    print(ct_summary.tail()[['mean_non_additivity', 'n_genes', 'mean_gini']].to_string())

    # 6. Cross-method correlations
    print(f"\n{'=' * 70}")
    print("CROSS-METHOD CORRELATIONS")
    print("=" * 70)
    correlations = cross_method_correlations(ct_summary, RESULTS_DIR)

    # 7. Figures
    print("\nGenerating figures...")
    plot_figures(gene_scores, ct_summary, beta_matrix, correlations, FIG_DIR)

    # 8. Save results
    gene_scores.to_csv(RESULTS_DIR / "dgsa_gene_scores.csv")
    ct_summary.to_csv(RESULTS_DIR / "dgsa_celltype_summary.csv")
    print(f"\nResults saved to {RESULTS_DIR}/")

    # 9. Summary
    print(f"\n{'=' * 70}")
    print("DGSA ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"  Genes analyzed: {len(gene_scores):,}")
    print(f"  Cell types: {len(ct_summary)}")
    print(f"  Mean non-additivity: {gene_scores['non_additivity'].mean():.3f}")
    if 'disease' in correlations:
        d = correlations['disease']
        print(f"  Non-additivity vs disease: rho={d['rho']:.3f}, P={d['p']:.2e}")
    if 'topple' in correlations:
        t = correlations['topple']
        print(f"  Non-additivity vs TOPPLE:  rho={t['rho']:.3f}, P={t['p']:.2e}")
    print(f"  Figures: fig3a-fig3d in {FIG_DIR}/")

    return gene_scores, ct_summary, correlations


if __name__ == "__main__":
    gene_scores, ct_summary, correlations = main()
