#!/usr/bin/env python3
"""
SICAI Coupling Complexity Analysis on CIMA Tables S8 + S15
==========================================================
Applies SICAI (Statistical Inference of Coupling Architecture Index)
to 4,692 eQTL sharing pairs (r_b) across 69 cell types.

Core algorithm:
  1. Build 69×69 r_b matrix from Table S8
  2. Per cell type: coupling complexity = Shannon entropy of its r_b profile
     (high entropy → broadly coupled; low entropy → few strong partners)
  3. Load Table S15: count disease associations per cell type
  4. Test: does coupling complexity predict disease pleiotropy?
     (hypothesis: more coupled cell types → more disease associations)

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr, entropy
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
import warnings
import sys
import io
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Configuration ──────────────────────────────────────────────────────
DATA_DIR = Path("data")
FIG_DIR = Path("figures")
RESULTS_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# ── Step 1: Load Table S8 (eQTL sharing) ──────────────────────────────
def load_table_s8(data_dir: Path) -> pd.DataFrame:
    """Load CIMA Table S8: 4,692 pairwise eQTL sharing (π₁, r_b)."""
    candidates = [
        "science.adt3130_table_s8.xlsx",
        "science_adt3130_table_s8.xlsx",
        "table_s8.xlsx",
        "Table_S8.xlsx",
        "science_adt3130_table_s8.csv",
        "table_s8.csv",
    ]
    
    for fname in candidates:
        fpath = data_dir / fname
        if fpath.exists():
            print(f"✓ Loading Table S8: {fpath}")
            if fpath.suffix == '.csv':
                df = pd.read_csv(fpath)
            else:
                df = pd.read_excel(fpath, sheet_name=0)
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns)}")
            return df
    
    raise FileNotFoundError(f"Table S8 not found in {data_dir}/")


def build_rb_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Build 69×69 r_b correlation matrix from pairwise data.
    
    Expected columns: cell_type_1, cell_type_2, r_b (or rb, r.b)
    4,692 pairs = 69 × 68 / 2 = 2,346 × 2 (if both directions) or C(69,2) = 2,346
    """
    cols = [c.lower().strip() for c in df.columns]
    col_map = {c: orig for c, orig in zip(cols, df.columns)}
    
    # Auto-detect columns
    ct1_col = ct2_col = rb_col = pi1_col = None
    
    for c in cols:
        if any(x in c for x in ['cell_type_1', 'celltype_1', 'ct1', 'reference_cell', 'reference cell']):
            ct1_col = col_map[c]
        elif any(x in c for x in ['cell_type_2', 'celltype_2', 'ct2', 'query_cell', 'query cell']):
            ct2_col = col_map[c]
        elif c in ['r_b', 'rb', 'r.b', 'r_b_mean']:
            rb_col = col_map[c]
        elif c in ['pi1', 'π1', 'pi_1', 'pi1_mean']:
            pi1_col = col_map[c]
    
    if any(x is None for x in [ct1_col, ct2_col, rb_col]):
        print("\n⚠ Column detection incomplete. Columns available:")
        for i, c in enumerate(df.columns):
            print(f"  [{i}] {c} — sample: {df[c].iloc[0]}")
        # Positional fallback
        ct1_col = ct1_col or df.columns[0]
        ct2_col = ct2_col or df.columns[1]
        rb_col = rb_col or df.columns[2]
    
    print(f"\n✓ Columns: ct1={ct1_col}, ct2={ct2_col}, r_b={rb_col}")
    
    # Get all unique cell types
    all_cts = sorted(set(df[ct1_col].unique()) | set(df[ct2_col].unique()))
    n = len(all_cts)
    print(f"✓ Cell types: {n} (expected 69)")
    
    # Build symmetric matrix
    rb_matrix = pd.DataFrame(np.eye(n), index=all_cts, columns=all_cts)
    
    for _, row in df.iterrows():
        ct1 = row[ct1_col]
        ct2 = row[ct2_col]
        rb = row[rb_col]
        if pd.notna(rb):
            rb_matrix.loc[ct1, ct2] = rb
            rb_matrix.loc[ct2, ct1] = rb
    
    print(f"✓ r_b matrix: {rb_matrix.shape}")
    print(f"  r_b range: [{rb_matrix.values[rb_matrix.values < 1].min():.3f}, "
          f"{rb_matrix.values[rb_matrix.values < 1].max():.3f}]")
    print(f"  Mean r_b: {rb_matrix.values[rb_matrix.values < 1].mean():.3f} (paper: 0.82)")
    
    return rb_matrix, pi1_col, df, ct1_col, ct2_col


# ── Step 2: Coupling Complexity ────────────────────────────────────────
def compute_coupling_complexity(rb_matrix: pd.DataFrame) -> pd.DataFrame:
    """Compute SICAI coupling complexity per cell type.
    
    Coupling complexity = Shannon entropy of the r_b profile.
    
    For cell type i:
      1. Extract r_b vector: [r_b(i,1), r_b(i,2), ..., r_b(i,n)] excluding self
      2. Normalize to probability distribution (shift to positive, L1 normalize)
      3. H(i) = -Σ p_j × log2(p_j)  [Shannon entropy]
      4. High H → broadly coupled to many cell types (complex coupling)
      5. Low H → concentrated coupling to few partners (simple coupling)
    
    Also compute:
      - Mean r_b per cell type (coupling strength)
      - Max r_b (strongest partner)
      - Coupling breadth = number of partners with r_b > 0.8
    """
    celltypes = rb_matrix.index.tolist()
    n = len(celltypes)
    
    results = []
    
    for ct in celltypes:
        # r_b profile (exclude self-comparison)
        profile = rb_matrix.loc[ct].drop(ct, errors='ignore').values.astype(float)
        profile = np.nan_to_num(profile, nan=0.0)
        
        # Shift to positive for entropy calculation
        profile_pos = profile - profile.min() + 1e-10
        
        # L1 normalize
        p = profile_pos / profile_pos.sum()
        
        # Shannon entropy
        H = entropy(p, base=2)
        
        # Additional metrics
        mean_rb = profile.mean()
        max_rb = profile.max()
        std_rb = profile.std()
        breadth = (profile > 0.8).sum()
        
        # Coupling concentration (Gini-like coefficient)
        sorted_p = np.sort(p)
        cum = np.cumsum(sorted_p)
        gini = 1 - 2 * np.trapezoid(cum, dx=1/len(cum))
        
        results.append({
            'cell_type': ct,
            'coupling_complexity': H,
            'mean_rb': mean_rb,
            'max_rb': max_rb,
            'std_rb': std_rb,
            'coupling_breadth': breadth,
            'concentration': gini,
        })
    
    cc_df = pd.DataFrame(results).set_index('cell_type')
    cc_df = cc_df.sort_values('coupling_complexity', ascending=False)
    
    print(f"\n✓ Coupling complexity computed for {len(cc_df)} cell types")
    print(f"  Complexity range: [{cc_df['coupling_complexity'].min():.3f}, "
          f"{cc_df['coupling_complexity'].max():.3f}]")
    print(f"\n  Top 5 most complex:")
    print(cc_df.head()[['coupling_complexity', 'mean_rb', 'coupling_breadth']].to_string())
    print(f"\n  Bottom 5 least complex:")
    print(cc_df.tail()[['coupling_complexity', 'mean_rb', 'coupling_breadth']].to_string())
    
    return cc_df


# ── Step 3: Load Table S15 (Disease associations) ─────────────────────
def load_table_s15(data_dir: Path) -> pd.DataFrame:
    """Load CIMA Table S15: 2,085 SMR pleiotropic associations."""
    candidates = [
        "science.adt3130_table_s15.xlsx",
        "science_adt3130_table_s15.xlsx",
        "table_s15.xlsx",
        "Table_S15.xlsx",
        "science_adt3130_table_s15.csv",
        "table_s15.csv",
    ]
    
    for fname in candidates:
        fpath = data_dir / fname
        if fpath.exists():
            print(f"\n✓ Loading Table S15: {fpath}")
            if fpath.suffix == '.csv':
                df = pd.read_csv(fpath)
            else:
                df = pd.read_excel(fpath, sheet_name=0)
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns)}")
            return df
    
    raise FileNotFoundError(f"Table S15 not found in {data_dir}/")


def count_disease_associations(s15_df: pd.DataFrame) -> pd.Series:
    """Count disease associations per cell type from Table S15.
    
    Returns Series: cell_type → number of unique trait associations
    """
    cols = [c.lower().strip() for c in s15_df.columns]
    col_map = {c: orig for c, orig in zip(cols, s15_df.columns)}
    
    # Find cell type column
    ct_col = None
    for c in cols:
        # Match 'celltype', 'cell_type', 'cell_type_l4', etc. but not partial matches like 'query_celltype'
        if c in ['celltype', 'cell_type', 'cell_type_l4', 'cell type', 'l4']:
            ct_col = col_map[c]
            break

    # Find trait column
    trait_col = None
    for c in cols:
        if c in ['trait', 'disease', 'phenotype', 'gwas_trait', 'trait_name']:
            trait_col = col_map[c]
            break
    
    if ct_col is None:
        print("  ⚠ Auto-detect failed, using positional fallback")
        ct_col = s15_df.columns[0]
    if trait_col is None:
        trait_col = s15_df.columns[1] if len(s15_df.columns) > 1 else ct_col
    
    print(f"  Using: cell_type={ct_col}, trait={trait_col}")
    
    # Count unique traits per cell type
    disease_counts = s15_df.groupby(ct_col)[trait_col].nunique()
    disease_counts.name = 'n_disease_associations'
    
    print(f"  Cell types with disease associations: {len(disease_counts)}")
    print(f"  Association range: [{disease_counts.min()}, {disease_counts.max()}]")
    
    return disease_counts


# ── Step 4: Correlation Test ───────────────────────────────────────────
def test_complexity_disease_correlation(cc_df, disease_counts):
    """Core test: does coupling complexity predict disease pleiotropy?
    
    Merge coupling complexity with disease association counts.
    Compute Spearman ρ and P-value.
    """
    # Align cell type names
    merged = cc_df.join(disease_counts, how='inner')
    
    if len(merged) == 0:
        # Try fuzzy matching
        print("\n⚠ No exact cell type name matches. Attempting fuzzy matching...")
        cc_cts = set(cc_df.index)
        disease_cts = set(disease_counts.index)
        print(f"  Coupling complexity cell types (sample): {list(cc_cts)[:5]}")
        print(f"  Disease table cell types (sample): {list(disease_cts)[:5]}")
        
        # Try stripping/normalizing
        cc_normalized = {ct.strip().replace('-', '_').lower(): ct for ct in cc_cts}
        disease_normalized = {ct.strip().replace('-', '_').lower(): ct for ct in disease_cts}
        
        common = set(cc_normalized.keys()) & set(disease_normalized.keys())
        if common:
            mapping = {cc_normalized[k]: disease_normalized[k] for k in common}
            disease_remapped = disease_counts.rename(index={v: k for k, v in mapping.items()})
            merged = cc_df.join(disease_remapped, how='inner')
            print(f"  Fuzzy match found {len(merged)} common cell types")
    
    if len(merged) < 5:
        print(f"\n❌ Only {len(merged)} overlapping cell types. Check naming conventions.")
        return None
    
    n = len(merged)
    rho, p_val = spearmanr(merged['coupling_complexity'], merged['n_disease_associations'])
    r_pearson, p_pearson = pearsonr(merged['coupling_complexity'], merged['n_disease_associations'])
    
    print(f"\n{'=' * 70}")
    print(f"COUPLING COMPLEXITY vs DISEASE PLEIOTROPY (n={n} cell types)")
    print(f"{'=' * 70}")
    print(f"  Spearman ρ = {rho:.3f}, P = {p_val:.2e}")
    print(f"  Pearson  r = {r_pearson:.3f}, P = {p_pearson:.2e}")
    print(f"\n  Comparison with psoriasis result:")
    print(f"  Psoriasis: ρ = 0.480, P = 0.020, n = 23")
    print(f"  CIMA:      ρ = {rho:.3f}, P = {p_val:.2e}, n = {n}")
    
    if p_val < 0.05:
        print(f"\n  ✅ SIGNIFICANT — coupling complexity predicts disease pleiotropy!")
    else:
        print(f"\n  ⚠ Not significant at α=0.05 — but direction matters")
    
    return merged, rho, p_val


# ── Step 5: Publication Figures ────────────────────────────────────────
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
    })


def plot_rb_heatmap(rb_matrix, fig_dir):
    """Figure 2a: r_b correlation heatmap across 69 cell types."""
    setup_figure_style()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Cluster for visual grouping
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    Z = linkage(1 - rb_matrix.values, method='ward')
    order = leaves_list(Z)
    rb_ordered = rb_matrix.iloc[order, order]
    
    sns.heatmap(
        rb_ordered,
        cmap='RdYlBu_r',
        vmin=0, vmax=1,
        xticklabels=True, yticklabels=True,
        cbar_kws={'label': 'r_b (genetic correlation)', 'shrink': 0.6},
        ax=ax
    )
    ax.set_title('eQTL Effect Size Correlation (r_b) — CIMA Atlas\n'
                 '69 cell types, hierarchically clustered', fontweight='bold')
    
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig2a_rb_heatmap.{ext}')
    plt.close()
    print(f"✓ Figure 2a saved")


def plot_complexity_vs_disease(merged, rho, p_val, fig_dir):
    """Figure 2b: Scatter — coupling complexity vs disease association count."""
    setup_figure_style()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    x = merged['coupling_complexity']
    y = merged['n_disease_associations']
    
    ax.scatter(x, y, s=40, alpha=0.7, edgecolors='white', linewidth=0.5,
              c='#2c3e50', zorder=3)
    
    # Regression line
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, p(x_line), '--', color='#e74c3c', linewidth=1.5, alpha=0.8)
    
    # Label outliers (top 5 by disease count)
    top5 = merged.nlargest(5, 'n_disease_associations')
    for idx, row in top5.iterrows():
        ax.annotate(idx, (row['coupling_complexity'], row['n_disease_associations']),
                   fontsize=6, ha='left', va='bottom',
                   xytext=(5, 5), textcoords='offset points')
    
    # Statistics box
    stats_text = (f'Spearman ρ = {rho:.3f}\n'
                  f'P = {p_val:.2e}\n'
                  f'n = {len(merged)} cell types')
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            fontsize=8, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))
    
    ax.set_xlabel('Coupling Complexity (Shannon entropy of r_b profile)')
    ax.set_ylabel('Number of Disease Associations (SMR)')
    ax.set_title('Coupling Complexity Predicts Disease Pleiotropy\n'
                 'SICAI Analysis — CIMA Atlas', fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig2b_complexity_vs_disease.{ext}')
    plt.close()
    print(f"✓ Figure 2b saved")


def plot_complexity_distribution(cc_df, fig_dir):
    """Figure 2c: Distribution of coupling complexity across cell types."""
    setup_figure_style()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax1.hist(cc_df['coupling_complexity'], bins=20, color='#3498db', 
             edgecolor='white', alpha=0.8)
    ax1.set_xlabel('Coupling Complexity')
    ax1.set_ylabel('Count')
    ax1.set_title('Distribution of Coupling Complexity', fontweight='bold')
    ax1.axvline(cc_df['coupling_complexity'].median(), color='red', 
                linestyle='--', label=f"Median: {cc_df['coupling_complexity'].median():.3f}")
    ax1.legend()
    
    # Complexity vs breadth
    ax2.scatter(cc_df['coupling_complexity'], cc_df['coupling_breadth'],
               s=40, alpha=0.7, c='#2c3e50', edgecolors='white')
    ax2.set_xlabel('Coupling Complexity')
    ax2.set_ylabel('Coupling Breadth (# partners with r_b > 0.8)')
    ax2.set_title('Complexity vs Breadth', fontweight='bold')
    
    # Label extremes
    for label_set, n in [(cc_df.nlargest(3, 'coupling_complexity'), 3),
                         (cc_df.nsmallest(3, 'coupling_complexity'), 3)]:
        for idx, row in label_set.iterrows():
            ax2.annotate(idx, (row['coupling_complexity'], row['coupling_breadth']),
                        fontsize=6, ha='left')
    
    plt.suptitle('SICAI Coupling Architecture — CIMA Atlas', fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig2c_complexity_distribution.{ext}')
    plt.close()
    print(f"✓ Figure 2c saved")


# ── Step 6: Alternative Coupling Metrics ──────────────────────────────
def compute_alternative_metrics(rb_matrix: pd.DataFrame) -> pd.DataFrame:
    """Compute alternative coupling metrics that capture more variance
    than Shannon entropy (which has a ceiling effect in healthy blood).

    Metrics:
      1. Coupling heterogeneity: CV (std/mean) of r_b profile
      2. Coupling specificity: # partners with r_b < 0.7
      3. Mean r_b (coupling strength)
      4. Range of r_b (max - min)
    """
    celltypes = rb_matrix.index.tolist()
    results = []

    for ct in celltypes:
        profile = rb_matrix.loc[ct].drop(ct, errors='ignore').values.astype(float)
        profile = np.nan_to_num(profile, nan=0.0)

        mean_rb = profile.mean()
        std_rb = profile.std()
        cv_rb = std_rb / (mean_rb + 1e-12)  # heterogeneity
        n_specific = (profile < 0.7).sum()   # specificity
        rb_range = profile.max() - profile.min()  # range

        results.append({
            'cell_type': ct,
            'heterogeneity_cv': cv_rb,
            'specificity_n_low': int(n_specific),
            'mean_rb': mean_rb,
            'rb_range': rb_range,
        })

    alt_df = pd.DataFrame(results).set_index('cell_type')
    return alt_df


def compare_metrics_vs_disease(alt_df: pd.DataFrame, disease_counts: pd.Series,
                               fig_dir: Path, results_dir: Path):
    """Test all 4 alternative metrics against disease association count.
    Print comparison table and generate scatter for the best metric.
    """
    merged = alt_df.join(disease_counts, how='inner')
    n = len(merged)

    if n < 5:
        print(f"\n  Only {n} overlapping cell types — cannot test.")
        return None

    metric_labels = {
        'heterogeneity_cv': 'Coupling Heterogeneity (CV)',
        'specificity_n_low': 'Coupling Specificity (# r_b < 0.7)',
        'mean_rb': 'Mean r_b (Coupling Strength)',
        'rb_range': 'r_b Range (max − min)',
    }

    print(f"\n{'=' * 70}")
    print(f"ALTERNATIVE COUPLING METRICS vs DISEASE PLEIOTROPY (n={n})")
    print(f"{'=' * 70}")
    print(f"{'Metric':<40} {'Spearman rho':>12} {'P-value':>12} {'n':>5}")
    print("-" * 70)

    best_rho_abs = 0
    best_metric = None
    best_rho = None
    best_p = None
    rows = []

    for metric, label in metric_labels.items():
        rho, p_val = spearmanr(merged[metric], merged['n_disease_associations'])
        rows.append({'metric': label, 'rho': rho, 'p_value': p_val, 'n': n})
        print(f"{label:<40} {rho:>12.3f} {p_val:>12.2e} {n:>5}")
        if abs(rho) > best_rho_abs:
            best_rho_abs = abs(rho)
            best_metric = metric
            best_rho = rho
            best_p = p_val

    sig_marker = "**" if best_p < 0.01 else "*" if best_p < 0.05 else "(n.s.)"
    print(f"\n  Best metric: {metric_labels[best_metric]}")
    print(f"  rho = {best_rho:.3f}, P = {best_p:.2e} {sig_marker}")
    print(f"\n  Comparison — Psoriasis: rho = 0.480, P = 0.020, n = 23")

    # Save comparison table
    comp_df = pd.DataFrame(rows)
    comp_df.to_csv(results_dir / "sicai_alternative_metrics.csv", index=False)
    print(f"\n  Saved: {results_dir / 'sicai_alternative_metrics.csv'}")

    # Save full merged data
    merged.to_csv(results_dir / "sicai_alt_metrics_disease_merged.csv")

    # Generate scatter for best metric
    setup_figure_style()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = merged[best_metric]
    y = merged['n_disease_associations']

    ax.scatter(x, y, s=40, alpha=0.7, edgecolors='white', linewidth=0.5,
              c='#2c3e50', zorder=3)

    # Regression line
    z = np.polyfit(x, y, 1)
    poly = np.poly1d(z)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, poly(x_line), '--', color='#e74c3c', linewidth=1.5, alpha=0.8)

    # Label top 5 by disease count
    top5 = merged.nlargest(5, 'n_disease_associations')
    for idx, row in top5.iterrows():
        ax.annotate(idx, (row[best_metric], row['n_disease_associations']),
                   fontsize=6, ha='left', va='bottom',
                   xytext=(5, 5), textcoords='offset points')

    stats_text = (f'Spearman rho = {best_rho:.3f}\n'
                  f'P = {best_p:.2e}\n'
                  f'n = {n} cell types')
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            fontsize=8, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))

    ax.set_xlabel(metric_labels[best_metric])
    ax.set_ylabel('Number of Disease Associations (SMR)')
    ax.set_title(f'{metric_labels[best_metric]} vs Disease Pleiotropy\n'
                 f'SICAI Alternative Metric — CIMA Atlas', fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig2d_best_alt_metric_vs_disease.{ext}')
    plt.close()
    print(f"  Figure 2d saved: fig2d_best_alt_metric_vs_disease.png/pdf")

    return merged, best_metric, best_rho, best_p


# ── Main ───────────────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("SICAI Coupling Complexity Analysis on CIMA Atlas (Science 2026)")
    print("4,692 eQTL sharing pairs × 69 cell types — Tables S8 + S15")
    print("=" * 70)
    
    # 1. Load Table S8
    s8_df = load_table_s8(DATA_DIR)
    
    # 2. Build r_b matrix
    rb_matrix, pi1_col, s8_raw, ct1_col, ct2_col = build_rb_matrix(s8_df)
    
    # 3. Compute coupling complexity
    cc_df = compute_coupling_complexity(rb_matrix)
    
    # 4. Load Table S15
    s15_df = load_table_s15(DATA_DIR)
    disease_counts = count_disease_associations(s15_df)
    
    # 5. Core test
    result = test_complexity_disease_correlation(cc_df, disease_counts)
    
    # 6. Figures
    print("\n⏳ Generating figures...")
    plot_rb_heatmap(rb_matrix, FIG_DIR)
    plot_complexity_distribution(cc_df, FIG_DIR)
    
    if result is not None:
        merged, rho, p_val = result
        plot_complexity_vs_disease(merged, rho, p_val, FIG_DIR)
        merged.to_csv(RESULTS_DIR / "sicai_complexity_disease_merged.csv")
    
    # 7. Alternative metrics analysis
    print("\n" + "=" * 70)
    print("ALTERNATIVE COUPLING METRICS (addressing entropy ceiling effect)")
    print("=" * 70)
    alt_df = compute_alternative_metrics(rb_matrix)

    print(f"\n  Metric variance comparison:")
    print(f"  {'Metric':<30} {'Min':>8} {'Max':>8} {'CV':>8}")
    print(f"  {'-'*56}")
    for col in ['heterogeneity_cv', 'specificity_n_low', 'mean_rb', 'rb_range']:
        vals = alt_df[col]
        cv = vals.std() / (vals.mean() + 1e-12)
        print(f"  {col:<30} {vals.min():>8.3f} {vals.max():>8.3f} {cv:>8.3f}")
    # Compare with entropy
    ent_vals = cc_df['coupling_complexity']
    ent_cv = ent_vals.std() / (ent_vals.mean() + 1e-12)
    print(f"  {'Shannon entropy':<30} {ent_vals.min():>8.3f} {ent_vals.max():>8.3f} {ent_cv:>8.3f}")

    alt_result = compare_metrics_vs_disease(alt_df, disease_counts, FIG_DIR, RESULTS_DIR)

    # 8. Save results
    cc_df.to_csv(RESULTS_DIR / "sicai_coupling_complexity.csv")
    rb_matrix.to_csv(RESULTS_DIR / "sicai_rb_matrix.csv")

    print(f"\n{'=' * 70}")
    print("DAY 2 COMPLETE — Ready for Day 3 (combined figures)")
    print(f"{'=' * 70}")

    return rb_matrix, cc_df, result


if __name__ == "__main__":
    rb_matrix, cc_df, result = main()
