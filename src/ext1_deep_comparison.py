#!/usr/bin/env python3
"""
Extension #1 — Deep distributional comparison: eQTL vs caQTL non-additivity.

Three analyses:
  1. Overlaid histograms + QQ-plot
  2. Paired comparison at matched loci (same variant_id + celltype)
  3. Per-cell-type caQTL non-additivity vs disease count (S15)

Input:  results/dgsa_eqtl.csv, results/ext1_caqtl_dgsa.csv,
        data/raw/CIMA_Table_S6.csv, data/raw/science.adt3130_table_s15.xlsx
Output: figures/ext1_distribution_comparison.pdf
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


# ---------- helpers ----------------------------------------------------------

def compute_gini(arr):
    """Gini coefficient for a 1-D array of non-negative values."""
    arr = np.sort(arr.copy())
    n = len(arr)
    if n == 0 or arr.sum() == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * arr) - (n + 1) * np.sum(arr)) / (n * np.sum(arr))


def setup_style():
    """Nature-style plot settings."""
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'lines.linewidth': 1.0,
    })


# ---------- Analysis 1: Distribution plots -----------------------------------

def plot_distributions(eqtl_na, caqtl_na, outpath):
    """Overlaid histograms and QQ-plot."""
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(180/25.4, 80/25.4))

    # Panel A: overlaid histograms
    ax = axes[0]
    bins = np.linspace(0.3, 1.0, 60)
    ax.hist(eqtl_na, bins=bins, alpha=0.6, density=True,
            color='#4878CF', label=f'eQTL (n={len(eqtl_na):,})', edgecolor='none')
    ax.hist(caqtl_na, bins=bins, alpha=0.6, density=True,
            color='#D65F5F', label=f'caQTL (n={len(caqtl_na):,})', edgecolor='none')
    ax.set_xlabel('Non-additivity (Gini of slope²)')
    ax.set_ylabel('Density')
    ax.set_title('A', fontsize=10, fontweight='bold', loc='left')
    ax.legend(frameon=False, fontsize=7)

    # Add vertical lines for means and medians
    ax.axvline(np.mean(eqtl_na), color='#4878CF', linestyle='--', linewidth=0.8)
    ax.axvline(np.mean(caqtl_na), color='#D65F5F', linestyle='--', linewidth=0.8)
    ax.axvline(np.median(eqtl_na), color='#4878CF', linestyle=':', linewidth=0.8)
    ax.axvline(np.median(caqtl_na), color='#D65F5F', linestyle=':', linewidth=0.8)

    # Panel B: QQ-plot
    ax = axes[1]
    # Quantile-quantile: caQTL quantiles vs eQTL quantiles
    n_q = 500
    probs = np.linspace(0.5, 99.5, n_q)
    eq_q = np.percentile(eqtl_na, probs)
    cq_q = np.percentile(caqtl_na, probs)
    ax.scatter(eq_q, cq_q, s=4, alpha=0.6, color='#333333', zorder=3)
    lims = [min(eq_q.min(), cq_q.min()) - 0.01, max(eq_q.max(), cq_q.max()) + 0.01]
    ax.plot(lims, lims, 'k--', linewidth=0.6, alpha=0.5, zorder=2)
    ax.set_xlabel('eQTL non-additivity quantiles')
    ax.set_ylabel('caQTL non-additivity quantiles')
    ax.set_title('B', fontsize=10, fontweight='bold', loc='left')
    ax.set_aspect('equal')

    # Annotate crossover
    crossover_idx = np.where(np.diff(np.sign(cq_q - eq_q)))[0]
    if len(crossover_idx) > 0:
        cx = crossover_idx[0]
        ax.annotate(f'crossover ~{eq_q[cx]:.2f}',
                    xy=(eq_q[cx], cq_q[cx]), fontsize=6,
                    xytext=(eq_q[cx] - 0.08, cq_q[cx] + 0.04),
                    arrowprops=dict(arrowstyle='->', lw=0.5, color='red'),
                    color='red')

    plt.tight_layout()
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    fig.savefig(str(outpath).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")
    plt.close(fig)


# ---------- Analysis 2: Matched loci paired comparison -----------------------

def matched_loci_analysis(s6_path, eqtl_pivot_celltypes, caqtl_pivot_celltypes):
    """Paired comparison of non-additivity for matched gene-peak pairs."""
    print("\n" + "=" * 60)
    print("ANALYSIS 2: MATCHED LOCI (same variant_id + celltype)")
    print("=" * 60)

    df = pd.read_csv(s6_path, usecols=[
        'phenotype_id', 'celltype', 'slope', 'analysis', 'variant_id'
    ])
    eqtl = df[df['analysis'] == 'cis-eQTL']
    caqtl = df[df['analysis'] == 'cis-caQTL']

    # Match on variant_id + celltype
    matched = eqtl.merge(caqtl, on=['variant_id', 'celltype'],
                         suffixes=('_eqtl', '_caqtl'))
    print(f"  Matched rows (variant+celltype): {len(matched):,}")
    print(f"  Unique genes: {matched['phenotype_id_eqtl'].nunique():,}")
    print(f"  Unique peaks: {matched['phenotype_id_caqtl'].nunique():,}")

    # For each matched gene: build its slope vector across all eQTL cell types
    # For each matched peak: build its slope vector across all caQTL cell types
    # Then pair them and compare non-additivity

    # Get unique gene-peak pairs (a gene can map to multiple peaks and vice versa)
    # Group by (gene, peak): collect slopes across cell types
    pairs = matched.groupby(['phenotype_id_eqtl', 'phenotype_id_caqtl']).agg(
        n_shared_celltypes=('celltype', 'count'),
    ).reset_index()
    pairs = pairs[pairs['n_shared_celltypes'] >= 3]
    print(f"  Gene-peak pairs with >= 3 shared cell types: {len(pairs):,}")

    if len(pairs) == 0:
        print("  Not enough matched pairs with >= 3 cell types.")
        # Fall back: compute per-gene and per-peak from their FULL profiles
        # then match on shared variant_id
        return _matched_full_profiles(matched, eqtl, caqtl,
                                      eqtl_pivot_celltypes, caqtl_pivot_celltypes)

    return _matched_full_profiles(matched, eqtl, caqtl,
                                  eqtl_pivot_celltypes, caqtl_pivot_celltypes)


def _matched_full_profiles(matched, eqtl_raw, caqtl_raw,
                           eqtl_celltypes, caqtl_celltypes):
    """Compute non-additivity from full cell-type profiles for matched genes/peaks."""
    # Get unique gene-peak pairs
    gene_peak_pairs = matched[['phenotype_id_eqtl', 'phenotype_id_caqtl']].drop_duplicates()
    print(f"  Unique gene-peak pairs: {len(gene_peak_pairs):,}")

    # Build full-space non-additivity for each gene and peak
    # eQTL: gene → slope across 69 cell types
    eqtl_pivot = eqtl_raw.pivot_table(
        index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
    ).reindex(columns=eqtl_celltypes, fill_value=np.nan).fillna(0)

    caqtl_pivot = caqtl_raw.pivot_table(
        index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
    ).reindex(columns=caqtl_celltypes, fill_value=np.nan).fillna(0)

    gene_na = {}
    for gene in eqtl_pivot.index:
        s = eqtl_pivot.loc[gene].values
        gene_na[gene] = compute_gini(s ** 2)

    peak_na = {}
    for peak in caqtl_pivot.index:
        s = caqtl_pivot.loc[peak].values
        peak_na[peak] = compute_gini(s ** 2)

    # Pair up
    paired_eqtl = []
    paired_caqtl = []
    for _, row in gene_peak_pairs.iterrows():
        gene = row['phenotype_id_eqtl']
        peak = row['phenotype_id_caqtl']
        if gene in gene_na and peak in peak_na:
            paired_eqtl.append(gene_na[gene])
            paired_caqtl.append(peak_na[peak])

    paired_eqtl = np.array(paired_eqtl)
    paired_caqtl = np.array(paired_caqtl)
    print(f"  Paired comparisons: {len(paired_eqtl):,}")

    # Wilcoxon signed-rank (paired)
    diff = paired_caqtl - paired_eqtl
    w_stat, w_p = stats.wilcoxon(diff, alternative='two-sided')
    _, w_p_greater = stats.wilcoxon(diff, alternative='greater')
    _, w_p_less = stats.wilcoxon(diff, alternative='less')

    print(f"\n  Paired statistics:")
    print(f"    eQTL  mean={np.mean(paired_eqtl):.4f}  median={np.median(paired_eqtl):.4f}")
    print(f"    caQTL mean={np.mean(paired_caqtl):.4f}  median={np.median(paired_caqtl):.4f}")
    print(f"    Mean diff (caQTL - eQTL): {np.mean(diff):.4f}")
    print(f"    Median diff:              {np.median(diff):.4f}")
    print(f"    Wilcoxon W:          {w_stat:.0f}")
    print(f"    P (two-sided):       {w_p:.2e}")
    print(f"    P (caQTL > eQTL):    {w_p_greater:.2e}")
    print(f"    P (caQTL < eQTL):    {w_p_less:.2e}")
    print(f"    Pairs with caQTL > eQTL: {np.sum(diff > 0):,} "
          f"({100*np.mean(diff > 0):.1f}%)")
    print(f"    Pairs with caQTL < eQTL: {np.sum(diff < 0):,} "
          f"({100*np.mean(diff < 0):.1f}%)")
    print(f"    Pairs with caQTL = eQTL: {np.sum(diff == 0):,}")

    # Spearman correlation between paired values
    rho, p_rho = stats.spearmanr(paired_eqtl, paired_caqtl)
    print(f"\n    Spearman (paired eQTL vs caQTL): rho={rho:.4f}, P={p_rho:.2e}")

    return paired_eqtl, paired_caqtl


# ---------- Analysis 3: Per-cell-type caQTL vs disease -----------------------

def celltype_disease_analysis(caqtl_df, caqtl_pivot_path, eqtl_df, eqtl_pivot_path,
                              s6_path, s15_path):
    """Per-cell-type mean caQTL non-additivity vs disease count."""
    print("\n" + "=" * 60)
    print("ANALYSIS 3: PER-CELL-TYPE caQTL vs DISEASE")
    print("=" * 60)

    s15 = pd.read_excel(s15_path, sheet_name='Sheet1')
    disease_ct = s15.groupby('celltype')['trait'].nunique().reset_index()
    disease_ct.columns = ['celltype', 'disease_count']

    # Rebuild pivots from raw data to get cell-type membership
    df = pd.read_csv(s6_path, usecols=['phenotype_id', 'celltype', 'slope', 'analysis'])

    results_both = {}
    for label, analysis_val, dgsa_df, id_col in [
        ('eQTL', 'cis-eQTL', eqtl_df, 'gene'),
        ('caQTL', 'cis-caQTL', caqtl_df, 'peak'),
    ]:
        sub = df[df['analysis'] == analysis_val]
        pivot = sub.pivot_table(
            index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
        )
        pivot = pivot[pivot.notna().sum(axis=1) >= 3]

        # Build lookup
        na_lookup = dict(zip(dgsa_df[id_col], dgsa_df['non_additivity']))

        ct_records = []
        for ct in pivot.columns:
            active = pivot[ct].dropna().index
            na_vals = [na_lookup[g] for g in active if g in na_lookup]
            if na_vals:
                ct_records.append({
                    'celltype': ct,
                    f'mean_na_{label}': np.mean(na_vals),
                    f'n_{label}': len(na_vals),
                })
        ct_df = pd.DataFrame(ct_records)
        merged = ct_df.merge(disease_ct, on='celltype', how='inner')
        results_both[label] = merged

        na_col = f'mean_na_{label}'
        rho, pval = stats.spearmanr(merged[na_col], merged['disease_count'])
        print(f"\n  {label}:")
        print(f"    Cell types with disease data: {len(merged)}")
        print(f"    Spearman rho: {rho:.4f}")
        print(f"    P-value:      {pval:.2e}")

    # Direct comparison: merge eQTL and caQTL cell-type tables
    e_ct = results_both['eQTL'][['celltype', 'mean_na_eQTL']].copy()
    c_ct = results_both['caQTL'][['celltype', 'mean_na_caQTL']].copy()
    both = e_ct.merge(c_ct, on='celltype', how='inner')
    if len(both) > 3:
        rho_ec, p_ec = stats.spearmanr(both['mean_na_eQTL'], both['mean_na_caQTL'])
        print(f"\n  eQTL vs caQTL cell-type non-additivity correlation:")
        print(f"    Cell types in both: {len(both)}")
        print(f"    Spearman rho: {rho_ec:.4f}")
        print(f"    P-value:      {p_ec:.2e}")


# ---------- main -------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Extension #1: deep distributional comparison')
    parser.add_argument('--s6', default='data/raw/CIMA_Table_S6.csv')
    parser.add_argument('--s15', default='data/raw/science.adt3130_table_s15.xlsx')
    parser.add_argument('--eqtl-results', default='results/dgsa_eqtl.csv')
    parser.add_argument('--caqtl-results', default='results/ext1_caqtl_dgsa.csv')
    parser.add_argument('--fig-output',
                        default='figures/ext1_distribution_comparison.pdf')
    args = parser.parse_args()

    eqtl_df = pd.read_csv(args.eqtl_results)
    caqtl_df = pd.read_csv(args.caqtl_results)

    eqtl_na = eqtl_df['non_additivity'].values
    caqtl_na = caqtl_df['non_additivity'].values

    # ---- Analysis 1: Distribution plots ----
    print("=" * 60)
    print("ANALYSIS 1: DISTRIBUTION COMPARISON")
    print("=" * 60)
    print(f"  eQTL:  mean={np.mean(eqtl_na):.4f}, median={np.median(eqtl_na):.4f}")
    print(f"  caQTL: mean={np.mean(caqtl_na):.4f}, median={np.median(caqtl_na):.4f}")

    # KS test
    ks_stat, ks_p = stats.ks_2samp(eqtl_na, caqtl_na)
    print(f"  KS test: D={ks_stat:.4f}, P={ks_p:.2e}")

    plot_distributions(eqtl_na, caqtl_na, args.fig_output)

    # ---- Analysis 2: Matched loci ----
    # Get cell type lists from the pivots
    raw = pd.read_csv(args.s6, usecols=['phenotype_id', 'celltype', 'analysis'])
    eqtl_cts = sorted(raw[raw['analysis'] == 'cis-eQTL']['celltype'].unique())
    caqtl_cts = sorted(raw[raw['analysis'] == 'cis-caQTL']['celltype'].unique())
    del raw

    matched_loci_analysis(args.s6, eqtl_cts, caqtl_cts)

    # ---- Analysis 3: Per-cell-type vs disease ----
    celltype_disease_analysis(
        caqtl_df, args.caqtl_results,
        eqtl_df, args.eqtl_results,
        args.s6, args.s15,
    )

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == '__main__':
    main()
