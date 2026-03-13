#!/usr/bin/env python3
"""
Extension #1: caQTL Geometric Architecture

Applies the same DGSA pipeline to cis-caQTLs (peaks instead of genes).
Tests the "epigenomic amplification" hypothesis: caQTL non-additivity > eQTL.

Input:  data/raw/CIMA_Table_S6.csv (analysis='cis-caQTL')
        results/dgsa_eqtl.csv (for comparison)
Output: results/ext1_caqtl_dgsa.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


def compute_gini(arr):
    """Gini coefficient for a 1-D array of non-negative values."""
    arr = np.sort(arr.copy())
    n = len(arr)
    if n == 0 or arr.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * arr) - (n + 1) * np.sum(arr)) / (n * np.sum(arr))


def compute_dgsa_metrics(slopes_full, n_active):
    """Compute DGSA metrics for a peak's slope vector in full cell-type space."""
    n = len(slopes_full)
    abs_slopes = np.abs(slopes_full)

    magnitude = np.linalg.norm(slopes_full)
    non_additivity = compute_gini(slopes_full ** 2)
    gini = compute_gini(abs_slopes)

    uniform = np.ones(n) / np.sqrt(n)
    norm = np.linalg.norm(abs_slopes)
    uniformity = np.dot(abs_slopes, uniform) / norm if norm > 0 else 0.0

    sparsity = 1.0 - (n_active / n)

    return {
        'n_celltypes': n_active,
        'magnitude': magnitude,
        'uniformity': uniformity,
        'non_additivity': non_additivity,
        'sparsity': sparsity,
        'gini': gini,
    }


def load_caqtl_data(s6_path, min_celltypes=3):
    """Load and pivot caQTL data from S6."""
    print(f"Reading {s6_path} ...")
    df = pd.read_csv(s6_path, usecols=['phenotype_id', 'celltype', 'slope', 'analysis'])
    print(f"  Total rows: {len(df):,}")

    caqtl = df[df['analysis'] == 'cis-caQTL'].copy()
    print(f"  cis-caQTL rows: {len(caqtl):,}")
    print(f"  Unique peaks: {caqtl['phenotype_id'].nunique():,}")
    print(f"  Unique cell types: {caqtl['celltype'].nunique()}")

    pivot = caqtl.pivot_table(
        index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
    )
    print(f"  Pivot shape (before filter): {pivot.shape}")

    gene_counts = pivot.notna().sum(axis=1)
    pivot_filtered = pivot[gene_counts >= min_celltypes].copy()
    print(f"  Peaks in >= {min_celltypes} cell types: {len(pivot_filtered):,}")

    return pivot_filtered


def run_dgsa(pivot):
    """Run DGSA on pivoted slope matrix (NaN → 0 for full-space metrics)."""
    filled = pivot.fillna(0)
    n_active_per_peak = pivot.notna().sum(axis=1)

    results = []
    for peak in filled.index:
        slopes_full = filled.loc[peak].values
        n_active = n_active_per_peak[peak]
        metrics = compute_dgsa_metrics(slopes_full, n_active)
        metrics['peak'] = peak
        results.append(metrics)

    dgsa_df = pd.DataFrame(results)
    dgsa_df = dgsa_df[['peak', 'n_celltypes', 'magnitude', 'uniformity',
                        'non_additivity', 'sparsity', 'gini']]
    dgsa_df = dgsa_df.sort_values('non_additivity', ascending=False).reset_index(drop=True)
    return dgsa_df


def compare_with_eqtl(caqtl_df, eqtl_path):
    """Mann-Whitney U test: caQTL non-additivity vs eQTL non-additivity."""
    print(f"\nReading eQTL results from {eqtl_path} ...")
    eqtl_df = pd.read_csv(eqtl_path)
    print(f"  eQTL genes: {len(eqtl_df):,}")

    eqtl_na = eqtl_df['non_additivity'].values
    caqtl_na = caqtl_df['non_additivity'].values

    # Mann-Whitney U tests
    u_stat, p_two = stats.mannwhitneyu(caqtl_na, eqtl_na, alternative='two-sided')
    _, p_greater = stats.mannwhitneyu(caqtl_na, eqtl_na, alternative='greater')
    _, p_less = stats.mannwhitneyu(caqtl_na, eqtl_na, alternative='less')

    # Effect size: rank-biserial correlation
    # r > 0 means caQTL tends to be larger; r < 0 means eQTL tends to be larger
    n1, n2 = len(caqtl_na), len(eqtl_na)
    rank_biserial = 1 - (2 * u_stat) / (n1 * n2)

    return {
        'eqtl_mean': np.mean(eqtl_na),
        'eqtl_median': np.median(eqtl_na),
        'eqtl_n': len(eqtl_na),
        'caqtl_mean': np.mean(caqtl_na),
        'caqtl_median': np.median(caqtl_na),
        'caqtl_n': len(caqtl_na),
        'U_statistic': u_stat,
        'p_twosided': p_two,
        'p_greater': p_greater,
        'p_less': p_less,
        'rank_biserial_r': rank_biserial,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Extension #1: caQTL geometric architecture')
    parser.add_argument('--s6', default='data/raw/CIMA_Table_S6.csv',
                        help='Path to S6 CSV file')
    parser.add_argument('--eqtl-results', default='results/dgsa_eqtl.csv',
                        help='Path to DGSA-eQTL results for comparison')
    parser.add_argument('--min-celltypes', type=int, default=3,
                        help='Minimum cell types per peak (default: 3)')
    parser.add_argument('--output', default='results/ext1_caqtl_dgsa.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # Load and pivot caQTL data
    pivot = load_caqtl_data(args.s6, args.min_celltypes)

    # Compute DGSA metrics
    print("\nComputing DGSA metrics for caQTLs ...")
    caqtl_df = run_dgsa(pivot)

    # Save
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    caqtl_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(caqtl_df):,} peaks to {args.output}")

    # Key statistics
    print("\n" + "=" * 60)
    print("EXTENSION #1: caQTL DGSA KEY STATISTICS")
    print("=" * 60)
    print(f"  Total peaks analyzed:     {len(caqtl_df):,}")
    print(f"  Mean non-additivity:      {caqtl_df['non_additivity'].mean():.4f}")
    print(f"  Median non-additivity:    {caqtl_df['non_additivity'].median():.4f}")
    print(f"  Std non-additivity:       {caqtl_df['non_additivity'].std():.4f}")
    print(f"  Mean magnitude:           {caqtl_df['magnitude'].mean():.4f}")
    print(f"  Mean Gini:                {caqtl_df['gini'].mean():.4f}")
    print(f"  Mean sparsity:            {caqtl_df['sparsity'].mean():.4f}")

    # Compare with eQTL
    if Path(args.eqtl_results).exists():
        comparison = compare_with_eqtl(caqtl_df, args.eqtl_results)
        print(f"\n  --- Mann-Whitney: caQTL vs eQTL non-additivity ---")
        print(f"  eQTL  mean={comparison['eqtl_mean']:.4f}  "
              f"median={comparison['eqtl_median']:.4f}  n={comparison['eqtl_n']:,}")
        print(f"  caQTL mean={comparison['caqtl_mean']:.4f}  "
              f"median={comparison['caqtl_median']:.4f}  n={comparison['caqtl_n']:,}")
        print(f"  U statistic:        {comparison['U_statistic']:.0f}")
        print(f"  P (two-sided):      {comparison['p_twosided']:.2e}")
        print(f"  P (caQTL > eQTL):   {comparison['p_greater']:.2e}")
        print(f"  P (caQTL < eQTL):   {comparison['p_less']:.2e}")
        print(f"  Rank-biserial r:    {comparison['rank_biserial_r']:.4f}")
        print(f"    (r > 0: caQTL tends larger; r < 0: eQTL tends larger)")
        if comparison['p_greater'] < 0.05:
            print(f"\n  >> CONFIRMED: caQTL non-additivity > eQTL "
                  f"(epigenomic amplification)")
        elif comparison['p_less'] < 0.05:
            print(f"\n  >> REVERSED: eQTL non-additivity > caQTL "
                  f"(P={comparison['p_less']:.2e})")
        else:
            print(f"\n  >> No significant difference")
    else:
        print(f"\n  WARNING: {args.eqtl_results} not found, skipping comparison.")

    print("=" * 60)


if __name__ == '__main__':
    main()
