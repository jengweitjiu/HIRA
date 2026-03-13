#!/usr/bin/env python3
"""
DGSA — Decomposed Geometric Signature Analysis (eQTL)

Decomposes eQTL effect-size vectors across cell types into geometric metrics.
Each gene's slope vector spans all 69 cell types (NaN → 0 for absent eQTLs).

Key metric: non_additivity = Gini(slope²) in full cell-type space.

Input:  data/raw/CIMA_Table_S6.csv (analysis='cis-eQTL')
Output: results/dgsa_eqtl.csv
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
    """Compute DGSA metrics for a gene's slope vector in full cell-type space.

    Parameters
    ----------
    slopes_full : array of shape (n_celltypes_total,), NaN filled with 0
    n_active : number of cell types with non-zero slope
    """
    n = len(slopes_full)
    abs_slopes = np.abs(slopes_full)

    magnitude = np.linalg.norm(slopes_full)
    non_additivity = compute_gini(slopes_full ** 2)
    gini = compute_gini(abs_slopes)

    # Uniformity: cosine similarity with uniform vector (full space, abs slopes)
    uniform = np.ones(n) / np.sqrt(n)
    norm = np.linalg.norm(abs_slopes)
    uniformity = np.dot(abs_slopes, uniform) / norm if norm > 0 else 0.0

    # Sparsity: fraction of cell types with zero slope
    sparsity = 1.0 - (n_active / n)

    return {
        'n_celltypes': n_active,
        'magnitude': magnitude,
        'uniformity': uniformity,
        'non_additivity': non_additivity,
        'sparsity': sparsity,
        'gini': gini,
    }


def load_eqtl_data(s6_path, min_celltypes=3):
    """Load and pivot eQTL data from S6."""
    print(f"Reading {s6_path} ...")
    df = pd.read_csv(s6_path, usecols=['phenotype_id', 'celltype', 'slope', 'analysis'])
    print(f"  Total rows: {len(df):,}")
    print(f"  Columns: {list(df.columns)}")

    eqtl = df[df['analysis'] == 'cis-eQTL'].copy()
    print(f"  cis-eQTL rows: {len(eqtl):,}")
    print(f"  Unique genes: {eqtl['phenotype_id'].nunique():,}")
    print(f"  Unique cell types: {eqtl['celltype'].nunique()}")

    pivot = eqtl.pivot_table(
        index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
    )
    print(f"  Pivot shape (before filter): {pivot.shape}")

    gene_counts = pivot.notna().sum(axis=1)
    pivot_filtered = pivot[gene_counts >= min_celltypes].copy()
    print(f"  Genes in >= {min_celltypes} cell types: {len(pivot_filtered):,}")

    return pivot_filtered


def run_dgsa(pivot):
    """Run DGSA on pivoted slope matrix (NaN → 0 for full-space metrics)."""
    filled = pivot.fillna(0)
    n_active_per_gene = pivot.notna().sum(axis=1)

    results = []
    for gene in filled.index:
        slopes_full = filled.loc[gene].values
        n_active = n_active_per_gene[gene]
        metrics = compute_dgsa_metrics(slopes_full, n_active)
        metrics['gene'] = gene
        results.append(metrics)

    dgsa_df = pd.DataFrame(results)
    dgsa_df = dgsa_df[['gene', 'n_celltypes', 'magnitude', 'uniformity',
                        'non_additivity', 'sparsity', 'gini']]
    dgsa_df = dgsa_df.sort_values('non_additivity', ascending=False).reset_index(drop=True)
    return dgsa_df


def correlate_with_disease(dgsa_df, pivot, s15_path):
    """Per-cell-type Spearman: mean non-additivity vs disease trait count."""
    print(f"\nReading {s15_path} for disease correlation ...")
    s15 = pd.read_excel(s15_path, sheet_name='Sheet1')
    print(f"  S15 columns: {list(s15.columns)}")
    print(f"  S15 rows: {len(s15):,}")

    # Build gene→non_additivity lookup
    gene_na = dict(zip(dgsa_df['gene'], dgsa_df['non_additivity']))

    # Per cell type: mean non-additivity of genes active in that cell type
    ct_records = []
    for ct in pivot.columns:
        active_genes = pivot[ct].dropna().index
        na_vals = [gene_na[g] for g in active_genes if g in gene_na]
        if na_vals:
            ct_records.append({
                'celltype': ct,
                'mean_non_additivity': np.mean(na_vals),
                'n_genes': len(na_vals),
            })
    ct_df = pd.DataFrame(ct_records)

    # Disease trait count per cell type
    disease_ct = s15.groupby('celltype')['trait'].nunique().reset_index()
    disease_ct.columns = ['celltype', 'disease_count']
    print(f"  Cell types with disease data: {len(disease_ct)}")

    merged = ct_df.merge(disease_ct, on='celltype', how='inner')
    print(f"  Cell types in both DGSA and S15: {len(merged)}")

    if len(merged) < 3:
        print("  WARNING: Too few overlapping cell types for correlation.")
        return None, None, len(merged)

    rho, pval = stats.spearmanr(merged['mean_non_additivity'], merged['disease_count'])
    return rho, pval, len(merged)


def main():
    parser = argparse.ArgumentParser(description='DGSA: eQTL geometric decomposition')
    parser.add_argument('--s6', default='data/raw/CIMA_Table_S6.csv',
                        help='Path to S6 CSV file')
    parser.add_argument('--s15', default='data/raw/science.adt3130_table_s15.xlsx',
                        help='Path to S15 Excel file')
    parser.add_argument('--min-celltypes', type=int, default=3,
                        help='Minimum cell types per gene (default: 3)')
    parser.add_argument('--output', default='results/dgsa_eqtl.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # Load and pivot
    pivot = load_eqtl_data(args.s6, args.min_celltypes)

    # Compute DGSA metrics
    print("\nComputing DGSA metrics ...")
    dgsa_df = run_dgsa(pivot)

    # Save
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    dgsa_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(dgsa_df):,} genes to {args.output}")

    # Key statistics
    print("\n" + "=" * 60)
    print("DGSA-eQTL KEY STATISTICS")
    print("=" * 60)
    print(f"  Total genes analyzed:     {len(dgsa_df):,}")
    print(f"  Mean non-additivity:      {dgsa_df['non_additivity'].mean():.4f}")
    print(f"  Median non-additivity:    {dgsa_df['non_additivity'].median():.4f}")
    print(f"  Std non-additivity:       {dgsa_df['non_additivity'].std():.4f}")
    print(f"  Mean magnitude:           {dgsa_df['magnitude'].mean():.4f}")
    print(f"  Mean Gini:                {dgsa_df['gini'].mean():.4f}")
    print(f"  Mean sparsity:            {dgsa_df['sparsity'].mean():.4f}")

    # Disease correlation (per cell type)
    rho, pval, n_ct = correlate_with_disease(dgsa_df, pivot, args.s15)
    if rho is not None:
        print(f"\n  Disease correlation (per cell type, n={n_ct}):")
        print(f"    Spearman rho: {rho:.4f}")
        print(f"    P-value:      {pval:.2e}")
    print("=" * 60)


if __name__ == '__main__':
    main()
