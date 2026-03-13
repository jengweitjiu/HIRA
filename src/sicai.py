#!/usr/bin/env python3
"""
SICAI — Shared Immune Cellular Architecture Index

Builds a 69x69 eQTL sharing matrix (r_b) and computes per-cell-type coupling
metrics: mean r_b, CV, Shannon entropy. Correlates with disease burden from S15.

Input:  data/raw/science.adt3130_table_s8.xlsx, sheet "cis_eQTL"
        data/raw/science.adt3130_table_s15.xlsx
Output: results/sicai_rb.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


def build_rb_matrix(s8_df):
    """Build symmetric r_b matrix from pairwise data."""
    celltypes = sorted(set(s8_df['reference_cell_type']) | set(s8_df['query_celltype']))
    n = len(celltypes)
    ct_idx = {ct: i for i, ct in enumerate(celltypes)}

    rb_mat = np.full((n, n), np.nan)
    np.fill_diagonal(rb_mat, 1.0)

    for _, row in s8_df.iterrows():
        i = ct_idx[row['reference_cell_type']]
        j = ct_idx[row['query_celltype']]
        rb_mat[j, i] = row['rb']  # query is row, reference is column
        if np.isnan(rb_mat[i, j]):
            rb_mat[i, j] = row['rb']  # fill symmetric if missing

    return pd.DataFrame(rb_mat, index=celltypes, columns=celltypes)


def compute_sicai_metrics(rb_matrix):
    """Per-cell-type coupling metrics from r_b matrix."""
    celltypes = rb_matrix.index.tolist()
    results = []

    for ct in celltypes:
        # Get all r_b values for this cell type (excluding self)
        rb_vals = rb_matrix.loc[ct].drop(ct).dropna().values

        if len(rb_vals) == 0:
            continue

        mean_rb = np.mean(rb_vals)
        std_rb = np.std(rb_vals)
        cv_rb = std_rb / mean_rb if mean_rb > 0 else 0.0

        # Shannon entropy of the r_b distribution (discretize into bins)
        # Normalize r_b to probability-like distribution
        rb_pos = np.clip(rb_vals, 0, None)
        rb_sum = rb_pos.sum()
        if rb_sum > 0:
            p = rb_pos / rb_sum
            entropy = -np.sum(p * np.log(p + 1e-30))
        else:
            entropy = 0.0

        results.append({
            'cell_type': ct,
            'mean_rb': mean_rb,
            'median_rb': np.median(rb_vals),
            'std_rb': std_rb,
            'cv_rb': cv_rb,
            'entropy': entropy,
            'min_rb': np.min(rb_vals),
            'max_rb': np.max(rb_vals),
            'n_pairs': len(rb_vals),
        })

    return pd.DataFrame(results).sort_values('mean_rb', ascending=False).reset_index(drop=True)


def correlate_with_disease(sicai_df, s15_path):
    """Spearman: per-cell-type mean r_b vs disease trait count."""
    print(f"\nReading {s15_path} for disease correlation ...")
    s15 = pd.read_excel(s15_path, sheet_name='Sheet1')
    print(f"  S15 columns: {list(s15.columns)}")

    disease_ct = s15.groupby('celltype')['trait'].nunique().reset_index()
    disease_ct.columns = ['cell_type', 'disease_count']
    print(f"  Cell types with disease data: {len(disease_ct)}")

    merged = sicai_df.merge(disease_ct, on='cell_type', how='inner')
    print(f"  Cell types in both SICAI and S15: {len(merged)}")

    if len(merged) < 3:
        return None, None, 0

    rho, pval = stats.spearmanr(merged['mean_rb'], merged['disease_count'])
    return rho, pval, len(merged)


def main():
    parser = argparse.ArgumentParser(description='SICAI: eQTL sharing topology')
    parser.add_argument('--s8', default='data/raw/science.adt3130_table_s8.xlsx',
                        help='Path to S8 Excel file')
    parser.add_argument('--s15', default='data/raw/science.adt3130_table_s15.xlsx',
                        help='Path to S15 Excel file')
    parser.add_argument('--output', default='results/sicai_rb.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # Load S8
    print(f"Reading {args.s8}, sheet 'cis_eQTL' ...")
    s8 = pd.read_excel(args.s8, sheet_name='cis_eQTL')
    print(f"  Columns: {list(s8.columns)}")
    print(f"  Shape: {s8.shape}")

    # Build r_b matrix
    print("\nBuilding r_b matrix ...")
    rb_matrix = build_rb_matrix(s8)
    print(f"  Matrix shape: {rb_matrix.shape}")
    print(f"  Global mean r_b (off-diagonal): {rb_matrix.values[~np.eye(69, dtype=bool)].mean():.4f}")

    # Compute per-cell-type metrics
    print("\nComputing SICAI metrics ...")
    sicai_df = compute_sicai_metrics(rb_matrix)

    # Save
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    sicai_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(sicai_df)} cell types to {args.output}")

    # Key statistics
    print("\n" + "=" * 60)
    print("SICAI KEY STATISTICS")
    print("=" * 60)
    print(f"  Cell types:       {len(sicai_df)}")
    print(f"  Mean r_b:         {sicai_df['mean_rb'].mean():.4f}")
    print(f"  Median r_b:       {sicai_df['median_rb'].mean():.4f}")
    print(f"  Mean CV:          {sicai_df['cv_rb'].mean():.4f}")
    print(f"  Mean entropy:     {sicai_df['entropy'].mean():.4f}")

    print(f"\n  Top 5 most coupled cell types:")
    for _, row in sicai_df.head(5).iterrows():
        print(f"    {row['cell_type']:30s}  mean_rb={row['mean_rb']:.4f}")

    print(f"\n  Bottom 5 least coupled cell types:")
    for _, row in sicai_df.tail(5).iterrows():
        print(f"    {row['cell_type']:30s}  mean_rb={row['mean_rb']:.4f}")

    # Disease correlation
    rho, pval, n = correlate_with_disease(sicai_df, args.s15)
    if rho is not None:
        print(f"\n  Disease correlation (n={n}):")
        print(f"    Spearman rho: {rho:.4f}")
        print(f"    P-value:      {pval:.2e}")
    print("=" * 60)


if __name__ == '__main__':
    main()
