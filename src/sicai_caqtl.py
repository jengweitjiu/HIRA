#!/usr/bin/env python3
"""
SICAI-caQTL — caQTL sharing topology + Mantel test vs eQTL rb.

Input:  data/raw/science.adt3130_table_s8.xlsx (sheets: cis_eQTL, cis_caQTL)
Output: results/sicai_caqtl_rb.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import squareform
from pathlib import Path


def build_rb_matrix(s8_df):
    """Build rb matrix from pairwise data. Returns DataFrame."""
    celltypes = sorted(set(s8_df['reference_cell_type']) | set(s8_df['query_celltype']))
    n = len(celltypes)
    ct_idx = {ct: i for i, ct in enumerate(celltypes)}

    rb_mat = np.full((n, n), np.nan)
    np.fill_diagonal(rb_mat, 1.0)

    for _, row in s8_df.iterrows():
        i = ct_idx[row['reference_cell_type']]
        j = ct_idx[row['query_celltype']]
        rb_mat[j, i] = row['rb']
        if np.isnan(rb_mat[i, j]):
            rb_mat[i, j] = row['rb']

    return pd.DataFrame(rb_mat, index=celltypes, columns=celltypes)


def compute_sicai_metrics(rb_matrix):
    """Per-cell-type coupling metrics."""
    celltypes = rb_matrix.index.tolist()
    results = []
    for ct in celltypes:
        rb_vals = rb_matrix.loc[ct].drop(ct).dropna().values
        if len(rb_vals) == 0:
            continue
        mean_rb = np.mean(rb_vals)
        std_rb = np.std(rb_vals)
        cv_rb = std_rb / mean_rb if mean_rb > 0 else 0.0
        rb_pos = np.clip(rb_vals, 0, None)
        rb_sum = rb_pos.sum()
        p = rb_pos / rb_sum if rb_sum > 0 else np.ones(len(rb_vals)) / len(rb_vals)
        entropy = -np.sum(p * np.log(p + 1e-30))
        results.append({
            'cell_type': ct, 'mean_rb': mean_rb, 'median_rb': np.median(rb_vals),
            'std_rb': std_rb, 'cv_rb': cv_rb, 'entropy': entropy,
            'n_pairs': len(rb_vals),
        })
    return pd.DataFrame(results).sort_values('mean_rb', ascending=False).reset_index(drop=True)


def mantel_test(mat1, mat2, n_perm=9999):
    """Mantel test between two distance/similarity matrices on shared cell types."""
    shared = sorted(set(mat1.index) & set(mat2.index))
    n = len(shared)
    print(f"  Shared cell types for Mantel: {n}")
    if n < 4:
        print("  Too few shared cell types.")
        return np.nan, np.nan, 0

    m1 = mat1.loc[shared, shared].values
    m2 = mat2.loc[shared, shared].values

    # Extract upper triangle (off-diagonal)
    triu_idx = np.triu_indices(n, k=1)
    v1 = m1[triu_idx]
    v2 = m2[triu_idx]

    # Drop pairs where either is NaN
    valid = ~(np.isnan(v1) | np.isnan(v2))
    v1 = v1[valid]
    v2 = v2[valid]
    print(f"  Valid pairs: {len(v1)}")

    # Observed Pearson correlation
    r_obs, _ = stats.pearsonr(v1, v2)

    # Permutation test
    count = 0
    for _ in range(n_perm):
        perm = np.random.permutation(n)
        m2_perm = m2[np.ix_(perm, perm)]
        v2_perm = m2_perm[triu_idx][valid]
        r_perm, _ = stats.pearsonr(v1, v2_perm)
        if r_perm >= r_obs:
            count += 1

    p_val = (count + 1) / (n_perm + 1)
    return r_obs, p_val, len(v1)


def main():
    parser = argparse.ArgumentParser(description='SICAI-caQTL + Mantel test')
    parser.add_argument('--s8', default='data/raw/science.adt3130_table_s8.xlsx')
    parser.add_argument('--output', default='results/sicai_caqtl_rb.csv')
    args = parser.parse_args()

    # Load caQTL
    print(f"Reading {args.s8}, sheet 'cis_caQTL' ...")
    caqtl = pd.read_excel(args.s8, sheet_name='cis_caQTL')
    print(f"  Columns: {list(caqtl.columns)}")
    print(f"  Rows: {len(caqtl)}, Cell types: {caqtl['reference_cell_type'].nunique()}")

    rb_caqtl = build_rb_matrix(caqtl)
    print(f"  caQTL rb matrix: {rb_caqtl.shape}")

    sicai_df = compute_sicai_metrics(rb_caqtl)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    sicai_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(sicai_df)} cell types to {args.output}")

    # Load eQTL for Mantel
    print(f"\nReading {args.s8}, sheet 'cis_eQTL' ...")
    eqtl = pd.read_excel(args.s8, sheet_name='cis_eQTL')
    rb_eqtl = build_rb_matrix(eqtl)
    print(f"  eQTL rb matrix: {rb_eqtl.shape}")

    # Print key stats
    print("\n" + "=" * 60)
    print("SICAI-caQTL KEY STATISTICS")
    print("=" * 60)
    n_ct = rb_caqtl.shape[0]
    offdiag = rb_caqtl.values[~np.eye(n_ct, dtype=bool)]
    offdiag = offdiag[~np.isnan(offdiag)]
    print(f"  Cell types:          {n_ct}")
    print(f"  Global mean rb:      {np.mean(offdiag):.4f}")
    print(f"  Global median rb:    {np.median(offdiag):.4f}")
    print(f"\n  Per-cell-type mean rb: {sicai_df['mean_rb'].mean():.4f}")

    print(f"\n  Top 5 most coupled:")
    for _, r in sicai_df.head(5).iterrows():
        print(f"    {r['cell_type']:30s}  mean_rb={r['mean_rb']:.4f}")
    print(f"  Bottom 5:")
    for _, r in sicai_df.tail(5).iterrows():
        print(f"    {r['cell_type']:30s}  mean_rb={r['mean_rb']:.4f}")

    # Mantel test
    print(f"\n  --- Mantel test: eQTL rb vs caQTL rb ---")
    np.random.seed(42)
    r_mantel, p_mantel, n_pairs = mantel_test(rb_eqtl, rb_caqtl)
    print(f"  Mantel r: {r_mantel:.4f}")
    print(f"  Mantel P: {p_mantel:.4f} (permutation, 9999 perms)")
    print(f"  N pairs:  {n_pairs}")
    print("=" * 60)


if __name__ == '__main__':
    main()
