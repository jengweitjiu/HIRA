#!/usr/bin/env python3
"""
TOPPLE — Transcriptional Orchestration via Perturbation-Probed Landscape Entropy

For each regulon, measure the redistribution of cell-type AUC profiles when that
regulon is removed (leave-one-out JSD perturbation).

Input:  data/raw/science.adt3130_table_s5.xlsx, sheet "eRegulons_Activators_Exp_AUC_RS"
Output: results/topple.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from pathlib import Path


def compute_ri(auc_matrix):
    """Compute redistribution index for each regulon via leave-one-out JSD.

    Parameters
    ----------
    auc_matrix : DataFrame, shape (n_regulons, n_celltypes), values = mean_AUC

    Returns
    -------
    DataFrame with regulon, mean_ri, per-celltype RI values
    """
    regulons = auc_matrix.index.tolist()
    celltypes = auc_matrix.columns.tolist()
    n_reg = len(regulons)

    # L1-normalize each cell type column (profile across regulons sums to 1)
    col_sums = auc_matrix.sum(axis=0)
    norm_matrix = auc_matrix.div(col_sums, axis=1)
    # norm_matrix shape: (203 regulons, 61 cell types)
    # Each column is a probability distribution over regulons

    results = []
    for idx_reg, reg in enumerate(regulons):
        ri_values = []
        for ct in celltypes:
            # Original L1-normalized profile for this cell type (203 regulons)
            p = norm_matrix[ct].values.copy()

            # Perturbed profile: zero out the removed regulon, renormalize
            q = p.copy()
            q[idx_reg] = 0.0
            q_sum = q.sum()
            if q_sum > 0:
                q = q / q_sum
            else:
                continue

            # JSD between original and perturbed profiles (both 203-dim)
            jsd_val = jensenshannon(p, q)
            if np.isfinite(jsd_val):
                ri_values.append(jsd_val ** 2)  # jensenshannon returns sqrt(JSD)

        mean_ri = np.mean(ri_values) if ri_values else 0.0
        results.append({
            'regulon': reg,
            'mean_RI': mean_ri,
            'std_RI': np.std(ri_values) if ri_values else 0.0,
            'max_RI': np.max(ri_values) if ri_values else 0.0,
            'min_RI': np.min(ri_values) if ri_values else 0.0,
        })

    df = pd.DataFrame(results)
    df = df.sort_values('mean_RI', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1
    return df


def classify_stability(topple_df):
    """Classify regulons into quartile-based stability classes."""
    q75 = topple_df['mean_RI'].quantile(0.75)
    q25 = topple_df['mean_RI'].quantile(0.25)
    topple_df['stability_class'] = 'intermediate'
    topple_df.loc[topple_df['mean_RI'] >= q75, 'stability_class'] = 'stabilizer'
    topple_df.loc[topple_df['mean_RI'] <= q25, 'stability_class'] = 'destabilizer'
    return topple_df


def main():
    parser = argparse.ArgumentParser(description='TOPPLE: regulon stability analysis')
    parser.add_argument('--s5', default='data/raw/science.adt3130_table_s5.xlsx',
                        help='Path to S5 Excel file')
    parser.add_argument('--output', default='results/topple.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # Load data
    print(f"Reading {args.s5}, sheet 'eRegulons_Activators_Exp_AUC_RS' ...")
    df = pd.read_excel(args.s5, sheet_name='eRegulons_Activators_Exp_AUC_RS')
    print(f"  Columns: {list(df.columns)}")
    print(f"  Shape: {df.shape}")
    print(f"  Unique eRegulons: {df['eRegulon'].nunique()}")
    print(f"  Unique cell types: {df['cell_type_l4'].nunique()}")

    # Pivot: regulon x celltype with mean_AUC
    auc_matrix = df.pivot_table(
        index='eRegulon', columns='cell_type_l4', values='mean_AUC', aggfunc='first'
    )
    print(f"  AUC matrix shape: {auc_matrix.shape}")

    # Compute TOPPLE
    print("\nComputing redistribution index (leave-one-out JSD) ...")
    topple_df = compute_ri(auc_matrix)
    topple_df = classify_stability(topple_df)

    # Save
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    topple_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(topple_df)} regulons to {args.output}")

    # Key statistics
    print("\n" + "=" * 60)
    print("TOPPLE KEY STATISTICS")
    print("=" * 60)
    print(f"  Total regulons:    {len(topple_df)}")
    print(f"  Mean RI:           {topple_df['mean_RI'].mean():.6f}")
    print(f"  Stabilizers:       {(topple_df['stability_class']=='stabilizer').sum()}")
    print(f"  Destabilizers:     {(topple_df['stability_class']=='destabilizer').sum()}")

    # Top 5 and HSF1/STAT5B
    print(f"\n  Top 10 regulons by RI:")
    for _, row in topple_df.head(10).iterrows():
        print(f"    #{int(row['rank']):3d}  {row['regulon']:25s}  RI={row['mean_RI']:.6f}  [{row['stability_class']}]")

    for target in ['HSF1_+', 'STAT5B_+']:
        match = topple_df[topple_df['regulon'] == target]
        if len(match) > 0:
            r = match.iloc[0]
            print(f"\n  {target}: rank #{int(r['rank'])}, RI={r['mean_RI']:.6f}")
        else:
            # Try without suffix
            for suffix in ['_+', '_-', '']:
                match = topple_df[topple_df['regulon'].str.startswith(target.replace('_+',''))]
                if len(match) > 0:
                    for _, r in match.iterrows():
                        print(f"\n  {r['regulon']}: rank #{int(r['rank'])}, RI={r['mean_RI']:.6f}")
                    break

    print("=" * 60)


if __name__ == '__main__':
    main()
