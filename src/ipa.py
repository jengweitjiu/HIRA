#!/usr/bin/env python3
"""
IPA — Immune Perturbation Architecture

Compares perturbation resistance (sex-difference |Log2_Fold_Change|) between
TOPPLE-classified stabilizers and destabilizers. Tests whether high-RI regulons
(stabilizers) are also more resistant to sex-based perturbation.

Input:  data/raw/science.adt3130_table_s5.xlsx, sheet "sex_difference_results"
        results/topple.csv (for stability classifications)
Output: results/ipa_sex.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='IPA: perturbation resistance')
    parser.add_argument('--s5', default='data/raw/science.adt3130_table_s5.xlsx',
                        help='Path to S5 Excel file')
    parser.add_argument('--topple', default='results/topple.csv',
                        help='Path to TOPPLE results')
    parser.add_argument('--output', default='results/ipa_sex.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # Load sex-difference data
    print(f"Reading {args.s5}, sheet 'sex_difference_results' ...")
    sex_df = pd.read_excel(args.s5, sheet_name='sex_difference_results')
    print(f"  Columns: {list(sex_df.columns)}")
    print(f"  Shape: {sex_df.shape}")
    print(f"  Unique regulons: {sex_df['eRegulon'].nunique()}")
    print(f"  Unique cell types: {sex_df['cell_type_l4'].nunique()}")

    sex_df['abs_log2FC'] = sex_df['Log2_Fold_Change'].abs()

    # Normalize regulon names: strip suffixes like "_(290g)" or "_extended"
    sex_df['regulon_clean'] = (sex_df['eRegulon']
                               .str.replace(r'_extended', '', regex=False)
                               .str.replace(r'_\(\d+g\)', '', regex=True))

    # Load TOPPLE classifications
    print(f"\nReading {args.topple} ...")
    topple = pd.read_csv(args.topple)
    print(f"  TOPPLE regulons: {len(topple)}")
    print(f"  Stability classes: {topple['stability_class'].value_counts().to_dict()}")

    # Merge on cleaned names
    merged = sex_df.merge(
        topple[['regulon', 'mean_RI', 'stability_class', 'rank']],
        left_on='regulon_clean', right_on='regulon', how='inner'
    )
    print(f"\n  Merged rows: {len(merged):,}")
    print(f"  Regulons matched: {merged['eRegulon'].nunique()}")

    # Split by stability class
    stab = merged[merged['stability_class'] == 'stabilizer']
    dest = merged[merged['stability_class'] == 'destabilizer']

    print(f"\n  Stabilizer rows:    {len(stab):,}  ({stab['eRegulon'].nunique()} regulons)")
    print(f"  Destabilizer rows:  {len(dest):,}  ({dest['eRegulon'].nunique()} regulons)")

    # Mann-Whitney: |log2FC| stabilizers vs destabilizers (drop NaN)
    stab_vals = stab['abs_log2FC'].dropna().values
    dest_vals = dest['abs_log2FC'].dropna().values
    print(f"  After dropping NaN: stab={len(stab_vals)}, dest={len(dest_vals)}")

    u_stat, p_two = stats.mannwhitneyu(
        stab_vals, dest_vals, alternative='two-sided', method='asymptotic'
    )
    _, p_less = stats.mannwhitneyu(
        stab_vals, dest_vals, alternative='less', method='asymptotic'
    )
    _, p_greater = stats.mannwhitneyu(
        stab_vals, dest_vals, alternative='greater', method='asymptotic'
    )

    n1, n2 = len(stab_vals), len(dest_vals)
    rank_biserial = 1 - (2 * u_stat) / (n1 * n2)

    # Per-regulon summary
    reg_summary = merged.groupby(['eRegulon', 'stability_class', 'mean_RI', 'rank']).agg(
        mean_abs_log2FC=('abs_log2FC', 'mean'),
        median_abs_log2FC=('abs_log2FC', 'median'),
        n_tests=('abs_log2FC', 'count'),
        n_significant=('adjust_P_Value', lambda x: (x < 0.05).sum()),
    ).reset_index()
    reg_summary = reg_summary.sort_values('rank')

    # Save
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    reg_summary.to_csv(args.output, index=False)
    print(f"\nSaved per-regulon summary to {args.output}")

    # Spearman: RI vs mean |log2FC| per regulon
    rho_ri, p_ri = stats.spearmanr(reg_summary['mean_RI'], reg_summary['mean_abs_log2FC'])

    # Key statistics
    print("\n" + "=" * 60)
    print("IPA KEY STATISTICS — Sex Perturbation")
    print("=" * 60)
    print(f"  Stabilizer |log2FC|:   mean={stab['abs_log2FC'].mean():.4f}  "
          f"median={stab['abs_log2FC'].median():.4f}")
    print(f"  Destabilizer |log2FC|: mean={dest['abs_log2FC'].mean():.4f}  "
          f"median={dest['abs_log2FC'].median():.4f}")
    print(f"\n  Mann-Whitney U test (stabilizer vs destabilizer):")
    print(f"    U statistic:      {u_stat:.0f}")
    print(f"    P (two-sided):    {p_two:.2e}")
    print(f"    P (stab < dest):  {p_less:.2e}")
    print(f"    P (stab > dest):  {p_greater:.2e}")
    print(f"    Rank-biserial r:  {rank_biserial:.4f}")
    print(f"      (r < 0 means stabilizers have smaller |log2FC|)")
    print(f"\n  Spearman (RI vs mean |log2FC| per regulon):")
    print(f"    rho={rho_ri:.4f}, P={p_ri:.2e}")
    print("=" * 60)


if __name__ == '__main__':
    main()
