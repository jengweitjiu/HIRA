#!/usr/bin/env python3
"""
IPA-Age — Immune Perturbation Architecture (Age Correlation)

Compares age-correlation magnitude (|Spearman rho|) between TOPPLE-classified
stabilizers and destabilizers. Tests whether high-RI regulons (stabilizers)
show stronger age-associated perturbation than destabilizers.

Input:  data/raw/science.adt3130_table_s5.xlsx, sheet "Age_Correlation"
        results/topple.csv (for stability classifications)
Output: results/ipa_age.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='IPA-Age: age perturbation resistance')
    parser.add_argument('--s5', default='data/raw/science.adt3130_table_s5.xlsx',
                        help='Path to S5 Excel file')
    parser.add_argument('--topple', default='results/topple.csv',
                        help='Path to TOPPLE results')
    parser.add_argument('--output', default='results/ipa_age.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # ── Load Age_Correlation sheet (wide format) ──────────────────────────
    print(f"Reading {args.s5}, sheet 'Age_Correlation' ...")
    wide = pd.read_excel(args.s5, sheet_name='Age_Correlation')
    print(f"  Columns: {list(wide.columns)}")
    print(f"  Shape: {wide.shape}")

    # Identify cell-type rho columns (non-Pvalue, non-index)
    rho_cols = [c for c in wide.columns
                if c != 'Unnamed: 0' and not c.endswith('_Pvalue')]
    pval_cols = [c + '_Pvalue' for c in rho_cols]
    print(f"  Cell types detected: {len(rho_cols)}")
    print(f"  Regulons: {wide.shape[0]}")

    # ── Melt to long format ───────────────────────────────────────────────
    # Melt rho values
    rho_long = wide.melt(
        id_vars='Unnamed: 0', value_vars=rho_cols,
        var_name='cell_type', value_name='age_rho'
    ).rename(columns={'Unnamed: 0': 'eRegulon'})

    # Melt pvalue values
    pval_long = wide.melt(
        id_vars='Unnamed: 0', value_vars=pval_cols,
        var_name='cell_type_pv', value_name='age_pvalue'
    ).rename(columns={'Unnamed: 0': 'eRegulon'})
    # Strip '_Pvalue' to align cell_type names
    pval_long['cell_type'] = pval_long['cell_type_pv'].str.replace('_Pvalue', '', regex=False)

    # Merge rho and pvalue
    age_df = rho_long.merge(pval_long[['eRegulon', 'cell_type', 'age_pvalue']],
                            on=['eRegulon', 'cell_type'], how='left')
    print(f"\n  Long-format rows: {len(age_df):,}")
    print(f"  Unique regulons: {age_df['eRegulon'].nunique()}")
    print(f"  Unique cell types: {age_df['cell_type'].nunique()}")

    # Compute absolute rho
    age_df['abs_age_rho'] = age_df['age_rho'].abs()

    # ── Clean regulon names (same as ipa.py) ──────────────────────────────
    age_df['regulon_clean'] = (age_df['eRegulon']
                               .str.replace(r'_extended', '', regex=False)
                               .str.replace(r'_\(\d+g\)', '', regex=True))

    # ── Load TOPPLE classifications ───────────────────────────────────────
    print(f"\nReading {args.topple} ...")
    topple = pd.read_csv(args.topple)
    print(f"  TOPPLE regulons: {len(topple)}")
    print(f"  Stability classes: {topple['stability_class'].value_counts().to_dict()}")

    # ── Merge on cleaned names ────────────────────────────────────────────
    merged = age_df.merge(
        topple[['regulon', 'mean_RI', 'stability_class', 'rank']],
        left_on='regulon_clean', right_on='regulon', how='inner'
    )
    print(f"\n  Merged rows: {len(merged):,}")
    print(f"  Regulons matched: {merged['eRegulon'].nunique()}")

    # ── Split by stability class ──────────────────────────────────────────
    stab = merged[merged['stability_class'] == 'stabilizer']
    dest = merged[merged['stability_class'] == 'destabilizer']

    print(f"\n  Stabilizer rows:    {len(stab):,}  ({stab['eRegulon'].nunique()} regulons)")
    print(f"  Destabilizer rows:  {len(dest):,}  ({dest['eRegulon'].nunique()} regulons)")

    # ── Mann-Whitney: |age_rho| stabilizers vs destabilizers ──────────────
    stab_vals = stab['abs_age_rho'].dropna().values
    dest_vals = dest['abs_age_rho'].dropna().values
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

    # ── Per-regulon summary ───────────────────────────────────────────────
    reg_summary = merged.groupby(['eRegulon', 'stability_class', 'mean_RI', 'rank']).agg(
        mean_abs_age_rho=('abs_age_rho', 'mean'),
        median_abs_age_rho=('abs_age_rho', 'median'),
        mean_age_rho=('age_rho', 'mean'),
        n_tests=('abs_age_rho', 'count'),
        n_significant=('age_pvalue', lambda x: (x < 0.05).sum()),
    ).reset_index()
    reg_summary = reg_summary.sort_values('rank')

    # ── Save ──────────────────────────────────────────────────────────────
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    reg_summary.to_csv(args.output, index=False)
    print(f"\nSaved per-regulon summary to {args.output}")

    # ── Spearman: RI vs mean |age_rho| per regulon ───────────────────────
    rho_ri, p_ri = stats.spearmanr(reg_summary['mean_RI'], reg_summary['mean_abs_age_rho'])

    # ── Key statistics ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("IPA KEY STATISTICS — Age Perturbation")
    print("=" * 60)
    print(f"  Stabilizer |age_rho|:   mean={stab['abs_age_rho'].mean():.4f}  "
          f"median={stab['abs_age_rho'].median():.4f}")
    print(f"  Destabilizer |age_rho|: mean={dest['abs_age_rho'].mean():.4f}  "
          f"median={dest['abs_age_rho'].median():.4f}")
    print(f"\n  Mann-Whitney U test (stabilizer vs destabilizer):")
    print(f"    U statistic:      {u_stat:.0f}")
    print(f"    P (two-sided):    {p_two:.2e}")
    print(f"    P (stab < dest):  {p_less:.2e}")
    print(f"    P (stab > dest):  {p_greater:.2e}")
    print(f"    Rank-biserial r:  {rank_biserial:.4f}")
    print(f"      (r > 0 means stabilizers have larger |age_rho|)")
    print(f"\n  Spearman (RI vs mean |age_rho| per regulon):")
    print(f"    rho={rho_ri:.4f}, P={p_ri:.2e}")
    print("=" * 60)


if __name__ == '__main__':
    main()
