#!/usr/bin/env python3
"""
Extension #1 remaining analyses:
A. Per-cell-type caQTL vs eQTL non-additivity correlation
B. SICAI metrics for caQTL coupling matrix
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def celltype_correlation(s6_path):
    """Compute per-cell-type mean non-additivity for caQTL and eQTL, then correlate."""
    print("=" * 60)
    print("EXT1 STEP 3c: Per-cell-type non-additivity correlation")
    print("=" * 60)

    s6 = pd.read_csv(s6_path)
    print(f"  S6 shape: {s6.shape}")

    s6_eqtl = s6[s6['analysis'] == 'cis-eQTL'].copy()
    s6_caqtl = s6[s6['analysis'] == 'cis-caQTL'].copy()

    # Compute Gini coefficient (non-additivity proxy) per cell type
    def gini(values):
        v = np.sort(np.abs(values))
        n = len(v)
        if n < 2 or v.sum() == 0:
            return np.nan
        idx = np.arange(1, n + 1)
        return (2 * np.sum(idx * v) - (n + 1) * np.sum(v)) / (n * np.sum(v))

    # For each cell type, gather all slopes and compute Gini
    eqtl_ct_gini = s6_eqtl.groupby('celltype')['slope'].apply(
        lambda x: gini(x.values) if len(x) >= 10 else np.nan
    ).dropna()

    caqtl_ct_gini = s6_caqtl.groupby('celltype')['slope'].apply(
        lambda x: gini(x.values) if len(x) >= 10 else np.nan
    ).dropna()

    print(f"  eQTL cell types with Gini: {len(eqtl_ct_gini)}")
    print(f"  caQTL cell types with Gini: {len(caqtl_ct_gini)}")

    # Alternative: use pre-computed DGSA scores aggregated per cell type
    # Load the raw S6 data and compute per-cell-type mean |slope| and variance
    eqtl_ct_stats = s6_eqtl.groupby('celltype').agg(
        n_genes=('phenotype_id', 'nunique'),
        mean_abs_slope=('slope', lambda x: x.abs().mean()),
        var_slope=('slope', 'var'),
        gini_slope=('slope', lambda x: gini(x.values) if len(x) >= 10 else np.nan),
    ).dropna()

    caqtl_ct_stats = s6_caqtl.groupby('celltype').agg(
        n_peaks=('phenotype_id', 'nunique'),
        mean_abs_slope=('slope', lambda x: x.abs().mean()),
        var_slope=('slope', 'var'),
        gini_slope=('slope', lambda x: gini(x.values) if len(x) >= 10 else np.nan),
    ).dropna()

    # Find shared cell types
    shared_cts = sorted(set(eqtl_ct_stats.index) & set(caqtl_ct_stats.index))
    print(f"  Shared cell types: {len(shared_cts)}")

    # Build correlation table
    results = []
    for ct in shared_cts:
        results.append({
            'cell_type': ct,
            'eqtl_n_genes': eqtl_ct_stats.loc[ct, 'n_genes'],
            'eqtl_mean_abs_slope': eqtl_ct_stats.loc[ct, 'mean_abs_slope'],
            'eqtl_gini': eqtl_ct_stats.loc[ct, 'gini_slope'],
            'caqtl_n_peaks': caqtl_ct_stats.loc[ct, 'n_peaks'],
            'caqtl_mean_abs_slope': caqtl_ct_stats.loc[ct, 'mean_abs_slope'],
            'caqtl_gini': caqtl_ct_stats.loc[ct, 'gini_slope'],
        })

    result_df = pd.DataFrame(results)

    # Correlate eQTL vs caQTL Gini
    valid = result_df.dropna(subset=['eqtl_gini', 'caqtl_gini'])
    rho, p = stats.spearmanr(valid['eqtl_gini'], valid['caqtl_gini'])
    print(f"\n  eQTL vs caQTL Gini correlation (n={len(valid)} CTs):")
    print(f"    Spearman rho = {rho:.3f}, P = {p:.2e}")

    # Mean comparison
    print(f"\n  Mean Gini (non-additivity proxy):")
    print(f"    eQTL:  {valid['eqtl_gini'].mean():.4f} +/- {valid['eqtl_gini'].std():.4f}")
    print(f"    caQTL: {valid['caqtl_gini'].mean():.4f} +/- {valid['caqtl_gini'].std():.4f}")

    w_stat, w_p = stats.wilcoxon(valid['caqtl_gini'], valid['eqtl_gini'])
    print(f"    Paired Wilcoxon: W={w_stat:.0f}, P={w_p:.2e}")
    direction = "caQTL > eQTL" if valid['caqtl_gini'].mean() > valid['eqtl_gini'].mean() else "eQTL > caQTL"
    print(f"    Direction: {direction}")

    # Slope correlation
    rho2, p2 = stats.spearmanr(valid['eqtl_mean_abs_slope'], valid['caqtl_mean_abs_slope'])
    print(f"\n  eQTL vs caQTL mean |slope| correlation:")
    print(f"    Spearman rho = {rho2:.3f}, P = {p2:.2e}")

    result_df.to_csv('results/ext1_celltype_correlation.csv', index=False)
    print(f"\nSaved to results/ext1_celltype_correlation.csv")
    return result_df


def sicai_caqtl_metrics(s8_path, s15_path):
    """Compute SICAI metrics for caQTL r_b matrix and correlate with disease."""
    print("\n" + "=" * 60)
    print("EXT1: SICAI metrics for caQTL coupling")
    print("=" * 60)

    # Load caQTL r_b from S8
    s8_caqtl = pd.read_excel(s8_path, sheet_name='cis_caQTL')
    print(f"  S8 caQTL r_b columns: {s8_caqtl.columns.tolist()}")
    print(f"  S8 caQTL r_b shape: {s8_caqtl.shape}")

    # Also load eQTL r_b for comparison
    s8_eqtl = pd.read_excel(s8_path, sheet_name='cis_eQTL')

    # Build caQTL metrics per cell type
    caqtl_cts = sorted(set(s8_caqtl['reference_cell_type']) | set(s8_caqtl['query_celltype']))
    print(f"  caQTL cell types: {len(caqtl_cts)}")

    caqtl_metrics = []
    for ct in caqtl_cts:
        rb_vals = pd.concat([
            s8_caqtl[s8_caqtl['reference_cell_type'] == ct]['rb'],
            s8_caqtl[s8_caqtl['query_celltype'] == ct]['rb']
        ]).dropna()
        if len(rb_vals) < 3:
            continue

        mean_rb = rb_vals.mean()
        cv = rb_vals.std() / mean_rb if mean_rb > 0 else np.nan

        # Shannon entropy of normalized r_b distribution
        rb_norm = rb_vals / rb_vals.sum()
        entropy = -np.sum(rb_norm * np.log2(rb_norm + 1e-10))

        caqtl_metrics.append({
            'cell_type': ct,
            'caqtl_mean_rb': mean_rb,
            'caqtl_cv': cv,
            'caqtl_entropy': entropy,
            'caqtl_n_pairs': len(rb_vals),
        })

    caqtl_df = pd.DataFrame(caqtl_metrics)

    # eQTL metrics for comparison
    eqtl_cts = sorted(set(s8_eqtl['reference_cell_type']) | set(s8_eqtl['query_celltype']))
    eqtl_metrics = []
    for ct in eqtl_cts:
        rb_vals = pd.concat([
            s8_eqtl[s8_eqtl['reference_cell_type'] == ct]['rb'],
            s8_eqtl[s8_eqtl['query_celltype'] == ct]['rb']
        ]).dropna()
        if len(rb_vals) < 3:
            continue
        mean_rb = rb_vals.mean()
        cv = rb_vals.std() / mean_rb if mean_rb > 0 else np.nan
        rb_norm = rb_vals / rb_vals.sum()
        entropy = -np.sum(rb_norm * np.log2(rb_norm + 1e-10))
        eqtl_metrics.append({
            'cell_type': ct,
            'eqtl_mean_rb': mean_rb,
            'eqtl_cv': cv,
            'eqtl_entropy': entropy,
        })
    eqtl_df = pd.DataFrame(eqtl_metrics)

    # Merge with disease counts from S15
    s15 = pd.read_excel(s15_path)
    disease_counts = s15.groupby('celltype').size().reset_index(name='n_disease')

    merged = caqtl_df.merge(disease_counts, left_on='cell_type', right_on='celltype', how='left')
    merged['n_disease'] = merged['n_disease'].fillna(0)
    merged = merged.merge(eqtl_df, on='cell_type', how='left')

    print(f"\n  caQTL SICAI metrics (n={len(merged)} cell types):")
    print(f"    Mean caQTL r_b: {merged['caqtl_mean_rb'].mean():.3f}")
    print(f"    Mean eQTL r_b:  {merged['eqtl_mean_rb'].mean():.3f}")

    # Correlate each metric with disease count
    valid = merged[merged['n_disease'] > 0]
    print(f"\n  Cell types with disease associations: {len(valid)}")

    for metric, label in [('caqtl_mean_rb', 'caQTL mean r_b'),
                           ('caqtl_cv', 'caQTL CV'),
                           ('caqtl_entropy', 'caQTL entropy'),
                           ('eqtl_mean_rb', 'eQTL mean r_b'),
                           ('eqtl_cv', 'eQTL CV'),
                           ('eqtl_entropy', 'eQTL entropy')]:
        v = merged.dropna(subset=[metric, 'n_disease'])
        if len(v) > 5:
            rho, p = stats.spearmanr(v[metric], v['n_disease'])
            print(f"    {label:20s} vs disease: rho={rho:.3f}, P={p:.2e} (n={len(v)})")

    # Compare caQTL vs eQTL r_b on shared cell types
    shared = merged.dropna(subset=['caqtl_mean_rb', 'eqtl_mean_rb'])
    if len(shared) > 3:
        rho, p = stats.spearmanr(shared['caqtl_mean_rb'], shared['eqtl_mean_rb'])
        print(f"\n  caQTL vs eQTL mean r_b correlation (n={len(shared)}):")
        print(f"    Spearman rho={rho:.3f}, P={p:.2e}")

        w, wp = stats.wilcoxon(shared['caqtl_mean_rb'], shared['eqtl_mean_rb'])
        print(f"    Paired Wilcoxon: W={w:.0f}, P={wp:.2e}")
        print(f"    caQTL mean: {shared['caqtl_mean_rb'].mean():.3f} vs eQTL: {shared['eqtl_mean_rb'].mean():.3f}")

    merged.to_csv('results/ext1_sicai_caqtl_metrics.csv', index=False)
    print(f"\nSaved to results/ext1_sicai_caqtl_metrics.csv")
    return merged


if __name__ == "__main__":
    celltype_correlation('data/raw/CIMA_Table_S6.csv')
    sicai_caqtl_metrics(
        'data/raw/science.adt3130_table_s8.xlsx',
        'data/raw/science.adt3130_table_s15.xlsx'
    )
