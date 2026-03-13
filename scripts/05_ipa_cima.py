#!/usr/bin/env python3
"""
IPA Perturbation Cascade Analysis on CIMA Table S5 (Sheets 4-5)
================================================================
Tests whether TOPPLE stabilizers resist natural perturbation (sex, age).

Core hypothesis: structurally critical regulons (high RI) should be
BUFFERED against demographic variation — lower |log2FC| for sex,
lower |age_corr| for age.

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import re
import warnings
warnings.filterwarnings('ignore')

from utils import (setup_stdout, ensure_dirs, setup_figure_style, save_figure,
                   find_file, scatter_panel, format_sig,
                   COLORS, DATA_DIR, FIG_DIR, RESULTS_DIR)
setup_stdout()
ensure_dirs()

S5_FILE = find_file(DATA_DIR, table_key='s5')


def strip_gene_count(name: str) -> str:
    """Remove gene-count suffix like '_(32g)' from eRegulon names."""
    return re.sub(r'_\(\d+g\)$', '', str(name).strip())


# ── Step 1: Load data ─────────────────────────────────────────────────
def load_sex_diff() -> pd.DataFrame:
    """Load sheet 4: sex_difference_results (21,289 rows)."""
    print("Loading sex_difference_results...")
    df = pd.read_excel(S5_FILE, sheet_name='sex_difference_results')
    print(f"  Shape: {df.shape}")
    # Clean regulon names
    df['eRegulon_clean'] = df['eRegulon'].apply(strip_gene_count)
    print(f"  Unique regulons (raw): {df['eRegulon'].nunique()}")
    print(f"  Unique regulons (clean): {df['eRegulon_clean'].nunique()}")
    print(f"  Unique cell types: {df['cell_type_l4'].nunique()}")
    return df


def load_age_corr() -> pd.DataFrame:
    """Load sheet 5: Age_Correlation (wide format) and reshape to long."""
    print("\nLoading Age_Correlation (wide -> long)...")
    df_wide = pd.read_excel(S5_FILE, sheet_name='Age_Correlation')
    print(f"  Wide shape: {df_wide.shape}")

    # First column is regulon name
    regulon_col = df_wide.columns[0]  # 'Unnamed: 0'
    df_wide = df_wide.rename(columns={regulon_col: 'eRegulon'})

    # Identify cell type columns: those WITHOUT '_Pvalue' suffix
    corr_cols = [c for c in df_wide.columns[1:] if not c.endswith('_Pvalue')]
    pval_cols = [c for c in df_wide.columns[1:] if c.endswith('_Pvalue')]

    rows = []
    for _, row in df_wide.iterrows():
        reg = row['eRegulon']
        for ct in corr_cols:
            pval_col = ct + '_Pvalue'
            corr_val = row[ct]
            pval_val = row.get(pval_col, np.nan)
            if pd.notna(corr_val):
                rows.append({
                    'eRegulon': reg,
                    'cell_type_l4': ct,
                    'age_corr': corr_val,
                    'age_pval': pval_val,
                })

    df_long = pd.DataFrame(rows)
    df_long['eRegulon_clean'] = df_long['eRegulon'].apply(strip_gene_count)
    print(f"  Long shape: {df_long.shape}")
    print(f"  Unique regulons: {df_long['eRegulon_clean'].nunique()}")
    print(f"  Unique cell types: {df_long['cell_type_l4'].nunique()}")
    return df_long


def load_topple_scores() -> pd.DataFrame:
    """Load TOPPLE stability scores."""
    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    print(f"\nTOPPLE scores: {len(scores)} regulons")
    print(f"  Roles: {scores['role'].value_counts().to_dict()}")
    return scores


# ── Step 2: Merge ─────────────────────────────────────────────────────
def merge_with_topple(sex_df, age_df, scores):
    """Merge sex/age data with TOPPLE stability scores."""
    # Build lookup from TOPPLE scores
    topple_lookup = scores[['mean_RI', 'role', 'stability_rank']].copy()
    topple_lookup.index.name = 'eRegulon_clean'

    # Merge sex
    sex_merged = sex_df.merge(
        topple_lookup, left_on='eRegulon_clean', right_index=True, how='inner')
    sex_merged['abs_log2FC'] = sex_merged['Log2_Fold_Change'].abs()

    # Merge age
    age_merged = age_df.merge(
        topple_lookup, left_on='eRegulon_clean', right_index=True, how='inner')
    age_merged['abs_age_corr'] = age_merged['age_corr'].abs()

    print(f"\n  Sex merged: {len(sex_merged)} rows "
          f"({sex_merged['eRegulon_clean'].nunique()} regulons)")
    print(f"  Age merged: {len(age_merged)} rows "
          f"({age_merged['eRegulon_clean'].nunique()} regulons)")

    return sex_merged, age_merged


# ── Step 3: Core tests ────────────────────────────────────────────────
def core_tests(sex_merged, age_merged):
    """Stabilizers vs destabilizers: Mann-Whitney U tests."""
    results = {}

    print(f"\n{'=' * 70}")
    print("CORE TEST: Do stabilizers resist natural perturbation?")
    print(f"{'=' * 70}")

    for label, df, val_col, sig_col, sig_thresh in [
        ('Sex (|log2FC|)', sex_merged, 'abs_log2FC', 'adjust_P_Value', 0.05),
        ('Age (|corr|)', age_merged, 'abs_age_corr', 'age_pval', 0.05),
    ]:
        stab = df[df['role'] == 'stabilizer'][val_col].dropna()
        dest = df[df['role'] == 'destabilizer'][val_col].dropna()
        neut = df[df['role'] == 'neutral'][val_col].dropna()

        U, p = mannwhitneyu(stab, dest, alternative='less')

        # Significance fractions
        stab_sig = (df[df['role'] == 'stabilizer'][sig_col] < sig_thresh).mean()
        dest_sig = (df[df['role'] == 'destabilizer'][sig_col] < sig_thresh).mean()
        neut_sig = (df[df['role'] == 'neutral'][sig_col] < sig_thresh).mean()

        sig_marker = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
        print(f"\n  {label}:")
        print(f"    Stabilizers:   median = {stab.median():.4f} (n={len(stab)})")
        print(f"    Neutral:       median = {neut.median():.4f} (n={len(neut)})")
        print(f"    Destabilizers: median = {dest.median():.4f} (n={len(dest)})")
        print(f"    Mann-Whitney U (stab < dest): U={U:,.0f}, P={p:.2e} {sig_marker}")
        print(f"    Frac significant: stab={stab_sig:.3f}, neut={neut_sig:.3f}, dest={dest_sig:.3f}")

        results[label] = {
            'stab_median': stab.median(), 'dest_median': dest.median(),
            'neut_median': neut.median(),
            'U': U, 'p': p,
            'stab_frac_sig': stab_sig, 'dest_frac_sig': dest_sig,
        }

    return results


# ── Step 4: Perturbation sensitivity per cell type ────────────────────
def perturbation_sensitivity(sex_merged, age_merged):
    """Mean |log2FC| + mean |age_corr| per cell type."""
    sex_ct = sex_merged.groupby('cell_type_l4')['abs_log2FC'].mean()
    sex_ct.name = 'mean_abs_log2FC'

    age_ct = age_merged.groupby('cell_type_l4')['abs_age_corr'].mean()
    age_ct.name = 'mean_abs_age_corr'

    sens = pd.concat([sex_ct, age_ct], axis=1).dropna()
    sens['perturbation_sensitivity'] = sens['mean_abs_log2FC'] + sens['mean_abs_age_corr']
    sens = sens.sort_values('perturbation_sensitivity', ascending=False)

    print(f"\n  Perturbation sensitivity computed for {len(sens)} cell types")
    return sens


def cross_correlations(sens):
    """Correlate perturbation sensitivity with disease, TOPPLE, DGSA."""
    results = {}

    print(f"\n{'=' * 70}")
    print("PERTURBATION SENSITIVITY vs OTHER METHODS")
    print(f"{'=' * 70}")

    # Disease
    alt_path = RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv"
    if alt_path.exists():
        alt = pd.read_csv(alt_path, index_col=0)
        merged = sens.join(alt['n_disease_associations'], how='inner')
        if len(merged) >= 5:
            rho, p = spearmanr(merged['perturbation_sensitivity'],
                               merged['n_disease_associations'])
            results['disease'] = {'rho': rho, 'p': p, 'n': len(merged), 'merged': merged}
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
            print(f"  vs Disease count: rho={rho:.3f}, P={p:.2e}, n={len(merged)} {sig}")

    # TOPPLE mean RI per cell type
    ri_path = RESULTS_DIR / "topple_ri_matrix.csv"
    if ri_path.exists():
        ri = pd.read_csv(ri_path, index_col=0)
        ct_ri = ri.mean(axis=0)
        ct_ri.name = 'mean_RI'
        merged = sens.join(ct_ri, how='inner')
        if len(merged) >= 5:
            rho, p = spearmanr(merged['perturbation_sensitivity'], merged['mean_RI'])
            results['topple'] = {'rho': rho, 'p': p, 'n': len(merged), 'merged': merged}
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
            print(f"  vs TOPPLE RI:     rho={rho:.3f}, P={p:.2e}, n={len(merged)} {sig}")

    # DGSA
    dgsa_path = RESULTS_DIR / "dgsa_celltype_summary.csv"
    if dgsa_path.exists():
        dgsa = pd.read_csv(dgsa_path, index_col=0)
        merged = sens.join(dgsa['mean_non_additivity'], how='inner')
        if len(merged) >= 5:
            rho, p = spearmanr(merged['perturbation_sensitivity'],
                               merged['mean_non_additivity'])
            results['dgsa'] = {'rho': rho, 'p': p, 'n': len(merged), 'merged': merged}
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
            print(f"  vs DGSA NonAdd:   rho={rho:.3f}, P={p:.2e}, n={len(merged)} {sig}")

    return results


# ── Step 5: Figures ───────────────────────────────────────────────────
def plot_figures(sex_merged, age_merged, sens, core_results, cross_results, scores):
    setup_figure_style()

    role_order = ['stabilizer', 'neutral', 'destabilizer']
    palette = {
        'stabilizer': COLORS['stabilizer'],
        'neutral': COLORS['neutral'],
        'destabilizer': COLORS['destabilizer'],
    }

    # ── Fig 4a: Violin — |log2FC| by role ──
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.violinplot(data=sex_merged, x='role', y='abs_log2FC', order=role_order,
                   palette=palette, inner='box', cut=0, ax=ax)
    r = core_results['Sex (|log2FC|)']
    sig = format_sig(r['p'])
    ax.set_title(f'A  Sex Perturbation Resistance by Stability Role\n'
                 f'Mann-Whitney U (stab < dest): P = {r["p"]:.2e} {sig}',
                 fontweight='bold', loc='left', fontsize=9)
    ax.set_xlabel('TOPPLE Stability Role')
    ax.set_ylabel('|log2 Fold Change| (sex difference)')
    plt.tight_layout()
    save_figure(fig, FIG_DIR, 'fig4a_sex_violin')
    print("  fig4a saved")

    # ── Fig 4b: Violin — |age_corr| by role ──
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.violinplot(data=age_merged, x='role', y='abs_age_corr', order=role_order,
                   palette=palette, inner='box', cut=0, ax=ax)
    r = core_results['Age (|corr|)']
    sig = format_sig(r['p'])
    ax.set_title(f'B  Age Perturbation Resistance by Stability Role\n'
                 f'Mann-Whitney U (stab < dest): P = {r["p"]:.2e} {sig}',
                 fontweight='bold', loc='left', fontsize=9)
    ax.set_xlabel('TOPPLE Stability Role')
    ax.set_ylabel('|Age Correlation| (Spearman)')
    plt.tight_layout()
    save_figure(fig, FIG_DIR, 'fig4b_age_violin')
    print("  fig4b saved")

    # ── Fig 4c: Scatter — perturbation sensitivity vs disease ──
    if 'disease' in cross_results:
        info = cross_results['disease']
        merged = info['merged']
        fig, ax = plt.subplots(figsize=(8, 6))
        scatter_panel(
            ax,
            merged['perturbation_sensitivity'],
            merged['n_disease_associations'],
            xlabel='Perturbation Sensitivity (|log2FC| + |age_corr|)',
            ylabel='Disease Associations (SMR)',
            title='C  Perturbation Sensitivity Predicts Disease Pleiotropy',
            label_top=5,
        )
        plt.tight_layout()
        save_figure(fig, FIG_DIR, 'fig4c_sensitivity_vs_disease')
        print("  fig4c saved")

    # ── Fig 4d: Heatmap — top 10 stabilizers × perturbation across cell types ──
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    top10 = scores.head(10).index.tolist()

    # Sex: top10 × cell types
    sex_pivot = sex_merged[sex_merged['eRegulon_clean'].isin(top10)].pivot_table(
        index='eRegulon_clean', columns='cell_type_l4', values='abs_log2FC', aggfunc='mean')
    # Reorder rows by TOPPLE rank
    sex_pivot = sex_pivot.reindex([r for r in top10 if r in sex_pivot.index])

    if not sex_pivot.empty:
        sns.heatmap(sex_pivot, cmap='YlOrRd', xticklabels=True, yticklabels=True,
                    cbar_kws={'label': '|log2FC|', 'shrink': 0.6}, ax=ax1)
        ax1.set_title('Top 10 Stabilizers: Sex Effect', fontweight='bold', fontsize=9)
        ax1.set_xlabel('Cell Types')
        ax1.set_ylabel('')

    # Age: top10 × cell types
    age_pivot = age_merged[age_merged['eRegulon_clean'].isin(top10)].pivot_table(
        index='eRegulon_clean', columns='cell_type_l4', values='abs_age_corr', aggfunc='mean')
    age_pivot = age_pivot.reindex([r for r in top10 if r in age_pivot.index])

    if not age_pivot.empty:
        sns.heatmap(age_pivot, cmap='YlOrRd', xticklabels=True, yticklabels=True,
                    cbar_kws={'label': '|age corr|', 'shrink': 0.6}, ax=ax2)
        ax2.set_title('Top 10 Stabilizers: Age Effect', fontweight='bold', fontsize=9)
        ax2.set_xlabel('Cell Types')
        ax2.set_ylabel('')

    plt.suptitle('D  Perturbation Response of Top 10 TOPPLE Stabilizers',
                 fontweight='bold', y=1.02)
    plt.tight_layout()
    save_figure(fig, FIG_DIR, 'fig4d_stabilizer_perturbation_heatmap')
    print("  fig4d saved")


# ── Main ───────────────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("IPA Perturbation Cascade Analysis on CIMA Atlas")
    print("Testing: Do TOPPLE stabilizers resist sex/age perturbation?")
    print("=" * 70)

    # 1. Load data
    sex_df = load_sex_diff()
    age_df = load_age_corr()
    scores = load_topple_scores()

    # 2. Merge
    sex_merged, age_merged = merge_with_topple(sex_df, age_df, scores)

    # 3. Core tests
    core_results = core_tests(sex_merged, age_merged)

    # 4. Perturbation sensitivity per cell type
    sens = perturbation_sensitivity(sex_merged, age_merged)
    cross_results = cross_correlations(sens)

    # 5. Figures
    print("\nGenerating figures...")
    plot_figures(sex_merged, age_merged, sens, core_results, cross_results, scores)

    # 6. Save
    sens.to_csv(RESULTS_DIR / "ipa_perturbation_sensitivity.csv")

    summary_rows = []
    for label, r in core_results.items():
        summary_rows.append({
            'test': label,
            'stabilizer_median': r['stab_median'],
            'destabilizer_median': r['dest_median'],
            'neutral_median': r['neut_median'],
            'U': r['U'],
            'p_value': r['p'],
            'stab_frac_sig': r['stab_frac_sig'],
            'dest_frac_sig': r['dest_frac_sig'],
        })
    pd.DataFrame(summary_rows).to_csv(RESULTS_DIR / "ipa_core_tests.csv", index=False)

    print(f"\nResults saved to {RESULTS_DIR}/")

    # 7. Summary
    print(f"\n{'=' * 70}")
    print("IPA ANALYSIS COMPLETE")
    print(f"{'=' * 70}")
    for label, r in core_results.items():
        direction = "CONFIRMED" if r['p'] < 0.05 else "not significant"
        print(f"  {label}: stab median={r['stab_median']:.4f} vs "
              f"dest={r['dest_median']:.4f}, P={r['p']:.2e} -> {direction}")
    for label, info in cross_results.items():
        print(f"  Sensitivity vs {label}: rho={info['rho']:.3f}, P={info['p']:.2e}")


if __name__ == "__main__":
    main()
