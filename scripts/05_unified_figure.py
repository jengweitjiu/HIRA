#!/usr/bin/env python3
"""
HIRA Complete 8-Panel Publication Figure (v1.0)
================================================
Row 1: TOPPLE (A: stabilizer ranking, B: stability landscape)
Row 2: CIMA disease (C: SICAI mean r_b, D: DGSA non-additivity)
Row 3: IPA + cross-method (E: sex perturbation violin, F: correlation matrix)
Row 4: Tissue validation (G: STRATA stabilizer expr, H: coupling heatmap)

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

from utils import (setup_stdout, setup_figure_style, save_figure, scatter_panel,
                   format_sig, check_prerequisites, find_file,
                   COLORS, RESULTS_DIR, FIG_DIR, DATA_DIR)
setup_stdout()


def load_ipa_sex_data():
    """Load sex perturbation data from Table S5 and merge with TOPPLE roles."""
    s5_file = find_file(DATA_DIR, table_key='s5')
    xls = pd.ExcelFile(s5_file)
    sex_sheet = [s for s in xls.sheet_names if 'sex' in s.lower()][0]
    sex_df = pd.read_excel(s5_file, sheet_name=sex_sheet)

    # Detect columns
    cols = sex_df.columns.tolist()
    reg_col = [c for c in cols if 'regulon' in c.lower() or 'eregulon' in c.lower()][0]
    fc_candidates = [c for c in cols if 'log2' in c.lower() or 'fc' in c.lower()]
    fc_col = fc_candidates[0] if fc_candidates else None

    if fc_col is None:
        return None

    sex_df['regulon_raw'] = sex_df[reg_col].astype(str).str.strip()
    sex_df['abs_log2FC'] = sex_df[fc_col].abs()

    # Load TOPPLE roles — match by prefix (e.g., "ADNP_+_(32g)" → "ADNP_+")
    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    role_map = scores['role'].to_dict()

    # Build prefix map: strip gene count suffix like "_(32g)"
    import re
    def clean_regulon(name):
        # "ADNP_+_(32g)" → "ADNP_+"
        return re.sub(r'_\(\d+g\)$', '', name)

    sex_df['regulon_clean'] = sex_df['regulon_raw'].apply(clean_regulon)
    sex_df['role'] = sex_df['regulon_clean'].map(role_map)
    merged = sex_df.dropna(subset=['role', 'abs_log2FC'])

    return merged[['regulon_clean', 'abs_log2FC', 'role']]


def create_complete_figure():
    """Create the definitive 8-panel HIRA figure."""
    setup_figure_style()

    # ── Load all results ──
    scores = pd.read_csv(RESULTS_DIR / "topple_stability_scores.csv", index_col=0)
    ri_matrix = pd.read_csv(RESULTS_DIR / "topple_ri_matrix.csv", index_col=0)
    alt_merged = pd.read_csv(RESULTS_DIR / "sicai_alt_metrics_disease_merged.csv", index_col=0)
    dgsa_ct = pd.read_csv(RESULTS_DIR / "dgsa_celltype_summary.csv", index_col=0)
    strata = pd.read_csv(RESULTS_DIR / "strata_sample_metrics.csv")
    atlas_mapping = pd.read_csv(RESULTS_DIR / "strata_atlas_mapping.csv", index_col=0)
    rb_matrix = pd.read_csv(RESULTS_DIR / "sicai_rb_matrix.csv", index_col=0)

    disease_counts = alt_merged['n_disease_associations']
    ct_mean_ri = ri_matrix.mean(axis=0)
    ct_mean_ri.name = 'mean_RI'

    # ── Figure layout: 4 rows × 2 cols ──
    fig = plt.figure(figsize=(16, 24))
    gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 1],
                           hspace=0.30, wspace=0.28)

    # ==================================================================
    # Panel A: TOPPLE Top 15 Stabilizer Bar Chart
    # ==================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    top15 = scores.head(15)
    colors_a = ['#b71c1c' if r <= 3 else '#c0392b' if r <= 5
                else '#e74c3c' if r <= 10 else '#f39c12'
                for r in range(1, 16)]
    ax_a.barh(range(len(top15)), top15['mean_RI'], color=colors_a,
              edgecolor='white', height=0.7)
    ax_a.set_yticks(range(len(top15)))
    ax_a.set_yticklabels(top15.index, fontsize=7)
    ax_a.set_xlabel('Mean Redistribution Index')
    ax_a.set_title('A  TOPPLE: Top 15 Stabilizer eRegulons',
                    fontweight='bold', loc='left')
    ax_a.invert_yaxis()
    ax_a.get_yticklabels()[0].set_fontweight('bold')
    ax_a.get_yticklabels()[0].set_color('#b71c1c')

    # STAT5B annotation
    stat5b_hits = [r for r in scores.index if r.upper().startswith('STAT5B')]
    if stat5b_hits:
        s5key = stat5b_hits[0]
        s5rank = int(scores.loc[s5key, 'stability_rank'])
        s5ri = scores.loc[s5key, 'mean_RI']
        ax_a.axvline(s5ri, color='#7f8c8d', linestyle=':', linewidth=1, alpha=0.7)
        ax_a.annotate(
            f'{s5key} (#{s5rank})\n#1 in psoriasis',
            xy=(s5ri, 14), fontsize=6.5, color='#7f8c8d', fontstyle='italic',
            xytext=(s5ri + (top15['mean_RI'].max() - s5ri) * 0.3, 13),
            arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=0.8),
            va='center')

    # ==================================================================
    # Panel B: TOPPLE Stability Landscape Heatmap
    # ==================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    top25 = scores.head(25).index.tolist()
    bot25 = scores.tail(25).index.tolist()
    ri_sel = ri_matrix.loc[top25 + bot25]
    ri_z = ri_sel.sub(ri_sel.mean(axis=1), axis=0).div(
        ri_sel.std(axis=1) + 1e-12, axis=0).fillna(0)

    sns.heatmap(ri_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'RI (z-score)', 'shrink': 0.5}, ax=ax_b)
    ax_b.set_title('B  Stability Landscape: Stabilizers vs Dispensable',
                    fontweight='bold', loc='left', fontsize=9)
    ax_b.set_xlabel('Cell Types (61 CIMA L4)')
    ax_b.set_ylabel('')
    ax_b.axhline(y=25, color='black', linewidth=1.5)
    top5_set = set(scores.head(5).index)
    for lbl in ax_b.get_yticklabels():
        if lbl.get_text() in top5_set:
            lbl.set_color('#c0392b')
            lbl.set_fontweight('bold')

    # ==================================================================
    # Panel C: SICAI Mean r_b vs Disease
    # ==================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    scatter_panel(
        ax_c,
        alt_merged['mean_rb'], alt_merged['n_disease_associations'],
        'Mean $r_b$ (Coupling Strength)',
        'Disease Associations (SMR)',
        'C  SICAI: Coupling Strength vs Disease Pleiotropy')

    # ==================================================================
    # Panel D: DGSA Non-additivity vs Disease
    # ==================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    dgsa_disease = dgsa_ct.join(disease_counts, how='inner')
    scatter_panel(
        ax_d,
        dgsa_disease['mean_non_additivity'], dgsa_disease['n_disease_associations'],
        'Mean Non-Additivity Index',
        'Disease Associations (SMR)',
        'D  DGSA: Genetic Non-Additivity vs Disease Pleiotropy')

    # ==================================================================
    # Panel E: IPA Sex Perturbation Violin
    # ==================================================================
    ax_e = fig.add_subplot(gs[2, 0])
    sex_data = load_ipa_sex_data()
    if sex_data is not None:
        sex_data = sex_data.dropna(subset=['abs_log2FC'])
        role_palette = {
            'stabilizer': COLORS['stabilizer'],
            'neutral': COLORS['neutral'],
            'destabilizer': COLORS['destabilizer'],
        }
        role_order = ['stabilizer', 'neutral', 'destabilizer']
        sns.violinplot(data=sex_data, x='role', y='abs_log2FC',
                       order=role_order, palette=role_palette,
                       inner='quartile', cut=0, ax=ax_e)
        ax_e.set_xlabel('')
        ax_e.set_ylabel('|log$_2$ Fold Change| (Sex Difference)')
        ax_e.set_title('E  IPA: Stabilizers Resist Sex Perturbation',
                        fontweight='bold', loc='left', fontsize=9)
        ax_e.set_ylim(0, sex_data['abs_log2FC'].quantile(0.99))

        # Add P-value annotation
        stab_v = sex_data[sex_data['role'] == 'stabilizer']['abs_log2FC']
        dest_v = sex_data[sex_data['role'] == 'destabilizer']['abs_log2FC']
        _, p_sex = mannwhitneyu(stab_v, dest_v, alternative='less')
        ymax_e = ax_e.get_ylim()[1]
        ax_e.plot([0, 0, 2, 2],
                  [ymax_e * 0.88, ymax_e * 0.91, ymax_e * 0.91, ymax_e * 0.88],
                  color='black', linewidth=1)
        ax_e.text(1, ymax_e * 0.93, f'P = {p_sex:.1e}',
                  ha='center', fontsize=8, fontweight='bold')
        # Median annotations
        ax_e.text(0, stab_v.median() + ymax_e * 0.02,
                  f'med={stab_v.median():.3f}', ha='center', fontsize=6.5,
                  color=COLORS['stabilizer'])
        ax_e.text(2, dest_v.median() + ymax_e * 0.02,
                  f'med={dest_v.median():.3f}', ha='center', fontsize=6.5,
                  color=COLORS['destabilizer'])

    # ==================================================================
    # Panel F: Cross-Method Correlation Matrix
    # ==================================================================
    ax_f = fig.add_subplot(gs[2, 1])

    topple_series = ct_mean_ri.rename('TOPPLE_RI')
    sicai_series = alt_merged['mean_rb'].rename('SICAI_mean_rb')
    dgsa_series = dgsa_ct['mean_non_additivity'].rename('DGSA_NonAdd')
    disease_series = disease_counts.rename('Disease_N')
    all_metrics = pd.concat([topple_series, sicai_series, dgsa_series,
                             disease_series], axis=1).dropna()

    labels = ['TOPPLE RI', 'SICAI\nmean $r_b$', 'DGSA\nNon-Add', 'Disease\nCount']
    n = len(labels)
    rho_mat = np.ones((n, n))
    p_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                r, p = spearmanr(all_metrics.iloc[:, i], all_metrics.iloc[:, j])
                rho_mat[i, j] = r
                p_mat[i, j] = p

    annot = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            if i == j:
                annot[i, j] = '1.000'
            else:
                sig = format_sig(p_mat[i, j])
                sig = '' if sig == 'n.s.' else sig
                annot[i, j] = f'{rho_mat[i, j]:.3f}{sig}'

    sns.heatmap(pd.DataFrame(rho_mat, index=labels, columns=labels),
                cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                annot=annot, fmt='', annot_kws={'fontsize': 9, 'fontweight': 'bold'},
                cbar_kws={'label': 'Spearman $\\rho$', 'shrink': 0.7},
                linewidths=1.5, linecolor='white', ax=ax_f)
    ax_f.set_title('F  Cross-Method Correlation Matrix (CIMA)',
                    fontweight='bold', loc='left')
    ax_f.set_xticklabels(labels, fontsize=8, rotation=0, ha='center')
    ax_f.set_yticklabels(labels, fontsize=8, rotation=0)
    ax_f.text(0.5, -0.08, f'n = {len(all_metrics)} cell types',
              transform=ax_f.transAxes, fontsize=7, ha='center', color='#7f8c8d')

    # ==================================================================
    # Panel G: STRATA — Stabilizer expression by condition
    # ==================================================================
    ax_g = fig.add_subplot(gs[3, 0])
    cond_palette = {'Healthy': '#27ae60', 'Non-Lesional': '#3498db',
                    'Lesional': '#e74c3c'}
    order = ['Healthy', 'Non-Lesional', 'Lesional']
    avail_order = [c for c in order if c in strata['condition'].values]

    sns.boxplot(data=strata, x='condition', y='mean_stabilizer_expr',
                order=avail_order, palette=cond_palette,
                width=0.5, fliersize=0, ax=ax_g)
    sns.stripplot(data=strata, x='condition', y='mean_stabilizer_expr',
                  order=avail_order, palette=cond_palette,
                  size=7, alpha=0.7, jitter=0.15, ax=ax_g)
    ax_g.set_xlabel('')
    ax_g.set_ylabel('Mean Stabilizer TF Expression')
    ax_g.set_title('G  STRATA: Stabilizer Expression (GSE202011 Visium)',
                    fontweight='bold', loc='left', fontsize=9)

    h_vals = strata[strata['condition'] == 'Healthy']['mean_stabilizer_expr']
    l_vals = strata[strata['condition'] == 'Lesional']['mean_stabilizer_expr']
    if len(h_vals) > 0 and len(l_vals) > 0:
        _, p_hl = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
        ymax_g = strata['mean_stabilizer_expr'].max()
        h_idx = avail_order.index('Healthy')
        l_idx = avail_order.index('Lesional')
        ax_g.plot([h_idx, h_idx, l_idx, l_idx],
                  [ymax_g * 1.05, ymax_g * 1.08, ymax_g * 1.08, ymax_g * 1.05],
                  color='black', linewidth=1)
        ax_g.text((h_idx + l_idx) / 2, ymax_g * 1.09,
                  f'P = {p_hl:.2e}', ha='center', fontsize=8)

    # ==================================================================
    # Panel H: STRATA-Atlas — Coupling comparison heatmap
    # ==================================================================
    ax_h = fig.add_subplot(gs[3, 1])

    # Import the mapping from 07_strata_atlas
    from importlib.util import spec_from_file_location, module_from_spec
    import os
    spec = spec_from_file_location("strata_atlas",
        os.path.join(os.path.dirname(__file__), "07_strata_atlas.py"))
    sa_mod = module_from_spec(spec)
    spec.loader.exec_module(sa_mod)

    cima_metrics, _ = sa_mod.load_cima_metrics()
    mapped = sa_mod.aggregate_cima_to_gse(cima_metrics)

    # Load GSE immune scores and compute coupling matrices
    sample_means = sa_mod.load_gse_immune_scores()
    cima_coupling, gse_coupling = sa_mod.compute_coupling_matrices(
        sample_means, rb_matrix, mapped)

    shared_cats = cima_coupling.index.intersection(gse_coupling.index)
    if len(shared_cats) >= 3:
        n_sc = len(shared_cats)
        cima_vals = cima_coupling.loc[shared_cats, shared_cats].values.astype(float)
        gse_vals = gse_coupling.loc[shared_cats, shared_cats].values.astype(float)

        combined = np.zeros((n_sc, n_sc))
        for i in range(n_sc):
            for j in range(n_sc):
                if i == j:
                    combined[i, j] = 1.0
                elif i < j:
                    combined[i, j] = cima_vals[i, j] if not np.isnan(cima_vals[i, j]) else 0
                else:
                    combined[i, j] = gse_vals[i, j] if not np.isnan(gse_vals[i, j]) else 0

        short_names = [c.replace('Macrophage_', 'Mac_') for c in shared_cats]
        comb_df = pd.DataFrame(combined, index=short_names, columns=short_names)

        sns.heatmap(comb_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={'label': 'Coupling', 'shrink': 0.6},
                    linewidths=0.5, linecolor='white', ax=ax_h)
        ax_h.set_title('H  Blood Genetic (upper) vs Skin Spatial (lower)',
                        fontweight='bold', loc='left', fontsize=9)
        ax_h.tick_params(labelsize=7)

        # Mantel test
        mantel_r, mantel_p, mantel_n = sa_mod.mantel_test(
            cima_coupling.loc[shared_cats, shared_cats],
            gse_coupling.loc[shared_cats, shared_cats], n_perm=9999)
        ax_h.text(0.5, -0.10,
                  f'Mantel r = {mantel_r:.3f}, P = {mantel_p:.3f} '
                  f'({mantel_n} pairs)',
                  transform=ax_h.transAxes, fontsize=8, ha='center',
                  color='#2c3e50', fontweight='bold')
        ax_h.text(0.82, 0.18, 'CIMA\n$r_b$', transform=ax_h.transAxes,
                  fontsize=7, ha='center', color='#7f8c8d', fontstyle='italic')
        ax_h.text(0.18, 0.82, 'GSE\nspatial', transform=ax_h.transAxes,
                  fontsize=7, ha='center', color='#7f8c8d', fontstyle='italic')

    # ── Suptitle ──
    fig.suptitle(
        'HIRA: Hierarchical Immune Regulatory Architecture\n'
        '6-Layer Framework Validated on CIMA Atlas + GSE202011 Visium',
        fontsize=14, fontweight='bold', y=0.995)

    save_figure(fig, FIG_DIR, 'fig_HIRA_complete')
    print("Complete figure saved: fig_HIRA_complete.png/pdf")

    return all_metrics, rho_mat, p_mat


def main():
    print("=" * 70)
    print("HIRA Complete 8-Panel Figure (v1.0)")
    print("=" * 70)

    required = [
        "topple_stability_scores.csv",
        "topple_ri_matrix.csv",
        "sicai_alt_metrics_disease_merged.csv",
        "sicai_rb_matrix.csv",
        "dgsa_celltype_summary.csv",
        "strata_sample_metrics.csv",
        "strata_atlas_mapping.csv",
    ]
    check_prerequisites(RESULTS_DIR, required)

    create_complete_figure()

    print(f"\n{'=' * 70}")
    print("DONE — fig_HIRA_complete.pdf ready for submission")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
