#!/usr/bin/env python3
"""
Lineage-Stratified TOPPLE — Fractal Stability Architecture Test.

Split 203 regulons × 61 cell types into 4 lineage groups and run TOPPLE
independently within each. Tests whether stability architecture is fractal:
does the same hierarchical pattern repeat at every scale of immune organization?

Input: CIMA Table S5 (eRegulon AUC matrix)
Output: results/lineage_stratified_topple.csv
        results/lineage_stratified_summary.csv
        figures/fig_lineage_stratified_topple.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


# ---- Lineage assignments for all 61 CIMA cell types ----
LINEAGE_MAP = {
    # T cell lineage (25 cell types)
    'CD4_CTL_GZMH': 'T_cell',
    'CD4_Tcm_CXCR5': 'T_cell',
    'CD4_Tcm_IFI44L': 'T_cell',
    'CD4_Tem_CCR5': 'T_cell',
    'CD4_Tem_CCR7neg': 'T_cell',
    'CD4_Tfh-like_CXCR5': 'T_cell',
    'CD4_Th1-like_GZMK': 'T_cell',
    'CD4_Th17-like_RORC': 'T_cell',
    'CD4_Th22-like_CCR10': 'T_cell',
    'CD4_Th_CCR4': 'T_cell',
    'CD4_Th_LMNA': 'T_cell',
    'CD4_Th_TNFRSF11A': 'T_cell',
    'CD4_Tn_CCR7': 'T_cell',
    'CD4_Tr1_IL10': 'T_cell',
    'CD4_Treg_FCRL3': 'T_cell',
    'CD4_Treg_FOXP3': 'T_cell',
    'CD8_CTL_GZMB': 'T_cell',
    'CD8_CTL_GZMK': 'T_cell',
    'CD8_CTL_IFI44L': 'T_cell',
    'CD8_Tcm_IFI44L': 'T_cell',
    'CD8_Tem_CCR7neg': 'T_cell',
    'CD8_Tn_CCR7': 'T_cell',
    'MAIT_SLC4A10': 'T_cell',
    'gdT1_TRDV1': 'T_cell',
    'gdT2_GZMH': 'T_cell',
    'gdT2_GZMK': 'T_cell',
    'gdT2_IL12RB2': 'T_cell',
    'Proliferative_T_MKI67': 'T_cell',

    # B cell lineage (10 cell types)
    'Bn_IFIT3': 'B_cell',
    'Bn_TCL1A': 'B_cell',
    'Atypical_Bm_ITGAX': 'B_cell',
    'Switched_Bm_IGHDneg': 'B_cell',
    'Switched_Bm_IGHE': 'B_cell',
    'Switched_activated_Bm_CD86': 'B_cell',
    'Unswitched_Bm_CD1C': 'B_cell',
    'Transitional_B_MME': 'B_cell',
    'pre-Switched_Bm_JAM3': 'B_cell',
    'Plasma_IGHA1': 'B_cell',
    'Plasma_IGHG1': 'B_cell',
    'Plasmablast_MKI67': 'B_cell',

    # Myeloid lineage (14 cell types)
    'cMono_CD14': 'Myeloid',
    'cMono_CXCL10': 'Myeloid',
    'cMono_IFI44L': 'Myeloid',
    'cMono_IL1B': 'Myeloid',
    'intMono_GFRA2': 'Myeloid',
    'ncMono_C1QA': 'Myeloid',
    'ncMono_FCGR3A': 'Myeloid',
    'ncMono_IFIT1': 'Myeloid',
    'DC1_CLEC9A': 'Myeloid',
    'DC2_CD1C': 'Myeloid',
    'DC_CSF2RA': 'Myeloid',
    'AS_DC': 'Myeloid',
    'pDC_IRF4': 'Myeloid',

    # NK/ILC lineage (7 cell types)
    'NK_bright_XCL1': 'NK_ILC',
    'Mature_NK_dim_FCGR3A': 'NK_ILC',
    'Terminal_NK_dim_CD160neg': 'NK_ILC',
    'Transitional_NK_GZMK': 'NK_ILC',
    'Proliferative_NK_MKI67': 'NK_ILC',
    'NKT_NCR1': 'NK_ILC',
    'ILC2_IL2RA': 'NK_ILC',

    # Progenitor (excluded from lineage analysis)
    'HSPC_CD34': 'Progenitor',
}


def run_topple(auc_matrix, regulon_names, cell_types, label=""):
    """
    Run TOPPLE on a subset of the AUC matrix.

    Parameters:
        auc_matrix: DataFrame (regulons × cell types) of mean AUC values
        regulon_names: list of regulon names
        cell_types: list of cell type names
        label: string label for this run

    Returns:
        DataFrame with per-regulon RI scores
    """
    n_regs = len(regulon_names)
    n_cts = len(cell_types)

    # Normalize each regulon to a probability distribution across cell types
    mat = auc_matrix.values.astype(float)
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # avoid division by zero
    dist = mat / row_sums

    results = []
    for i in range(n_regs):
        full_dist = dist[i, :]
        if full_dist.sum() == 0:
            continue

        jsd_values = []
        for j in range(n_cts):
            # Remove cell type j
            reduced = np.delete(full_dist, j)
            if reduced.sum() > 0:
                reduced_norm = reduced / reduced.sum()
                # Full distribution without cell type j (renormalized)
                full_without = np.delete(full_dist, j)
                full_without_norm = full_without / full_without.sum() if full_without.sum() > 0 else full_without

                # Compare original distribution (with j zeroed) to reduced
                modified = full_dist.copy()
                modified[j] = 0
                if modified.sum() > 0:
                    modified_norm = modified / modified.sum()
                    jsd = jensenshannon(full_dist, modified_norm) ** 2
                    jsd_values.append(jsd)

        if jsd_values:
            ri = np.mean(jsd_values)
            max_auc = float(auc_matrix.iloc[i].max())
            mean_auc = float(auc_matrix.iloc[i].mean())
            n_active = int((auc_matrix.iloc[i] > 0.01).sum())

            results.append({
                'regulon': regulon_names[i],
                'mean_RI': ri,
                'max_RI': np.max(jsd_values),
                'std_RI': np.std(jsd_values),
                'mean_AUC': mean_auc,
                'max_AUC': max_auc,
                'n_active_cell_types': n_active,
            })

    results_df = pd.DataFrame(results).sort_values('mean_RI', ascending=False)
    results_df['rank'] = range(1, len(results_df) + 1)
    results_df['lineage'] = label
    return results_df


def classify_regulon_type(regulon_name, lineage_results):
    """Classify a regulon as universal or lineage-specific based on its ranking."""
    # Universal stress-response TFs (from CIMA TOPPLE)
    universal_tfs = {
        'HSF1', 'EGR1', 'KLF9', 'JUN', 'JUNB', 'JUND', 'FOS', 'FOSB', 'FOSL1', 'FOSL2',
        'NFKB1', 'NFKB2', 'RELA', 'REL',
        'STAT1', 'STAT3',
        'NR4A1', 'NR4A2', 'NR4A3',
        'NFE2L2', 'BACH1',
        'ATF3', 'ATF4',
        'KLF2', 'KLF4', 'KLF6',
    }
    tf = regulon_name.replace('_+', '').replace('_-', '')
    return 'universal' if tf in universal_tfs else 'lineage-specific'


def main():
    print("=" * 60)
    print("Lineage-Stratified TOPPLE — Fractal Stability Test")
    print("=" * 60)

    # ---- Load S5 AUC data ----
    s5 = pd.read_excel(
        'data/raw/science.adt3130_table_s5.xlsx',
        sheet_name='eRegulons_Activators_Exp_AUC_RS'
    )
    print(f"  S5 loaded: {len(s5)} rows, {s5['eRegulon'].nunique()} regulons, "
          f"{s5['cell_type_l4'].nunique()} cell types")

    # Pivot to regulon × cell_type matrix
    auc_matrix = s5.pivot_table(
        index='eRegulon', columns='cell_type_l4', values='mean_AUC', aggfunc='first'
    ).fillna(0)
    print(f"  AUC matrix: {auc_matrix.shape}")

    # ---- Assign lineages ----
    all_cts = set(auc_matrix.columns)
    mapped_cts = set(LINEAGE_MAP.keys())
    unmapped = all_cts - mapped_cts
    if unmapped:
        print(f"  WARNING: unmapped cell types: {unmapped}")

    lineages = {}
    for ct, lin in LINEAGE_MAP.items():
        if ct in auc_matrix.columns and lin != 'Progenitor':
            lineages.setdefault(lin, []).append(ct)

    for lin, cts in sorted(lineages.items()):
        print(f"  {lin}: {len(cts)} cell types")

    # ---- Run TOPPLE globally (for comparison) ----
    print("\n--- Global TOPPLE (all 61 cell types) ---")
    all_cell_types = list(auc_matrix.columns)
    global_results = run_topple(
        auc_matrix[all_cell_types],
        list(auc_matrix.index),
        all_cell_types,
        label='Global'
    )
    print(f"  {len(global_results)} regulons scored")
    print(f"  Top 5: {', '.join(global_results.head(5)['regulon'].tolist())}")

    # ---- Run TOPPLE within each lineage ----
    all_lineage_results = [global_results]

    for lineage_name in ['T_cell', 'B_cell', 'Myeloid', 'NK_ILC']:
        cts = lineages[lineage_name]
        print(f"\n--- {lineage_name} TOPPLE ({len(cts)} cell types) ---")

        # Subset matrix
        sub_matrix = auc_matrix[cts]

        # Filter regulons with any activity in this lineage
        active_mask = sub_matrix.max(axis=1) > 0.01
        sub_matrix_active = sub_matrix[active_mask]
        print(f"  Active regulons: {active_mask.sum()}/{len(auc_matrix)}")

        results = run_topple(
            sub_matrix_active,
            list(sub_matrix_active.index),
            cts,
            label=lineage_name
        )
        all_lineage_results.append(results)

        print(f"  Top 10 stabilizers:")
        for _, row in results.head(10).iterrows():
            reg_type = classify_regulon_type(row['regulon'], results)
            print(f"    #{int(row['rank']):3d} {row['regulon']:20s} RI={row['mean_RI']:.5f}  "
                  f"AUC={row['mean_AUC']:.3f}  n_CT={int(row['n_active_cell_types'])}  [{reg_type}]")

        # Count universal vs lineage-specific in top 20
        top20 = results.head(20)
        n_universal = sum(1 for r in top20['regulon']
                         if classify_regulon_type(r, results) == 'universal')
        print(f"  Top 20: {n_universal}/20 universal, {20-n_universal}/20 lineage-specific")

    # ---- Combine and save ----
    combined = pd.concat(all_lineage_results, ignore_index=True)
    combined.to_csv('results/lineage_stratified_topple.csv', index=False)
    print(f"\nSaved results/lineage_stratified_topple.csv ({len(combined)} rows)")

    # ---- Cross-lineage comparison ----
    print("\n" + "=" * 60)
    print("Cross-Lineage Stability Comparison")
    print("=" * 60)

    # For each regulon, get its rank in each lineage
    lineage_names = ['T_cell', 'B_cell', 'Myeloid', 'NK_ILC']
    rank_matrix = {}

    for lin in lineage_names:
        lin_data = combined[combined['lineage'] == lin][['regulon', 'rank']].copy()
        lin_data = lin_data.set_index('regulon')['rank']
        rank_matrix[lin] = lin_data

    rank_df = pd.DataFrame(rank_matrix)
    common = rank_df.dropna()
    print(f"  Regulons present in all 4 lineages: {len(common)}")

    # Pairwise Spearman correlation of rankings
    print(f"\n  Pairwise rank correlation (Spearman rho):")
    corr_results = []
    for i, lin1 in enumerate(lineage_names):
        for j, lin2 in enumerate(lineage_names):
            if i >= j:
                continue
            shared = rank_df[[lin1, lin2]].dropna()
            if len(shared) > 5:
                rho, p = stats.spearmanr(shared[lin1], shared[lin2])
                print(f"    {lin1:10s} vs {lin2:10s}: rho={rho:.3f}, P={p:.2e}, n={len(shared)}")
                corr_results.append({
                    'lineage_1': lin1, 'lineage_2': lin2,
                    'rho': rho, 'p_value': p, 'n_shared': len(shared)
                })

    # Global vs each lineage
    global_ranks = combined[combined['lineage'] == 'Global'][['regulon', 'rank']].set_index('regulon')['rank']
    print(f"\n  Global vs lineage rank correlation:")
    for lin in lineage_names:
        lin_ranks = combined[combined['lineage'] == lin][['regulon', 'rank']].set_index('regulon')['rank']
        shared_regs = set(global_ranks.index) & set(lin_ranks.index)
        if len(shared_regs) > 5:
            g = global_ranks.loc[list(shared_regs)]
            l = lin_ranks.loc[list(shared_regs)]
            rho, p = stats.spearmanr(g, l)
            print(f"    Global vs {lin:10s}: rho={rho:.3f}, P={p:.2e}, n={len(shared_regs)}")
            corr_results.append({
                'lineage_1': 'Global', 'lineage_2': lin,
                'rho': rho, 'p_value': p, 'n_shared': len(shared_regs)
            })

    # ---- Identify conserved vs lineage-specific stabilizers ----
    print(f"\n  Conserved stabilizers (top 20 in >=3 lineages):")
    top20_sets = {}
    for lin in lineage_names:
        lin_top = set(combined[combined['lineage'] == lin].nlargest(20, 'mean_RI')['regulon'])
        top20_sets[lin] = lin_top

    all_top20 = set()
    for s in top20_sets.values():
        all_top20 |= s

    conserved = []
    for reg in all_top20:
        count = sum(1 for lin in lineage_names if reg in top20_sets[lin])
        if count >= 3:
            conserved.append((reg, count))
            reg_type = classify_regulon_type(reg, None)
            print(f"    {reg:20s} in {count}/4 lineages  [{reg_type}]")

    lineage_unique = {}
    for lin in lineage_names:
        unique = top20_sets[lin] - set(r for r, c in conserved)
        for other_lin in lineage_names:
            if other_lin != lin:
                unique -= top20_sets[other_lin]
        if unique:
            lineage_unique[lin] = unique
            print(f"\n  {lin}-specific stabilizers (top 20, unique to this lineage):")
            for reg in sorted(unique):
                rank_in_lin = combined[(combined['lineage'] == lin) &
                                       (combined['regulon'] == reg)]['rank'].values[0]
                print(f"    {reg:20s} rank=#{int(rank_in_lin)}")

    # ---- Summary statistics ----
    summary_rows = []
    for lin in ['Global'] + lineage_names:
        lin_data = combined[combined['lineage'] == lin]
        summary_rows.append({
            'lineage': lin,
            'n_cell_types': len(lineages.get(lin, all_cell_types)),
            'n_regulons': len(lin_data),
            'mean_RI': lin_data['mean_RI'].mean(),
            'std_RI': lin_data['mean_RI'].std(),
            'top1_regulon': lin_data.iloc[0]['regulon'] if len(lin_data) > 0 else '',
            'top1_RI': lin_data.iloc[0]['mean_RI'] if len(lin_data) > 0 else 0,
            'n_universal_in_top20': sum(1 for r in lin_data.head(20)['regulon']
                                        if classify_regulon_type(r, None) == 'universal'),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv('results/lineage_stratified_summary.csv', index=False)
    print(f"\nSaved results/lineage_stratified_summary.csv")

    # Save correlation results
    if corr_results:
        corr_df = pd.DataFrame(corr_results)
        corr_df.to_csv('results/lineage_stratified_correlations.csv', index=False)
        print(f"Saved results/lineage_stratified_correlations.csv")

    # ---- Figure ----
    make_figure(combined, lineage_names, rank_df, corr_results)

    print("\n" + "=" * 60)
    print("Lineage-stratified TOPPLE complete!")
    print("=" * 60)


def make_figure(combined, lineage_names, rank_df, corr_results):
    """Generate 4-panel figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})
    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 160/25.4))

    lineage_colors = {
        'T_cell': '#E74C3C',
        'B_cell': '#3498DB',
        'Myeloid': '#F39C12',
        'NK_ILC': '#2ECC71',
        'Global': '#95A5A6',
    }

    # Panel A: Top 15 per lineage (horizontal grouped bars)
    ax = axes[0, 0]
    n_show = 10
    y_offset = 0
    y_ticks = []
    y_labels = []
    for lin in lineage_names:
        lin_data = combined[combined['lineage'] == lin].nlargest(n_show, 'mean_RI')
        color = lineage_colors[lin]
        for idx, (_, row) in enumerate(lin_data.iterrows()):
            ax.barh(y_offset, row['mean_RI'], height=0.7, color=color, alpha=0.8)
            y_ticks.append(y_offset)
            y_labels.append(row['regulon'])
            y_offset += 1
        y_offset += 1  # gap between lineages

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=4)
    ax.invert_yaxis()
    ax.set_xlabel('Redistribution Index')
    ax.set_title('Top Stabilizers per Lineage', fontsize=9, fontweight='bold')
    ax.legend(handles=[Patch(color=lineage_colors[l], label=l) for l in lineage_names],
              fontsize=5, loc='lower right')

    # Panel B: Rank correlation heatmap
    ax = axes[0, 1]
    all_lins = lineage_names
    n = len(all_lins)
    corr_matrix = np.ones((n, n))
    for cr in corr_results:
        if cr['lineage_1'] in all_lins and cr['lineage_2'] in all_lins:
            i = all_lins.index(cr['lineage_1'])
            j = all_lins.index(cr['lineage_2'])
            corr_matrix[i, j] = cr['rho']
            corr_matrix[j, i] = cr['rho']

    im = ax.imshow(corr_matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=1.0)
    ax.set_xticks(range(n))
    ax.set_xticklabels(all_lins, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(all_lins, fontsize=7)
    # Annotate
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f'{corr_matrix[i,j]:.2f}', ha='center', va='center', fontsize=6)
    plt.colorbar(im, ax=ax, shrink=0.7, label='Spearman rho')
    ax.set_title('Rank Correlation Across Lineages', fontsize=9, fontweight='bold')

    # Panel C: RI distribution per lineage (violin/box)
    ax = axes[1, 0]
    data_for_box = []
    labels_for_box = []
    colors_for_box = []
    for lin in ['Global'] + lineage_names:
        ri_vals = combined[combined['lineage'] == lin]['mean_RI'].values
        data_for_box.append(ri_vals)
        labels_for_box.append(lin)
        colors_for_box.append(lineage_colors[lin])

    bp = ax.boxplot(data_for_box, labels=labels_for_box, patch_artist=True, widths=0.6)
    for patch, color in zip(bp['boxes'], colors_for_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_xticklabels(labels_for_box, rotation=30, ha='right', fontsize=7)
    ax.set_ylabel('Redistribution Index')
    ax.set_title('RI Distribution by Lineage', fontsize=9, fontweight='bold')

    # Panel D: Universal vs lineage-specific in top 20
    ax = axes[1, 1]
    universal_counts = []
    specific_counts = []
    for lin in lineage_names:
        top20 = combined[combined['lineage'] == lin].nlargest(20, 'mean_RI')
        n_univ = sum(1 for r in top20['regulon']
                     if classify_regulon_type(r, None) == 'universal')
        universal_counts.append(n_univ)
        specific_counts.append(20 - n_univ)

    x = np.arange(len(lineage_names))
    ax.bar(x, universal_counts, 0.6, label='Universal (stress-response)',
           color='#E74C3C', alpha=0.7)
    ax.bar(x, specific_counts, 0.6, bottom=universal_counts,
           label='Lineage-specific', color='#3498DB', alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(lineage_names, rotation=30, ha='right', fontsize=7)
    ax.set_ylabel('Count (top 20)')
    ax.set_title('Universal vs Lineage-Specific\nin Top 20 Stabilizers', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6)

    # Panel labels
    for ax, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Fractal Stability Architecture — Lineage-Stratified TOPPLE',
                 fontsize=10, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_lineage_stratified_topple.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_lineage_stratified_topple.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("\nSaved figures/fig_lineage_stratified_topple.pdf + .png")


if __name__ == "__main__":
    main()
