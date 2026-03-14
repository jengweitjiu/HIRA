#!/usr/bin/env python3
"""
DECODE-AD gap analyses:
G. LD expansion verification (already done -- confirm 500kb window)
H. IPA kill experiment on AD regulons
I. Permutation test for coupling disruption
J. Type 2 axis analysis
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def ipa_kill_experiment(s5_path, topple_path, eqtl_path, s4_path):
    """IPA: compute RI for AD TF regulons vs non-AD regulons."""
    print("=" * 60)
    print("AIM 2 GAP: IPA kill experiment on AD regulons")
    print("=" * 60)

    # Load S5 AUC matrix
    s5 = pd.read_excel(s5_path, sheet_name='eRegulons_Activators_Exp_AUC_RS')
    print(f"  S5 columns (first 5): {s5.columns[:5].tolist()}")

    # Identify regulon and cell_type columns
    if 'eRegulon' in s5.columns:
        s5 = s5.rename(columns={'eRegulon': 'regulon', 'cell_type_l4': 'cell_type'})
    elif 'regulon' not in s5.columns:
        # Try to find the right columns
        print(f"  All columns: {s5.columns.tolist()}")
        return None

    # Pivot to matrix
    if 'mean_AUC' in s5.columns:
        auc_col = 'mean_AUC'
    elif 'AUC' in s5.columns:
        auc_col = 'AUC'
    else:
        # Find AUC-like column
        auc_cols = [c for c in s5.columns if 'auc' in c.lower()]
        auc_col = auc_cols[0] if auc_cols else s5.columns[2]
    print(f"  Using AUC column: {auc_col}")

    mat = s5.pivot_table(index='regulon', columns='cell_type', values=auc_col, fill_value=0)
    print(f"  AUC matrix: {mat.shape}")

    # L1 normalize each cell type
    mat_norm = mat.div(mat.sum(axis=0), axis=1)

    # Compute RI per regulon (leave-one-out JSD)
    def compute_ri(regulon_name):
        original = mat_norm.copy()
        perturbed = original.copy()
        perturbed.loc[regulon_name] = 0
        perturbed = perturbed.div(perturbed.sum(axis=0), axis=1)

        ri_values = []
        for ct in mat_norm.columns:
            p = original[ct].values
            q = perturbed[ct].values
            if np.any(np.isnan(p)) or np.any(np.isnan(q)):
                continue
            jsd = jensenshannon(p, q) ** 2
            ri_values.append(jsd)
        return np.mean(ri_values) if ri_values else np.nan

    # Load AD eGenes and S4 to identify AD-related TFs
    ad_eqtl = pd.read_csv(eqtl_path)
    ad_egenes = set(ad_eqtl['phenotype_id'].unique())

    s4 = pd.read_csv(s4_path)
    s4_tfs = set(s4['TF'].unique())

    # AD TFs: genes that are both AD eGenes and TFs in S4
    ad_tfs = ad_egenes & s4_tfs
    print(f"  AD eGene TFs: {sorted(ad_tfs)}")

    # Map TF names to regulon names (add _+ suffix)
    regulon_names = set(mat.index)
    ad_regulons = set()
    non_ad_regulons = set()
    for reg in regulon_names:
        tf_name = reg.replace('_+', '').replace('_-', '')
        if tf_name in ad_tfs:
            ad_regulons.add(reg)
        else:
            non_ad_regulons.add(reg)

    print(f"  AD regulons in matrix: {sorted(ad_regulons)}")
    print(f"  Non-AD regulons: {len(non_ad_regulons)}")

    # Compute RI for all regulons (reuse from topple.csv if available)
    topple = pd.read_csv(topple_path)
    topple_ri = topple.set_index('regulon')['mean_RI'].to_dict()

    results = []
    for reg in regulon_names:
        ri = topple_ri.get(reg, np.nan)
        tf_name = reg.replace('_+', '').replace('_-', '')
        results.append({
            'regulon': reg,
            'tf_name': tf_name,
            'mean_RI': ri,
            'is_ad_regulon': reg in ad_regulons,
            'stability_class': topple.set_index('regulon').get('stability_class', pd.Series()).get(reg, 'unknown'),
        })

    result_df = pd.DataFrame(results)

    # Mann-Whitney: AD regulon RI vs non-AD regulon RI
    ad_ri = result_df[result_df['is_ad_regulon']]['mean_RI'].dropna()
    nonad_ri = result_df[~result_df['is_ad_regulon']]['mean_RI'].dropna()

    print(f"\n  AD regulon RI (n={len(ad_ri)}):")
    for _, row in result_df[result_df['is_ad_regulon']].iterrows():
        print(f"    {row['regulon']:15s}  RI={row['mean_RI']:.6f}  [{row['stability_class']}]")

    if len(ad_ri) >= 2 and len(nonad_ri) >= 2:
        u_stat, u_p = stats.mannwhitneyu(ad_ri, nonad_ri, alternative='two-sided')
        print(f"\n  Mann-Whitney U test:")
        print(f"    AD mean RI:     {ad_ri.mean():.6f} (median {ad_ri.median():.6f})")
        print(f"    Non-AD mean RI: {nonad_ri.mean():.6f} (median {nonad_ri.median():.6f})")
        print(f"    U={u_stat:.0f}, P={u_p:.2e}")
        direction = "AD > non-AD" if ad_ri.mean() > nonad_ri.mean() else "non-AD > AD"
        print(f"    Direction: {direction}")

    result_df.to_csv('results/decode_ad_ipa_kill.csv', index=False)
    print(f"\nSaved to results/decode_ad_ipa_kill.csv")
    return result_df


def permutation_test(s6_path, s8_path, eqtl_path, n_perm=10000):
    """Permutation test for AD coupling disruption significance.

    Strategy: For each pair of cell types, compute correlation of eQTL slopes
    across AD genes (AD-specific coupling). Then permute gene labels to get null.
    Compare mean AD coupling to null distribution.
    """
    print("\n" + "=" * 60)
    print("AIM 3 GAP: Permutation test for coupling disruption")
    print("=" * 60)

    # Load S8 r_b (background coupling)
    s8 = pd.read_excel(s8_path, sheet_name='cis_eQTL')
    cts_s8 = sorted(set(s8['reference_cell_type']) | set(s8['query_celltype']))
    rb_bg = {}
    for _, row in s8.iterrows():
        key = tuple(sorted([row['reference_cell_type'], row['query_celltype']]))
        rb_bg[key] = row['rb']

    # Load S6 eQTL data
    s6 = pd.read_csv(s6_path)
    s6_eqtl = s6[s6['analysis'] == 'cis-eQTL']

    # Load AD eGenes
    ad_eqtl = pd.read_csv(eqtl_path)
    ad_genes = list(ad_eqtl['phenotype_id'].unique())
    all_genes = list(s6_eqtl['phenotype_id'].unique())
    print(f"  AD eGenes: {len(ad_genes)}")
    print(f"  Total eQTL genes: {len(all_genes)}")

    # Build gene x cell-type slope matrix for efficient computation
    print("  Building gene x cell-type slope matrix...")
    pivot = s6_eqtl.pivot_table(index='phenotype_id', columns='celltype',
                                 values='slope', aggfunc='first')
    print(f"  Slope matrix: {pivot.shape}")

    def compute_mean_pairwise_corr(gene_list):
        """Compute mean pairwise Pearson correlation of slopes across gene set."""
        sub = pivot.loc[pivot.index.isin(gene_list)]
        if sub.shape[0] < 5:
            return np.nan
        # Compute cell-type pairwise correlations using gene slopes
        ct_cols = [c for c in sub.columns if c in cts_s8]
        sub_valid = sub[ct_cols].dropna(axis=1, how='all')
        if sub_valid.shape[1] < 2:
            return np.nan
        # Mean pairwise correlation between cell types
        corr_mat = sub_valid.corr()
        mask = np.triu(np.ones(corr_mat.shape, dtype=bool), k=1)
        vals = corr_mat.values[mask]
        return np.nanmean(vals)

    # Observed
    observed = compute_mean_pairwise_corr(ad_genes)
    global_val = compute_mean_pairwise_corr(all_genes)
    print(f"  Observed AD pairwise correlation: {observed:.4f}")
    print(f"  Global pairwise correlation: {global_val:.4f}")
    print(f"  Disruption (global - AD): {global_val - observed:.4f}")

    # Permutation
    print(f"\n  Running {n_perm} permutations...")
    np.random.seed(42)
    null_vals = []
    n_ad = len(ad_genes)
    genes_array = np.array(all_genes)

    for i in range(n_perm):
        perm_genes = np.random.choice(genes_array, size=n_ad, replace=False)
        perm_val = compute_mean_pairwise_corr(perm_genes)
        null_vals.append(perm_val)
        if (i + 1) % 2000 == 0:
            print(f"    {i+1}/{n_perm} done (last={perm_val:.4f})")

    null_array = np.array(null_vals)
    null_clean = null_array[~np.isnan(null_array)]

    if len(null_clean) > 0 and not np.isnan(observed):
        p_lower = np.mean(null_clean <= observed)
        p_upper = np.mean(null_clean >= observed)
        p_two = 2 * min(p_lower, p_upper)
        z_score = (observed - null_clean.mean()) / null_clean.std() if null_clean.std() > 0 else 0

        print(f"\n  Permutation results ({len(null_clean)} valid / {n_perm}):")
        print(f"    Null mean: {null_clean.mean():.4f} +/- {null_clean.std():.4f}")
        print(f"    Observed:  {observed:.4f}")
        print(f"    P (AD weaker coupling): {p_lower:.4f}")
        print(f"    P (two-tailed): {p_two:.4f}")
        print(f"    Z-score: {z_score:.2f}")
    else:
        p_lower = p_two = z_score = np.nan
        print("  Could not compute permutation P-value")

    perm_df = pd.DataFrame({
        'metric': ['observed_coupling', 'global_coupling', 'disruption',
                    'null_mean', 'null_std', 'p_one_tailed', 'p_two_tailed',
                    'z_score', 'n_permutations', 'n_ad_genes'],
        'value': [observed, global_val, global_val - observed if not np.isnan(global_val) else np.nan,
                  null_clean.mean() if len(null_clean) > 0 else np.nan,
                  null_clean.std() if len(null_clean) > 0 else np.nan,
                  p_lower, p_two, z_score, n_perm, n_ad],
    })
    perm_df.to_csv('results/decode_ad_permutation.csv', index=False)
    print(f"\nSaved to results/decode_ad_permutation.csv")
    return perm_df


def type2_axis(s8_path, eqtl_path, coupling_path):
    """Analyze Th2-related cell type coupling in AD context."""
    print("\n" + "=" * 60)
    print("AIM 3 GAP: Type 2 axis analysis")
    print("=" * 60)

    # Th2-related cell types in CIMA
    th2_types = [
        'CD4_Th_CCR4',           # Th2 proxy
        'ILC2_IL2RA',            # ILC2
        'pDC_IRF4',              # pDC
        'DC1_CLEC9A',            # cDC1
        'DC2_CD1C',              # cDC2
        'CD4_Treg_FOXP3',       # Treg (tolerance)
        'CD4_Th22-like_CCR10',  # Th22 (skin-relevant)
        'CD4_Th17-like_RORC',   # Th17
    ]

    # Load eQTL r_b
    s8 = pd.read_excel(s8_path, sheet_name='cis_eQTL')
    cts = sorted(set(s8['reference_cell_type']) | set(s8['query_celltype']))
    rb_matrix = pd.DataFrame(0.0, index=cts, columns=cts)
    for _, row in s8.iterrows():
        rb_matrix.loc[row['reference_cell_type'], row['query_celltype']] = row['rb']
        rb_matrix.loc[row['query_celltype'], row['reference_cell_type']] = row['rb']
    np.fill_diagonal(rb_matrix.values, 1.0)

    # Filter to available Th2 types
    available = [ct for ct in th2_types if ct in cts]
    print(f"  Th2-related cell types available: {len(available)}")
    for ct in available:
        print(f"    {ct}")

    # Background r_b for Th2 pairs
    th2_pairs = []
    for i in range(len(available)):
        for j in range(i+1, len(available)):
            rb = rb_matrix.loc[available[i], available[j]]
            th2_pairs.append({
                'cell_type_1': available[i],
                'cell_type_2': available[j],
                'global_rb': rb,
            })
    th2_bg = pd.DataFrame(th2_pairs)
    print(f"\n  Th2 pairwise background r_b (n={len(th2_bg)}):")
    print(f"    Mean: {th2_bg['global_rb'].mean():.3f}")

    # Load AD coupling
    coupling = pd.read_csv(coupling_path)

    # Get AD coupling for Th2 pairs
    th2_ad = []
    for _, pair in th2_bg.iterrows():
        ct1, ct2 = pair['cell_type_1'], pair['cell_type_2']
        match = coupling[
            ((coupling['cell_type_1'] == ct1) & (coupling['cell_type_2'] == ct2)) |
            ((coupling['cell_type_1'] == ct2) & (coupling['cell_type_2'] == ct1))
        ]
        if len(match) > 0:
            ad_r = match.iloc[0]['ad_coupling_r']
            th2_ad.append({
                'cell_type_1': ct1,
                'cell_type_2': ct2,
                'global_rb': pair['global_rb'],
                'ad_coupling': ad_r,
                'disruption': pair['global_rb'] - ad_r,
            })

    th2_ad_df = pd.DataFrame(th2_ad)
    if len(th2_ad_df) > 0:
        print(f"\n  Th2 pairs with AD coupling data: {len(th2_ad_df)}")
        print(f"    Global mean r_b: {th2_ad_df['global_rb'].mean():.3f}")
        print(f"    AD mean coupling: {th2_ad_df['ad_coupling'].mean():.3f}")
        print(f"    Mean disruption: {th2_ad_df['disruption'].mean():.3f}")

        # Compare with non-Th2 pairs
        all_pairs = coupling[coupling['global_rb'] > 0].copy()
        all_pairs['disruption'] = all_pairs['global_rb'] - all_pairs['ad_coupling_r']
        non_th2 = all_pairs[
            ~(all_pairs['cell_type_1'].isin(available) & all_pairs['cell_type_2'].isin(available))
        ]

        print(f"\n  Non-Th2 pairs: mean disruption = {non_th2['disruption'].mean():.3f} (n={len(non_th2)})")

        if len(th2_ad_df) >= 3:
            u, p = stats.mannwhitneyu(
                th2_ad_df['disruption'], non_th2['disruption'], alternative='two-sided'
            )
            print(f"  Th2 vs non-Th2 disruption: U={u:.0f}, P={p:.2e}")

        # Detail each pair
        print(f"\n  Th2 pair details:")
        for _, row in th2_ad_df.sort_values('disruption', ascending=False).iterrows():
            print(f"    {row['cell_type_1']:25s} <-> {row['cell_type_2']:25s}  "
                  f"global={row['global_rb']:.3f}  AD={row['ad_coupling']:.3f}  "
                  f"disruption={row['disruption']:.3f}")

    th2_ad_df.to_csv('results/decode_ad_type2_axis.csv', index=False)
    print(f"\nSaved to results/decode_ad_type2_axis.csv")
    return th2_ad_df


if __name__ == "__main__":
    ipa_kill_experiment(
        'data/raw/science.adt3130_table_s5.xlsx',
        'results/topple.csv',
        'results/decode_ad_eqtl.csv',
        'data/raw/CIMA_Table_S4.csv',
    )

    permutation_test(
        'data/raw/CIMA_Table_S6.csv',
        'data/raw/science.adt3130_table_s8.xlsx',
        'results/decode_ad_eqtl.csv',
        n_perm=10000,
    )

    type2_axis(
        'data/raw/science.adt3130_table_s8.xlsx',
        'results/decode_ad_eqtl.csv',
        'results/decode_ad_coupling.csv',
    )
