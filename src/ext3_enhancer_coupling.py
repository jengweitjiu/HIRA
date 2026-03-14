#!/usr/bin/env python3
"""
Extension #3 — Enhancer-Driven Coupling Network

Builds a regulatory coupling matrix from peak-gene linkages (S4 GRN)
and peak activity (S3 binary matrix), then performs three-way Mantel tests
against eQTL and caQTL coupling matrices.

Pipeline:
  1. Load S4 -> unique Region-Gene pairs (133,574 links, 84,625 regions, 13,645 genes)
  2. Load S3 -> binary peak x cell-type activity (338,036 peaks x 65 cell types)
  3. For each Region-Gene pair, mark it active in cell types where the Region is active
  4. Build Jaccard coupling matrix between cell types (shared active R-G pairs)
  5. Three-way Mantel: eQTL r_b x caQTL r_b x regulatory Jaccard
  6. Per-regulon: regulatory breadth vs TOPPLE RI

Input:  data/raw/CIMA_Table_S4.csv, data/raw/CIMA_Table_S3.csv,
        results/sicai_rb.csv, results/sicai_caqtl_rb.csv, results/topple.csv
Output: results/ext3_regulatory_coupling.csv, results/ext3_regulon_breadth.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


def load_region_gene_pairs(s4_path):
    """Load unique Region-Gene pairs from S4 GRN table."""
    print(f"Loading {s4_path} ...")
    s4 = pd.read_csv(s4_path, usecols=['Region', 'Gene', 'TF', 'Consensus_name'])
    print(f"  Total rows: {len(s4):,}")
    rg = s4[['Region', 'Gene']].drop_duplicates()
    print(f"  Unique Region-Gene pairs: {len(rg):,}")
    print(f"  Unique Regions: {rg['Region'].nunique():,}")
    print(f"  Unique Genes: {rg['Gene'].nunique():,}")
    return rg, s4


def load_peak_activity(s3_path):
    """Load binary peak x cell-type matrix from S3."""
    print(f"\nLoading {s3_path} ...")
    s3 = pd.read_csv(s3_path, index_col='Peaks')
    print(f"  Shape: {s3.shape} (peaks x cell types)")
    # Convert to bool
    s3 = s3.astype(bool)
    celltypes = list(s3.columns)
    print(f"  Cell types: {len(celltypes)}")
    active_per_ct = s3.sum(axis=0)
    print(f"  Active peaks per cell type: min={active_per_ct.min()}, "
          f"max={active_per_ct.max()}, mean={active_per_ct.mean():.0f}")
    return s3, celltypes


def build_regulatory_coupling(rg_pairs, peak_activity, celltypes):
    """Build Jaccard coupling matrix between cell types based on shared active R-G pairs.

    For each R-G pair, it is 'active' in a cell type if the Region is active there.
    For each pair of cell types, Jaccard = |intersection| / |union| of active R-G sets.
    """
    print("\nBuilding cell-type active R-G pair sets ...")
    regions_in_s3 = set(peak_activity.index)

    # For each cell type, collect set of (Region, Gene) pairs where Region is active
    ct_rg_sets = {}
    # Pre-filter rg_pairs to regions in S3
    rg_in_s3 = rg_pairs[rg_pairs['Region'].isin(regions_in_s3)]
    print(f"  R-G pairs with Region in S3: {len(rg_in_s3):,}")

    # Build lookup: region -> list of genes
    region_to_genes = rg_in_s3.groupby('Region')['Gene'].apply(set).to_dict()

    for ct in celltypes:
        active_peaks = peak_activity.index[peak_activity[ct]].tolist()
        rg_set = set()
        for peak in active_peaks:
            if peak in region_to_genes:
                for gene in region_to_genes[peak]:
                    rg_set.add((peak, gene))
        ct_rg_sets[ct] = rg_set

    # Print stats
    sizes = {ct: len(s) for ct, s in ct_rg_sets.items()}
    print(f"  Active R-G pairs per cell type: min={min(sizes.values()):,}, "
          f"max={max(sizes.values()):,}, mean={np.mean(list(sizes.values())):.0f}")

    # Build Jaccard matrix
    print("\nComputing Jaccard coupling matrix ...")
    n = len(celltypes)
    jaccard_mat = np.zeros((n, n))

    for i in range(n):
        jaccard_mat[i, i] = 1.0
        set_i = ct_rg_sets[celltypes[i]]
        for j in range(i + 1, n):
            set_j = ct_rg_sets[celltypes[j]]
            intersection = len(set_i & set_j)
            union = len(set_i | set_j)
            jac = intersection / union if union > 0 else 0.0
            jaccard_mat[i, j] = jac
            jaccard_mat[j, i] = jac

    jaccard_df = pd.DataFrame(jaccard_mat, index=celltypes, columns=celltypes)
    print(f"  Jaccard matrix shape: {jaccard_df.shape}")
    offdiag = jaccard_mat[np.triu_indices(n, k=1)]
    print(f"  Mean Jaccard (off-diagonal): {offdiag.mean():.4f}")
    print(f"  Median Jaccard: {np.median(offdiag):.4f}")

    return jaccard_df, ct_rg_sets


def load_rb_matrix(sicai_path, s8_path, sheet):
    """Reconstruct r_b matrix from S8 data."""
    s8 = pd.read_excel(s8_path, sheet_name=sheet)
    celltypes = sorted(set(s8['reference_cell_type']) | set(s8['query_celltype']))
    n = len(celltypes)
    ct_idx = {ct: i for i, ct in enumerate(celltypes)}
    rb_mat = np.full((n, n), np.nan)
    np.fill_diagonal(rb_mat, 1.0)
    for _, row in s8.iterrows():
        i = ct_idx[row['reference_cell_type']]
        j = ct_idx[row['query_celltype']]
        rb_mat[j, i] = row['rb']
        if np.isnan(rb_mat[i, j]):
            rb_mat[i, j] = row['rb']
    return pd.DataFrame(rb_mat, index=celltypes, columns=celltypes)


def mantel_test(mat1, mat2, label, n_perm=9999):
    """Mantel test between two matrices on shared cell types."""
    shared = sorted(set(mat1.index) & set(mat2.index))
    n = len(shared)
    if n < 4:
        print(f"  {label}: too few shared cell types ({n})")
        return np.nan, np.nan, 0

    m1 = mat1.loc[shared, shared].values
    m2 = mat2.loc[shared, shared].values
    triu = np.triu_indices(n, k=1)
    v1 = m1[triu]
    v2 = m2[triu]
    valid = ~(np.isnan(v1) | np.isnan(v2))
    v1 = v1[valid]
    v2 = v2[valid]

    r_obs, _ = stats.pearsonr(v1, v2)
    count = 0
    for _ in range(n_perm):
        perm = np.random.permutation(n)
        m2p = m2[np.ix_(perm, perm)]
        v2p = m2p[triu][valid]
        rp, _ = stats.pearsonr(v1, v2p)
        if rp >= r_obs:
            count += 1
    p_val = (count + 1) / (n_perm + 1)

    print(f"  {label}: r={r_obs:.4f}, P={p_val:.4f} ({n} cell types, {len(v1)} pairs)")
    return r_obs, p_val, len(v1)


def regulon_breadth_analysis(s4_full, peak_activity, celltypes, topple_path):
    """Per-regulon: count R-G pairs, compute regulatory breadth, correlate with TOPPLE RI."""
    print("\n" + "=" * 60)
    print("PER-REGULON ANALYSIS")
    print("=" * 60)

    topple = pd.read_csv(topple_path)
    print(f"  TOPPLE regulons: {len(topple)}")

    # Clean consensus name to match TOPPLE regulon names
    s4_full = s4_full.copy()
    # Consensus_name like 'ADNP_+_+' -> take first two parts -> 'ADNP_+'
    s4_full['regulon_clean'] = s4_full['Consensus_name'].str.extract(r'^([A-Za-z0-9]+_[+-])')[0]

    # Per regulon: unique R-G pairs and regions
    reg_stats = s4_full.groupby('regulon_clean').agg(
        n_rg_pairs=('Gene', 'count'),
        n_unique_regions=('Region', 'nunique'),
        n_unique_genes=('Gene', 'nunique'),
    ).reset_index()

    # Regulatory breadth: for each regulon, across how many cell types
    # are its regions active? (mean fraction of cell types per region)
    regions_in_s3 = set(peak_activity.index)
    breadth_records = []
    for regulon in reg_stats['regulon_clean'].unique():
        reg_rows = s4_full[s4_full['regulon_clean'] == regulon]
        regions = reg_rows['Region'].unique()
        regions_active = [r for r in regions if r in regions_in_s3]
        if not regions_active:
            breadth_records.append({
                'regulon_clean': regulon,
                'mean_ct_breadth': 0.0,
                'max_ct_breadth': 0,
            })
            continue

        ct_counts = []
        for region in regions_active:
            n_active = peak_activity.loc[region].sum()
            ct_counts.append(n_active)

        breadth_records.append({
            'regulon_clean': regulon,
            'mean_ct_breadth': np.mean(ct_counts),
            'max_ct_breadth': np.max(ct_counts),
        })

    breadth_df = pd.DataFrame(breadth_records)
    reg_stats = reg_stats.merge(breadth_df, on='regulon_clean', how='left')

    # Merge with TOPPLE
    merged = reg_stats.merge(topple[['regulon', 'mean_RI', 'rank', 'stability_class']],
                             left_on='regulon_clean', right_on='regulon', how='inner')
    print(f"  Regulons matched to TOPPLE: {len(merged)}")

    # Correlations
    rho_rg, p_rg = stats.spearmanr(merged['n_rg_pairs'], merged['mean_RI'])
    rho_br, p_br = stats.spearmanr(merged['mean_ct_breadth'], merged['mean_RI'])
    rho_genes, p_genes = stats.spearmanr(merged['n_unique_genes'], merged['mean_RI'])

    print(f"\n  Spearman correlations with TOPPLE RI:")
    print(f"    n_rg_pairs vs RI:        rho={rho_rg:.4f}, P={p_rg:.2e}")
    print(f"    mean_ct_breadth vs RI:   rho={rho_br:.4f}, P={p_br:.2e}")
    print(f"    n_unique_genes vs RI:    rho={rho_genes:.4f}, P={p_genes:.2e}")

    # Stabilizer vs destabilizer comparison
    stab = merged[merged['stability_class'] == 'stabilizer']
    dest = merged[merged['stability_class'] == 'destabilizer']
    for metric in ['n_rg_pairs', 'n_unique_regions', 'n_unique_genes', 'mean_ct_breadth']:
        u, p = stats.mannwhitneyu(stab[metric], dest[metric],
                                  alternative='two-sided', method='asymptotic')
        print(f"\n    {metric}:")
        print(f"      Stabilizer mean={stab[metric].mean():.1f}, Destabilizer mean={dest[metric].mean():.1f}")
        print(f"      MW P={p:.2e}")

    return merged


def main():
    parser = argparse.ArgumentParser(description='Extension #3: enhancer-driven coupling')
    parser.add_argument('--s4', default='data/raw/CIMA_Table_S4.csv')
    parser.add_argument('--s3', default='data/raw/CIMA_Table_S3.csv')
    parser.add_argument('--s8', default='data/raw/science.adt3130_table_s8.xlsx')
    parser.add_argument('--topple', default='results/topple.csv')
    parser.add_argument('--output-coupling', default='results/ext3_regulatory_coupling.csv')
    parser.add_argument('--output-breadth', default='results/ext3_regulon_breadth.csv')
    args = parser.parse_args()

    # Step 1: Load S4 Region-Gene pairs
    rg_pairs, s4_full = load_region_gene_pairs(args.s4)

    # Step 2: Load S3 peak activity
    peak_activity, celltypes_s3 = load_peak_activity(args.s3)

    # Step 3: Build Jaccard coupling matrix
    jaccard_df, ct_rg_sets = build_regulatory_coupling(rg_pairs, peak_activity, celltypes_s3)

    # Per-cell-type stats
    ct_stats = []
    for ct in celltypes_s3:
        ct_stats.append({
            'cell_type': ct,
            'n_active_rg_pairs': len(ct_rg_sets[ct]),
            'mean_jaccard': jaccard_df.loc[ct].drop(ct).mean(),
        })
    ct_df = pd.DataFrame(ct_stats).sort_values('mean_jaccard', ascending=False).reset_index(drop=True)

    Path(args.output_coupling).parent.mkdir(parents=True, exist_ok=True)
    ct_df.to_csv(args.output_coupling, index=False)
    print(f"\nSaved per-cell-type regulatory coupling to {args.output_coupling}")

    # Step 4-5: Three-way Mantel tests
    print("\n" + "=" * 60)
    print("THREE-WAY MANTEL TESTS")
    print("=" * 60)

    rb_eqtl = load_rb_matrix(None, args.s8, 'cis_eQTL')
    rb_caqtl = load_rb_matrix(None, args.s8, 'cis_caQTL')
    print(f"  eQTL rb: {rb_eqtl.shape}, caQTL rb: {rb_caqtl.shape}, Jaccard: {jaccard_df.shape}")

    np.random.seed(42)
    r_ec, p_ec, n_ec = mantel_test(rb_eqtl, rb_caqtl, "eQTL rb vs caQTL rb")
    r_ej, p_ej, n_ej = mantel_test(rb_eqtl, jaccard_df, "eQTL rb vs Regulatory Jaccard")
    r_cj, p_cj, n_cj = mantel_test(rb_caqtl, jaccard_df, "caQTL rb vs Regulatory Jaccard")

    print(f"\n  Summary:")
    print(f"    eQTL r_b   vs caQTL r_b:   r={r_ec:.4f}, P={p_ec:.4f}")
    print(f"    eQTL r_b   vs Reg Jaccard: r={r_ej:.4f}, P={p_ej:.4f}")
    print(f"    caQTL r_b  vs Reg Jaccard: r={r_cj:.4f}, P={p_cj:.4f}")

    # Step 6: Per-regulon analysis
    regulon_df = regulon_breadth_analysis(s4_full, peak_activity, celltypes_s3, args.topple)
    regulon_df.to_csv(args.output_breadth, index=False)
    print(f"\nSaved regulon breadth to {args.output_breadth}")

    # Final summary
    print("\n" + "=" * 60)
    print("EXTENSION #3 KEY STATISTICS")
    print("=" * 60)
    print(f"  GRN: {rg_pairs['Region'].nunique():,} regions -> {rg_pairs['Gene'].nunique():,} genes "
          f"({len(rg_pairs):,} links)")
    print(f"  Peak activity: {peak_activity.shape[0]:,} peaks x {peak_activity.shape[1]} cell types")
    print(f"  Regulatory Jaccard matrix: {jaccard_df.shape}")
    print(f"  Mean Jaccard: {jaccard_df.values[np.triu_indices(len(celltypes_s3), k=1)].mean():.4f}")
    print(f"\n  Mantel tests:")
    print(f"    eQTL vs caQTL:   r={r_ec:.4f}")
    print(f"    eQTL vs Jaccard: r={r_ej:.4f}")
    print(f"    caQTL vs Jaccard: r={r_cj:.4f}")
    print("=" * 60)


if __name__ == '__main__':
    main()
