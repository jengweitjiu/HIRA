#!/usr/bin/env python3
"""
DECODE-AD — Decoding Atopic Dermatitis through CIMA Cell-Type-Resolved xQTLs

Maps AD GWAS lead loci (Budu-Aggrey et al., Nat Commun 2023; GCST90027161)
onto CIMA cis-eQTL and cis-caQTL catalogues to resolve disease genetics at
cell-type resolution.

Pipeline:
  1. Load AD GWAS summary stats, extract genome-wide significant SNPs (P < 5e-8)
  2. Clump into independent lead loci (1 Mb window, keep lowest P)
  3. Load CIMA S6 eQTLs/caQTLs, match by chr:pos coordinate
  4. Build AD-locus x cell-type matrices for eQTL and caQTL
  5. Apply DGSA: compare non-additivity of AD genes vs non-AD genes (Mann-Whitney)

Input:  data/raw/AD_GWAS_Budu-Aggrey_2023.tsv.gz
        data/raw/CIMA_Table_S6.csv (from science.adt3130_table_s6.zip)
        results/dgsa_eqtl.csv (for non-AD gene comparison)
Output: results/decode_ad_eqtl.csv, results/decode_ad_caqtl.csv,
        results/decode_ad_locus_celltype.csv, results/decode_ad_summary.csv
"""

import argparse
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
from collections import defaultdict


def load_gwas(gwas_path, p_threshold=5e-8):
    """Load AD GWAS, extract genome-wide significant SNPs."""
    print(f"Loading {gwas_path} ...")
    chunks = []
    with gzip.open(gwas_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        print(f"  Columns: {header}")
        for line in f:
            fields = line.strip().split('\t')
            try:
                pval = float(fields[1])
            except (ValueError, IndexError):
                continue
            if pval < p_threshold:
                chunks.append(fields)

    gwas = pd.DataFrame(chunks, columns=header)
    for col in ['p_value', 'beta', 'standard_error', 'effect_allele_frequency']:
        if col in gwas.columns:
            gwas[col] = pd.to_numeric(gwas[col], errors='coerce')
    gwas['chromosome'] = gwas['chromosome'].astype(str)
    gwas['base_pair_location'] = pd.to_numeric(gwas['base_pair_location'], errors='coerce').astype(int)

    print(f"  Total genome-wide significant (P < {p_threshold}): {len(gwas):,}")
    print(f"  Chromosomes: {sorted(gwas['chromosome'].unique(), key=lambda x: int(x) if x.isdigit() else 99)}")
    return gwas


def clump_loci(gwas, window=1_000_000):
    """Simple distance-based clumping: within each chromosome, group SNPs
    within `window` bp and keep the one with lowest P as the lead."""
    print(f"\nClumping loci (window = {window/1e6:.0f} Mb) ...")
    leads = []
    for chrom in gwas['chromosome'].unique():
        sub = gwas[gwas['chromosome'] == chrom].sort_values('base_pair_location').reset_index(drop=True)
        clusters = []
        current_cluster = [0]
        for i in range(1, len(sub)):
            if sub.loc[i, 'base_pair_location'] - sub.loc[current_cluster[-1], 'base_pair_location'] <= window:
                current_cluster.append(i)
            else:
                clusters.append(current_cluster)
                current_cluster = [i]
        clusters.append(current_cluster)

        for cluster in clusters:
            cluster_df = sub.iloc[cluster]
            lead_idx = cluster_df['p_value'].idxmin()
            lead = cluster_df.loc[lead_idx].copy()
            lead['n_snps_in_locus'] = len(cluster)
            lead['locus_start'] = cluster_df['base_pair_location'].min()
            lead['locus_end'] = cluster_df['base_pair_location'].max()
            leads.append(lead)

    leads_df = pd.DataFrame(leads).sort_values('p_value').reset_index(drop=True)
    leads_df['locus_id'] = [f"AD_locus_{i+1}" for i in range(len(leads_df))]
    print(f"  Independent loci: {len(leads_df)}")
    print(f"  Top 10 lead SNPs:")
    for _, row in leads_df.head(10).iterrows():
        print(f"    {row['variant_id']:20s}  chr{row['chromosome']}:{row['base_pair_location']:>12,}  "
              f"P={row['p_value']:.2e}  beta={row['beta']:.4f}")
    return leads_df


def load_cima_xqtl(s6_path, analysis_type):
    """Load CIMA S6 xQTLs of specified type, parse chr:pos from variant_id."""
    print(f"\nLoading CIMA {analysis_type} from {s6_path} ...")
    df = pd.read_csv(s6_path, usecols=['phenotype_id', 'variant_id', 'celltype',
                                        'slope', 'slope_se', 'pval_nominal', 'analysis'])
    df = df[df['analysis'] == analysis_type].copy()
    print(f"  {analysis_type} rows: {len(df):,}")

    # Parse variant_id: chr1_814733 -> (1, 814733)
    parts = df['variant_id'].str.extract(r'chr(\w+)_(\d+)')
    df['chrom'] = parts[0].astype(str)
    df['pos'] = pd.to_numeric(parts[1], errors='coerce')
    df = df.dropna(subset=['pos'])
    df['pos'] = df['pos'].astype(int)

    print(f"  Unique variants: {df['variant_id'].nunique():,}")
    print(f"  Unique genes/peaks: {df['phenotype_id'].nunique():,}")
    print(f"  Cell types: {df['celltype'].nunique()}")
    return df


def match_gwas_to_xqtl(gwas_leads, xqtl_df, window=500_000, label='eQTL'):
    """Match GWAS lead loci to xQTLs by proximity (within window bp)."""
    print(f"\nMatching AD loci to CIMA {label} (window = +/-{window/1e3:.0f} kb) ...")

    # Build lookup: for each chromosome, list of (pos, locus_id, rsid)
    gwas_by_chr = defaultdict(list)
    for _, row in gwas_leads.iterrows():
        gwas_by_chr[str(row['chromosome'])].append(
            (int(row['base_pair_location']), row['locus_id'], row['variant_id']))

    # For each xQTL, check if it falls within window of any GWAS lead
    matches = []
    for chrom, loci in gwas_by_chr.items():
        xqtl_chr = xqtl_df[xqtl_df['chrom'] == chrom]
        if len(xqtl_chr) == 0:
            continue
        for gwas_pos, locus_id, rsid in loci:
            nearby = xqtl_chr[
                (xqtl_chr['pos'] >= gwas_pos - window) &
                (xqtl_chr['pos'] <= gwas_pos + window)
            ].copy()
            if len(nearby) > 0:
                nearby['locus_id'] = locus_id
                nearby['gwas_rsid'] = rsid
                nearby['gwas_pos'] = gwas_pos
                nearby['distance'] = (nearby['pos'] - gwas_pos).abs()
                matches.append(nearby)

    if not matches:
        print(f"  No matches found!")
        return pd.DataFrame()

    matched = pd.concat(matches, ignore_index=True)
    print(f"  Matched xQTL entries: {len(matched):,}")
    print(f"  AD loci with xQTL overlap: {matched['locus_id'].nunique()}")
    print(f"  Unique genes/peaks: {matched['phenotype_id'].nunique()}")
    print(f"  Cell types represented: {matched['celltype'].nunique()}")
    return matched


def build_locus_celltype_matrix(matched_df, label='eQTL'):
    """Build AD locus x cell type matrix (count of significant xQTLs per cell)."""
    print(f"\nBuilding AD locus x cell type matrix ({label}) ...")

    # Count unique genes/peaks per locus x cell type
    pivot = matched_df.groupby(['locus_id', 'celltype'])['phenotype_id'].nunique().reset_index()
    pivot.columns = ['locus_id', 'celltype', 'n_genes']
    matrix = pivot.pivot_table(index='locus_id', columns='celltype',
                               values='n_genes', fill_value=0)
    print(f"  Matrix shape: {matrix.shape} (loci x cell types)")
    print(f"  Mean genes per locus-celltype: {matrix.values[matrix.values > 0].mean():.2f}")

    # Top cell types by total genes
    ct_totals = matrix.sum(axis=0).sort_values(ascending=False)
    print(f"  Top 5 cell types by total xQTL genes:")
    for ct, count in ct_totals.head(5).items():
        print(f"    {ct:35s}  {int(count)} genes")

    return matrix


def compute_gini(arr):
    """Gini coefficient of an array."""
    arr = np.sort(np.abs(arr))
    n = len(arr)
    if n == 0 or arr.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * arr) - (n + 1) * np.sum(arr)) / (n * np.sum(arr))


def dgsa_ad_comparison(matched_eqtl, xqtl_full, dgsa_path, all_celltypes):
    """Compare non-additivity of AD-overlapping genes vs non-AD genes."""
    print("\n" + "=" * 60)
    print("DGSA: AD vs non-AD gene non-additivity comparison")
    print("=" * 60)

    # Get set of AD-overlapping genes
    ad_genes = set(matched_eqtl['phenotype_id'].unique())
    print(f"  AD-overlapping genes (from locus matching): {len(ad_genes)}")

    # Load pre-computed DGSA results
    dgsa = pd.read_csv(dgsa_path)
    print(f"  DGSA total genes: {len(dgsa)}")
    print(f"  Columns: {list(dgsa.columns)}")

    # Split into AD and non-AD
    dgsa['is_AD'] = dgsa['gene'].isin(ad_genes)
    ad_dgsa = dgsa[dgsa['is_AD']]
    nonad_dgsa = dgsa[~dgsa['is_AD']]
    print(f"  AD genes in DGSA: {len(ad_dgsa)}")
    print(f"  Non-AD genes in DGSA: {len(nonad_dgsa)}")

    if len(ad_dgsa) < 5:
        print("  Too few AD genes for comparison.")
        return dgsa

    # Compare non-additivity
    u_stat, p_val = stats.mannwhitneyu(
        ad_dgsa['non_additivity'], nonad_dgsa['non_additivity'],
        alternative='two-sided')
    _, p_less = stats.mannwhitneyu(
        ad_dgsa['non_additivity'], nonad_dgsa['non_additivity'],
        alternative='less')
    _, p_greater = stats.mannwhitneyu(
        ad_dgsa['non_additivity'], nonad_dgsa['non_additivity'],
        alternative='greater')

    print(f"\n  Non-additivity comparison:")
    print(f"    AD genes:      mean={ad_dgsa['non_additivity'].mean():.4f}, "
          f"median={ad_dgsa['non_additivity'].median():.4f} (n={len(ad_dgsa)})")
    print(f"    Non-AD genes:  mean={nonad_dgsa['non_additivity'].mean():.4f}, "
          f"median={nonad_dgsa['non_additivity'].median():.4f} (n={len(nonad_dgsa)})")
    print(f"    MW U = {u_stat:.0f}")
    print(f"    P (two-sided):    {p_val:.2e}")
    print(f"    P (AD < non-AD):  {p_less:.2e}")
    print(f"    P (AD > non-AD):  {p_greater:.2e}")

    # Also compare other metrics if available
    for metric in ['gini', 'sparsity', 'magnitude']:
        if metric in dgsa.columns:
            u, p = stats.mannwhitneyu(
                ad_dgsa[metric], nonad_dgsa[metric], alternative='two-sided')
            print(f"\n    {metric}: AD mean={ad_dgsa[metric].mean():.4f} vs "
                  f"non-AD mean={nonad_dgsa[metric].mean():.4f}, P={p:.2e}")

    # Also compute non-additivity directly for AD-matched eQTLs
    # using the full slope profiles
    print("\n  Computing non-additivity for AD eQTLs de novo ...")
    eqtl_sub = xqtl_full[xqtl_full['analysis'] == 'cis-eQTL']
    ad_eqtl = eqtl_sub[eqtl_sub['phenotype_id'].isin(ad_genes)]
    ad_na_records = []
    for gene, grp in ad_eqtl.groupby('phenotype_id'):
        if grp['celltype'].nunique() < 3:
            continue
        slopes_full = np.zeros(len(all_celltypes))
        ct_idx = {ct: i for i, ct in enumerate(all_celltypes)}
        for _, row in grp.iterrows():
            if row['celltype'] in ct_idx:
                slopes_full[ct_idx[row['celltype']]] = row['slope']
        na_val = compute_gini(slopes_full ** 2)
        ad_na_records.append({
            'gene': gene,
            'non_additivity_denovo': na_val,
            'n_celltypes': grp['celltype'].nunique(),
        })

    if ad_na_records:
        ad_na_df = pd.DataFrame(ad_na_records)
        print(f"    De novo AD genes (>=3 CTs): {len(ad_na_df)}")
        print(f"    Mean non-additivity: {ad_na_df['non_additivity_denovo'].mean():.4f}")

    return dgsa


def main():
    parser = argparse.ArgumentParser(description='DECODE-AD: AD GWAS x CIMA xQTL integration')
    parser.add_argument('--gwas', default='data/raw/AD_GWAS_Budu-Aggrey_2023.tsv.gz')
    parser.add_argument('--s6', default='data/raw/CIMA_Table_S6.csv')
    parser.add_argument('--dgsa', default='results/dgsa_eqtl.csv')
    parser.add_argument('--p-threshold', type=float, default=5e-8)
    parser.add_argument('--clump-window', type=int, default=1_000_000)
    parser.add_argument('--match-window', type=int, default=500_000)
    parser.add_argument('--output-dir', default='results')
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load GWAS and extract significant SNPs
    gwas = load_gwas(args.gwas, args.p_threshold)

    # Step 2: Clump into independent loci
    leads = clump_loci(gwas, window=args.clump_window)

    # Step 3: Load CIMA xQTLs
    eqtl = load_cima_xqtl(args.s6, 'cis-eQTL')
    caqtl = load_cima_xqtl(args.s6, 'cis-caQTL')

    all_celltypes_eqtl = sorted(eqtl['celltype'].unique())
    all_celltypes_caqtl = sorted(caqtl['celltype'].unique())

    # Step 4: Match GWAS loci to xQTLs
    matched_eqtl = match_gwas_to_xqtl(leads, eqtl, window=args.match_window, label='eQTL')
    matched_caqtl = match_gwas_to_xqtl(leads, caqtl, window=args.match_window, label='caQTL')

    # Save matched results
    if len(matched_eqtl) > 0:
        matched_eqtl.to_csv(outdir / 'decode_ad_eqtl.csv', index=False)
        print(f"\nSaved eQTL matches to {outdir / 'decode_ad_eqtl.csv'}")

    if len(matched_caqtl) > 0:
        matched_caqtl.to_csv(outdir / 'decode_ad_caqtl.csv', index=False)
        print(f"Saved caQTL matches to {outdir / 'decode_ad_caqtl.csv'}")

    # Step 5: Build locus x cell type matrices
    if len(matched_eqtl) > 0:
        eqtl_matrix = build_locus_celltype_matrix(matched_eqtl, label='eQTL')
        eqtl_matrix.to_csv(outdir / 'decode_ad_locus_celltype_eqtl.csv')

    if len(matched_caqtl) > 0:
        caqtl_matrix = build_locus_celltype_matrix(matched_caqtl, label='caQTL')
        caqtl_matrix.to_csv(outdir / 'decode_ad_locus_celltype_caqtl.csv')

    # Step 6: DGSA comparison — AD vs non-AD genes
    xqtl_full = pd.read_csv(args.s6, usecols=['phenotype_id', 'variant_id', 'celltype',
                                                'slope', 'slope_se', 'pval_nominal', 'analysis'])
    dgsa = dgsa_ad_comparison(matched_eqtl, xqtl_full, args.dgsa, all_celltypes_eqtl)

    # Step 7: Summary
    print("\n" + "=" * 60)
    print("DECODE-AD SUMMARY")
    print("=" * 60)
    print(f"  GWAS SNPs (P < {args.p_threshold}): {len(gwas):,}")
    print(f"  Independent lead loci: {len(leads)}")
    print(f"  eQTL overlap: {matched_eqtl['locus_id'].nunique() if len(matched_eqtl) > 0 else 0} loci, "
          f"{matched_eqtl['phenotype_id'].nunique() if len(matched_eqtl) > 0 else 0} genes")
    print(f"  caQTL overlap: {matched_caqtl['locus_id'].nunique() if len(matched_caqtl) > 0 else 0} loci, "
          f"{matched_caqtl['phenotype_id'].nunique() if len(matched_caqtl) > 0 else 0} peaks")

    if len(matched_eqtl) > 0:
        # Cell type enrichment summary
        ct_counts = matched_eqtl.groupby('celltype')['phenotype_id'].nunique().sort_values(ascending=False)
        print(f"\n  Top 10 cell types by AD eQTL gene count:")
        for ct, count in ct_counts.head(10).items():
            print(f"    {ct:40s}  {count} genes")

    if len(matched_caqtl) > 0:
        ct_counts_ca = matched_caqtl.groupby('celltype')['phenotype_id'].nunique().sort_values(ascending=False)
        print(f"\n  Top 10 cell types by AD caQTL peak count:")
        for ct, count in ct_counts_ca.head(10).items():
            print(f"    {ct:40s}  {count} peaks")

    # Save lead loci with overlap annotation
    leads['has_eqtl'] = leads['locus_id'].isin(matched_eqtl['locus_id'].unique() if len(matched_eqtl) > 0 else [])
    leads['has_caqtl'] = leads['locus_id'].isin(matched_caqtl['locus_id'].unique() if len(matched_caqtl) > 0 else [])
    leads.to_csv(outdir / 'decode_ad_summary.csv', index=False)
    print(f"\nSaved lead loci summary to {outdir / 'decode_ad_summary.csv'}")
    print("=" * 60)


if __name__ == '__main__':
    main()
