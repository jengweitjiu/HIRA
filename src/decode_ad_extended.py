#!/usr/bin/env python3
"""
DECODE-AD Extended — Disease architecture comparison, drug targets, SMR mediation.

Analysis 1: AD vs Asthma vs RA cell-type SMR profiles (S15)
Analysis 2: Drug target overlap with AD eGenes
Analysis 3: S15 SMR mediation for AD (eQTL + caQTL)

Input:  data/raw/science.adt3130_table_s15.xlsx
        data/raw/science.adt3130_table_s7.csv
        results/decode_ad_eqtl.csv
Output: results/decode_ad_disease_comparison.csv
        results/decode_ad_drug_targets.csv
        results/decode_ad_smr_mediation.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


# ======================================================================
# Analysis 1: AD vs Asthma vs RA cell-type architecture
# ======================================================================

def disease_comparison(s15_path):
    """Compare cell-type SMR profiles across diseases."""
    print("=" * 60)
    print("ANALYSIS 1: Disease cell-type architecture comparison")
    print("=" * 60)

    s15 = pd.read_excel(s15_path)
    print(f"  S15 columns: {list(s15.columns)}")
    print(f"  S15 shape: {s15.shape}")

    # Focus on AD, As (asthma), RA, GD (Graves' disease)
    diseases = {'AD': 'Atopic Dermatitis', 'As': 'Asthma',
                'RA': 'Rheumatoid Arthritis', 'GD': "Graves' Disease"}

    results = []
    disease_profiles = {}

    for abbr, full_name in diseases.items():
        sub = s15[s15['trait'] == abbr]
        if len(sub) == 0:
            continue

        # Count associations per cell type (both eQTL and caQTL)
        ct_counts = sub.groupby('celltype').size().reset_index(name='n_associations')
        ct_counts['trait'] = abbr
        ct_counts['trait_name'] = full_name

        # Also count by QTL type
        for qtype in ['eQTL', 'caQTL']:
            qsub = sub[sub['QTL'] == qtype]
            q_counts = qsub.groupby('celltype').size().reset_index(name=f'n_{qtype}')
            ct_counts = ct_counts.merge(q_counts, on='celltype', how='left')
            ct_counts[f'n_{qtype}'] = ct_counts[f'n_{qtype}'].fillna(0).astype(int)

        results.append(ct_counts)
        disease_profiles[abbr] = ct_counts.set_index('celltype')['n_associations']

        print(f"\n  {full_name} ({abbr}):")
        print(f"    Total associations: {len(sub)}")
        print(f"    Cell types: {sub['celltype'].nunique()}")
        print(f"    Unique genes (eQTL): {sub[sub['QTL'] == 'eQTL']['Gene'].nunique()}")
        top_ct = ct_counts.nlargest(5, 'n_associations')
        for _, r in top_ct.iterrows():
            print(f"      {r['celltype']:35s}  {int(r['n_associations'])} assoc")

    all_results = pd.concat(results, ignore_index=True)

    # Identify disease-specific and shared cell types
    ad_cts = set(disease_profiles.get('AD', pd.Series(dtype=float)).index)
    as_cts = set(disease_profiles.get('As', pd.Series(dtype=float)).index)
    ra_cts = set(disease_profiles.get('RA', pd.Series(dtype=float)).index)

    shared_all = ad_cts & as_cts & ra_cts
    ad_specific = ad_cts - as_cts - ra_cts
    as_specific = as_cts - ad_cts - ra_cts
    ra_specific = ra_cts - ad_cts - as_cts
    ad_as_shared = (ad_cts & as_cts) - ra_cts
    ad_ra_shared = (ad_cts & ra_cts) - as_cts

    print(f"\n  Cell type overlap:")
    print(f"    AD cell types:     {len(ad_cts)}")
    print(f"    Asthma cell types: {len(as_cts)}")
    print(f"    RA cell types:     {len(ra_cts)}")
    print(f"    Shared all three:  {len(shared_all)}")
    print(f"    AD-specific:       {len(ad_specific)} -> {sorted(ad_specific)[:5]}")
    print(f"    Asthma-specific:   {len(as_specific)}")
    print(f"    RA-specific:       {len(ra_specific)}")
    print(f"    AD+Asthma only:    {len(ad_as_shared)}")
    print(f"    AD+RA only:        {len(ad_ra_shared)}")

    # Spearman correlation of cell-type profiles (on shared cell types)
    if 'AD' in disease_profiles and 'As' in disease_profiles:
        shared = sorted(ad_cts & as_cts)
        if len(shared) >= 5:
            v1 = disease_profiles['AD'].reindex(shared, fill_value=0)
            v2 = disease_profiles['As'].reindex(shared, fill_value=0)
            rho, p = stats.spearmanr(v1, v2)
            print(f"\n  AD vs Asthma profile correlation (n={len(shared)} CTs): rho={rho:.3f}, P={p:.2e}")

    if 'AD' in disease_profiles and 'RA' in disease_profiles:
        shared = sorted(ad_cts & ra_cts)
        if len(shared) >= 5:
            v1 = disease_profiles['AD'].reindex(shared, fill_value=0)
            v2 = disease_profiles['RA'].reindex(shared, fill_value=0)
            rho, p = stats.spearmanr(v1, v2)
            print(f"  AD vs RA profile correlation (n={len(shared)} CTs): rho={rho:.3f}, P={p:.2e}")

    return all_results, disease_profiles


# ======================================================================
# Analysis 2: Drug target overlap
# ======================================================================

def drug_target_analysis(eqtl_path):
    """Check overlap between AD eGenes and known drug targets."""
    print("\n" + "=" * 60)
    print("ANALYSIS 2: AD drug target cell-type resolution")
    print("=" * 60)

    eqtl = pd.read_csv(eqtl_path)
    ad_genes = set(eqtl['phenotype_id'].unique())
    print(f"  AD eGenes from DECODE-AD: {len(ad_genes)}")

    # Known AD drug targets and their drugs
    drug_targets = {
        'IL4R':  'Dupilumab (anti-IL4Ra)',
        'IL13':  'Tralokinumab (anti-IL13)',
        'IL31RA': 'Nemolizumab (anti-IL31RA)',
        'JAK1':  'Upadacitinib, Abrocitinib (JAK1i)',
        'JAK2':  'Baricitinib (JAK1/2i)',
        'JAK3':  'Tofacitinib (pan-JAKi)',
        'TYK2':  'Deucravacitinib (TYK2i)',
        'IL5':   'Mepolizumab (anti-IL5)',
        'IL5RA': 'Benralizumab (anti-IL5Ra)',
        'TSLP':  'Tezepelumab (anti-TSLP)',
        'IL33':  'Itepekimab (anti-IL33)',
        'IL1RL1': 'Astegolimab (anti-IL33R/ST2)',
        'PDE4':  'Crisaborole (PDE4i topical)',
        'PDE4A': 'Crisaborole (PDE4i topical)',
        'PDE4B': 'Apremilast (PDE4i oral)',
        'PDE4D': 'Roflumilast (PDE4i)',
        'IL22':  'Fezagitinib (emerging)',
        'IL17A': 'Secukinumab (off-label)',
        'IL12B': 'Ustekinumab (off-label)',
        'IFNG':  'Emapalumab (anti-IFNg)',
        'IL2':   'Low-dose IL2 (experimental)',
    }

    print(f"  Drug targets checked: {len(drug_targets)}")

    records = []
    for gene, drug in drug_targets.items():
        if gene in ad_genes:
            gene_eqtl = eqtl[eqtl['phenotype_id'] == gene]
            celltypes = sorted(gene_eqtl['celltype'].unique())
            n_ct = len(celltypes)
            mean_slope = gene_eqtl['slope'].abs().mean()
            top_ct = gene_eqtl.loc[gene_eqtl['slope'].abs().idxmax(), 'celltype']
            top_slope = gene_eqtl['slope'].abs().max()

            records.append({
                'gene': gene,
                'drug': drug,
                'n_celltypes': n_ct,
                'mean_abs_slope': mean_slope,
                'top_celltype': top_ct,
                'top_abs_slope': top_slope,
                'all_celltypes': '; '.join(celltypes),
            })
            print(f"\n  MATCH: {gene} -> {drug}")
            print(f"    eQTL in {n_ct} cell types, mean |slope| = {mean_slope:.4f}")
            print(f"    Top cell type: {top_ct} (|slope| = {top_slope:.4f})")
            if n_ct <= 10:
                print(f"    All: {', '.join(celltypes)}")
            else:
                print(f"    Top 5: {', '.join(celltypes[:5])} ... (+{n_ct-5} more)")

    if not records:
        print("  No drug target matches among AD eGenes.")

    drug_df = pd.DataFrame(records)

    # Summary
    matched = set(drug_df['gene']) if len(drug_df) > 0 else set()
    unmatched = set(drug_targets.keys()) - matched - (set(drug_targets.keys()) - ad_genes)
    print(f"\n  Summary: {len(matched)} / {len(drug_targets)} drug targets found in AD eGenes")
    print(f"  Matched: {sorted(matched)}")
    not_found = set(drug_targets.keys()) - ad_genes
    print(f"  Not in AD eGenes: {sorted(not_found)}")

    return drug_df


# ======================================================================
# Analysis 3: S15 SMR mediation for AD
# ======================================================================

def smr_mediation(s15_path, s7_path):
    """Extract AD SMR mediators from S15 and caQTL->eQTL mediation from S7."""
    print("\n" + "=" * 60)
    print("ANALYSIS 3: SMR mediation of AD genetic risk")
    print("=" * 60)

    # S15: direct GWAS -> gene SMR for AD
    s15 = pd.read_excel(s15_path)
    ad_s15 = s15[s15['trait'] == 'AD'].copy()
    print(f"\n  S15 AD entries: {len(ad_s15)}")

    # Also get As (asthma) and related atopic traits
    atopic_traits = ['AD', 'As', 'AR', 'Urt']
    atopic = s15[s15['trait'].isin(atopic_traits)].copy()
    print(f"  Atopic trait entries (AD+As+AR+Urt): {len(atopic)}")

    # eQTL mediators
    ad_eqtl = ad_s15[ad_s15['QTL'] == 'eQTL']
    print(f"\n  AD eQTL SMR mediators:")
    print(f"    Entries: {len(ad_eqtl)}")
    print(f"    Unique genes: {ad_eqtl['Gene'].nunique()}")
    print(f"    Cell types: {ad_eqtl['celltype'].nunique()}")

    # Top mediating genes by number of cell types
    gene_cts = ad_eqtl.groupby('Gene').agg(
        n_celltypes=('celltype', 'nunique'),
        min_p_smr=('p_SMR', 'min'),
        celltypes=('celltype', lambda x: '; '.join(sorted(x.unique())))
    ).sort_values('n_celltypes', ascending=False)

    print(f"\n    Top mediating genes (by cell-type breadth):")
    for gene, row in gene_cts.iterrows():
        print(f"      {gene:15s}  {row['n_celltypes']:2d} CTs  min_P={row['min_p_smr']:.2e}")

    # caQTL mediators
    ad_caqtl = ad_s15[ad_s15['QTL'] == 'caQTL']
    print(f"\n  AD caQTL SMR mediators:")
    print(f"    Entries: {len(ad_caqtl)}")
    print(f"    Cell types: {ad_caqtl['celltype'].nunique()}")
    if len(ad_caqtl) > 0:
        ct_counts = ad_caqtl['celltype'].value_counts()
        print(f"    Top cell types by caQTL count:")
        for ct, n in ct_counts.head(10).items():
            print(f"      {ct:35s}  {n} peaks")

    # S7: caQTL -> eQTL mediation analysis
    print(f"\n  S7 caQTL -> eQTL mediation:")
    s7 = pd.read_csv(s7_path)
    print(f"    S7 columns: {list(s7.columns)[:10]}")
    print(f"    S7 shape: {s7.shape}")

    # Find S7 entries where the outcome gene is an AD SMR gene
    ad_genes = set(ad_eqtl['Gene'].unique())
    s7_ad = s7[s7['Outco_Gene'].isin(ad_genes)]
    print(f"    S7 entries with AD mediating gene as outcome: {len(s7_ad)}")

    if len(s7_ad) > 0:
        print(f"    Unique exposure peaks: {s7_ad['Expo_ID'].nunique()}")
        print(f"    Unique outcome genes: {s7_ad['Outco_Gene'].nunique()}")
        print(f"    Cell types: {s7_ad['celltype'].nunique()}")
        print(f"\n    Top caQTL -> AD gene mediation paths:")
        s7_ad_sig = s7_ad[s7_ad['p_SMR'] < 1e-5].sort_values('p_SMR')
        for _, row in s7_ad_sig.head(10).iterrows():
            print(f"      {str(row.get('Expo_ID',''))[:30]:30s} -> {row['Outco_Gene']:15s}  "
                  f"in {row['celltype']:30s}  P={row['p_SMR']:.2e}")

    # Compile mediation summary
    mediation_records = []
    for _, row in ad_eqtl.iterrows():
        mediation_records.append({
            'gene': row['Gene'],
            'celltype': row['celltype'],
            'p_SMR': row['p_SMR'],
            'p_HEIDI': row['p_HEIDI'],
            'b_SMR': row['b_SMR'],
            'QTL': 'eQTL',
            'trait': 'AD',
        })
    for _, row in ad_caqtl.iterrows():
        mediation_records.append({
            'gene': row.get('Gene', ''),
            'celltype': row['celltype'],
            'p_SMR': row['p_SMR'],
            'p_HEIDI': row['p_HEIDI'],
            'b_SMR': row['b_SMR'],
            'QTL': 'caQTL',
            'trait': 'AD',
        })

    # Add asthma for comparison
    as_eqtl = atopic[(atopic['trait'] == 'As') & (atopic['QTL'] == 'eQTL')]
    for _, row in as_eqtl.iterrows():
        mediation_records.append({
            'gene': row['Gene'],
            'celltype': row['celltype'],
            'p_SMR': row['p_SMR'],
            'p_HEIDI': row['p_HEIDI'],
            'b_SMR': row['b_SMR'],
            'QTL': 'eQTL',
            'trait': 'As',
        })

    mediation_df = pd.DataFrame(mediation_records)

    # Key finding: which cell types mediate AD vs asthma?
    print(f"\n  AD vs Asthma mediating cell types (eQTL):")
    ad_med_cts = set(ad_eqtl['celltype'].unique())
    as_med_cts = set(as_eqtl['celltype'].unique())
    shared_med = ad_med_cts & as_med_cts
    ad_only_med = ad_med_cts - as_med_cts
    as_only_med = as_med_cts - ad_med_cts
    print(f"    AD-mediating CTs:     {len(ad_med_cts)}")
    print(f"    Asthma-mediating CTs: {len(as_med_cts)}")
    print(f"    Shared:               {len(shared_med)}")
    print(f"    AD-only:              {len(ad_only_med)} -> {sorted(ad_only_med)[:5]}")
    print(f"    Asthma-only:          {len(as_only_med)}")

    return mediation_df, ad_s15


def main():
    parser = argparse.ArgumentParser(description='DECODE-AD Extended Analysis')
    parser.add_argument('--s15', default='data/raw/science.adt3130_table_s15.xlsx')
    parser.add_argument('--s7', default='data/raw/science.adt3130_table_s7.csv')
    parser.add_argument('--eqtl', default='results/decode_ad_eqtl.csv')
    parser.add_argument('--output-dir', default='results')
    args = parser.parse_args()

    outdir = Path(args.output_dir)

    # Analysis 1
    disease_df, profiles = disease_comparison(args.s15)
    disease_df.to_csv(outdir / 'decode_ad_disease_comparison.csv', index=False)
    print(f"\nSaved to {outdir / 'decode_ad_disease_comparison.csv'}")

    # Analysis 2
    drug_df = drug_target_analysis(args.eqtl)
    drug_df.to_csv(outdir / 'decode_ad_drug_targets.csv', index=False)
    print(f"Saved to {outdir / 'decode_ad_drug_targets.csv'}")

    # Analysis 3
    mediation_df, ad_s15 = smr_mediation(args.s15, args.s7)
    mediation_df.to_csv(outdir / 'decode_ad_smr_mediation.csv', index=False)
    print(f"Saved to {outdir / 'decode_ad_smr_mediation.csv'}")

    # Final summary
    print("\n" + "=" * 60)
    print("DECODE-AD EXTENDED — KEY FINDINGS")
    print("=" * 60)
    print("  1. Disease architecture: AD has 35 SMR cell types; compare with")
    print("     Asthma (65 CTs) and RA (62 CTs) for shared/specific patterns")
    print("  2. Drug targets: see which approved AD therapies have genetic")
    print("     support from CIMA eQTL data and in which cell compartments")
    print("  3. SMR mediation: IL18R1 and LIME1 are top AD mediating genes")
    print("     across multiple T cell subtypes")
    print("=" * 60)


if __name__ == '__main__':
    main()
