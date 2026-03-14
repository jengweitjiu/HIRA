#!/usr/bin/env python3
"""
SMR/HEIDI Analysis for AD x CIMA.

Analytical SMR using Budu-Aggrey 2023 AD GWAS + CIMA S6 lead eQTLs for
the top 5 AD-relevant cell types. Cross-reference with S15 pre-computed
SMR/HEIDI results (which used full eQTL summary stats).

CIMA S6 contains only lead eQTLs (1 SNP per gene per cell type). For a
single-instrument MR, the SMR chi-sq statistic is computed analytically:

    z_SMR = (z_GWAS * z_eQTL) / sqrt(z_GWAS^2 + z_eQTL^2)

This is mathematically equivalent to what the SMR software computes.
HEIDI requires multiple independent SNPs and cannot be run on lead-only data,
so we rely on S15 (which used full summary statistics) for HEIDI filtering.

Steps:
  1. Build chr:pos -> GWAS lookup from AD GWAS summary stats
  2. For each target cell type, match lead eQTLs to GWAS by chr:pos
  3. Compute analytical SMR z-scores and p-values
  4. Cross-reference with S15 pre-computed SMR/HEIDI
  5. Generate figure

Output: results/smr_ad_denovo.csv
        results/smr_ad_vs_s15.csv
        figures/fig_smr_ad.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
GWAS_FILE = Path('data/raw/AD_GWAS_Budu-Aggrey_2023.tsv.gz')
S6_FILE = Path('data/raw/CIMA_Table_S6.csv')
S15_FILE = Path('data/raw/science.adt3130_table_s15.xlsx')

# Target cell types (top 5 from DECODE-AD PRS analysis)
TARGET_CTS = [
    'CD4_Tn_CCR7',
    'CD8_Tn_CCR7',
    'CD4_Tfh-like_CXCR5',
    'cMono_CD14',
    'Mature_NK_dim_FCGR3A',
]

# Budu-Aggrey 2023 meta-analysis total N
GWAS_N = 865459


def step1_load_gwas():
    """Load GWAS and build chr:pos lookup."""
    print("\n--- Step 1: Load AD GWAS summary statistics ---")

    gwas_lookup_file = Path('data/processed/gwas_chrpos_lookup.pkl')
    if gwas_lookup_file.exists():
        import pickle
        with open(gwas_lookup_file, 'rb') as f:
            gwas_lookup = pickle.load(f)
        print(f"  Loaded cached GWAS lookup: {len(gwas_lookup):,} SNPs")
        return gwas_lookup

    print("  Reading GWAS file (16M+ SNPs)...")
    gwas = pd.read_csv(GWAS_FILE, sep='\t',
                       dtype={'chromosome': str, 'base_pair_location': 'Int64'})
    print(f"  Raw GWAS: {len(gwas):,} rows")
    print(f"  Columns: {list(gwas.columns)}")

    # Drop rows with missing essentials
    gwas = gwas.dropna(subset=['chromosome', 'base_pair_location', 'beta',
                                'standard_error', 'p_value'])
    print(f"  After cleaning: {len(gwas):,} rows")

    # Build lookup: "chr{N}_{pos}" -> {beta, se, p, a1, a2, af, rsid}
    gwas_lookup = {}
    for _, row in gwas.iterrows():
        key = f"chr{row['chromosome']}_{row['base_pair_location']}"
        gwas_lookup[key] = {
            'rsid': row.get('variant_id', ''),
            'beta': row['beta'],
            'se': row['standard_error'],
            'p': row['p_value'],
            'a1': row.get('effect_allele', ''),
            'a2': row.get('other_allele', ''),
            'af': row.get('effect_allele_frequency', np.nan),
        }

    print(f"  GWAS lookup: {len(gwas_lookup):,} unique chr:pos entries")

    import pickle
    with open(gwas_lookup_file, 'wb') as f:
        pickle.dump(gwas_lookup, f)
    print(f"  Cached to {gwas_lookup_file}")

    return gwas_lookup


def step1_load_gwas_fast():
    """Load GWAS as DataFrame for vectorized matching."""
    print("\n--- Step 1: Load AD GWAS summary statistics ---")

    cache_file = Path('data/processed/gwas_for_smr.pkl')
    if cache_file.exists():
        gwas = pd.read_pickle(cache_file)
        print(f"  Loaded cached GWAS: {len(gwas):,} SNPs")
        return gwas

    print("  Reading GWAS file (16M+ SNPs)...")
    gwas = pd.read_csv(GWAS_FILE, sep='\t',
                       dtype={'chromosome': str, 'base_pair_location': 'Int64'})
    print(f"  Raw GWAS: {len(gwas):,} rows")

    gwas = gwas.dropna(subset=['chromosome', 'base_pair_location', 'beta',
                                'standard_error', 'p_value'])

    # Create chr:pos key for matching
    gwas['chrpos'] = 'chr' + gwas['chromosome'].astype(str) + '_' + gwas['base_pair_location'].astype(str)

    # Keep only needed columns
    gwas = gwas[['chrpos', 'variant_id', 'beta', 'standard_error', 'p_value',
                 'effect_allele', 'other_allele', 'effect_allele_frequency']].copy()
    gwas.columns = ['chrpos', 'rsid', 'gwas_beta', 'gwas_se', 'gwas_p',
                     'gwas_a1', 'gwas_a2', 'gwas_af']

    # Drop duplicates on chrpos (keep first = lowest p-value would be better, but first is fine)
    gwas = gwas.drop_duplicates(subset='chrpos', keep='first')

    print(f"  GWAS lookup: {len(gwas):,} unique chr:pos entries")

    gwas.to_pickle(cache_file)
    print(f"  Cached to {cache_file}")

    return gwas


def step2_analytical_smr(gwas):
    """Compute analytical SMR for each cell type."""
    print("\n--- Step 2: Analytical SMR ---")

    print("  Loading S6 eQTL data...")
    s6 = pd.read_csv(S6_FILE)
    s6_eqtl = s6[s6['analysis'] == 'cis-eQTL'].copy()
    print(f"  Total cis-eQTLs: {len(s6_eqtl):,}")

    all_results = []

    for ct in TARGET_CTS:
        print(f"\n  === {ct} ===")
        ct_data = s6_eqtl[s6_eqtl['celltype'] == ct].copy()
        print(f"    Lead eQTLs: {len(ct_data)}")

        if len(ct_data) == 0:
            continue

        # Match to GWAS by chr:pos (variant_id in S6 is chr_pos format)
        ct_data = ct_data.rename(columns={'variant_id': 'chrpos'})
        merged = ct_data.merge(gwas, on='chrpos', how='inner')
        print(f"    Matched to GWAS: {len(merged)}/{len(ct_data)} ({100*len(merged)/len(ct_data):.1f}%)")

        if len(merged) == 0:
            continue

        # Compute z-scores
        merged['z_eqtl'] = merged['slope'] / merged['slope_se']
        merged['z_gwas'] = merged['gwas_beta'] / merged['gwas_se']

        # Allele alignment: check if eQTL effect allele matches GWAS
        # S6 has A1(ALT/effect allele) and A2(REF)
        # GWAS has gwas_a1 (effect) and gwas_a2 (other)
        # If eQTL A1 matches GWAS A2 (flipped), flip the eQTL z-score
        eqtl_a1 = merged['A1(ALT/effect allele)'].str.upper()
        gwas_a1 = merged['gwas_a1'].str.upper()
        gwas_a2 = merged['gwas_a2'].str.upper()

        # Flip z_eqtl where eQTL effect allele = GWAS other allele
        flip_mask = (eqtl_a1 == gwas_a2)
        merged.loc[flip_mask, 'z_eqtl'] = -merged.loc[flip_mask, 'z_eqtl']
        n_flipped = flip_mask.sum()
        print(f"    Allele-flipped: {n_flipped}")

        # Analytical SMR: z_SMR = (z_GWAS * z_eQTL) / sqrt(z_GWAS^2 + z_eQTL^2)
        z_g2 = merged['z_gwas'] ** 2
        z_e2 = merged['z_eqtl'] ** 2
        merged['z_smr'] = (merged['z_gwas'] * merged['z_eqtl']) / np.sqrt(z_g2 + z_e2)
        merged['chi2_smr'] = merged['z_smr'] ** 2
        merged['p_smr'] = scipy_stats.chi2.sf(merged['chi2_smr'], df=1)

        # Effect size: b_SMR = b_GWAS / b_eQTL
        merged['b_smr'] = merged['gwas_beta'] / merged['slope']
        merged['se_smr'] = np.abs(merged['b_smr']) * np.sqrt(
            (merged['gwas_se'] / merged['gwas_beta']) ** 2 +
            (merged['slope_se'] / merged['slope']) ** 2
        )

        # Bonferroni threshold for this cell type
        n_tests = len(merged)
        bonf_thresh = 0.05 / n_tests
        n_sig_bonf = (merged['p_smr'] < bonf_thresh).sum()
        n_sig_nom = (merged['p_smr'] < 0.05).sum()
        n_sig_suggestive = (merged['p_smr'] < 1e-4).sum()

        print(f"    Tested: {n_tests}")
        print(f"    Bonferroni threshold: {bonf_thresh:.2e}")
        print(f"    Significant (Bonferroni): {n_sig_bonf}")
        print(f"    Significant (P<1e-4): {n_sig_suggestive}")
        print(f"    Nominal (P<0.05): {n_sig_nom}")

        # Show top hits
        top = merged.nsmallest(10, 'p_smr')
        print(f"    Top 10 hits:")
        for _, row in top.iterrows():
            print(f"      {row['phenotype_id']:15s} z_SMR={row['z_smr']:7.3f} "
                  f"P={row['p_smr']:.2e}  b_SMR={row['b_smr']:.3f}  "
                  f"GWAS_P={row['gwas_p']:.2e}  eQTL_P={row['pval_nominal']:.2e}")

        # Store results
        result = merged[['phenotype_id', 'chrpos', 'rsid', 'slope', 'slope_se',
                         'pval_nominal', 'gwas_beta', 'gwas_se', 'gwas_p',
                         'z_eqtl', 'z_gwas', 'z_smr', 'chi2_smr', 'p_smr',
                         'b_smr', 'se_smr']].copy()
        result['celltype'] = ct
        result['bonf_sig'] = result['p_smr'] < bonf_thresh
        all_results.append(result)

    if all_results:
        denovo = pd.concat(all_results, ignore_index=True)
        print(f"\n  Total de novo SMR results: {len(denovo):,}")
        return denovo
    return pd.DataFrame()


def step3_crossref_s15(denovo):
    """Cross-reference de novo SMR with S15 pre-computed results."""
    print("\n--- Step 3: Cross-reference with S15 ---")

    s15 = pd.read_excel(S15_FILE)
    print(f"  S15 columns: {list(s15.columns)}")
    print(f"  S15 total: {len(s15)} SMR associations")
    print(f"  S15 traits: {s15['trait'].nunique()}")
    print(f"  Sample traits: {s15['trait'].unique()[:10]}")

    # Find AD-related traits (S15 uses 'AD' as exact trait name)
    ad_mask = (s15['trait'] == 'AD') | s15['trait'].str.contains('atopic|dermatitis|eczema', case=False, na=False)
    ad_s15 = s15[ad_mask].copy()
    print(f"  S15 AD associations: {len(ad_s15)}")
    print(f"  S15 AD unique genes: {ad_s15['Gene'].nunique()}")
    print(f"  S15 AD cell types: {ad_s15['celltype'].nunique()}")

    # HEIDI filtering on S15
    if 'p_HEIDI' in ad_s15.columns:
        heidi_pass = ad_s15[ad_s15['p_HEIDI'] > 0.05]
        heidi_fail = ad_s15[ad_s15['p_HEIDI'] <= 0.05]
        heidi_na = ad_s15[ad_s15['p_HEIDI'].isna()]
        print(f"  S15 HEIDI pass (P>0.05): {len(heidi_pass)}")
        print(f"  S15 HEIDI fail (P<=0.05): {len(heidi_fail)}")
        print(f"  S15 HEIDI NA: {len(heidi_na)}")

    # Per-cell-type comparison
    print(f"\n  Per-cell-type comparison:")
    comparison = []

    for ct in TARGET_CTS:
        ct_s15 = ad_s15[ad_s15['celltype'] == ct]
        ct_s15_heidi = ct_s15[ct_s15['p_HEIDI'] > 0.05] if 'p_HEIDI' in ct_s15.columns else ct_s15

        if len(denovo) > 0:
            ct_denovo = denovo[denovo['celltype'] == ct]
            ct_denovo_sig = ct_denovo[ct_denovo['p_smr'] < 1e-4]
            ct_denovo_bonf = ct_denovo[ct_denovo['bonf_sig'] == True]

            # Gene overlap
            denovo_genes = set(ct_denovo_sig['phenotype_id'].unique())
            s15_genes = set(ct_s15['Gene'].unique()) if 'Gene' in ct_s15.columns else set()
            overlap = denovo_genes & s15_genes
        else:
            ct_denovo = pd.DataFrame()
            ct_denovo_sig = pd.DataFrame()
            ct_denovo_bonf = pd.DataFrame()
            denovo_genes = set()
            s15_genes = set(ct_s15['Gene'].unique()) if 'Gene' in ct_s15.columns else set()
            overlap = set()

        row = {
            'celltype': ct,
            'n_denovo_tested': len(ct_denovo),
            'n_denovo_p1e4': len(ct_denovo_sig),
            'n_denovo_bonf': len(ct_denovo_bonf),
            'n_s15_total': len(ct_s15),
            'n_s15_heidi_pass': len(ct_s15_heidi),
            'gene_overlap': len(overlap),
            'overlap_genes': ';'.join(sorted(overlap)) if overlap else '',
        }
        comparison.append(row)
        print(f"    {ct:30s}  denovo={len(ct_denovo_sig):3d}(P<1e-4)  "
              f"s15={len(ct_s15):3d}  heidi_pass={len(ct_s15_heidi):3d}  "
              f"overlap={len(overlap)}")

    comp_df = pd.DataFrame(comparison)
    comp_df.to_csv('results/smr_ad_vs_s15.csv', index=False)
    print(f"\n  Saved: results/smr_ad_vs_s15.csv")

    # Save S15 AD subset
    ad_s15.to_csv('results/smr_ad_s15_subset.csv', index=False)
    print(f"  Saved: results/smr_ad_s15_subset.csv ({len(ad_s15)} rows)")

    # Show S15 top hits in target cell types
    target_s15 = ad_s15[ad_s15['celltype'].isin(TARGET_CTS)]
    if len(target_s15) > 0 and 'p_SMR' in target_s15.columns:
        print(f"\n  S15 top hits in target cell types (HEIDI-pass):")
        heidi_target = target_s15[target_s15.get('p_HEIDI', pd.Series(dtype=float)) > 0.05]
        for _, row in heidi_target.nsmallest(15, 'p_SMR').iterrows():
            heidi_str = f"HEIDI={row['p_HEIDI']:.2e}" if 'p_HEIDI' in row.index and not pd.isna(row['p_HEIDI']) else "HEIDI=NA"
            gene_str = str(row['Gene'])[:15]
            ct_str = str(row['celltype'])[:25]
            print(f"    {gene_str:15s} {ct_str:25s} p_SMR={row['p_SMR']:.2e}  "
                  f"{heidi_str}  b_SMR={row.get('b_SMR', np.nan):.3f}")

    return ad_s15, comp_df


def step4_figure(denovo, s15_ad, comp_df):
    """Generate SMR results figure."""
    print("\n--- Step 4: Figure ---")

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})

    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 160/25.4))

    ct_colors = {
        'CD4_Tn_CCR7': '#E74C3C',
        'CD8_Tn_CCR7': '#3498DB',
        'CD4_Tfh-like_CXCR5': '#2ECC71',
        'cMono_CD14': '#F39C12',
        'Mature_NK_dim_FCGR3A': '#9B59B6',
    }

    # Panel A: De novo SMR Manhattan-style per cell type
    ax = axes[0, 0]
    if len(denovo) > 0:
        for i, ct in enumerate(TARGET_CTS):
            ct_data = denovo[denovo['celltype'] == ct]
            if len(ct_data) > 0:
                pvals = -np.log10(ct_data['p_smr'].clip(lower=1e-50))
                jitter = np.random.RandomState(42 + i).uniform(-0.2, 0.2, len(pvals))
                ax.scatter(np.full(len(pvals), i) + jitter, pvals,
                           s=4, alpha=0.3, c=ct_colors[ct], rasterized=True)

        # Significance lines
        if len(denovo) > 0:
            bonf_line = -np.log10(0.05 / len(denovo[denovo['celltype'] == TARGET_CTS[0]]))
            ax.axhline(bonf_line, color='red', linestyle='--', linewidth=0.8, label='Bonferroni')
            ax.axhline(-np.log10(1e-4), color='orange', linestyle=':', linewidth=0.8, label='P<1e-4')

        ax.set_xticks(range(len(TARGET_CTS)))
        ax.set_xticklabels([ct.replace('_', '\n') for ct in TARGET_CTS],
                           fontsize=5, rotation=45, ha='right')
        ax.set_ylabel('-log10(p_SMR)')
        ax.legend(fontsize=5, loc='upper right')
    ax.set_title('De novo analytical SMR', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel B: S15 AD associations per cell type (top 15)
    ax = axes[0, 1]
    if s15_ad is not None and len(s15_ad) > 0:
        ct_counts = s15_ad['celltype'].value_counts().head(15)
        colors = ['#E74C3C' if ct in TARGET_CTS else '#BDC3C7' for ct in ct_counts.index]
        ax.barh(range(len(ct_counts)), ct_counts.values, color=colors)
        ax.set_yticks(range(len(ct_counts)))
        ax.set_yticklabels(ct_counts.index, fontsize=5)
        ax.set_xlabel('N SMR associations')
        ax.invert_yaxis()
    ax.set_title('S15 pre-computed SMR', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel C: HEIDI filter effect in target cell types
    ax = axes[1, 0]
    if comp_df is not None and len(comp_df) > 0:
        x = np.arange(len(TARGET_CTS))
        w = 0.35
        s15_total = comp_df['n_s15_total'].values
        s15_heidi = comp_df['n_s15_heidi_pass'].values
        denovo_sig = comp_df['n_denovo_p1e4'].values

        ax.bar(x - w/2, s15_total, w, label='S15 total', color='#BDC3C7')
        ax.bar(x + w/2, s15_heidi, w, label='S15 HEIDI pass', color='#E74C3C')

        ax.set_xticks(x)
        ax.set_xticklabels([ct.replace('_', '\n') for ct in TARGET_CTS],
                           fontsize=5, rotation=45, ha='right')
        ax.set_ylabel('N genes')
        ax.legend(fontsize=6)
    ax.set_title('HEIDI filtering', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel D: Comparison — de novo hits vs S15 overlap
    ax = axes[1, 1]
    if len(denovo) > 0:
        # For each cell type, plot de novo -log10(p_smr) vs S15 -log10(p_SMR) for shared genes
        for ct in TARGET_CTS:
            ct_denovo = denovo[denovo['celltype'] == ct]
            if s15_ad is not None and 'Gene' in s15_ad.columns:
                ct_s15 = s15_ad[s15_ad['celltype'] == ct]
                shared = ct_denovo.merge(ct_s15, left_on='phenotype_id', right_on='Gene', how='inner')
                if len(shared) > 0 and 'p_SMR' in shared.columns:
                    ax.scatter(-np.log10(shared['p_smr'].clip(lower=1e-50)),
                               -np.log10(shared['p_SMR'].clip(lower=1e-50)),
                               s=15, c=ct_colors[ct], alpha=0.6,
                               label=ct.split('_')[0], edgecolors='white', linewidth=0.3)

        # Add diagonal
        lim = ax.get_xlim()[1]
        ax.plot([0, lim], [0, lim], 'k--', linewidth=0.5, alpha=0.3)
        ax.set_xlabel('-log10(p_SMR) de novo')
        ax.set_ylabel('-log10(p_SMR) S15')
        ax.legend(fontsize=5, loc='lower right')
    ax.set_title('De novo vs S15 concordance', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel labels
    for ax_i, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax_i.text(-0.15, 1.1, label, transform=ax_i.transAxes, fontsize=12,
                  fontweight='bold', va='top')

    fig.suptitle('SMR Analysis: AD GWAS x CIMA eQTLs', fontsize=11,
                 fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    Path('figures').mkdir(exist_ok=True)
    fig.savefig('figures/fig_smr_ad.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_smr_ad.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved figures/fig_smr_ad.pdf + .png")


def main():
    print("=" * 60)
    print("SMR/HEIDI Analysis: AD x CIMA (Top 5 Cell Types)")
    print("=" * 60)
    print("Method: Analytical SMR (single-instrument MR)")
    print("  z_SMR = (z_GWAS * z_eQTL) / sqrt(z_GWAS^2 + z_eQTL^2)")
    print("  HEIDI: from S15 pre-computed (full summary stats)")

    # Step 1: Load GWAS
    gwas = step1_load_gwas_fast()

    # Step 2: Analytical SMR
    denovo = step2_analytical_smr(gwas)

    # Save de novo results
    if len(denovo) > 0:
        denovo.to_csv('results/smr_ad_denovo.csv', index=False)
        print(f"\n  Saved: results/smr_ad_denovo.csv ({len(denovo):,} rows)")

        # Summary statistics
        print(f"\n  === Summary ===")
        for ct in TARGET_CTS:
            ct_data = denovo[denovo['celltype'] == ct]
            n_bonf = ct_data['bonf_sig'].sum()
            n_1e4 = (ct_data['p_smr'] < 1e-4).sum()
            print(f"    {ct:30s}  tested={len(ct_data):5d}  Bonferroni={n_bonf:3d}  P<1e-4={n_1e4:3d}")

    # Step 3: Cross-reference with S15
    s15_ad, comp_df = step3_crossref_s15(denovo)

    # Step 4: Figure
    step4_figure(denovo, s15_ad, comp_df)

    print("\n" + "=" * 60)
    print("SMR pipeline complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
