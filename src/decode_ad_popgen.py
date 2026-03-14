#!/usr/bin/env python3
"""
DECODE-AD Population Genetics & Enhancer Architecture.

Analysis 1: Population allele frequency divergence (S9) + positive selection (S10)
Analysis 2: AD caQTL peak accessibility classification (S3)
Analysis 3: AD coupling disruption network
Analysis 4: Combined network figure

Input:  data/raw/science.adt3130_table_s9.xlsx
        data/raw/science.adt3130_table_s10.xlsx
        data/raw/CIMA_Table_S3.csv
        results/decode_ad_summary.csv
        results/decode_ad_caqtl.csv
        results/decode_ad_coupling.csv
Output: results/decode_ad_popgen.csv
        results/decode_ad_peak_accessibility.csv
        results/decode_ad_disruption_network.csv
        figures/fig_decode_ad_popgen.pdf + .png
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path


# ======================================================================
# Analysis 1: Population genetics -- allele frequency divergence
# ======================================================================

def population_genetics(s9_path, s10_path, lead_snps_path,
                        eqtl_path='results/decode_ad_eqtl.csv',
                        caqtl_path='results/decode_ad_caqtl.csv'):
    """Compare allele frequencies across populations for AD xQTL SNPs."""
    print("=" * 60)
    print("ANALYSIS 1: Population allele frequency divergence")
    print("=" * 60)

    # Load S9 (contains eQTL/caQTL SNPs with population AF)
    s9 = pd.read_excel(s9_path)
    print(f"  S9 columns: {s9.columns.tolist()}")
    print(f"  S9 shape: {s9.shape}")

    # Load lead SNPs and AD xQTL associations
    lead = pd.read_csv(lead_snps_path)
    eqtl = pd.read_csv(eqtl_path)
    caqtl = pd.read_csv(caqtl_path)
    print(f"  AD lead SNPs: {len(lead)}")

    # S9 contains eQTL/caQTL variants, not GWAS lead SNPs
    # Match via AD eQTL/caQTL variant_ids
    ad_variants = pd.concat([eqtl[['variant_id', 'locus_id']],
                             caqtl[['variant_id', 'locus_id']]]).drop_duplicates('variant_id')
    print(f"  Unique AD xQTL variants: {ad_variants['variant_id'].nunique()}")

    merged = ad_variants.merge(s9, on='variant_id', how='inner')
    print(f"  Matched in S9: {len(merged)} variants across "
          f"{merged['locus_id'].nunique()} loci")

    # Compute AF divergence
    merged['AF_diff_EUR'] = (merged['CIMA'] - merged['EUR']).abs()
    merged['AF_diff_AFR'] = (merged['CIMA'] - merged['AFR']).abs()
    merged['divergent_EUR'] = merged['AF_diff_EUR'] > 0.1

    print(f"\n  Population AF comparison (n={len(merged)} xQTL SNPs):")
    print(f"    Mean |AF_CIMA - AF_EUR| = {merged['AF_diff_EUR'].mean():.4f}")
    print(f"    Mean |AF_CIMA - AF_AFR| = {merged['AF_diff_AFR'].mean():.4f}")
    print(f"    SNPs with |CIMA-EUR| > 0.1: {merged['divergent_EUR'].sum()} "
          f"({100*merged['divergent_EUR'].mean():.1f}%)")

    if merged['divergent_EUR'].any():
        print(f"\n  Divergent SNPs (|CIMA-EUR| > 0.1, top 10):")
        div = merged[merged['divergent_EUR']].sort_values('AF_diff_EUR', ascending=False)
        for _, row in div.head(10).iterrows():
            print(f"    {row['variant_id']:25s}  {row['locus_id']}"
                  f"  CIMA={row['CIMA']:.3f}  EUR={row['EUR']:.3f}  AFR={row['AFR']:.3f}"
                  f"  |diff|={row['AF_diff_EUR']:.3f}")

    # Load S10 -- positive selection loci
    s10 = pd.read_excel(s10_path, sheet_name=0)
    print(f"\n  S10 columns: {s10.columns.tolist()}")
    print(f"  S10 shape: {s10.shape} (Chinese selection loci)")

    # Match AD SNPs to selection loci (within 500kb)
    selection_hits = []
    for _, sel_row in s10.iterrows():
        site = str(sel_row.get('site_38', ''))
        if '_' not in site:
            continue
        parts = site.split('_')
        sel_chr = parts[0].replace('chr', '')
        try:
            sel_pos = int(parts[1])
        except (ValueError, IndexError):
            continue

        for _, snp_row in lead.iterrows():
            if str(snp_row['chromosome']) == sel_chr:
                dist = abs(snp_row['base_pair_location'] - sel_pos)
                if dist <= 500000:
                    selection_hits.append({
                        'ad_snp': snp_row['variant_id'],
                        'ad_chr': snp_row['chromosome'],
                        'ad_pos': snp_row['base_pair_location'],
                        'selection_site': site,
                        'selection_gene': sel_row.get('gene_reported', ''),
                        'selection_study': sel_row.get('study', ''),
                        'distance_bp': dist,
                    })

    print(f"\n  AD SNPs near selection loci (±500kb): {len(selection_hits)}")
    if selection_hits:
        sel_df = pd.DataFrame(selection_hits)
        for _, row in sel_df.iterrows():
            print(f"    {row['ad_snp']} (chr{row['ad_chr']}:{row['ad_pos']}) "
                  f"<-> {row['selection_site']} ({row['selection_gene']}) "
                  f"dist={row['distance_bp']/1000:.0f}kb [{row['selection_study']}]")
    else:
        print("    None found -- AD GWAS loci do not overlap Chinese selection signals")

    # Save results
    out = merged[['variant_id', 'locus_id', 'ID', 'CIMA', 'EUR', 'AFR', 'TOT',
                   'AF_diff_EUR', 'AF_diff_AFR', 'divergent_EUR']].copy()
    out.to_csv('results/decode_ad_popgen.csv', index=False)
    print(f"\nSaved to results/decode_ad_popgen.csv")

    return merged, selection_hits


# ======================================================================
# Analysis 2: AD caQTL peak accessibility classification
# ======================================================================

def peak_accessibility(s3_path, caqtl_path):
    """Classify AD caQTL peaks by breadth of chromatin accessibility."""
    print("\n" + "=" * 60)
    print("ANALYSIS 2: AD caQTL peak accessibility classification")
    print("=" * 60)

    # Load S3 (large file -- read in chunks if needed)
    print("  Loading S3 peak accessibility matrix...")
    s3 = pd.read_csv(s3_path, index_col=0)
    print(f"  S3 shape: {s3.shape}")
    print(f"  S3 cell types: {s3.shape[1]}")

    # Load AD caQTL peaks
    caqtl = pd.read_csv(caqtl_path)
    ad_peaks = caqtl['phenotype_id'].unique()
    print(f"  Unique AD caQTL peaks: {len(ad_peaks)}")

    # Match peaks to S3
    matched_peaks = [p for p in ad_peaks if p in s3.index]
    print(f"  Matched in S3: {len(matched_peaks)} / {len(ad_peaks)}")

    if len(matched_peaks) == 0:
        # Try normalizing peak format
        print("  Trying format normalization...")
        s3_peaks_set = set(s3.index)
        # caQTL format: chr1:151898585-151899086
        # S3 format: chr1:9928-10429
        matched_peaks = [p for p in ad_peaks if p in s3_peaks_set]
        print(f"  After normalization: {len(matched_peaks)}")

    if len(matched_peaks) == 0:
        print("  WARNING: No peaks matched. Checking format...")
        print(f"    AD caQTL sample: {ad_peaks[:3]}")
        print(f"    S3 index sample: {list(s3.index[:3])}")
        return None

    # Count cell types per peak
    peak_accessibility = s3.loc[matched_peaks].sum(axis=1).astype(int)

    # Classify
    results = pd.DataFrame({
        'peak': matched_peaks,
        'n_celltypes_accessible': peak_accessibility.values,
    })
    results['category'] = pd.cut(
        results['n_celltypes_accessible'],
        bins=[-1, 3, 10, 100],
        labels=['specific (1-3 CTs)', 'moderate (4-10 CTs)', 'broad (>10 CTs)']
    )

    cat_counts = results['category'].value_counts()
    print(f"\n  Peak accessibility distribution (n={len(results)}):")
    for cat in ['specific (1-3 CTs)', 'moderate (4-10 CTs)', 'broad (>10 CTs)']:
        n = cat_counts.get(cat, 0)
        pct = 100 * n / len(results)
        print(f"    {cat}: {n} ({pct:.1f}%)")

    print(f"\n  Mean cell types per peak: {results['n_celltypes_accessible'].mean():.1f}")
    print(f"  Median: {results['n_celltypes_accessible'].median():.0f}")

    # Cell-type-specific peaks in cMono_CD14
    if 'cMono_CD14' in s3.columns:
        specific = results[results['category'] == 'specific (1-3 CTs)']['peak'].values
        cmono_specific = [p for p in specific if s3.loc[p, 'cMono_CD14']]
        print(f"\n  Cell-type-specific peaks active in cMono_CD14: {len(cmono_specific)}")
        for p in cmono_specific[:10]:
            active_cts = s3.columns[s3.loc[p].astype(bool)].tolist()
            print(f"    {p} -> active in: {', '.join(active_cts)}")

    # Per-cell-type count of AD peaks
    ct_peak_counts = s3.loc[matched_peaks].sum(axis=0).sort_values(ascending=False)
    print(f"\n  Top 10 cell types by AD peak accessibility:")
    for ct, count in ct_peak_counts.head(10).items():
        print(f"    {ct}: {int(count)} peaks")

    results.to_csv('results/decode_ad_peak_accessibility.csv', index=False)
    print(f"\nSaved to results/decode_ad_peak_accessibility.csv")

    return results, ct_peak_counts


# ======================================================================
# Analysis 3: AD coupling disruption network
# ======================================================================

def coupling_disruption(coupling_path):
    """Build edge list of coupling disruption between cell types."""
    print("\n" + "=" * 60)
    print("ANALYSIS 3: AD coupling disruption network")
    print("=" * 60)

    coupling = pd.read_csv(coupling_path)
    print(f"  Coupling edges: {len(coupling)}")
    print(f"  Columns: {coupling.columns.tolist()}")

    # Filter to edges with real global_rb (>0) for meaningful comparison
    valid = coupling[coupling['global_rb'] > 0].copy()
    print(f"  Edges with global_rb > 0: {len(valid)}")

    # Compute disruption score: (global - AD) / global, higher = more disrupted
    valid['disruption_score'] = (valid['global_rb'] - valid['ad_coupling_r']) / valid['global_rb']
    valid = valid.sort_values('disruption_score', ascending=False)

    print(f"\n  Disruption score distribution:")
    print(f"    Mean: {valid['disruption_score'].mean():.3f}")
    print(f"    Median: {valid['disruption_score'].median():.3f}")
    print(f"    Positive (AD weaker): {(valid['disruption_score'] > 0).sum()} "
          f"({100*(valid['disruption_score'] > 0).mean():.1f}%)")
    print(f"    Negative (AD stronger): {(valid['disruption_score'] < 0).sum()}")

    print(f"\n  Top 10 most disrupted edges (AD coupling << global):")
    for i, (_, row) in enumerate(valid.head(10).iterrows()):
        print(f"    {i+1}. {row['cell_type_1']:30s} <-> {row['cell_type_2']:30s}  "
              f"global={row['global_rb']:.3f}  AD={row['ad_coupling_r']:.3f}  "
              f"disruption={row['disruption_score']:.3f}")

    print(f"\n  Top 5 strengthened edges (AD coupling >> global):")
    for i, (_, row) in enumerate(valid.tail(5).iloc[::-1].iterrows()):
        print(f"    {i+1}. {row['cell_type_1']:30s} <-> {row['cell_type_2']:30s}  "
              f"global={row['global_rb']:.3f}  AD={row['ad_coupling_r']:.3f}  "
              f"disruption={row['disruption_score']:.3f}")

    # Node-level disruption summary
    nodes = set(valid['cell_type_1']) | set(valid['cell_type_2'])
    node_scores = []
    for ct in nodes:
        edges = valid[(valid['cell_type_1'] == ct) | (valid['cell_type_2'] == ct)]
        node_scores.append({
            'cell_type': ct,
            'mean_disruption': edges['disruption_score'].mean(),
            'n_edges': len(edges),
            'n_disrupted': (edges['disruption_score'] > 0.2).sum(),
            'mean_ad_coupling': edges['ad_coupling_r'].mean(),
            'mean_global_rb': edges['global_rb'].mean(),
        })
    node_df = pd.DataFrame(node_scores).sort_values('mean_disruption', ascending=False)

    print(f"\n  Top 10 most disrupted cell types (mean disruption score):")
    for _, row in node_df.head(10).iterrows():
        print(f"    {row['cell_type']:35s}  mean_disruption={row['mean_disruption']:.3f}  "
              f"edges={row['n_edges']:.0f}")

    valid.to_csv('results/decode_ad_disruption_network.csv', index=False)
    print(f"\nSaved to results/decode_ad_disruption_network.csv")

    return valid, node_df


# ======================================================================
# Analysis 4: Combined figure
# ======================================================================

def make_figure(popgen_df, selection_hits, peak_results, ct_peak_counts,
                disruption_df, node_df):
    """4-panel figure: AF divergence, peak classification, disruption network, node disruption."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.linewidth': 0.8,
    })

    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 160/25.4))

    # Panel A: Allele frequency comparison (CIMA vs EUR)
    ax = axes[0, 0]
    if popgen_df is not None and len(popgen_df) > 0:
        ax.scatter(popgen_df['EUR'], popgen_df['CIMA'], s=30, c='#3498DB',
                   alpha=0.7, edgecolors='white', linewidth=0.5, zorder=3)
        # Highlight divergent
        div = popgen_df[popgen_df['divergent_EUR']]
        if len(div) > 0:
            ax.scatter(div['EUR'], div['CIMA'], s=50, c='#E74C3C',
                       alpha=0.9, edgecolors='black', linewidth=0.5, zorder=4)
            for _, row in div.head(5).iterrows():
                ax.annotate(row.get('ID', row['variant_id']), (row['EUR'], row['CIMA']),
                            fontsize=5, xytext=(3, 3), textcoords='offset points')

        ax.plot([0, 0.5], [0, 0.5], 'k--', alpha=0.3, linewidth=0.5)
        ax.fill_between([0, 0.5], [0.1, 0.6], [0, 0.5], alpha=0.05, color='red')
        ax.fill_between([0, 0.5], [0, 0.5], [-0.1, 0.4], alpha=0.05, color='red')
        ax.set_xlabel('EUR allele frequency')
        ax.set_ylabel('CIMA allele frequency')
        ax.set_xlim(-0.02, 0.52)
        ax.set_ylim(-0.02, 0.52)
    ax.set_title('AD SNP population AF divergence', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel B: Peak accessibility classification
    ax = axes[0, 1]
    if peak_results is not None:
        cats = ['specific (1-3 CTs)', 'moderate (4-10 CTs)', 'broad (>10 CTs)']
        cat_counts = peak_results['category'].value_counts()
        counts = [cat_counts.get(c, 0) for c in cats]
        colors = ['#E74C3C', '#F39C12', '#3498DB']
        bars = ax.bar(range(3), counts, color=colors, edgecolor='white', linewidth=0.5)
        ax.set_xticks(range(3))
        ax.set_xticklabels(['Specific\n(1-3 CTs)', 'Moderate\n(4-10 CTs)', 'Broad\n(>10 CTs)'],
                           fontsize=7)
        ax.set_ylabel('Number of peaks')
        for bar, count in zip(bars, counts):
            pct = 100 * count / sum(counts) if sum(counts) > 0 else 0
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                    f'{count}\n({pct:.0f}%)', ha='center', fontsize=7)
    ax.set_title('AD caQTL peak accessibility', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel C: Disruption score distribution
    ax = axes[1, 0]
    if disruption_df is not None:
        ax.hist(disruption_df['disruption_score'], bins=40, color='#8E44AD',
                alpha=0.7, edgecolor='white', linewidth=0.3)
        ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        mean_d = disruption_df['disruption_score'].mean()
        ax.axvline(mean_d, color='#E74C3C', linestyle='-', linewidth=1)
        ax.text(mean_d + 0.02, ax.get_ylim()[1]*0.9, f'mean={mean_d:.3f}',
                fontsize=7, color='#E74C3C')
        ax.set_xlabel('Disruption score (positive = AD weaker)')
        ax.set_ylabel('Edge count')
    ax.set_title('Coupling disruption distribution', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel D: Top disrupted cell types (horizontal bar)
    ax = axes[1, 1]
    if node_df is not None:
        top_nodes = node_df.head(12).iloc[::-1]
        y = np.arange(len(top_nodes))
        colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(top_nodes)))
        ax.barh(y, top_nodes['mean_disruption'], color=colors,
                edgecolor='white', linewidth=0.3)
        ax.set_yticks(y)
        labels = [ct.split('_')[0] + '_' + '_'.join(ct.split('_')[1:])
                  if len(ct.split('_')) > 1 else ct
                  for ct in top_nodes['cell_type']]
        ax.set_yticklabels(labels, fontsize=5.5)
        ax.set_xlabel('Mean disruption score')
        ax.axvline(0, color='black', linestyle='--', linewidth=0.5, alpha=0.3)
    ax.set_title('Most disrupted cell types', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel labels
    for ax, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.suptitle('DECODE-AD Population Genetics & Enhancer Architecture',
                 fontsize=11, fontweight='bold', y=0.98)

    fig_dir = Path("figures")
    fig_dir.mkdir(exist_ok=True)
    fig.savefig(fig_dir / "fig_decode_ad_popgen.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(fig_dir / "fig_decode_ad_popgen.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("\nSaved figures/fig_decode_ad_popgen.pdf + .png")


# ======================================================================
# Main
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description="DECODE-AD Population Genetics & Enhancer Architecture")
    parser.add_argument('--s9', default='data/raw/science.adt3130_table_s9.xlsx')
    parser.add_argument('--s10', default='data/raw/science.adt3130_table_s10.xlsx')
    parser.add_argument('--s3', default='data/raw/CIMA_Table_S3.csv')
    parser.add_argument('--lead', default='results/decode_ad_summary.csv')
    parser.add_argument('--caqtl', default='results/decode_ad_caqtl.csv')
    parser.add_argument('--coupling', default='results/decode_ad_coupling.csv')
    args = parser.parse_args()

    # Analysis 1
    popgen_df, selection_hits = population_genetics(
        args.s9, args.s10, args.lead,
        eqtl_path=args.caqtl.replace('caqtl', 'eqtl'),
        caqtl_path=args.caqtl)

    # Analysis 2
    result = peak_accessibility(args.s3, args.caqtl)
    if result is not None:
        peak_results, ct_peak_counts = result
    else:
        peak_results, ct_peak_counts = None, None

    # Analysis 3
    disruption_df, node_df = coupling_disruption(args.coupling)

    # Analysis 4: Figure
    make_figure(popgen_df, selection_hits, peak_results, ct_peak_counts,
                disruption_df, node_df)

    print("\n" + "=" * 60)
    print("DECODE-AD POPGEN & ENHANCER -- KEY FINDINGS")
    print("=" * 60)
    if popgen_df is not None:
        n_div = popgen_df['divergent_EUR'].sum()
        print(f"  1. Population genetics: {n_div} AD SNPs with |CIMA-EUR| > 0.1")
        print(f"     {len(selection_hits)} AD loci overlap Chinese selection signals")
    if peak_results is not None:
        cats = peak_results['category'].value_counts()
        print(f"  2. Peak accessibility: {cats.get('specific (1-3 CTs)', 0)} specific, "
              f"{cats.get('moderate (4-10 CTs)', 0)} moderate, "
              f"{cats.get('broad (>10 CTs)', 0)} broad")
    if disruption_df is not None:
        print(f"  3. Coupling disruption: {(disruption_df['disruption_score'] > 0).sum()} / "
              f"{len(disruption_df)} edges weakened in AD context")
    print("=" * 60)


if __name__ == "__main__":
    main()
