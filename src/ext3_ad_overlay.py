#!/usr/bin/env python3
"""
Extension #3 remaining:
C. AD GWAS overlay on peak-gene pairs from S4
D. Alluvial/chord visualization of AD peak-gene pairs
"""

import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def ad_peak_gene_overlay(s4_path, s3_path, summary_path, topple_path):
    """Check which S4 Region-Gene pairs overlap AD GWAS loci."""
    print("=" * 60)
    print("EXT3 STEP 5: AD GWAS overlay on peak-gene pairs")
    print("=" * 60)

    # Load AD loci
    summary = pd.read_csv(summary_path)
    print(f"  AD loci: {len(summary)}")

    # Build locus windows
    loci = []
    for _, row in summary.iterrows():
        loci.append({
            'locus_id': row['locus_id'],
            'chr': int(row['chromosome']),
            'start': int(row['base_pair_location']) - 500000,
            'end': int(row['base_pair_location']) + 500000,
        })

    # Load S4 GRN
    print("  Loading S4 GRN...")
    s4 = pd.read_csv(s4_path)
    print(f"  S4 shape: {s4.shape}")

    # Parse Region coordinates
    s4['region_chr'] = s4['Region'].str.extract(r'chr(\d+):')[0].astype(float)
    region_coords = s4['Region'].str.extract(r'chr\d+:(\d+)-(\d+)')
    s4['region_start'] = region_coords[0].astype(float)
    s4['region_end'] = region_coords[1].astype(float)

    # Get unique Region-Gene pairs
    rg_pairs = s4[['Region', 'Gene', 'TF', 'region_chr', 'region_start', 'region_end']].drop_duplicates(
        subset=['Region', 'Gene'])
    print(f"  Unique Region-Gene pairs: {len(rg_pairs)}")

    # Match to AD loci
    hits = []
    for locus in loci:
        mask = (
            (rg_pairs['region_chr'] == locus['chr']) &
            (rg_pairs['region_start'] >= locus['start']) &
            (rg_pairs['region_end'] <= locus['end'])
        )
        locus_hits = rg_pairs[mask].copy()
        if len(locus_hits) > 0:
            locus_hits['locus_id'] = locus['locus_id']
            hits.append(locus_hits)

    if hits:
        hit_df = pd.concat(hits, ignore_index=True)
    else:
        hit_df = pd.DataFrame()

    print(f"\n  AD-overlapping Region-Gene pairs: {len(hit_df)}")
    if len(hit_df) > 0:
        print(f"  Unique regions: {hit_df['Region'].nunique()}")
        print(f"  Unique genes: {hit_df['Gene'].nunique()}")
        print(f"  Unique TFs: {hit_df['TF'].nunique()}")
        print(f"  AD loci with hits: {hit_df['locus_id'].nunique()}")

    # Load S3 for cell-type activity
    print("\n  Loading S3 peak accessibility...")
    s3 = pd.read_csv(s3_path, index_col=0)

    # Check which AD peaks are in S3
    if len(hit_df) > 0:
        ad_regions = hit_df['Region'].unique()
        matched_in_s3 = [r for r in ad_regions if r in s3.index]
        print(f"  AD regions in S3: {len(matched_in_s3)} / {len(ad_regions)}")

        if matched_in_s3:
            s3_subset = s3.loc[matched_in_s3]
            ct_activity = s3_subset.sum(axis=0).sort_values(ascending=False)
            print(f"\n  Top 10 cell types by AD peak-gene activity:")
            for ct, n in ct_activity.head(10).items():
                print(f"    {ct}: {int(n)} peaks")

            # Annotate each hit with active cell types count
            hit_df['n_active_celltypes'] = hit_df['Region'].apply(
                lambda r: int(s3.loc[r].sum()) if r in s3.index else 0
            )

        # Load TOPPLE for regulon membership
        topple = pd.read_csv(topple_path)
        tf_class = topple.set_index(
            topple['regulon'].str.replace(r'_\+$', '', regex=True)
        )['stability_class'].to_dict()

        hit_df['tf_stability'] = hit_df['TF'].map(tf_class).fillna('not_in_TOPPLE')

        # Summary by TF
        tf_summary = hit_df.groupby('TF').agg(
            n_regions=('Region', 'nunique'),
            n_genes=('Gene', 'nunique'),
            n_loci=('locus_id', 'nunique'),
            stability=('tf_stability', 'first'),
        ).sort_values('n_regions', ascending=False)

        print(f"\n  Top 10 TFs (regulons) by AD peak-gene coverage:")
        for tf, row in tf_summary.head(10).iterrows():
            print(f"    {tf:15s}: {row['n_regions']:3.0f} regions, "
                  f"{row['n_genes']:.0f} genes, {row['n_loci']:.0f} loci  [{row['stability']}]")

        # By locus
        locus_summary = hit_df.groupby('locus_id').agg(
            n_pairs=('Region', 'size'),
            n_regions=('Region', 'nunique'),
            n_genes=('Gene', 'nunique'),
            n_tfs=('TF', 'nunique'),
        ).sort_values('n_pairs', ascending=False)
        print(f"\n  AD loci by peak-gene pair count:")
        for loc, row in locus_summary.head(10).iterrows():
            print(f"    {loc}: {row['n_pairs']:.0f} pairs, "
                  f"{row['n_regions']:.0f} regions, {row['n_genes']:.0f} genes")

    hit_df.to_csv('results/ext3_ad_peak_gene_overlay.csv', index=False)
    print(f"\nSaved to results/ext3_ad_peak_gene_overlay.csv")
    return hit_df


def alluvial_figure(overlay_path):
    """Create alluvial-style plot: AD locus -> peak -> gene -> TF/regulon."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.sankey import Sankey
    import matplotlib.patches as mpatches

    print("\n" + "=" * 60)
    print("EXT3: AD overlay figure")
    print("=" * 60)

    overlay = pd.read_csv(overlay_path)
    if len(overlay) == 0:
        print("  No overlay data, skipping figure.")
        return

    # Create a grouped bar chart: top loci by peak-gene pairs, colored by TF stability
    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 160/25.4))

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})

    # Panel A: Locus x peak-gene pair count
    ax = axes[0, 0]
    locus_counts = overlay.groupby('locus_id').size().sort_values(ascending=False).head(15)
    y = np.arange(len(locus_counts))
    ax.barh(y, locus_counts.values, color='#3498DB', edgecolor='white', linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(locus_counts.index, fontsize=6)
    ax.set_xlabel('Peak-gene pairs')
    ax.set_title('AD loci: regulatory pairs', fontsize=9, fontweight='bold')
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel B: TF regulon distribution (stabilizer vs destabilizer)
    ax = axes[0, 1]
    if 'tf_stability' in overlay.columns:
        stab_counts = overlay.groupby('tf_stability').size()
        cats = ['stabilizer', 'intermediate', 'destabilizer', 'not_in_TOPPLE']
        colors_map = {'stabilizer': '#2ECC71', 'intermediate': '#F39C12',
                      'destabilizer': '#E74C3C', 'not_in_TOPPLE': '#BDC3C7'}
        vals = [stab_counts.get(c, 0) for c in cats]
        bars = ax.bar(range(len(cats)), vals,
                       color=[colors_map[c] for c in cats],
                       edgecolor='white', linewidth=0.5)
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(['Stab', 'Inter', 'Destab', 'N/A'], fontsize=7)
        ax.set_ylabel('Peak-gene pairs')
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                        str(v), ha='center', fontsize=7)
    ax.set_title('Regulatory stability class', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel C: Top genes by peak count
    ax = axes[1, 0]
    gene_counts = overlay.groupby('Gene').agg(
        n_peaks=('Region', 'nunique'),
        n_loci=('locus_id', 'nunique'),
    ).sort_values('n_peaks', ascending=False).head(15)
    y = np.arange(len(gene_counts))
    ax.barh(y, gene_counts['n_peaks'], color='#E74C3C', edgecolor='white', linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(gene_counts.index, fontsize=6)
    ax.set_xlabel('Enhancer peaks')
    ax.set_title('Top AD target genes', fontsize=9, fontweight='bold')
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel D: Cell-type accessibility of AD peaks
    ax = axes[1, 1]
    if 'n_active_celltypes' in overlay.columns:
        ax.hist(overlay['n_active_celltypes'], bins=30, color='#8E44AD',
                alpha=0.7, edgecolor='white', linewidth=0.3)
        mean_ct = overlay['n_active_celltypes'].mean()
        ax.axvline(mean_ct, color='#E74C3C', linestyle='--', linewidth=1)
        ax.text(mean_ct + 1, ax.get_ylim()[1]*0.9 if ax.get_ylim()[1] > 0 else 10,
                f'mean={mean_ct:.1f}', fontsize=7, color='#E74C3C')
        ax.set_xlabel('Cell types accessible')
        ax.set_ylabel('Peak-gene pairs')
    ax.set_title('Peak accessibility breadth', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for ax, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')

    fig.suptitle('Extension #3: AD GWAS Regulatory Overlay', fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_ext3_ad_overlay.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_ext3_ad_overlay.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_ext3_ad_overlay.pdf + .png")


if __name__ == "__main__":
    hit_df = ad_peak_gene_overlay(
        'data/raw/CIMA_Table_S4.csv',
        'data/raw/CIMA_Table_S3.csv',
        'results/decode_ad_summary.csv',
        'results/topple.csv',
    )
    alluvial_figure('results/ext3_ad_peak_gene_overlay.csv')
