#!/usr/bin/env python3
"""
DECODE-AD Final Analyses -- six remaining analyses.

Analysis 1: Cell-type PRS decomposition
Analysis 2: Regulatory cascade analysis (S4 GRN + S11 trans-eQTL)
Analysis 3: Network visualization (networkx)
Analysis 4: Drug target deep analysis (eQTL + TOPPLE)
Analysis 5: S12 interaction QTLs for AD
Analysis 6: East Asian (Tanaka 2021) sensitivity analysis

Output: results/decode_ad_celltype_prs.csv
        results/decode_ad_cascades.csv
        results/decode_ad_drug_targets_deep.csv
        results/decode_ad_interaction_qtl.csv
        results/decode_ad_tanaka_comparison.csv
        figures/fig_decode_ad_network.pdf
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


# ======================================================================
# Analysis 1: Cell-type PRS decomposition
# ======================================================================

def celltype_prs(s6_path, eqtl_path, summary_path):
    """Decompose AD polygenic risk by cell type using eQTL slopes."""
    print("=" * 60)
    print("ANALYSIS 1: Cell-type PRS decomposition")
    print("=" * 60)

    # Load AD eQTL results (already matched to AD loci)
    ad_eqtl = pd.read_csv(eqtl_path)
    summary = pd.read_csv(summary_path)
    print(f"  AD eQTL associations: {len(ad_eqtl)}")
    print(f"  Unique AD eGenes: {ad_eqtl['phenotype_id'].nunique()}")

    # Get GWAS beta per locus
    locus_beta = summary.set_index('locus_id')['beta'].to_dict()

    # For each association: PRS contribution = |GWAS_beta| * |eQTL_slope|
    ad_eqtl['gwas_beta'] = ad_eqtl['locus_id'].map(locus_beta)
    ad_eqtl['prs_weight'] = ad_eqtl['gwas_beta'].abs() * ad_eqtl['slope'].abs()

    # Sum per cell type
    ct_prs = ad_eqtl.groupby('celltype').agg(
        total_prs_weight=('prs_weight', 'sum'),
        n_genes=('phenotype_id', 'nunique'),
        n_loci=('locus_id', 'nunique'),
        mean_abs_slope=('slope', lambda x: x.abs().mean()),
    ).sort_values('total_prs_weight', ascending=False)

    print(f"\n  Top 15 cell types by PRS contribution:")
    for i, (ct, row) in enumerate(ct_prs.head(15).iterrows()):
        print(f"    {i+1:2d}. {ct:35s}  PRS={row['total_prs_weight']:.4f}  "
              f"genes={row['n_genes']:.0f}  loci={row['n_loci']:.0f}")

    # Normalize to percentage
    ct_prs['prs_pct'] = 100 * ct_prs['total_prs_weight'] / ct_prs['total_prs_weight'].sum()

    # Lineage grouping
    lineage_map = {}
    for ct in ct_prs.index:
        if any(x in ct for x in ['Mono', 'DC', 'pDC']):
            lineage_map[ct] = 'Myeloid'
        elif any(x in ct for x in ['CD4', 'CD8', 'MAIT', 'gdT', 'NKT', 'Treg']):
            lineage_map[ct] = 'T cell'
        elif any(x in ct for x in ['NK', 'ILC']):
            lineage_map[ct] = 'NK/ILC'
        elif any(x in ct for x in ['Bn', 'Bm', 'Plasma', 'Switched', 'Unswitched', 'Transitional_B']):
            lineage_map[ct] = 'B cell'
        else:
            lineage_map[ct] = 'Other'

    ct_prs['lineage'] = ct_prs.index.map(lineage_map)
    lineage_prs = ct_prs.groupby('lineage')['total_prs_weight'].sum()
    lineage_prs = 100 * lineage_prs / lineage_prs.sum()
    print(f"\n  PRS by lineage:")
    for lin, pct in lineage_prs.sort_values(ascending=False).items():
        print(f"    {lin:15s}: {pct:.1f}%")

    ct_prs.to_csv('results/decode_ad_celltype_prs.csv')
    print(f"\nSaved to results/decode_ad_celltype_prs.csv")
    return ct_prs


# ======================================================================
# Analysis 2: Regulatory cascade analysis
# ======================================================================

def regulatory_cascades(s4_path, s11_path, eqtl_path):
    """Find regulatory cascades: AD SNP -> cis-eQTL on TF -> TF regulates targets."""
    print("\n" + "=" * 60)
    print("ANALYSIS 2: Regulatory cascade analysis")
    print("=" * 60)

    # Load data
    print("  Loading S4 GRN...")
    s4 = pd.read_csv(s4_path)
    print(f"  S4 columns: {s4.columns.tolist()}")
    print(f"  S4 shape: {s4.shape}")

    print("  Loading S11 trans-eQTLs...")
    s11 = pd.read_excel(s11_path)
    print(f"  S11 columns: {s11.columns.tolist()}")
    print(f"  S11 shape: {s11.shape}")

    ad_eqtl = pd.read_csv(eqtl_path)
    ad_egenes = set(ad_eqtl['phenotype_id'].unique())
    print(f"  AD eGenes: {len(ad_egenes)}")

    # Get unique TFs in S4
    s4_tfs = set(s4['TF'].unique())
    print(f"  S4 TFs: {len(s4_tfs)}")

    # Find AD eGenes that are TFs in S4
    ad_tf_genes = ad_egenes & s4_tfs
    print(f"  AD eGenes that are TFs in S4: {len(ad_tf_genes)}")
    if ad_tf_genes:
        print(f"    {sorted(ad_tf_genes)}")

    # Build cascade chains
    cascades = []

    # Chain type 1: AD SNP -> cis-eQTL on TF -> TF regulates targets via S4
    for tf in ad_tf_genes:
        # Get eQTL info for this TF
        tf_eqtls = ad_eqtl[ad_eqtl['phenotype_id'] == tf]

        # Get S4 targets of this TF
        tf_targets = s4[s4['TF'] == tf][['Gene', 'Region']].drop_duplicates()

        for _, eqtl_row in tf_eqtls.iterrows():
            for _, target_row in tf_targets.iterrows():
                cascades.append({
                    'chain_type': 'cis-eQTL->TF->target',
                    'ad_locus': eqtl_row['locus_id'],
                    'ad_variant': eqtl_row['variant_id'],
                    'celltype': eqtl_row['celltype'],
                    'cis_gene_TF': tf,
                    'eqtl_slope': eqtl_row['slope'],
                    'eqtl_pval': eqtl_row['pval_nominal'],
                    'target_gene': target_row['Gene'],
                    'target_region': target_row['Region'],
                    'trans_eqtl_confirmed': False,
                    'trans_pval': np.nan,
                })

    print(f"\n  Type 1 cascades (cis-eQTL->TF->S4 target): {len(cascades)}")

    # Chain type 2: Check if any trans-eQTLs confirm the cascade
    # S11 has cis_eGene and trans_eGene columns
    if 'cis_eGene' in s11.columns and 'trans_eGene' in s11.columns:
        for _, trans_row in s11.iterrows():
            cis_gene = trans_row.get('cis_eGene', '')
            trans_gene = trans_row.get('trans_eGene', '')

            if cis_gene in ad_tf_genes:
                # Mark matching cascades as trans-confirmed
                mask = [(c['cis_gene_TF'] == cis_gene and c['target_gene'] == trans_gene)
                        for c in cascades]
                n_confirmed = sum(mask)
                for i, m in enumerate(mask):
                    if m:
                        cascades[i]['trans_eqtl_confirmed'] = True
                        cascades[i]['trans_pval'] = trans_row.get('pval', np.nan)

                if n_confirmed == 0 and cis_gene in ad_tf_genes:
                    # Add as new cascade even if target not in S4
                    cascades.append({
                        'chain_type': 'cis-eQTL->TF->trans-eQTL',
                        'ad_locus': 'multiple',
                        'ad_variant': trans_row.get('variant_id', ''),
                        'celltype': trans_row.get('celltype', ''),
                        'cis_gene_TF': cis_gene,
                        'eqtl_slope': np.nan,
                        'eqtl_pval': np.nan,
                        'target_gene': trans_gene,
                        'target_region': '',
                        'trans_eqtl_confirmed': True,
                        'trans_pval': trans_row.get('pval', np.nan),
                    })

    # Also check trans-eQTLs where the trans_eGene is an AD eGene
    trans_ad_hits = s11[s11['trans_eGene'].isin(ad_egenes)] if 'trans_eGene' in s11.columns else pd.DataFrame()
    print(f"  Trans-eQTLs targeting AD eGenes: {len(trans_ad_hits)}")
    if len(trans_ad_hits) > 0:
        for _, row in trans_ad_hits.iterrows():
            print(f"    {row.get('cis_eGene','')} -> {row['trans_eGene']} "
                  f"in {row.get('celltype','')} (P={row.get('pval','?')})")

    cascade_df = pd.DataFrame(cascades)
    n_confirmed = cascade_df['trans_eqtl_confirmed'].sum() if len(cascade_df) > 0 else 0
    print(f"\n  Total cascades: {len(cascade_df)}")
    print(f"  Trans-eQTL confirmed: {n_confirmed}")

    if len(cascade_df) > 0:
        # Summary per TF
        tf_summary = cascade_df.groupby('cis_gene_TF').agg(
            n_targets=('target_gene', 'nunique'),
            n_celltypes=('celltype', 'nunique'),
            n_confirmed=('trans_eqtl_confirmed', 'sum'),
        ).sort_values('n_targets', ascending=False)
        print(f"\n  Cascade TF summary:")
        for tf, row in tf_summary.iterrows():
            print(f"    {tf:15s}: {row['n_targets']:5.0f} targets, "
                  f"{row['n_celltypes']:.0f} CTs, {row['n_confirmed']:.0f} trans-confirmed")

        cascade_df.to_csv('results/decode_ad_cascades.csv', index=False)
        print(f"\nSaved to results/decode_ad_cascades.csv")

    return cascade_df


# ======================================================================
# Analysis 3: Network visualization
# ======================================================================

def network_figure(disruption_path):
    """Create network figure of coupling disruption."""
    print("\n" + "=" * 60)
    print("ANALYSIS 3: Network visualization")
    print("=" * 60)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    try:
        import networkx as nx
    except ImportError:
        print("  WARNING: networkx not installed. Installing...")
        import subprocess
        subprocess.check_call(['pip', 'install', 'networkx'])
        import networkx as nx

    disruption = pd.read_csv(disruption_path)
    print(f"  Disruption edges: {len(disruption)}")

    # Node-level scores
    nodes_all = set(disruption['cell_type_1']) | set(disruption['cell_type_2'])
    node_scores = {}
    for ct in nodes_all:
        edges = disruption[(disruption['cell_type_1'] == ct) |
                           (disruption['cell_type_2'] == ct)]
        node_scores[ct] = edges['disruption_score'].mean()

    # Top 20 most affected nodes (by absolute disruption)
    sorted_nodes = sorted(node_scores.items(), key=lambda x: abs(x[1]), reverse=True)
    top_nodes = [n for n, _ in sorted_nodes[:20]]

    # Build network
    G = nx.Graph()
    for ct in top_nodes:
        G.add_node(ct, disruption=node_scores[ct])

    # Add edges between top nodes
    for _, row in disruption.iterrows():
        ct1, ct2 = row['cell_type_1'], row['cell_type_2']
        if ct1 in top_nodes and ct2 in top_nodes:
            G.add_edge(ct1, ct2,
                        disruption=row['disruption_score'],
                        global_rb=row['global_rb'],
                        ad_coupling=row['ad_coupling_r'])

    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Layout
    pos = nx.spring_layout(G, k=2.5, iterations=100, seed=42)

    fig, ax = plt.subplots(1, 1, figsize=(180/25.4, 180/25.4))

    # Draw edges
    for u, v, data in G.edges(data=True):
        d = data['disruption']
        if d > 0.2:
            color = '#E74C3C'  # red = disrupted
            alpha = min(0.8, 0.3 + abs(d) * 0.3)
        elif d < -0.2:
            color = '#3498DB'  # blue = strengthened
            alpha = min(0.8, 0.3 + abs(d) * 0.3)
        else:
            color = '#BDC3C7'  # gray = unchanged
            alpha = 0.15
        width = max(0.3, min(3.0, abs(d) * 2))
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width,
                                edge_color=color, alpha=alpha, ax=ax)

    # Draw nodes
    node_sizes = [max(100, abs(node_scores[n]) * 800) for n in G.nodes()]
    node_colors = []
    for n in G.nodes():
        if n == 'CD4_Tr1_IL10':
            node_colors.append('#E74C3C')
        elif n == 'MK_GP9':
            node_colors.append('#3498DB')
        elif node_scores[n] > 0.2:
            node_colors.append('#F1948A')
        elif node_scores[n] < -0.2:
            node_colors.append('#85C1E9')
        else:
            node_colors.append('#D5D8DC')

    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                            edgecolors='black', linewidths=0.5, ax=ax)

    # Labels
    labels = {n: n.split('_')[0] + '_' + '_'.join(n.split('_')[1:2])
              for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=5, ax=ax)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#E74C3C', label='CD4_Tr1_IL10 (most disrupted)'),
        mpatches.Patch(facecolor='#3498DB', label='MK_GP9 (strengthened)'),
        plt.Line2D([0], [0], color='#E74C3C', linewidth=2, label='Disrupted edge'),
        plt.Line2D([0], [0], color='#3498DB', linewidth=2, label='Strengthened edge'),
        plt.Line2D([0], [0], color='#BDC3C7', linewidth=1, label='Unchanged'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=6, frameon=True,
              fancybox=True, framealpha=0.9)

    ax.set_title('AD Coupling Disruption Network', fontsize=11, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()

    fig_dir = Path("figures")
    fig_dir.mkdir(exist_ok=True)
    fig.savefig(fig_dir / "fig_decode_ad_network.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(fig_dir / "fig_decode_ad_network.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_decode_ad_network.pdf + .png")


# ======================================================================
# Analysis 4: Drug target deep analysis
# ======================================================================

def drug_target_deep(s6_path, topple_path):
    """Deep analysis of drug target eQTLs + TOPPLE stability."""
    print("\n" + "=" * 60)
    print("ANALYSIS 4: Drug target deep analysis")
    print("=" * 60)

    # Drug targets
    targets = ['IL4R', 'IL13', 'IL13RA1', 'IL13RA2', 'IL31RA', 'JAK1', 'JAK2',
               'JAK3', 'TYK2', 'IL5', 'IL5RA', 'TSLP', 'IL33', 'IL1RL1',
               'PDE4A', 'PDE4B', 'PDE4D']

    # Load S6 (large file)
    print("  Loading S6 cis-eQTLs...")
    s6 = pd.read_csv(s6_path)
    s6_eqtl = s6[s6['analysis'] == 'cis-eQTL']
    print(f"  S6 eQTL rows: {len(s6_eqtl)}")

    # Load TOPPLE
    topple = pd.read_csv(topple_path)
    print(f"  TOPPLE regulons: {len(topple)}")
    # Map TF name to RI (strip _+ suffix)
    topple['tf_name'] = topple['regulon'].str.replace(r'_\+$', '', regex=True)
    tf_ri = topple.set_index('tf_name')['mean_RI'].to_dict()
    tf_class = topple.set_index('tf_name')['stability_class'].to_dict()

    results = []
    for gene in targets:
        gene_eqtls = s6_eqtl[s6_eqtl['phenotype_id'] == gene]
        if len(gene_eqtls) == 0:
            results.append({
                'gene': gene, 'found_in_s6': False, 'n_celltypes': 0,
                'mean_abs_slope': np.nan, 'celltypes': '',
                'topple_ri': tf_ri.get(gene, np.nan),
                'stability_class': tf_class.get(gene, 'not_in_TOPPLE'),
                'druggability_score': 0,
            })
            continue

        celltypes = gene_eqtls['celltype'].unique()
        mean_slope = gene_eqtls['slope'].abs().mean()

        # Get mean TOPPLE RI across cell types that have eQTL for this gene
        # (RI is per regulon, not per cell type, so use the gene's own RI if it's a TF)
        gene_ri = tf_ri.get(gene, np.nan)

        # Druggability = n_celltypes * mean_RI (if available)
        drug_score = len(celltypes) * (gene_ri if not np.isnan(gene_ri) else mean_slope)

        results.append({
            'gene': gene,
            'found_in_s6': True,
            'n_celltypes': len(celltypes),
            'mean_abs_slope': mean_slope,
            'celltypes': ';'.join(sorted(celltypes)),
            'topple_ri': gene_ri,
            'stability_class': tf_class.get(gene, 'not_in_TOPPLE'),
            'druggability_score': drug_score,
        })

    result_df = pd.DataFrame(results).sort_values('druggability_score', ascending=False)

    print(f"\n  Drug target eQTL summary:")
    found = result_df[result_df['found_in_s6']]
    not_found = result_df[~result_df['found_in_s6']]
    print(f"    Found in S6 eQTLs: {len(found)} / {len(targets)}")
    print(f"    Not found: {', '.join(not_found['gene'].values)}")

    print(f"\n  Drug targets with eQTLs (ranked by druggability):")
    for _, row in found.iterrows():
        print(f"    {row['gene']:12s}  CTs={row['n_celltypes']:3.0f}  "
              f"|slope|={row['mean_abs_slope']:.3f}  "
              f"RI={row['topple_ri']:.5f}" if not np.isnan(row['topple_ri'])
              else f"    {row['gene']:12s}  CTs={row['n_celltypes']:3.0f}  "
              f"|slope|={row['mean_abs_slope']:.3f}  RI=N/A  "
              f"drug_score={row['druggability_score']:.3f}")

    result_df.to_csv('results/decode_ad_drug_targets_deep.csv', index=False)
    print(f"\nSaved to results/decode_ad_drug_targets_deep.csv")
    return result_df


# ======================================================================
# Analysis 5: S12 interaction QTLs
# ======================================================================

def interaction_qtls(s12_path, summary_path, eqtl_path, caqtl_path):
    """Check if AD-associated variants are context-dependent (interaction QTLs)."""
    print("\n" + "=" * 60)
    print("ANALYSIS 5: S12 interaction QTLs for AD")
    print("=" * 60)

    s12 = pd.read_csv(s12_path)
    print(f"  S12 columns: {s12.columns.tolist()}")
    print(f"  S12 shape: {s12.shape}")

    summary = pd.read_csv(summary_path)
    ad_eqtl = pd.read_csv(eqtl_path)
    ad_caqtl = pd.read_csv(caqtl_path)

    # S12 variant_id is in chr_pos format; AD eQTL/caQTL variant_ids also chr_pos
    ad_xqtl_variants = set(ad_eqtl['variant_id'].unique()) | set(ad_caqtl['variant_id'].unique())
    print(f"  Unique AD xQTL variants: {len(ad_xqtl_variants)}")
    print(f"  S12 unique variants: {s12['variant_id'].nunique()}")

    # Direct match
    direct = s12[s12['variant_id'].isin(ad_xqtl_variants)].copy()
    print(f"  Direct variant overlap: {len(direct)}")

    if len(direct) > 0:
        # Reshape direct matches to standard format
        overlap_rows = []
        for _, hit in direct.iterrows():
            overlap_rows.append({
                'ad_locus': 'direct_match',
                'ad_snp': hit['variant_id'],
                'ad_chr': '',
                'ad_pos': '',
                'iqtl_variant': hit['variant_id'],
                'iqtl_gene': hit['phenotype_id'],
                'iqtl_celltype': hit['celltype'],
                'iqtl_slope': hit['slope'],
                'iqtl_pval': hit['pval_nominal'],
                'dynamic_LRT_p': hit.get('dynamic_LRT_p', np.nan),
                'distance': 0,
            })
        overlap = pd.DataFrame(overlap_rows)
    else:
        overlap = pd.DataFrame()

    # Also try proximity matching
    if True:
        print("  Trying proximity matching (10kb window around AD lead SNPs)...")
        # Parse S12 variant positions
        s12['chrom'] = s12['variant_id'].str.extract(r'chr(\d+)_')[0].astype(float)
        s12['pos'] = s12['variant_id'].str.extract(r'chr\d+_(\d+)')[0].astype(float)

        prox_hits = []
        for _, snp in summary.iterrows():
            s12_chr = s12[s12['chrom'] == snp['chromosome']]
            nearby = s12_chr[
                (s12_chr['pos'] >= snp['base_pair_location'] - 500000) &
                (s12_chr['pos'] <= snp['base_pair_location'] + 500000)
            ]
            if len(nearby) > 0:
                for _, hit in nearby.iterrows():
                    prox_hits.append({
                        'ad_locus': snp['locus_id'],
                        'ad_snp': snp['variant_id'],
                        'ad_chr': snp['chromosome'],
                        'ad_pos': snp['base_pair_location'],
                        'iqtl_variant': hit['variant_id'],
                        'iqtl_gene': hit['phenotype_id'],
                        'iqtl_celltype': hit['celltype'],
                        'iqtl_slope': hit['slope'],
                        'iqtl_pval': hit['pval_nominal'],
                        'dynamic_LRT_p': hit.get('dynamic_LRT_p', np.nan),
                        'distance': abs(hit['pos'] - snp['base_pair_location']),
                    })

        prox_df = pd.DataFrame(prox_hits)
        print(f"  Proximity matches (500kb): {len(prox_df)}")
        if len(prox_df) > 0:
            overlap = pd.concat([overlap, prox_df], ignore_index=True).drop_duplicates(
                subset=['iqtl_variant', 'iqtl_gene', 'iqtl_celltype'])

    if len(overlap) > 0:
        print(f"\n  Interaction QTL hits:")
        print(f"    Unique AD loci: {overlap['ad_locus'].nunique()}")
        print(f"    Unique genes: {overlap['iqtl_gene'].nunique()}")
        print(f"    Unique cell types: {overlap['iqtl_celltype'].nunique()}")

        if 'iqtl_celltype' in overlap.columns:
            ct_counts = overlap['iqtl_celltype'].value_counts().head(10)
            print(f"\n  Top cell types with AD interaction QTLs:")
            for ct, n in ct_counts.items():
                print(f"    {ct}: {n}")

        if 'iqtl_gene' in overlap.columns:
            gene_counts = overlap['iqtl_gene'].value_counts().head(10)
            print(f"\n  Top genes with AD interaction QTLs:")
            for gene, n in gene_counts.items():
                print(f"    {gene}: {n}")

        overlap.to_csv('results/decode_ad_interaction_qtl.csv', index=False)
        print(f"\nSaved to results/decode_ad_interaction_qtl.csv")
    else:
        print("  No interaction QTL overlaps found.")
        pd.DataFrame(columns=['note']).to_csv('results/decode_ad_interaction_qtl.csv', index=False)

    return overlap


# ======================================================================
# Analysis 6: East Asian (Tanaka 2021) sensitivity analysis
# ======================================================================

def tanaka_comparison(s6_path, s9_path, summary_path):
    """Compare European vs Japanese AD GWAS overlap with CIMA xQTLs."""
    print("\n" + "=" * 60)
    print("ANALYSIS 6: East Asian (Tanaka 2021) sensitivity analysis")
    print("=" * 60)

    # Tanaka 2021 Japanese AD GWAS -- 17 lead SNPs
    # Known loci with approximate hg38 positions from literature
    tanaka_snps = [
        {'rsid': 'rs4713555', 'chr': 2, 'pos': 102300000},   # IL18R1 region
        {'rsid': 'rs2228145', 'chr': 1, 'pos': 154426264},   # IL6R
        {'rsid': 'rs12188917', 'chr': 5, 'pos': 132660151},  # IL13/IRF1
        {'rsid': 'rs2897442', 'chr': 11, 'pos': 76297081},   # LRRC32
        {'rsid': 'rs7127307', 'chr': 11, 'pos': 76369609},   # EMSY/LRRC32
        {'rsid': 'rs2041733', 'chr': 12, 'pos': 56346850},   # STAT6
        {'rsid': 'rs11923593', 'chr': 3, 'pos': 188390580},  # LPP
        {'rsid': 'rs12951971', 'chr': 16, 'pos': 27324424},  # IL4R region
        {'rsid': 'rs10214237', 'chr': 10, 'pos': 6076785},   # IL15RA
        {'rsid': 'rs12295535', 'chr': 11, 'pos': 36381082},  # FLG/TSHZ2
        {'rsid': 'rs1250546', 'chr': 2, 'pos': 60966955},    # BCL11A
        {'rsid': 'rs2592555', 'chr': 5, 'pos': 110410850},   # CAMK4
        {'rsid': 'rs6419573', 'chr': 20, 'pos': 48563555},   # PTGIS
        {'rsid': 'rs2143950', 'chr': 7, 'pos': 37375270},    # ELMO1
        {'rsid': 'rs7512552', 'chr': 1, 'pos': 67846927},    # C1orf141
        {'rsid': 'rs12634229', 'chr': 3, 'pos': 71258849},   # FOXP1
        {'rsid': 'rs4312054', 'chr': 8, 'pos': 130695218},   # DOCK8
    ]
    tanaka_df = pd.DataFrame(tanaka_snps)
    print(f"  Tanaka 2021 lead SNPs: {len(tanaka_df)}")

    # Try to get exact positions from S9
    s9 = pd.read_excel(s9_path)
    for i, row in tanaka_df.iterrows():
        s9_match = s9[s9['ID'] == row['rsid']]
        if len(s9_match) > 0:
            # Parse chr_pos from S9 variant_id
            vid = s9_match.iloc[0]['variant_id']
            parts = vid.split('_')
            tanaka_df.loc[i, 'chr'] = int(parts[0].replace('chr', ''))
            tanaka_df.loc[i, 'pos'] = int(parts[1])
            tanaka_df.loc[i, 's9_matched'] = True
        else:
            tanaka_df.loc[i, 's9_matched'] = False
    print(f"  Matched in S9: {tanaka_df['s9_matched'].sum()}")

    # Load S6
    print("  Loading S6...")
    s6 = pd.read_csv(s6_path)
    s6['chrom'] = s6['variant_id'].str.extract(r'chr(\d+)_')[0].astype(float)
    s6['pos_s6'] = s6['variant_id'].str.extract(r'chr\d+_(\d+)')[0].astype(float)

    s6_eqtl = s6[s6['analysis'] == 'cis-eQTL']
    s6_caqtl = s6[s6['analysis'] == 'cis-caQTL']

    # Match Tanaka SNPs to S6 (500kb window)
    tanaka_eqtl_hits = []
    tanaka_caqtl_hits = []

    for _, snp in tanaka_df.iterrows():
        chr_val = int(snp['chr'])
        pos_val = int(snp['pos'])

        # eQTL
        eqtl_chr = s6_eqtl[s6_eqtl['chrom'] == chr_val]
        eqtl_nearby = eqtl_chr[
            (eqtl_chr['pos_s6'] >= pos_val - 500000) &
            (eqtl_chr['pos_s6'] <= pos_val + 500000)
        ]
        for _, hit in eqtl_nearby.iterrows():
            tanaka_eqtl_hits.append({
                'tanaka_rsid': snp['rsid'],
                'tanaka_chr': chr_val,
                'tanaka_pos': pos_val,
                'phenotype_id': hit['phenotype_id'],
                'variant_id': hit['variant_id'],
                'celltype': hit['celltype'],
                'slope': hit['slope'],
                'pval_nominal': hit['pval_nominal'],
            })

        # caQTL
        caqtl_chr = s6_caqtl[s6_caqtl['chrom'] == chr_val]
        caqtl_nearby = caqtl_chr[
            (caqtl_chr['pos_s6'] >= pos_val - 500000) &
            (caqtl_chr['pos_s6'] <= pos_val + 500000)
        ]
        for _, hit in caqtl_nearby.iterrows():
            tanaka_caqtl_hits.append({
                'tanaka_rsid': snp['rsid'],
                'tanaka_chr': chr_val,
                'tanaka_pos': pos_val,
                'phenotype_id': hit['phenotype_id'],
                'variant_id': hit['variant_id'],
                'celltype': hit['celltype'],
                'slope': hit['slope'],
                'pval_nominal': hit['pval_nominal'],
            })

    tanaka_eqtl = pd.DataFrame(tanaka_eqtl_hits)
    tanaka_caqtl = pd.DataFrame(tanaka_caqtl_hits)

    n_eqtl_loci = tanaka_eqtl['tanaka_rsid'].nunique() if len(tanaka_eqtl) > 0 else 0
    n_caqtl_loci = tanaka_caqtl['tanaka_rsid'].nunique() if len(tanaka_caqtl) > 0 else 0

    print(f"\n  Tanaka SNP -> CIMA overlap:")
    print(f"    eQTL: {n_eqtl_loci}/17 loci ({len(tanaka_eqtl)} associations)")
    print(f"    caQTL: {n_caqtl_loci}/17 loci ({len(tanaka_caqtl)} associations)")

    # Compare with European GWAS
    eu_summary = pd.read_csv(summary_path)
    eu_eqtl_rate = eu_summary['has_eqtl'].mean()
    eu_caqtl_rate = eu_summary['has_caqtl'].mean()
    jp_eqtl_rate = n_eqtl_loci / 17
    jp_caqtl_rate = n_caqtl_loci / 17

    print(f"\n  Overlap rate comparison:")
    print(f"    European GWAS: eQTL {eu_eqtl_rate:.1%} ({eu_summary['has_eqtl'].sum()}/25), "
          f"caQTL {eu_caqtl_rate:.1%} ({eu_summary['has_caqtl'].sum()}/25)")
    print(f"    Japanese GWAS: eQTL {jp_eqtl_rate:.1%} ({n_eqtl_loci}/17), "
          f"caQTL {jp_caqtl_rate:.1%} ({n_caqtl_loci}/17)")

    # Cell-type enrichment comparison
    if len(tanaka_eqtl) > 0:
        jp_ct_counts = tanaka_eqtl['celltype'].value_counts()
        print(f"\n  Tanaka top cell types (eQTL):")
        for ct, n in jp_ct_counts.head(10).items():
            print(f"    {ct}: {n}")

        jp_genes = set(tanaka_eqtl['phenotype_id'].unique())
        eu_eqtl = pd.read_csv('results/decode_ad_eqtl.csv')
        eu_genes = set(eu_eqtl['phenotype_id'].unique())
        shared_genes = jp_genes & eu_genes
        print(f"\n  Gene overlap: {len(shared_genes)} shared / "
              f"{len(jp_genes)} Japanese / {len(eu_genes)} European")
        if shared_genes:
            print(f"    Shared: {sorted(shared_genes)[:20]}")

    # Save
    comparison = pd.DataFrame({
        'metric': ['n_loci', 'eqtl_overlap', 'caqtl_overlap', 'eqtl_rate', 'caqtl_rate',
                    'n_eqtl_genes', 'n_eqtl_celltypes'],
        'european': [25, eu_summary['has_eqtl'].sum(), eu_summary['has_caqtl'].sum(),
                     eu_eqtl_rate, eu_caqtl_rate,
                     eu_eqtl['phenotype_id'].nunique() if 'eu_eqtl' in dir() else 195,
                     eu_eqtl['celltype'].nunique() if 'eu_eqtl' in dir() else 70],
        'japanese': [17, n_eqtl_loci, n_caqtl_loci, jp_eqtl_rate, jp_caqtl_rate,
                     tanaka_eqtl['phenotype_id'].nunique() if len(tanaka_eqtl) > 0 else 0,
                     tanaka_eqtl['celltype'].nunique() if len(tanaka_eqtl) > 0 else 0],
    })
    comparison.to_csv('results/decode_ad_tanaka_comparison.csv', index=False)

    # Also save full hit lists
    if len(tanaka_eqtl) > 0:
        tanaka_eqtl.to_csv('results/decode_ad_tanaka_eqtl.csv', index=False)
    if len(tanaka_caqtl) > 0:
        tanaka_caqtl.to_csv('results/decode_ad_tanaka_caqtl.csv', index=False)

    print(f"\nSaved to results/decode_ad_tanaka_comparison.csv")
    return comparison


# ======================================================================
# Main
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description="DECODE-AD Final Analyses")
    parser.add_argument('--s6', default='data/raw/CIMA_Table_S6.csv')
    parser.add_argument('--s4', default='data/raw/CIMA_Table_S4.csv')
    parser.add_argument('--s9', default='data/raw/science.adt3130_table_s9.xlsx')
    parser.add_argument('--s11', default='data/raw/science.adt3130_table_s11.xlsx')
    parser.add_argument('--s12', default='data/raw/science.adt3130_table_s12.csv')
    parser.add_argument('--topple', default='results/topple.csv')
    parser.add_argument('--eqtl', default='results/decode_ad_eqtl.csv')
    parser.add_argument('--caqtl', default='results/decode_ad_caqtl.csv')
    parser.add_argument('--summary', default='results/decode_ad_summary.csv')
    parser.add_argument('--disruption', default='results/decode_ad_disruption_network.csv')
    args = parser.parse_args()

    # Analysis 1: Cell-type PRS
    ct_prs = celltype_prs(args.s6, args.eqtl, args.summary)

    # Analysis 2: Regulatory cascades
    cascades = regulatory_cascades(args.s4, args.s11, args.eqtl)

    # Analysis 3: Network figure
    network_figure(args.disruption)

    # Analysis 4: Drug target deep
    drug_df = drug_target_deep(args.s6, args.topple)

    # Analysis 5: Interaction QTLs
    iqtl = interaction_qtls(args.s12, args.summary, args.eqtl, args.caqtl)

    # Analysis 6: Tanaka comparison
    tanaka = tanaka_comparison(args.s6, args.s9, args.summary)

    # Final summary
    print("\n" + "=" * 60)
    print("DECODE-AD FINAL ANALYSES -- COMPLETE")
    print("=" * 60)
    print(f"  1. Cell-type PRS: {len(ct_prs)} cell types ranked")
    print(f"  2. Regulatory cascades: {len(cascades)} chains found")
    print(f"  3. Network figure: fig_decode_ad_network.pdf")
    print(f"  4. Drug targets: {(drug_df['found_in_s6']).sum()}/{len(drug_df)} found in S6")
    print(f"  5. Interaction QTLs: {len(iqtl)} hits")
    print(f"  6. Tanaka comparison: saved")
    print("=" * 60)


if __name__ == "__main__":
    main()
