#!/usr/bin/env python3
"""
DECODE-AD Deep — Extended analyses for AD GWAS x CIMA integration.

Three analyses:
  1. Top AD eGenes cross-referenced with TOPPLE regulons via S4 GRN
  2. cMono_CD14 caQTL pathway analysis via peak-gene linkages in S4
  3. AD-specific SICAI coupling matrix from eQTL slope profiles

Input:  results/decode_ad_eqtl.csv, results/decode_ad_caqtl.csv,
        results/topple.csv, data/raw/CIMA_Table_S4.csv,
        results/sicai_rb_matrix.csv
Output: results/decode_ad_egenes_regulons.csv,
        results/decode_ad_cmono_pathways.csv,
        results/decode_ad_coupling.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
from itertools import combinations


# ---------------------------------------------------------------------------
# Hardcoded immune pathway gene sets for pathway enrichment
# ---------------------------------------------------------------------------
IMMUNE_PATHWAYS = {
    'JAK-STAT signaling': [
        'JAK1', 'JAK2', 'JAK3', 'TYK2', 'STAT1', 'STAT2', 'STAT3', 'STAT4',
        'STAT5A', 'STAT5B', 'STAT6', 'SOCS1', 'SOCS2', 'SOCS3', 'SOCS4',
        'SOCS5', 'SOCS6', 'SOCS7', 'CISH', 'PIAS1', 'PIAS2', 'PIAS3',
        'PIAS4', 'IL6', 'IL6R', 'IL6ST', 'IL2', 'IL2RA', 'IL2RB', 'IL2RG',
        'IL7', 'IL7R', 'IL15', 'IL15RA', 'IL21', 'IL21R', 'IL12A', 'IL12B',
        'IL12RB1', 'IL12RB2', 'IL23A', 'IL23R', 'IFNG', 'IFNGR1', 'IFNGR2',
        'IL4', 'IL4R', 'IL13', 'IL13RA1', 'IL10', 'IL10RA', 'IL10RB',
    ],
    'NF-kB signaling': [
        'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL', 'NFKBIA', 'NFKBIB',
        'NFKBIE', 'IKBKA', 'IKBKB', 'IKBKG', 'CHUK', 'MAP3K7', 'MAP3K14',
        'TRAF1', 'TRAF2', 'TRAF3', 'TRAF5', 'TRAF6', 'RIPK1', 'BIRC2',
        'BIRC3', 'TNFAIP3', 'CYLD', 'BCL3', 'LTA', 'LTB', 'LTBR',
        'CD40', 'CD40LG', 'TNFRSF1A', 'TNFRSF1B', 'TICAM1', 'MYD88',
        'IRAK1', 'IRAK2', 'IRAK4', 'TAB1', 'TAB2', 'TAB3',
    ],
    'TNF signaling': [
        'TNF', 'TNFRSF1A', 'TNFRSF1B', 'TRADD', 'FADD', 'RIPK1', 'RIPK3',
        'TRAF1', 'TRAF2', 'TRAF3', 'TRAF5', 'CASP3', 'CASP7', 'CASP8',
        'CASP10', 'BIRC2', 'BIRC3', 'CFLAR', 'TNFAIP3', 'MAP3K7', 'MAPK8',
        'MAPK9', 'MAPK14', 'JUNB', 'JUN', 'FOS', 'ICAM1', 'VCAM1',
        'SELE', 'CCL2', 'CCL5', 'CCL20', 'CXCL1', 'CXCL2', 'CXCL3',
        'CXCL10', 'IL6', 'IL1B', 'CSF2', 'MMP9', 'MMP14',
    ],
    'T cell receptor signaling': [
        'CD3D', 'CD3E', 'CD3G', 'CD247', 'ZAP70', 'LCK', 'FYN', 'LAT',
        'SLP76', 'ITK', 'PLCG1', 'RASGRP1', 'VAV1', 'NCK1', 'GRAP2',
        'CARD11', 'BCL10', 'MALT1', 'NFATC1', 'NFATC2', 'NFATC3',
        'NFKB1', 'RELA', 'AP1', 'JUN', 'FOS', 'MAPK1', 'MAPK3',
        'MAP2K1', 'MAP2K2', 'RAF1', 'KRAS', 'HRAS', 'NRAS', 'CTLA4',
        'PDCD1', 'ICOS', 'CD28', 'PIK3CA', 'PIK3CD', 'AKT1', 'AKT2',
    ],
    'Cytokine-cytokine receptor interaction': [
        'IL1A', 'IL1B', 'IL1R1', 'IL1R2', 'IL1RN', 'IL2', 'IL2RA', 'IL2RB',
        'IL3', 'IL3RA', 'IL4', 'IL4R', 'IL5', 'IL5RA', 'IL6', 'IL6R',
        'IL7', 'IL7R', 'IL9', 'IL9R', 'IL10', 'IL10RA', 'IL12A', 'IL12B',
        'IL13', 'IL13RA1', 'IL15', 'IL15RA', 'IL17A', 'IL17F', 'IL17RA',
        'IL18', 'IL18R1', 'IL21', 'IL21R', 'IL22', 'IL22RA1', 'IL23A',
        'IL23R', 'IL33', 'IL33R', 'IFNG', 'IFNGR1', 'TGFB1', 'TGFBR1',
        'CSF1', 'CSF1R', 'CSF2', 'CSF2RA', 'CSF3', 'CSF3R',
    ],
    'Antigen processing and presentation': [
        'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1',
        'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'B2M', 'TAP1', 'TAP2',
        'TAPBP', 'CALR', 'CANX', 'PDIA3', 'PSME1', 'PSME2', 'PSME3',
        'PSMB8', 'PSMB9', 'PSMB10', 'CIITA', 'RFXANK', 'RFX5', 'RFXAP',
        'CD74', 'CTSS', 'CTSL', 'LGMN', 'GILT', 'CD1A', 'CD1B', 'CD1C',
        'CD1D', 'MR1',
    ],
    'Th1/Th2 differentiation': [
        'TBX21', 'GATA3', 'STAT1', 'STAT4', 'STAT6', 'IL12RB1', 'IL12RB2',
        'IL4R', 'IFNG', 'IL4', 'IL5', 'IL13', 'IL2', 'IL12A', 'IL12B',
        'IL18', 'IL18R1', 'RUNX3', 'EOMES', 'IRF1', 'IRF4', 'MAF',
        'BATF', 'BCL6', 'NFATC1', 'NFATC2', 'GATA1', 'IL27', 'IL27RA',
        'IL21', 'CD28', 'ICOS', 'CTLA4',
    ],
    'Th17 differentiation': [
        'RORC', 'RORA', 'STAT3', 'IL6', 'IL6R', 'IL21', 'IL21R', 'IL23A',
        'IL23R', 'TGFB1', 'TGFBR1', 'TGFBR2', 'IL17A', 'IL17F', 'IL22',
        'CCR6', 'CCL20', 'BATF', 'IRF4', 'AHR', 'HIF1A', 'RUNX1',
        'FOXP3', 'SMAD2', 'SMAD3', 'SMAD7', 'IL1B', 'IL1R1', 'IL2',
    ],
    'Chemokine signaling': [
        'CCR1', 'CCR2', 'CCR3', 'CCR4', 'CCR5', 'CCR6', 'CCR7', 'CCR8',
        'CCR9', 'CCR10', 'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CXCR5',
        'CXCR6', 'CX3CR1', 'XCR1', 'CCL1', 'CCL2', 'CCL3', 'CCL4',
        'CCL5', 'CCL7', 'CCL8', 'CCL11', 'CCL13', 'CCL17', 'CCL19',
        'CCL20', 'CCL21', 'CCL22', 'CCL25', 'CXCL1', 'CXCL2', 'CXCL3',
        'CXCL8', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL12', 'CXCL13',
        'CXCL16', 'CX3CL1',
    ],
    'Toll-like receptor signaling': [
        'TLR1', 'TLR2', 'TLR3', 'TLR4', 'TLR5', 'TLR6', 'TLR7', 'TLR8',
        'TLR9', 'TLR10', 'MYD88', 'TIRAP', 'TICAM1', 'TICAM2', 'IRAK1',
        'IRAK2', 'IRAK4', 'TRAF3', 'TRAF6', 'TAB1', 'TAB2', 'MAP3K7',
        'MAPK8', 'MAPK14', 'IRF3', 'IRF5', 'IRF7', 'NFKB1', 'RELA',
        'IFNB1', 'IFNA1', 'TNF', 'IL6', 'IL1B', 'IL12A', 'IL12B',
        'CXCL10', 'CD14', 'LBP', 'CD86',
    ],
    'IL-17 signaling': [
        'IL17A', 'IL17B', 'IL17C', 'IL17D', 'IL17E', 'IL17F', 'IL17RA',
        'IL17RB', 'IL17RC', 'IL17RD', 'IL17RE', 'ACT1', 'TRAF3IP2',
        'TRAF6', 'NFKB1', 'RELA', 'MAPK1', 'MAPK3', 'MAPK8', 'MAPK14',
        'JUN', 'FOS', 'CEBPB', 'CEBPD', 'CXCL1', 'CXCL2', 'CXCL3',
        'CXCL5', 'CXCL8', 'CCL2', 'CCL7', 'CCL20', 'CSF2', 'CSF3',
        'MMP1', 'MMP3', 'MMP9', 'MMP13', 'S100A7', 'S100A8', 'S100A9',
        'DEFB4A', 'LCN2',
    ],
}


# ---------------------------------------------------------------------------
# Analysis 1: Top AD eGenes x TOPPLE regulons
# ---------------------------------------------------------------------------
def analysis1_egenes_regulons(eqtl_path, s4_path, topple_path):
    """Find top AD eGenes, map to TFs via GRN, cross-ref with TOPPLE."""
    print("\n" + "=" * 70)
    print("ANALYSIS 1: Top AD eGenes cross-referenced with TOPPLE regulons")
    print("=" * 70)

    # Load AD eQTL matches
    eqtl = pd.read_csv(eqtl_path)
    print(f"\nLoaded {eqtl_path}")
    print(f"  Columns: {list(eqtl.columns)}")
    print(f"  Rows: {len(eqtl):,}")

    # Find top 10 eGenes by number of unique cell types
    gene_ct_counts = eqtl.groupby('phenotype_id')['celltype'].nunique().sort_values(ascending=False)
    top10_genes = gene_ct_counts.head(10)
    print(f"\nTop 10 AD eGenes by cell-type breadth:")
    for gene, n_ct in top10_genes.items():
        print(f"  {gene:20s}  {n_ct} cell types")

    top10_gene_list = list(top10_genes.index)

    # Load S4 GRN (only Gene and TF columns needed)
    print(f"\nLoading S4 GRN from {s4_path} ...")
    s4 = pd.read_csv(s4_path, usecols=['TF', 'Gene'])
    print(f"  Columns available: TF, Gene")
    print(f"  Rows: {len(s4):,}")
    print(f"  Unique TFs: {s4['TF'].nunique()}")
    print(f"  Unique Genes: {s4['Gene'].nunique()}")

    # Find TFs that regulate each top AD eGene
    gene_tf_map = {}
    for gene in top10_gene_list:
        tfs = sorted(s4[s4['Gene'] == gene]['TF'].unique())
        gene_tf_map[gene] = tfs

    # Load TOPPLE results
    topple = pd.read_csv(topple_path)
    print(f"\nLoaded {topple_path}")
    print(f"  Columns: {list(topple.columns)}")
    print(f"  Rows: {len(topple)}")

    # Create TF -> stability_class lookup (strip _+ suffix from regulon)
    topple['tf_name'] = topple['regulon'].str.replace(r'_\+$', '', regex=True)
    tf_stability = dict(zip(topple['tf_name'], topple['stability_class']))

    # Build results table
    results_rows = []
    stabilizer_count = 0
    destabilizer_count = 0
    total_mapped = 0

    print(f"\n{'Gene':<20s} {'#CTs':<6s} {'#TFs':<6s} {'TFs (stability_class)'}")
    print("-" * 90)

    for gene in top10_gene_list:
        tfs = gene_tf_map[gene]
        n_cts = int(top10_genes[gene])
        tf_info_list = []
        for tf in tfs:
            sc = tf_stability.get(tf, 'not_in_TOPPLE')
            tf_info_list.append(f"{tf}({sc})")
            if sc == 'stabilizer':
                stabilizer_count += 1
                total_mapped += 1
            elif sc == 'destabilizer':
                destabilizer_count += 1
                total_mapped += 1

        tf_info_str = ', '.join(tf_info_list) if tf_info_list else 'none'
        print(f"  {gene:<20s} {n_cts:<6d} {len(tfs):<6d} {tf_info_str}")

        results_rows.append({
            'gene': gene,
            'n_celltypes': n_cts,
            'n_regulating_TFs': len(tfs),
            'TFs': ';'.join(tfs) if tfs else 'none',
            'TF_stability': ';'.join(
                [tf_stability.get(tf, 'not_in_TOPPLE') for tf in tfs]
            ) if tfs else 'none',
            'n_stabilizer_TFs': sum(1 for tf in tfs if tf_stability.get(tf) == 'stabilizer'),
            'n_destabilizer_TFs': sum(1 for tf in tfs if tf_stability.get(tf) == 'destabilizer'),
            'n_unmapped_TFs': sum(1 for tf in tfs if tf not in tf_stability),
        })

    results_df = pd.DataFrame(results_rows)

    # Summary statistics
    print(f"\n--- Summary ---")
    print(f"  Total TF-gene links mapped to TOPPLE: {total_mapped}")
    print(f"  Stabilizer TFs regulating top AD eGenes: {stabilizer_count}")
    print(f"  Destabilizer TFs regulating top AD eGenes: {destabilizer_count}")
    if total_mapped > 0:
        stab_frac = stabilizer_count / total_mapped
        print(f"  Stabilizer fraction: {stab_frac:.3f} ({stabilizer_count}/{total_mapped})")
        # Binomial test: is stabilizer fraction significantly different from 50%?
        binom_result = stats.binomtest(stabilizer_count, total_mapped, 0.5)
        binom_p = binom_result.pvalue
        print(f"  Binomial test (H0: 50/50): P = {binom_p:.4e}")

        # Also compare to background ratio in TOPPLE
        n_stab_bg = sum(1 for v in tf_stability.values() if v == 'stabilizer')
        n_destab_bg = sum(1 for v in tf_stability.values() if v == 'destabilizer')
        bg_frac = n_stab_bg / (n_stab_bg + n_destab_bg) if (n_stab_bg + n_destab_bg) > 0 else 0.5
        print(f"  Background stabilizer fraction in TOPPLE: {bg_frac:.3f} ({n_stab_bg}/{n_stab_bg + n_destab_bg})")
        binom_result2 = stats.binomtest(stabilizer_count, total_mapped, bg_frac)
        binom_p2 = binom_result2.pvalue
        print(f"  Binomial test (H0: background ratio): P = {binom_p2:.4e}")

    return results_df


# ---------------------------------------------------------------------------
# Analysis 2: cMono_CD14 caQTL pathway analysis
# ---------------------------------------------------------------------------
def analysis2_cmono_pathways(caqtl_path, s4_path):
    """Map cMono_CD14 AD peaks to target genes, run pathway enrichment."""
    print("\n" + "=" * 70)
    print("ANALYSIS 2: cMono_CD14 caQTL pathway analysis")
    print("=" * 70)

    # Load AD caQTL matches
    caqtl = pd.read_csv(caqtl_path)
    print(f"\nLoaded {caqtl_path}")
    print(f"  Columns: {list(caqtl.columns)}")
    print(f"  Rows: {len(caqtl):,}")

    # Filter for cMono_CD14
    cmono = caqtl[caqtl['celltype'] == 'cMono_CD14'].copy()
    print(f"\n  cMono_CD14 caQTL entries: {len(cmono):,}")
    cmono_peaks = cmono['phenotype_id'].unique()
    print(f"  Unique peaks in cMono_CD14: {len(cmono_peaks)}")

    # Load S4 GRN
    print(f"\nLoading S4 GRN from {s4_path} ...")
    s4 = pd.read_csv(s4_path, usecols=['Region', 'Gene'])
    print(f"  Rows: {len(s4):,}")
    print(f"  Unique Regions: {s4['Region'].nunique()}")
    print(f"  Unique Genes: {s4['Gene'].nunique()}")

    # Match peaks to S4 Regions
    # Peak format in caQTL: chr1:151857303-151857804
    # Region format in S4:  chr1:107132199-107132700
    # They should match directly
    matched_s4 = s4[s4['Region'].isin(cmono_peaks)]
    target_genes = sorted(matched_s4['Gene'].unique())
    matched_peaks = matched_s4['Region'].nunique()

    print(f"\n  Peaks matched in S4 GRN: {matched_peaks} / {len(cmono_peaks)}")
    print(f"  Target genes linked to cMono_CD14 AD peaks: {len(target_genes)}")
    if len(target_genes) <= 30:
        print(f"  Genes: {', '.join(target_genes)}")
    else:
        print(f"  First 30 genes: {', '.join(target_genes[:30])}")

    # Background: all unique genes in S4
    all_genes_s4 = s4['Gene'].unique()
    N_background = len(all_genes_s4)
    print(f"  Background gene universe (all S4 genes): {N_background}")

    # Pathway enrichment using hypergeometric test
    target_gene_set = set(target_genes)
    n_target = len(target_gene_set)

    pathway_results = []
    print(f"\n  Pathway enrichment (hypergeometric test):")
    print(f"  {'Pathway':<45s} {'Overlap':<10s} {'PathSize':<10s} {'P-value':<12s} {'Genes'}")
    print("  " + "-" * 110)

    for pathway_name, pathway_genes in IMMUNE_PATHWAYS.items():
        pathway_gene_set = set(pathway_genes)
        # Genes in pathway that are also in S4 background
        pathway_in_bg = pathway_gene_set & set(all_genes_s4)
        K = len(pathway_in_bg)  # pathway genes in background
        if K == 0:
            continue

        overlap = target_gene_set & pathway_in_bg
        k = len(overlap)  # hits

        # Hypergeometric test: P(X >= k)
        # scipy.stats.hypergeom.sf(k-1, N, K, n)
        pval = stats.hypergeom.sf(k - 1, N_background, K, n_target) if k > 0 else 1.0

        overlap_genes = sorted(overlap)
        pathway_results.append({
            'pathway': pathway_name,
            'overlap_count': k,
            'pathway_size_in_bg': K,
            'target_genes': n_target,
            'background_size': N_background,
            'p_value': pval,
            'overlap_genes': ';'.join(overlap_genes) if overlap_genes else '',
        })
        print(f"  {pathway_name:<45s} {k:<10d} {K:<10d} {pval:<12.4e} {', '.join(overlap_genes) if overlap_genes else '-'}")

    pathway_df = pd.DataFrame(pathway_results).sort_values('p_value').reset_index(drop=True)

    # Apply Bonferroni correction
    n_tests = len(pathway_df)
    pathway_df['p_bonferroni'] = np.minimum(pathway_df['p_value'] * n_tests, 1.0)

    print(f"\n  Results sorted by P-value:")
    for _, row in pathway_df.iterrows():
        sig = '*' if row['p_bonferroni'] < 0.05 else ''
        print(f"    {row['pathway']:<45s}  overlap={row['overlap_count']}/{row['pathway_size_in_bg']}  "
              f"P={row['p_value']:.4e}  P_bonf={row['p_bonferroni']:.4e} {sig}")

    return pathway_df


# ---------------------------------------------------------------------------
# Analysis 3: AD-specific SICAI coupling
# ---------------------------------------------------------------------------
def analysis3_ad_coupling(eqtl_path, sicai_matrix_path):
    """Build AD-specific cell-type coupling from eQTL slope profiles."""
    print("\n" + "=" * 70)
    print("ANALYSIS 3: AD-specific SICAI coupling")
    print("=" * 70)

    # Load AD eQTL matches
    eqtl = pd.read_csv(eqtl_path)
    print(f"\nLoaded {eqtl_path}")
    print(f"  Columns: {list(eqtl.columns)}")
    print(f"  Rows: {len(eqtl):,}")

    # For each AD locus x gene, get slope across cell types
    # Use locus_id + phenotype_id as the unit (each gene at each locus)
    eqtl['locus_gene'] = eqtl['locus_id'] + '::' + eqtl['phenotype_id']
    n_locus_genes = eqtl['locus_gene'].nunique()
    print(f"  Unique locus-gene combinations: {n_locus_genes}")

    # Get all cell types present
    all_celltypes = sorted(eqtl['celltype'].unique())
    print(f"  Cell types in AD eQTL data: {len(all_celltypes)}")

    # Build slope matrix: rows = locus_gene, columns = cell types
    # For each locus-gene, take the slope per cell type (if multiple variants,
    # take the one with smallest P-value)
    best_per_ct = eqtl.sort_values('pval_nominal').drop_duplicates(
        subset=['locus_gene', 'celltype'], keep='first'
    )

    slope_pivot = best_per_ct.pivot_table(
        index='locus_gene', columns='celltype', values='slope', fill_value=0
    )
    print(f"  Slope matrix shape: {slope_pivot.shape}")

    # Filter: keep locus-genes present in >= 3 cell types (non-zero)
    nonzero_counts = (slope_pivot != 0).sum(axis=1)
    slope_filtered = slope_pivot[nonzero_counts >= 3]
    print(f"  Locus-genes with >= 3 cell types: {len(slope_filtered)}")

    if len(slope_filtered) < 5:
        print("  Too few locus-genes for coupling analysis. Skipping.")
        return pd.DataFrame()

    # Compute cell-type pair coupling: Pearson correlation of slope profiles
    # For each pair of cell types, correlate their slope vectors across locus-genes
    ct_list = list(slope_filtered.columns)
    n_ct = len(ct_list)
    print(f"  Computing pairwise coupling for {n_ct} cell types ...")

    coupling_records = []
    coupling_matrix = pd.DataFrame(np.nan, index=ct_list, columns=ct_list)

    for i, ct1 in enumerate(ct_list):
        for j, ct2 in enumerate(ct_list):
            if i >= j:
                continue
            v1 = slope_filtered[ct1].values
            v2 = slope_filtered[ct2].values
            # Only use locus-genes where both cell types have non-zero slopes
            mask = (v1 != 0) & (v2 != 0)
            n_shared = mask.sum()
            if n_shared >= 5:
                r, p = stats.pearsonr(v1[mask], v2[mask])
            else:
                r, p = np.nan, np.nan
                n_shared = 0

            coupling_matrix.loc[ct1, ct2] = r
            coupling_matrix.loc[ct2, ct1] = r
            coupling_records.append({
                'cell_type_1': ct1,
                'cell_type_2': ct2,
                'ad_coupling_r': r,
                'ad_coupling_p': p,
                'n_shared_locus_genes': n_shared,
            })

    np.fill_diagonal(coupling_matrix.values, 1.0)
    coupling_df = pd.DataFrame(coupling_records).dropna(subset=['ad_coupling_r'])
    coupling_df = coupling_df.sort_values('ad_coupling_r', ascending=False).reset_index(drop=True)

    print(f"\n  Valid cell-type pairs: {len(coupling_df)}")
    print(f"  Mean AD coupling (Pearson r): {coupling_df['ad_coupling_r'].mean():.4f}")
    print(f"  Median AD coupling: {coupling_df['ad_coupling_r'].median():.4f}")

    # Top 10 most strongly coupled pairs
    print(f"\n  Top 10 most strongly coupled cell-type pairs (AD-specific):")
    for _, row in coupling_df.head(10).iterrows():
        print(f"    {row['cell_type_1']:35s} x {row['cell_type_2']:35s}  "
              f"r={row['ad_coupling_r']:.4f}  P={row['ad_coupling_p']:.2e}  "
              f"n={int(row['n_shared_locus_genes'])}")

    # Compare with global SICAI coupling (r_b matrix)
    print(f"\nLoading global SICAI r_b matrix from {sicai_matrix_path} ...")
    sicai_rb = pd.read_csv(sicai_matrix_path, index_col=0)
    print(f"  SICAI r_b matrix shape: {sicai_rb.shape}")
    print(f"  Columns (first 5): {list(sicai_rb.columns[:5])}")

    # For each AD coupling pair, look up global r_b
    enrichment_records = []
    for _, row in coupling_df.iterrows():
        ct1 = row['cell_type_1']
        ct2 = row['cell_type_2']
        ad_r = row['ad_coupling_r']

        global_rb = np.nan
        if ct1 in sicai_rb.index and ct2 in sicai_rb.columns:
            global_rb = sicai_rb.loc[ct1, ct2]
        elif ct2 in sicai_rb.index and ct1 in sicai_rb.columns:
            global_rb = sicai_rb.loc[ct2, ct1]

        diff = ad_r - global_rb if not np.isnan(global_rb) else np.nan
        enrichment_records.append({
            'cell_type_1': ct1,
            'cell_type_2': ct2,
            'ad_coupling_r': ad_r,
            'global_rb': global_rb,
            'ad_minus_global': diff,
            'n_shared_locus_genes': int(row['n_shared_locus_genes']),
        })

    enrichment_df = pd.DataFrame(enrichment_records).dropna(subset=['ad_minus_global'])
    enrichment_df = enrichment_df.sort_values('ad_minus_global', ascending=False).reset_index(drop=True)

    print(f"\n  Cell-type pairs with both AD and global coupling: {len(enrichment_df)}")

    # Top 10 AD-enriched pairs (AD coupling >> global coupling)
    print(f"\n  Top 10 cell-type pairs with AD-SPECIFIC coupling enrichment:")
    print(f"  (Largest positive AD_r - global_r_b)")
    for _, row in enrichment_df.head(10).iterrows():
        print(f"    {row['cell_type_1']:35s} x {row['cell_type_2']:35s}  "
              f"AD_r={row['ad_coupling_r']:.4f}  global_rb={row['global_rb']:.4f}  "
              f"diff={row['ad_minus_global']:+.4f}")

    # Bottom 10 (AD-depleted coupling)
    print(f"\n  Top 10 cell-type pairs with AD-DEPLETED coupling:")
    for _, row in enrichment_df.tail(10).iterrows():
        print(f"    {row['cell_type_1']:35s} x {row['cell_type_2']:35s}  "
              f"AD_r={row['ad_coupling_r']:.4f}  global_rb={row['global_rb']:.4f}  "
              f"diff={row['ad_minus_global']:+.4f}")

    # Overall correlation between AD coupling and global coupling
    valid = enrichment_df.dropna(subset=['ad_coupling_r', 'global_rb'])
    if len(valid) >= 10:
        corr_r, corr_p = stats.pearsonr(valid['ad_coupling_r'], valid['global_rb'])
        print(f"\n  Correlation between AD coupling and global r_b:")
        print(f"    Pearson r = {corr_r:.4f}, P = {corr_p:.2e} (n = {len(valid)} pairs)")
        rho, rho_p = stats.spearmanr(valid['ad_coupling_r'], valid['global_rb'])
        print(f"    Spearman rho = {rho:.4f}, P = {rho_p:.2e}")

    # Mean difference
    mean_diff = enrichment_df['ad_minus_global'].mean()
    median_diff = enrichment_df['ad_minus_global'].median()
    print(f"\n  Mean (AD - global) difference: {mean_diff:+.4f}")
    print(f"  Median (AD - global) difference: {median_diff:+.4f}")

    # One-sample t-test: is difference != 0?
    t_stat, t_p = stats.ttest_1samp(enrichment_df['ad_minus_global'].dropna(), 0)
    print(f"  One-sample t-test (H0: diff = 0): t = {t_stat:.3f}, P = {t_p:.2e}")

    return enrichment_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='DECODE-AD Deep: extended AD GWAS x CIMA analyses')
    parser.add_argument('--eqtl', default='results/decode_ad_eqtl.csv',
                        help='Path to AD eQTL matches')
    parser.add_argument('--caqtl', default='results/decode_ad_caqtl.csv',
                        help='Path to AD caQTL matches')
    parser.add_argument('--s4', default='data/raw/CIMA_Table_S4.csv',
                        help='Path to S4 GRN table')
    parser.add_argument('--topple', default='results/topple.csv',
                        help='Path to TOPPLE results')
    parser.add_argument('--sicai-matrix', default='results/sicai_rb_matrix.csv',
                        help='Path to SICAI r_b matrix')
    parser.add_argument('--output-dir', default='results',
                        help='Output directory')
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Analysis 1 ---
    egenes_df = analysis1_egenes_regulons(args.eqtl, args.s4, args.topple)
    out1 = outdir / 'decode_ad_egenes_regulons.csv'
    egenes_df.to_csv(out1, index=False)
    print(f"\nSaved -> {out1}")

    # --- Analysis 2 ---
    pathways_df = analysis2_cmono_pathways(args.caqtl, args.s4)
    out2 = outdir / 'decode_ad_cmono_pathways.csv'
    pathways_df.to_csv(out2, index=False)
    print(f"\nSaved -> {out2}")

    # --- Analysis 3 ---
    coupling_df = analysis3_ad_coupling(args.eqtl, args.sicai_matrix)
    out3 = outdir / 'decode_ad_coupling.csv'
    coupling_df.to_csv(out3, index=False)
    print(f"\nSaved -> {out3}")

    print("\n" + "=" * 70)
    print("DECODE-AD Deep analysis complete.")
    print("=" * 70)


if __name__ == '__main__':
    main()
