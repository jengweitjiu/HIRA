#!/usr/bin/env python3
"""
Thompson et al. (Nature 2026) Rete Ridge Validation.

Apply HIRA framework (TOPPLE, SICAI, STRATA) to pig skin scRNA-seq and
Stereo-seq data from GSE305111 to validate regulon stability principles
in an independent tissue and species.

Input: data/thompson/ — RData and h5ad files from GEO GSE305111
Output: results/thompson_*.csv, figures/fig_thompson_*.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


# ---- Cell type annotations ----
# P3 all clusters (19 clusters from scRNA-seq)
ALLCLUSTER_ANNOTATIONS = {
    '0': 'Krt_Progenitor_KRT19',
    '1': 'Immune_Mixed',
    '2': 'Krt_Differentiated_KRT10',
    '3': 'Krt_HF_Bulge_KRT15',
    '4': 'Krt_Basal_POSTN',
    '5': 'Endothelial',
    '6': 'Krt_Glandular_KRT19',
    '7': 'Krt_HF_Inner_KRT25',
    '8': 'Fibroblast_COL1A1',
    '9': 'Pericyte_RGS5',
    '10': 'Krt_Proliferating',
    '11': 'Krt_KRT17_Appendageal',
    '12': 'Fibroblast_Dermal_SFRP2',
    '13': 'Krt_Cornified_KRT35',
    '14': 'Mesenchymal_VIM',
    '15': 'Krt_Glandular_CHIA',
    '16': 'Melanocyte_KIT',
    '17': 'Lymphatic_CCL21',
    '18': 'Neural_ADAM23',
}

# P3 keratinocyte subset (9 clusters)
KRTNO_ANNOTATIONS = {
    '0': 'Differentiated_KRT10',
    '1': 'Basal_Progenitor',
    '2': 'Proliferating',
    '3': 'Basal_KRT14_POSTN',
    '4': 'KRT17_Appendageal',
    '5': 'HF_Bulge_KRT15',
    '6': 'LGR5_Stem',
    '7': 'EMT_like',
    '8': 'Ribosomal_High',
}


def load_seurat_to_anndata(rdata_path, key=None):
    """Load a Seurat RData file and convert to AnnData."""
    import rdata as rd
    import anndata as ad
    from scipy import sparse
    import gzip

    # Decompress if needed
    if str(rdata_path).endswith('.gz'):
        raw_path = str(rdata_path).replace('.gz', '')
        if not Path(raw_path).exists():
            with gzip.open(rdata_path, 'rb') as f:
                data = f.read()
            with open(raw_path, 'wb') as f:
                f.write(data)
        rdata_path = raw_path

    parsed = rd.parser.parse_file(str(rdata_path))
    converted = rd.conversion.convert(parsed)
    obj_key = key or list(converted.keys())[0]
    seurat = converted[obj_key]

    # Extract SCT assay expression
    sct = seurat.assays.get('SCT', list(seurat.assays.values())[0])
    sp_data = sct.data
    i = np.array(sp_data.i)
    p = np.array(sp_data.p)
    x = np.array(sp_data.x)
    dims = np.array(sp_data.Dim)
    dimnames = sp_data.Dimnames
    mat = sparse.csc_matrix((x, i, p), shape=tuple(dims))

    adata = ad.AnnData(
        X=mat.T,
        obs=seurat.__getattribute__('meta.data').copy(),
        var=pd.DataFrame(index=dimnames[0])
    )
    adata.obs.index = dimnames[1]
    return adata


def step1_topple_tf(adata, cell_type_col='cell_type', min_cells=20):
    """
    TOPPLE-like analysis using TF expression across cell types.

    Instead of eRegulon AUC (CIMA-specific), we use mean TF expression
    per cell type as a proxy for regulon activity, then compute
    leave-one-out JSD perturbation (Redistribution Index).
    """
    import scanpy as sc

    print("\n" + "=" * 60)
    print("STEP 1: TOPPLE — TF Stability in Pig Skin")
    print("=" * 60)

    # Normalize
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # Get known TFs (use a curated list of key skin TFs + CIMA top TFs)
    tf_list = [
        # CIMA top stabilizers
        'HSF1', 'EGR1', 'JUN', 'JUNB', 'JUND', 'FOS', 'FOSB',
        'NFKB1', 'NFKB2', 'RELA', 'REL', 'RELB',
        'STAT1', 'STAT3', 'STAT5A', 'STAT5B',
        'IRF1', 'IRF4', 'IRF8',
        'SPI1', 'CEBPA', 'CEBPB', 'CEBPD',
        # BMP pathway TFs (Thompson key pathway)
        'SMAD1', 'SMAD5', 'SMAD9', 'SMAD4',
        'ID1', 'ID2', 'ID3', 'ID4',
        'MSX1', 'MSX2', 'DLX3',
        # Skin-specific TFs
        'TP63', 'KLF4', 'KLF5', 'GRHL1', 'GRHL3',
        'SOX9', 'SOX2', 'LEF1', 'TCF7L2',
        'GATA3', 'GATA6', 'OVOL1', 'OVOL2',
        'FOXN1', 'FOXO1', 'FOSL1', 'FOSL2',
        'ETS1', 'ETS2', 'ELF3', 'ELF5',
        'NOTCH1', 'NOTCH2', 'HES1', 'HEY1',
        'MYC', 'MYCN', 'MAX',
        'RUNX1', 'RUNX2', 'RUNX3',
        'ZNF750', 'MAF', 'MAFB',
        'NR4A1', 'NR4A2', 'NR4A3',
        'BACH1', 'BACH2', 'NFE2L2',
        'YAP1', 'TEAD1', 'TEAD4',
        'RBPJ', 'CTNNB1',
    ]

    # Filter to TFs present in the data
    present_tfs = [tf for tf in tf_list if tf in adata_norm.var_names]
    print(f"  TFs present in pig data: {len(present_tfs)}/{len(tf_list)}")

    # Get cell types with sufficient cells
    ct_counts = adata_norm.obs[cell_type_col].value_counts()
    valid_cts = ct_counts[ct_counts >= min_cells].index.tolist()
    print(f"  Valid cell types (>={min_cells} cells): {len(valid_cts)}")

    # Build TF expression matrix (TFs × cell types)
    tf_matrix = pd.DataFrame(index=present_tfs, columns=valid_cts, dtype=float)
    for ct in valid_cts:
        mask = adata_norm.obs[cell_type_col] == ct
        for tf in present_tfs:
            vals = adata_norm[mask, tf].X
            if hasattr(vals, 'toarray'):
                vals = vals.toarray()
            tf_matrix.loc[tf, ct] = float(np.mean(vals))

    # Normalize each TF to a distribution across cell types
    tf_dist = tf_matrix.div(tf_matrix.sum(axis=1), axis=0).fillna(0)

    # TOPPLE: Leave-one-out JSD perturbation
    # Only include TFs expressed in ≥5 cell types (otherwise RI is trivially high)
    n_cts = len(valid_cts)
    results = []

    for tf in present_tfs:
        full_dist = tf_dist.loc[tf].values.astype(float)
        if full_dist.sum() == 0:
            continue
        # Require expression in at least 5 cell types
        n_expressed = int((tf_matrix.loc[tf] > 0.1).sum())
        if n_expressed < 5:
            continue

        jsd_values = []
        for j in range(n_cts):
            # Remove cell type j, renormalize
            reduced = np.delete(full_dist, j)
            if reduced.sum() > 0:
                reduced = reduced / reduced.sum()
                # Compare to full distribution without cell type j
                full_reduced = np.delete(full_dist, j)
                full_reduced = full_reduced / full_reduced.sum() if full_reduced.sum() > 0 else full_reduced
                # Actually compare full vs reduced
                full_normed = full_dist / full_dist.sum()
                # Reinsert 0 at position j for alignment
                padded = np.insert(reduced * (1 - full_normed[j]), j, 0)
                padded = padded / padded.sum() if padded.sum() > 0 else padded
                jsd = jensenshannon(full_normed, padded) ** 2
                jsd_values.append(jsd)

        if jsd_values:
            ri = np.mean(jsd_values)
            results.append({
                'TF': tf,
                'mean_RI': ri,
                'max_RI': np.max(jsd_values),
                'mean_expression': float(tf_matrix.loc[tf].mean()),
                'n_cell_types_expressed': int((tf_matrix.loc[tf] > 0.1).sum()),
            })

    results_df = pd.DataFrame(results).sort_values('mean_RI', ascending=False)
    results_df['rank'] = range(1, len(results_df) + 1)

    print(f"\n  Top 10 most stable TFs (highest RI):")
    for _, row in results_df.head(10).iterrows():
        print(f"    #{int(row['rank']):2d} {row['TF']:10s} RI={row['mean_RI']:.5f}  "
              f"expr={row['mean_expression']:.2f}  n_CTs={int(row['n_cell_types_expressed'])}")

    # Check CIMA top stabilizers
    print(f"\n  CIMA stabilizer TF ranks in pig skin:")
    cima_top = ['HSF1', 'EGR1', 'JUN', 'JUNB', 'FOS', 'NFKB1', 'STAT3']
    for tf in cima_top:
        match = results_df[results_df['TF'] == tf]
        if len(match) > 0:
            row = match.iloc[0]
            print(f"    {tf:10s} rank=#{int(row['rank']):2d}/{len(results_df)}  "
                  f"RI={row['mean_RI']:.5f}")
        else:
            print(f"    {tf:10s} not found in data")

    # BMP pathway TFs
    print(f"\n  BMP pathway TF ranks:")
    bmp_tfs = ['SMAD1', 'SMAD5', 'SMAD4', 'ID1', 'ID2', 'ID3', 'MSX2']
    for tf in bmp_tfs:
        match = results_df[results_df['TF'] == tf]
        if len(match) > 0:
            row = match.iloc[0]
            print(f"    {tf:10s} rank=#{int(row['rank']):2d}/{len(results_df)}  "
                  f"RI={row['mean_RI']:.5f}")
        else:
            print(f"    {tf:10s} not found in data")

    return results_df, tf_matrix


def step2_sicai_coupling(adata, cell_type_col='cell_type', min_cells=20, n_hvg=2000):
    """
    SICAI-like coupling analysis between keratinocyte subtypes.

    Compute pairwise gene expression correlation (r_b proxy) between
    cell types using highly variable genes.
    """
    import scanpy as sc

    print("\n" + "=" * 60)
    print("STEP 2: SICAI — Cell-Type Coupling in Pig Skin")
    print("=" * 60)

    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # Get valid cell types
    ct_counts = adata_norm.obs[cell_type_col].value_counts()
    valid_cts = sorted(ct_counts[ct_counts >= min_cells].index.tolist())
    print(f"  Cell types: {len(valid_cts)}")

    # Select HVGs
    sc.pp.highly_variable_genes(adata_norm, n_top_genes=min(n_hvg, adata_norm.shape[1]))
    hvg_mask = adata_norm.var['highly_variable']
    print(f"  HVGs: {hvg_mask.sum()}")

    # Compute mean expression per cell type for HVGs
    mean_expr = pd.DataFrame(index=adata_norm.var_names[hvg_mask], columns=valid_cts, dtype=float)
    for ct in valid_cts:
        mask = adata_norm.obs[cell_type_col] == ct
        X_sub = adata_norm[mask][:, hvg_mask].X
        if hasattr(X_sub, 'toarray'):
            X_sub = X_sub.toarray()
        mean_expr[ct] = X_sub.mean(axis=0)

    # Pairwise Spearman correlation (proxy for r_b)
    n = len(valid_cts)
    r_b_matrix = pd.DataFrame(np.ones((n, n)), index=valid_cts, columns=valid_cts)

    for i in range(n):
        for j in range(i + 1, n):
            rho, _ = stats.spearmanr(mean_expr[valid_cts[i]], mean_expr[valid_cts[j]])
            r_b_matrix.iloc[i, j] = rho
            r_b_matrix.iloc[j, i] = rho

    # Summary statistics (like SICAI)
    print(f"\n  Pairwise coupling (Spearman rho proxy for r_b):")
    mean_rb = r_b_matrix.values[np.triu_indices(n, k=1)].mean()
    print(f"    Mean coupling: {mean_rb:.3f}")
    print(f"    Min coupling: {r_b_matrix.values[np.triu_indices(n, k=1)].min():.3f}")
    print(f"    Max coupling: {r_b_matrix.values[np.triu_indices(n, k=1)].max():.3f}")

    # Per cell-type mean coupling
    ct_coupling = {}
    for ct in valid_cts:
        others = [c for c in valid_cts if c != ct]
        ct_coupling[ct] = r_b_matrix.loc[ct, others].mean()

    ct_coupling_df = pd.DataFrame([
        {'cell_type': ct, 'mean_coupling': v}
        for ct, v in ct_coupling.items()
    ]).sort_values('mean_coupling', ascending=False)

    print(f"\n  Per-cell-type mean coupling:")
    for _, row in ct_coupling_df.iterrows():
        print(f"    {row['cell_type']:30s} r_b={row['mean_coupling']:.3f}")

    # Identify keratinocyte vs non-keratinocyte coupling
    krt_cts = [ct for ct in valid_cts if ct.startswith('Krt_')]
    nonkrt_cts = [ct for ct in valid_cts if not ct.startswith('Krt_')]
    if krt_cts and nonkrt_cts:
        krt_krt = []
        krt_nonkrt = []
        for i, ct1 in enumerate(valid_cts):
            for j, ct2 in enumerate(valid_cts):
                if i >= j:
                    continue
                val = r_b_matrix.loc[ct1, ct2]
                if ct1 in krt_cts and ct2 in krt_cts:
                    krt_krt.append(val)
                elif (ct1 in krt_cts) != (ct2 in krt_cts):
                    krt_nonkrt.append(val)

        if krt_krt and krt_nonkrt:
            stat, p = stats.mannwhitneyu(krt_krt, krt_nonkrt, alternative='greater')
            print(f"\n  Keratinocyte-Keratinocyte coupling: {np.mean(krt_krt):.3f} (n={len(krt_krt)})")
            print(f"  Keratinocyte-NonKrt coupling:       {np.mean(krt_nonkrt):.3f} (n={len(krt_nonkrt)})")
            print(f"  Mann-Whitney P (Krt-Krt > Krt-NonKrt): {p:.2e}")

    return r_b_matrix, ct_coupling_df


def step3_strata_spatial(spatial_path, cell_type_col=None):
    """
    STRATA-like spatial analysis on Stereo-seq data.

    Compute spatial autocorrelation (Moran's I) of BMP pathway genes
    to test whether rete ridge-associated signals show spatial structure.
    """
    import scanpy as sc
    import h5py

    print("\n" + "=" * 60)
    print("STEP 3: STRATA — Spatial Validation (Stereo-seq)")
    print("=" * 60)

    # Read spatial data
    print(f"  Loading {spatial_path}...")
    f = h5py.File(spatial_path, 'r')

    # Extract expression matrix
    data_grp = f['exp_matrix']
    from scipy import sparse
    X = sparse.csc_matrix((
        np.array(data_grp['data']),
        np.array(data_grp['indices']),
        np.array(data_grp['indptr'])
    ))

    # Get gene names and positions
    genes = [g.decode() if isinstance(g, bytes) else str(g)
             for g in f['genes']['var']['_index'][:]]
    positions = np.array(f['position'][:])
    n_spots = X.shape[1]
    n_genes = X.shape[0]

    print(f"  Spots: {n_spots}, Genes: {n_genes}")
    print(f"  Spatial range: x=[{positions[:,0].min()}-{positions[:,0].max()}], "
          f"y=[{positions[:,1].min()}-{positions[:,1].max()}]")

    # Create AnnData
    import anndata as ad
    adata = ad.AnnData(X=X.T)  # spots × genes
    adata.var_names = genes
    adata.obsm['spatial'] = positions.astype(float)

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Test BMP pathway genes for spatial autocorrelation
    bmp_genes = ['BMP2', 'BMP4', 'BMP7', 'BMPR1A', 'BMPR2',
                 'ID1', 'ID2', 'ID3', 'SMAD1', 'SMAD5']
    krt_genes = ['KRT14', 'KRT10', 'KRT1', 'KRT17', 'KRT5', 'KRT15']
    test_genes = bmp_genes + krt_genes

    present = [g for g in test_genes if g in adata.var_names]
    print(f"  Test genes present: {len(present)}/{len(test_genes)}")

    if len(present) > 0:
        # Compute spatial neighbors
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=min(20, n_spots - 1))
        nn.fit(positions)
        distances, indices = nn.kneighbors(positions)

        # Compute Moran's I for each gene
        results = []
        for gene in present:
            gene_idx = list(adata.var_names).index(gene)
            expr = adata[:, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = np.array(expr).flatten()

            if expr.std() == 0:
                continue

            # Moran's I
            expr_centered = expr - expr.mean()
            n = len(expr)
            W = 0  # sum of weights
            numerator = 0
            denominator = np.sum(expr_centered ** 2)

            if denominator == 0:
                continue

            # Use k nearest neighbors as spatial weights
            for i_spot in range(n):
                for j_nb in indices[i_spot, 1:7]:  # 6 nearest neighbors
                    numerator += expr_centered[i_spot] * expr_centered[j_nb]
                    W += 1

            morans_i = (n / W) * (numerator / denominator) if W > 0 else 0

            results.append({
                'gene': gene,
                'morans_i': morans_i,
                'mean_expression': float(expr.mean()),
                'pct_expressed': float((expr > 0).mean() * 100),
                'category': 'BMP_pathway' if gene in bmp_genes else 'Keratin',
            })

        results_df = pd.DataFrame(results).sort_values('morans_i', ascending=False)
        print(f"\n  Spatial autocorrelation (Moran's I):")
        for _, row in results_df.iterrows():
            print(f"    {row['gene']:10s} I={row['morans_i']:.4f}  "
                  f"mean={row['mean_expression']:.2f}  "
                  f"pct={row['pct_expressed']:.1f}%  [{row['category']}]")

        # Compare BMP vs keratin spatial structure
        bmp_morans = results_df[results_df['category'] == 'BMP_pathway']['morans_i'].values
        krt_morans = results_df[results_df['category'] == 'Keratin']['morans_i'].values
        if len(bmp_morans) > 1 and len(krt_morans) > 1:
            stat, p = stats.mannwhitneyu(bmp_morans, krt_morans, alternative='two-sided')
            print(f"\n  BMP vs Keratin Moran's I: MW P={p:.3f}")
            print(f"    BMP mean I={np.mean(bmp_morans):.4f}, Keratin mean I={np.mean(krt_morans):.4f}")

        f.close()
        return results_df
    else:
        print("  No test genes found in spatial data")
        f.close()
        return pd.DataFrame()


def step4_cross_species_comparison(topple_results):
    """
    Compare TOPPLE stability rankings between pig skin (Thompson)
    and human blood (CIMA).
    """
    print("\n" + "=" * 60)
    print("STEP 4: Cross-Species Comparison (Pig Skin vs CIMA Blood)")
    print("=" * 60)

    # Load CIMA TOPPLE results
    cima_path = Path('results/topple_stability_scores.csv')
    if not cima_path.exists():
        print("  CIMA TOPPLE results not found, skipping comparison")
        return None

    cima = pd.read_csv(cima_path, index_col=0)
    cima = cima.reset_index().rename(columns={'index': 'regulon'})
    print(f"  CIMA TOPPLE: {len(cima)} regulons")
    print(f"  Thompson TF-TOPPLE: {len(topple_results)} TFs")

    # Extract TF names from CIMA regulon names (format: "TF_+", "TF_-")
    cima['TF'] = cima['regulon'].str.replace(r'_[+-]$', '', regex=True)

    # Find overlapping TFs
    overlap = set(topple_results['TF']) & set(cima['TF'])
    print(f"  Overlapping TFs: {len(overlap)}")

    if len(overlap) < 5:
        print("  Too few overlapping TFs for comparison")
        return None

    # Merge on TF name (use max RI for CIMA if TF appears as both + and -)
    cima_agg = cima.groupby('TF')['mean_RI'].max().reset_index()
    merged = topple_results[['TF', 'mean_RI']].merge(
        cima_agg, on='TF', how='inner', suffixes=('_pig', '_cima')
    )

    # Spearman correlation of stability rankings
    if merged['mean_RI_pig'].std() == 0 or merged['mean_RI_cima'].std() == 0:
        print(f"\n  Stability correlation: insufficient variance for correlation")
        rho, p = np.nan, np.nan
    else:
        rho, p = stats.spearmanr(merged['mean_RI_pig'], merged['mean_RI_cima'])
        print(f"\n  Stability correlation (n={len(merged)} TFs):")
        print(f"    Spearman rho = {rho:.3f}, P = {p:.2e}")

    # Top stabilizers comparison
    print(f"\n  Top 5 pig stabilizers → CIMA rank:")
    pig_top5 = topple_results.nlargest(5, 'mean_RI')
    for _, row in pig_top5.iterrows():
        tf = row['TF']
        cima_match = cima_agg[cima_agg['TF'] == tf]
        if len(cima_match) > 0:
            cima_ri = cima_match.iloc[0]['mean_RI']
            cima_rank = int((cima_agg['mean_RI'] >= cima_ri).sum())
            total = len(cima_agg)
            print(f"    {tf:10s}  pig_RI={row['mean_RI']:.5f}  cima_rank=#{cima_rank}/{total}")
        else:
            print(f"    {tf:10s}  pig_RI={row['mean_RI']:.5f}  not in CIMA")

    return merged


def step5_figure(topple_results, coupling_matrix, spatial_results, comparison):
    """Generate Thompson validation figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})

    n_panels = 4
    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 150/25.4))

    # Panel A: TF stability ranking (bar chart, top 20)
    ax = axes[0, 0]
    top20 = topple_results.nlargest(20, 'mean_RI')
    colors = []
    cima_tfs = {'HSF1', 'EGR1', 'JUN', 'JUNB', 'FOS', 'NFKB1', 'STAT3', 'IRF1', 'CEBPB'}
    bmp_tfs = {'SMAD1', 'SMAD5', 'SMAD4', 'ID1', 'ID2', 'ID3', 'MSX2'}
    for tf in top20['TF']:
        if tf in cima_tfs:
            colors.append('#E74C3C')
        elif tf in bmp_tfs:
            colors.append('#3498DB')
        else:
            colors.append('#95A5A6')
    ax.barh(range(len(top20)), top20['mean_RI'].values, color=colors)
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20['TF'].values, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel('Redistribution Index')
    ax.set_title('TF Stability (TOPPLE)', fontsize=9, fontweight='bold')
    # Legend
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color='#E74C3C', label='CIMA stabilizer'),
        Patch(color='#3498DB', label='BMP pathway'),
        Patch(color='#95A5A6', label='Other'),
    ], fontsize=5, loc='lower right')

    # Panel B: Cell-type coupling heatmap
    ax = axes[0, 1]
    if coupling_matrix is not None and len(coupling_matrix) > 0:
        n = len(coupling_matrix)
        # Abbreviate labels
        labels = [c[:15] for c in coupling_matrix.index]
        im = ax.imshow(coupling_matrix.values, cmap='RdYlBu_r', vmin=0.5, vmax=1.0)
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=90, fontsize=4)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels, fontsize=4)
        plt.colorbar(im, ax=ax, shrink=0.7, label='Coupling (rho)')
    ax.set_title('Cell-Type Coupling (SICAI)', fontsize=9, fontweight='bold')

    # Panel C: Spatial autocorrelation
    ax = axes[1, 0]
    if spatial_results is not None and len(spatial_results) > 0:
        bmp_mask = spatial_results['category'] == 'BMP_pathway'
        krt_mask = spatial_results['category'] == 'Keratin'
        ax.barh(range(len(spatial_results)), spatial_results['morans_i'].values,
                color=['#3498DB' if b else '#2ECC71' for b in bmp_mask])
        ax.set_yticks(range(len(spatial_results)))
        ax.set_yticklabels(spatial_results['gene'].values, fontsize=6)
        ax.set_xlabel("Moran's I")
        ax.legend(handles=[
            Patch(color='#3498DB', label='BMP pathway'),
            Patch(color='#2ECC71', label='Keratin'),
        ], fontsize=5, loc='lower right')
    ax.set_title('Spatial Structure (STRATA)', fontsize=9, fontweight='bold')

    # Panel D: Cross-species comparison
    ax = axes[1, 1]
    if comparison is not None and len(comparison) > 0:
        ax.scatter(comparison['mean_RI_cima'], comparison['mean_RI_pig'],
                   s=20, c='#2C3E50', alpha=0.6, edgecolors='white', linewidth=0.5)
        # Label top TFs
        for _, row in comparison.nlargest(5, 'mean_RI_pig').iterrows():
            ax.annotate(row['TF'], (row['mean_RI_cima'], row['mean_RI_pig']),
                       fontsize=5, ha='left')
        rho, p = stats.spearmanr(comparison['mean_RI_cima'], comparison['mean_RI_pig'])
        ax.text(0.05, 0.95, f'rho={rho:.3f}\nP={p:.2e}\nn={len(comparison)}',
                transform=ax.transAxes, fontsize=6, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.set_xlabel('CIMA Blood RI')
        ax.set_ylabel('Pig Skin RI')
    ax.set_title('Cross-Species (CIMA vs Thompson)', fontsize=9, fontweight='bold')

    # Panel labels
    for ax, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('HIRA Validation in Pig Skin (Thompson et al. 2026, GSE305111)',
                 fontsize=10, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_thompson_validation.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_thompson_validation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("\nSaved figures/fig_thompson_validation.pdf + .png")


def main():
    print("=" * 60)
    print("Thompson et al. (Nature 2026) — HIRA Rete Ridge Validation")
    print("GSE305111 — Pig Skin scRNA-seq + Stereo-seq")
    print("=" * 60)

    # ---- Load P3 allclusters data ----
    print("\n--- Loading P3 pig skin allclusters ---")
    rdata_path = Path('data/thompson/P3_pig_krtno_sub.RData.gz')
    allclusters_path = Path('data/thompson/p3_pig_allclusters_v3.RData.gz')

    if allclusters_path.exists() and allclusters_path.stat().st_size > 1000:
        adata = load_seurat_to_anndata(allclusters_path, key='P3')
        annotations = ALLCLUSTER_ANNOTATIONS
        print(f"  Loaded allclusters: {adata.shape}")
    elif rdata_path.exists():
        adata = load_seurat_to_anndata(rdata_path, key='P3_krtno_sub')
        annotations = KRTNO_ANNOTATIONS
        print(f"  Loaded keratinocyte subset: {adata.shape}")
    else:
        print("ERROR: No data files found")
        return

    # Annotate cell types
    adata.obs['cluster'] = adata.obs['seurat_clusters'].astype(str)
    adata.obs['cell_type'] = adata.obs['cluster'].map(annotations).fillna('Unknown')
    print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

    # ---- Step 1: TOPPLE ----
    topple_results, tf_matrix = step1_topple_tf(adata, cell_type_col='cell_type')
    topple_results.to_csv('results/thompson_topple_tf.csv', index=False)
    print(f"  Saved results/thompson_topple_tf.csv")

    # ---- Step 2: SICAI ----
    coupling_matrix, coupling_df = step2_sicai_coupling(adata, cell_type_col='cell_type')
    coupling_matrix.to_csv('results/thompson_sicai_coupling.csv')
    coupling_df.to_csv('results/thompson_sicai_per_ct.csv', index=False)
    print(f"  Saved results/thompson_sicai_*.csv")

    # ---- Step 3: STRATA (Stereo-seq spatial) ----
    spatial_path = Path('data/thompson/PigP3_bin20.stereo.h5ad')
    spatial_results = None
    if spatial_path.exists():
        spatial_results = step3_strata_spatial(str(spatial_path))
        if len(spatial_results) > 0:
            spatial_results.to_csv('results/thompson_strata_spatial.csv', index=False)
            print(f"  Saved results/thompson_strata_spatial.csv")

    # ---- Step 4: Cross-species comparison ----
    comparison = step4_cross_species_comparison(topple_results)
    if comparison is not None:
        comparison.to_csv('results/thompson_cross_species.csv', index=False)
        print(f"  Saved results/thompson_cross_species.csv")

    # ---- Step 5: Figure ----
    step5_figure(topple_results, coupling_matrix, spatial_results, comparison)

    print("\n" + "=" * 60)
    print("Thompson validation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
