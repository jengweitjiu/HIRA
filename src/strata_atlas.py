#!/usr/bin/env python3
"""
STRATA-Atlas — Cross-Tissue Mantel Test for Coupling Conservation

Tests whether inter-cellular coupling patterns observed in blood (CIMA eQTL r_b)
are preserved in skin tissue (psoriasis Visium spatial co-expression).

Pipeline:
  1. Load S5 AUC matrix → select top 5 marker regulons per cell type (by AUC)
  2. Extract TF name from each marker regulon as proxy for cell-type activity
  3. Score each Visium spot for each cell type using TF expression
  4. Compute Pearson correlation of TF expression profiles across all spots
     (pooled across samples) → spatial co-expression coupling matrix
  5. Load eQTL r_b from S8 (same build_rb_matrix as sicai.py)
  6. Mantel test between the two matrices on shared cell types
  7. Print key statistics and save results

Input:  data/raw/CIMA_Table_S5.xlsx (sheet 1: regulon × cell_type AUC)
        data/raw/science.adt3130_table_s8.xlsx (sheet cis_eQTL: pairwise r_b)
        data/visium/*.h5 (10x Visium samples)
Output: results/strata_atlas.csv
"""

import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import squareform
from pathlib import Path

warnings.filterwarnings('ignore', category=UserWarning)


# ---------------------------------------------------------------------------
# Step 1: Select top marker regulons per cell type from S5
# ---------------------------------------------------------------------------

def select_marker_regulons(s5_df, top_n=5):
    """Select top N marker regulons per cell type by mean_AUC.

    Returns dict: cell_type -> list of TF names.
    """
    # Use mean_AUC as the ranking metric
    if 'mean_AUC' in s5_df.columns:
        score_col = 'mean_AUC'
    elif 'RSS' in s5_df.columns:
        score_col = 'RSS'
    else:
        raise ValueError(f"No AUC or RSS column found. Columns: {list(s5_df.columns)}")

    print(f"  Ranking regulons by '{score_col}'")

    markers = {}
    for ct, grp in s5_df.groupby('cell_type'):
        top = grp.nlargest(top_n, score_col)
        # Extract TF name: strip trailing _+ or _- from regulon name
        tfs = top['regulon'].str.replace(r'_[+-]$', '', regex=True).tolist()
        markers[ct] = tfs

    print(f"  Cell types with markers: {len(markers)}")
    print(f"  Example — {list(markers.keys())[0]}: {markers[list(markers.keys())[0]]}")
    return markers


# ---------------------------------------------------------------------------
# Step 2-3: Score Visium spots and build co-expression matrix
# ---------------------------------------------------------------------------

def score_visium_spots(visium_dir, markers):
    """For each cell type, compute mean TF expression per spot across all samples.

    Returns DataFrame: rows = spots (pooled), columns = cell types.
    """
    import scanpy as sc

    h5_files = sorted(Path(visium_dir).glob('*.h5'))
    print(f"\n  Found {len(h5_files)} Visium .h5 files")

    # Collect all unique TFs we need
    all_tfs = set()
    for tfs in markers.values():
        all_tfs.update(tfs)
    print(f"  Unique TFs needed: {len(all_tfs)}")

    # Pool spot-level expression across all samples
    all_spot_scores = []

    for h5_path in h5_files:
        print(f"    Loading {h5_path.name} ... ", end='')
        adata = sc.read_10x_h5(str(h5_path))
        adata.var_names_make_unique()

        # Convert to dense if sparse
        try:
            X = adata.X.toarray()
        except AttributeError:
            X = adata.X

        genes = list(adata.var_names)
        gene_idx = {g: i for i, g in enumerate(genes)}

        # Compute per-spot score for each cell type
        n_spots = adata.n_obs
        ct_scores = {}
        for ct, tfs in markers.items():
            valid_tfs = [tf for tf in tfs if tf in gene_idx]
            if len(valid_tfs) == 0:
                ct_scores[ct] = np.zeros(n_spots)
            else:
                cols = [gene_idx[tf] for tf in valid_tfs]
                ct_scores[ct] = X[:, cols].mean(axis=1)

        spot_df = pd.DataFrame(ct_scores)
        all_spot_scores.append(spot_df)
        print(f"{n_spots} spots")

    pooled = pd.concat(all_spot_scores, ignore_index=True)
    print(f"\n  Total pooled spots: {len(pooled)}")
    return pooled


def build_spatial_coupling_matrix(spot_scores):
    """Pearson correlation between cell-type expression profiles across spots.

    Returns symmetric DataFrame (cell_type x cell_type).
    """
    # Drop cell types with zero variance (no expression detected)
    variances = spot_scores.var()
    active_cts = variances[variances > 0].index.tolist()
    dropped = variances[variances == 0].index.tolist()
    if dropped:
        print(f"  Dropped {len(dropped)} cell types with zero variance")

    spot_active = spot_scores[active_cts]
    corr_mat = spot_active.corr(method='pearson')
    print(f"  Spatial coupling matrix: {corr_mat.shape[0]} x {corr_mat.shape[1]} cell types")
    return corr_mat


# ---------------------------------------------------------------------------
# Step 5: Load eQTL r_b matrix (reused from sicai.py)
# ---------------------------------------------------------------------------

def build_rb_matrix(s8_df):
    """Build symmetric r_b matrix from pairwise data."""
    celltypes = sorted(set(s8_df['reference_cell_type']) | set(s8_df['query_celltype']))
    n = len(celltypes)
    ct_idx = {ct: i for i, ct in enumerate(celltypes)}

    rb_mat = np.full((n, n), np.nan)
    np.fill_diagonal(rb_mat, 1.0)

    for _, row in s8_df.iterrows():
        i = ct_idx[row['reference_cell_type']]
        j = ct_idx[row['query_celltype']]
        rb_mat[j, i] = row['rb']  # query is row, reference is column
        if np.isnan(rb_mat[i, j]):
            rb_mat[i, j] = row['rb']  # fill symmetric if missing

    return pd.DataFrame(rb_mat, index=celltypes, columns=celltypes)


# ---------------------------------------------------------------------------
# Step 6: Mantel test
# ---------------------------------------------------------------------------

def mantel_test(mat_a, mat_b, permutations=9999):
    """Mantel test between two symmetric distance/similarity matrices.

    Uses upper-triangle elements. Returns observed Pearson r, P-value.
    """
    n = mat_a.shape[0]
    assert mat_a.shape == mat_b.shape, "Matrices must be same size"

    # Extract upper triangle (excluding diagonal)
    iu = np.triu_indices(n, k=1)
    vec_a = mat_a[iu]
    vec_b = mat_b[iu]

    # Remove NaN pairs
    valid = ~(np.isnan(vec_a) | np.isnan(vec_b))
    vec_a = vec_a[valid]
    vec_b = vec_b[valid]
    n_pairs = len(vec_a)

    if n_pairs < 3:
        return np.nan, np.nan, 0

    # Observed correlation
    r_obs, _ = stats.pearsonr(vec_a, vec_b)

    # Permutation test: shuffle rows/columns of one matrix
    n_greater = 0
    for _ in range(permutations):
        perm = np.random.permutation(n)
        mat_b_perm = mat_b[np.ix_(perm, perm)]
        vec_b_perm = mat_b_perm[iu][valid]
        r_perm, _ = stats.pearsonr(vec_a, vec_b_perm)
        if r_perm >= r_obs:
            n_greater += 1

    p_val = (n_greater + 1) / (permutations + 1)
    return r_obs, p_val, n_pairs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='STRATA-Atlas: cross-tissue Mantel test (blood r_b vs skin spatial coupling)')
    parser.add_argument('--s5', default='data/raw/science.adt3130_table_s5.xlsx',
                        help='Path to S5 Excel file (regulon AUC)')
    parser.add_argument('--s8', default='data/raw/science.adt3130_table_s8.xlsx',
                        help='Path to S8 Excel file (eQTL r_b)')
    parser.add_argument('--visium-dir', default='data/visium',
                        help='Directory containing Visium .h5 files')
    parser.add_argument('--top-n', type=int, default=5,
                        help='Number of top marker regulons per cell type')
    parser.add_argument('--permutations', type=int, default=9999,
                        help='Number of permutations for Mantel test')
    parser.add_argument('--output', default='results/strata_atlas.csv',
                        help='Output CSV path')
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Load S5 regulon AUC scores
    # ------------------------------------------------------------------
    print(f"Reading {args.s5}, sheet 'eRegulons_Activators_Exp_AUC_RS' ...")
    s5 = pd.read_excel(args.s5, sheet_name='eRegulons_Activators_Exp_AUC_RS')
    print(f"  Columns: {list(s5.columns)}")
    print(f"  Shape: {s5.shape}")
    # Rename columns to standard names used by select_marker_regulons
    s5 = s5.rename(columns={'cell_type_l4': 'cell_type', 'eRegulon': 'regulon'})

    markers = select_marker_regulons(s5, top_n=args.top_n)

    # ------------------------------------------------------------------
    # 2-4. Score Visium spots and build spatial coupling matrix
    # ------------------------------------------------------------------
    print(f"\nScoring Visium spots for {len(markers)} cell types ...")
    spot_scores = score_visium_spots(args.visium_dir, markers)

    print("\nBuilding spatial co-expression coupling matrix ...")
    spatial_mat = build_spatial_coupling_matrix(spot_scores)

    # ------------------------------------------------------------------
    # 5. Load eQTL r_b from S8
    # ------------------------------------------------------------------
    print(f"\nReading {args.s8}, sheet 'cis_eQTL' ...")
    s8 = pd.read_excel(args.s8, sheet_name='cis_eQTL')
    print(f"  Columns: {list(s8.columns)}")
    print(f"  Shape: {s8.shape}")

    rb_matrix = build_rb_matrix(s8)
    print(f"  r_b matrix: {rb_matrix.shape[0]} x {rb_matrix.shape[1]} cell types")

    # ------------------------------------------------------------------
    # 6. Find shared cell types and run Mantel test
    # ------------------------------------------------------------------
    shared_cts = sorted(set(spatial_mat.index) & set(rb_matrix.index))
    print(f"\nShared cell types between spatial and r_b: {len(shared_cts)}")
    if len(shared_cts) < 3:
        print("ERROR: Too few shared cell types for Mantel test. Exiting.")
        return

    # Subset both matrices to shared cell types
    spatial_shared = spatial_mat.loc[shared_cts, shared_cts].values
    rb_shared = rb_matrix.loc[shared_cts, shared_cts].values

    print(f"\nRunning Mantel test ({args.permutations} permutations) ...")
    r_mantel, p_mantel, n_pairs = mantel_test(
        spatial_shared, rb_shared, permutations=args.permutations)

    # ------------------------------------------------------------------
    # 7. Save results and print statistics
    # ------------------------------------------------------------------
    # Build per-cell-type summary
    results_records = []
    for ct in shared_cts:
        # Mean spatial coupling for this cell type
        spatial_vals = spatial_mat.loc[ct].drop(ct, errors='ignore').dropna()
        rb_vals = rb_matrix.loc[ct].drop(ct, errors='ignore').dropna()

        results_records.append({
            'cell_type': ct,
            'mean_spatial_coupling': spatial_vals.mean() if len(spatial_vals) > 0 else np.nan,
            'mean_rb_coupling': rb_vals.mean() if len(rb_vals) > 0 else np.nan,
            'n_spatial_partners': len(spatial_vals),
            'n_rb_partners': len(rb_vals),
            'marker_tfs': ', '.join(markers.get(ct, [])),
        })

    results_df = pd.DataFrame(results_records)
    results_df = results_df.sort_values('mean_spatial_coupling', ascending=False).reset_index(drop=True)

    # Append Mantel test result as metadata row
    mantel_row = pd.DataFrame([{
        'cell_type': '__MANTEL_TEST__',
        'mean_spatial_coupling': r_mantel,
        'mean_rb_coupling': p_mantel,
        'n_spatial_partners': n_pairs,
        'n_rb_partners': len(shared_cts),
        'marker_tfs': f'permutations={args.permutations}',
    }])
    results_df = pd.concat([results_df, mantel_row], ignore_index=True)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(results_df)} rows to {args.output}")

    # Correlation between mean spatial coupling and mean r_b coupling
    valid_df = results_df[results_df['cell_type'] != '__MANTEL_TEST__'].dropna(
        subset=['mean_spatial_coupling', 'mean_rb_coupling'])
    if len(valid_df) >= 3:
        rho_ct, p_ct = stats.spearmanr(
            valid_df['mean_spatial_coupling'], valid_df['mean_rb_coupling'])
    else:
        rho_ct, p_ct = np.nan, np.nan

    # ------------------------------------------------------------------
    # Print key statistics
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STRATA-ATLAS KEY STATISTICS")
    print("=" * 60)
    print(f"  Shared cell types:        {len(shared_cts)}")
    print(f"  Matrix pairs tested:      {n_pairs}")
    print(f"  Mantel r (Pearson):       {r_mantel:.4f}")
    print(f"  Mantel P-value:           {p_mantel:.4e}")
    print(f"  Permutations:             {args.permutations}")
    print(f"")
    print(f"  Per-cell-type Spearman:")
    print(f"    rho (spatial vs r_b):   {rho_ct:.4f}")
    print(f"    P-value:                {p_ct:.2e}")

    print(f"\n  Top 5 spatially coupled cell types:")
    for _, row in valid_df.head(5).iterrows():
        print(f"    {row['cell_type']:30s}  spatial={row['mean_spatial_coupling']:.4f}  "
              f"r_b={row['mean_rb_coupling']:.4f}")

    print(f"\n  Bottom 5 spatially coupled cell types:")
    for _, row in valid_df.tail(5).iterrows():
        print(f"    {row['cell_type']:30s}  spatial={row['mean_spatial_coupling']:.4f}  "
              f"r_b={row['mean_rb_coupling']:.4f}")
    print("=" * 60)


if __name__ == '__main__':
    main()
