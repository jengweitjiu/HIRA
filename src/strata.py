#!/usr/bin/env python3
"""
STRATA — Spatial Transcriptomic Regulatory Architecture Testing & Analysis

Validates TOPPLE stability classifications using independent spatial
transcriptomics data from psoriasis Visium (GSE202011, Nat Commun 2022).

Pipeline:
  1. Parse sample metadata from filenames (L=Lesional, NL=Non-Lesional, H=Healthy)
  2. For each .h5 file, extract mean expression of stabilizer/destabilizer TFs
  3. Compare stabilizer vs destabilizer expression ratios across conditions
  4. Test lesional vs healthy differential expression (Mann-Whitney)

Input:  data/visium/*.h5, results/topple.csv
Output: results/strata_spatial.csv, results/strata_summary.csv
"""

import argparse
import re
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

warnings.filterwarnings('ignore', category=UserWarning)


def parse_metadata(visium_dir):
    """Infer sample metadata from filenames."""
    records = []
    for f in sorted(Path(visium_dir).glob('*.h5')):
        name = f.stem
        # Extract GSM ID
        gsm = name.split('_')[0]
        # Extract sample name (everything between GSM and _filtered)
        sample = name.replace(f'{gsm}_', '').replace('_filtered_feature_bc_matrix', '')

        # Determine condition
        if '_HM_' in name or '_HF_' in name:
            condition = 'healthy'
            # Patient ID from e.g. ST_HM_1 or ST_HF_2E
            patient = sample
        elif '_NL' in name:
            condition = 'non-lesional'
            patient = re.search(r'ST_(\d+)', name).group(1)
        elif '_L' in name:
            condition = 'lesional'
            patient = re.search(r'ST_(\d+|L_RP\d)', name).group(1)
        else:
            condition = 'unknown'
            patient = sample

        records.append({
            'gsm': gsm,
            'sample': sample,
            'condition': condition,
            'patient': patient,
            'filepath': str(f),
        })

    meta = pd.DataFrame(records)
    print(f"Sample metadata: {len(meta)} samples")
    print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")
    return meta


def extract_tf_expression(filepath, stabilizer_tfs, destabilizer_tfs):
    """Extract per-spot mean expression of stabilizer and destabilizer TFs."""
    import scanpy as sc
    adata = sc.read_10x_h5(filepath)
    adata.var_names_make_unique()

    # Convert to dense if sparse
    try:
        X = adata.X.toarray()
    except AttributeError:
        X = adata.X

    genes = list(adata.var_names)
    gene_idx = {g: i for i, g in enumerate(genes)}

    def get_mean_expression(tf_list):
        """Mean expression across TFs (per spot), then mean across spots."""
        valid_tfs = [tf for tf in tf_list if tf in gene_idx]
        if not valid_tfs:
            return 0.0, 0.0, 0, []
        cols = [gene_idx[tf] for tf in valid_tfs]
        expr = X[:, cols]  # spots x TFs
        per_spot_mean = expr.mean(axis=1)  # mean across TFs per spot
        return per_spot_mean.mean(), per_spot_mean.std(), len(valid_tfs), per_spot_mean.tolist()

    stab_mean, stab_std, n_stab, stab_spots = get_mean_expression(stabilizer_tfs)
    dest_mean, dest_std, n_dest, dest_spots = get_mean_expression(destabilizer_tfs)

    ratio = stab_mean / dest_mean if dest_mean > 0 else np.nan
    n_spots = adata.n_obs

    return {
        'n_spots': n_spots,
        'n_stab_tfs': n_stab,
        'n_dest_tfs': n_dest,
        'stab_mean_expr': stab_mean,
        'stab_std_expr': stab_std,
        'dest_mean_expr': dest_mean,
        'dest_std_expr': dest_std,
        'stab_dest_ratio': ratio,
        'stab_spot_values': stab_spots,
        'dest_spot_values': dest_spots,
    }


def main():
    parser = argparse.ArgumentParser(description='STRATA: spatial validation of TOPPLE')
    parser.add_argument('--visium-dir', default='data/visium')
    parser.add_argument('--topple', default='results/topple.csv')
    parser.add_argument('--output', default='results/strata_spatial.csv')
    parser.add_argument('--output-summary', default='results/strata_summary.csv')
    args = parser.parse_args()

    # Load TOPPLE classifications
    topple = pd.read_csv(args.topple)
    topple['tf_name'] = topple['regulon'].str.replace(r'_[+-]$', '', regex=True)
    stabilizers = topple[topple['stability_class'] == 'stabilizer']['tf_name'].tolist()
    destabilizers = topple[topple['stability_class'] == 'destabilizer']['tf_name'].tolist()
    print(f"TOPPLE: {len(stabilizers)} stabilizer TFs, {len(destabilizers)} destabilizer TFs")
    print(f"  Top stabilizers: {stabilizers[:10]}")
    print(f"  Top destabilizers: {destabilizers[:10]}")

    # Parse metadata
    meta = parse_metadata(args.visium_dir)

    # Extract TF expression from each sample
    print(f"\nExtracting TF expression from {len(meta)} samples ...")
    results = []
    for _, row in meta.iterrows():
        print(f"  {row['sample']} ({row['condition']}) ... ", end='')
        expr = extract_tf_expression(row['filepath'], stabilizers, destabilizers)
        print(f"{expr['n_spots']} spots, stab={expr['stab_mean_expr']:.4f}, "
              f"dest={expr['dest_mean_expr']:.4f}, ratio={expr['stab_dest_ratio']:.3f}")
        record = {**row.to_dict(), **{k: v for k, v in expr.items()
                                       if k not in ('stab_spot_values', 'dest_spot_values')}}
        record['_stab_spots'] = expr['stab_spot_values']
        record['_dest_spots'] = expr['dest_spot_values']
        results.append(record)

    df = pd.DataFrame(results)

    # Save per-sample results (without spot-level data)
    save_cols = [c for c in df.columns if not c.startswith('_')]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df[save_cols].to_csv(args.output, index=False)
    print(f"\nSaved per-sample results to {args.output}")

    # ============================================================
    # Statistical tests
    # ============================================================
    print("\n" + "=" * 60)
    print("STRATA KEY STATISTICS")
    print("=" * 60)

    # 1. Stabilizer vs destabilizer expression (all spots pooled)
    all_stab = []
    all_dest = []
    for r in results:
        all_stab.extend(r['_stab_spots'])
        all_dest.extend(r['_dest_spots'])

    u_all, p_all = stats.mannwhitneyu(all_stab, all_dest, alternative='two-sided')
    print(f"\n  1. Stabilizer vs Destabilizer (all {len(all_stab)} spots):")
    print(f"     Stabilizer mean: {np.mean(all_stab):.4f}")
    print(f"     Destabilizer mean: {np.mean(all_dest):.4f}")
    print(f"     MW P = {p_all:.2e}")

    # 2. Lesional vs Healthy comparison
    les = df[df['condition'] == 'lesional']
    heal = df[df['condition'] == 'healthy']
    nl = df[df['condition'] == 'non-lesional']

    print(f"\n  2. Per-sample stabilizer/destabilizer ratio:")
    for cond, sub in [('lesional', les), ('non-lesional', nl), ('healthy', heal)]:
        print(f"     {cond:15s}: mean ratio={sub['stab_dest_ratio'].mean():.3f} "
              f"(n={len(sub)})")

    # Lesional vs Healthy: stabilizer expression
    if len(les) > 0 and len(heal) > 0:
        u_lh_stab, p_lh_stab = stats.mannwhitneyu(
            les['stab_mean_expr'], heal['stab_mean_expr'], alternative='two-sided')
        u_lh_dest, p_lh_dest = stats.mannwhitneyu(
            les['dest_mean_expr'], heal['dest_mean_expr'], alternative='two-sided')
        u_lh_ratio, p_lh_ratio = stats.mannwhitneyu(
            les['stab_dest_ratio'], heal['stab_dest_ratio'], alternative='two-sided')

        print(f"\n  3. Lesional vs Healthy (Mann-Whitney):")
        print(f"     Stabilizer expr:  les={les['stab_mean_expr'].mean():.4f} vs "
              f"heal={heal['stab_mean_expr'].mean():.4f}, P={p_lh_stab:.2e}")
        print(f"     Destabilizer expr: les={les['dest_mean_expr'].mean():.4f} vs "
              f"heal={heal['dest_mean_expr'].mean():.4f}, P={p_lh_dest:.2e}")
        print(f"     Ratio:            les={les['stab_dest_ratio'].mean():.3f} vs "
              f"heal={heal['stab_dest_ratio'].mean():.3f}, P={p_lh_ratio:.2e}")

    # 3. Spot-level: lesional vs healthy for stabilizer expression
    les_stab_spots = []
    heal_stab_spots = []
    les_dest_spots = []
    heal_dest_spots = []
    for r in results:
        if r['condition'] == 'lesional':
            les_stab_spots.extend(r['_stab_spots'])
            les_dest_spots.extend(r['_dest_spots'])
        elif r['condition'] == 'healthy':
            heal_stab_spots.extend(r['_stab_spots'])
            heal_dest_spots.extend(r['_dest_spots'])

    if les_stab_spots and heal_stab_spots:
        u_spot, p_spot = stats.mannwhitneyu(
            les_stab_spots, heal_stab_spots, alternative='two-sided')
        print(f"\n  4. Spot-level lesional vs healthy (stabilizer expression):")
        print(f"     Lesional: {len(les_stab_spots)} spots, mean={np.mean(les_stab_spots):.4f}")
        print(f"     Healthy:  {len(heal_stab_spots)} spots, mean={np.mean(heal_stab_spots):.4f}")
        print(f"     MW P = {p_spot:.2e}")

        u_spot_d, p_spot_d = stats.mannwhitneyu(
            les_dest_spots, heal_dest_spots, alternative='two-sided')
        print(f"\n  5. Spot-level lesional vs healthy (destabilizer expression):")
        print(f"     Lesional: {len(les_dest_spots)} spots, mean={np.mean(les_dest_spots):.4f}")
        print(f"     Healthy:  {len(heal_dest_spots)} spots, mean={np.mean(heal_dest_spots):.4f}")
        print(f"     MW P = {p_spot_d:.2e}")

    # 4. Non-lesional vs healthy
    if len(nl) > 0 and len(heal) > 0:
        nl_stab_spots = []
        heal_stab_spots2 = []
        for r in results:
            if r['condition'] == 'non-lesional':
                nl_stab_spots.extend(r['_stab_spots'])
            elif r['condition'] == 'healthy':
                heal_stab_spots2.extend(r['_stab_spots'])
        if nl_stab_spots and heal_stab_spots2:
            u_nl, p_nl = stats.mannwhitneyu(
                nl_stab_spots, heal_stab_spots2, alternative='two-sided')
            print(f"\n  6. Spot-level non-lesional vs healthy (stabilizer expression):")
            print(f"     Non-lesional: {len(nl_stab_spots)} spots, mean={np.mean(nl_stab_spots):.4f}")
            print(f"     Healthy:      {len(heal_stab_spots2)} spots, mean={np.mean(heal_stab_spots2):.4f}")
            print(f"     MW P = {p_nl:.2e}")

    # Summary table
    summary_records = []
    for cond in ['lesional', 'non-lesional', 'healthy']:
        sub = df[df['condition'] == cond]
        if len(sub) == 0:
            continue
        summary_records.append({
            'condition': cond,
            'n_samples': len(sub),
            'total_spots': sub['n_spots'].sum(),
            'mean_stab_expr': sub['stab_mean_expr'].mean(),
            'mean_dest_expr': sub['dest_mean_expr'].mean(),
            'mean_ratio': sub['stab_dest_ratio'].mean(),
            'std_ratio': sub['stab_dest_ratio'].std(),
        })
    summary = pd.DataFrame(summary_records)
    summary.to_csv(args.output_summary, index=False)
    print(f"\nSaved summary to {args.output_summary}")
    print("=" * 60)


if __name__ == '__main__':
    main()
