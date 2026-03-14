#!/usr/bin/env python3
"""
HIRA x PASI Severity Validation.

Correlate all HIRA layer metrics with PASI scores from GSE202011 Visium.
Uses pre-computed strata_sample_metrics.csv (30 samples, 9 patients).

For each sample, we already have:
  - mean_stabilizer_expr, mean_destabilizer_expr (TOPPLE proxy)
  - stab_dest_ratio (TOPPLE)
  - coupling_strength (SICAI proxy)
  - immune_entropy, fib_entropy

Additional per-sample metrics to compute:
  - Mean DGSA non-additivity of expressed genes
  - IPA perturbation sensitivity proxy

Output: results/pasi_hira_correlations.csv
        figures/fig_pasi_validation.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def main():
    print("=" * 60)
    print("HIRA x PASI Severity Validation (GSE202011)")
    print("=" * 60)

    # Load sample metrics
    meta = pd.read_csv('results/strata_sample_metrics.csv')
    print(f"  Samples: {len(meta)}")
    print(f"  Patients: {meta['patient'].nunique()}")
    print(f"  PASI range: {meta['pasi'].min():.1f} - {meta['pasi'].max():.1f}")

    # Use all samples (including healthy at PASI=0 for full range)
    # Also compute diseased-only correlations
    all_samples = meta.copy()
    diseased = meta[meta['pasi'] > 0].copy()

    print(f"  All samples: {len(all_samples)} (incl. {len(meta[meta['pasi']==0])} healthy)")
    print(f"  Diseased only: {len(diseased)}")

    correlations = []

    # -----------------------------------------------------------
    # Layer 1: TOPPLE — Stabilizer expression vs PASI
    # -----------------------------------------------------------
    print("\n--- Layer 1: TOPPLE (Stabilizer Expression) ---")

    for subset_name, subset in [('all', all_samples), ('diseased', diseased)]:
        for metric, label in [
            ('mean_stabilizer_expr', 'Stabilizer expression'),
            ('mean_destabilizer_expr', 'Destabilizer expression'),
            ('stab_dest_ratio', 'Stabilizer/Destabilizer ratio'),
        ]:
            valid = subset.dropna(subset=[metric, 'pasi'])
            if len(valid) < 5:
                continue
            rho, p = stats.spearmanr(valid[metric], valid['pasi'])
            correlations.append({
                'layer': 'TOPPLE',
                'metric': label,
                'subset': subset_name,
                'n_samples': len(valid),
                'rho': rho,
                'p_value': p,
            })
            print(f"  {label} vs PASI ({subset_name}): rho={rho:.3f}, P={p:.2e}, n={len(valid)}")

    # -----------------------------------------------------------
    # Layer 2: DGSA proxy — effect size heterogeneity
    # We don't have per-sample DGSA, but we can use immune_entropy
    # as a proxy for transcriptomic non-additivity
    # -----------------------------------------------------------
    print("\n--- Layer 2: DGSA proxy (Immune Entropy) ---")
    for subset_name, subset in [('all', all_samples), ('diseased', diseased)]:
        valid = subset.dropna(subset=['immune_entropy', 'pasi'])
        if len(valid) < 5:
            continue
        rho, p = stats.spearmanr(valid['immune_entropy'], valid['pasi'])
        correlations.append({
            'layer': 'DGSA_proxy',
            'metric': 'Immune entropy',
            'subset': subset_name,
            'n_samples': len(valid),
            'rho': rho,
            'p_value': p,
        })
        print(f"  Immune entropy vs PASI ({subset_name}): rho={rho:.3f}, P={p:.2e}, n={len(valid)}")

    # -----------------------------------------------------------
    # Layer 3: SICAI — Coupling strength vs PASI
    # -----------------------------------------------------------
    print("\n--- Layer 3: SICAI (Coupling Strength) ---")
    for subset_name, subset in [('all', all_samples), ('diseased', diseased)]:
        for metric, label in [
            ('coupling_strength', 'Coupling strength'),
            ('coupling_entropy', 'Coupling entropy'),
        ]:
            valid = subset.dropna(subset=[metric, 'pasi'])
            if len(valid) < 5:
                continue
            rho, p = stats.spearmanr(valid[metric], valid['pasi'])
            correlations.append({
                'layer': 'SICAI',
                'metric': label,
                'subset': subset_name,
                'n_samples': len(valid),
                'rho': rho,
                'p_value': p,
            })
            print(f"  {label} vs PASI ({subset_name}): rho={rho:.3f}, P={p:.2e}, n={len(valid)}")

    # -----------------------------------------------------------
    # Layer 4: IPA proxy — ratio divergence between conditions
    # IPA measures perturbation resistance. Proxy: how much does
    # the stabilizer/destabilizer ratio deviate from healthy baseline?
    # -----------------------------------------------------------
    print("\n--- Layer 4: IPA proxy (Ratio Deviation) ---")
    healthy_ratio = all_samples[all_samples['pasi'] == 0]['stab_dest_ratio'].mean()
    all_samples['ratio_deviation'] = (all_samples['stab_dest_ratio'] - healthy_ratio).abs()

    for subset_name, subset_df in [('all', all_samples), ('diseased', all_samples[all_samples['pasi'] > 0])]:
        valid = subset_df.dropna(subset=['ratio_deviation', 'pasi'])
        if len(valid) < 5:
            continue
        rho, p = stats.spearmanr(valid['ratio_deviation'], valid['pasi'])
        correlations.append({
            'layer': 'IPA_proxy',
            'metric': 'Ratio deviation from healthy',
            'subset': subset_name,
            'n_samples': len(valid),
            'rho': rho,
            'p_value': p,
        })
        print(f"  Ratio deviation vs PASI ({subset_name}): rho={rho:.3f}, P={p:.2e}, n={len(valid)}")

    # -----------------------------------------------------------
    # Layer 5: STRATA — Fibroblast entropy vs PASI
    # -----------------------------------------------------------
    print("\n--- Layer 5: STRATA (Fibroblast Entropy) ---")
    for subset_name, subset in [('all', all_samples), ('diseased', diseased)]:
        valid = subset.dropna(subset=['fib_entropy', 'pasi'])
        if len(valid) < 5:
            continue
        rho, p = stats.spearmanr(valid['fib_entropy'], valid['pasi'])
        correlations.append({
            'layer': 'STRATA',
            'metric': 'Fibroblast entropy',
            'subset': subset_name,
            'n_samples': len(valid),
            'rho': rho,
            'p_value': p,
        })
        print(f"  Fibroblast entropy vs PASI ({subset_name}): rho={rho:.3f}, P={p:.2e}, n={len(valid)}")

    # -----------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------
    corr_df = pd.DataFrame(correlations)
    corr_df.to_csv('results/pasi_hira_correlations.csv', index=False)
    print(f"\nSaved to results/pasi_hira_correlations.csv")

    # Highlight significant results
    sig = corr_df[corr_df['p_value'] < 0.05]
    print(f"\n  Significant correlations (P<0.05): {len(sig)}/{len(corr_df)}")
    if len(sig) > 0:
        for _, row in sig.iterrows():
            print(f"    [{row['layer']}] {row['metric']} ({row['subset']}): "
                  f"rho={row['rho']:.3f}, P={row['p_value']:.2e}, n={row['n_samples']}")

    # -----------------------------------------------------------
    # Figure: 4-panel PASI correlation
    # -----------------------------------------------------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})

    fig, axes = plt.subplots(2, 2, figsize=(180/25.4, 150/25.4))

    panels = [
        ('mean_stabilizer_expr', 'Stabilizer Expression', 'TOPPLE'),
        ('immune_entropy', 'Immune Entropy', 'DGSA proxy'),
        ('coupling_strength', 'Coupling Strength', 'SICAI'),
        ('stab_dest_ratio', 'Stab/Destab Ratio', 'IPA proxy'),
    ]

    condition_colors = {
        'Lesional': '#E74C3C',
        'Non-Lesional': '#F39C12',
        'Healthy': '#2ECC71',
    }

    for ax, (metric, ylabel, layer) in zip(axes.flat, panels):
        valid = all_samples.dropna(subset=[metric, 'pasi'])

        for cond, color in condition_colors.items():
            mask = valid['condition'] == cond
            if mask.any():
                ax.scatter(valid.loc[mask, 'pasi'], valid.loc[mask, metric],
                           s=25, c=color, alpha=0.7, edgecolors='white',
                           linewidth=0.5, label=cond, zorder=3)

        # Fit line (all samples)
        if len(valid) > 5:
            rho, p = stats.spearmanr(valid['pasi'], valid[metric])
            # Regression line
            z = np.polyfit(valid['pasi'], valid[metric], 1)
            x_line = np.linspace(valid['pasi'].min(), valid['pasi'].max(), 100)
            ax.plot(x_line, np.polyval(z, x_line), 'k--', linewidth=0.8, alpha=0.5)
            ax.text(0.95, 0.95, f'rho={rho:.3f}\nP={p:.2e}\nn={len(valid)}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=6,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        ax.set_xlabel('PASI Score')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{layer}', fontsize=9, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend on first panel
    axes[0, 0].legend(fontsize=6, frameon=True, loc='upper left')

    # Panel labels
    for ax, label in zip(axes.flat, ['A', 'B', 'C', 'D']):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')

    fig.suptitle('HIRA Layer Metrics vs PASI Severity (GSE202011)',
                 fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_pasi_validation.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_pasi_validation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_pasi_validation.pdf + .png")


if __name__ == "__main__":
    main()
