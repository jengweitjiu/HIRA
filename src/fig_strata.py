#!/usr/bin/env python3
"""
STRATA Figure — Spatial validation of TOPPLE stability classifications.

4 panels:
  A. Stabilizer vs destabilizer mean expression by condition (grouped bar)
  B. Stabilizer/destabilizer ratio by condition (box + strip)
  C. Spot-level stabilizer expression distributions (lesional vs healthy)
  D. Per-sample scatter: stabilizer vs destabilizer expression, colored by condition

Output: figures/fig_strata.pdf + .png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from pathlib import Path


def main():
    # Style
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
    })

    df = pd.read_csv('results/strata_spatial.csv')

    # Color scheme
    colors = {'lesional': '#D62728', 'non-lesional': '#FF7F0E', 'healthy': '#2CA02C'}
    cond_order = ['healthy', 'non-lesional', 'lesional']

    fig, axes = plt.subplots(2, 2, figsize=(7.08, 6.5))

    # Panel A: Grouped bar — stabilizer vs destabilizer by condition
    ax = axes[0, 0]
    x = np.arange(len(cond_order))
    width = 0.35
    stab_means = [df[df['condition'] == c]['stab_mean_expr'].mean() for c in cond_order]
    dest_means = [df[df['condition'] == c]['dest_mean_expr'].mean() for c in cond_order]
    stab_sems = [df[df['condition'] == c]['stab_mean_expr'].sem() for c in cond_order]
    dest_sems = [df[df['condition'] == c]['dest_mean_expr'].sem() for c in cond_order]

    bars1 = ax.bar(x - width/2, stab_means, width, yerr=stab_sems, capsize=3,
                   color='#4393C3', edgecolor='black', linewidth=0.5, label='Stabilizer TFs')
    bars2 = ax.bar(x + width/2, dest_means, width, yerr=dest_sems, capsize=3,
                   color='#D6604D', edgecolor='black', linewidth=0.5, label='Destabilizer TFs')
    ax.set_xticks(x)
    ax.set_xticklabels(['Healthy', 'Non-lesional', 'Lesional'], fontsize=7)
    ax.set_ylabel('Mean expression')
    ax.legend(fontsize=6, frameon=False)
    ax.set_title('TF expression by condition', fontsize=9, fontweight='bold')
    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold')

    # Add significance bracket for stabilizer les vs healthy
    y_max = max(stab_means) + max(stab_sems) + 0.02
    ax.plot([0 - width/2, 2 - width/2], [y_max, y_max], 'k-', lw=0.8)
    ax.text(1 - width/2, y_max + 0.005, 'P = 7.7e-4', ha='center', fontsize=6)

    # Panel B: Box plot — ratio by condition
    ax = axes[0, 1]
    data_by_cond = [df[df['condition'] == c]['stab_dest_ratio'].values for c in cond_order]
    bp = ax.boxplot(data_by_cond, positions=range(len(cond_order)), widths=0.5,
                    patch_artist=True, showfliers=False)
    for patch, cond in zip(bp['boxes'], cond_order):
        patch.set_facecolor(colors[cond])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
        patch.set_linewidth(0.8)
    for cond, pos in zip(cond_order, range(len(cond_order))):
        vals = df[df['condition'] == cond]['stab_dest_ratio'].values
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(pos + jitter, vals, c=colors[cond], s=20, edgecolors='black',
                   linewidth=0.5, zorder=3, alpha=0.8)
    ax.set_xticks(range(len(cond_order)))
    ax.set_xticklabels(['Healthy', 'Non-lesional', 'Lesional'], fontsize=7)
    ax.set_ylabel('Stabilizer / Destabilizer ratio')
    ax.set_title('Expression ratio by condition', fontsize=9, fontweight='bold')
    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold')

    # Significance bracket
    y_top = df['stab_dest_ratio'].max() + 0.15
    ax.plot([0, 2], [y_top, y_top], 'k-', lw=0.8)
    ax.text(1, y_top + 0.03, 'P = 1.1e-3', ha='center', fontsize=6)

    # Panel C: Sample-level stabilizer expression by condition
    ax = axes[1, 0]
    for cond in cond_order:
        sub = df[df['condition'] == cond]
        vals = sub['stab_mean_expr'].values
        jitter = np.random.default_rng(hash(cond) % 2**32).uniform(-0.12, 0.12, len(vals))
        pos = cond_order.index(cond)
        ax.scatter(pos + jitter, vals, c=colors[cond], s=25, edgecolors='black',
                   linewidth=0.5, alpha=0.8, label=cond.capitalize())
        ax.plot([pos - 0.2, pos + 0.2], [vals.mean(), vals.mean()],
                color='black', lw=2, zorder=4)
    ax.set_xticks(range(len(cond_order)))
    ax.set_xticklabels(['Healthy', 'Non-lesional', 'Lesional'], fontsize=7)
    ax.set_ylabel('Mean stabilizer TF expression')
    ax.set_title('Stabilizer expression per sample', fontsize=9, fontweight='bold')
    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=12, fontweight='bold')

    # Panel D: Scatter — stabilizer vs destabilizer, colored by condition
    ax = axes[1, 1]
    for cond in cond_order:
        sub = df[df['condition'] == cond]
        ax.scatter(sub['dest_mean_expr'], sub['stab_mean_expr'],
                   c=colors[cond], s=30, edgecolors='black', linewidth=0.5,
                   label=cond.capitalize(), alpha=0.8)

    # Fit line
    from scipy import stats
    slope, intercept, r, p, se = stats.linregress(df['dest_mean_expr'], df['stab_mean_expr'])
    x_line = np.linspace(0, df['dest_mean_expr'].max() * 1.05, 100)
    ax.plot(x_line, slope * x_line + intercept, 'k--', lw=0.8, alpha=0.6)
    ax.text(0.05, 0.95, f'r = {r:.3f}\nP = {p:.1e}', transform=ax.transAxes,
            fontsize=6, va='top')
    ax.set_xlabel('Mean destabilizer TF expression')
    ax.set_ylabel('Mean stabilizer TF expression')
    ax.set_title('Stabilizer vs destabilizer', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6, frameon=False, loc='lower right')
    ax.text(-0.15, 1.05, 'D', transform=ax.transAxes, fontsize=12, fontweight='bold')

    plt.tight_layout()
    Path('figures').mkdir(exist_ok=True)
    fig.savefig('figures/fig_strata.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_strata.png', dpi=300, bbox_inches='tight')
    print("Saved figures/fig_strata.pdf and .png")
    plt.close()


if __name__ == '__main__':
    main()
