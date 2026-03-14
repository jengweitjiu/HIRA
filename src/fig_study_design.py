#!/usr/bin/env python3
"""
DECODE-AD Study Design Overview Figure (Fig 1) and Circos-style summary (Fig 6).

Fig 1: Flowchart — AD GWAS -> CIMA -> HIRA layers -> key findings
Fig 6: Summary radial plot — AD loci x cell types x regulatory layers
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
from pathlib import Path

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8})


def fig_study_design():
    """Create study design flowchart."""
    fig, ax = plt.subplots(1, 1, figsize=(180/25.4, 220/25.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis('off')

    # Color scheme
    c_input = '#3498DB'
    c_hira = '#2ECC71'
    c_decode = '#E74C3C'
    c_finding = '#9B59B6'

    def box(x, y, w, h, text, color, fontsize=7):
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                               facecolor=color, edgecolor='black', linewidth=0.8, alpha=0.85)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='white',
                wrap=True)

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#555', lw=1.2))

    # Title
    ax.text(5, 11.6, 'DECODE-AD Study Design', ha='center', fontsize=14,
            fontweight='bold')

    # Row 1: Input data
    box(0.5, 10.2, 3.5, 0.9, 'AD GWAS\n(Budu-Aggrey 2023)\n25 lead loci', c_input, 7)
    box(6, 10.2, 3.5, 0.9, 'CIMA Atlas\n(Yin et al. 2026)\n69 cell types', c_input, 7)

    # Arrows to HIRA
    arrow(2.25, 10.2, 2.25, 9.5)
    arrow(7.75, 10.2, 7.75, 9.5)

    # Row 2: HIRA Framework
    box(0.2, 8.6, 2.0, 0.8, 'TOPPLE\n203 regulons', c_hira, 6)
    box(2.5, 8.6, 1.8, 0.8, 'DGSA\n5,253 genes', c_hira, 6)
    box(4.6, 8.6, 1.8, 0.8, 'SICAI\n69x69 r_b', c_hira, 6)
    box(6.7, 8.6, 1.5, 0.8, 'IPA\nSex/Age', c_hira, 6)
    box(8.5, 8.6, 1.3, 0.8, 'STRATA\nVisium', c_hira, 6)
    ax.text(5, 9.55, 'HIRA Six-Layer Framework', ha='center', fontsize=9,
            fontweight='bold', color=c_hira)

    # Row 3: DECODE-AD analyses
    arrow(5, 8.6, 5, 7.8)
    ax.text(5, 7.9, 'DECODE-AD Integration', ha='center', fontsize=9,
            fontweight='bold', color=c_decode)

    box(0.2, 6.5, 2.2, 1.0, 'xQTL Mapping\n23 eQTL loci\n24 caQTL loci\n195 eGenes', c_decode, 6)
    box(2.7, 6.5, 2.2, 1.0, 'Cell-type PRS\nT cells: 68.7%\nMyeloid: 11.3%\nB cells: 10.3%', c_decode, 6)
    box(5.2, 6.5, 2.2, 1.0, 'Regulatory\nCascades\n5 AD TF regulons\nETS1, IRF1, RELA', c_decode, 6)
    box(7.7, 6.5, 2.2, 1.0, 'Coupling\nDisruption\n79% edges weaker\nCd4_Tr1 most hit', c_decode, 6)

    # Row 4: Extended analyses
    arrow(5, 6.5, 5, 5.7)

    box(0.2, 4.5, 2.8, 1.0, 'Population Genetics\n78.7% AF divergent\n5 selection overlaps\n(all MHC region)', c_decode, 6)
    box(3.3, 4.5, 2.0, 1.0, 'Drug Targets\nJAK3: 27 CTs\nIL4R: 18 CTs\nIFNG: AD eGene', c_decode, 6)
    box(5.6, 4.5, 2.2, 1.0, 'East Asian\nReplication\n17/17 loci overlap\n31 shared genes', c_decode, 6)
    box(8.1, 4.5, 1.8, 1.0, 'Interaction\nQTLs (S12)\n60 hits\nMono + B cells', c_decode, 6)

    # Row 5: Key findings
    arrow(5, 4.5, 5, 3.7)
    ax.text(5, 3.8, 'Key Findings', ha='center', fontsize=9,
            fontweight='bold', color=c_finding)

    findings = [
        'T cells carry 68.7% of AD genetic risk weight',
        'cMono_CD14 dominates caQTL landscape (712 peaks)',
        'IRF1 cascade: 998 target genes from AD cis-eQTL',
        'CD4_Tr1_IL10 coupling most disrupted in AD context',
        'Japanese GWAS shows 100% CIMA overlap (vs 92% EUR)',
    ]
    for i, f in enumerate(findings):
        y = 3.0 - i * 0.5
        ax.text(0.5, y, f'{i+1}. {f}', fontsize=6.5, va='center',
                color='#333', style='italic')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=c_input, label='Input Data'),
        mpatches.Patch(facecolor=c_hira, label='HIRA Framework'),
        mpatches.Patch(facecolor=c_decode, label='DECODE-AD Analysis'),
        mpatches.Patch(facecolor=c_finding, label='Key Findings'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=6, frameon=True)

    plt.tight_layout()
    Path("figures").mkdir(exist_ok=True)
    fig.savefig('figures/fig_study_design.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_study_design.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_study_design.pdf + .png")


def fig_circos_summary():
    """Create a radial summary plot (Circos-style) of AD loci x cell types."""
    fig, ax = plt.subplots(1, 1, figsize=(180/25.4, 180/25.4), subplot_kw={'projection': 'polar'})

    # Load data
    eqtl_matrix = pd.read_csv('results/decode_ad_locus_celltype_eqtl.csv', index_col=0)
    caqtl_matrix = pd.read_csv('results/decode_ad_locus_celltype_caqtl.csv', index_col=0)

    # Sum across loci per cell type
    eqtl_ct = eqtl_matrix.sum(axis=0).sort_values(ascending=False)
    caqtl_ct = caqtl_matrix.sum(axis=0).sort_values(ascending=False)

    # Top 20 cell types
    top_cts = eqtl_ct.head(20).index.tolist()

    n = len(top_cts)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    width = 2 * np.pi / n * 0.8

    # eQTL bars (inner)
    eqtl_vals = [eqtl_ct.get(ct, 0) for ct in top_cts]
    bars1 = ax.bar(theta, eqtl_vals, width=width, bottom=0,
                    color='#3498DB', alpha=0.7, label='eQTL')

    # caQTL bars (stacked)
    caqtl_vals = [caqtl_ct.get(ct, 0) for ct in top_cts]
    bars2 = ax.bar(theta, caqtl_vals, width=width, bottom=eqtl_vals,
                    color='#E74C3C', alpha=0.7, label='caQTL')

    # Labels
    ax.set_xticks(theta)
    labels = [ct.split('_')[0] + '_' + '_'.join(ct.split('_')[1:2])
              if len(ct.split('_')) > 1 else ct for ct in top_cts]
    ax.set_xticklabels(labels, fontsize=5, rotation=0)

    # Rotate labels to match angle
    for label, angle in zip(ax.get_xticklabels(), theta):
        if angle > np.pi/2 and angle < 3*np.pi/2:
            label.set_rotation(np.degrees(angle) + 180)
            label.set_ha('right')
        else:
            label.set_rotation(np.degrees(angle))
            label.set_ha('left')

    ax.set_title('DECODE-AD: Cell-Type xQTL Landscape\n(Top 20 Cell Types)',
                  fontsize=10, fontweight='bold', pad=20)
    ax.legend(loc='lower right', fontsize=7, bbox_to_anchor=(1.15, -0.05))

    plt.tight_layout()
    fig.savefig('figures/fig_circos_summary.pdf', dpi=300, bbox_inches='tight')
    fig.savefig('figures/fig_circos_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved figures/fig_circos_summary.pdf + .png")


if __name__ == "__main__":
    fig_study_design()
    fig_circos_summary()
