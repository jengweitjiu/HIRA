#!/usr/bin/env python3
"""
Publication figure: Extension #1 — Epigenomic amplification of non-additivity.

Four panels:
  (A) Overlaid density histograms: eQTL vs caQTL non-additivity
  (B) QQ-plot of the two distributions with crossover annotation
  (C) Paired scatter: matched-locus eQTL vs caQTL non-additivity
  (D) Bar chart: Spearman rho with disease count (eQTL vs caQTL)

Input:  results/dgsa_eqtl.csv, results/ext1_caqtl_dgsa.csv,
        data/raw/CIMA_Table_S6.csv, data/raw/science.adt3130_table_s15.xlsx
Output: figures/fig_ext1_epigenomic_amplification.pdf + .png
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


# ---------- helpers ----------------------------------------------------------

def compute_gini(arr):
    arr = np.sort(arr.copy())
    n = len(arr)
    if n == 0 or arr.sum() == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * arr) - (n + 1) * np.sum(arr)) / (n * np.sum(arr))


def setup_style():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 7,
        'axes.linewidth': 0.6,
        'axes.labelsize': 7.5,
        'axes.titlesize': 9,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'lines.linewidth': 0.8,
        'savefig.dpi': 300,
    })


BLUE = '#4878CF'
RED = '#D65F5F'
GREY = '#555555'


# ---------- data loading -----------------------------------------------------

def load_matched_pairs(s6_path):
    """Return (paired_eqtl_na, paired_caqtl_na) for matched variant_id+celltype."""
    df = pd.read_csv(s6_path, usecols=[
        'phenotype_id', 'celltype', 'slope', 'analysis', 'variant_id',
    ])
    eqtl = df[df['analysis'] == 'cis-eQTL']
    caqtl = df[df['analysis'] == 'cis-caQTL']

    matched = eqtl.merge(caqtl, on=['variant_id', 'celltype'],
                         suffixes=('_eqtl', '_caqtl'))
    gene_peak_pairs = matched[
        ['phenotype_id_eqtl', 'phenotype_id_caqtl']
    ].drop_duplicates()

    # Full-space non-additivity for every gene and peak
    eqtl_cts = sorted(eqtl['celltype'].unique())
    caqtl_cts = sorted(caqtl['celltype'].unique())

    eqtl_pivot = eqtl.pivot_table(
        index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
    ).reindex(columns=eqtl_cts).fillna(0)

    caqtl_pivot = caqtl.pivot_table(
        index='phenotype_id', columns='celltype', values='slope', aggfunc='first'
    ).reindex(columns=caqtl_cts).fillna(0)

    gene_na = {g: compute_gini(eqtl_pivot.loc[g].values ** 2) for g in eqtl_pivot.index}
    peak_na = {p: compute_gini(caqtl_pivot.loc[p].values ** 2) for p in caqtl_pivot.index}

    pe, pc = [], []
    for _, row in gene_peak_pairs.iterrows():
        g, p = row['phenotype_id_eqtl'], row['phenotype_id_caqtl']
        if g in gene_na and p in peak_na:
            pe.append(gene_na[g])
            pc.append(peak_na[p])
    return np.array(pe), np.array(pc)


def load_disease_rhos(s6_path, s15_path, eqtl_df, caqtl_df):
    """Return dict with rho and pval for eQTL and caQTL."""
    s15 = pd.read_excel(s15_path, sheet_name='Sheet1')
    disease_ct = s15.groupby('celltype')['trait'].nunique().reset_index()
    disease_ct.columns = ['celltype', 'disease_count']

    df = pd.read_csv(s6_path, usecols=['phenotype_id', 'celltype', 'slope', 'analysis'])

    out = {}
    for label, analysis_val, dgsa_df, id_col in [
        ('eQTL', 'cis-eQTL', eqtl_df, 'gene'),
        ('caQTL', 'cis-caQTL', caqtl_df, 'peak'),
    ]:
        sub = df[df['analysis'] == analysis_val]
        pivot = sub.pivot_table(
            index='phenotype_id', columns='celltype', values='slope', aggfunc='first')
        pivot = pivot[pivot.notna().sum(axis=1) >= 3]
        na_lookup = dict(zip(dgsa_df[id_col], dgsa_df['non_additivity']))

        ct_records = []
        for ct in pivot.columns:
            active = pivot[ct].dropna().index
            vals = [na_lookup[g] for g in active if g in na_lookup]
            if vals:
                ct_records.append({'celltype': ct, 'mean_na': np.mean(vals)})
        ct_df = pd.DataFrame(ct_records)
        merged = ct_df.merge(disease_ct, on='celltype', how='inner')
        rho, pval = stats.spearmanr(merged['mean_na'], merged['disease_count'])
        out[label] = {'rho': rho, 'pval': pval, 'n': len(merged)}
    return out


# ---------- figure -----------------------------------------------------------

def make_figure(eqtl_na, caqtl_na, paired_e, paired_c, disease_rhos, outpath):
    setup_style()
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 150 / 25.4))

    # ---- Panel A: Overlaid histograms ----
    ax = axes[0, 0]
    bins = np.linspace(0.3, 1.0, 60)
    ax.hist(eqtl_na, bins=bins, alpha=0.55, density=True, color=BLUE,
            label=f'eQTL (n = {len(eqtl_na):,})', edgecolor='none')
    ax.hist(caqtl_na, bins=bins, alpha=0.55, density=True, color=RED,
            label=f'caQTL (n = {len(caqtl_na):,})', edgecolor='none')
    # Mean lines
    ax.axvline(np.mean(eqtl_na), color=BLUE, ls='--', lw=0.7,
               label=f'eQTL mean = {np.mean(eqtl_na):.3f}')
    ax.axvline(np.mean(caqtl_na), color=RED, ls='--', lw=0.7,
               label=f'caQTL mean = {np.mean(caqtl_na):.3f}')
    ax.set_xlabel('Non-additivity')
    ax.set_ylabel('Density')
    ax.legend(frameon=False, loc='upper left')
    ks_d, ks_p = stats.ks_2samp(eqtl_na, caqtl_na)
    ax.text(0.97, 0.95, f'KS D = {ks_d:.3f}\nP = {ks_p:.1e}',
            transform=ax.transAxes, ha='right', va='top', fontsize=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
    ax.set_title('A', fontweight='bold', loc='left', pad=4)

    # ---- Panel B: QQ-plot ----
    ax = axes[0, 1]
    probs = np.linspace(0.5, 99.5, 500)
    eq_q = np.percentile(eqtl_na, probs)
    cq_q = np.percentile(caqtl_na, probs)
    ax.scatter(eq_q, cq_q, s=3, alpha=0.5, color=GREY, zorder=3, linewidths=0)
    lims = [min(eq_q.min(), cq_q.min()) - 0.02, max(eq_q.max(), cq_q.max()) + 0.02]
    ax.plot(lims, lims, 'k--', lw=0.5, alpha=0.4, zorder=2)
    ax.set_xlabel('eQTL quantiles')
    ax.set_ylabel('caQTL quantiles')
    ax.set_aspect('equal')
    # Shade regions
    ax.fill_between(lims, lims, [lims[1], lims[1]], alpha=0.03, color=RED,
                    label='caQTL > eQTL')
    ax.fill_between(lims, [lims[0], lims[0]], lims, alpha=0.03, color=BLUE,
                    label='eQTL > caQTL')
    # Crossover annotation
    cross_idx = np.where(np.diff(np.sign(cq_q - eq_q)))[0]
    if len(cross_idx) > 0:
        cx = cross_idx[0]
        ax.plot(eq_q[cx], cq_q[cx], 'o', ms=4, mfc='none', mec='red', mew=0.8, zorder=5)
        ax.annotate(f'crossover\n~{eq_q[cx]:.2f}',
                    xy=(eq_q[cx], cq_q[cx]),
                    xytext=(eq_q[cx] - 0.12, cq_q[cx] + 0.02),
                    fontsize=6, color='red',
                    arrowprops=dict(arrowstyle='->', color='red', lw=0.5))
    ax.legend(frameon=False, loc='lower right', fontsize=5.5)
    ax.set_title('B', fontweight='bold', loc='left', pad=4)

    # ---- Panel C: Paired scatter ----
    ax = axes[1, 0]
    ax.scatter(paired_e, paired_c, s=2, alpha=0.15, color=GREY, linewidths=0,
               rasterized=True)
    lims_c = [0.3, 1.02]
    ax.plot(lims_c, lims_c, 'k--', lw=0.5, alpha=0.4)
    ax.set_xlabel('eQTL non-additivity (matched gene)')
    ax.set_ylabel('caQTL non-additivity (matched peak)')
    ax.set_xlim(lims_c)
    ax.set_ylim(lims_c)
    ax.set_aspect('equal')
    # Stats
    diff = paired_c - paired_e
    w_stat, w_p = stats.wilcoxon(diff, alternative='greater')
    n_above = np.sum(diff > 0)
    pct_above = 100 * n_above / len(diff)
    stats_text = (f'n = {len(paired_e):,} pairs\n'
                  f'caQTL > eQTL: {pct_above:.1f}%\n'
                  f'Wilcoxon P = {w_p:.1e}')
    ax.text(0.03, 0.97, stats_text, transform=ax.transAxes,
            ha='left', va='top', fontsize=6,
            bbox=dict(facecolor='white', edgecolor='#cccccc', alpha=0.9, lw=0.4))
    ax.set_title('C', fontweight='bold', loc='left', pad=4)

    # ---- Panel D: Bar chart — rho with disease ----
    ax = axes[1, 1]
    labels = ['eQTL', 'caQTL']
    rhos = [disease_rhos['eQTL']['rho'], disease_rhos['caQTL']['rho']]
    pvals = [disease_rhos['eQTL']['pval'], disease_rhos['caQTL']['pval']]
    ns = [disease_rhos['eQTL']['n'], disease_rhos['caQTL']['n']]
    colors = [BLUE, RED]

    bars = ax.bar(labels, rhos, width=0.5, color=colors, alpha=0.75, edgecolor='none')
    ax.set_ylabel('Spearman $\\rho$ with disease count')
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_locator(MultipleLocator(0.2))

    for i, (bar, rho, pv, n) in enumerate(zip(bars, rhos, pvals, ns)):
        y = bar.get_height()
        sig = _pval_stars(pv)
        ax.text(bar.get_x() + bar.get_width() / 2, y + 0.02,
                f'{rho:.3f}{sig}\n(n = {n})',
                ha='center', va='bottom', fontsize=6)
    ax.set_title('D', fontweight='bold', loc='left', pad=4)

    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    fig.savefig(str(outpath).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close(fig)


def _pval_stars(p):
    if p < 1e-10:
        return '***'
    elif p < 1e-4:
        return '**'
    elif p < 0.05:
        return '*'
    return ' n.s.'


# ---------- main -------------------------------------------------------------

def main():
    s6_path = 'data/raw/CIMA_Table_S6.csv'
    s15_path = 'data/raw/science.adt3130_table_s15.xlsx'
    eqtl_path = 'results/dgsa_eqtl.csv'
    caqtl_path = 'results/ext1_caqtl_dgsa.csv'
    outpath = 'figures/fig_ext1_epigenomic_amplification.pdf'

    print("Loading data ...")
    eqtl_df = pd.read_csv(eqtl_path)
    caqtl_df = pd.read_csv(caqtl_path)

    print("Computing matched pairs ...")
    paired_e, paired_c = load_matched_pairs(s6_path)
    print(f"  {len(paired_e):,} matched gene-peak pairs")

    print("Computing disease correlations ...")
    disease_rhos = load_disease_rhos(s6_path, s15_path, eqtl_df, caqtl_df)
    for k, v in disease_rhos.items():
        print(f"  {k}: rho={v['rho']:.4f}, P={v['pval']:.2e}, n={v['n']}")

    print("Generating figure ...")
    make_figure(
        eqtl_df['non_additivity'].values,
        caqtl_df['non_additivity'].values,
        paired_e, paired_c,
        disease_rhos,
        outpath,
    )
    print("Done.")


if __name__ == '__main__':
    main()
