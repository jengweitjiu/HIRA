#!/usr/bin/env python3
"""
Figure 1 — Complete HIRA 8-panel overview.

(A) TOPPLE: stability ranking bar (top 15 regulons)
(B) TOPPLE: heatmap (top/bottom 25 regulons × 61 cell types)
(C) SICAI: mean rb vs disease count scatter
(D) DGSA: non-additivity vs disease count scatter (per cell type)
(E) IPA: violin of |log2FC| by stability class
(F) Cross-method correlation matrix (4 per-cell-type metrics)
(G) Ext #1: paired scatter eQTL vs caQTL non-additivity
(H) Ext #1: overlaid distribution histogram

Output: figures/figure1_complete.pdf + .png
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator


# ---------- style ------------------------------------------------------------

BLUE = '#4878CF'
RED = '#D65F5F'
GREY = '#555555'
GREEN = '#6ACC65'
ORANGE = '#D5A021'


def setup_style():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 6,
        'axes.linewidth': 0.5,
        'axes.labelsize': 6.5,
        'axes.titlesize': 8,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'xtick.labelsize': 5.5,
        'ytick.labelsize': 5.5,
        'legend.fontsize': 5.5,
        'lines.linewidth': 0.7,
    })


def compute_gini(arr):
    arr = np.sort(arr.copy())
    n = len(arr)
    if n == 0 or arr.sum() == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * arr) - (n + 1) * np.sum(arr)) / (n * np.sum(arr))


def pval_annot(p):
    if p < 1e-10: return '***'
    elif p < 1e-4: return '**'
    elif p < 0.05: return '*'
    return 'n.s.'


# ---------- data loading -----------------------------------------------------

def load_all():
    """Load all results needed for the figure."""
    data = {}

    # TOPPLE
    data['topple'] = pd.read_csv('results/topple.csv')

    # AUC matrix for heatmap
    s5 = pd.read_excel('data/raw/science.adt3130_table_s5.xlsx',
                        sheet_name='eRegulons_Activators_Exp_AUC_RS')
    data['auc_matrix'] = s5.pivot_table(
        index='eRegulon', columns='cell_type_l4', values='mean_AUC', aggfunc='first')

    # SICAI
    data['sicai'] = pd.read_csv('results/sicai_rb.csv')

    # DGSA
    data['dgsa'] = pd.read_csv('results/dgsa_eqtl.csv')

    # IPA sex data
    sex = pd.read_excel('data/raw/science.adt3130_table_s5.xlsx',
                         sheet_name='sex_difference_results')
    sex['abs_log2FC'] = sex['Log2_Fold_Change'].abs()
    sex['regulon_clean'] = (sex['eRegulon']
                            .str.replace(r'_extended', '', regex=False)
                            .str.replace(r'_\(\d+g\)', '', regex=True))
    topple = data['topple']
    sex_merged = sex.merge(topple[['regulon', 'stability_class']],
                           left_on='regulon_clean', right_on='regulon', how='inner')
    data['ipa_merged'] = sex_merged

    # S15 disease counts
    s15 = pd.read_excel('data/raw/science.adt3130_table_s15.xlsx', sheet_name='Sheet1')
    data['disease_ct'] = s15.groupby('celltype')['trait'].nunique().reset_index()
    data['disease_ct'].columns = ['cell_type', 'disease_count']

    # DGSA per-celltype (from S6)
    s6 = pd.read_csv('data/raw/CIMA_Table_S6.csv',
                      usecols=['phenotype_id', 'celltype', 'slope', 'analysis'])
    eqtl = s6[s6['analysis'] == 'cis-eQTL']
    pivot = eqtl.pivot_table(index='phenotype_id', columns='celltype',
                              values='slope', aggfunc='first')
    pivot = pivot[pivot.notna().sum(axis=1) >= 3]
    filled = pivot.fillna(0)
    gene_na = {g: compute_gini(filled.loc[g].values ** 2) for g in filled.index}
    ct_dgsa = []
    for ct in pivot.columns:
        genes = pivot[ct].dropna().index
        vals = [gene_na[g] for g in genes if g in gene_na]
        if vals:
            ct_dgsa.append({'cell_type': ct, 'mean_na': np.mean(vals)})
    data['dgsa_ct'] = pd.DataFrame(ct_dgsa)

    # Extension #1 matched pairs
    s6_full = pd.read_csv('data/raw/CIMA_Table_S6.csv',
                           usecols=['phenotype_id', 'celltype', 'slope', 'analysis', 'variant_id'])
    eqtl_full = s6_full[s6_full['analysis'] == 'cis-eQTL']
    caqtl_full = s6_full[s6_full['analysis'] == 'cis-caQTL']
    matched = eqtl_full.merge(caqtl_full, on=['variant_id', 'celltype'],
                               suffixes=('_eqtl', '_caqtl'))
    pairs = matched[['phenotype_id_eqtl', 'phenotype_id_caqtl']].drop_duplicates()

    eqtl_cts = sorted(eqtl_full['celltype'].unique())
    caqtl_cts = sorted(caqtl_full['celltype'].unique())
    ep = eqtl_full.pivot_table(index='phenotype_id', columns='celltype',
                                values='slope', aggfunc='first').reindex(columns=eqtl_cts).fillna(0)
    cp = caqtl_full.pivot_table(index='phenotype_id', columns='celltype',
                                 values='slope', aggfunc='first').reindex(columns=caqtl_cts).fillna(0)
    g_na = {g: compute_gini(ep.loc[g].values ** 2) for g in ep.index}
    p_na = {p: compute_gini(cp.loc[p].values ** 2) for p in cp.index}
    pe, pc = [], []
    for _, row in pairs.iterrows():
        g, p = row['phenotype_id_eqtl'], row['phenotype_id_caqtl']
        if g in g_na and p in p_na:
            pe.append(g_na[g]); pc.append(p_na[p])
    data['paired_eqtl'] = np.array(pe)
    data['paired_caqtl'] = np.array(pc)

    # Full distributions
    data['eqtl_na'] = pd.read_csv('results/dgsa_eqtl.csv')['non_additivity'].values
    data['caqtl_na'] = pd.read_csv('results/ext1_caqtl_dgsa.csv')['non_additivity'].values

    return data


# ---------- per-celltype cross-method metrics --------------------------------

def build_cross_method(data):
    """Build per-cell-type DataFrame with 4 HIRA metrics."""
    sicai = data['sicai'][['cell_type', 'mean_rb']].copy()

    # TOPPLE per cell type: mean RI across all regulons (all have 61 ct)
    topple = data['topple']
    auc = data['auc_matrix']
    # For each cell type, mean RI weighted by that cell type's AUC rank
    # Simpler: just mean RI of top-15-RSS regulons per cell type
    # Or: since all regulons appear in all cell types, just use per-regulon RI
    # mapped to cell types by their RSS rank. Simplest: mean RI.
    # All regulons appear in all cell types, so mean RI is same for all.
    # Better: weight by AUC
    ct_topple = []
    for ct in auc.columns:
        aucs = auc[ct]
        ris = topple.set_index('regulon')['mean_RI'].reindex(auc.index)
        valid = aucs.notna() & ris.notna()
        if valid.sum() > 0:
            # Weighted mean RI by AUC
            w = aucs[valid].values
            r = ris[valid].values
            wm = np.average(r, weights=w)
            ct_topple.append({'cell_type': ct, 'weighted_RI': wm})
    topple_ct = pd.DataFrame(ct_topple)

    # DGSA per cell type
    dgsa_ct = data['dgsa_ct'].copy()

    # IPA per cell type: mean |log2FC|
    ipa = data['ipa_merged']
    ipa_ct = ipa.groupby('cell_type_l4')['abs_log2FC'].mean().reset_index()
    ipa_ct.columns = ['cell_type', 'mean_abs_log2FC']

    # Merge all
    merged = sicai.merge(topple_ct, on='cell_type', how='outer')
    merged = merged.merge(dgsa_ct, on='cell_type', how='outer')
    merged = merged.merge(ipa_ct, on='cell_type', how='outer')
    return merged


# ---------- panels -----------------------------------------------------------

def panel_a(ax, topple):
    """TOPPLE bar chart: top 15 regulons by RI."""
    top15 = topple.head(15).copy()
    top15 = top15.iloc[::-1]  # reverse for horizontal bar
    names = [r.replace('_+', '') for r in top15['regulon']]
    bars = ax.barh(range(len(names)), top15['mean_RI'].values,
                   color=BLUE, alpha=0.8, edgecolor='none', height=0.7)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=5)
    ax.set_xlabel('Redistribution Index (RI)')
    ax.set_title('A', fontweight='bold', loc='left', pad=3)
    # Annotate rank
    for i, (_, row) in enumerate(top15.iterrows()):
        ax.text(row['mean_RI'] + 0.0001, i, f"#{int(row['rank'])}",
                va='center', fontsize=4.5, color=GREY)


def panel_b(ax, topple, auc_matrix):
    """TOPPLE heatmap: top/bottom 25 regulons × 61 cell types."""
    top25 = topple.head(25)['regulon'].tolist()
    bot25 = topple.tail(25)['regulon'].tolist()
    selected = top25 + bot25

    sub = auc_matrix.loc[auc_matrix.index.isin(selected)].copy()
    # Reorder: top25 first, then bottom25
    order = [r for r in selected if r in sub.index]
    sub = sub.loc[order]

    # Z-score each row for visualization
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-10, axis=0)

    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                   interpolation='nearest')
    ax.set_xlabel('Cell types (61)')
    ax.set_ylabel('')
    ax.set_xticks([])
    # Show regulon names on y-axis
    names = [r.replace('_+', '') for r in order]
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=3)
    # Divider line between stabilizers and destabilizers
    ax.axhline(len(top25) - 0.5, color='black', lw=0.5, ls='--')
    ax.set_title('B', fontweight='bold', loc='left', pad=3)

    # Colorbar
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.ax.tick_params(labelsize=4.5)
    cb.set_label('Z-score', fontsize=5)


def panel_c(ax, sicai, disease_ct):
    """SICAI scatter: mean rb vs disease count."""
    merged = sicai.merge(disease_ct, on='cell_type', how='inner')
    ax.scatter(merged['mean_rb'], merged['disease_count'], s=8, alpha=0.6,
               color=BLUE, edgecolors='none')
    rho, pval = stats.spearmanr(merged['mean_rb'], merged['disease_count'])
    # Regression line
    z = np.polyfit(merged['mean_rb'], merged['disease_count'], 1)
    x_line = np.linspace(merged['mean_rb'].min(), merged['mean_rb'].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), '--', color=GREY, lw=0.6)
    ax.set_xlabel('Mean $r_b$ (SICAI)')
    ax.set_ylabel('Disease trait count')
    ax.text(0.03, 0.97, f'$\\rho$ = {rho:.3f}\nP = {pval:.1e}\nn = {len(merged)}',
            transform=ax.transAxes, va='top', fontsize=5,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
    ax.set_title('C', fontweight='bold', loc='left', pad=3)


def panel_d(ax, dgsa_ct, disease_ct):
    """DGSA scatter: mean non-additivity vs disease count per cell type."""
    merged = dgsa_ct.merge(disease_ct, on='cell_type', how='inner')
    ax.scatter(merged['mean_na'], merged['disease_count'], s=8, alpha=0.6,
               color=RED, edgecolors='none')
    rho, pval = stats.spearmanr(merged['mean_na'], merged['disease_count'])
    z = np.polyfit(merged['mean_na'], merged['disease_count'], 1)
    x_line = np.linspace(merged['mean_na'].min(), merged['mean_na'].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), '--', color=GREY, lw=0.6)
    ax.set_xlabel('Mean non-additivity (DGSA)')
    ax.set_ylabel('Disease trait count')
    ax.text(0.03, 0.97, f'$\\rho$ = {rho:.3f}\nP = {pval:.1e}\nn = {len(merged)}',
            transform=ax.transAxes, va='top', fontsize=5,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
    ax.set_title('D', fontweight='bold', loc='left', pad=3)


def panel_e(ax, ipa_merged):
    """IPA violin: |log2FC| by stability class."""
    classes = ['stabilizer', 'intermediate', 'destabilizer']
    colors_map = {'stabilizer': BLUE, 'intermediate': GREY, 'destabilizer': RED}
    plot_data = []
    positions = []
    for i, cls in enumerate(classes):
        vals = ipa_merged[ipa_merged['stability_class'] == cls]['abs_log2FC'].dropna().values
        plot_data.append(vals)
        positions.append(i)

    vp = ax.violinplot(plot_data, positions=positions, showmeans=True,
                       showmedians=True, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors_map[classes[i]])
        body.set_alpha(0.6)
    vp['cmeans'].set_color('black')
    vp['cmeans'].set_linewidth(0.5)
    vp['cmedians'].set_color('black')
    vp['cmedians'].set_linewidth(0.5)
    vp['cmedians'].set_linestyle('--')

    ax.set_xticks(positions)
    ax.set_xticklabels(['Stab.', 'Inter.', 'Destab.'], fontsize=5.5)
    ax.set_ylabel('|Log$_2$FC| (sex difference)')
    ax.set_ylim(-0.05, 1.5)

    # Annotate P-value
    s = plot_data[0]
    d = plot_data[2]
    _, p = stats.mannwhitneyu(s, d, alternative='less', method='asymptotic')
    ax.text(0.5, 0.97, f'MW P = {p:.1e}', transform=ax.transAxes,
            ha='center', va='top', fontsize=5)
    ax.set_title('E', fontweight='bold', loc='left', pad=3)


def panel_f(ax, cross):
    """Cross-method correlation matrix."""
    cols = ['mean_rb', 'weighted_RI', 'mean_na', 'mean_abs_log2FC']
    labels = ['SICAI\n$r_b$', 'TOPPLE\nRI', 'DGSA\nNA', 'IPA\n|FC|']
    sub = cross[cols].dropna()
    n = len(cols)
    corr_mat = np.zeros((n, n))
    pval_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                corr_mat[i, j] = 1.0
                pval_mat[i, j] = 0.0
            else:
                r, p = stats.spearmanr(sub[cols[i]], sub[cols[j]])
                corr_mat[i, j] = r
                pval_mat[i, j] = p

    im = ax.imshow(corr_mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, fontsize=5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=5)

    for i in range(n):
        for j in range(n):
            r = corr_mat[i, j]
            sig = pval_annot(pval_mat[i, j]) if i != j else ''
            color = 'white' if abs(r) > 0.5 else 'black'
            ax.text(j, i, f'{r:.2f}\n{sig}', ha='center', va='center',
                    fontsize=4.5, color=color)

    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=4.5)
    cb.set_label('Spearman $\\rho$', fontsize=5)
    ax.set_title('F', fontweight='bold', loc='left', pad=3)


def panel_g(ax, paired_e, paired_c):
    """Paired scatter: eQTL vs caQTL non-additivity at matched loci."""
    ax.scatter(paired_e, paired_c, s=1.5, alpha=0.1, color=GREY,
               linewidths=0, rasterized=True)
    lims = [0.3, 1.02]
    ax.plot(lims, lims, 'k--', lw=0.4, alpha=0.4)
    ax.set_xlabel('eQTL non-additivity')
    ax.set_ylabel('caQTL non-additivity')
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_aspect('equal')

    diff = paired_c - paired_e
    _, w_p = stats.wilcoxon(diff, alternative='greater')
    pct = 100 * np.mean(diff > 0)
    ax.text(0.03, 0.97,
            f'n = {len(paired_e):,}\ncaQTL > eQTL: {pct:.1f}%\nP = {w_p:.1e}',
            transform=ax.transAxes, va='top', fontsize=5,
            bbox=dict(facecolor='white', edgecolor='#ccc', alpha=0.9, lw=0.3))
    ax.set_title('G', fontweight='bold', loc='left', pad=3)


def panel_h(ax, eqtl_na, caqtl_na):
    """Overlaid histograms eQTL vs caQTL non-additivity."""
    bins = np.linspace(0.3, 1.0, 50)
    ax.hist(eqtl_na, bins=bins, alpha=0.55, density=True, color=BLUE,
            label=f'eQTL (n={len(eqtl_na):,})', edgecolor='none')
    ax.hist(caqtl_na, bins=bins, alpha=0.55, density=True, color=RED,
            label=f'caQTL (n={len(caqtl_na):,})', edgecolor='none')
    ax.axvline(np.mean(eqtl_na), color=BLUE, ls='--', lw=0.5)
    ax.axvline(np.mean(caqtl_na), color=RED, ls='--', lw=0.5)
    ax.set_xlabel('Non-additivity')
    ax.set_ylabel('Density')
    ax.legend(frameon=False, loc='upper left', fontsize=5)
    ks_d, ks_p = stats.ks_2samp(eqtl_na, caqtl_na)
    ax.text(0.97, 0.97, f'KS D={ks_d:.3f}\nP={ks_p:.1e}',
            transform=ax.transAxes, ha='right', va='top', fontsize=5,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
    ax.set_title('H', fontweight='bold', loc='left', pad=3)


# ---------- main -------------------------------------------------------------

def main():
    setup_style()

    print("Loading all data ...")
    data = load_all()

    print("Building cross-method metrics ...")
    cross = build_cross_method(data)
    print(f"  Cross-method table: {len(cross)} cell types, cols: {list(cross.columns)}")

    print("Generating figure ...")
    fig = plt.figure(figsize=(180 / 25.4, 220 / 25.4))

    # Layout: 4 rows × 2 cols, but panel B (heatmap) is tall
    # Row 1: A (bar) | B (heatmap, spans rows 1-2)
    # Row 2: C       | B continued
    # Row 3: D | E
    # Row 4: F | G | H  (bottom row has 3 panels)
    gs = fig.add_gridspec(4, 4, hspace=0.45, wspace=0.45,
                          height_ratios=[1, 1, 1, 1])

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0:2, 2:4])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1:3])
    ax_g = fig.add_subplot(gs[3, 0:2])
    ax_h = fig.add_subplot(gs[3, 2:4])

    panel_a(ax_a, data['topple'])
    panel_b(ax_b, data['topple'], data['auc_matrix'])
    panel_c(ax_c, data['sicai'], data['disease_ct'])
    panel_d(ax_d, data['dgsa_ct'], data['disease_ct'])
    panel_e(ax_e, data['ipa_merged'])
    panel_f(ax_f, cross)
    panel_g(ax_g, data['paired_eqtl'], data['paired_caqtl'])
    panel_h(ax_h, data['eqtl_na'], data['caqtl_na'])

    outpath = 'figures/figure1_complete.pdf'
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    fig.savefig(outpath.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close(fig)


if __name__ == '__main__':
    main()
