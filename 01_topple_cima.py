#!/usr/bin/env python3
"""
TOPPLE Stability Analysis on CIMA Table S5
==========================================
Applies TOPPLE (Transcriptomic Perturbation-based Phenotype Landscape Explorer)
to 203 eRegulons × 61 cell types from Yin et al., Science 2026.

Core algorithm:
  1. Build regulon × cell_type AUC matrix from Table S5
  2. Leave-one-out perturbation: remove each regulon, measure redistribution
  3. Redistribution Index (RI) = Jensen-Shannon divergence of renormalized profile
  4. Cross-cell-type stability score = mean RI across all cell types
  5. Rank: high RI → stabilizer (system depends on it); low RI → dispensable

Author: Jeng-Wei Tjiu, M.D.
Date: 2026-03
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr, zscore
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
import warnings
import sys
import io
warnings.filterwarnings('ignore')
# Fix Windows encoding for Unicode output
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Configuration ──────────────────────────────────────────────────────
DATA_DIR = Path("data")
FIG_DIR = Path("figures")
RESULTS_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# Nature-style color palette
COLORS = {
    'stabilizer': '#c0392b',
    'destabilizer': '#2980b9',
    'neutral': '#95a5a6',
    'accent': '#27ae60',
}

# ── Step 1: Load Table S5 ─────────────────────────────────────────────
def load_table_s5(data_dir: Path) -> pd.DataFrame:
    """Load CIMA Table S5 (203 eRegulons × 61 cell types).
    
    Tries multiple possible filenames and formats.
    Returns raw dataframe for inspection.
    """
    candidates = [
        "science_adt3130_table_s5.xlsx",
        "table_s5.xlsx",
        "Table_S5.xlsx",
        "science.adt3130_table_s5.xlsx",
        "adt3130_table_s5.xlsx",
    ]
    
    for fname in candidates:
        fpath = data_dir / fname
        if fpath.exists():
            print(f"✓ Loading: {fpath}")
            # Try multiple sheets — prefer the AUC/RSS data sheet
            xls = pd.ExcelFile(fpath)
            print(f"  Sheets: {xls.sheet_names}")
            # Look for the main data sheet with AUC/RSS per regulon × cell type
            target_sheet = 0
            for idx, sname in enumerate(xls.sheet_names):
                if 'auc' in sname.lower() or 'activator' in sname.lower():
                    target_sheet = idx
                    break
            df = pd.read_excel(fpath, sheet_name=target_sheet)
            print(f"  Using sheet: {xls.sheet_names[target_sheet]}")
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns[:10])}...")
            return df
    
    # Also try CSV
    for fname in ["science_adt3130_table_s5.csv", "table_s5.csv"]:
        fpath = data_dir / fname
        if fpath.exists():
            print(f"✓ Loading CSV: {fpath}")
            df = pd.read_csv(fpath)
            print(f"  Shape: {df.shape}")
            return df
    
    raise FileNotFoundError(
        f"Table S5 not found in {data_dir}/. "
        f"Expected one of: {candidates}\n"
        f"Available files: {list(data_dir.glob('*'))}"
    )


def build_auc_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Extract regulon × cell_type AUC matrix from Table S5.
    
    Table S5 has 12,383 entries = 203 regulons × 61 cell types.
    Expected columns include: regulon/TF name, cell type, AUC score, RSS.
    
    This function auto-detects column names.
    """
    cols = [c.lower().strip() for c in df.columns]
    col_map = {c: orig for c, orig in zip(cols, df.columns)}
    
    # Auto-detect regulon column
    regulon_col = None
    for candidate in ['eregulon', 'regulon', 'tf', 'tf_name', 'gene', 'eregulon_name']:
        if candidate in cols:
            regulon_col = col_map[candidate]
            break
    
    # Auto-detect cell type column
    celltype_col = None
    for candidate in ['cell_type_l4', 'cell_type', 'celltype', 'cell type', 'l4', 'l4_celltype', 'cluster']:
        if candidate in cols:
            celltype_col = col_map[candidate]
            break
    
    # Auto-detect AUC column
    auc_col = None
    for candidate in ['auc', 'auc_score', 'mean_auc', 'regulon_auc']:
        if candidate in cols:
            auc_col = col_map[candidate]
            break
    
    # Auto-detect RSS column
    rss_col = None
    for candidate in ['rss', 'rss_score', 'regulon_specificity_score']:
        if candidate in cols:
            rss_col = col_map[candidate]
            break
    
    # If auto-detect fails, print columns and let user adjust
    if any(x is None for x in [regulon_col, celltype_col, auc_col]):
        print("\n⚠ Column auto-detection incomplete. Available columns:")
        for i, c in enumerate(df.columns):
            print(f"  [{i}] {c} — dtype={df[c].dtype}, sample: {df[c].iloc[0]}")
        print("\n→ Please update column mapping in build_auc_matrix()")
        
        # Fallback: assume first few columns
        print("\n  Attempting positional fallback...")
        regulon_col = regulon_col or df.columns[0]
        celltype_col = celltype_col or df.columns[1]
        auc_col = auc_col or df.columns[2]
    
    print(f"\n✓ Column mapping:")
    print(f"  Regulon:   {regulon_col}")
    print(f"  Cell type: {celltype_col}")
    print(f"  AUC:       {auc_col}")
    if rss_col:
        print(f"  RSS:       {rss_col}")
    
    # Pivot to matrix
    auc_matrix = df.pivot_table(
        index=regulon_col, 
        columns=celltype_col, 
        values=auc_col,
        aggfunc='mean'
    )
    
    print(f"\n✓ AUC matrix: {auc_matrix.shape[0]} regulons × {auc_matrix.shape[1]} cell types")
    print(f"  Expected: 203 × 61")
    
    # Fill NaN with 0 (regulon not active in that cell type)
    auc_matrix = auc_matrix.fillna(0)
    
    return auc_matrix, rss_col, df, regulon_col, celltype_col


# ── Step 2: TOPPLE Core Algorithm ─────────────────────────────────────
def compute_redistribution_index(auc_matrix: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-out perturbation → redistribution index per regulon per cell type.
    
    For each cell type j:
      1. Original profile: p_j = AUC vector across all regulons (L1-normalized)
      2. For each regulon i:
         - Remove regulon i → renormalize remaining → q_j^{-i}
         - RI(i, j) = JSD(p_j, q_j^{-i})  [Jensen-Shannon divergence]
      3. High RI(i,j) = removing regulon i strongly disrupts cell type j's landscape
    
    Returns DataFrame: regulons × cell_types of RI values
    """
    n_reg, n_ct = auc_matrix.shape
    regulons = auc_matrix.index.tolist()
    celltypes = auc_matrix.columns.tolist()
    
    ri_matrix = np.zeros((n_reg, n_ct))
    
    print(f"\n⏳ Computing redistribution index ({n_reg} regulons × {n_ct} cell types)...")
    
    for j, ct in enumerate(celltypes):
        # Original profile (L1-normalized)
        profile = auc_matrix[ct].values.copy()
        profile_sum = profile.sum()
        if profile_sum == 0:
            continue
        p = profile / profile_sum
        
        for i in range(n_reg):
            # Leave-one-out: mask regulon i
            perturbed = profile.copy()
            perturbed[i] = 0
            perturbed_sum = perturbed.sum()
            
            if perturbed_sum == 0:
                ri_matrix[i, j] = 1.0  # Total collapse
                continue
            
            q = perturbed / perturbed_sum
            
            # Jensen-Shannon divergence (symmetric, bounded [0, 1])
            # Add small epsilon to avoid log(0)
            eps = 1e-12
            p_safe = p + eps
            q_safe = q + eps
            p_safe /= p_safe.sum()
            q_safe /= q_safe.sum()
            
            ri_matrix[i, j] = jensenshannon(p_safe, q_safe) ** 2  # squared for JSD
    
    ri_df = pd.DataFrame(ri_matrix, index=regulons, columns=celltypes)
    
    print(f"✓ RI matrix computed: {ri_df.shape}")
    print(f"  RI range: [{ri_df.values.min():.6f}, {ri_df.values.max():.6f}]")
    print(f"  Mean RI: {ri_df.values.mean():.6f}")
    
    return ri_df


def compute_stability_scores(ri_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-regulon stability metrics across cell types.
    
    Returns DataFrame with columns:
    - mean_RI: average redistribution index (higher = more stabilizing)
    - max_RI: maximum RI across any cell type
    - std_RI: variability (high = lineage-specific stabilizer)
    - n_critical: number of cell types where RI > 95th percentile
    - stability_rank: rank by mean_RI (1 = top stabilizer)
    - role: 'stabilizer' / 'destabilizer' / 'neutral'
    """
    threshold_95 = np.percentile(ri_df.values, 95)
    
    scores = pd.DataFrame({
        'mean_RI': ri_df.mean(axis=1),
        'median_RI': ri_df.median(axis=1),
        'max_RI': ri_df.max(axis=1),
        'std_RI': ri_df.std(axis=1),
        'cv_RI': ri_df.std(axis=1) / (ri_df.mean(axis=1) + 1e-12),
        'n_critical': (ri_df > threshold_95).sum(axis=1),
    })
    
    scores = scores.sort_values('mean_RI', ascending=False)
    scores['stability_rank'] = range(1, len(scores) + 1)
    
    # Classify roles
    q75 = scores['mean_RI'].quantile(0.75)
    q25 = scores['mean_RI'].quantile(0.25)
    scores['role'] = 'neutral'
    scores.loc[scores['mean_RI'] >= q75, 'role'] = 'stabilizer'
    scores.loc[scores['mean_RI'] <= q25, 'role'] = 'destabilizer'
    
    return scores


# ── Step 3: Cross-lineage analysis ────────────────────────────────────
def lineage_analysis(ri_df: pd.DataFrame, df_raw: pd.DataFrame) -> dict:
    """Group cell types by immune lineage (L1) and analyze per-lineage stability.
    
    CIMA L1 categories: T, B, NK, Mono, DC, Mega, Ery, Plasma, etc.
    """
    # Infer lineage from cell type names
    lineage_map = {}
    for ct in ri_df.columns:
        ct_lower = ct.lower()
        if any(x in ct_lower for x in ['cd4', 'cd8', 'treg', 'mait', 'proliferative_t', 'gdt', 'gd_t', 'gdt1', 'gdt2', 'inkt', 'isg_t']):
            lineage_map[ct] = 'T cell'
        elif any(x in ct_lower for x in ['nk', 'ilc']):
            lineage_map[ct] = 'NK/ILC'
        elif any(x in ct_lower for x in ['bn_', 'bm_', 'transitional_b', 'switched_', 'unswitched_', 'pre-switched', 'b_mme', 'b_cd1c']):
            lineage_map[ct] = 'B cell'
        elif any(x in ct_lower for x in ['plasma', 'plasmablast']):
            lineage_map[ct] = 'Plasma'
        elif any(x in ct_lower for x in ['mono', 'cmono', 'ncmono', 'intmono']):
            lineage_map[ct] = 'Monocyte'
        elif any(x in ct_lower for x in ['dc', 'pdc']):
            lineage_map[ct] = 'DC'
        elif 'mega' in ct_lower or 'plt' in ct_lower:
            lineage_map[ct] = 'Megakaryocyte'
        elif 'ery' in ct_lower or 'rbc' in ct_lower:
            lineage_map[ct] = 'Erythroid'
        elif 'hsc' in ct_lower or 'hspc' in ct_lower or 'prog' in ct_lower:
            lineage_map[ct] = 'Progenitor'
        else:
            lineage_map[ct] = 'Other'
    
    # Per-lineage mean RI
    lineage_ri = {}
    for lineage in set(lineage_map.values()):
        cts = [ct for ct, lin in lineage_map.items() if lin == lineage]
        if cts:
            lineage_ri[lineage] = ri_df[cts].mean(axis=1)
    
    return lineage_map, lineage_ri


# ── Step 4: Publication Figures ────────────────────────────────────────
def setup_figure_style():
    """Nature-style figure formatting."""
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })


def plot_stability_landscape(ri_df, scores, fig_dir):
    """Figure 1: TOPPLE Stability Landscape Heatmap.
    
    Rows = regulons (sorted by mean RI), columns = cell types (grouped by lineage).
    Top stabilizers and destabilizers labeled.
    """
    setup_figure_style()
    
    # Sort regulons by mean RI
    order = scores.sort_values('mean_RI', ascending=False).index
    top_n = min(50, len(order))  # Show top 50 for readability
    
    # Select top stabilizers and destabilizers
    top_stab = order[:top_n // 2]
    top_dest = order[-top_n // 2:]
    selected = list(top_stab) + list(top_dest)
    
    ri_plot = ri_df.loc[selected]
    
    # Z-score for visualization (row-wise)
    ri_z = ri_plot.sub(ri_plot.mean(axis=1), axis=0).div(ri_plot.std(axis=1) + 1e-12, axis=0)
    ri_z = ri_z.fillna(0)
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    sns.heatmap(
        ri_z,
        cmap='RdBu_r',
        center=0,
        vmin=-2, vmax=2,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Redistribution Index (z-score)', 'shrink': 0.6},
        ax=ax
    )
    
    ax.set_xlabel('Cell Types (CIMA L4)', fontsize=10)
    ax.set_ylabel('eRegulons', fontsize=10)
    ax.set_title(
        'TOPPLE Stability Landscape — CIMA Atlas\n'
        f'Top {top_n//2} stabilizers (top) vs destabilizers (bottom)',
        fontsize=11, fontweight='bold'
    )
    
    # Highlight top 5 stabilizer names
    top5 = list(scores.head(5).index)
    for i, label in enumerate(ax.get_yticklabels()):
        if label.get_text() in top5:
            label.set_color(COLORS['stabilizer'])
            label.set_fontweight('bold')
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig1_topple_stability_landscape.{ext}')
    plt.close()
    print(f"✓ Figure 1 saved: fig1_topple_stability_landscape.png/pdf")


def plot_stabilizer_ranking(scores, fig_dir):
    """Figure 1b: Bar chart of top 20 stabilizers and bottom 20 destabilizers."""
    setup_figure_style()
    
    top20 = scores.head(20)
    bot20 = scores.tail(20).iloc[::-1]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Top stabilizers
    colors1 = [COLORS['stabilizer']] * len(top20)
    ax1.barh(range(len(top20)), top20['mean_RI'], color=colors1, edgecolor='white', height=0.7)
    ax1.set_yticks(range(len(top20)))
    ax1.set_yticklabels(top20.index, fontsize=7)
    ax1.set_xlabel('Mean Redistribution Index')
    ax1.set_title('Top 20 Stabilizers', fontweight='bold', color=COLORS['stabilizer'])
    ax1.invert_yaxis()
    
    # Check for STAT5B (regulon names may have _+ or _- suffix)
    stat5b_hits = [r for r in top20.index if r.upper().startswith('STAT5B')]
    for stat5b_key in stat5b_hits:
        idx = list(top20.index).index(stat5b_key)
        ax1.get_yticklabels()[idx].set_fontweight('bold')
        ax1.get_yticklabels()[idx].set_fontsize(9)
        ax1.annotate('★ #1 in psoriasis',
                     xy=(top20.loc[stat5b_key, 'mean_RI'], idx),
                     xytext=(top20['mean_RI'].max() * 0.7, idx + 2),
                     fontsize=7, color=COLORS['stabilizer'],
                     arrowprops=dict(arrowstyle='->', color=COLORS['stabilizer']))
    
    # Bottom destabilizers
    colors2 = [COLORS['destabilizer']] * len(bot20)
    ax2.barh(range(len(bot20)), bot20['mean_RI'], color=colors2, edgecolor='white', height=0.7)
    ax2.set_yticks(range(len(bot20)))
    ax2.set_yticklabels(bot20.index, fontsize=7)
    ax2.set_xlabel('Mean Redistribution Index')
    ax2.set_title('Bottom 20 (Dispensable)', fontweight='bold', color=COLORS['destabilizer'])
    ax2.invert_yaxis()
    
    plt.suptitle('TOPPLE Regulon Stability Ranking — CIMA Atlas (203 eRegulons × 61 cell types)',
                 fontsize=11, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig1b_stabilizer_ranking.{ext}')
    plt.close()
    print(f"✓ Figure 1b saved: fig1b_stabilizer_ranking.png/pdf")


def plot_lineage_comparison(ri_df, lineage_map, scores, fig_dir):
    """Figure 1c: Per-lineage stability profiles of top regulons."""
    setup_figure_style()
    
    # Get unique lineages
    lineages = sorted(set(lineage_map.values()))
    top10 = scores.head(10).index.tolist()
    
    # Compute per-lineage mean RI for top regulons
    lineage_data = {}
    for lin in lineages:
        cts = [ct for ct, l in lineage_map.items() if l == lin]
        if cts:
            lineage_data[lin] = ri_df.loc[top10, cts].mean(axis=1)
    
    lin_df = pd.DataFrame(lineage_data)
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.heatmap(
        lin_df,
        cmap='YlOrRd',
        annot=True, fmt='.4f',
        xticklabels=True, yticklabels=True,
        cbar_kws={'label': 'Mean RI', 'shrink': 0.6},
        ax=ax
    )
    ax.set_xlabel('Immune Lineage')
    ax.set_ylabel('Top 10 Stabilizer eRegulons')
    ax.set_title('Lineage-Specific Stability — Top 10 TOPPLE Stabilizers', fontweight='bold')
    
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig1c_lineage_stability.{ext}')
    plt.close()
    print(f"✓ Figure 1c saved: fig1c_lineage_stability.png/pdf")


# ── Main ───────────────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("TOPPLE Stability Analysis on CIMA Atlas (Science 2026)")
    print("203 eRegulons × 61 cell types — Table S5")
    print("=" * 70)
    
    # 1. Load data
    df_raw = load_table_s5(DATA_DIR)
    
    # 2. Build AUC matrix
    auc_matrix, rss_col, df, reg_col, ct_col = build_auc_matrix(df_raw)
    
    # 3. Compute redistribution index
    ri_df = compute_redistribution_index(auc_matrix)
    
    # 4. Stability scores
    scores = compute_stability_scores(ri_df)
    
    # 5. Print top results
    print("\n" + "=" * 70)
    print("TOP 10 STABILIZERS (highest mean RI — system depends on them)")
    print("=" * 70)
    print(scores.head(10)[['mean_RI', 'max_RI', 'std_RI', 'n_critical', 'stability_rank']].to_string())
    
    print("\n" + "=" * 70)
    print("BOTTOM 10 (lowest mean RI — dispensable)")
    print("=" * 70)
    print(scores.tail(10)[['mean_RI', 'max_RI', 'std_RI', 'n_critical', 'stability_rank']].to_string())
    
    # 6. Check STAT5B
    # Try exact and prefix match for STAT5B (regulon names have _+ or _- suffix)
    stat5b_matches = [r for r in scores.index if r.upper().startswith('STAT5B')]
    if stat5b_matches:
        stat5b_key = stat5b_matches[0]
        rank = scores.loc[stat5b_key, 'stability_rank']
        ri = scores.loc[stat5b_key, 'mean_RI']
        print(f"\n★ {stat5b_key}: rank #{rank}, mean RI = {ri:.6f}")
        print(f"  (In psoriasis GSE173706: #1 stabilizer, RI 0.03→0.61)")
        if len(stat5b_matches) > 1:
            for m in stat5b_matches[1:]:
                print(f"  {m}: rank #{scores.loc[m, 'stability_rank']}, RI = {scores.loc[m, 'mean_RI']:.6f}")
    elif 'STAT5B' in scores.index:
        rank = scores.loc['STAT5B', 'stability_rank']
        ri = scores.loc['STAT5B', 'mean_RI']
        print(f"\n★ STAT5B: rank #{rank}, mean RI = {ri:.6f}")
        print(f"  (In psoriasis GSE173706: #1 stabilizer, RI 0.03→0.61)")
    else:
        # Try partial match
        matches = [r for r in scores.index if 'STAT5' in r.upper()]
        if matches:
            print(f"\n★ STAT5-related regulons found: {matches}")
            for m in matches:
                print(f"  {m}: rank #{scores.loc[m, 'stability_rank']}, RI = {scores.loc[m, 'mean_RI']:.6f}")
        else:
            print("\n⚠ STAT5B not found in 203 eRegulons. Regulon naming may differ.")
    
    # 7. Lineage analysis
    lineage_map, lineage_ri = lineage_analysis(ri_df, df_raw)
    lineage_counts = pd.Series(lineage_map).value_counts()
    print(f"\n✓ Lineage breakdown:")
    for lin, count in lineage_counts.items():
        print(f"  {lin}: {count} cell types")
    
    # 8. Generate figures
    print("\n⏳ Generating figures...")
    plot_stability_landscape(ri_df, scores, FIG_DIR)
    plot_stabilizer_ranking(scores, FIG_DIR)
    plot_lineage_comparison(ri_df, lineage_map, scores, FIG_DIR)
    
    # 9. Save results
    scores.to_csv(RESULTS_DIR / "topple_stability_scores.csv")
    ri_df.to_csv(RESULTS_DIR / "topple_ri_matrix.csv")
    print(f"\n✓ Results saved to {RESULTS_DIR}/")
    
    # 10. Summary for Day 2
    print("\n" + "=" * 70)
    print("DAY 1 COMPLETE — Ready for Day 2 (SICAI on Table S8)")
    print("=" * 70)
    print(f"  Stabilizers identified: {(scores['role'] == 'stabilizer').sum()}")
    print(f"  Destabilizers: {(scores['role'] == 'destabilizer').sum()}")
    print(f"  Neutral: {(scores['role'] == 'neutral').sum()}")
    print(f"  Figures: {list(FIG_DIR.glob('fig1*'))}")
    
    return ri_df, scores, lineage_map


if __name__ == "__main__":
    ri_df, scores, lineage_map = main()
