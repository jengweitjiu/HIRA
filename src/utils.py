#!/usr/bin/env python3
"""
HIRA shared utility functions.

Provides common data loading and metric computation used across scripts.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

# Default data paths
DATA_RAW = Path("data/raw")
RESULTS = Path("results")
FIGURES = Path("figures")


def load_s5(path=None, sheet='eRegulons_Activators_Exp_AUC_RS'):
    """Load S5 AUC matrix (203 regulons x 61 cell types)."""
    path = path or DATA_RAW / "science.adt3130_table_s5.xlsx"
    df = pd.read_excel(path, sheet_name=sheet)
    if 'eRegulon' in df.columns:
        df = df.rename(columns={'eRegulon': 'regulon', 'cell_type_l4': 'cell_type'})
    return df


def load_s6(path=None):
    """Load S6 cis-eQTL + cis-caQTL (223K rows)."""
    path = path or DATA_RAW / "CIMA_Table_S6.csv"
    return pd.read_csv(path)


def load_s8(path=None, sheet='cis_eQTL'):
    """Load S8 pairwise r_b sharing matrix."""
    path = path or DATA_RAW / "science.adt3130_table_s8.xlsx"
    return pd.read_excel(path, sheet_name=sheet)


def load_s15(path=None):
    """Load S15 SMR pleiotropic associations (2,085 rows)."""
    path = path or DATA_RAW / "science.adt3130_table_s15.xlsx"
    return pd.read_excel(path)


def load_topple_classes(path=None):
    """Load TOPPLE regulon stability classes. Returns dict: regulon -> class."""
    path = path or RESULTS / "topple.csv"
    topple = pd.read_csv(path)
    return topple.set_index('regulon')['stability_class'].to_dict()


def load_topple_ri(path=None):
    """Load TOPPLE RI scores. Returns dict: regulon -> mean_RI."""
    path = path or RESULTS / "topple.csv"
    topple = pd.read_csv(path)
    return topple.set_index('regulon')['mean_RI'].to_dict()


def gini(values):
    """Compute Gini coefficient of an array."""
    v = np.sort(np.abs(np.asarray(values, dtype=float)))
    n = len(v)
    if n < 2 or v.sum() == 0:
        return np.nan
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * v) - (n + 1) * np.sum(v)) / (n * np.sum(v))


def dgsa_metrics(slopes):
    """Compute DGSA geometric metrics from a vector of effect sizes.

    Parameters
    ----------
    slopes : array-like
        Cell-type-specific eQTL/caQTL slopes for one gene/peak.

    Returns
    -------
    dict with keys: magnitude, uniformity, non_additivity, sparsity, gini
    """
    s = np.asarray(slopes, dtype=float)
    n_total = 69  # total CIMA cell types

    magnitude = np.sqrt(np.sum(s ** 2))
    if magnitude == 0:
        return {'magnitude': 0, 'uniformity': 0, 'non_additivity': 0, 'sparsity': 1, 'gini': 0}

    uniform = np.ones(len(s)) / np.sqrt(len(s))
    s_norm = s / magnitude
    uniformity = np.dot(s_norm, uniform)

    non_additivity = gini(s ** 2)
    sparsity = 1 - len(s) / n_total
    gini_val = gini(np.abs(s))

    return {
        'magnitude': magnitude,
        'uniformity': uniformity,
        'non_additivity': non_additivity,
        'sparsity': sparsity,
        'gini': gini_val,
    }


def build_rb_matrix(s8_df):
    """Build symmetric r_b matrix from S8 pairwise data.

    Parameters
    ----------
    s8_df : DataFrame with columns: reference_cell_type, query_celltype, rb

    Returns
    -------
    DataFrame (n_ct x n_ct) symmetric matrix with 1.0 diagonal
    """
    cts = sorted(set(s8_df['reference_cell_type']) | set(s8_df['query_celltype']))
    mat = pd.DataFrame(0.0, index=cts, columns=cts)
    for _, row in s8_df.iterrows():
        mat.loc[row['reference_cell_type'], row['query_celltype']] = row['rb']
        mat.loc[row['query_celltype'], row['reference_cell_type']] = row['rb']
    np.fill_diagonal(mat.values, 1.0)
    return mat


def nature_style_figure(nrows=1, ncols=1, width_mm=180, height_mm=None, **kwargs):
    """Create a figure with Nature-style formatting.

    Parameters
    ----------
    nrows, ncols : int
        Subplot grid dimensions.
    width_mm : float
        Figure width in mm (Nature single column = 89mm, double = 180mm).
    height_mm : float, optional
        Figure height in mm. If None, auto-calculated.

    Returns
    -------
    fig, axes
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
    })

    w = width_mm / 25.4
    h = (height_mm / 25.4) if height_mm else w * nrows / ncols * 0.8

    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h), **kwargs)
    return fig, axes


def add_panel_labels(axes, labels=None):
    """Add A, B, C... panel labels to axes array."""
    if labels is None:
        labels = [chr(65 + i) for i in range(len(axes.flat))]
    for ax, label in zip(axes.flat, labels):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')


def clean_spines(ax, keep=('bottom', 'left')):
    """Remove unwanted spines from axis."""
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(spine in keep)
