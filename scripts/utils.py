"""
HIRA shared utilities — figure styling, column detection, data loading,
scatter plot helpers, and Windows encoding fix.
"""

import sys
import io
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from pathlib import Path


# ── Paths ─────────────────────────────────────────────────────────────
DATA_DIR = Path("data")
FIG_DIR = Path("figures")
RESULTS_DIR = Path("results")

# ── Color palette (Nature-style) ─────────────────────────────────────
COLORS = {
    'stabilizer': '#c0392b',
    'destabilizer': '#2980b9',
    'neutral': '#95a5a6',
    'accent': '#27ae60',
    'primary': '#2c3e50',
    'regression': '#e74c3c',
    'secondary': '#3498db',
}


def setup_stdout():
    """Fix UTF-8 encoding on Windows (cp950/cp1252 can't handle Unicode)."""
    if hasattr(sys.stdout, 'buffer'):
        sys.stdout = io.TextIOWrapper(
            sys.stdout.buffer, encoding='utf-8', errors='replace')


def ensure_dirs():
    """Create output directories if they don't exist."""
    FIG_DIR.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)


# ── Figure styling ────────────────────────────────────────────────────
def setup_figure_style():
    """Set Nature-style matplotlib rcParams (Arial, 300 DPI, compact labels)."""
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
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
    })


def save_figure(fig, fig_dir, name):
    """Save figure as both PNG and PDF."""
    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'{name}.{ext}')
    plt.close(fig)


# ── Column detection ──────────────────────────────────────────────────
def detect_column(df, candidates, label="column"):
    """Auto-detect a column name from a list of candidates (case-insensitive).

    Parameters
    ----------
    df : DataFrame
    candidates : list of str
        Lowercase candidate names to search for.
    label : str
        Human-readable label for error messages.

    Returns
    -------
    str or None
        The original column name if found, else None.
    """
    cols_lower = [c.lower().strip() for c in df.columns]
    col_map = {cl: orig for cl, orig in zip(cols_lower, df.columns)}
    for candidate in candidates:
        if candidate in col_map:
            return col_map[candidate]
    return None


# ── File finding ──────────────────────────────────────────────────────
# Canonical filename candidates per table
TABLE_CANDIDATES = {
    's5': [
        "science.adt3130_table_s5.xlsx",
        "science_adt3130_table_s5.xlsx",
        "table_s5.xlsx",
        "Table_S5.xlsx",
    ],
    's6': [
        "CIMA_Table_S6.csv",
        "science.adt3130_table_s6.csv",
        "table_s6.csv",
    ],
    's8': [
        "science.adt3130_table_s8.xlsx",
        "science_adt3130_table_s8.xlsx",
        "table_s8.xlsx",
        "Table_S8.xlsx",
    ],
    's15': [
        "science.adt3130_table_s15.xlsx",
        "science_adt3130_table_s15.xlsx",
        "table_s15.xlsx",
        "Table_S15.xlsx",
    ],
}


def find_file(data_dir, table_key=None, candidates=None):
    """Find a data file from a list of candidate filenames.

    Parameters
    ----------
    data_dir : Path
    table_key : str, optional
        Key into TABLE_CANDIDATES (e.g. 's5', 's8').
    candidates : list of str, optional
        Explicit filename candidates (overrides table_key).

    Returns
    -------
    Path
        Path to the found file.

    Raises
    ------
    FileNotFoundError
        If no candidate file exists.
    """
    if candidates is None:
        candidates = TABLE_CANDIDATES.get(table_key, [])
    for fname in candidates:
        fpath = Path(data_dir) / fname
        if fpath.exists():
            return fpath
    available = list(Path(data_dir).glob('*'))
    raise FileNotFoundError(
        f"No file found for {table_key or candidates}. "
        f"Available: {[f.name for f in available]}"
    )


# ── Scatter + regression panel ────────────────────────────────────────
def scatter_panel(ax, x, y, xlabel, ylabel, title,
                  label_top=5, label_bot=0):
    """Reusable scatter plot with regression line and Spearman stats box.

    Parameters
    ----------
    ax : matplotlib Axes
    x, y : array-like (with index for labeling)
    xlabel, ylabel, title : str
    label_top : int
        Number of top-y points to label.
    label_bot : int
        Number of bottom-x points to label.

    Returns
    -------
    rho, p : float
        Spearman correlation and p-value.
    """
    rho, p = spearmanr(x, y)
    ax.scatter(x, y, s=40, alpha=0.7, c=COLORS['primary'],
               edgecolors='white', linewidth=0.5, zorder=3)

    # Regression line
    z = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, np.poly1d(z)(x_line), '--',
            color=COLORS['regression'], linewidth=1.5)

    # Labels
    df_tmp = pd.DataFrame({'x': x, 'y': y})
    if label_top > 0:
        for idx, row in df_tmp.nlargest(label_top, 'y').iterrows():
            ax.annotate(idx, (row['x'], row['y']), fontsize=5.5,
                        ha='left', va='bottom', xytext=(4, 3),
                        textcoords='offset points')
    if label_bot > 0:
        for idx, row in df_tmp.nsmallest(label_bot, 'x').iterrows():
            ax.annotate(idx, (row['x'], row['y']), fontsize=5.5,
                        ha='right', va='top', xytext=(-4, -3),
                        textcoords='offset points', color='#7f8c8d')

    # Stats box
    sig = format_sig(p)
    stats = f'$\\rho$ = {rho:.3f} {sig}\nP = {p:.2e}\nn = {len(x)}'
    ax.text(0.05, 0.95, stats, transform=ax.transAxes, fontsize=8, va='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.8))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold', loc='left')
    ax.grid(True, alpha=0.2)
    return rho, p


def format_sig(p):
    """Return significance stars for a p-value."""
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'n.s.'


def gini_coefficient(x):
    """Gini coefficient of absolute values in array x.

    Takes |x| before computing, so mixed-sign effect sizes are handled
    as magnitudes. Returns 0 for empty or all-zero arrays.

    Parameters
    ----------
    x : array-like
        Values (may be negative; absolute values are used).

    Returns
    -------
    float
        Gini coefficient in [0, 1]. 0 = perfect equality, 1 = max inequality.
    """
    x = np.abs(np.asarray(x, dtype=float))
    x = np.sort(x)
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * x) / (n * np.sum(x))) - (n + 1) / n


def check_prerequisites(results_dir, required_files):
    """Check that all required result files exist before proceeding.

    Parameters
    ----------
    results_dir : Path
    required_files : list of str

    Raises
    ------
    SystemExit
        If any file is missing.
    """
    missing = [f for f in required_files if not (results_dir / f).exists()]
    if missing:
        for f in missing:
            print(f"  Missing: {results_dir / f}")
        print("  Run prerequisite scripts first.")
        raise SystemExit(1)
