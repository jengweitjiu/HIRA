# 🚀 Quick Start — 72-Hour Execution Checklist

## Pre-flight (10 min)

```bash
# 1. Create project
mkdir -p ~/CIMA_Breakthrough/{data,scripts,results,figures}

# 2. Copy these files into the project
cp CLAUDE.md ~/CIMA_Breakthrough/
cp scripts/*.py ~/CIMA_Breakthrough/scripts/

# 3. Copy CIMA tables into data/
cp science_adt3130_table_s5.xlsx ~/CIMA_Breakthrough/data/
cp science_adt3130_table_s8.xlsx ~/CIMA_Breakthrough/data/
cp science_adt3130_table_s15.xlsx ~/CIMA_Breakthrough/data/

# 4. Install dependencies
pip install pandas numpy scipy matplotlib seaborn openpyxl scikit-learn
```

## Option A: Run with Claude Code (recommended)

```bash
cd ~/CIMA_Breakthrough
claude

# Then in Claude Code session:
> Run scripts/01_topple_cima.py. If any column detection fails, 
  inspect the Excel file and fix the column mapping automatically.
```

Claude Code will execute → debug → fix → re-run automatically.

## Option B: Run directly / on Colab

```bash
cd ~/CIMA_Breakthrough
python scripts/01_topple_cima.py   # Day 1
python scripts/02_sicai_cima.py    # Day 2
python scripts/03_combined_figure.py  # Day 3
```

## Expected Outputs

```
results/
├── topple_stability_scores.csv     # 203 regulons ranked
├── topple_ri_matrix.csv            # 203 × 61 RI values
├── sicai_coupling_complexity.csv   # 69 cell types scored
├── sicai_rb_matrix.csv             # 69 × 69 correlation
├── sicai_complexity_disease_merged.csv  # The money result
└── genome_biology_framing.md       # 500-word intro draft

figures/
├── fig1_topple_stability_landscape.pdf
├── fig1b_stabilizer_ranking.pdf
├── fig1c_lineage_stability.pdf
├── fig2a_rb_heatmap.pdf
├── fig2b_complexity_vs_disease.pdf  ← THE key figure
├── fig2c_complexity_distribution.pdf
└── fig_combined_proof_of_concept.pdf  ← For submission
```

## 🔴 Don't Forget

- [ ] **DGSA transfer → Bioinformatics Advances** (30-day expiry from 7-Mar-2026)
- [ ] TOPPLE: add references, Zenodo deposit, push GitHub, submit to NatMethods
- [ ] SICAI: pre-submit checks (Fig3c-d, Fig5, Ref11) before Genome Biology
