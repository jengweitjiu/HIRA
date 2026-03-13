# CLAUDE.md — CIMA × Methods Breakthrough Project

## Project Goal
Apply TOPPLE, SICAI, DGSA frameworks to CIMA (Yin et al., Science 2026, 10.1126/science.adt3130)
published supplementary tables. Generate proof-of-concept figures for Genome Biology / NatComms.

## Key Constraint
- NO raw data access — all analyses use published supplementary tables ONLY
- Environment: Colab Pro+ (T4 GPU, <64GB RAM) or local Python 3.10+
- Single author: Dr. Jeng-Wei Tjiu, Dept. Dermatology, NTUH

## Data Files (in data/)
| File | Contents | Rows | Use |
|------|----------|------|-----|
| table_s5 | 203 eRegulons × 61 cell types (AUC, RSS, age_corr, sex_diff) | 12,383 | TOPPLE |
| table_s6 | 223,405 lead xQTLs (eQTL + caQTL) | 223,405 | DGSA |
| table_s8 | 4,692 eQTL sharing pairs (π₁, r_b) | 4,692 | SICAI |
| table_s15 | 2,085 SMR pleiotropic associations (68 traits × 68 cell types) | 2,085 | Disease correlation |

## Methods
- **TOPPLE**: Regulon stability via leave-one-out perturbation → redistribution index
  - Psoriasis result: STAT5B = #1 stabilizer, RI 0.03→0.61
- **SICAI**: Coupling complexity = entropy of pairwise r_b profile per cell type
  - Psoriasis result: complexity–PASI ρ=0.480, P=0.020, n=23
- **DGSA**: Geometric decomposition of effect-size vectors across cell types

## Coding Standards
- Python 3.10+, pandas, numpy, scipy, matplotlib, seaborn, openpyxl
- Figures: 300 DPI, PDF + PNG, publication quality, Nature-style palette
- All scripts self-contained and reproducible
- Use `if __name__ == "__main__":` pattern

## Build
```bash
pip install pandas numpy scipy matplotlib seaborn openpyxl xlrd scikit-learn
```
