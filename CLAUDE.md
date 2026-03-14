# CLAUDE.md — HIRA Project (Hierarchical Immune Regulatory Architecture)

## What is this project?
HIRA is a six-layer analytical framework that extracts architectural insights from
the Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026, DOI: 10.1126/science.adt3130).
Single author: Dr. Jeng-Wei Tjiu, M.D., Dept. of Dermatology, NTUH, Taipei, Taiwan.

## Key Constraint
- **NO raw data access** — all analyses use published supplementary tables ONLY
- Environment: Windows + Claude Code (primary) or Colab Pro+ (GPU tasks only)
- Python 3.10+, all scripts must be self-contained and runnable standalone

## Data Files (in data/raw/)
### Actual file names: science.adt3130_table_sN.xlsx (or .zip for S6)

| File | Sheet | Contents | Key Columns | Rows | Use |
|------|-------|----------|-------------|------|-----|
| table_s5.xlsx | `eRegulons_Activators_Exp_AUC_RS` | 203 eRegulons x 61 cell types | eRegulon, cell_type_l4, mean_AUC, RSS, Expression, Repressor_Activator | 12,383 | TOPPLE, IPA |
| table_s5.xlsx | `sex_difference_results` | Sex-difference data | cell_type_l4, eRegulon, Log2_Fold_Change, P_Value, adjust_P_Value | 21,289 | IPA |
| table_s5.xlsx | `Age_Correlation` | Age-correlation (wide format) | row=eRegulon, cols=CellType_value/CellType_Pvalue pairs | 404 | IPA |
| table_s6.zip -> CIMA_Table_S6.csv | analysis=`cis-eQTL` | Lead cis-eQTLs | phenotype_id, celltype, slope, slope_se, pval_nominal, af, analysis | 71,530 | DGSA |
| table_s6.zip -> CIMA_Table_S6.csv | analysis=`cis-caQTL` | Lead cis-caQTLs | phenotype_id, celltype, slope, slope_se, pval_nominal, af, analysis | 151,875 | **Extension #1** |
| table_s8.xlsx | `cis_eQTL` | Pairwise eQTL sharing | reference_cell_type, query_celltype, rb, pi1 | 4,692 | SICAI |
| table_s8.xlsx | `cis_caQTL` | Pairwise caQTL sharing | reference_cell_type, query_celltype, rb, pi1 | 1,722 | **Extension #1** |
| table_s15.xlsx | `Sheet1` | SMR pleiotropic associations | Gene, celltype, trait, trait_category, p_SMR, p_HEIDI, b_SMR | 2,085 | Disease correlation |
| table_s3.zip -> CIMA_Table_S3.csv | — | Peak x cell-type binary activity | Peaks (index), 65 cell-type columns (TRUE/FALSE) | 338,036 | **Extension #3** |
| table_s4.zip -> CIMA_Table_S4.csv | — | SCENIC+ GRN TF->Region->Gene | TF, Region, Gene, R2G_importance, R2G_rho, Consensus_name | 469,157 | **Extension #3** |

**NOTE:** S5 Age_Correlation is wide-format (each column is CellType_eRegulon + CellType_eRegulon_Pvalue).
S6 is a single CSV inside a zip; filter by `analysis` column for eQTL vs caQTL.
S8 caQTL has only 1,722 rows (42 cell types vs 69 for eQTL).
S4 has 133,574 unique Region-Gene pairs (84,625 regions -> 13,645 genes).
S3 has 65 cell types (vs S5's 61 for eRegulons and S8's 69 for eQTL).

## Six HIRA Layers
1. **TOPPLE** (src/topple.py): Regulon stability via leave-one-out JSD perturbation
   - Input: S5 AUC matrix (203 regulons x 61 cell types)
   - Output: Redistribution index (RI) per regulon, stability ranking
2. **DGSA** (src/dgsa.py): eQTL effect-size geometry decomposition
   - Input: S6 cis-eQTL (genes in >=3 cell types -> 5,253 genes)
   - Output: Non-additivity, magnitude, uniformity, sparsity, Gini per gene
3. **SICAI** (src/sicai.py): Inter-cellular coupling topology via r_b
   - Input: S8 pairwise r_b (69 x 69 matrix)
   - Output: Mean r_b, CV, entropy per cell type
4. **IPA** (src/ipa.py + src/ipa_age.py): Perturbation resistance (sex + age)
   - Input: S5 sheets (sex log2FC, age Spearman rho)
   - Output: Stabilizer vs destabilizer comparison (Mann-Whitney)
5. **STRATA** (src/strata.py): Spatial validation on psoriasis Visium GSE202011
6. **STRATA-Atlas** (src/strata_atlas.py): Cross-tissue Mantel test (blood vs skin coupling)

## Extensions
### Extension #1: caQTL Geometric Architecture (src/ext1_caqtl_dgsa.py) -- DONE
- caQTL mean non-additivity = 0.871 (> eQTL 0.864)
- Paired Wilcoxon at matched loci: P = 7e-35 (caQTL > eQTL)
- Mantel eQTL vs caQTL r_b: r = 0.482, P = 0.0001

### Extension #3: Enhancer-Driven Coupling Network (src/ext3_enhancer_coupling.py) -- DONE
- Source: S4 (SCENIC+ GRN, 469K rows) + S3 (338K peaks x 65 cell types binary)
- 133,574 unique Region-Gene pairs, mean ~43K active per cell type
- Jaccard coupling matrix: 65 x 65, mean Jaccard = 0.406
- Three-way Mantel: eQTL rb vs caQTL rb (0.482), eQTL vs Jaccard (0.388), caQTL vs Jaccard (0.566)
- Stabilizers have 2.4x more target genes (417 vs 175, P=7.9e-4)
- Regulatory breadth vs TOPPLE RI: rho=0.283 (P=4.3e-5)

## Expected Key Numbers (for validation)
| Metric | Expected Value | Source |
|--------|---------------|--------|
| TOPPLE: HSF1 RI | ~0.0063 (rank #1) | HIRA manuscript |
| TOPPLE: STAT5B RI | ~0.0031 (rank #18) | HIRA manuscript |
| DGSA-eQTL: mean non-additivity | ~0.867 | HIRA manuscript |
| DGSA-eQTL: rho with disease | ~0.690 (P = 7.4e-11) | HIRA manuscript |
| SICAI: mean r_b | ~0.818 | HIRA manuscript, matches CIMA published 0.82 |
| SICAI: rho with disease | ~0.353 (P = 0.003) | HIRA manuscript |
| IPA: sex perturbation | Mann-Whitney P ~ 4.4e-218 | HIRA manuscript |
| Extension #1: caQTL non-additivity | 0.871 (mean); paired Wilcoxon P=7e-35 | confirmed |
| Extension #1: Mantel eQTL vs caQTL rb | r=0.482, P=0.0001 | confirmed |
| Extension #3: Mantel eQTL vs Jaccard | r=0.388, P<0.0001 | confirmed |
| Extension #3: Mantel caQTL vs Jaccard | r=0.566, P<0.0001 | confirmed |
| Extension #3: Breadth vs RI | rho=0.283, P=4.3e-5 | confirmed |
| STRATA: lesional vs healthy ratio | 3.42 vs 2.88, P=1.1e-3 | confirmed |
| STRATA: stabilizer expr les vs heal | 0.257 vs 0.078, P=7.7e-4 | confirmed |
| STRATA-Atlas: Mantel blood vs skin | r=0.143, P=0.0001 (59 shared CTs) | confirmed |
| IPA-age: stab vs dest |age_rho| | 0.139 vs 0.086, P=1.3e-112 | confirmed |
| IPA-age: RI vs |age_rho| | rho=0.435, P=8.8e-11 | confirmed |

## Coding Conventions
- Python 3.10+, pandas, numpy, scipy, matplotlib, seaborn, openpyxl
- All scripts runnable standalone: `python src/xxx.py`
- Results -> results/ as .csv
- Figures -> figures/ as .pdf + .png (300 DPI, Nature-style, Arial font, 180mm wide)
- Use argparse for configurable parameters
- Print key statistics to stdout for quick validation
- Every script must have `if __name__ == "__main__":` block
- **ALWAYS print column names** when first reading any Excel file

## Git Conventions
- Remote: github.com/jengweitjiu/HIRA
- Branch: main (stable), dev (work-in-progress)
- Commit messages: `feat:`, `fix:`, `data:`, `fig:` prefixes
- Run pytest tests/ before every commit

## Build
```bash
pip install pandas numpy scipy matplotlib seaborn openpyxl xlrd scikit-learn scikit-bio
```

## Project Structure
```
D:\HIRA\
├── CLAUDE.md              <- this file
├── data\raw\              <- CIMA xlsx files go here
├── data\visium\           <- GSE202011 Visium .h5 files
├── data\processed\        <- intermediate outputs
├── src\                   <- all Python scripts
├── tests\                 <- pytest test files
├── figures\               <- publication figures
├── results\               <- csv result tables
├── drafts\                <- manuscript results paragraphs
└── requirements.txt
```
