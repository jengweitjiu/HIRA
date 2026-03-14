# CLAUDE.md — HIRA Project (Hierarchical Immune Regulatory Architecture)

## What is this project?
HIRA is a six-layer analytical framework that extracts architectural insights from
the Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026, DOI: 10.1126/science.adt3130).
Single author: Dr. Jeng-Wei Tjiu, M.D., Dept. of Dermatology, NTUH, Taipei, Taiwan.

## Key Constraint
- **NO raw data access** — all analyses use published supplementary tables ONLY
- Environment: Windows + Claude Code (primary) or Colab Pro+ (GPU tasks only)
- Python 3.10+, all scripts must be self-contained and runnable standalone

## Git Branch Status
- **main**: Core HIRA pipeline complete (TOPPLE, DGSA, SICAI, IPA, STRATA, Extensions)
  - Latest: `f3b86c7` — manuscript-ready HIRA pipeline
- **decode-ad** (active): DECODE-AD extension (AD GWAS overlay on CIMA)
  - Latest: `5d10a9c` — PASI severity validation
  - Contains all main branch content plus 15+ DECODE-AD analyses
- Remote: github.com/jengweitjiu/HIRA

---

## Data Files (in data/raw/)

### Verified Column Names

| File | Format | Sheet/Filter | Rows | Key Columns |
|------|--------|-------------|------|-------------|
| `science.adt3130_table_s3.zip` -> `CIMA_Table_S3.csv` | CSV | — | 338,036 | Peaks (index), 65 cell-type columns (TRUE/FALSE) |
| `science.adt3130_table_s4.zip` -> `CIMA_Table_S4.csv` | CSV | — | 469,157 | `TF`, `Region`, `Gene`, `R2G_importance`, `R2G_rho`, `Consensus_name` |
| `science.adt3130_table_s5.xlsx` | XLSX | `eRegulons_Activators_Exp_AUC_RS` | 12,383 | `eRegulon`, `cell_type_l4`, `mean_AUC`, `RSS`, `Expression`, `Repressor_Activator` |
| `science.adt3130_table_s5.xlsx` | XLSX | `sex_difference_results` | 21,289 | `cell_type_l4`, `eRegulon`, `Log2_Fold_Change`, `P_Value`, `adjust_P_Value` |
| `science.adt3130_table_s5.xlsx` | XLSX | `Age_Correlation` | 404 | Wide format: row=eRegulon, cols=`CellType_value`/`CellType_Pvalue` pairs |
| `science.adt3130_table_s6.zip` -> `CIMA_Table_S6.csv` | CSV | `analysis=cis-eQTL` | 71,530 | `phenotype_id`, `celltype`, `slope`, `slope_se`, `pval_nominal`, `af`, `analysis`, `variant_id` |
| `science.adt3130_table_s6.zip` -> `CIMA_Table_S6.csv` | CSV | `analysis=cis-caQTL` | 151,875 | Same as eQTL columns |
| `science.adt3130_table_s8.xlsx` | XLSX | `cis_eQTL` | 4,692 | `reference_cell_type`, `query_celltype`, `rb`, `pi1` |
| `science.adt3130_table_s8.xlsx` | XLSX | `cis_caQTL` | 1,722 | `reference_cell_type`, `query_celltype`, `rb`, `pi1` |
| `science.adt3130_table_s9.xlsx` | XLSX | — | 13,176 | `ID`, `REF`, `ALT`, `EUR`, `AFR`, `TOT`, `pair`, `variant_id`, `CIMA` |
| `science.adt3130_table_s10.xlsx` | XLSX | — | 98 | `study_pmid`, `study`, `original_site`, `site_38`, `gene_reported` |
| `science.adt3130_table_s11.xlsx` | XLSX | — | 271 | `variant_id`, `pval`, `b`, `b_se`, `af`, `variant_chr`, `phenotype_chr`, `celltype`, `padj`, `cis_eGene`, `trans_eGene`, `A1(ALT/effect allele)`, `A2(REF)` |
| `science.adt3130_table_s12.csv` | CSV | — | 3,126 | `phenotype_id`, `variant_id`, `slope`, `slope_se`, `pval_nominal`, `celltype`, `af`, `dynamic_LRT_p`, `dynamic_LRT_q` |
| `science.adt3130_table_s15.xlsx` | XLSX | `Sheet1` | 2,085 | `Gene`, `celltype`, `trait`, `trait_category`, `p_SMR`, `p_HEIDI`, `b_SMR` |
| `AD_GWAS_Budu-Aggrey_2023.tsv.gz` | TSV.GZ | — | ~8M | `rsid`, `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `beta`, `standard_error`, `p_value`, `effect_allele_frequency` |

**Notes:**
- S5 Age_Correlation is wide-format (each column is CellType_eRegulon + CellType_eRegulon_Pvalue)
- S6 is a single CSV inside a zip; filter by `analysis` column for eQTL vs caQTL
- S8 caQTL has only 1,722 rows (42 cell types vs 69 for eQTL)
- S9 contains eQTL/caQTL variant SNPs (NOT GWAS lead SNPs); match via variant_id
- S10 `site_38` column is in chr_pos format (e.g., `chr1_12345678`)
- S12 `variant_id` is in chr_pos format; 3,126 interaction QTLs

---

## HIRA Results — Actual Measured Values

### Layer 1: TOPPLE (src/topple.py)
- 203 regulons x 61 cell types, leave-one-out JSD perturbation
- 51 stabilizers, 51 destabilizers, 101 neutral
- HSF1_+ RI = 0.00632 (rank #1), EGR1_+ = 0.00574 (#2), KLF9_+ = 0.00569 (#3)
- Output: `results/topple.csv`

### Layer 2: DGSA (src/dgsa.py)
- 5,253 genes in >=3 cell types
- Mean non-additivity = 0.8673, mean Gini = 0.1151
- Non-additivity vs disease (S15 n_traits): rho = 0.690, P = 7.4e-11
- Output: `results/dgsa_gene_scores.csv`, `results/dgsa_celltype_summary.csv`

### Layer 3: SICAI (src/sicai.py)
- 69 x 69 eQTL r_b matrix, mean r_b = 0.8182
- Coupling complexity vs disease burden: rho = 0.353, P = 0.003
- Output: `results/sicai_coupling_complexity.csv`, `results/sicai_rb_matrix.csv`

### Layer 4: IPA — Sex (src/ipa.py)
- Stabilizers have LOWER sex perturbation: median |log2FC| = 0.018 vs destab 0.057
- Mann-Whitney P = 4.4e-218
- Output: `results/ipa_sex.csv`

### Layer 4: IPA — Age (src/ipa_age.py)
- Stabilizers have HIGHER age correlation: mean |rho| = 0.139 vs destab 0.086
- Mann-Whitney P = 1.3e-112 (opposite direction: stabilizers MORE age-sensitive)
- RI vs |age_rho|: rho = 0.435, P = 8.8e-11
- Output: `results/ipa_age.csv`

### Layer 5: STRATA (src/strata.py)
- GSE202011 Visium psoriasis, 30 samples (14 lesional, 9 non-lesional, 7 healthy)
- Stabilizer expression: lesional 0.257 vs healthy 0.078, P = 7.7e-4
- Stab/Destab ratio: lesional 3.42 vs healthy 2.88, P = 1.1e-3
- Output: `results/strata_summary.csv`, `results/strata_sample_metrics.csv`

### Layer 6: STRATA-Atlas (src/strata_atlas.py)
- Blood vs skin coupling Mantel test: r = 0.143, P = 0.0001 (59 shared cell types)
- Output: `results/strata_atlas.csv`

### Extension #1: caQTL Geometric Architecture (src/ext1_caqtl_dgsa.py, src/ext1_remaining.py)
- caQTL mean non-additivity = 0.871 (> eQTL 0.864) — "epigenomic amplification"
- Paired Wilcoxon at matched loci: P = 7e-35 (caQTL > eQTL)
- Mantel eQTL vs caQTL r_b: r = 0.482, P = 0.0001
- Per-cell-type eQTL vs caQTL Gini: rho = 0.746, P = 1.5e-8 (42 cell types)
- caQTL mean r_b = 0.788 vs eQTL 0.825, P = 3.9e-7 (caQTL less coupled)
- caQTL mean r_b best disease predictor: rho = 0.511, P = 5.5e-4
- Output: `results/ext1_caqtl_dgsa.csv`, `results/ext1_celltype_correlation.csv`, `results/ext1_sicai_caqtl_metrics.csv`

### Extension #3: Enhancer-Driven Coupling Network (src/ext3_enhancer_coupling.py, src/ext3_ad_overlay.py)
- S4 GRN: 469K rows, 133,574 unique Region-Gene pairs (84,625 regions -> 13,645 genes)
- Jaccard coupling matrix: 65 x 65, mean Jaccard = 0.406
- Three-way Mantel: eQTL vs caQTL (0.482), eQTL vs Jaccard (0.388), caQTL vs Jaccard (0.566)
- Stabilizers have 2.4x more target genes (417 vs 175, P = 7.9e-4)
- Regulatory breadth vs TOPPLE RI: rho = 0.283, P = 4.3e-5
- AD overlay: 2,661 Region-Gene pairs in AD loci (1,510 regions, 277 genes, 196 TFs)
  - Top TFs: ELF1 (144 regions), ETS1 (114), FLI1 (92)
- Output: `results/ext3_regulatory_coupling.csv`, `results/ext3_ad_peak_gene_overlay.csv`

### PASI Severity Validation (src/decode_pasi_validation.py)
- GSE202011 Visium, 30 samples, 9 patients (6 psoriasis PASI 1.8-32.0, 3 healthy PASI=0)
- **Significant correlations (P<0.05):**

| Layer | Metric | Subset | rho | P |
|-------|--------|--------|-----|---|
| TOPPLE | Stabilizer expression | all (n=30) | 0.446 | 1.3e-2 |
| TOPPLE | Destabilizer expression | all (n=30) | 0.373 | 4.2e-2 |
| TOPPLE | Stab/Destab ratio | all (n=30) | 0.404 | 2.7e-2 |
| DGSA proxy | Immune entropy | all (n=30) | 0.384 | 3.6e-2 |
| DGSA proxy | Immune entropy | diseased (n=23) | **0.479** | **2.1e-2** |

- SICAI coupling strength vs PASI: rho = -0.148, P = 0.43 (NOT significant)
- Key insight: Immune entropy is the only metric significant in diseased-only subset
- Output: `results/pasi_hira_correlations.csv`, `figures/fig_pasi_validation.pdf`

---

## DECODE-AD Results (decode-ad branch)

### Core DECODE-AD Pipeline (src/decode_ad.py)
- AD GWAS: 25 loci (Budu-Aggrey 2023), 500kb window matching
- 1,497 eQTL hits (195 genes, 69 cell types), 2,880 caQTL hits (898 peaks, 40 cell types)
- Output: `results/decode_ad_eqtl.csv`, `results/decode_ad_caqtl.csv`, `results/decode_ad_summary.csv`

### Analysis Results Table

| # | Analysis | Script | Key Result | Output |
|---|----------|--------|------------|--------|
| 1 | Disease comparison | `decode_ad_extended.py` | AD 93 rows/35 CTs, As 245/65, RA 298/62 | `decode_ad_disease_comparison.csv` |
| 2 | Drug targets | `decode_ad_extended.py` | IFNG sole S15 match; JAK3 highest (27 CTs) | `decode_ad_drug_targets.csv` |
| 3 | SMR mediation | `decode_ad_extended.py` | 258 genes; KDELR2 (11 CTs), LIME1 (9 CTs) top mediators | `decode_ad_smr_mediation.csv` |
| 4 | Population genetics | `decode_ad_popgen.py` | 78.7% AF divergent (|CIMA-EUR|>0.1); 5 MHC selection overlaps | `decode_ad_popgen.csv` |
| 5 | Peak accessibility | `decode_ad_popgen.py` | 76.6% AD peaks are broad (>=5 CTs) | `decode_ad_peak_accessibility.csv` |
| 6 | Coupling disruption | `decode_ad_popgen.py` | 79.1% edges weakened; CD4_Tr1_IL10 most disrupted | `decode_ad_disruption_network.csv` |
| 7 | Cell-type PRS | `decode_ad_final.py` | T cells 68.7% PRS; CD4_Tn_CCR7 top (6.9%) | `decode_ad_celltype_prs.csv` |
| 8 | Regulatory cascades | `decode_ad_final.py` | 5 AD TF regulons; ETS1 2135 targets, RELA 1466, IRF1 998 | `decode_ad_cascades.csv` |
| 9 | Network visualization | `decode_ad_final.py` | — | `fig_decode_ad_network.pdf` |
| 10 | Drug target deep | `decode_ad_final.py` | JAK3 27 CTs, IL4R 18 CTs, TYK2 12 CTs | `decode_ad_drug_targets_deep.csv` |
| 11 | Interaction QTLs (S12) | `decode_ad_final.py` | 60 iQTL hits (Monocyte + B cells) | `decode_ad_interaction_qtl.csv` |
| 12 | Tanaka sensitivity | `decode_ad_final.py` | 100% overlap (17/17 GWAS loci with xQTL) | `decode_ad_tanaka_comparison.csv` |
| 13 | IPA kill experiment | `decode_ad_gaps.py` | 3 AD TF regulons (ETS1, IRF1, ZBTB10); RI not different P=0.59 | `decode_ad_ipa_kill.csv` |
| 14 | Permutation test | `decode_ad_gaps.py` | AD coupling r=0.579 vs null 0.619, P=0.38 — disruption is locus-specific | `decode_ad_permutation.csv` |
| 15 | Type 2 axis | `decode_ad_gaps.py` | CD4_Th_CCR4<->CD4_Th22 most disrupted (0.456) | `decode_ad_type2_axis.csv` |
| 16 | Rebound analysis | `decode_ad_rebound.py` | Drug targets mean disruption 0.122 vs AD eGenes 0.193, P=0.091 | `decode_ad_rebound_analysis.csv` |
| 17 | PASI validation | `decode_pasi_validation.py` | Immune entropy rho=0.479 P=0.021 (diseased only) | `pasi_hira_correlations.csv` |

### Most Novel Finding: Chromatin-Expression Dissociation
- **Chromatin accessibility (caQTL)**: Most disrupted in CD14_Mono (monocytes)
- **Gene expression (eQTL)**: Most PRS weight in T cells (68.7%)
- This dissociation explains **rebound** in AD: drugs targeting expression layer (T cell cytokines like IL-4R/JAK) miss the chromatin layer (monocyte epigenomic rewiring), which persists and drives relapse
- Drug targets have lower coupling disruption (0.122 vs 0.193) but P=0.091 — they operate in less disrupted network regions
- PDE4B uniquely shows negative disruption (-0.173) — may paradoxically strengthen coupling

---

## Remaining Items

### Not Yet Done
1. **SMR/HEIDI filtering**: Current SMR mediation uses raw p_SMR; should filter p_HEIDI > 0.05 for pleiotropic (not linkage) associations
2. **Thompson et al. 2024 validation**: Independent AD eQTL dataset for replication of DECODE-AD cell-type enrichment
3. **Manuscript submission**: Drafts exist in `drafts/` but need final formatting for journal submission
4. **Merge decode-ad -> main**: After review, merge the decode-ad branch into main

### Previously Investigated but Not Reproducible
- The "n=23, rho=0.480" SICAI-PASI result was hardcoded from earlier Colab work (`02_sicai_cima.py`); the source cohort/metric is unknown and does not match current GSE202011 data (n=30, rho=-0.148)

---

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
- Branches: main (stable), decode-ad (DECODE-AD extension)
- Commit messages: `feat:`, `fix:`, `data:`, `fig:` prefixes

## Build
```bash
pip install pandas numpy scipy matplotlib seaborn openpyxl xlrd scikit-learn scikit-bio networkx
```

---

## Project Structure
```
D:\HIRA\
├── CLAUDE.md                          <- this file
├── .claude\commands\                  <- custom slash commands (7 total)
│   ├── figures.md                     <- /figures — regenerate all figures
│   ├── run-topple.md                  <- /run-topple
│   ├── run-dgsa.md                    <- /run-dgsa
│   ├── run-sicai.md                   <- /run-sicai
│   ├── run-ext1.md                    <- /run-ext1
│   ├── run-ext3.md                    <- /run-ext3
│   └── validate.md                    <- /validate — run all validation tests
├── data\
│   ├── raw\                           <- CIMA xlsx/csv/zip + AD GWAS tsv.gz
│   ├── visium\                        <- GSE202011 Visium .h5 files
│   └── processed\                     <- intermediate outputs
├── src\                               <- all Python scripts (29 files)
│   ├── topple.py                      <- Layer 1: TOPPLE
│   ├── dgsa.py                        <- Layer 2: DGSA
│   ├── sicai.py                       <- Layer 3: SICAI
│   ├── ipa.py                         <- Layer 4: IPA sex
│   ├── ipa_age.py                     <- Layer 4: IPA age
│   ├── strata.py                      <- Layer 5: STRATA
│   ├── strata_atlas.py                <- Layer 6: STRATA-Atlas
│   ├── ext1_caqtl_dgsa.py             <- Extension #1: caQTL DGSA
│   ├── ext1_remaining.py              <- Extension #1: cell-type corr + SICAI caQTL
│   ├── ext1_deep_comparison.py        <- Extension #1: deep eQTL vs caQTL
│   ├── ext3_enhancer_coupling.py      <- Extension #3: enhancer coupling
│   ├── ext3_ad_overlay.py             <- Extension #3: AD GWAS overlay on S4
│   ├── decode_ad.py                   <- DECODE-AD core pipeline
│   ├── decode_ad_deep.py              <- DECODE-AD deep analysis
│   ├── decode_ad_extended.py          <- DECODE-AD disease comparison + drugs + SMR
│   ├── decode_ad_popgen.py            <- DECODE-AD population genetics
│   ├── decode_ad_final.py             <- DECODE-AD 6 final analyses
│   ├── decode_ad_gaps.py              <- DECODE-AD IPA kill + permutation + Type 2
│   ├── decode_ad_rebound.py           <- DECODE-AD rebound analysis
│   ├── decode_pasi_validation.py      <- PASI severity validation
│   ├── utils.py                       <- Shared utility functions
│   ├── fig_decode_ad.py               <- DECODE-AD main figure
│   ├── fig_decode_ad_extended.py      <- DECODE-AD extended figure
│   ├── fig_ext1_epigenomic_amplification.py
│   ├── fig_ext3_coupling.py
│   ├── fig_strata.py
│   ├── fig_study_design.py            <- Study design + circos figures
│   ├── figure1_complete.py            <- Complete Figure 1
│   └── sicai_caqtl.py                 <- SICAI caQTL analysis
├── results\                           <- 55 CSV files
├── figures\                           <- 68 PDF/PNG figure files
├── drafts\                            <- Manuscript text (8 .md files)
│   ├── HIRA_manuscript_v2_complete.md
│   ├── topple_results.md
│   ├── dgsa_results.md
│   ├── sicai_results.md
│   ├── ipa_results.md
│   ├── strata_results_text.md
│   ├── ext1_results_text.md
│   └── ext3_results_text.md
├── tests\                             <- pytest test files
└── requirements.txt
```

---

## Key Numbers for Validation (Actual Measured Values)

| Metric | Value | Status |
|--------|-------|--------|
| TOPPLE: HSF1_+ RI | 0.00632 (rank #1) | confirmed |
| TOPPLE: regulon counts | 51 stab / 51 destab / 101 neutral | confirmed |
| DGSA: mean non-additivity | 0.8673 (5,253 genes) | confirmed |
| DGSA: non-additivity vs disease | rho=0.690, P=7.4e-11 | confirmed |
| SICAI: mean r_b | 0.8182 (69 cell types) | confirmed |
| SICAI: coupling vs disease | rho=0.353, P=0.003 | confirmed |
| IPA sex: stab vs destab |log2FC| | 0.018 vs 0.057, P=4.4e-218 | confirmed |
| IPA age: stab vs destab |rho| | 0.139 vs 0.086, P=1.3e-112 | confirmed |
| STRATA: lesional vs healthy ratio | 3.42 vs 2.88, P=1.1e-3 | confirmed |
| STRATA-Atlas: Mantel blood vs skin | r=0.143, P=0.0001 | confirmed |
| Ext #1: caQTL non-additivity | 0.871 (> eQTL 0.864), Wilcoxon P=7e-35 | confirmed |
| Ext #1: caQTL mean r_b | 0.788 (< eQTL 0.825), P=3.9e-7 | confirmed |
| Ext #1: Mantel eQTL vs caQTL rb | r=0.482, P=0.0001 | confirmed |
| Ext #3: Mantel caQTL vs Jaccard | r=0.566, P<0.0001 | confirmed |
| Ext #3: stab vs destab target genes | 417 vs 175, P=7.9e-4 | confirmed |
| PASI: immune entropy (diseased) | rho=0.479, P=0.021, n=23 | confirmed |
| PASI: SICAI coupling vs PASI | rho=-0.148, P=0.43 (NS) | confirmed |
| DECODE-AD: eQTL hits | 1,497 (195 genes, 69 CTs) | confirmed |
| DECODE-AD: caQTL hits | 2,880 (898 peaks, 40 CTs) | confirmed |
| DECODE-AD: T cell PRS fraction | 68.7% | confirmed |
| DECODE-AD: Tanaka overlap | 100% (17/17 loci) | confirmed |
| DECODE-AD: permutation P | 0.38 (disruption is locus-specific) | confirmed |
| DECODE-AD: drug disruption vs AD | 0.122 vs 0.193, P=0.091 | confirmed |
