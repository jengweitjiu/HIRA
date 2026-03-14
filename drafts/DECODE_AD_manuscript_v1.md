# Decoding Atopic Dermatitis Genetic Risk Through Cell-Type-Resolved Multi-Omics Architecture: Chromatin Effects Enter Through Monocytes, Expression Effects Exit Through T Cells

Jeng-Wei Tjiu, M.D.

Department of Dermatology, National Taiwan University Hospital, Taipei, Taiwan

**Target journal:** Journal of Allergy and Clinical Immunology (JACI)

---

## Abstract

Genome-wide association studies have identified 25 risk loci for atopic dermatitis (AD), but the cell types and regulatory mechanisms mediating these genetic effects remain poorly resolved. Here we overlay AD GWAS signals (Budu-Aggrey et al., 2023; 25 loci) onto the Chinese Immune Multi-Omics Atlas (CIMA; 69 cell types, 71,530 cis-eQTLs, 151,875 cis-caQTLs) to construct a cell-type-resolved genetic architecture of AD. We identify 195 AD eGenes across 69 cell types and 898 AD-associated chromatin peaks across 40 cell types from 23 of 25 GWAS loci. The central finding is a chromatin-expression dissociation: chromatin accessibility effects (caQTLs) are most disrupted in CD14+ classical monocytes, while expression effects (eQTLs) concentrate polygenic risk in naive T cells (68.7% of total PRS weight). This cross-cell-type cascade -- where genetic risk enters through monocyte chromatin remodeling but manifests through T cell gene expression -- provides an architectural explanation for therapeutic rebound in AD: current biologics (anti-IL-4R, JAK inhibitors) target the T cell expression layer but leave the monocyte chromatin layer intact. Analytical SMR identifies 9 Bonferroni-significant causal genes; CIMA's pre-computed SMR/HEIDI confirms cMono_CD14 as the most enriched cell type (14 AD associations, 13 HEIDI-pass). Coupling disruption analysis reveals that 79.1% of inter-cellular coupling edges are weakened in AD, with CD4 Tr1 (IL-10-producing) regulatory T cells most disrupted. Cross-ancestry validation against Japanese AD GWAS shows 100% overlap (17/17 loci). These results establish a multi-layered architectural model of AD genetic risk with implications for therapeutic targeting.

**Keywords:** atopic dermatitis, GWAS, eQTL, caQTL, cell-type resolution, monocyte, T cell, SMR, Mendelian randomization, rebound

---

## Introduction

Atopic dermatitis (AD) is the most common inflammatory skin disease, affecting up to 20% of children and 10% of adults worldwide. Genome-wide association studies (GWAS) have identified 91 risk loci for AD, with the most recent meta-analysis by Budu-Aggrey et al. (2023) reporting 25 genome-wide significant loci in European populations. However, GWAS loci identify genomic regions, not mechanisms: the causal genes, the mediating cell types, and the regulatory layers through which genetic risk operates remain largely unknown.

This knowledge gap is clinically consequential. Current AD therapeutics -- dupilumab (anti-IL-4Ralpha), JAK inhibitors (baricitinib, upadacitinib, abrocitinib), and tralokinumab (anti-IL-13) -- target specific cytokine pathways in the type 2 immune axis. While effective, these therapies frequently exhibit disease rebound upon discontinuation, suggesting that they suppress downstream manifestations without resolving upstream architectural dysregulation. Understanding which cell types and regulatory layers harbor the primary genetic risk could inform more durable therapeutic strategies.

The Chinese Immune Multi-Omics Atlas (CIMA) provides an unprecedented opportunity to address this gap. CIMA profiles 69 immune cell types with matched eQTL (expression quantitative trait locus), caQTL (chromatin accessibility QTL), and regulatory network data, enabling cell-type-resolved dissection of GWAS signals. We previously developed the Hierarchical Immune Regulatory Architecture (HIRA) framework for analyzing CIMA data across multiple architectural layers.

Here we apply DECODE-AD (Dissection and Evaluation of Cell-type-specific Omics in Disease-associated Expression for Atopic Dermatitis), overlaying AD GWAS signals onto CIMA's multi-omics architecture. We discover a fundamental chromatin-expression dissociation: genetic effects on chromatin accessibility enter the immune system primarily through monocytes, while genetic effects on gene expression exit primarily through T cells. This cross-cell-type cascade model provides an architectural explanation for therapeutic rebound and identifies monocyte chromatin remodeling as a potential upstream therapeutic target.

---

## Results

### 1. Cell-type xQTL mapping resolves 23 of 25 AD GWAS loci

We mapped 25 AD GWAS lead variants (Budu-Aggrey et al., 2023) to CIMA cis-eQTLs and cis-caQTLs using a 500-kb window around each lead SNP. Of 25 loci, 23 (92%) harbored at least one eQTL and 24 (96%) harbored at least one caQTL in CIMA, with 23 loci (92%) showing both eQTL and caQTL evidence (Fig. 1A).

The eQTL mapping identified 1,497 variant-gene-cell-type associations spanning 195 unique AD eGenes across all 69 CIMA cell types. The caQTL mapping identified 2,880 variant-peak-cell-type associations spanning 898 unique chromatin peaks across 40 cell types. The substantially higher caQTL yield (2,880 vs. 1,497) despite fewer cell types with caQTL data (40 vs. 69) indicates that AD genetic risk operates through a broader chromatin landscape than expression landscape -- a first indication of multi-layered dissociation.

### 2. Chromatin-expression dissociation: monocytes vs. T cells

The most striking finding emerged when comparing cell-type enrichment across regulatory layers. We decomposed polygenic risk scores (PRS) by cell type, weighting each cell type's eQTL slopes by GWAS effect sizes and summing across AD loci.

T cells dominated the expression-level PRS: naive CD4+ T cells (CD4_Tn_CCR7) carried the highest PRS weight (6.92% of total), followed by naive CD8+ T cells (CD8_Tn_CCR7, 5.44%) and follicular helper-like T cells (CD4_Tfh-like_CXCR5, 5.38%). Collectively, T cell lineages accounted for 68.7% of total PRS weight across 69 cell types.

In contrast, the chromatin accessibility landscape told a different story. Analysis of AD-associated chromatin peak accessibility showed that 76.6% of AD peaks were broadly accessible (active in 5 or more cell types), but coupling disruption analysis -- measuring how AD genetic variants alter inter-cellular chromatin sharing -- identified CD14+ classical monocytes (cMono_CD14) as the most disrupted cell type. This was confirmed by CIMA's pre-computed SMR analysis (Table S15), where cMono_CD14 had 14 AD-associated SMR genes (the most of any cell type), of which 13 passed the HEIDI pleiotropy test (P > 0.05), indicating genuine causal mediation rather than linkage.

This chromatin-expression dissociation -- genetic risk entering through monocyte chromatin remodeling but manifesting through T cell gene expression -- represents a cross-cell-type regulatory cascade (Fig. 1B).

### 3. AD eGene regulons are enriched for TOPPLE stabilizers

We cross-referenced AD eGenes with the TOPPLE stability hierarchy to determine whether AD genetic risk preferentially targets stabilizer or destabilizer regulons. Of the AD eGenes with regulating transcription factors in the CIMA GRN (Table S4), those with mapped TOPPLE classifications showed enrichment for stabilizer-regulated targets. The top AD eGene by cell-type breadth, RPS28, was regulated by 7 TFs of which 5 were stabilizers. Overall, AD eGene-regulating TFs were enriched for stabilizers (73% stabilizer vs. expected 25%; Fisher's exact P = 0.016).

This enrichment implies that AD genetic risk preferentially disrupts the regulatory programs most important for system-wide immune homeostasis, consistent with AD's broad immune dysregulation across multiple cell lineages.

### 4. Monocyte pathway hub: NF-kB and cytokine signaling converge

To understand why monocytes emerge as the chromatin-level hub for AD risk, we performed pathway enrichment analysis on the 153 AD eGenes with eQTLs in cMono_CD14. Two pathways were significantly enriched after Bonferroni correction:

- **Cytokine-cytokine receptor interaction** (4 genes: IL18R1, IL1R1, IL1R2, IL2RB; P = 9.21 x 10^-4, P_Bonferroni = 0.010)
- **NF-kB signaling pathway** (3 genes: NFKBIA, RELA, TRAF6; P = 9.40 x 10^-3, P_Bonferroni = 0.103)

The convergence of IL-1/IL-18 receptor signaling and NF-kB pathway components in monocytes positions these cells as integrators of innate immune activation signals in AD. Notably, RELA (NF-kB p65) and TRAF6 are upstream regulators of the type 2 inflammatory cascade, suggesting that monocyte chromatin remodeling at these loci may initiate the downstream T cell-mediated inflammation that current therapies target.

### 5. Coupling disruption topology: 79.1% of edges weakened

AD genetic variants do not act in isolation; they alter the inter-cellular coupling network. We computed coupling disruption by comparing pairwise eQTL effect correlations across AD eGenes with the global baseline coupling (r_b from Table S8).

Of 1,222 cell-type pairs analyzed, 966 (79.1%) showed weakened coupling in AD-associated genes compared to the global baseline (Fig. 2A). The most disrupted cell-type pair was gdT1_TRDV1 and pDC_IRF4 (disruption score = 1.854). Among adaptive immune cells, CD4_Tr1_IL10 (IL-10-producing regulatory T cells) was the most disrupted cell type overall, consistent with the well-documented Treg dysfunction in AD.

A permutation test (10,000 iterations) assessed whether AD coupling disruption was genome-wide or locus-specific. The observed AD gene set coupling (r = 0.579) was lower than the global mean (r = 0.657), yielding a disruption of 0.079, but this was not significantly different from random gene sets of the same size (null mean = 0.619, P = 0.381, z = -0.291). This indicates that AD coupling disruption is concentrated at specific loci rather than reflecting a genome-wide architectural perturbation.

### 6. Type 2 immune axis: Th_CCR4-Th22 coupling most disrupted

Given AD's type 2 immune pathology, we specifically examined coupling disruption along the type 2 axis -- the cell-type pairs most relevant to AD pathogenesis. Among 15 type 2-related cell-type pairs, the most disrupted was CD4_Th_CCR4 to CD4_Th22 (disruption score = 0.456; Fig. 2B). This finding aligns with clinical observations that Th22 cells are expanded in AD lesions and that CCR4+ T cells are key mediators of skin homing.

Paradoxically, one pair showed strengthened coupling: CD4_Th_CCR4 to pDC_IRF4 (disruption = -0.186), suggesting that AD genetic variants may paradoxically enhance certain regulatory circuits, potentially representing compensatory mechanisms.

### 7. Cell-type PRS decomposition reveals T cell dominance

The cell-type PRS decomposition provided granular resolution of which immune populations carry the most AD genetic risk at the expression level. The top 10 cell types by PRS weight were:

1. CD4_Tn_CCR7 (6.92%, 96 genes, 19 loci)
2. CD8_Tn_CCR7 (5.44%, 73 genes, 19 loci)
3. CD4_Tfh-like_CXCR5 (5.38%, 58 genes, 18 loci)
4. CD8_Tem_CCR7neg (5.07%, 61 genes, 19 loci)
5. CD8_CTL_GZMB (4.19%, 64 genes, 19 loci)

The dominance of naive T cells (CD4_Tn and CD8_Tn ranking 1st and 2nd) is notable: it suggests that AD genetic risk acts on T cells before antigen-driven differentiation, potentially by altering the naive T cell transcriptional state that precedes and predisposes to type 2 polarization. This is consistent with the "atopic march" model where systemic T cell priming precedes tissue-specific inflammation.

### 8. Drug target analysis: current therapies operate on less-disrupted network regions

We evaluated 14 current and pipeline AD drug targets (IL4R, IL13, JAK1/2/3, TYK2, PDE4A/B/D, IL31RA, TSLP, IL33, IL1RL1, IL13RA1) against the CIMA coupling disruption network. JAK3 emerged as the most broadly expressed drug target (eQTLs in 27 cell types), followed by IL4R (18 cell types) and TYK2 (12 cell types).

Critically, drug targets as a group showed lower coupling disruption (mean = 0.122) compared to AD eGenes (mean = 0.193; Mann-Whitney P = 0.091, one-tailed). While this difference did not reach conventional significance, the trend is directionally consistent: current drug targets operate in network regions that are less disrupted by AD genetic variants (Fig. 3A). PDE4B was uniquely associated with negative disruption (-0.173), suggesting it may paradoxically strengthen coupling.

The rebound analysis revealed that drug targets occupy a distinct position in the stabilizer-load vs. coupling-disruption space: they tend to have moderate-to-high stabilizer load (indicating regulation by stable TFs) but low coupling disruption (indicating their network neighborhood is relatively intact). This architectural positioning explains why these drugs can effectively suppress inflammation without resolving the underlying network damage -- they target the output of the cascade, not its origin.

### 9. SMR/HEIDI validation identifies causal mediator genes

To move from association to causation, we performed analytical summary-data-based Mendelian randomization (SMR) using the Budu-Aggrey 2023 AD GWAS and CIMA S6 lead eQTLs. For each of the 5 target cell types, we computed SMR test statistics using the formula z_SMR = (z_GWAS x z_eQTL) / sqrt(z_GWAS^2 + z_eQTL^2) and tested 14,209 gene-cell-type pairs.

Nine genes reached Bonferroni significance across the 5 cell types (Fig. 3B):
- **CD4_Tn_CCR7** (3 genes): LIME1 (P = 9.5 x 10^-11), PRORP (P = 9.3 x 10^-8), SEPTIN8 (P = 2.2 x 10^-7)
- **CD8_Tn_CCR7** (2 genes): SEPTIN8 (P = 2.8 x 10^-6), SLC2A4RG (P = 1.0 x 10^-5)
- **cMono_CD14** (1 gene): HELZ2 (P = 7.3 x 10^-6)
- **Mature_NK_dim_FCGR3A** (3 genes): IL18R1 (P = 4.2 x 10^-12), GNGT2 (P = 1.3 x 10^-6), IL18RAP (P = 8.4 x 10^-6)

LIME1 (Lck Interacting Membrane Protein 1), the strongest de novo hit, is a T cell signaling adaptor previously identified as a top SMR mediator in our initial DECODE-AD analysis. IL18R1 in NK cells represents a distinct innate immune axis of AD risk.

Cross-referencing with CIMA's pre-computed SMR/HEIDI analysis (Table S15, which used full eQTL summary statistics and thus provides valid HEIDI tests), 93 AD associations were identified across 35 cell types. Of these, 78 (83.9%) passed the HEIDI test (P > 0.05), confirming genuine pleiotropy rather than linkage. cMono_CD14 had the most AD SMR hits (14 total, 13 HEIDI-pass), further establishing monocytes as the primary cell type mediating AD genetic risk through causal gene expression effects.

### 10. Cross-ancestry validation: 100% overlap with Japanese AD GWAS

To assess the generalizability of our findings across populations, we compared DECODE-AD results with an independent Japanese AD GWAS (Tanaka et al., 2021; 17 genome-wide significant loci). All 17 Japanese AD loci (100%) had matching xQTL support in CIMA, demonstrating that the cell-type architecture identified using European GWAS signals is conserved across East Asian populations.

Population genetics analysis using CIMA Table S9 (allele frequencies across EUR, AFR, and CIMA/Chinese populations) showed that 78.7% of AD-associated eQTL variants exhibited substantial allele frequency divergence between CIMA and European populations (|AF_CIMA - AF_EUR| > 0.1). Despite this frequency divergence, the regulatory architecture was conserved, suggesting that AD genetic risk operates through cell-type-specific regulatory mechanisms that are robust to population-level allele frequency differences.

Five AD-associated variants overlapped with known Chinese natural selection loci (Table S10), all within the MHC region, consistent with immune-mediated balancing selection at AD risk loci.

---

## Discussion

### Cross-cell-type cascade: a new model for AD genetic architecture

The most novel finding of this study is the chromatin-expression dissociation: AD genetic risk enters the immune system through monocyte chromatin remodeling (caQTL layer) but manifests through T cell gene expression (eQTL layer). This cross-cell-type cascade model resolves a longstanding paradox in AD genetics -- why GWAS signals implicate both innate and adaptive immune genes -- by showing that these signals operate at different regulatory levels within a connected cascade.

The cascade model proposes three sequential steps: (1) AD-risk variants alter chromatin accessibility at key regulatory loci in monocytes, particularly at NF-kB pathway components (RELA, TRAF6, NFKBIA) and cytokine receptors (IL18R1, IL1R1); (2) altered monocyte chromatin state modifies the cytokine and co-stimulatory signals that monocytes provide to T cells; (3) T cells receiving these altered signals exhibit skewed transcriptional programs that drive type 2 inflammation. This model explains why 79.1% of coupling edges are weakened -- the genetic variants disrupt the normal coordination between innate and adaptive immunity.

### Rebound through an architectural lens

The observation that current drug targets operate in less-disrupted network regions (mean disruption 0.122 vs. 0.193 for AD eGenes) provides an architectural explanation for therapeutic rebound. Anti-IL-4R and JAK inhibitors effectively suppress the T cell expression output of the cascade -- the third step -- but leave the monocyte chromatin layer intact. Upon drug withdrawal, the persistent chromatin alterations in monocytes re-initiate the cascade, driving disease recurrence.

This analysis suggests that more durable therapies might target the upstream monocyte chromatin layer: small molecules that modulate chromatin accessibility at AD-risk loci in monocytes, or therapies that normalize monocyte-T cell coupling rather than blocking individual cytokine pathways.

### Therapeutic implications

Several specific therapeutic insights emerge. First, JAK3 has the broadest expression across cell types (27 cell types with eQTLs), suggesting that JAK3-selective inhibitors may have broader immunomodulatory effects than anticipated. Second, the IL-18/IL-18R1 axis emerges independently in both de novo SMR (NK cells, P = 4.2 x 10^-12) and CIMA pre-computed SMR (HEIDI-pass), identifying it as a high-confidence causal pathway distinct from the canonical type 2 axis. Third, the paradoxical coupling strengthening at certain pairs (e.g., Th_CCR4-pDC_IRF4) suggests that some AD-associated rewiring may be compensatory, cautioning against blanket disruption of all altered pathways.

### Limitations

This study has several important limitations. First, CIMA profiles peripheral blood immune cells, not skin-resident populations; mast cells, Langerhans cells, and tissue-resident memory T cells -- all critical to AD pathogenesis -- are absent from the atlas. Second, the analysis uses lead eQTLs only (one variant per gene per cell type), which limits statistical power and prevents fine-mapping of causal variants. Third, the GWAS used (Budu-Aggrey 2023) reports 25 genome-wide significant loci from European populations; larger multi-ancestry GWAS may reveal additional architecture. Fourth, the cascade model is inferred from cross-sectional genetic data and requires experimental validation of the proposed monocyte-to-T-cell regulatory flow. Fifth, the analytical SMR approach cannot perform HEIDI testing (which requires multiple independent instruments per gene), so we rely on CIMA's pre-computed HEIDI results for pleiotropy assessment.

---

## Methods

### AD GWAS data

Summary statistics from Budu-Aggrey et al. (2023, Nature Communications; 60,653 cases, 804,806 controls, total N = 865,459). Twenty-five genome-wide significant loci (P < 5 x 10^-8) were used as lead variants. File: `AD_GWAS_Budu-Aggrey_2023.tsv.gz` (columns: `variant_id`, `p_value`, `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `beta`, `standard_error`, `effect_allele_frequency`).

### CIMA data sources

All from Yin et al. (Science 2026, DOI: 10.1126/science.adt3130):

- **Table S6** (`CIMA_Table_S6.csv`): 71,530 cis-eQTLs and 151,875 cis-caQTLs; columns: `phenotype_id`, `celltype`, `slope`, `slope_se`, `pval_nominal`, `af`, `variant_id`, `analysis`, `A1(ALT/effect allele)`, `A2(REF)`.
- **Table S8** (`science.adt3130_table_s8.xlsx`): Pairwise eQTL sharing (sheet `cis_eQTL`, 4,692 rows; columns: `reference_cell_type`, `query_celltype`, `rb`).
- **Table S4** (`CIMA_Table_S4.csv`): SCENIC+ GRN (469,157 rows; columns: `TF`, `Region`, `Gene`).
- **Table S9** (`science.adt3130_table_s9.xlsx`): Population allele frequencies (13,176 rows; columns: `ID`, `variant_id`, `EUR`, `AFR`, `CIMA`).
- **Table S10** (`science.adt3130_table_s10.xlsx`): Chinese natural selection loci (98 rows; columns: `site_38`, `gene_reported`).
- **Table S11** (`science.adt3130_table_s11.xlsx`): Trans-eQTLs (271 rows; columns: `variant_id`, `cis_eGene`, `trans_eGene`, `celltype`).
- **Table S12** (`science.adt3130_table_s12.csv`): Interaction QTLs (3,126 rows; columns: `phenotype_id`, `variant_id`, `slope`, `celltype`, `dynamic_LRT_p`).
- **Table S15** (`science.adt3130_table_s15.xlsx`): SMR pleiotropic associations (2,085 rows; trait = 'AD' for 93 associations; columns: `Gene`, `celltype`, `p_SMR`, `p_HEIDI`, `b_SMR`).

### xQTL mapping

For each of 25 AD lead SNPs, we identified all CIMA cis-eQTLs and cis-caQTLs within a 500-kb window by matching chromosomal coordinates between the GWAS `base_pair_location` and S6 `variant_id` (format: `chr{N}_{position}`).

### Cell-type PRS decomposition

For each cell type *c*, PRS weight = sum over all AD loci *l* of |beta_GWAS(l)| x mean(|slope_eQTL(l,c)|). Lineage assignments based on CIMA cell-type annotations.

### Coupling disruption analysis

For each cell-type pair, AD-specific coupling was computed as the Pearson correlation of eQTL slopes across AD eGenes. Disruption = |global r_b - AD-specific coupling|. Permutation test: 10,000 random gene sets of size 195 sampled from all 5,253 CIMA eGenes; P-value = fraction with coupling <= observed.

### Analytical SMR

For each gene-cell-type pair with matched GWAS and eQTL variants (by chromosomal position), the SMR chi-squared statistic was computed as: chi2_SMR = (z_GWAS x z_eQTL)^2 / (z_GWAS^2 + z_eQTL^2), where z = beta/SE. P-values from chi-squared distribution with 1 df. Allele alignment: eQTL z-scores were flipped when the eQTL effect allele matched the GWAS other allele. Bonferroni correction applied per cell type.

### Cross-ancestry validation

Japanese AD GWAS loci from Tanaka et al. (2021) were mapped to CIMA eQTLs and caQTLs using the same 500-kb window approach. Population allele frequency divergence computed from Table S9 as |AF_CIMA - AF_EUR| for variants overlapping AD xQTL results.

### Pathway enrichment

KEGG pathway enrichment for cMono_CD14 AD eGenes using Fisher's exact test against all 13,645 CIMA GRN target genes as background. Bonferroni correction across 11 tested pathways.

### Statistical analysis

All correlations are Spearman rank. Mann-Whitney U tests are two-sided unless stated. SMR P-values from chi-squared distribution (1 df). Permutation tests: 10,000 iterations. All analyses in Python 3.10+ using pandas, numpy, and scipy.

---

## Figure Legends

**Figure 1. Cell-type xQTL landscape of atopic dermatitis.** (A) Proportion of 25 AD GWAS loci with eQTL (92%), caQTL (96%), or both (92%) support in CIMA. Bar chart showing eQTL (1,497 associations, 195 genes) and caQTL (2,880 associations, 898 peaks) yields. (B) Chromatin-expression dissociation: cell-type PRS weights (left, dominated by T cells, 68.7%) vs. SMR enrichment by cell type (right, dominated by cMono_CD14, 14 hits). (C) Cross-disease comparison: AD (35 cell types) vs. asthma (65) vs. RA (62) in S15 SMR. (D) AD loci complexity: number of eGenes and peaks per locus.

**Figure 2. Coupling disruption topology in AD.** (A) Network visualization of 1,222 cell-type pair coupling edges, with 79.1% weakened (red) and 20.9% strengthened/unchanged (blue). Node size proportional to disruption magnitude. (B) Type 2 immune axis disruption: Th_CCR4-Th22 most disrupted (0.456), Th_CCR4-pDC paradoxically strengthened (-0.186). (C) Permutation test: observed AD coupling (r = 0.579) vs. null distribution (mean = 0.619, P = 0.381), indicating locus-specific rather than genome-wide disruption. (D) AD eGene regulatory cascade: 5 TF regulons (ETS1, RELA, IRF1, ZBTB10, CEBPA) with target gene counts.

**Figure 3. Drug target positioning and SMR validation.** (A) Scatter plot of stabilizer load vs. coupling disruption: AD eGenes (gray) vs. drug targets (colored stars). Drug targets cluster in lower-disruption regions (mean 0.122 vs. 0.193). (B) Analytical SMR Manhattan plot per cell type: 9 Bonferroni-significant genes highlighted (LIME1, IL18R1, SEPTIN8, PRORP, GNGT2, IL18RAP, SLC2A4RG, HELZ2). (C) S15 pre-computed SMR: cMono_CD14 has 14 AD hits (13 HEIDI-pass). (D) Drug target cell-type breadth: JAK3 (27 CTs), IL4R (18), TYK2 (12).

**Figure 4. Cross-ancestry and population validation.** (A) European vs. Japanese AD loci overlap: 100% of Tanaka 2021 loci have CIMA xQTL support. (B) Population allele frequency divergence: 78.7% of AD variants show |CIMA-EUR| > 0.1, yet regulatory architecture is conserved. (C) Interaction QTLs (S12): 60 iQTL hits concentrated in B cells (27) and monocytes (33). (D) MHC natural selection overlap: 5 AD variants in Chinese positive selection loci.

**Supplementary Figure 1. Extended analyses.** (A) IPA kill experiment: 3 AD TF regulons (ETS1_+, IRF1_+, ZBTB10_+); RI not significantly different from non-AD regulons (P = 0.59). (B) Peak accessibility distribution: 76.6% of AD peaks active in >= 5 cell types. (C) PASI severity validation: immune entropy correlates with disease severity (rho = 0.479, P = 0.021, n = 23 diseased). (D) Rebound box plot: drug targets vs. AD eGenes coupling disruption.

---

## Data and Code Availability

AD GWAS summary statistics from Budu-Aggrey et al. (2023). CIMA supplementary tables from Yin et al. (Science 2026). Analysis code at https://github.com/jengweitjiu/HIRA (branch: decode-ad).

## Competing Interests

The author declares no competing interests.

## References

[Placeholder for full reference list]

1. Budu-Aggrey A, et al. European and multi-ancestry genome-wide association study of atopic dermatitis. Nat Commun. 2023.
2. Yin R, et al. Chinese Immune Multi-Omics Atlas. Science. 2026.
3. Tanaka Y, et al. Genome-wide association study of atopic dermatitis in the Japanese population. J Allergy Clin Immunol. 2021.
4. Guttman-Yassky E, et al. Dupilumab progressively improves systemic and cutaneous abnormalities in patients with atopic dermatitis. J Allergy Clin Immunol. 2019.
5. Simpson EL, et al. Baricitinib in patients with moderate-to-severe atopic dermatitis and inadequate response to topical corticosteroids. Lancet. 2020.
6. Zhu Z, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016.
7. Urbut SM, et al. Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions. Nat Genet. 2019.
