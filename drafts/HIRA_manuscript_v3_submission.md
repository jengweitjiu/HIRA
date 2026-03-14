# Hierarchical Immune Regulatory Architecture: Multi-Scale Geometric Analysis Reveals Stability Landscapes, Epigenomic Amplification, and Multi-Layered Coupling Topology of the Human Immune System

Jeng-Wei Tjiu, M.D.

Department of Dermatology, National Taiwan University Hospital, Taipei, Taiwan

**Target journal:** Genome Biology

---

## Abstract

Large-scale single-cell atlases have catalogued the transcriptional and epigenomic diversity of the human immune system, yet the architectural principles governing how hundreds of regulatory programs maintain system-wide homeostasis remain uncharacterized. Here we introduce the Hierarchical Immune Regulatory Architecture (HIRA), a six-layer analytical framework applied to the Chinese Immune Multi-Omics Atlas (CIMA; 203 eRegulons, 71,530 cis-eQTLs, 151,875 cis-caQTLs across up to 69 cell types). Layer 1 (TOPPLE) identifies 51 stabilizer regulons whose removal maximally redistributes system activity; HSF1 ranks first (redistribution index = 0.00632). Layer 2 (DGSA) decomposes eQTL effect-size geometry across 5,253 genes, revealing near-universal non-additivity (mean = 0.864) that predicts disease pleiotropy (rho = 0.796, P = 5.2 x 10^-16). Layer 3 (SICAI) maps a 69 x 69 inter-cellular coupling topology (mean r_b = 0.819) where coupling complexity predicts disease burden (rho = 0.509, P = 9.4 x 10^-6). Layer 4 (IPA) uncovers a striking dissociation: stabilizers resist sex-based perturbation (P = 1.96 x 10^-243) yet are preferentially sensitive to aging (P = 1.27 x 10^-112). Two extensions demonstrate that chromatin-level non-additivity exceeds expression-level non-additivity at matched loci (paired Wilcoxon P = 7.0 x 10^-35; "epigenomic amplification") and that three coupling modalities -- genetic, epigenomic, and regulatory -- form a convergent hierarchy (all Mantel P < 0.0001). Spatial validation in psoriasis Visium data confirms that stabilizer-to-destabilizer ratios scale with disease severity (PASI rho = 0.479, P = 0.021). HIRA establishes that immune homeostasis is maintained through hierarchically organized stability, coupling, and perturbation-resistance landscapes.

**Keywords:** systems immunology, regulatory architecture, eQTL geometry, inter-cellular coupling, epigenomic amplification, single-cell atlas

---

## Introduction

The past decade has witnessed an explosion of single-cell and multi-omics atlases profiling the human immune system at unprecedented resolution. The Human Cell Atlas, Tabula Sapiens, and most recently the Chinese Immune Multi-Omics Atlas (CIMA) have catalogued hundreds of cell types, thousands of regulatory programs, and millions of genetic associations across immune lineages. Yet a fundamental gap remains: while these atlases describe *what* exists, they have not addressed *how* the system is architecturally organized to maintain homeostasis.

This gap matters because immune dysregulation underlies autoimmune diseases, chronic inflammation, and immunodeficiency. Understanding which regulatory programs are structurally essential versus dispensable, how genetic effects propagate across cell types, and what topological features predict disease susceptibility requires moving beyond cataloguing to architecture.

Several conceptual precedents exist. Network biology has long emphasized that biological systems exhibit scale-free and hierarchical organization. Perturbation biology has shown that system responses to knockouts reveal functional hierarchy. Information-theoretic approaches have quantified cellular diversity and communication. However, no framework has integrated these perspectives across multiple data modalities -- transcription factor activity, genetic effect sizes, inter-cellular sharing, and chromatin accessibility -- within a single atlas.

Here we introduce HIRA (Hierarchical Immune Regulatory Architecture), a six-layer analytical framework that extracts architectural principles from the CIMA dataset without requiring raw data access. HIRA operates entirely on published summary statistics: eRegulon activity scores (Table S5), cis-eQTL and cis-caQTL effect sizes (Table S6), pairwise genetic sharing coefficients (Table S8), and regulatory network linkages (Tables S3-S4). Each layer addresses a distinct architectural question: Which regulons are structurally essential (TOPPLE)? How do genetic effects distribute across cell types (DGSA)? What is the coupling topology (SICAI)? How do different perturbations interact with stability (IPA)? We further extend the framework with two analyses of chromatin-level architecture and multi-modal coupling convergence, and validate predictions in independent spatial transcriptomics data from psoriasis tissue.

---

## Results

### 1. TOPPLE identifies a stability hierarchy among 203 immune regulons

To identify which regulatory programs are structurally essential for immune homeostasis, we developed TOPPLE (Transcriptomic Organization via Perturbation of Program-Level Expression), a leave-one-out perturbation framework. For each of the 203 eRegulons profiled across 61 cell types in CIMA Table S5, we computed the redistribution index (RI): the mean Jensen-Shannon divergence between the original activity profile and all profiles obtained by removing one regulon at a time.

RI values spanned two orders of magnitude, from 0.00632 (HSF1) to 3.55 x 10^-5 (GATA1), revealing a continuous stability landscape rather than a binary classification (Fig. 1A). HSF1 ranked first (RI = 0.00632), followed by EGR1 (0.00574), KLF9 (0.00569), JUNB (0.00566), and JUN (0.00563). We classified the top quartile (n = 51) as stabilizers and the bottom quartile (n = 51) as destabilizers. The AP-1 transcription factor family was strikingly enriched among stabilizers (JUN rank 5, JUNB rank 4, FOSB rank 10), consistent with its known role as a master regulator of immune cell activation. Conversely, lineage-specifying factors occupied destabilizer positions: GATA1 (rank 201, RI = 3.55 x 10^-5), GATA2 (rank 200, RI = 3.64 x 10^-5), and FOXP3 (rank 192, RI = 1.50 x 10^-4).

This pattern reveals a design principle: factors that operate broadly across the immune system (stabilizers) contribute most to system-wide architecture, while factors that define specific lineages (destabilizers) are locally important but globally dispensable.

### 2. DGSA reveals near-universal non-additivity of eQTL effects that predicts disease pleiotropy

To characterize how genetic effects on gene expression are distributed across immune cell types, we developed DGSA (Decomposition of Genetic effect-Size Architecture). We analyzed 71,530 cis-eQTLs from CIMA Table S6, restricting to 5,253 genes with lead eQTLs in three or more cell types. For each gene, we decomposed the vector of effect sizes (slopes) into four geometric components: magnitude (L2 norm), uniformity (cosine similarity to a uniform vector), sparsity (fraction of zero entries), and non-additivity (1 minus the Gini coefficient of squared slopes).

The eQTL landscape was dominated by non-additivity: mean = 0.864, median = 0.911, s.d. = 0.126 (Fig. 2A). This indicates that for most genes, eQTL effects are highly unequal across cell types -- a minority of cell types carry the majority of genetic regulatory impact. Mean sparsity was 0.819, indicating that a typical gene has detectable eQTLs in approximately 12.5 of the 69 profiled cell types.

Critically, per-cell-type mean non-additivity predicted disease pleiotropy. We correlated each cell type's mean non-additivity with the number of disease-trait associations from CIMA's SMR analysis (Table S15, 68 cell types with both metrics). The correlation was strongly positive (Spearman rho = 0.796, P = 5.2 x 10^-16; Fig. 2B), indicating that cell types where eQTL effects are most unevenly distributed across the genome are precisely those most implicated in disease.

### 3. SICAI maps a coupling topology where complexity predicts disease burden

Inter-cellular coupling -- the degree to which genetic regulatory effects are shared across cell types -- is a fundamental but poorly characterized dimension of immune architecture. We developed SICAI (Systematic Inter-Cellular Architecture Index) to analyze the 69 x 69 pairwise eQTL sharing matrix (r_b coefficients from CIMA Table S8).

The global mean r_b was 0.819, consistent with the published CIMA estimate of 0.82, indicating extensive sharing of genetic regulatory architecture across immune cell types (Fig. 3A). However, coupling was not uniform. The most coupled cell type was CD4 Th22-like CCR10 (mean r_b = 0.882), while megakaryocytes (MK GP9) were most isolated (mean r_b = 0.603). We computed coupling complexity as Shannon entropy of each cell type's r_b distribution, with higher entropy indicating more uniformly distributed coupling.

Coupling complexity predicted disease burden. Correlating per-cell-type coupling complexity with the number of SMR disease associations (Table S15) yielded a significant positive association (Spearman rho = 0.509, P = 9.4 x 10^-6, n = 69; Fig. 3B). Cell types that are more uniformly coupled to the rest of the immune system -- those embedded in the densest regions of the coupling network -- carry the heaviest disease burden.

### 4. IPA reveals a striking dissociation between sex perturbation resistance and age sensitivity

If stabilizer regulons truly maintain system-wide homeostasis, they should resist environmental perturbations that could destabilize the immune system. We tested this prediction using IPA (Immune Perturbation Architecture), analyzing two natural perturbation axes: sex (Table S5, sex_difference_results sheet, 21,289 tests) and age (Table S5, Age_Correlation sheet, 404 regulons across up to 61 cell types).

For sex-based perturbation, the prediction held dramatically. Stabilizers exhibited significantly smaller sex-differential activity (median |log2FC| = 0.018) compared to destabilizers (median |log2FC| = 0.057; Mann-Whitney P = 1.96 x 10^-243; Fig. 4A). The correlation between RI and mean |log2FC| was strongly negative (Spearman rho = -0.574, P = 3.71 x 10^-19): the more stabilizing a regulon, the less it is perturbed by sex. The most perturbation-resistant regulons included PATZ1 (mean |log2FC| = 0.0069), PHF1 (0.0070), and KLF9 (0.0098), while extreme destabilizers showed over 100-fold greater perturbation: GATA1 (1.209), ZBTB8A (1.258), TCF7L1 (1.660).

For age, however, the pattern reversed. Stabilizers showed significantly *higher* age-correlated activity changes (mean |age rho| = 0.139) than destabilizers (mean |age rho| = 0.086; Mann-Whitney P = 1.27 x 10^-112; Fig. 4B). The correlation between RI and mean |age rho| was positive (rho = 0.435, P = 8.83 x 10^-11). A larger fraction of stabilizer-cell-type pairs showed significant age correlation (61.1%) than destabilizer pairs (39.6%).

This sex-age dissociation reveals that immune stability is not a monolithic property. Stabilizers resist acute, categorical perturbations (sex) while remaining responsive to gradual, continuous perturbations (aging) -- a pattern consistent with adaptive homeostasis rather than rigid buffering.

### 5. Epigenomic amplification: chromatin-level non-additivity exceeds expression-level non-additivity at matched loci

Having established the eQTL geometric landscape, we asked whether the same architectural principles apply to chromatin accessibility QTLs (caQTLs). We applied the identical DGSA pipeline to 151,875 cis-caQTLs from CIMA Table S6 (52,361 peaks across 42 cell types) and compared the geometric distributions.

At the aggregate level, caQTL and eQTL non-additivity distributions were statistically distinct (Kolmogorov-Smirnov D = 0.235, P = 4.8 x 10^-190) but with similar means (caQTL = 0.871, eQTL = 0.864). The critical test required matching: at 4,087 loci where the same genetic variant acts as both an eQTL and a caQTL in the same cell-type context, we performed paired comparisons. The Wilcoxon signed-rank test revealed that caQTL non-additivity was significantly higher (P = 7.0 x 10^-35), with 58.4% of matched pairs showing higher chromatin-level non-additivity (mean difference = +0.028).

We term this phenomenon "epigenomic amplification": genetic effects on chromatin accessibility are more cell-type-specific than the same variants' effects on gene expression. This has implications for disease genetics, as it suggests that GWAS signals operating through chromatin may have more cell-type-restricted functional consequences than those operating through expression.

The caQTL coupling matrix (42 x 42, from Table S8 cis_caQTL sheet) showed lower mean r_b (0.788) than the eQTL matrix (0.825; Mann-Whitney P = 3.9 x 10^-7), confirming that chromatin-level genetic architecture is less shared across cell types. Critically, per-cell-type caQTL mean r_b was a stronger predictor of disease burden (Spearman rho = 0.511, P = 5.5 x 10^-4) than any eQTL-derived metric, suggesting that chromatin coupling topology is more directly relevant to disease susceptibility.

### 6. Three coupling modalities converge into a hierarchical architecture

The immune system couples cell types through at least three modalities: shared genetic effects on expression (eQTL r_b), shared genetic effects on chromatin (caQTL r_b), and shared regulatory programs (enhancer-gene linkages). We asked whether these modalities converge or diverge.

Using CIMA Table S4 (SCENIC+ GRN, 469,157 TF-Region-Gene linkages comprising 133,574 unique region-gene pairs across 84,625 enhancer regions and 13,645 target genes) and Table S3 (338,036 peaks with binary activity across 65 cell types), we constructed a 65 x 65 regulatory coupling matrix based on Jaccard similarity of active region-gene pairs between cell types (mean Jaccard = 0.406).

Three-way Mantel tests revealed significant convergence across all modality pairs (all P < 0.0001; Fig. 3C): caQTL r_b vs. regulatory Jaccard showed the strongest correlation (r = 0.566), followed by eQTL r_b vs. caQTL r_b (r = 0.482), and eQTL r_b vs. regulatory Jaccard (r = 0.388). This ordering -- chromatin-regulatory > genetic-chromatin > genetic-expression -- suggests a hierarchy in which regulatory coupling is most directly mediated through shared chromatin architecture.

We further linked regulatory architecture to TOPPLE stability. Stabilizer regulons regulated 2.4-fold more target genes than destabilizers (mean 417 vs. 175; Mann-Whitney P = 7.9 x 10^-4) and exhibited broader cell-type activity (mean 29.0 vs. 23.1 cell types; P = 3.3 x 10^-4). Regulatory breadth (number of target region-gene pairs) correlated positively with TOPPLE RI (Spearman rho = 0.283, P = 4.3 x 10^-5), confirming that stabilizers achieve their architectural importance through extensive regulatory footprints.

### 7. Spatial validation confirms HIRA predictions in psoriasis tissue

To validate HIRA predictions in an independent tissue context, we analyzed spatial transcriptomics data from psoriasis (GSE202011 Visium, 30 samples: 14 lesional, 9 non-lesional, 7 healthy; 24,227 total spots).

TOPPLE predictions were confirmed: stabilizer expression was 3.3-fold higher in lesional tissue compared to healthy (0.257 vs. 0.078; P = 7.7 x 10^-4), while destabilizer expression showed a 2.7-fold increase (0.074 vs. 0.027; P = 1.6 x 10^-3). The stabilizer-to-destabilizer expression ratio was significantly elevated in lesional tissue (3.42 vs. 2.88 for healthy; Mann-Whitney P = 1.1 x 10^-3; Fig. 5), indicating that disease states amplify the activity of system-stabilizing programs.

Severity validation using PASI scores (6 psoriasis patients, PASI range 1.8-32.0, plus 3 healthy controls at PASI = 0) identified immune entropy -- a proxy for transcriptomic non-additivity -- as the metric most robustly correlated with disease severity in diseased patients alone (Spearman rho = 0.479, P = 0.021, n = 23). TOPPLE metrics (stabilizer expression rho = 0.446, P = 0.013; ratio rho = 0.404, P = 0.027) were significant when including healthy controls but lost significance within the diseased subset, suggesting that TOPPLE captures the healthy-to-disease transition while non-additivity scales with severity among affected individuals.

### 8. Cross-tissue coupling architecture persists from blood to skin

Finally, we tested whether the coupling topology identified in blood-derived CIMA data generalizes to tissue. Using STRATA-Atlas, we compared eQTL-derived r_b coupling (from CIMA, representing blood) with spatial co-expression coupling (from psoriasis Visium, representing skin) across 59 shared cell types (1,696 pairwise comparisons).

The Mantel test revealed a significant positive correlation (r = 0.143, P = 0.0001; Fig. 6), demonstrating that the genetic coupling architecture identified in blood partially persists in tissue. The moderate effect size is expected given the fundamental biological differences between circulating and tissue-resident immune populations, and suggests that a substantial fraction of inter-cellular coupling is intrinsic to genetic regulatory architecture rather than environmentally determined.

---

## Discussion

HIRA establishes five principal findings about the architectural organization of the human immune system.

**First, immune regulons are organized into a continuous stability hierarchy.** The TOPPLE analysis reveals that a small number of broadly expressed, system-stabilizing regulons (exemplified by HSF1 and the AP-1 family) exert disproportionate influence on immune architecture, while lineage-specifying factors (GATA1/2, FOXP3) are locally important but globally dispensable. This is reminiscent of the distinction between housekeeping and tissue-specific genes, but operates at the regulatory program level and is quantified by a perturbation-based metric rather than expression breadth alone.

**Second, eQTL effects are near-universally non-additive across cell types, and this non-additivity predicts disease.** The strong correlation between per-cell-type non-additivity and disease pleiotropy (rho = 0.796) suggests that cell types where genetic effects are most concentrated -- where a few variants carry outsized regulatory impact -- are precisely those most vulnerable to disease-associated genetic variation.

**Third, inter-cellular coupling topology predicts disease burden.** The SICAI finding that coupling complexity correlates with disease associations (rho = 0.509) implies that cell types embedded in dense coupling networks serve as conduits through which genetic perturbations propagate across the immune system. This has implications for understanding why certain cell types are disproportionately implicated in immune-mediated diseases.

**Fourth, the sex-age perturbation dissociation reveals adaptive homeostasis.** Stabilizers resist sex-based perturbation (an acute, categorical axis) yet are preferentially responsive to aging (a gradual, continuous axis). This suggests that immune stability is not rigid buffering but adaptive maintenance -- the system holds firm against acute perturbations while remaining responsive to gradual environmental change.

**Fifth, epigenomic amplification positions chromatin as a more cell-type-specific regulatory layer.** The finding that caQTL effects are more non-additive than eQTL effects at matched loci, combined with the observation that caQTL coupling is a stronger disease predictor than eQTL coupling, positions chromatin accessibility as a more cell-type-resolved -- and potentially more disease-relevant -- layer of genetic regulation. The three-way Mantel convergence further shows that regulatory coupling is most directly mediated through chromatin rather than expression.

**Limitations.** HIRA operates on summary statistics from a single atlas (CIMA), which may not capture the full diversity of human immune variation. The lead-eQTL design (one variant per gene per cell type) limits the resolution of DGSA and prevents multi-variant analyses. The psoriasis validation is limited to one disease in one tissue; broader validation across diseases and tissues is needed. The proxy metrics used for PASI severity validation (immune entropy for DGSA, coupling strength for SICAI) are approximations of the atlas-derived quantities and should be interpreted accordingly.

**Future directions.** HIRA could be applied to other multi-omics atlases as they become available, including the Tabula Sapiens and tissue-specific atlases. The framework could be extended to incorporate protein-level QTLs (pQTLs), splicing QTLs (sQTLs), and single-cell-resolution eQTL data. The disease predictions generated by DGSA and SICAI could be tested prospectively using genetic risk scores stratified by cell-type-specific eQTL architecture.

---

## Methods

### Data sources

All analyses used published supplementary tables from the Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026, DOI: 10.1126/science.adt3130). No raw sequencing data were accessed.

- **Table S5** (`science.adt3130_table_s5.xlsx`): eRegulon activity (sheet `eRegulons_Activators_Exp_AUC_RS`; columns: `eRegulon`, `cell_type_l4`, `mean_AUC`, `RSS`; 12,383 rows, 203 regulons x 61 cell types); sex difference (sheet `sex_difference_results`; columns: `cell_type_l4`, `eRegulon`, `Log2_Fold_Change`, `P_Value`, `adjust_P_Value`; 21,289 rows); age correlation (sheet `Age_Correlation`; wide format, 404 rows).
- **Table S6** (`CIMA_Table_S6.csv` from `science.adt3130_table_s6.zip`): cis-eQTLs (71,530 rows) and cis-caQTLs (151,875 rows); columns: `phenotype_id`, `celltype`, `slope`, `slope_se`, `pval_nominal`, `af`, `variant_id`, `analysis`.
- **Table S8** (`science.adt3130_table_s8.xlsx`): pairwise eQTL sharing (sheet `cis_eQTL`, 4,692 rows, columns: `reference_cell_type`, `query_celltype`, `rb`, `pi1`) and caQTL sharing (sheet `cis_caQTL`, 1,722 rows).
- **Tables S3/S4** (from zip archives): S3 provides binary peak activity across 65 cell types (338,036 peaks); S4 provides SCENIC+ GRN linkages (469,157 rows; columns: `TF`, `Region`, `Gene`, `R2G_importance`, `R2G_rho`).
- **Table S15** (`science.adt3130_table_s15.xlsx`): SMR pleiotropic associations (2,085 rows; columns: `Gene`, `celltype`, `trait`, `p_SMR`, `p_HEIDI`, `b_SMR`).

### TOPPLE: regulon stability via leave-one-out perturbation

For each regulon *r* among 203 eRegulons, we computed the activity profile *P_r(c)* = mean_AUC for regulon *r* in cell type *c*, normalized to sum to 1 across all 61 cell types. The redistribution index was computed as:

RI(r) = (1 / (N-1)) * sum_{r' != r} JSD(P_full, P_{-r'})

where JSD is the Jensen-Shannon divergence and P_{-r'} is the renormalized profile after removing regulon r'. Regulons were ranked by RI (higher = more stabilizing) and classified into stabilizers (top quartile, ranks 1-51), intermediates (ranks 52-152), and destabilizers (bottom quartile, ranks 153-203).

### DGSA: eQTL effect-size geometry decomposition

For each gene *g* with cis-eQTLs in >= 3 cell types, let **s_g** be the vector of signed slopes. Non-additivity was defined as 1 - Gini(|s_g|^2). Magnitude = L2 norm of |s_g|. Uniformity = cosine similarity of |s_g| to a uniform vector. Sparsity = fraction of zero entries. Disease correlation: per-cell-type mean non-additivity was correlated (Spearman) with number of disease-associated traits from Table S15.

### SICAI: inter-cellular coupling topology

The 69 x 69 pairwise r_b matrix was constructed from Table S8 cis_eQTL sheet by averaging r_b values for each ordered cell-type pair. Per-cell-type metrics: mean r_b, coefficient of variation (CV), Shannon entropy of the r_b distribution. Disease correlation: coupling complexity (Shannon entropy) vs. number of S15 disease traits.

### IPA: perturbation resistance

For sex perturbation: per-regulon mean |Log2_Fold_Change| across all cell types from Table S5 sex_difference_results. For age perturbation: per-regulon mean |Spearman rho| with age across cell types from Table S5 Age_Correlation (wide format, extracting value columns). Stabilizer vs. destabilizer comparison: two-sided Mann-Whitney U test.

### Extension 1: caQTL DGSA and coupling

Identical DGSA pipeline applied to 151,875 cis-caQTLs. Matched analysis: identified 4,087 variant-cell-type pairs where the same variant acts as both eQTL and caQTL; paired Wilcoxon signed-rank test on non-additivity. caQTL coupling: 42 x 42 r_b matrix from Table S8 cis_caQTL sheet. Mantel test between eQTL and caQTL r_b matrices (Pearson, 10,000 permutations).

### Extension 3: regulatory coupling

Regulatory coupling matrix: for each cell-type pair, Jaccard index of active region-gene pairs (a pair is active if the region is accessible in Table S3 AND the TF-Region-Gene linkage exists in Table S4). Three-way Mantel tests between eQTL r_b, caQTL r_b, and regulatory Jaccard matrices (restricted to shared cell types). Regulon breadth: number of target region-gene pairs per TF from Table S4, correlated with TOPPLE RI.

### STRATA: spatial validation

Visium spatial transcriptomics from GSE202011 (psoriasis, 30 samples). Stabilizer and destabilizer marker genes defined from TOPPLE top/bottom 20 regulons. Per-spot expression computed; sample-level means compared by condition (lesional vs. non-lesional vs. healthy) using Mann-Whitney U tests. PASI correlation: Spearman correlation of sample-level HIRA metrics with clinical PASI scores.

### STRATA-Atlas: cross-tissue coupling

For 59 cell types shared between CIMA (blood) and Visium (skin), computed pairwise spatial co-expression coupling (Pearson correlation of marker gene expression across spots) and compared with eQTL r_b coupling using Mantel test (10,000 permutations).

### Statistical analysis

All correlations are Spearman rank unless stated otherwise. Mann-Whitney U tests are two-sided unless stated. Multiple testing correction: Bonferroni where applicable. All analyses performed in Python 3.10+ using pandas, numpy, scipy, and scikit-bio (for Mantel tests).

---

## Figure Legends

**Figure 1. TOPPLE stability landscape of immune regulons.** (A) Redistribution index (RI) ranked across 203 eRegulons. Top stabilizers (red) and bottom destabilizers (blue) are highlighted. HSF1 ranks first (RI = 0.00632). (B) Stabilizer ranking showing the top 20 regulons by RI. AP-1 family members (JUN, JUNB, FOSB) are marked. (C) Lineage-specific destabilizers: GATA1, GATA2, and FOXP3 occupy the lowest RI positions. (D) Heatmap of RI values across cell-type lineages.

**Figure 2. DGSA reveals non-additive eQTL architecture that predicts disease.** (A) Distribution of non-additivity scores across 5,253 genes (mean = 0.864). (B) Per-cell-type mean non-additivity vs. number of disease associations (Spearman rho = 0.796, P = 5.2 x 10^-16, n = 68). (C) Non-additivity vs. TOPPLE stability class. (D) Top 20 genes by non-additivity with disease association counts.

**Figure 3. SICAI coupling topology and multi-modal convergence.** (A) Heatmap of 69 x 69 eQTL r_b sharing matrix (mean = 0.819). (B) Coupling complexity vs. disease burden (rho = 0.509, P = 9.4 x 10^-6). (C) Three-way Mantel correlations: caQTL-Jaccard (r = 0.566), eQTL-caQTL (r = 0.482), eQTL-Jaccard (r = 0.388); all P < 0.0001. (D) Stabilizers vs. destabilizers: target gene count (417 vs. 175, P = 7.9 x 10^-4).

**Figure 4. IPA sex-age perturbation dissociation.** (A) Stabilizers show lower sex-differential activity than destabilizers (median |log2FC| = 0.018 vs. 0.057; P = 1.96 x 10^-243). (B) Stabilizers show higher age sensitivity than destabilizers (mean |rho| = 0.139 vs. 0.086; P = 1.27 x 10^-112). (C) RI vs. perturbation sensitivity scatter, colored by stability class. (D) Heatmap of stabilizer perturbation magnitude across cell types.

**Supplementary Figure 1. Epigenomic amplification (Extension 1).** (A) caQTL vs. eQTL non-additivity distributions (KS D = 0.235, P = 4.8 x 10^-190). (B) Matched-locus paired comparison (Wilcoxon P = 7.0 x 10^-35, 58.4% favor caQTL). (C) caQTL coupling matrix (42 x 42, mean r_b = 0.788). (D) caQTL mean r_b vs. disease burden (rho = 0.511, P = 5.5 x 10^-4).

**Supplementary Figure 2. STRATA spatial validation.** (A) Stabilizer expression across conditions (lesional 0.257 vs. healthy 0.078, P = 7.7 x 10^-4). (B) Stabilizer/destabilizer ratio (3.42 vs. 2.88, P = 1.1 x 10^-3). (C) PASI severity correlations for TOPPLE (rho = 0.446), immune entropy (rho = 0.479), and coupling metrics. (D) STRATA-Atlas Mantel test (r = 0.143, P = 0.0001).

---

## Data and Code Availability

All source data are from CIMA supplementary tables (Yin et al., Science 2026). Analysis code is available at https://github.com/jengweitjiu/HIRA.

## Competing Interests

The author declares no competing interests.

## References

[Placeholder for full reference list]

1. Yin R, et al. Chinese Immune Multi-Omics Atlas. Science. 2026.
2. Regev A, et al. The Human Cell Atlas. eLife. 2017.
3. Tabula Sapiens Consortium. The Tabula Sapiens. Science. 2022.
4. Barabasi AL, Oltvai ZN. Network biology. Nat Rev Genet. 2004.
5. Consortium GTEx. Genetic effects on gene expression across human tissues. Nature. 2017.
6. Urbut SM, et al. Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions. Nat Genet. 2019.
