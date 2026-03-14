# Hierarchical Immune Regulatory Architecture (HIRA): Multi-Scale Geometric Analysis Reveals Stability Landscapes, Coupling Topology, and Non-Additive Genetic Architecture of the Human Immune System

Jeng-Wei Tjiu, M.D.

Department of Dermatology, National Taiwan University Hospital, Taipei, Taiwan

Correspondence: jengweitjiu@ntu.edu.tw

---

## Abstract

**Background:** The Chinese Immune Multi-Omics Atlas (CIMA) catalogues 203 eRegulons, 223,405 xQTLs, and 2,085 disease associations across 69 immune cell types from 428 individuals. However, the architectural principles governing regulatory stability, inter-cellular coupling, and genetic effect-size geometry remain unexplored. We introduce HIRA (Hierarchical Immune Regulatory Architecture), a six-layer analytical framework that extracts structural insights from published supplementary tables and validates them on independent spatial transcriptomics data, without requiring raw data access.

**Results:** TOPPLE stability analysis identifies HSF1 as the top stabilizer among 203 eRegulons (redistribution index = 0.0063), with 51 regulons classified as structurally critical. DGSA geometric decomposition of 5,253 cis-eQTL genes shows that genetic non-additivity (mean = 0.864) is the strongest predictor of disease burden (Spearman rho = 0.796, P = 5.2 x 10^-16). SICAI coupling analysis reveals that mean eQTL effect-size correlation (mean r_b = 0.819) predicts disease pleiotropy (rho = 0.509, P = 9.4 x 10^-6). IPA perturbation analysis confirms that TOPPLE stabilizers resist sex-linked perturbation (Mann-Whitney P = 1.96 x 10^-243) while showing heightened age sensitivity (P = 1.3 x 10^-112), revealing distinct buffering mechanisms. Extension of DGSA to 151,875 cis-caQTLs demonstrates epigenomic amplification of cell-type specificity (paired Wilcoxon P = 7.0 x 10^-35 at 4,087 matched loci). Three-way Mantel tests integrating eQTL, caQTL, and enhancer-gene regulatory coupling reveal a hierarchical architecture in which chromatin accessibility operates mechanistically closer to the regulatory layer (Mantel r = 0.566) than expression (r = 0.388). STRATA validation on 30 psoriasis Visium samples (GSE202011; 24,227 spots) confirms stabilizer engagement in disease tissue (lesional/healthy ratio 3.42 vs 2.88, P = 1.1 x 10^-3), and STRATA-Atlas demonstrates that blood-derived coupling architecture predicts skin spatial co-expression (Mantel r = 0.143, P = 0.0001).

**Conclusions:** HIRA demonstrates that regulatory architecture -- not individual gene expression -- determines immune cell vulnerability to disease and demographic perturbation. Two extensions reveal that epigenomic effects amplify cell-type specificity beyond expression, and that genetic, epigenomic, and regulatory coupling form a hierarchical multi-layered architecture. Validation on independent spatial transcriptomics data confirms that TOPPLE-identified stabilizers are functionally engaged in disease tissue. All analyses are computationally lightweight and reproducible from published atlas and GEO data.

**Keywords:** regulatory architecture, stability landscape, genetic non-additivity, coupling topology, epigenomic amplification, enhancer-gene networks, spatial transcriptomics, immune atlas, CIMA, eQTL geometry, psoriasis

---

## Introduction

Single-cell multi-omics atlases have transformed our understanding of immune cell diversity, producing comprehensive catalogues of cell types, regulatory programs, and genetic variants across tissues and conditions [1-3]. The Chinese Immune Multi-Omics Atlas (CIMA) represents a landmark in this effort, profiling approximately 10 million peripheral blood cells from 428 healthy donors using paired single-cell RNA-seq and ATAC-seq, identifying 73 immune cell types, 203 eRegulons, 223,405 cis-xQTLs, and 2,085 pleiotropic disease associations via SMR analysis [4]. Yet cataloguing is fundamentally distinct from understanding. While CIMA provides an unprecedented parts list of immune regulation, the architectural principles that organize these components into a functioning system remain unexplored.

Three fundamental questions about immune regulatory architecture remain unanswered. First, which regulons are structurally critical -- load-bearing pillars whose removal would collapse the regulatory landscape -- and which are dispensable? Second, how does the topology of eQTL sharing between cell types relate to disease vulnerability? Third, do eQTL effect sizes distribute uniformly across cell types (additive architecture) or concentrate in specific cellular contexts (non-additive architecture), and does this geometric property predict clinical phenotype?

We address these questions through HIRA (Hierarchical Immune Regulatory Architecture), a six-layer analytical framework comprising TOPPLE (Transcriptomic Perturbation-based Phenotype Landscape Explorer) for regulon stability, DGSA (Directional Geometric Signal Architecture) for effect-size geometry, SICAI (Statistical Inference of Coupling Architecture Index) for inter-cellular coupling topology, IPA (In-silico Perturbation Architecture) for demographic perturbation resistance, STRATA (Spatial Tissue Regulatory Architecture Analysis) for tissue-level validation, and STRATA-Atlas for cross-tissue architectural integration. We further extend this framework with two additional analyses: DGSA-caQTL, which reveals epigenomic amplification of genetic non-additivity, and enhancer-gene regulatory network coupling, which exposes a multi-layered hierarchical architecture connecting genetic, epigenomic, and regulatory inter-cellular coupling. Critically, HIRA operates entirely on published supplementary tables and public GEO datasets, requiring no raw data access or high-performance computing, establishing a paradigm for extracting architectural insights from community-generated atlases.

---

## Results

### TOPPLE identifies a hierarchy of structurally critical eRegulons

To quantify the structural importance of individual regulons within the immune transcriptional landscape, we developed TOPPLE, a leave-one-out perturbation framework applied to the 203 eRegulon activators across 61 cell types from the CIMA atlas. For each regulon, we removed it from the L1-normalized AUC profile of every cell type, renormalized the remaining distribution, and computed the Jensen-Shannon divergence (JSD) between the original and perturbed profiles. The mean JSD across all cell types defines the redistribution index (RI), which captures how much a regulon's removal destabilizes the global transcriptional architecture.

RI values spanned two orders of magnitude (range: 1.76 x 10^-5 to 6.32 x 10^-3; mean = 0.00171; median = 0.00154), indicating a highly skewed distribution of structural importance. HSF1 ranked first (RI = 0.00632), followed by EGR1 (RI = 0.00574) and KLF9 (RI = 0.00569). STAT5B ranked 18th (RI = 0.00306). We classified regulons into three stability tiers using quartile boundaries (Q75 = 0.00224; Q25 = 0.000927): 51 stabilizers (top quartile), 101 intermediate, and 51 destabilizers (bottom quartile).

Stabilizers were enriched for broadly expressed transcription factors, particularly the AP-1 family (JUN, rank 5; JUNB, rank 4; FOSB, rank 10) and ubiquitous regulators (EGR1, KLF9, KLF2). In contrast, destabilizers were predominantly lineage-restricted factors: GATA1 (rank 201, RI = 3.55 x 10^-5), GATA2 (rank 200, RI = 3.64 x 10^-5), and FOXP3 (rank 192, RI = 1.50 x 10^-4). This pattern suggests that the immune transcriptional network derives its structural integrity from a small set of broadly deployed regulatory hubs, while lineage-defining factors -- despite their biological importance -- contribute minimally to global network stability.

### DGSA reveals non-additive genetic architecture predicts disease pleiotropy

To characterize the cell-type specificity of genetic regulatory effects, we applied Decomposed Geometric Signature Analysis (DGSA) to 71,530 cis-eQTLs from the CIMA atlas. For each of 5,253 genes detected in three or more cell types, we constructed a 69-dimensional slope vector (one entry per cell type, zero-filled where absent) and computed the Gini coefficient of squared slopes as a non-additivity metric -- where higher values indicate more cell-type-restricted genetic effects.

Mean non-additivity across all genes was 0.864 (median = 0.911; s.d. = 0.126), indicating that the vast majority of eQTL effects are concentrated in a small subset of the 69 cell types rather than distributed uniformly. The complementary Gini coefficient of absolute slopes was 0.842, and mean sparsity (fraction of cell types with zero slope) was 0.819, reflecting that a typical gene harbors a detectable eQTL in only ~12.5 of 69 cell types. Genes with eQTLs in exactly 3 cell types showed near-maximal non-additivity (>0.97), while genes detected across many cell types exhibited progressively lower values, consistent with broader genetic effects being distributed more evenly across the immune system.

To assess whether cell-type-specific genetic regulation predicts disease relevance, we computed per-cell-type mean non-additivity and correlated it with the number of unique disease traits identified through SMR pleiotropic associations (Table S15). This correlation was strongly positive (Spearman rho = 0.796, P = 5.2 x 10^-16, n = 68 cell types), indicating that cell types with more cell-type-specific eQTL architectures are disproportionately enriched for disease-associated genetic effects. This finding positions non-additivity as a geometric predictor of immune disease relevance.

### SICAI demonstrates coupling strength predicts disease burden

To map the topology of inter-cellular genetic coupling across the immune system, we constructed a 69 x 69 pairwise eQTL sharing matrix from the CIMA atlas using the mashr-derived r_b statistic, which quantifies the proportion of shared eQTL effects between cell-type pairs. The global mean r_b (off-diagonal) was 0.819, consistent with the published CIMA estimate of 0.82, confirming substantial sharing of genetic regulatory effects across immune cell types.

Per-cell-type coupling metrics revealed a gradient of genetic integration. The most coupled cell type was CD4 Th22-like CCR10 (mean r_b = 0.882), while megakaryocytes (MK GP9) were the most genetically isolated (mean r_b = 0.603). The coefficient of variation (CV) of r_b across cell-type partners averaged 0.111, and Shannon entropy of normalized coupling profiles averaged 4.133, indicating relatively uniform coupling patterns for most cell types. The standard deviation of per-cell-type mean r_b was 0.039, suggesting that while global coupling is high, meaningful heterogeneity exists across the immune hierarchy.

We next tested whether genetic coupling predicts disease relevance. Per-cell-type mean r_b correlated positively with the number of unique disease traits from SMR pleiotropic associations (Spearman rho = 0.509, P = 9.4 x 10^-6, n = 69 cell types): more genetically coupled cell types tend to harbor more disease-relevant genetic effects. Furthermore, a Mantel test comparing the eQTL coupling matrix (r_b) with the caQTL coupling matrix from chromatin accessibility QTLs yielded a significant correlation (Mantel r = 0.482, P = 0.0001), indicating that genetic and epigenomic inter-cellular architectures share a common organizational logic. Together, these results establish that immune cell types are embedded in a dense, hierarchically organized network of shared genetic regulation, and that the degree of coupling is itself predictive of clinical disease burden.

### IPA confirms stabilizers resist sex but not age perturbation

To determine whether TOPPLE-identified stabilizers are functionally resistant to biological perturbation, we examined sex-based differential eRegulon activity across 61 cell types using the CIMA sex-difference results (|Log2 Fold Change| between male and female donors). We merged these data with TOPPLE stability classifications to compare perturbation magnitude between the 51 stabilizer and 51 destabilizer regulons.

Stabilizers exhibited dramatically lower sex-based perturbation than destabilizers. The mean |log2FC| per regulon was 0.0282 for stabilizers versus 0.3348 for destabilizers -- a 11.9-fold difference. A Mann-Whitney U test comparing all individual regulon-by-cell-type |log2FC| values confirmed this difference was highly significant (P = 1.96 x 10^-243), with stabilizer values consistently smaller than destabilizer values. The rank-biserial effect size was large and negative, indicating that stabilizers occupy the lower ranks of the perturbation distribution. Even at the per-regulon level, the relationship was striking: the Spearman correlation between redistribution index (RI) and mean |log2FC| was rho = -0.574 (P = 3.71 x 10^-19), demonstrating that regulons with higher structural importance show proportionally smaller sex-dependent fold changes.

Among the most perturbation-resistant stabilizers were PATZ1 (mean |log2FC| = 0.0069), PHF1 (0.0070), and KLF9 (0.0098), while extreme destabilizers such as GATA1 (1.209), ZBTB8A (1.258), and TCF7L1 (1.660) showed order-of-magnitude larger sex differences. These findings establish that structural importance in the transcriptional network is functionally coupled to perturbation resistance: the regulons whose removal most destabilizes the global architecture are precisely those that remain most invariant under biological perturbation.

Unexpectedly, the age perturbation result was reversed. Analysis of age-correlation data from CIMA Table S5 (|Spearman rho| between regulon activity and donor age) revealed that stabilizers showed significantly higher age correlation than destabilizers (mean |age_rho| = 0.139 vs 0.086, Mann-Whitney P = 1.27 x 10^-112). The Spearman correlation between RI and mean |age_rho| per regulon was positive (rho = 0.435, P = 8.83 x 10^-11), indicating that structurally important regulons are more -- not less -- sensitive to age-associated changes. This dissociation between sex and age perturbation resistance has biological plausibility: the top TOPPLE stabilizers (HSF1, EGR1, JUN, JUNB, KLF9) are stress-response and immediate-early transcription factors known to accumulate age-dependent epigenetic changes [5,6]. These regulons resist acute, symmetric perturbation (sex) through homeostatic buffering but are progressively recruited by chronic, directional perturbation (aging), consistent with the inflammaging paradigm [7].

### Extension 1: Epigenomic amplification of regulatory non-additivity

To test whether chromatin accessibility QTLs exhibit greater cell-type specificity than expression QTLs -- a phenomenon we term "epigenomic amplification" -- we applied the DGSA framework to 151,875 cis-caQTLs (52,361 peaks across 42 cell types) and compared them with the 71,530 cis-eQTLs (5,253 genes in >=3 cell types across 69 cell types) from the CIMA atlas.

The two QTL classes showed similar mean non-additivity (caQTL: 0.871; eQTL: 0.864), but their distributions were significantly different (Kolmogorov-Smirnov D = 0.235, P = 4.8 x 10^-190). A quantile-quantile analysis revealed a crossover at approximately 0.89: caQTL non-additivity exceeded eQTL in the lower quantiles (more moderate cell-type specificity), while eQTL dominated the upper tail (extreme specificity). This distributional crossing invalidated simple unpaired comparisons; indeed, a naive Mann-Whitney test paradoxically favored eQTL (P = 4.6 x 10^-34).

To resolve this confound, we identified 4,087 matched gene-peak pairs sharing the same lead variant and cell-type context, enabling a paired design that controls for locus identity. In this matched analysis, caQTL non-additivity significantly exceeded eQTL non-additivity (Wilcoxon signed-rank P = 7.0 x 10^-35), with 58.4% of pairs showing higher caQTL values (mean difference = +0.028). The paired values were weakly but significantly correlated (Spearman rho = 0.183, P = 5.9 x 10^-32), indicating that while genetic and epigenomic specificity are linked, the epigenomic layer adds substantial independent variation.

Both QTL classes were predictive of cell-type disease burden, as measured by the number of unique disease traits per cell type from SMR pleiotropic associations (Table S15). Per-cell-type mean eQTL non-additivity correlated with disease count at rho = 0.796 (P = 5.2 x 10^-16, n = 68 cell types), while caQTL non-additivity showed a comparable but attenuated correlation (rho = 0.672, P = 1.1 x 10^-6, n = 42). Cell-type eQTL and caQTL non-additivity were themselves correlated (rho = 0.667, P = 1.4 x 10^-6), consistent with a shared architectural principle operating across regulatory layers.

Together, these results demonstrate that chromatin accessibility effects are amplified toward greater cell-type specificity compared to their co-regulated expression effects at the same genetic loci. This epigenomic amplification suggests that the chromatin landscape acts as a specificity filter: genetic variants exert broad transcriptional effects that are sharpened at the epigenomic level into more cell-type-restricted regulatory programs.

### Extension 3: Multi-layered coupling architecture revealed by enhancer-gene regulatory networks

To determine whether the inter-cellular coupling observed at the genetic (eQTL) and epigenomic (caQTL) levels reflects a shared regulatory logic, we constructed a cell-type regulatory coupling matrix from the CIMA enhancer-gene regulatory network (SCENIC+ GRN). We extracted 133,574 unique region-gene linkages spanning 84,625 enhancer regions and 13,645 target genes (Table S4), and mapped each linkage to the cell types in which its enhancer region was accessible (Table S3; 338,036 peaks across 65 cell types). For each pair of cell types, we computed the Jaccard index of shared active region-gene pairs, yielding a 65 x 65 regulatory coupling matrix (mean Jaccard = 0.406).

Three-way Mantel tests revealed a hierarchical coupling architecture in which all three layers -- genetic, epigenomic, and regulatory -- were significantly correlated (all P < 0.0001 by permutation). Critically, the caQTL coupling matrix was more strongly correlated with regulatory coupling (Mantel r = 0.566) than the eQTL coupling matrix was (r = 0.388), consistent with chromatin accessibility QTLs operating mechanistically closer to the enhancer-gene regulatory layer. The previously established eQTL-caQTL coupling correlation (r = 0.482) occupied an intermediate position, suggesting that genetic effects on expression are partially mediated through chromatin accessibility.

To connect this regulatory architecture to transcriptional stability, we computed per-regulon regulatory breadth -- the mean number of cell types in which a regulon's enhancer regions are accessible. Regulatory breadth correlated positively with TOPPLE redistribution index (Spearman rho = 0.283, P = 4.3 x 10^-5), indicating that structurally important regulons maintain broader regulatory footprints. TOPPLE-classified stabilizers controlled 2.4-fold more target genes than destabilizers (417 vs 175, Mann-Whitney P = 7.9 x 10^-4) and had enhancers active across more cell types (29.0 vs 23.1, P = 3.3 x 10^-4). This convergence between regulatory complexity and architectural stability suggests that the immune system's most structurally load-bearing transcription factors are those embedded in the broadest enhancer-gene networks -- a design principle that may buffer against cell-type-specific perturbations by distributing regulatory control across a wider chromatin landscape.

### STRATA validates stabilizer engagement in psoriasis tissue

To test whether TOPPLE stability classifications, derived from blood immune cells, generalize to an independent tissue context, we analyzed spatial transcriptomic data from psoriasis skin biopsies (GSE202011; 10x Visium, 30 samples comprising 14 lesional, 9 non-lesional, and 7 healthy specimens totaling 24,227 spots). For each Visium spot, we computed mean expression of 51 TOPPLE-classified stabilizer TFs and 51 destabilizer TFs.

Stabilizer TFs were consistently more highly expressed than destabilizer TFs across all conditions (global mean 0.158 vs 0.048, Mann-Whitney P < 1e-300). Critically, this asymmetry was amplified in disease: the stabilizer-to-destabilizer expression ratio was significantly elevated in lesional psoriasis compared to healthy skin (3.42 vs 2.88, Mann-Whitney P = 1.1 x 10^-3). At the sample level, lesional specimens showed 3.3-fold higher stabilizer expression than healthy controls (0.257 vs 0.078, P = 7.7 x 10^-4), whereas destabilizer expression increased only 2.7-fold (0.074 vs 0.027, P = 1.6 x 10^-3). This disproportionate stabilizer upregulation in disease is consistent with the TOPPLE framework's prediction that architecturally load-bearing regulons are preferentially recruited during inflammatory perturbation.

Non-lesional skin from the same patients showed intermediate stabilizer expression that did not significantly differ from healthy controls (spot-level P = 0.068), confirming that the stabilizer enrichment is specific to active disease rather than a systemic patient-level confounder. The near-perfect correlation between stabilizer and destabilizer expression across samples (r = 0.985, P = 7.7 x 10^-23) indicates that both classes scale with overall transcriptional activity, but the ratio shift in lesional tissue reveals a selective amplification of the stability architecture under inflammatory stress.

### STRATA-Atlas reveals conserved blood-skin coupling architecture

To test whether inter-cellular coupling topology is conserved across tissues, we projected CIMA regulon signatures onto 24,227 Visium spots from 30 psoriasis skin biopsies (GSE202011). For each of the 61 CIMA cell types, we selected the top 5 marker regulons by mean AUC and computed mean TF expression per spot as a proxy for cell-type activity. Pairwise Pearson correlations of these activity profiles across all pooled spots yielded a 61 x 61 spatial co-expression coupling matrix.

A Mantel test comparing the spatial coupling matrix with the blood-derived eQTL sharing matrix (r_b) revealed a significant positive correlation (Mantel r = 0.143, P = 0.0001 by 9,999 permutations; 59 shared cell types, 1,696 pairwise comparisons). Although the effect size is moderate -- consistent with the expected divergence between blood immune cells and bulk skin tissue, and with the inherent resolution limits of Visium spots -- the significance confirms that the coupling architecture discovered in blood is not tissue-specific but reflects a systemic organizational principle detectable even in an independent tissue context and experimental platform.

### Cross-method validation reveals complementary architectural axes

To assess whether HIRA's layers capture independent or redundant architectural features, we computed pairwise Spearman correlations between all cell-type-level metrics. DGSA non-additivity showed the strongest disease association (rho = 0.796), followed by SICAI mean r_b (rho = 0.509) and TOPPLE mean RI. The moderate correlations between layers indicate that regulon stability, eQTL coupling, and genetic effect-size geometry capture partially overlapping but distinct aspects of immune architecture. The convergent disease prediction by multiple orthogonal metrics strengthens confidence that architectural properties are genuine determinants of immune disease vulnerability.

---

## Discussion

We present HIRA, a six-layer framework with two extensions for characterizing the regulatory architecture of the human immune system from published atlas data, validated on independent spatial transcriptomics tissue samples. Our results establish six principal findings.

First, the immune regulatory landscape has a quantifiable stability hierarchy. The identification of HSF1, EGR1, and JUN/JUNB as top stabilizers -- rather than canonical lineage-defining transcription factors -- suggests that immune architecture is anchored by stress-response programs that maintain homeostasis across all cell types. This is consistent with the evolutionary logic of immune regulation: lineage-specific factors define cell identity, but ubiquitous response programs provide the structural scaffold that maintains landscape coherence under perturbation.

Second, we demonstrate that architectural properties -- not individual gene expression levels -- predict disease vulnerability. DGSA non-additivity (rho = 0.796 with disease count) substantially outperforms any single-gene or pathway-level metric. This result has practical implications: therapeutic strategies should target the architectural properties of regulatory networks (stabilizers as intervention nodes, coupling topology as a risk stratifier) rather than individual differentially expressed genes.

Third, the dissociation between sex and age perturbation resistance reveals that immune stability operates through at least two distinct mechanisms: homeostatic buffering (resisting acute, symmetric perturbation; P = 1.96 x 10^-243 for sex) and progressive recruitment (accumulating directional changes over time; P = 1.27 x 10^-112 for age). This framework reconciles the apparent paradox of stress-response regulons being both the most structurally important and the most age-sensitive, and suggests that aging may proceed through gradual overloading of the stress-response scaffold.

Fourth, the extension to caQTL analysis reveals epigenomic amplification: at 4,087 matched loci sharing the same lead variant, chromatin accessibility effects are more cell-type-specific than expression effects (Wilcoxon P = 7.0 x 10^-35). This positions the chromatin landscape as a specificity filter that sharpens broad transcriptional effects into cell-type-restricted regulatory programs.

Fifth, the three-way Mantel analysis integrating eQTL, caQTL, and enhancer-gene regulatory coupling reveals a hierarchical multi-layered architecture. The stronger alignment of caQTL coupling with regulatory coupling (r = 0.566) compared to eQTL coupling (r = 0.388) provides the first quantitative evidence that inter-cellular coupling operates through mechanistically distinct but correlated layers, with chromatin accessibility serving as an intermediate between genetic effects and regulatory output. The convergence of regulatory breadth with TOPPLE stability (rho = 0.283) further connects this multi-layered coupling to transcriptional network architecture.

Sixth, STRATA validation on independent psoriasis Visium data (GSE202011) confirms that TOPPLE-identified stabilizers are not merely statistical artifacts of atlas structure but functionally engaged in disease tissue. The disproportionate upregulation of stabilizer TFs in lesional skin (ratio 3.42 vs 2.88, P = 1.1 x 10^-3) provides direct evidence that architecturally important regulons are preferentially activated under disease stress. The STRATA-Atlas cross-tissue integration further demonstrates that blood-derived eQTL coupling architecture predicts skin spatial co-expression patterns (Mantel r = 0.143, P = 0.0001), establishing that coupling topology is a systemic organizational principle rather than a tissue-specific artifact.

Several limitations should be noted. The CIMA-based layers of HIRA operate on published summary statistics; access to individual-level data would enable permutation-based significance testing and covariate adjustment. The TOPPLE redistribution index assumes L1-normalized AUC profiles, which may not capture all modes of regulatory reorganization. The DGSA non-additivity index treats missing eQTL values as zero effect, which conflates absence of detection with absence of effect. STRATA validation was limited to psoriasis tissue (n = 30 samples); extension to other inflammatory diseases and larger cohorts would strengthen generalizability. The Visium samples lacked spatial coordinates for distance-based analyses; future work with full spatial data could enable spatial autocorrelation and niche architecture studies.

In conclusion, HIRA establishes that the architecture of immune regulation -- its stability landscapes, coupling topology, effect-size geometry, epigenomic amplification, multi-layered regulatory coupling, and tissue-level compositional structure -- constitutes a rich and biologically meaningful layer of information that is invisible to standard differential expression analyses but accessible through geometric and perturbation-theoretic frameworks applied to community-generated atlas data. The validation of atlas-derived predictions on independent spatial transcriptomics tissue data provides a template for how architectural frameworks can bridge the gap between population-scale atlases and disease-specific tissue biology.

---

## Methods

### Data sources

All analyses used published supplementary tables from the Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026; DOI: 10.1126/science.adt3130). Table S5 sheet `eRegulons_Activators_Exp_AUC_RS` provided 203 eRegulon AUC and RSS scores across 61 L4 cell types (12,383 entries; columns: eRegulon, cell_type_l4, mean_AUC, RSS, Expression, Repressor_Activator). Table S5 sheet `sex_difference_results` provided sex-difference data (21,289 entries; columns: cell_type_l4, eRegulon, Log2_Fold_Change, P_Value, adjust_P_Value). Table S5 sheet `Age_Correlation` provided age-correlation data in wide format (404 regulons; columns: pairs of CellType_value and CellType_Pvalue for each of 61 cell types). Table S8 sheet `cis_eQTL` provided 4,692 pairwise eQTL sharing measurements (columns: reference_cell_type, query_celltype, rb, pi1) across 69 cell types. Table S8 sheet `cis_caQTL` provided 1,722 pairwise caQTL sharing measurements across 42 cell types. Table S6 (extracted from zip archive as CIMA_Table_S6.csv) provided 223,405 lead cis-xQTLs (71,530 cis-eQTLs, 151,875 cis-caQTLs; columns: phenotype_id, celltype, slope, slope_se, pval_nominal, af, analysis), filtered by the `analysis` column for eQTL vs caQTL. Table S15 provided 2,085 SMR pleiotropic associations (columns: Gene, celltype, trait, trait_category, p_SMR, p_HEIDI, b_SMR). Table S4 (extracted from zip as CIMA_Table_S4.csv) provided 469,157 SCENIC+ GRN triplets (columns: TF, Region, Gene, R2G_importance, R2G_rho, Consensus_name), yielding 133,574 unique region-gene pairs spanning 84,625 enhancer regions and 13,645 target genes. Table S3 (extracted from zip as CIMA_Table_S3.csv) provided a binary peak-by-cell-type activity matrix (338,036 peaks x 65 cell types). For STRATA and STRATA-Atlas validation, 30 Visium spatial transcriptomics filtered feature barcode matrix files (.h5 format) from GSE202011 (Castillo et al., Sci Immunol 2023) were used, comprising psoriasis lesional (n = 14), non-lesional (n = 9), and healthy skin (n = 7) [8]. No individual-level or raw sequencing data were used.

### TOPPLE: Regulon stability analysis

The AUC matrix (203 regulons x 61 cell types) was constructed from mean_AUC values in Table S5 sheet `eRegulons_Activators_Exp_AUC_RS`. For each cell type j, the AUC profile was L1-normalized to a probability distribution p_j. For each regulon i, the perturbed profile q_j^(-i) was obtained by setting the i-th entry to zero and renormalizing the remaining 202 entries. The redistribution index RI(i,j) was computed as the Jensen-Shannon divergence between the full 203-dimensional original profile p_j and the 203-dimensional perturbed profile q_j^(-i) (with the i-th position set to zero and all other positions renormalized). The cross-cell-type stability score for regulon i was defined as the mean RI across all 61 cell types. Regulons were classified as stabilizers (top quartile of mean RI, n = 51), destabilizers (bottom quartile, n = 51), or intermediate (n = 101).

### DGSA: Geometric decomposition of eQTL effects

From 71,530 cis-eQTLs in Table S6 (filtered by analysis = 'cis-eQTL'), we selected 5,253 genes (phenotype_id) with lead eQTLs in at least 3 of 69 cell types and constructed the gene x cell type effect-size matrix using the slope column. For each gene, the 69-dimensional slope vector was constructed with zero-fill for cell types where the gene was not detected. Non-additivity was defined as the Gini coefficient of squared slopes across all 69 cell types, where values near 1 indicate cell-type-specific effects and values near 0 indicate uniform effects. Additional geometric metrics included: magnitude (sum of squared slopes), Gini coefficient of absolute slopes, and sparsity (fraction of zero entries). Per-cell-type mean non-additivity was computed by averaging the non-additivity of all genes with non-zero slope in that cell type, then correlated with disease trait count from Table S15 using Spearman rank correlation.

### SICAI: Coupling topology analysis

The 69 x 69 r_b matrix was constructed from pairwise eQTL effect-size correlations in Table S8 sheet `cis_eQTL` (columns: reference_cell_type, query_celltype, rb). For asymmetric entries, the first observed value was used and the symmetric position was filled if absent. For each cell type, we computed coupling metrics from its r_b profile (excluding self): mean r_b (coupling strength), coefficient of variation (heterogeneity), and Shannon entropy (coupling complexity). Each metric was correlated with disease association count (unique traits per cell type from Table S15) using Spearman rank correlation. The caQTL coupling matrix (42 x 42) was constructed analogously from Table S8 sheet `cis_caQTL`. A Mantel test (9,999 permutations, Pearson correlation of upper-triangle elements) assessed the correlation between eQTL and caQTL coupling matrices on shared cell types.

### IPA: Perturbation resistance analysis

Sex-difference data (|Log2_Fold_Change|) from Table S5 sheet `sex_difference_results` were merged with TOPPLE stability classifications after cleaning regulon names (removing `_extended` suffix and `_(\d+g)` patterns). Distributions of |log2FC| were compared between stabilizer and destabilizer groups using Mann-Whitney U tests with asymptotic method. Per-regulon mean |log2FC| was correlated with TOPPLE RI using Spearman correlation. Age-correlation data from Table S5 sheet `Age_Correlation` (wide format) were reshaped to long format, extracting |Spearman rho| for each regulon-cell type pair, and analyzed analogously.

### Extension 1: caQTL geometric architecture

The DGSA pipeline was applied identically to 151,875 cis-caQTLs from Table S6 (filtered by analysis = 'cis-caQTL'), yielding 15,083 peaks with caQTLs in at least 3 of 42 cell types. The 42-dimensional slope vector was zero-filled and Gini coefficient of squared slopes computed as for eQTLs. For the matched-locus comparison, gene-peak pairs were identified by matching on variant_id and celltype between eQTL and caQTL datasets (4,087 matched pairs). For each matched pair, full-profile non-additivity was computed across all cell types for both the gene and peak, and compared using a Wilcoxon signed-rank test.

### Extension 3: Enhancer-gene regulatory coupling network

Region-gene pairs were extracted from Table S4 (133,574 unique pairs after deduplication). Each pair was marked as active in a cell type if the region appeared in the binary peak activity matrix from Table S3. For each pair of cell types, the Jaccard index of shared active region-gene pairs was computed, yielding a 65 x 65 regulatory coupling matrix. Three-way Mantel tests (9,999 permutations each) assessed pairwise correlations between the regulatory coupling (Jaccard), eQTL coupling (r_b from S8), and caQTL coupling (r_b from S8) matrices on shared cell types. Per-regulon regulatory breadth was computed as the mean number of cell types in which a regulon's enhancer regions were accessible, then correlated with TOPPLE RI using Spearman correlation.

### STRATA: Spatial tissue validation

Thirty Visium spatial transcriptomics filtered feature barcode matrices (.h5 files) from GSE202011 were read using scanpy. Sample conditions were inferred from filenames (L = lesional, NL = non-lesional, HM/HF = healthy male/female). For each spot, mean expression was computed across the 51 TOPPLE-classified stabilizer TFs and 51 destabilizer TFs from the full CIMA regulon set. Per-sample metrics included mean stabilizer expression, mean destabilizer expression, and the stabilizer-to-destabilizer ratio. Group comparisons (lesional vs healthy, non-lesional vs healthy) used Mann-Whitney U tests at both the sample and spot levels.

### STRATA-Atlas: Cross-tissue coupling integration

For each of the 61 CIMA L4 cell types, the top 5 marker regulons by mean AUC (from Table S5) were selected and TF names extracted. Mean TF expression per Visium spot was computed as a proxy for cell-type activity, pooled across all 30 samples (24,227 spots). Pairwise Pearson correlations of cell-type activity profiles across spots yielded a 61 x 61 spatial co-expression coupling matrix. Cell types with zero variance were excluded. A Mantel test (9,999 permutations) compared the spatial coupling matrix with the blood-derived eQTL r_b matrix on 59 shared cell types.

### Statistical analysis

All correlations were assessed using Spearman rank correlation. Group comparisons used Mann-Whitney U tests with asymptotic P-value computation. Mantel tests used 9,999 permutations with Pearson correlation of upper-triangle matrix elements. No multiple testing correction was applied to primary cross-method correlations, as each represents a distinct pre-specified hypothesis. Analyses were performed in Python 3.13 using pandas, scipy, numpy, matplotlib, seaborn, scanpy, and openpyxl. All code is available at https://github.com/jengweitjiu/HIRA.

---

## Figure Legends

### Figure 1. HIRA: Multi-scale architectural analysis of the human immune system.

**(A)** TOPPLE stability ranking of the top 15 stabilizer eRegulons by mean redistribution index (RI) across 61 cell types. HSF1_+ is the top stabilizer (RI = 0.0063). Dashed line indicates the position of STAT5B_+ (rank #18, RI = 0.0031).

**(B)** TOPPLE stability landscape heatmap showing redistribution index (z-scored per regulon) for the top 25 stabilizers and bottom 25 destabilizers. Columns represent 61 L4 cell types. Red indicates high RI (removal disrupts landscape); blue indicates low RI (dispensable).

**(C)** SICAI: mean eQTL effect-size correlation (mean r_b) per cell type versus number of disease associations (SMR). Spearman rho = 0.509, P = 9.4 x 10^-6, n = 69 cell types.

**(D)** DGSA: mean genetic non-additivity per cell type versus number of disease associations (SMR). Spearman rho = 0.796, P = 5.2 x 10^-16, n = 68 cell types.

**(E)** IPA: Violin plots of |log2 fold change| (sex difference) for stabilizer, intermediate, and destabilizer regulons. Stabilizers show 11.9-fold lower sex-linked variation (mean 0.028 vs 0.335, Mann-Whitney P = 1.96 x 10^-243).

**(F)** Cross-method Spearman correlation matrix between cell-type-level metrics: TOPPLE RI, SICAI mean r_b, DGSA non-additivity, and disease association count.

**(G)** STRATA: Per-sample stabilizer vs destabilizer mean TF expression across psoriasis Visium samples, colored by condition. Scatter shows near-perfect correlation (r = 0.985) with lesional samples shifted toward higher expression.

**(H)** Distribution of TOPPLE RI values across all 203 regulons, showing the long-tailed stability hierarchy.

### Figure 2. Epigenomic amplification of regulatory non-additivity (Extension 1).

**(A)** Overlaid histograms of non-additivity distributions for eQTL (5,253 genes) and caQTL (15,083 peaks), showing distributional crossover at ~0.89.

**(B)** Quantile-quantile plot comparing eQTL and caQTL non-additivity distributions. Departure from the diagonal confirms distinct distributional shapes.

**(C)** Paired scatter plot of non-additivity at 4,087 matched loci (same variant and cell type). Points above the diagonal indicate caQTL > eQTL. Wilcoxon signed-rank P = 7.0 x 10^-35; 58.4% of pairs favor caQTL.

**(D)** Per-cell-type mean non-additivity versus disease trait count for both eQTL (rho = 0.796) and caQTL (rho = 0.672), demonstrating convergent disease prediction across regulatory layers.

### Figure 3. Multi-layered coupling architecture (Extension 3).

**(A)** Bar chart of per-cell-type mean Jaccard regulatory coupling, showing top and bottom 10 cell types. Overall mean Jaccard = 0.406.

**(B)** Three-way Mantel correlation comparison: eQTL-caQTL (r = 0.482), eQTL-Jaccard (r = 0.388), caQTL-Jaccard (r = 0.566). All P < 0.0001. The stronger caQTL-Jaccard alignment indicates chromatin accessibility operates closer to the regulatory layer.

**(C)** Regulon breadth versus TOPPLE RI scatter plot, colored by stability class. Spearman rho = 0.283, P = 4.3 x 10^-5.

**(D)** Stabilizer vs destabilizer comparison of regulatory complexity: unique target genes (417 vs 175, P = 7.9 x 10^-4) and mean cell-type breadth (29.0 vs 23.1, P = 3.3 x 10^-4).

### Figure 4. Spatial validation of TOPPLE stability classifications (STRATA).

**(A)** Grouped bar chart of mean stabilizer and destabilizer TF expression by condition (healthy, non-lesional, lesional). Stabilizer expression in lesional vs healthy: P = 7.7 x 10^-4.

**(B)** Box plot of stabilizer-to-destabilizer expression ratio by condition. Lesional ratio (3.42) significantly exceeds healthy (2.88); P = 1.1 x 10^-3.

**(C)** Per-sample stabilizer TF expression by condition, showing clear separation between lesional and healthy with non-lesional intermediate.

**(D)** Scatter plot of mean stabilizer vs destabilizer expression per sample, colored by condition (r = 0.985, P = 7.7 x 10^-23), demonstrating that the ratio shift in lesional tissue reflects selective amplification rather than global expression changes.

---

## References

1. Tabula Sapiens Consortium. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. Science 376, eabl4896 (2022).

2. Dominguez Conde C, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022).

3. Suo C, et al. Mapping the developing human immune system across organs. Science 376, eabo0510 (2022).

4. Yin Y, et al. A multi-omics atlas of the human immune system. Science (2026). DOI: 10.1126/science.adt3130.

5. Lopez-Otin C, et al. Hallmarks of aging: an expanding universe. Cell 186, 243-278 (2023).

6. Benayoun BA, et al. Remodeling of epigenome and transcriptome landscapes with aging in mice reveals widespread induction of inflammatory responses. Genome Res 29, 697-709 (2019).

7. Franceschi C, et al. Inflammaging: a new immune-metabolic viewpoint for age-related diseases. Nat Rev Endocrinol 14, 576-590 (2018).

8. Castillo RL, et al. Spatial transcriptomics stratifies psoriasis disease severity by emergent cellular ecosystems. Sci Immunol 8, eabq7991 (2023). GSE202011.
