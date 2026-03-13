# Multi-Scale Architectural Analysis of the Human Immune System
## Revealed by Geometric Decomposition of the CIMA Atlas

### Framing for Genome Biology / Nature Communications

The Chinese Immune Multi-Omics Atlas (CIMA; Yin et al., Science 2026) represents
a landmark cataloging effort -- 10 million cells from 428 individuals across 73 cell
types, with paired single-cell RNA-seq and ATAC-seq enabling the most comprehensive
map of human immune regulation to date. However, cataloging is not understanding.
CIMA's 203 eRegulons, 223,405 xQTLs, and 2,085 disease associations remain
uninterpreted as architectural entities. What stabilizes the regulatory landscape?
Which regulons are load-bearing, and which are dispensable? How does a cell type's
coupling and effect-size architecture determine its disease vulnerability?

We address these questions by applying three novel mathematical frameworks to CIMA's
published supplementary data:

**TOPPLE stability analysis** of 203 eRegulons across 61
cell types reveals a striking hierarchy: HSF1_+ emerges as the top stabilizer, with
51 regulons classified as structurally critical (top quartile
redistribution index). Notably, STAT5B ¡X the #1 stabilizer in our prior psoriasis analysis (GSE173706, RI 0.03 -> 0.61) ¡X ranks #18 in healthy blood, consistent with disease-specific amplification of regulatory importance.

**SICAI coupling analysis** of 4,692 eQTL sharing pairs across 69 cell types
demonstrates that mean eQTL effect-size correlation (mean r_b) significantly
predicts disease pleiotropy (Spearman rho = 0.353, P = 3.13e-03,
n = 68 cell types). Shannon entropy of the r_b profile exhibits
a ceiling effect in healthy blood (CV = 0.009), while mean r_b captures the
biologically relevant variance (CV = 0.070).

**DGSA geometric decomposition** of 5,253 cis-eQTL genes reveals
that genetic non-additivity -- a measure of cell-type-specific effect-size geometry
-- is the strongest predictor of disease pleiotropy across all three methods
(Spearman rho = 0.690, P = 7.43e-11, n = 68 cell
types). Cell types whose eQTL effects are concentrated in specific directions of
effect-size space (high non-additivity) carry substantially more disease associations
than those with uniform, shared effects.

Cross-method validation confirms complementary architectural axes: DGSA non-additivity
and SICAI coupling strength capture distinct but convergent aspects of disease
vulnerability, while TOPPLE stability identifies the regulons whose perturbation
most disrupts cell-type identity.

### Key Statistics Summary

| Method | Metric vs Disease | rho | P-value | n |
|--------|-------------------|-----|---------|---|
| SICAI  | Mean r_b          | 0.353 | 3.13e-03 | 68 |
| DGSA   | Non-additivity    | 0.690 | 7.43e-11 | 68 |

Our framework is computationally lightweight -- all analyses use only CIMA's published
supplementary tables, requiring no raw data access or high-performance computing. This
establishes a paradigm where mathematical frameworks, applied to community-generated
atlases, can extract fundamentally new biological insights from existing data.

**Keywords:** regulatory architecture, stability landscape, genetic non-additivity,
coupling strength, single-cell transcriptomics, immune atlas, CIMA, eQTL geometry
