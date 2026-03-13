# Extension 1: Epigenomic Amplification of Regulatory Non-Additivity

## Results

To test whether chromatin accessibility QTLs exhibit greater cell-type specificity than
expression QTLs — a phenomenon we term "epigenomic amplification" — we applied the
Decomposed Geometric Signature Analysis (DGSA) framework to 151,875 cis-caQTLs
(52,361 peaks across 42 cell types) and compared them with the 71,530 cis-eQTLs
(5,253 genes in ≥3 cell types across 69 cell types) from the CIMA atlas.

The two QTL classes showed similar mean non-additivity (caQTL: 0.871; eQTL: 0.864),
but their distributions were significantly different (Kolmogorov–Smirnov D = 0.235,
P = 4.8 × 10⁻¹⁹⁰). A quantile–quantile analysis revealed a crossover at approximately
0.89: caQTL non-additivity exceeded eQTL in the lower quantiles (more moderate
cell-type specificity), while eQTL dominated the upper tail (extreme specificity)
(Fig. A–B). This distributional crossing invalidated simple unpaired comparisons;
indeed, a naive Mann–Whitney test paradoxically favored eQTL (P = 4.6 × 10⁻³⁴).

To resolve this confound, we identified 4,087 matched gene–peak pairs sharing the
same lead variant and cell-type context, enabling a paired design that controls for
locus identity. In this matched analysis, caQTL non-additivity significantly exceeded
eQTL non-additivity (Wilcoxon signed-rank P = 7.0 × 10⁻³⁵), with 58.4% of pairs
showing higher caQTL values (mean difference = +0.028) (Fig. C). The paired values
were weakly but significantly correlated (Spearman ρ = 0.183, P = 5.9 × 10⁻³²),
indicating that while genetic and epigenomic specificity are linked, the epigenomic
layer adds substantial independent variation.

Both QTL classes were predictive of cell-type disease burden, as measured by the
number of unique disease traits per cell type from SMR pleiotropic associations
(Table S15). Per-cell-type mean eQTL non-additivity correlated with disease count
at ρ = 0.796 (P = 5.2 × 10⁻¹⁶, n = 68 cell types), while caQTL non-additivity
showed a comparable but attenuated correlation (ρ = 0.672, P = 1.1 × 10⁻⁶, n = 42)
(Fig. D). Cell-type eQTL and caQTL non-additivity were themselves correlated
(ρ = 0.667, P = 1.4 × 10⁻⁶), consistent with a shared architectural principle
operating across regulatory layers.

Together, these results demonstrate that chromatin accessibility effects are
amplified toward greater cell-type specificity compared to their co-regulated
expression effects at the same genetic loci. This epigenomic amplification suggests
that the chromatin landscape acts as a specificity filter: genetic variants exert
broad transcriptional effects that are sharpened at the epigenomic level into more
cell-type-restricted regulatory programs. The convergent disease prediction by both
layers — with eQTL providing stronger signal — implies that disease relevance is
established at the expression level but refined through epigenomic architecture.
