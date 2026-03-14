# IPA — Immune Perturbation Architecture

## Results

To determine whether TOPPLE-identified stabilizers are functionally resistant to
biological perturbation, we examined sex-based differential eRegulon activity across
61 cell types using the CIMA sex-difference results (|Log2 Fold Change| between male
and female donors). We merged these data with TOPPLE stability classifications to
compare perturbation magnitude between the 51 stabilizer and 51 destabilizer regulons.

Stabilizers exhibited dramatically lower sex-based perturbation than destabilizers.
The mean |log2FC| per regulon was 0.0282 for stabilizers versus 0.3348 for
destabilizers — a 11.9-fold difference. A Mann-Whitney U test comparing all individual
regulon-by-cell-type |log2FC| values confirmed this difference was highly significant
(P = 1.96 x 10^-243), with stabilizer values consistently smaller than destabilizer
values. The rank-biserial effect size was large and negative, indicating that
stabilizers occupy the lower ranks of the perturbation distribution. Even at the
per-regulon level, the relationship was striking: the Spearman correlation between
redistribution index (RI) and mean |log2FC| was rho = -0.574 (P = 3.71 x 10^-19),
demonstrating that regulons with higher structural importance show proportionally
smaller sex-dependent fold changes.

Among the most perturbation-resistant stabilizers were PATZ1 (mean |log2FC| = 0.0069),
PHF1 (0.0070), and KLF9 (0.0098), while extreme destabilizers such as GATA1 (1.209),
ZBTB8A (1.258), and TCF7L1 (1.660) showed order-of-magnitude larger sex differences.
These findings establish that structural importance in the transcriptional network is
functionally coupled to perturbation resistance: the regulons whose removal most
destabilizes the global architecture are precisely those that remain most invariant
under biological perturbation, consistent with a design principle in which
load-bearing transcriptional elements are buffered against environmental variation.
