# DGSA — Decomposed Geometric Signature Analysis

## Results

To characterize the cell-type specificity of genetic regulatory effects, we applied
Decomposed Geometric Signature Analysis (DGSA) to 71,530 cis-eQTLs from the CIMA atlas.
For each of 5,253 genes detected in three or more cell types, we constructed a
69-dimensional slope vector (one entry per cell type, zero-filled where absent) and
computed the Gini coefficient of squared slopes as a non-additivity metric — where
higher values indicate more cell-type-restricted genetic effects.

Mean non-additivity across all genes was 0.864 (median = 0.911; s.d. = 0.126), indicating
that the vast majority of eQTL effects are concentrated in a small subset of the 69 cell
types rather than distributed uniformly. The complementary Gini coefficient of absolute
slopes was 0.842, and mean sparsity (fraction of cell types with zero slope) was 0.819,
reflecting that a typical gene harbors a detectable eQTL in only ~12.5 of 69 cell types.
Genes with eQTLs in exactly 3 cell types showed near-maximal non-additivity (>0.97),
while genes detected across many cell types exhibited progressively lower values,
consistent with broader genetic effects being distributed more evenly across the immune
system.

To assess whether cell-type-specific genetic regulation predicts disease relevance, we
computed per-cell-type mean non-additivity and correlated it with the number of unique
disease traits identified through SMR pleiotropic associations (Table S15). This
correlation was strongly positive (Spearman rho = 0.796, P = 5.2 x 10^-16, n = 68
cell types), indicating that cell types with more cell-type-specific eQTL architectures
are disproportionately enriched for disease-associated genetic effects. This finding
positions non-additivity as a geometric predictor of immune disease relevance.
