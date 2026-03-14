# SICAI — Shared Immune Cellular Architecture Index

## Results

To map the topology of inter-cellular genetic coupling across the immune system, we
constructed a 69 x 69 pairwise eQTL sharing matrix from the CIMA atlas using the
mashr-derived r_b statistic, which quantifies the proportion of shared eQTL effects
between cell-type pairs. The global mean r_b (off-diagonal) was 0.819, consistent with
the published CIMA estimate of 0.82, confirming substantial sharing of genetic regulatory
effects across immune cell types.

Per-cell-type coupling metrics revealed a gradient of genetic integration. The most
coupled cell type was CD4 Th22-like CCR10 (mean r_b = 0.882), while megakaryocytes
(MK GP9) were the most genetically isolated (mean r_b = 0.603). The coefficient of
variation (CV) of r_b across cell-type partners averaged 0.111, and Shannon entropy of
normalized coupling profiles averaged 4.133, indicating relatively uniform coupling
patterns for most cell types. The standard deviation of per-cell-type mean r_b was 0.039,
suggesting that while global coupling is high, meaningful heterogeneity exists across
the immune hierarchy.

We next tested whether genetic coupling predicts disease relevance. Per-cell-type mean
r_b correlated positively with the number of unique disease traits from SMR pleiotropic
associations (Spearman rho = 0.509, P = 9.4 x 10^-6, n = 69 cell types): more
genetically coupled cell types tend to harbor more disease-relevant genetic effects.
Furthermore, a Mantel test comparing the eQTL coupling matrix (r_b) with the caQTL
coupling matrix from chromatin accessibility QTLs yielded a significant correlation
(Mantel r = 0.482, P = 0.0001), indicating that genetic and epigenomic inter-cellular
architectures share a common organizational logic. Together, these results establish
that immune cell types are embedded in a dense, hierarchically organized network of
shared genetic regulation, and that the degree of coupling is itself predictive of
clinical disease burden.
