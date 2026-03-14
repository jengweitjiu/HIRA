# Extension 3: Multi-Layered Coupling Architecture Revealed by Enhancer-Gene Regulatory Networks

## Results

To determine whether the inter-cellular coupling observed at the genetic (eQTL) and
epigenomic (caQTL) levels reflects a shared regulatory logic, we constructed a
cell-type regulatory coupling matrix from the CIMA enhancer-gene regulatory network
(SCENIC+ GRN). We extracted 133,574 unique region-gene linkages spanning 84,625
enhancer regions and 13,645 target genes (Table S4), and mapped each linkage to the
cell types in which its enhancer region was accessible (Table S3; 338,036 peaks across
65 cell types). For each pair of cell types, we computed the Jaccard index of shared
active region-gene pairs, yielding a 65 x 65 regulatory coupling matrix (mean
Jaccard = 0.406).

Three-way Mantel tests revealed a hierarchical coupling architecture in which all
three layers -- genetic, epigenomic, and regulatory -- were significantly correlated
(all P < 0.0001 by permutation). Critically, the caQTL coupling matrix was more
strongly correlated with regulatory coupling (Mantel r = 0.566) than the eQTL
coupling matrix was (r = 0.388), consistent with chromatin accessibility QTLs
operating mechanistically closer to the enhancer-gene regulatory layer. The previously
established eQTL-caQTL coupling correlation (r = 0.482) occupied an intermediate
position, suggesting that genetic effects on expression are partially mediated through
chromatin accessibility.

To connect this regulatory architecture to transcriptional stability, we computed
per-regulon regulatory breadth -- the mean number of cell types in which a regulon's
enhancer regions are accessible. Regulatory breadth correlated positively with TOPPLE
redistribution index (Spearman rho = 0.283, P = 4.3 x 10^-5), indicating that
structurally important regulons maintain broader regulatory footprints. TOPPLE-classified
stabilizers controlled 2.4-fold more target genes than destabilizers (417 vs 175,
Mann-Whitney P = 7.9 x 10^-4) and had enhancers active across more cell types (29.0
vs 23.1, P = 3.3 x 10^-4). This convergence between regulatory complexity and
architectural stability suggests that the immune system's most structurally load-bearing
transcription factors are those embedded in the broadest enhancer-gene networks --
a design principle that may buffer against cell-type-specific perturbations by
distributing regulatory control across a wider chromatin landscape.
