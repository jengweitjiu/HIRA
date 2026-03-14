# TOPPLE — Transcriptional Orchestration via Perturbation-Probed Landscape Entropy

## Results

To quantify the structural importance of individual regulons within the immune
transcriptional landscape, we developed TOPPLE, a leave-one-out perturbation framework
applied to the 203 eRegulon activators across 61 cell types from the CIMA atlas.
For each regulon, we removed it from the L1-normalized AUC profile of every cell type,
renormalized the remaining distribution, and computed the Jensen-Shannon divergence (JSD)
between the original and perturbed profiles. The mean JSD across all cell types defines
the redistribution index (RI), which captures how much a regulon's removal destabilizes
the global transcriptional architecture.

RI values spanned two orders of magnitude (range: 1.76 x 10^-5 to 6.32 x 10^-3; mean =
0.00171; median = 0.00154), indicating a highly skewed distribution of structural
importance. HSF1 ranked first (RI = 0.00632), followed by EGR1 (RI = 0.00574) and KLF9
(RI = 0.00569). STAT5B ranked 18th (RI = 0.00306). We classified regulons into three
stability tiers using quartile boundaries (Q75 = 0.00224; Q25 = 0.000927): 51 stabilizers
(top quartile), 101 intermediate, and 51 destabilizers (bottom quartile).

Stabilizers were enriched for broadly expressed transcription factors, particularly
the AP-1 family (JUN, rank 5; JUNB, rank 4; FOSB, rank 10) and ubiquitous regulators
(EGR1, KLF9, KLF2). In contrast, destabilizers were predominantly lineage-restricted
factors: GATA1 (rank 201, RI = 3.55 x 10^-5), GATA2 (rank 200, RI = 3.64 x 10^-5),
and FOXP3 (rank 192, RI = 1.50 x 10^-4). This pattern suggests that the immune
transcriptional network derives its structural integrity from a small set of broadly
deployed regulatory hubs, while lineage-defining factors — despite their biological
importance — contribute minimally to global network stability.
