# STRATA: Spatial Transcriptomic Validation of Regulatory Stability Architecture

## Results

To test whether TOPPLE stability classifications, derived from blood immune cells,
generalize to an independent tissue context, we analyzed spatial transcriptomic data
from psoriasis skin biopsies (GSE202011; 10x Visium, 30 samples comprising 14
lesional, 9 non-lesional, and 7 healthy specimens totaling 24,227 spots). For each
Visium spot, we computed mean expression of 51 TOPPLE-classified stabilizer TFs and
51 destabilizer TFs.

Stabilizer TFs were consistently more highly expressed than destabilizer TFs across
all conditions (global mean 0.158 vs 0.048, Mann-Whitney P < 1e-300). Critically,
this asymmetry was amplified in disease: the stabilizer-to-destabilizer expression
ratio was significantly elevated in lesional psoriasis compared to healthy skin
(3.42 vs 2.88, Mann-Whitney P = 1.1 x 10^-3). At the sample level, lesional
specimens showed 3.3-fold higher stabilizer expression than healthy controls
(0.257 vs 0.078, P = 7.7 x 10^-4), whereas destabilizer expression increased
only 2.7-fold (0.074 vs 0.027, P = 1.6 x 10^-3). This disproportionate
stabilizer upregulation in disease is consistent with the TOPPLE framework's
prediction that architecturally load-bearing regulons are preferentially recruited
during inflammatory perturbation.

Non-lesional skin from the same patients showed intermediate stabilizer expression
that did not significantly differ from healthy controls (spot-level P = 0.068),
confirming that the stabilizer enrichment is specific to active disease rather
than a systemic patient-level confounder. The near-perfect correlation between
stabilizer and destabilizer expression across samples (r = 0.985, P = 7.7 x 10^-23)
indicates that both classes scale with overall transcriptional activity, but the
ratio shift in lesional tissue reveals a selective amplification of the stability
architecture under inflammatory stress.

## STRATA-Atlas: Cross-Tissue Coupling Conservation

To test whether inter-cellular coupling topology is conserved across tissues, we
projected CIMA regulon signatures onto 24,227 Visium spots from 30 psoriasis skin
biopsies (GSE202011). For each of the 61 CIMA cell types, we selected the top 5
marker regulons by mean AUC and computed mean TF expression per spot as a proxy
for cell-type activity. Pairwise Pearson correlations of these activity profiles
across all pooled spots yielded a 61 x 61 spatial co-expression coupling matrix.

A Mantel test comparing the spatial coupling matrix with the blood-derived eQTL
sharing matrix (r_b) revealed a significant positive correlation (Mantel r = 0.143,
P = 0.0001 by 9,999 permutations; 59 shared cell types, 1,696 pairwise comparisons).
Although the effect size is moderate -- consistent with the expected divergence between
blood immune cells and bulk skin tissue, and with the inherent resolution limits of
Visium spots -- the significance confirms that the coupling architecture discovered
in blood is not tissue-specific but reflects a systemic organizational principle
detectable even in an independent tissue context and experimental platform.
