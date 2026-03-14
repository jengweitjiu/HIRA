Run Extension #3: Enhancer-Driven Coupling Network analysis.

Steps:
1. Run `python src/ext3_enhancer_coupling.py`
   - Builds regulatory coupling from S4 GRN peak-gene linkages
   - Computes Jaccard index between cell types
   - Three-way Mantel test: eQTL r_b vs caQTL r_b vs regulatory Jaccard
2. Run `python src/ext3_ad_overlay.py`
   - Overlays AD GWAS loci onto S4 Region-Gene pairs (500kb window)
   - Classifies by cell-type accessibility (S3) and regulon membership
   - Generates 4-panel figure: fig_ext3_ad_overlay.pdf
3. Verify outputs:
   - results/ext3_regulatory_coupling.csv
   - results/ext3_regulon_breadth.csv
   - results/ext3_ad_peak_gene_overlay.csv
   - figures/fig_ext3_ad_overlay.pdf

Input: data/raw/CIMA_Table_S4.csv, data/raw/CIMA_Table_S3.csv, data/raw/science.adt3130_table_s8.xlsx
