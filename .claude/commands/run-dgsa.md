Run the DGSA analysis (Layer 2: eQTL effect-size geometry decomposition).

Steps:
1. Run `python src/dgsa.py`
2. Verify output: results/dgsa_eqtl.csv should have ~5,253 genes (those in >=3 cell types)
3. Check metrics: non_additivity, magnitude, uniformity, sparsity, gini per gene
4. Verify mean non-additivity ~0.867 (HIRA manuscript value)
5. Print top 10 genes by non-additivity

Input: data/raw/CIMA_Table_S6.csv (filter analysis == 'cis-eQTL')
Output: results/dgsa_eqtl.csv, results/dgsa_gene_scores.csv
