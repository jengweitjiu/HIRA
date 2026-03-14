Run the TOPPLE analysis (Layer 1: regulon stability via leave-one-out JSD perturbation).

Steps:
1. Run `python src/topple.py`
2. Verify output: results/topple.csv should have 203 regulons
3. Check key values: HSF1_+ RI should be ~0.006322 (rank #1)
4. Verify stability_class column: top 25% = stabilizer, bottom 25% = destabilizer
5. Print top 5 stabilizers and top 5 destabilizers

Input: data/raw/science.adt3130_table_s5.xlsx (sheet: eRegulons_Activators_Exp_AUC_RS)
Output: results/topple.csv, results/topple_stability_scores.csv
