Run the SICAI analysis (Layer 3: inter-cellular coupling topology via r_b).

Steps:
1. Run `python src/sicai.py`
2. Verify output: results/sicai_rb.csv should have 69 cell types with mean_rb, cv, entropy
3. Check: mean r_b across all cell types should be ~0.818 (matches CIMA published 0.82)
4. Verify disease correlation: rho ~0.353 (P = 0.003)
5. Print top 5 most coupled and top 5 least coupled cell types

Input: data/raw/science.adt3130_table_s8.xlsx (sheet: cis_eQTL, columns: reference_cell_type, query_celltype, rb)
Output: results/sicai_rb.csv, results/sicai_rb_matrix.csv, results/sicai_coupling_complexity.csv
