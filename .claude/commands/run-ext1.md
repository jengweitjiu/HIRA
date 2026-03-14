Run Extension #1: caQTL Geometric Architecture analysis.

Steps:
1. Run `python src/ext1_caqtl_dgsa.py`
   - Applies DGSA to 151,875 cis-caQTLs from S6 caQTL sheet
   - Filters to peaks in >=3 cell types -> ~15,083 peaks
   - KEY TEST: Mann-Whitney comparing caQTL vs eQTL non-additivity
   - Expected: caQTL non-additivity > eQTL non-additivity (epigenomic amplification)
2. Run `python src/ext1_remaining.py`
   - Per-cell-type caQTL vs eQTL non-additivity correlation
   - SICAI metrics for caQTL r_b coupling (S8 cis_caQTL sheet, 1,722 rows, 42 cell types)
3. Verify outputs:
   - results/ext1_caqtl_dgsa.csv
   - results/ext1_celltype_correlation.csv
   - results/ext1_sicai_caqtl_metrics.csv

Input: data/raw/CIMA_Table_S6.csv, data/raw/science.adt3130_table_s8.xlsx
