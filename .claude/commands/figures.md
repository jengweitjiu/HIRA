Regenerate all HIRA publication figures.

Steps:
1. Core HIRA figures:
   - `python src/fig_topple.py` -> figures/fig_topple.pdf (if exists)
   - `python src/fig_ext3_coupling.py` -> figures/fig_ext3_coupling.pdf
   - `python src/fig_strata.py` -> figures/fig_strata.pdf
2. DECODE-AD figures:
   - `python src/fig_decode_ad.py` -> figures/fig_decode_ad.pdf
   - `python src/fig_decode_ad_extended.py` -> figures/fig_decode_ad_extended.pdf
   - `python src/fig_decode_ad_popgen.py` -> figures/fig_decode_ad_popgen.pdf (analysis 4 only)
   - `python src/fig_study_design.py` -> figures/fig_study_design.pdf + fig_circos_summary.pdf
   - `python src/ext3_ad_overlay.py` -> figures/fig_ext3_ad_overlay.pdf (analysis D only)
3. Verify all PDFs exist in figures/ directory
4. All figures: 300 DPI, Nature-style, Arial font, 180mm wide

Output directory: figures/
Format: .pdf (vector) + .png (raster, 300 DPI)
