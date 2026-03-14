Run all validation tests for the HIRA pipeline.

Steps:
1. Run `python -m pytest tests/ -v`
2. Check that all 45 tests pass across 10 test classes
3. Verify key expected values:
   - TOPPLE: HSF1 RI ~0.0063 (rank #1), STAT5B RI ~0.0031 (rank #18)
   - DGSA: mean non-additivity ~0.867
   - SICAI: mean r_b ~0.818, disease rho ~0.353
   - IPA: sex perturbation P ~4.4e-218
   - Extension #1: caQTL non-additivity > eQTL (paired Wilcoxon)
4. If any test fails, read the test file and the source script to diagnose

Test files: tests/test_validation.py, tests/conftest.py
Config: pytest.ini
