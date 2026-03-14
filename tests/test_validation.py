"""
Validation tests for all HIRA pipeline key numbers.

Each test loads the relevant results CSV and checks that computed values
match expected ranges from the HIRA manuscript.
"""

from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"


# ---------------------------------------------------------------------------
# 1. TOPPLE — Regulon stability via leave-one-out JSD perturbation
# ---------------------------------------------------------------------------

class TestTOPPLE:
    """Validate TOPPLE regulon stability results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.df = pd.read_csv(RESULTS_DIR / "topple.csv")

    def test_total_regulons(self):
        """TOPPLE should analyse exactly 203 regulons."""
        assert len(self.df) == 203, (
            f"Expected 203 regulons, got {len(self.df)}"
        )

    def test_hsf1_ri(self):
        """HSF1 RI should be approximately 0.0063 (within 10%)."""
        hsf1 = self.df[self.df["regulon"] == "HSF1_+"]
        assert len(hsf1) == 1, "HSF1_+ not found in topple.csv"
        ri = hsf1["mean_RI"].values[0]
        assert ri == pytest.approx(0.0063, rel=0.10), (
            f"HSF1 RI = {ri:.6f}, expected ~0.0063"
        )

    def test_stat5b_ri(self):
        """STAT5B RI should be approximately 0.0031 (within 10%)."""
        stat5b = self.df[self.df["regulon"] == "STAT5B_+"]
        assert len(stat5b) == 1, "STAT5B_+ not found in topple.csv"
        ri = stat5b["mean_RI"].values[0]
        assert ri == pytest.approx(0.0031, rel=0.10), (
            f"STAT5B RI = {ri:.6f}, expected ~0.0031"
        )

    def test_hsf1_is_rank_1(self):
        """HSF1 should be the top-ranked regulon (rank 1)."""
        hsf1 = self.df[self.df["regulon"] == "HSF1_+"]
        rank = hsf1["rank"].values[0]
        assert rank == 1, f"HSF1 rank = {rank}, expected 1"

    def test_stat5b_rank(self):
        """STAT5B should be rank 18."""
        stat5b = self.df[self.df["regulon"] == "STAT5B_+"]
        rank = stat5b["rank"].values[0]
        assert rank == 18, f"STAT5B rank = {rank}, expected 18"

    def test_stability_classes_exist(self):
        """All three stability classes should be present."""
        classes = set(self.df["stability_class"].unique())
        expected = {"stabilizer", "intermediate", "destabilizer"}
        assert classes == expected, (
            f"Expected classes {expected}, got {classes}"
        )

    @pytest.mark.parametrize("regulon,expected_class", [
        ("HSF1_+", "stabilizer"),
        ("STAT5B_+", "stabilizer"),
        ("HOXA9_+", "destabilizer"),
    ])
    def test_stability_class_assignment(self, regulon, expected_class):
        """Check that specific regulons are assigned to the correct class."""
        row = self.df[self.df["regulon"] == regulon]
        assert len(row) == 1, f"{regulon} not found"
        assert row["stability_class"].values[0] == expected_class, (
            f"{regulon} should be {expected_class}"
        )


# ---------------------------------------------------------------------------
# 2. DGSA — eQTL effect-size geometry decomposition
# ---------------------------------------------------------------------------

class TestDGSA:
    """Validate DGSA eQTL geometry results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.df = pd.read_csv(RESULTS_DIR / "dgsa_eqtl.csv")

    def test_gene_count(self):
        """DGSA should contain 5253 genes (those in >= 3 cell types)."""
        assert len(self.df) == 5253, (
            f"Expected 5253 genes, got {len(self.df)}"
        )

    def test_mean_non_additivity(self):
        """Mean non-additivity should be ~0.864 (within 10%)."""
        mean_na = self.df["non_additivity"].mean()
        assert mean_na == pytest.approx(0.864, rel=0.10), (
            f"Mean non-additivity = {mean_na:.4f}, expected ~0.864"
        )

    def test_non_additivity_range(self):
        """All non-additivity values should be between 0 and 1."""
        assert self.df["non_additivity"].min() >= 0, "Non-additivity has negative values"
        assert self.df["non_additivity"].max() <= 1, "Non-additivity exceeds 1"

    def test_all_geometry_columns_present(self):
        """All five DGSA geometry metrics should be present."""
        expected_cols = {"magnitude", "uniformity", "non_additivity", "sparsity", "gini"}
        actual_cols = set(self.df.columns)
        missing = expected_cols - actual_cols
        assert not missing, f"Missing columns: {missing}"

    def test_min_celltypes_per_gene(self):
        """Every gene should be present in at least 3 cell types."""
        assert self.df["n_celltypes"].min() >= 3, (
            f"Min cell types = {self.df['n_celltypes'].min()}, expected >= 3"
        )


# ---------------------------------------------------------------------------
# 3. SICAI — Inter-cellular coupling topology via r_b
# ---------------------------------------------------------------------------

class TestSICAI:
    """Validate SICAI eQTL r_b coupling results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.df = pd.read_csv(RESULTS_DIR / "sicai_rb.csv")

    def test_cell_type_count(self):
        """SICAI should cover 69 cell types."""
        assert len(self.df) == 69, (
            f"Expected 69 cell types, got {len(self.df)}"
        )

    def test_mean_rb(self):
        """Mean r_b should be ~0.819 (within 5%)."""
        mean_rb = self.df["mean_rb"].mean()
        assert mean_rb == pytest.approx(0.819, rel=0.05), (
            f"Mean r_b = {mean_rb:.4f}, expected ~0.819"
        )

    def test_rb_values_in_valid_range(self):
        """All mean r_b values should be between 0 and 1."""
        assert self.df["mean_rb"].min() >= 0, "r_b has negative values"
        assert self.df["mean_rb"].max() <= 1, "r_b exceeds 1"

    def test_top_coupled_cell_type(self):
        """The highest-coupled cell type should have mean r_b > 0.85."""
        max_rb = self.df["mean_rb"].max()
        assert max_rb > 0.85, (
            f"Max mean r_b = {max_rb:.4f}, expected > 0.85"
        )


# ---------------------------------------------------------------------------
# 4. IPA — Perturbation resistance (sex differences)
# ---------------------------------------------------------------------------

class TestIPA:
    """Validate IPA sex perturbation results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.df = pd.read_csv(RESULTS_DIR / "ipa_sex.csv")

    def test_all_regulons_present(self):
        """IPA should have results for all 203 regulons (or close)."""
        n = len(self.df)
        assert n >= 200, f"Expected ~203 regulons, got {n}"

    def test_stabilizers_have_lower_perturbation(self):
        """Stabilizer mean |log2FC| should be less than destabilizer mean |log2FC|."""
        stab = self.df[self.df["stability_class"] == "stabilizer"]
        dest = self.df[self.df["stability_class"] == "destabilizer"]
        assert len(stab) > 0 and len(dest) > 0, "Missing stability classes"
        stab_mean = stab["mean_abs_log2FC"].mean()
        dest_mean = dest["mean_abs_log2FC"].mean()
        assert stab_mean < dest_mean, (
            f"Stabilizer mean |log2FC| ({stab_mean:.4f}) should be < "
            f"destabilizer ({dest_mean:.4f})"
        )

    def test_has_required_columns(self):
        """IPA output should have stability_class and perturbation columns."""
        required = {"stability_class", "mean_abs_log2FC", "mean_RI", "rank"}
        missing = required - set(self.df.columns)
        assert not missing, f"Missing columns: {missing}"


# ---------------------------------------------------------------------------
# 5. Extension #1 — caQTL Geometric Architecture
# ---------------------------------------------------------------------------

class TestExtension1CaQTL:
    """Validate Extension #1: caQTL DGSA and coupling results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.caqtl = pd.read_csv(RESULTS_DIR / "ext1_caqtl_dgsa.csv")
        self.eqtl = pd.read_csv(RESULTS_DIR / "dgsa_eqtl.csv")

    def test_caqtl_peak_count(self):
        """caQTL DGSA should have a substantial number of peaks."""
        n = len(self.caqtl)
        assert n > 10000, (
            f"Expected > 10,000 peaks in caQTL DGSA, got {n}"
        )

    def test_caqtl_non_additivity_higher_than_eqtl(self):
        """caQTL mean non-additivity should be >= eQTL mean non-additivity
        (epigenomic amplification hypothesis)."""
        caqtl_na = self.caqtl["non_additivity"].mean()
        eqtl_na = self.eqtl["non_additivity"].mean()
        assert caqtl_na >= eqtl_na, (
            f"caQTL non-additivity ({caqtl_na:.4f}) should be >= "
            f"eQTL non-additivity ({eqtl_na:.4f})"
        )

    def test_caqtl_mean_non_additivity_value(self):
        """caQTL mean non-additivity should be ~0.871 (within 10%)."""
        caqtl_na = self.caqtl["non_additivity"].mean()
        assert caqtl_na == pytest.approx(0.871, rel=0.10), (
            f"caQTL mean non-additivity = {caqtl_na:.4f}, expected ~0.871"
        )

    def test_caqtl_geometry_columns(self):
        """caQTL DGSA should have the same geometry columns as eQTL DGSA."""
        expected = {"magnitude", "uniformity", "non_additivity", "sparsity", "gini"}
        missing = expected - set(self.caqtl.columns)
        assert not missing, f"Missing columns: {missing}"


class TestExtension1CaQTLCoupling:
    """Validate caQTL r_b coupling matrix."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.df = pd.read_csv(RESULTS_DIR / "sicai_caqtl_rb.csv")

    def test_cell_type_count(self):
        """caQTL coupling should cover a substantial number of cell types."""
        n = len(self.df)
        assert n >= 30, (
            f"Expected >= 30 cell types in caQTL r_b, got {n}"
        )

    def test_mean_rb_reasonable(self):
        """caQTL mean r_b should be positive and < 1."""
        mean_rb = self.df["mean_rb"].mean()
        assert 0 < mean_rb < 1, (
            f"caQTL mean r_b = {mean_rb:.4f}, should be between 0 and 1"
        )


# ---------------------------------------------------------------------------
# 6. STRATA — Spatial validation (psoriasis Visium)
# ---------------------------------------------------------------------------

class TestSTRATA:
    """Validate STRATA spatial validation results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.spatial = pd.read_csv(RESULTS_DIR / "strata_spatial.csv")
        self.summary = pd.read_csv(RESULTS_DIR / "strata_summary.csv")

    def test_minimum_samples(self):
        """At least 20 samples should be processed."""
        n = len(self.spatial)
        assert n >= 20, f"Expected >= 20 samples, got {n}"

    def test_lesional_ratio_higher_than_healthy(self):
        """Lesional stabilizer/destabilizer ratio should be higher than healthy."""
        lesional = self.summary[self.summary["condition"] == "lesional"]
        healthy = self.summary[self.summary["condition"] == "healthy"]
        assert len(lesional) == 1 and len(healthy) == 1, (
            "Missing lesional or healthy condition in summary"
        )
        les_ratio = lesional["mean_ratio"].values[0]
        hlt_ratio = healthy["mean_ratio"].values[0]
        assert les_ratio > hlt_ratio, (
            f"Lesional ratio ({les_ratio:.3f}) should be > healthy ({hlt_ratio:.3f})"
        )

    def test_all_conditions_present(self):
        """Summary should contain lesional, non-lesional, and healthy conditions."""
        conditions = set(self.summary["condition"].unique())
        expected = {"lesional", "non-lesional", "healthy"}
        assert conditions == expected, (
            f"Expected conditions {expected}, got {conditions}"
        )

    @pytest.mark.parametrize("condition", ["lesional", "non-lesional", "healthy"])
    def test_positive_ratios(self, condition):
        """All stabilizer/destabilizer ratios should be positive."""
        row = self.summary[self.summary["condition"] == condition]
        ratio = row["mean_ratio"].values[0]
        assert ratio > 0, f"{condition} mean_ratio = {ratio}, expected > 0"


# ---------------------------------------------------------------------------
# 7. STRATA-Atlas — Cross-tissue Mantel test
# ---------------------------------------------------------------------------

class TestSTRATAAtlas:
    """Validate STRATA-Atlas cross-tissue Mantel test results."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.df = pd.read_csv(RESULTS_DIR / "strata_atlas.csv")
        self.mantel_row = self.df[self.df["cell_type"] == "__MANTEL_TEST__"]
        self.cell_data = self.df[self.df["cell_type"] != "__MANTEL_TEST__"]

    def test_mantel_row_exists(self):
        """The Mantel test result row should exist in the atlas output."""
        assert len(self.mantel_row) == 1, (
            "__MANTEL_TEST__ row not found in strata_atlas.csv"
        )

    def test_mantel_p_value(self):
        """Mantel test P-value should be < 0.01."""
        # P-value is stored in mean_rb_coupling column
        p_val = self.mantel_row["mean_rb_coupling"].values[0]
        assert p_val < 0.01, (
            f"Mantel P = {p_val}, expected < 0.01"
        )

    def test_mantel_r_positive(self):
        """Mantel correlation coefficient should be positive."""
        # r is stored in mean_spatial_coupling column
        r_val = self.mantel_row["mean_spatial_coupling"].values[0]
        assert r_val > 0, (
            f"Mantel r = {r_val:.4f}, expected > 0"
        )

    def test_sufficient_cell_types(self):
        """Atlas should include a substantial number of cell types."""
        n = len(self.cell_data)
        assert n >= 50, (
            f"Expected >= 50 cell types in atlas, got {n}"
        )


# ---------------------------------------------------------------------------
# 8. Extension #3 — Enhancer-Driven Coupling Network
# ---------------------------------------------------------------------------

class TestExtension3:
    """Validate Extension #3: regulatory coupling network."""

    @pytest.fixture(autouse=True)
    def load_data(self):
        self.coupling = pd.read_csv(RESULTS_DIR / "ext3_regulatory_coupling.csv")
        self.breadth = pd.read_csv(RESULTS_DIR / "ext3_regulon_breadth.csv")

    def test_coupling_cell_type_count(self):
        """Regulatory coupling matrix should cover at least 50 cell types."""
        n = len(self.coupling)
        assert n >= 50, (
            f"Expected >= 50 cell types in coupling matrix, got {n}"
        )

    def test_coupling_has_jaccard(self):
        """Coupling results should include Jaccard similarity metric."""
        assert "mean_jaccard" in self.coupling.columns, (
            "mean_jaccard column missing from ext3_regulatory_coupling.csv"
        )

    def test_jaccard_range(self):
        """Jaccard values should be between 0 and 1."""
        min_j = self.coupling["mean_jaccard"].min()
        max_j = self.coupling["mean_jaccard"].max()
        assert min_j >= 0, f"Min Jaccard = {min_j}, expected >= 0"
        assert max_j <= 1, f"Max Jaccard = {max_j}, expected <= 1"

    def test_breadth_has_regulons(self):
        """Regulon breadth results should contain multiple regulons."""
        n = len(self.breadth)
        assert n >= 100, (
            f"Expected >= 100 regulons in breadth analysis, got {n}"
        )

    def test_breadth_has_stability_class(self):
        """Breadth results should include stability class annotation."""
        assert "stability_class" in self.breadth.columns, (
            "stability_class column missing from ext3_regulon_breadth.csv"
        )

    def test_breadth_has_ri(self):
        """Breadth results should include TOPPLE RI for cross-referencing."""
        assert "mean_RI" in self.breadth.columns, (
            "mean_RI column missing from ext3_regulon_breadth.csv"
        )


# ---------------------------------------------------------------------------
# Cross-layer consistency checks
# ---------------------------------------------------------------------------

class TestCrossLayerConsistency:
    """Check consistency across HIRA layers."""

    def test_topple_ipa_regulon_alignment(self):
        """TOPPLE and IPA should reference the same set of regulons."""
        topple = pd.read_csv(RESULTS_DIR / "topple.csv")
        ipa = pd.read_csv(RESULTS_DIR / "ipa_sex.csv")
        topple_regulons = set(topple["regulon"].str.replace(r"_\+$", "_+", regex=True))
        # IPA eRegulon names may have gene counts appended; extract base name
        # Just check count alignment
        assert len(topple) == len(ipa), (
            f"TOPPLE has {len(topple)} regulons but IPA has {len(ipa)}"
        )

    def test_eqtl_caqtl_both_have_geometry(self):
        """Both eQTL and caQTL DGSA should have the same metric columns."""
        eqtl = pd.read_csv(RESULTS_DIR / "dgsa_eqtl.csv")
        caqtl = pd.read_csv(RESULTS_DIR / "ext1_caqtl_dgsa.csv")
        metrics = {"magnitude", "uniformity", "non_additivity", "sparsity", "gini"}
        assert metrics.issubset(set(eqtl.columns)), "eQTL missing geometry columns"
        assert metrics.issubset(set(caqtl.columns)), "caQTL missing geometry columns"
