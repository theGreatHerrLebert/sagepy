"""Tests for sagepy.core.mass module."""
import warnings
import pytest
from sagepy.core.mass import Tolerance, Composition, CONSTANTS, monoisotopic


class TestTolerance:
    """Tests for the Tolerance class."""

    def test_ppm_creation(self):
        """Test creating a tolerance with ppm values."""
        tol = Tolerance(ppm=(-10.0, 10.0))
        assert tol.ppm == (-10.0, 10.0)
        assert tol.da is None

    def test_da_creation(self):
        """Test creating a tolerance with Da values."""
        tol = Tolerance(da=(-0.5, 0.5))
        assert tol.da == (-0.5, 0.5)
        assert tol.ppm is None

    def test_both_raises_error(self):
        """Test that providing both da and ppm raises an error."""
        with pytest.raises(ValueError):
            Tolerance(da=(-0.5, 0.5), ppm=(-10.0, 10.0))

    def test_neither_raises_error(self):
        """Test that providing neither da nor ppm raises an error."""
        with pytest.raises(ValueError):
            Tolerance()

    def test_zero_width_ppm_warns(self):
        """A zero-width ppm window is almost always a mistake (the user
        meant the symmetric form). It should warn but still construct."""
        with pytest.warns(UserWarning, match="zero-width"):
            tol = Tolerance(ppm=(20.0, 20.0))
        assert tol.ppm == (20.0, 20.0)

    def test_zero_width_da_warns(self):
        """Same check for the Da branch."""
        with pytest.warns(UserWarning, match="zero-width"):
            tol = Tolerance(da=(0.5, 0.5))
        assert tol.da == (0.5, 0.5)

    def test_reversed_ppm_raises(self):
        """A reversed (low > high) interval is always a bug."""
        with pytest.raises(ValueError, match="reversed"):
            Tolerance(ppm=(10.0, -10.0))

    def test_reversed_da_raises(self):
        """A reversed Da interval is always a bug."""
        with pytest.raises(ValueError, match="reversed"):
            Tolerance(da=(0.5, -0.5))

    def test_symmetric_no_warning(self):
        """The conventional symmetric usage must NOT warn."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # any warning becomes an exception
            Tolerance(ppm=(-10.0, 10.0))
            Tolerance(da=(-0.5, 0.5))

    def test_asymmetric_no_warning(self):
        """Legitimately asymmetric tolerances (e.g. systematic offset) must
        not warn — only equal/reversed bounds are flagged."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            Tolerance(ppm=(-5.0, 15.0))

    def test_ppm_bounds(self):
        """Test PPM tolerance bounds calculation.

        Based on sage-core test: Tolerance::Ppm(-10.0, 20.0).bounds(1000.0) == (999.99, 1000.02)
        """
        tol = Tolerance(ppm=(-10.0, 20.0))
        lo, hi = tol.bounds(1000.0)
        assert abs(lo - 999.99) < 0.001
        assert abs(hi - 1000.02) < 0.001

    def test_ppm_bounds_smaller_mass(self):
        """Test PPM tolerance bounds with smaller mass.

        Based on sage-core test: Tolerance::Ppm(-10.0, 10.0).bounds(487.0) == (486.99513, 487.00487)
        """
        tol = Tolerance(ppm=(-10.0, 10.0))
        lo, hi = tol.bounds(487.0)
        assert abs(lo - 486.99513) < 0.0001
        assert abs(hi - 487.00487) < 0.0001

    def test_da_bounds(self):
        """Test Da tolerance bounds calculation."""
        tol = Tolerance(da=(-0.5, 0.5))
        lo, hi = tol.bounds(500.0)
        assert abs(lo - 499.5) < 0.001
        assert abs(hi - 500.5) < 0.001

    def test_contains(self):
        """Test tolerance contains method."""
        tol = Tolerance(ppm=(-10.0, 10.0))
        # 1000.0 +/- 10 ppm = 999.99 to 1000.01
        assert tol.contains(1000.0, 1000.005)
        assert tol.contains(1000.0, 999.995)
        assert not tol.contains(1000.0, 1000.02)
        assert not tol.contains(1000.0, 999.97)

    def test_multiplication_float(self):
        """Test tolerance multiplication with float."""
        tol = Tolerance(ppm=(-10.0, 10.0))
        tol2 = tol * 2.0
        assert tol2.ppm == (-20.0, 20.0)

    def test_multiplication_int(self):
        """Test tolerance multiplication with int."""
        tol = Tolerance(da=(-0.5, 0.5))
        tol2 = tol * 2
        assert tol2.da == (-1.0, 1.0)

    def test_repr_ppm(self):
        """Test string representation for ppm tolerance."""
        tol = Tolerance(ppm=(-10.0, 10.0))
        assert "ppm" in repr(tol)

    def test_repr_da(self):
        """Test string representation for Da tolerance."""
        tol = Tolerance(da=(-0.5, 0.5))
        assert "da" in repr(tol)


class TestComposition:
    """Tests for the Composition class."""

    def test_creation(self):
        """Test creating a composition."""
        comp = Composition(carbon=10, sulfur=2)
        assert comp.carbon == 10
        assert comp.sulfur == 2

    def test_aa_composition(self):
        """Test getting amino acid composition."""
        # Methionine contains sulfur
        comp = Composition.aa_composition('M')
        assert comp.sulfur == 1
        assert comp.carbon > 0

    def test_repr(self):
        """Test string representation."""
        comp = Composition(carbon=5, sulfur=1)
        assert "carbon" in repr(comp)
        assert "sulfur" in repr(comp)


class TestConstants:
    """Tests for mass constants."""

    def test_h2o(self):
        """Test H2O constant is in expected range."""
        assert 18.0 < CONSTANTS.H2O < 18.02

    def test_proton(self):
        """Test PROTON constant is in expected range."""
        assert 1.0 < CONSTANTS.PROTON < 1.01

    def test_neutron(self):
        """Test NEUTRON constant is in expected range."""
        assert 1.0 < CONSTANTS.NEUTRON < 1.01

    def test_nh3(self):
        """Test NH3 constant is in expected range."""
        assert 17.0 < CONSTANTS.NH3 < 17.03


class TestMonoisotopic:
    """Tests for the monoisotopic mass function."""

    def test_valid_amino_acid(self):
        """Test getting monoisotopic mass for valid amino acids."""
        # Glycine is the smallest amino acid
        mass_g = monoisotopic('G')
        assert 57.0 < mass_g < 58.0  # Glycine ~57.02

    def test_all_standard_amino_acids(self):
        """Test that all standard amino acids return positive masses."""
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        for aa in amino_acids:
            mass = monoisotopic(aa)
            assert mass > 0, f"Amino acid {aa} should have positive mass"

    def test_invalid_input_lowercase(self):
        """Test that lowercase amino acids raise an error."""
        with pytest.raises(Exception):
            monoisotopic('a')

    def test_invalid_input_multiple_chars(self):
        """Test that multiple characters raise an error."""
        with pytest.raises(Exception):
            monoisotopic('AA')
