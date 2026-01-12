"""Tests for sagepy.core.ion_series module."""
import pytest
from sagepy.core.ion_series import IonType, Ion, IonSeries
from sagepy.core.peptide import Peptide
from sagepy.core.enzyme import Position


class TestIonType:
    """Tests for the IonType class."""

    def test_valid_ion_types(self):
        """Test creating valid ion types."""
        valid_types = ['a', 'b', 'c', 'x', 'y', 'z']
        for ion_type in valid_types:
            it = IonType(ion_type)
            assert it.to_string().lower() == ion_type.lower()

    def test_case_insensitive(self):
        """Test that ion types are case insensitive."""
        it_lower = IonType('b')
        it_upper = IonType('B')
        assert it_lower.to_string() == it_upper.to_string()

    def test_invalid_ion_type(self):
        """Test that invalid ion types raise an error."""
        with pytest.raises(ValueError):
            IonType('invalid')

    def test_y_factory(self):
        """Test the y() factory method."""
        y = IonType.y()
        assert y.to_string().lower() == 'y'

    def test_b_factory(self):
        """Test the b() factory method."""
        b = IonType.b()
        assert b.to_string().lower() == 'b'

    def test_equality(self):
        """Test ion type equality."""
        b1 = IonType('b')
        b2 = IonType('b')
        y = IonType('y')
        assert b1 == b2
        assert b1 != y

    def test_hash(self):
        """Test that ion types are hashable."""
        b = IonType('b')
        y = IonType('y')
        ion_set = {b, y}
        assert len(ion_set) == 2

    def test_repr(self):
        """Test string representation."""
        b = IonType('b')
        assert 'B' in repr(b) or 'b' in repr(b).lower()


class TestIon:
    """Tests for the Ion class."""

    def test_creation(self):
        """Test creating an ion."""
        ion_type = IonType('b')
        ion = Ion(ion_type, mass=500.0)
        assert ion.mono_isotopic_mass == 500.0

    def test_ion_type_property(self):
        """Test getting ion type from ion."""
        ion_type = IonType('y')
        ion = Ion(ion_type, mass=600.0)
        assert ion.ion_type.to_string().lower() == 'y'

    def test_repr(self):
        """Test string representation."""
        ion_type = IonType('b')
        ion = Ion(ion_type, mass=500.0)
        repr_str = repr(ion)
        assert '500' in repr_str


def make_peptide(sequence: str) -> Peptide:
    """Helper to create a simple peptide for testing."""
    return Peptide(
        decoy=False,
        sequence=sequence,
        modifications=[0.0] * len(sequence),
        mono_isotopic=sum(57.0 for _ in sequence),  # Approximate
        missed_cleavages=0,
        position=Position('internal'),
        proteins=['TEST_PROTEIN']
    )


class TestIonSeries:
    """Tests for the IonSeries class."""

    def test_b_ion_series_creation(self):
        """Test creating a b-ion series."""
        peptide = make_peptide('PEPTIDE')
        b_type = IonType.b()
        series = IonSeries(peptide, b_type)
        assert series.ion_type.to_string().lower() == 'b'

    def test_y_ion_series_creation(self):
        """Test creating a y-ion series."""
        peptide = make_peptide('PEPTIDE')
        y_type = IonType.y()
        series = IonSeries(peptide, y_type)
        assert series.ion_type.to_string().lower() == 'y'

    def test_ion_series_length(self):
        """Test that ion series has correct length (n-1 fragments)."""
        peptide = make_peptide('PEPTIDE')  # 7 residues
        b_type = IonType.b()
        series = IonSeries(peptide, b_type)
        ions = series.get_ion_series()
        assert len(ions) == 6  # n-1 = 6

    def test_b_ion_masses_increase(self):
        """Test that b-ion masses increase along the series."""
        peptide = make_peptide('PEPTIDE')
        b_type = IonType.b()
        series = IonSeries(peptide, b_type)
        ions = series.get_ion_series()

        for i in range(1, len(ions)):
            assert ions[i].mono_isotopic_mass > ions[i - 1].mono_isotopic_mass, \
                "B-ion masses should increase"

    def test_y_ion_masses_decrease(self):
        """Test that y-ion masses decrease along the series."""
        peptide = make_peptide('PEPTIDE')
        y_type = IonType.y()
        series = IonSeries(peptide, y_type)
        ions = series.get_ion_series()

        for i in range(1, len(ions)):
            assert ions[i].mono_isotopic_mass < ions[i - 1].mono_isotopic_mass, \
                "Y-ion masses should decrease"

    def test_peptide_property(self):
        """Test getting peptide from ion series."""
        peptide = make_peptide('PEPTIDE')
        b_type = IonType.b()
        series = IonSeries(peptide, b_type)
        assert series.peptide.sequence == 'PEPTIDE'

    def test_cumulative_mass_property(self):
        """Test cumulative mass property exists."""
        peptide = make_peptide('PEPTIDE')
        b_type = IonType.b()
        series = IonSeries(peptide, b_type)
        # Just verify it's a number
        assert isinstance(series.cumulative_mass, float)

    def test_repr(self):
        """Test string representation."""
        peptide = make_peptide('PEP')
        b_type = IonType.b()
        series = IonSeries(peptide, b_type)
        repr_str = repr(series)
        assert 'IonSeries' in repr_str
