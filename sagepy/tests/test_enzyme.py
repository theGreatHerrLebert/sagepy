"""Tests for sagepy.core.enzyme module."""
import pytest
from sagepy.core.enzyme import Position, EnzymeParameters


class TestPosition:
    """Tests for the Position class."""

    def test_valid_positions(self):
        """Test creating valid positions."""
        # Position returns formatted strings like 'Cterm', 'Nterm', 'Full', 'Internal'
        valid_positions = ['c_term', 'n_term', 'full', 'internal']
        expected = ['Cterm', 'Nterm', 'Full', 'Internal']
        for pos, exp in zip(valid_positions, expected):
            p = Position(pos)
            assert p.position == exp

    def test_default_is_internal(self):
        """Test that default position is internal."""
        p = Position()
        assert p.position == 'Internal'

    def test_invalid_position(self):
        """Test that invalid position raises an error."""
        with pytest.raises(ValueError):
            Position('invalid')

    def test_repr(self):
        """Test string representation."""
        p = Position('n_term')
        assert 'Nterm' in repr(p)


class TestEnzymeParameters:
    """Tests for the EnzymeParameters class."""

    def test_defaults(self):
        """Test default parameters."""
        params = EnzymeParameters()
        assert params.missed_cleavages == 0
        assert params.min_len == 5
        assert params.max_len == 50

    def test_custom_parameters(self):
        """Test custom parameters."""
        params = EnzymeParameters(
            missed_cleavages=2,
            min_len=7,
            max_len=30
        )
        assert params.missed_cleavages == 2
        assert params.min_len == 7
        assert params.max_len == 30

    def test_enzyme_property(self):
        """Test enzyme property is accessible through EnzymeParameters."""
        params = EnzymeParameters()
        # Verify enzyme is accessible through params (may be None for default)
        enzyme = params.enzyme
        assert enzyme is not None

    def test_digest_method(self):
        """Test digest method."""
        params = EnzymeParameters(missed_cleavages=0, min_len=3, max_len=20)
        digests = params.digest('PEPTIDEKPEPTIDER', 'TEST_PROTEIN')
        assert len(digests) > 0
        for d in digests:
            assert d.protein == 'TEST_PROTEIN'

    def test_cleavage_sites(self):
        """Test cleavage sites method."""
        params = EnzymeParameters()
        sites = params.cleavage_sites('PEPTIDEKR')
        assert len(sites) > 0

    def test_repr(self):
        """Test string representation."""
        params = EnzymeParameters(missed_cleavages=1, min_len=6)
        repr_str = repr(params)
        assert 'missed_cleavages' in repr_str
