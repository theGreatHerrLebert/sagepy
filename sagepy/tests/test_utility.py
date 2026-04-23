"""Tests for sagepy utility functions."""
import pytest
import numpy as np
from sagepy.utility import calculate_ppm_error, calculate_ppms, mean_ppm, median_ppm

# Import directly from sagepy_connector for tests that don't need full utility module
import sagepy_connector
psc = sagepy_connector.py_utility


def prosit_intensities_to_fragments_map(intensities):
    """Convert Prosit intensity array to fragment map."""
    return psc.flat_prosit_array_to_fragments_map(list(intensities))


def sage_sequence_to_unimod(sequence, modifications, expected_modifications):
    """Convert sage sequence to UNIMOD format."""
    return psc.sage_sequence_to_unimod(sequence, list(modifications), expected_modifications)


class TestPpmCalculations:
    """Tests for PPM error calculation functions."""

    def test_calculate_ppm_error_zero(self):
        """Test PPM error is zero when values are equal."""
        error = calculate_ppm_error(1000.0, 1000.0)
        assert error == pytest.approx(0.0)

    def test_calculate_ppm_error_positive(self):
        """Test positive PPM error."""
        # 1000.001 vs 1000.0 should be 1 ppm
        error = calculate_ppm_error(1000.001, 1000.0)
        assert error == pytest.approx(1.0, rel=0.01)

    def test_calculate_ppm_error_negative(self):
        """Test negative PPM error."""
        error = calculate_ppm_error(999.999, 1000.0)
        assert error == pytest.approx(-1.0, rel=0.01)

    def test_calculate_ppm_error_larger_mass(self):
        """Test PPM error with larger mass values."""
        # At 2000 Da, 0.002 difference should still be 1 ppm
        error = calculate_ppm_error(2000.002, 2000.0)
        assert error == pytest.approx(1.0, rel=0.01)

    def test_calculate_ppms_array(self):
        """Test PPM calculation for arrays."""
        measured = np.array([1000.001, 500.0005, 2000.002])
        reference = np.array([1000.0, 500.0, 2000.0])
        errors = calculate_ppms(measured, reference)
        assert len(errors) == 3
        assert errors[0] == pytest.approx(1.0, rel=0.01)
        assert errors[1] == pytest.approx(1.0, rel=0.01)
        assert errors[2] == pytest.approx(1.0, rel=0.01)

    def test_mean_ppm(self):
        """Test mean PPM calculation."""
        measured = np.array([1000.001, 1000.002, 1000.003])
        reference = np.array([1000.0, 1000.0, 1000.0])
        mean = mean_ppm(measured, reference)
        assert mean == pytest.approx(2.0, rel=0.01)  # Average of 1, 2, 3

    def test_median_ppm(self):
        """Test median PPM calculation."""
        measured = np.array([1000.001, 1000.002, 1000.003])
        reference = np.array([1000.0, 1000.0, 1000.0])
        median = median_ppm(measured, reference)
        assert median == pytest.approx(2.0, rel=0.01)  # Median of 1, 2, 3


class TestPrositIntensities:
    """Tests for Prosit intensity conversion functions."""

    def test_prosit_intensities_to_fragments_map_shape(self):
        """Test that prosit intensities are converted to the correct shape."""
        # Prosit outputs 174 values (29 * 2 * 3)
        intensities = np.zeros(174, dtype=np.float32)
        intensities[0] = 1.0  # Set one intensity

        fragment_map = prosit_intensities_to_fragments_map(intensities)

        # Should return a dictionary
        assert isinstance(fragment_map, dict)

    def test_prosit_intensities_keys(self):
        """Test that keys are (ion_type, charge, ordinal) tuples."""
        intensities = np.ones(174, dtype=np.float32)
        fragment_map = prosit_intensities_to_fragments_map(intensities)

        # Check that keys are tuples of length 3
        for key in fragment_map.keys():
            assert len(key) == 3
            ion_type, charge, ordinal = key
            assert ion_type in (0, 1)  # B or Y
            assert 1 <= charge <= 3
            assert 1 <= ordinal <= 29


class TestSageSequenceToUnimod:
    """Tests for sage sequence to UNIMOD conversion."""

    def test_unmodified_sequence(self):
        """Test that unmodified sequence is returned unchanged."""
        result = sage_sequence_to_unimod('PEPTIDE', [0.0] * 7, set())
        assert result == 'PEPTIDE'

    def test_with_known_modification(self):
        """Test sequence with known UNIMOD modification."""
        mods = [0.0, 0.0, 0.0, 15.9949, 0.0, 0.0, 0.0]  # Oxidation on position 4
        expected_mods = {'[UNIMOD:35]'}  # Oxidation
        result = sage_sequence_to_unimod('PEPTMDE', mods, expected_mods)
        assert 'UNIMOD:35' in result

    def test_with_unknown_modification(self):
        """Test sequence with unknown modification shows delta mass."""
        mods = [0.0, 0.0, 123.456, 0.0, 0.0, 0.0, 0.0]  # Unknown mod
        result = sage_sequence_to_unimod('PEPTIDE', mods, set())
        # Should contain the delta mass
        assert '123' in result or '+123' in result

    def test_empty_expected_mods(self):
        """Test with empty expected modifications set."""
        mods = [0.0, 15.9949, 0.0]
        result = sage_sequence_to_unimod('PEP', mods, set())
        # Should show delta mass since no expected mods match
        assert '+' in result or '15.99' in result

    def test_sequence_length_matches_mods(self):
        """Test that sequence length must match modifications length."""
        # The Rust code panics with mismatched lengths, which raises BaseException
        with pytest.raises(BaseException):
            sage_sequence_to_unimod('PEPTIDE', [0.0, 0.0], set())
