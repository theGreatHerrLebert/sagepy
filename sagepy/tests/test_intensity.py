"""Tests for the intensity prediction interface layer."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from sagepy.core.intensity import (
    PredictionRequest,
    PredictionResult,
    validate_prediction_result,
    aggregate_predictions_by_peptide,
    write_intensity_file,
    read_intensity_file,
    get_peptide_intensities,
    get_sequence_length,
    create_uniform_intensity,
    ION_KIND_B,
    ION_KIND_Y,
)


class TestPredictionRequest:
    def test_creation(self):
        request = PredictionRequest(
            sequences=np.array(["PEPTIDEK", "ANOTHERPEPTIDE"]),
            charges=np.array([2, 3]),
            peptide_indices=np.array([0, 1]),
        )
        assert len(request) == 2
        assert request.charges.dtype == np.int32
        assert request.peptide_indices.dtype == np.int64

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="Array length mismatch"):
            PredictionRequest(
                sequences=np.array(["PEPTIDEK", "ANOTHER"]),
                charges=np.array([2]),  # Wrong length
                peptide_indices=np.array([0, 1]),
            )

    def test_to_dataframe(self):
        request = PredictionRequest(
            sequences=np.array(["PEPTIDEK", "ANOTHERPEPTIDE"]),
            charges=np.array([2, 3]),
            peptide_indices=np.array([0, 1]),
        )
        df = request.to_dataframe()
        assert list(df.columns) == ["sequence", "charge", "peptide_idx"]
        assert len(df) == 2

    def test_from_dataframe(self):
        import pandas as pd

        df = pd.DataFrame({
            "sequence": ["PEPTIDEK", "ANOTHERPEPTIDE"],
            "charge": [2, 3],
            "peptide_idx": [0, 1],
        })
        request = PredictionRequest.from_dataframe(df)
        assert len(request) == 2
        assert request.sequences[0] == "PEPTIDEK"

    def test_repr(self):
        request = PredictionRequest(
            sequences=np.array(["PEPTIDEK", "PEPTIDEK", "ANOTHER"]),
            charges=np.array([2, 3, 2]),
            peptide_indices=np.array([0, 0, 1]),
        )
        repr_str = repr(request)
        assert "n_entries=3" in repr_str
        assert "unique_peptides=2" in repr_str


class TestPredictionResult:
    def test_creation(self):
        result = PredictionResult(
            peptide_indices=np.array([0, 1]),
            charges=np.array([2, 3]),
            intensities=[
                np.random.rand(2, 7, 2).astype(np.float32),  # 8 AA peptide
                np.random.rand(2, 13, 2).astype(np.float32),  # 14 AA peptide
            ],
        )
        assert len(result) == 2
        assert result.ion_kinds == [ION_KIND_B, ION_KIND_Y]
        assert result.max_fragment_charge == 2

    def test_custom_ion_kinds(self):
        result = PredictionResult(
            peptide_indices=np.array([0]),
            charges=np.array([2]),
            intensities=[np.random.rand(3, 7, 3).astype(np.float32)],
            ion_kinds=[ION_KIND_B, ION_KIND_Y, 0],  # A ion too
            max_fragment_charge=3,
        )
        assert len(result.ion_kinds) == 3
        assert result.max_fragment_charge == 3


class TestGetSequenceLength:
    def test_simple_sequence(self):
        assert get_sequence_length("PEPTIDEK") == 8

    def test_with_unimod(self):
        assert get_sequence_length("PEPTC[UNIMOD:4]IDEK") == 9

    def test_multiple_mods(self):
        assert get_sequence_length("M[UNIMOD:35]PEPTC[UNIMOD:4]IDEK") == 10

    def test_n_terminal_mod(self):
        # N-terminal acetylation
        assert get_sequence_length("[UNIMOD:1]PEPTIDEK") == 8


class TestValidatePredictionResult:
    def test_valid_result(self):
        request = PredictionRequest(
            sequences=np.array(["PEPTIDEK"]),
            charges=np.array([2]),
            peptide_indices=np.array([0]),
        )
        result = PredictionResult(
            peptide_indices=np.array([0]),
            charges=np.array([2]),
            intensities=[np.random.rand(2, 7, 2).astype(np.float32)],
        )
        assert validate_prediction_result(request, result) is True

    def test_length_mismatch(self):
        request = PredictionRequest(
            sequences=np.array(["PEPTIDEK", "ANOTHER"]),
            charges=np.array([2, 3]),
            peptide_indices=np.array([0, 1]),
        )
        result = PredictionResult(
            peptide_indices=np.array([0]),
            charges=np.array([2]),
            intensities=[np.random.rand(2, 7, 2).astype(np.float32)],
        )
        with pytest.raises(ValueError, match="Length mismatch"):
            validate_prediction_result(request, result)

    def test_indices_mismatch(self):
        request = PredictionRequest(
            sequences=np.array(["PEPTIDEK"]),
            charges=np.array([2]),
            peptide_indices=np.array([0]),
        )
        result = PredictionResult(
            peptide_indices=np.array([1]),  # Wrong index
            charges=np.array([2]),
            intensities=[np.random.rand(2, 7, 2).astype(np.float32)],
        )
        with pytest.raises(ValueError, match="Peptide indices don't match"):
            validate_prediction_result(request, result)


class TestAggregatePredictions:
    def test_max_charge_aggregation(self):
        # Same peptide at charges 2 and 3
        result = PredictionResult(
            peptide_indices=np.array([0, 0]),
            charges=np.array([2, 3]),
            intensities=[
                np.full((2, 7, 2), 0.5, dtype=np.float32),  # charge 2
                np.full((2, 7, 2), 0.9, dtype=np.float32),  # charge 3
            ],
        )
        aggregated = aggregate_predictions_by_peptide(result, "max_charge")
        assert len(aggregated) == 1
        assert 0 in aggregated
        # Should have kept the charge 3 prediction (0.9 values)
        assert np.allclose(aggregated[0], 0.9)

    def test_min_charge_aggregation(self):
        result = PredictionResult(
            peptide_indices=np.array([0, 0]),
            charges=np.array([2, 3]),
            intensities=[
                np.full((2, 7, 2), 0.5, dtype=np.float32),  # charge 2
                np.full((2, 7, 2), 0.9, dtype=np.float32),  # charge 3
            ],
        )
        aggregated = aggregate_predictions_by_peptide(result, "min_charge")
        assert np.allclose(aggregated[0], 0.5)

    def test_mean_aggregation(self):
        result = PredictionResult(
            peptide_indices=np.array([0, 0]),
            charges=np.array([2, 3]),
            intensities=[
                np.full((2, 7, 2), 0.4, dtype=np.float32),
                np.full((2, 7, 2), 0.6, dtype=np.float32),
            ],
        )
        aggregated = aggregate_predictions_by_peptide(result, "mean")
        assert np.allclose(aggregated[0], 0.5)

    def test_multiple_peptides(self):
        result = PredictionResult(
            peptide_indices=np.array([0, 1, 0]),
            charges=np.array([2, 2, 3]),
            intensities=[
                np.full((2, 7, 2), 0.1, dtype=np.float32),
                np.full((2, 9, 2), 0.2, dtype=np.float32),
                np.full((2, 7, 2), 0.3, dtype=np.float32),
            ],
        )
        aggregated = aggregate_predictions_by_peptide(result, "max_charge")
        assert len(aggregated) == 2
        assert np.allclose(aggregated[0], 0.3)  # charge 3 for peptide 0
        assert np.allclose(aggregated[1], 0.2)  # only charge 2 for peptide 1


class TestSagiFileIO:
    def test_write_and_read_roundtrip(self):
        # Create test predictions
        predictions = [
            np.random.rand(2, 7, 2).astype(np.float32),  # 8 AA peptide
            np.random.rand(2, 9, 2).astype(np.float32),  # 10 AA peptide
            np.random.rand(2, 5, 2).astype(np.float32),  # 6 AA peptide
        ]
        peptide_lengths = [8, 10, 6]

        with tempfile.NamedTemporaryFile(suffix=".sagi", delete=False) as f:
            path = f.name

        try:
            write_intensity_file(path, predictions, peptide_lengths)
            data = read_intensity_file(path)

            assert data["peptide_count"] == 3
            assert data["max_charge"] == 2
            assert data["ion_kinds"] == [ION_KIND_B, ION_KIND_Y]
            assert len(data["offsets"]) == 3
        finally:
            Path(path).unlink()

    def test_custom_ion_kinds(self):
        predictions = [np.random.rand(3, 7, 3).astype(np.float32)]
        peptide_lengths = [8]

        with tempfile.NamedTemporaryFile(suffix=".sagi", delete=False) as f:
            path = f.name

        try:
            write_intensity_file(
                path, predictions, peptide_lengths,
                max_charge=3,
                ion_kinds=[0, 1, 4],  # A, B, Y
            )
            data = read_intensity_file(path)

            assert data["max_charge"] == 3
            assert data["ion_kinds"] == [0, 1, 4]
        finally:
            Path(path).unlink()

    def test_get_peptide_intensities(self):
        # Create predictable test data
        predictions = [
            np.arange(28, dtype=np.float32).reshape(2, 7, 2),  # 8 AA peptide
            np.arange(36, dtype=np.float32).reshape(2, 9, 2) + 100,  # 10 AA peptide
        ]
        peptide_lengths = [8, 10]

        with tempfile.NamedTemporaryFile(suffix=".sagi", delete=False) as f:
            path = f.name

        try:
            write_intensity_file(path, predictions, peptide_lengths)
            data = read_intensity_file(path)

            # Retrieve first peptide
            retrieved_0 = get_peptide_intensities(data, 0, 8)
            assert retrieved_0.shape == (2, 7, 2)
            assert np.allclose(retrieved_0, predictions[0])

            # Retrieve second peptide
            retrieved_1 = get_peptide_intensities(data, 1, 10)
            assert retrieved_1.shape == (2, 9, 2)
            assert np.allclose(retrieved_1, predictions[1])
        finally:
            Path(path).unlink()

    def test_invalid_magic_raises(self):
        with tempfile.NamedTemporaryFile(suffix=".sagi", delete=False) as f:
            f.write(b"XXXX")  # Invalid magic
            path = f.name

        try:
            with pytest.raises(ValueError, match="Invalid magic"):
                read_intensity_file(path)
        finally:
            Path(path).unlink()


class TestUniformIntensity:
    def test_default_parameters(self):
        uniform = create_uniform_intensity(10)
        assert uniform.shape == (2, 9, 2)  # 2 ion kinds, 9 positions, 2 charges
        assert np.allclose(uniform, 1.0)

    def test_custom_parameters(self):
        uniform = create_uniform_intensity(
            peptide_length=8,
            ion_kinds=[0, 1, 4],  # A, B, Y
            max_charge=3,
            value=0.5,
        )
        assert uniform.shape == (3, 7, 3)
        assert np.allclose(uniform, 0.5)
