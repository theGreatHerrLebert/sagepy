"""Tests for sagepy.core.ml.pep."""

import numpy as np
import pytest

from sagepy.core.ml.pep import calculate_pep, calculate_pep_single, posterior_error


def test_calculate_pep_single_shape_and_metadata():
    scores = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    decoys = np.array([True, True, False, False], dtype=bool)

    pep_bins, min_score, score_step = calculate_pep_single(scores, decoys, bins=8)

    assert pep_bins.shape == (8,)
    assert min_score == pytest.approx(1.0)
    assert score_step == pytest.approx(3.0 / 7.0)
    assert np.all(np.diff(pep_bins) <= 1e-12)
    assert pep_bins[0] > pep_bins[-1]


def test_posterior_error_interpolates_between_bins():
    scores = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    decoys = np.array([True, True, False, False], dtype=bool)
    pep_bins, min_score, score_step = calculate_pep_single(scores, decoys, bins=8)

    pep = posterior_error(pep_bins, min_score, score_step, 2.5)
    assert pep == pytest.approx(0.5)


def test_calculate_pep_matches_known_values():
    scores = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    decoys = np.array([True, True, False, False], dtype=bool)

    pep = calculate_pep(scores, decoys, bins=8)

    expected = np.array(
        [0.9999251184, 0.8860137784, 0.1139862216, 0.0000748816],
        dtype=np.float64,
    )
    assert np.allclose(pep, expected, rtol=1e-6, atol=1e-8)


def test_calculate_pep_single_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        calculate_pep_single(np.array([], dtype=np.float64), np.array([], dtype=bool))

    with pytest.raises(ValueError):
        calculate_pep_single(
            np.array([1.0, 2.0], dtype=np.float64),
            np.array([True], dtype=bool),
        )

    with pytest.raises(ValueError):
        calculate_pep_single(
            np.array([1.0, 2.0], dtype=np.float64),
            np.array([True, False], dtype=bool),
            bins=1,
        )


def test_calculate_pep_single_requires_targets_and_decoys():
    scores = np.array([1.0, 2.0, 3.0], dtype=np.float64)

    with pytest.raises(ValueError):
        calculate_pep_single(scores, np.array([True, True, True], dtype=bool))

    with pytest.raises(ValueError):
        calculate_pep_single(scores, np.array([False, False, False], dtype=bool))
