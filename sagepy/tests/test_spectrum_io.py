"""Tests for native spectrum I/O helpers."""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sagepy import utility
from sagepy.core.spectrum import read_spectra


DATA_DIR = Path(__file__).resolve().parent / "data"


class _FakePrecursor:
    def __init__(self, mz, charge, intensity, collision_energy):
        self.mz = mz
        self.charge = charge
        self.intensity = intensity
        self.collision_energy = collision_energy


class _FakeSpectrum:
    def __init__(self, spec_id, precursors):
        self.id = spec_id
        self.precursors = precursors
        self.scan_start_time = 12.5
        self.ion_injection_time = 8.0
        self.total_ion_current = 1234.0


class _FakeSpectrumProcessor:
    def __init__(self, take_top_n, min_deisotope_mz, deisotope):
        self.take_top_n = take_top_n
        self.min_deisotope_mz = min_deisotope_mz
        self.deisotope = deisotope

    def process(self, spectrum):
        return {"processed": spectrum.id}


def test_extract_ms_data_uses_native_reader(monkeypatch):
    spectra = [
        _FakeSpectrum(
            "scan=1",
            [_FakePrecursor(500.2, 2, 1000.0, 27.5)],
        )
    ]

    monkeypatch.setattr(utility, "read_spectra", lambda **_: spectra)
    monkeypatch.setattr(utility, "SpectrumProcessor", _FakeSpectrumProcessor)

    result = utility.extract_ms_data("example.mzML")

    assert isinstance(result, pd.DataFrame)
    assert list(result["spec_id"]) == ["scan=1"]
    assert list(result["precursor_mz"]) == [500.2]
    assert list(result["collision_energy"]) == [27.5]
    assert list(result["processed_spec"]) == [{"processed": "scan=1"}]


def test_extract_mgf_data_delegates_to_generic_loader(monkeypatch):
    monkeypatch.setattr(
        utility,
        "extract_ms_data",
        lambda file_path, **kwargs: {"file_path": file_path, "kwargs": kwargs},
    )

    result = utility.extract_mgf_data("example.mgf", take_top_n_peaks=50)

    assert result == {
        "file_path": "example.mgf",
        "kwargs": {"take_top_n_peaks": 50},
    }


def test_matteo_first10_mgf_and_pmsms_are_equivalent():
    fixture_dir = DATA_DIR / "matteo_first10"
    mgf_spectra = read_spectra(str(fixture_dir / "first10_corresponding.mgf"))
    pmsms_spectra = read_spectra(
        str(fixture_dir / "first10_sage_input_mgf_equivalent.pmsms")
    )

    assert len(mgf_spectra) == 10
    assert len(pmsms_spectra) == len(mgf_spectra)

    for mgf, pmsms in zip(mgf_spectra, pmsms_spectra):
        mgf_precursor_idx = re.search(r"precursor_idx=(\d+)", mgf.id).group(1)
        assert pmsms.id == f"precursor_idx={mgf_precursor_idx}"
        assert mgf.ms_level == pmsms.ms_level == 2
        assert mgf.total_ion_current == pmsms.total_ion_current
        assert mgf.scan_start_time == pytest.approx(pmsms.scan_start_time, abs=1e-5)

        assert len(mgf.precursors) == len(pmsms.precursors) == 1
        mgf_precursor = mgf.precursors[0]
        pmsms_precursor = pmsms.precursors[0]
        assert mgf_precursor.charge == pmsms_precursor.charge
        assert mgf_precursor.mz == pytest.approx(pmsms_precursor.mz, abs=5e-4)

        assert len(mgf.mz) == len(pmsms.mz)
        assert len(mgf.intensity) == len(pmsms.intensity)
        np.testing.assert_allclose(mgf.mz, pmsms.mz, rtol=0, atol=1e-6)
        np.testing.assert_allclose(mgf.intensity, pmsms.intensity, rtol=0, atol=1e-6)
