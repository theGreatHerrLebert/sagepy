"""Tests for native spectrum I/O helpers."""
import pandas as pd

from sagepy import utility


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
