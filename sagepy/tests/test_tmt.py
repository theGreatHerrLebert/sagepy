"""Tests for sagepy.core.tmt."""

import pytest

import sagepy.core.tmt as tmt


class _PtrWrapper:
    def __init__(self, ptr):
        self._ptr = ptr

    def get_py_ptr(self):
        return self._ptr


def test_isobaric_types_and_modification_masses():
    expected = {
        "tmt6": 229.16293334960938,
        "tmt10": 229.16293334960938,
        "tmt11": 229.16293334960938,
        "tmt16": 304.20709228515625,
        "tmt18": 304.2135009765625,
    }

    for type_name, mass in expected.items():
        isobaric = tmt.Isobaric(type_name)
        assert isobaric.type_name == type_name
        assert isobaric.modification_mass() == pytest.approx(mass)


def test_isobaric_rejects_invalid_type():
    with pytest.raises(ValueError):
        tmt.Isobaric("nope")


def test_purity_properties_and_repr():
    purity = tmt.Purity(0.8, 2, 1)

    assert purity.ratio == pytest.approx(0.8)
    assert purity.correct_precursors == 2
    assert purity.incorrect_precursors == 1
    assert "Purity(" in repr(purity)


def test_quant_wraps_backend_objects(monkeypatch):
    class _FakePyQuant:
        def __init__(self, hit, hit_purity, spectrum, chimera, chimera_purity, intensities):
            self._hit = hit
            self._hit_purity = hit_purity
            self._spectrum = spectrum
            self._chimera = chimera
            self._chimera_purity = chimera_purity
            self._intensities = intensities

        def hit(self):
            return self._hit

        def hit_purity(self):
            return self._hit_purity

        def spectrum(self):
            return self._spectrum

        def chimera(self):
            return self._chimera

        def chimera_purity(self):
            return self._chimera_purity

        def intensities(self):
            return self._intensities

    class _DummyBackend:
        PyQuant = _FakePyQuant

    monkeypatch.setattr(tmt, "psc", _DummyBackend)
    monkeypatch.setattr(
        tmt.Feature,
        "from_py_feature",
        staticmethod(lambda value: ("feature", value)),
    )
    monkeypatch.setattr(
        tmt.ProcessedSpectrum,
        "from_py_processed_spectrum",
        staticmethod(lambda value: ("spectrum", value)),
    )
    monkeypatch.setattr(
        tmt.Purity,
        "from_py_purity",
        staticmethod(lambda value: ("purity", value)),
    )
    monkeypatch.setattr(
        tmt.Peak,
        "from_py_peak",
        staticmethod(lambda value: ("peak", value)),
    )

    quant = tmt.Quant(
        hit=_PtrWrapper("hit"),
        hit_purity=_PtrWrapper("hit_purity"),
        spectrum=_PtrWrapper("spectrum"),
        chimera=_PtrWrapper("chimera"),
        chimera_purity=_PtrWrapper("chimera_purity"),
        intensities=[_PtrWrapper("i1"), _PtrWrapper("i2")],
    )

    assert quant.hit == ("feature", "hit")
    assert quant.hit_purity == ("purity", "hit_purity")
    assert quant.spectrum == ("spectrum", "spectrum")
    assert quant.chimera == ("feature", "chimera")
    assert quant.chimera_purity == ("purity", "chimera_purity")
    assert quant.intensities == [("peak", "i1"), ("peak", "i2")]
    assert "Quant(" in repr(quant)
