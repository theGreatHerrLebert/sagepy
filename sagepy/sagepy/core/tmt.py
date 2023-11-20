from typing import Optional, List

import sagepy_connector

from sagepy.core import ProcessedSpectrum
from sagepy.core.scoring import Feature
from sagepy.core.spectrum import Peak

psc = sagepy_connector.py_tmt


class Isobaric:
    def __init__(self, type_name: str):
        types = ["tmt6", "tmt10", "tmt11", "tmt16", "tmt18"]
        if type_name in types:
            self.__isobaric_ptr = psc.PyIsobaric(type_name)
        else:
            raise ValueError(f"Invalid isobaric type, allowed values are: {types}")

    @classmethod
    def from_py_isobaric(cls, isobaric: psc.PyIsobaric):
        instance = cls.__new__(cls)
        instance.__isobaric_ptr = isobaric
        return instance

    @property
    def type_name(self):
        return self.__isobaric_ptr.type_name

    def __repr__(self):
        return f"Isobaric({self.__isobaric_ptr.type_name})"

    def get_py_ptr(self):
        return self.__isobaric_ptr

    def modification_mass(self) -> Optional[float]:
        maybe_mass = self.__isobaric_ptr.modification_mass()
        if maybe_mass is None:
            return None
        return maybe_mass


class Purity:
    def __init__(self, ratio: float, correct_precursors: int, incorrect_precursors: int):
        self.__purity_ptr = psc.PyPurity(ratio, correct_precursors, incorrect_precursors)

    @classmethod
    def from_py_purity(cls, purity: psc.PyPurity):
        instance = cls.__new__(cls)
        instance.__purity_ptr = purity
        return instance

    @property
    def ratio(self):
        return self.__purity_ptr.ratio

    @property
    def correct_precursors(self):
        return self.__purity_ptr.correct_precursors

    @property
    def incorrect_precursors(self):
        return self.__purity_ptr.incorrect_precursors

    def __repr__(self):
        return (f"Purity(ratio={self.ratio}, correct_precursors={self.correct_precursors}, "
                f"incorrect_precursors={self.incorrect_precursors})")

    def get_py_ptr(self):
        return self.__purity_ptr


class Quant:
    def __init__(self, hit: Feature, hit_purity: Purity, spectrum: ProcessedSpectrum,
                 chimera: Optional[Feature] = None, chimera_purity: Optional[Purity] = None,
                 intensities: Optional[List[Peak]] = None):

        if chimera is not None:
            chimera = chimera.get_py_ptr()

        if chimera_purity is not None:
            chimera_purity = chimera_purity.get_py_ptr()

        if intensities is not None:
            intensities = [peak.get_py_ptr() for peak in intensities]

        self.__quant_ptr = psc.PyQuant(hit.get_py_ptr(), hit_purity.get_py_ptr(),
                                       spectrum.get_py_ptr(), chimera, chimera_purity, intensities)

    @classmethod
    def from_py_quant(cls, quant: psc.PyQuant):
        instance = cls.__new__(cls)
        instance.__quant_ptr = quant
        return instance

    @property
    def hit(self):
        return Feature.from_py_feature(self.__quant_ptr.hit())

    @property
    def hit_purity(self):
        return Purity.from_py_purity(self.__quant_ptr.hit_purity())

    @property
    def spectrum(self):
        return ProcessedSpectrum.from_py_processed_spectrum(self.__quant_ptr.spectrum())

    @property
    def chimera(self):
        maybe_chimera = self.__quant_ptr.chimera()
        if maybe_chimera is None:
            return None
        return Feature.from_py_feature(maybe_chimera)

    @property
    def chimera_purity(self):
        maybe_chimera_purity = self.__quant_ptr.chimera_purity()
        if maybe_chimera_purity is None:
            return None
        return Purity.from_py_purity(maybe_chimera_purity)

    @property
    def intensities(self):
        intensities = self.__quant_ptr.intensities()
        if intensities is None:
            return None
        return [Peak.from_py_peak(peak) for peak in intensities]

    def __repr__(self):
        return (f"Quant(hit={self.hit}, hit_purity={self.hit_purity}, spectrum={self.spectrum}, "
                f"chimera={self.chimera}, chimera_purity={self.chimera_purity}, intensities={self.intensities})")

    def get_py_ptr(self):
        return self.__quant_ptr
