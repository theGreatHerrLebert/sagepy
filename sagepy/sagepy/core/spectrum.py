import numpy as np

from typing import List, Optional

import sagepy_connector
from numpy.typing import NDArray

from sagepy.core.mass import Tolerance
psc = sagepy_connector.py_spectrum


class Peak:
    def __init__(self, mass: float, intensity: float):
        """Peak class

        Args:
            mass (float): The mass of the peak
            intensity (float): The intensity of the peak
        """
        self.__peak_ptr = psc.PyPeak(mass, intensity)

    @classmethod
    def from_py_peak(cls, peak: psc.PyPeak):
        instance = cls.__new__(cls)
        instance.__peak_ptr = peak
        return instance

    @property
    def mass(self):
        return self.__peak_ptr.mass

    @property
    def intensity(self):
        return self.__peak_ptr.intensity

    def get_py_ptr(self):
        return self.__peak_ptr

    def __repr__(self):
        return f"Peak(mass: {self.mass}, intensity: {self.intensity})"


class Deisotoped:
    def __init__(self, mz: float, intensity: float, charge: int = None, envelope: int = None):
        """Deisotoped class

        Args:
            mz (float): The mz of the peak
            intensity (float): The intensity of the peak
            charge (int, optional): The charge of the peak. Defaults to None.
            envelope (int, optional): The envelope of the peak. Defaults to None.
        """
        self.__deisotoped_ptr = psc.PyDeisotoped(mz, intensity, charge, envelope)

    @classmethod
    def from_py_deisotoped(cls, deisotoped: psc.PyDeisotoped):
        instance = cls.__new__(cls)
        instance.__deisotoped_ptr = deisotoped
        return instance

    @property
    def mz(self):
        return self.__deisotoped_ptr.mz

    @property
    def intensity(self):
        return self.__deisotoped_ptr.intensity

    @property
    def charge(self):
        return self.__deisotoped_ptr.charge

    @property
    def envelope(self):
        return self.__deisotoped_ptr.envelope

    def get_py_ptr(self):
        return self.__deisotoped_ptr

    def __repr__(self):
        return f"Deisotoped(mz: {self.mz}, intensity: {self.intensity}, charge: {self.charge}, envelope: {self.envelope})"


class Precursor:
    def __init__(self, mz: float, intensity: float = None, charge: int = None,
                 spectrum_ref: str = None, isolation_window: Tolerance = None,
                 inverse_ion_mobility: float = None, collision_energy: Optional[float] = None):
        """Precursor class

        Args:
            mz (float): The mz of the precursor
            intensity (float, optional): The intensity of the precursor. Defaults to None.
            charge (int, optional): The charge of the precursor. Defaults to None.
            spectrum_ref (str, optional): The spectrum reference of the precursor. Defaults to None.
            isolation_window (Tolerance, optional): The isolation window of the precursor. Defaults to None.
            inverse_ion_mobility (float, optional): The inverse ion mobility of the precursor. Defaults to None.
        """
        if isolation_window is not None:
            self.__precursor_ptr = psc.PyPrecursor(mz, intensity, charge,
                                                   spectrum_ref, isolation_window.get_py_ptr(), inverse_ion_mobility,
                                                   collision_energy)
        else:
            self.__precursor_ptr = psc.PyPrecursor(mz, intensity, charge, spectrum_ref, None, inverse_ion_mobility,
                                                   collision_energy)

    @classmethod
    def from_py_precursor(cls, precursor: psc.PyPrecursor):
        instance = cls.__new__(cls)
        instance.__precursor_ptr = precursor
        return instance

    @property
    def mz(self):
        return self.__precursor_ptr.mz

    def calibrate_mz_ppm(self, ppm: float):
        self.__precursor_ptr.calibrate_mz_ppm(ppm)

    @property
    def intensity(self):
        return self.__precursor_ptr.intensity

    @property
    def charge(self):
        return self.__precursor_ptr.charge

    @property
    def spectrum_ref(self):
        return self.__precursor_ptr.spectrum_ref

    @property
    def isolation_window(self):
        is_window_ptr = self.__precursor_ptr.isolation_window
        if is_window_ptr is not None:
            return Tolerance.from_py_tolerance(is_window_ptr)
        else:
            return None

    @property
    def inverse_ion_mobility(self):
        return self.__precursor_ptr.inverse_ion_mobility

    @property
    def collision_energy(self):
        return self.__precursor_ptr.collision_energy

    def get_py_ptr(self):
        return self.__precursor_ptr

    def __repr__(self):
        return (f"Precursor(mz: {np.round(self.mz, 2)}, "
                f"intensity: {self.intensity}, "
                f"charge: {self.charge}, "
                f"spectrum_ref: {self.spectrum_ref}, "
                f"isolation_window: {self.isolation_window}), "
                f"inverse_ion_mobility: {np.round(self.inverse_ion_mobility, 2) if self.inverse_ion_mobility is not None else self.inverse_ion_mobility}, "
                f"collision_energy: {np.round(self.collision_energy, 2) if self.collision_energy is not None else self.collision_energy})")


class Representation:
    def __init__(self, representation: str = 'centroid'):
        """Representation class

        Args:
            representation (str, optional): The representation of the spectrum. Defaults to 'centroid'.
        """
        self.__representation_ptr = psc.PyRepresentation(representation)

    @classmethod
    def from_py_representation(cls, representation: psc.PyRepresentation):
        instance = cls.__new__(cls)
        instance.__representation_ptr = representation
        return instance

    @property
    def representation(self):
        return self.__representation_ptr.representation_as_string

    def get_py_ptr(self):
        return self.__representation_ptr

    def __repr__(self):
        return f"Representation(representation: {self.representation})"


class ProcessedSpectrum:
    def __init__(self,
                 id: str,
                 level: int,
                 file_id: int,
                 scan_start_time: float,
                 ion_injection_time: float,
                 precursors: List[Precursor],
                 peaks: List[Peak],
                 total_ion_current: float):
        """ProcessedSpectrum class

        Args:
            id (str): The id of the spectrum
            level (int): The level of the spectrum
            file_id (int): The file id of the spectrum
            scan_start_time (float): The scan start time of the spectrum
            ion_injection_time (float): The ion injection time of the spectrum
            precursors (List[Precursor]): The precursors of the spectrum
            peaks (List[Peak]): The peaks of the spectrum
            total_ion_current (float): The total ion current of the spectrum
        """
        self.__processed_spectrum_ptr = psc.PyProcessedSpectrum(
            level, id, file_id, scan_start_time,
            ion_injection_time, [p.get_py_ptr() for p in precursors],
            [p.get_py_ptr() for p in peaks], total_ion_current)

    @classmethod
    def from_py_processed_spectrum(cls, processed_spectrum: psc.PyProcessedSpectrum):
        instance = cls.__new__(cls)
        instance.__processed_spectrum_ptr = processed_spectrum
        return instance

    @property
    def level(self):
        return self.__processed_spectrum_ptr.level

    @property
    def id(self):
        return self.__processed_spectrum_ptr.id

    @property
    def file_id(self):
        return self.__processed_spectrum_ptr.file_id

    @property
    def scan_start_time(self):
        return self.__processed_spectrum_ptr.scan_start_time

    @property
    def ion_injection_time(self):
        return self.__processed_spectrum_ptr.ion_injection_time

    @property
    def precursors(self):
        return [Precursor.from_py_precursor(p) for p in self.__processed_spectrum_ptr.precursors]

    @property
    def peaks(self):
        return [Peak.from_py_peak(p) for p in self.__processed_spectrum_ptr.peaks]

    @property
    def total_ion_current(self):
        return self.__processed_spectrum_ptr.total_ion_current

    @property
    def collision_energies(self):
        return self.__processed_spectrum_ptr.collision_energies

    def calibrate_mz_ppm(self, ppm: float, calibrate_precursors: bool = True):
        if calibrate_precursors:
            for precursor in self.precursors:
                precursor.calibrate_mz_ppm(ppm)
        self.__processed_spectrum_ptr.calibrate_mz_ppm(ppm)

    def get_py_ptr(self):
        return self.__processed_spectrum_ptr

    def __repr__(self):
        return (f"ProcessedSpectrum(level: {self.level}, "
                f"id: {self.id}, "
                f"file_id: {self.file_id}, "
                f"scan_start_time: {np.round(self.scan_start_time, 2)}, "
                f"ion_injection_time: {np.round(self.ion_injection_time, 2)}, "
                f"num_precursors: {len(self.precursors)}, "
                f"num_peaks: {len(self.peaks)}, "
                f"collision_energies: {self.collision_energies},"
                f"total_ion_current: {self.total_ion_current})")


class RawSpectrum:
    def __init__(self,
                 file_id: int,
                 spec_id: str,
                 total_ion_current: float,
                 precursors: List[Precursor],
                 mz: NDArray[float],
                 intensity: NDArray[float],
                 representation: Representation = Representation(),
                 scan_start_time: float = 0.0,
                 ion_injection_time: float = 0.0,
                 ms_level: int = 2,
                 ):
        """RawSpectrum class

        Args:
            file_id (int): The file id of the spectrum
            spec_id (str): The id of the spectrum
            total_ion_current (float): The total ion current of the spectrum
            precursors (List[Precursor]): The precursors of the spectrum
            mz (NDArray[float]): The mz of the peaks of the spectrum
            intensity (NDArray[float]): The intensity of the peaks of the spectrum
            representation (Representation, optional): The representation of the spectrum. Defaults to Representation().
            scan_start_time (float, optional): The scan start time of the spectrum. Defaults to 0.0.
            ion_injection_time (float, optional): The ion injection time of the spectrum. Defaults to 0.0.
            ms_level (int, optional): The ms level of the spectrum. Defaults to 2.
        """
        self.__raw_spectrum_ptr = psc.PyRawSpectrum(file_id, ms_level, spec_id, [p.get_py_ptr() for p in precursors],
                                                    representation.get_py_ptr(),
                                                    scan_start_time, ion_injection_time, total_ion_current,
                                                    mz, intensity)
    @classmethod
    def from_py_raw_spectrum(cls, raw_spectrum: psc.PyRawSpectrum):
        instance = cls.__new__(cls)
        instance.__raw_spectrum_ptr = raw_spectrum
        return instance

    @property
    def ms_level(self):
        return self.__raw_spectrum_ptr.ms_level

    @property
    def id(self):
        return self.__raw_spectrum_ptr.id

    @property
    def precursors(self):
        return [Precursor.from_py_precursor(p) for p in self.__raw_spectrum_ptr.precursors]

    @property
    def representation(self):
        return Representation.from_py_representation(self.__raw_spectrum_ptr.representation)

    @property
    def scan_start_time(self):
        return self.__raw_spectrum_ptr.scan_start_time

    @property
    def ion_injection_time(self):
        return self.__raw_spectrum_ptr.ion_injection_time

    @property
    def total_ion_current(self):
        return self.__raw_spectrum_ptr.total_ion_current

    @property
    def mz(self) -> NDArray[float]:
        return self.__raw_spectrum_ptr.mz

    @property
    def intensity(self) -> NDArray[float]:
        return self.__raw_spectrum_ptr.intensity

    def get_py_ptr(self):
        return self.__raw_spectrum_ptr

    def filter_top_n_peaks(self, n: int) -> 'RawSpectrum':
        return RawSpectrum.from_py_raw_spectrum(self.__raw_spectrum_ptr.filter_top_n_peaks(n))

    def __repr__(self):
        return (f"RawSpectrum(ms_level: {self.ms_level}, "
                f"id: {self.id}, "
                f"precursors: {self.precursors}, "
                f"representation: {self.representation}, "
                f"scan_start_time: {self.scan_start_time}, "
                f"ion_injection_time: {self.ion_injection_time}, "
                f"total_ion_current: {self.total_ion_current}, "
                f"num_peaks: {len(self.mz)})")


class SpectrumProcessor:
    def __init__(
            self, take_top_n: int = 150,
            min_fragment_mz: float = 150,
            max_fragment_mz: float = 2000,
            min_deisotope_mz: float = 0.0,
            deisotope: bool = False,
    ):
        """SpectrumProcessor class

        Args:
            take_top_n (int, optional): The number of peaks to take. Defaults to 150.
            min_fragment_mz (float, optional): The minimum fragment mz. Defaults to 150.
            max_fragment_mz (float, optional): The maximum fragment mz. Defaults to 2000.
            min_deisotope_mz (float, optional): The minimum deisotope mz. Defaults to 0.0.
            deisotope (bool, optional): Whether to deisotope the spectrum. Defaults to False.
        """
        self.__spectrum_processor_ptr = psc.PySpectrumProcessor(
            take_top_n, max_fragment_mz, min_fragment_mz, min_deisotope_mz, deisotope)

    @classmethod
    def from_py_spectrum_processor(cls, spectrum_processor: psc.PySpectrumProcessor):
        instance = cls.__new__(cls)
        instance.__spectrum_processor_ptr = spectrum_processor
        return instance

    def get_py_ptr(self):
        return self.__spectrum_processor_ptr

    @property
    def take_top_n(self):
        return self.__spectrum_processor_ptr.take_top_n

    @property
    def max_fragment_mz(self):
        return self.__spectrum_processor_ptr.max_fragment_mz

    @property
    def min_fragment_mz(self):
        return self.__spectrum_processor_ptr.min_fragment_mz

    @property
    def deisotope(self):
        return self.__spectrum_processor_ptr.deisotope

    def __repr__(self):
        return f"SpectrumProcessor(take_top_n: {self.take_top_n}, max_fragment_mz: {self.max_fragment_mz}, " \
               f"min_fragment_mz: {self.min_fragment_mz}, deisotope: {self.deisotope})"

    def process(self, raw_spectrum: RawSpectrum) -> ProcessedSpectrum:
        return ProcessedSpectrum.from_py_processed_spectrum(self.__spectrum_processor_ptr.process(raw_spectrum.get_py_ptr()))
