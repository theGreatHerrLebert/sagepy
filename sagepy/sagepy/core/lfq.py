from typing import List, Tuple, Dict

from sagepy.core import Feature, Psm, IndexedDatabase, ProcessedSpectrum
from sagepy.core.database import PeptideIx
import sagepy_connector

from sagepy.core.ml.retention_alignment import Alignment
from sagepy.core.spectrum import ProcessedIMSpectrum

psc = sagepy_connector.py_lfq


class Peak:
    def __init__(self, rt: int, spectral_angle: float, score: float, q_value: float):
        self.__peak_ptr = psc.PyPeak(rt, spectral_angle, score, q_value)

    @classmethod
    def from_py_ptr(cls, peak: psc.PyPeak):
        instance = cls.__new__(cls)
        instance.__peak_ptr = peak
        return instance

    @property
    def rt(self) -> int:
        return self.__peak_ptr.rt

    @property
    def spectral_angle(self) -> float:
        return self.__peak_ptr.spectral_angle

    @property
    def score(self) -> float:
        return self.__peak_ptr.score

    @property
    def q_value(self) -> float:
        return self.__peak_ptr.q_value

    def __repr__(self):
        return f"Peak(rt: {self.rt}, spectral_angle: {self.spectral_angle}, score: {self.score}, q_value: {self.q_value})"

    def get_py_ptr(self):
        return self.__peak_ptr


class PeakScoringStrategy:
    """PeakScoringStrategy class

    Args:
        strategy (str): The peak scoring strategy, allowed values are: retention_time, spectral_angle, intensity, hybrid
    """
    def __init__(self, strategy: str = "hybrid"):
        strategies = ["retention_time", "spectral_angle", "intensity", "hybrid"]
        if strategy in strategies:
            self.__peak_scoring_strategy_ptr = psc.PyPeakScoringStrategy(strategy)
        else:
            raise ValueError(f"Invalid peak scoring strategy, allowed values are: {strategies}")

    @classmethod
    def from_py_ptr(cls, peak_scoring_strategy: psc.PyPeakScoringStrategy):
        instance = cls.__new__(cls)
        instance.__peak_scoring_strategy_ptr = peak_scoring_strategy
        return instance

    def get_py_ptr(self):
        return self.__peak_scoring_strategy_ptr

    @property
    def strategy(self):
        return self.__peak_scoring_strategy_ptr.strategy

    def __repr__(self):
        return f"PeakScoringStrategy({self.__peak_scoring_strategy_ptr.strategy})"


class IntegrationStrategy:
    def __init__(self, strategy: str = "sum"):
        strategies = ['apex', 'sum']
        if strategy in strategies:
            self.__integration_strategy_ptr = psc.PyIntegrationStrategy(strategy)
        else:
            raise ValueError(f"Invalid integration strategy, allowed values are: {strategies}")

    @classmethod
    def from_py_integration_strategy(cls, integration_strategy: psc.PyIntegrationStrategy):
        instance = cls.__new__(cls)
        instance.__integration_strategy_ptr = integration_strategy
        return instance

    @property
    def strategy(self):
        return self.__integration_strategy_ptr.strategy

    def __repr__(self):
        return f"IntegrationStrategy({self.__integration_strategy_ptr.strategy})"

    def get_py_ptr(self):
        return self.__integration_strategy_ptr


class PrecursorId:
    def __init__(self, peptide_id: PeptideIx):
        self.__precursor_id_ptr = psc.PyPrecursorId(peptide_id.get_py_ptr())

    @classmethod
    def from_py_precursor_id(cls, precursor_id: psc.PyPrecursorId):
        instance = cls.__new__(cls)
        instance.__precursor_id_ptr = precursor_id
        return instance

    @classmethod
    def from_charged(cls, peptide_id: PeptideIx, charge: int) -> 'PrecursorId':
        return cls.from_py_precursor_id(psc.PyPrecursorId.from_charged(peptide_id.get_py_ptr(), charge))

    def get_ptr(self):
        return self.__precursor_id_ptr

    @property
    def peptide_id(self):
        return PeptideIx.from_py_peptide_ix(self.__precursor_id_ptr.peptide_id)

    @property
    def charge(self):
        return self.__precursor_id_ptr.charge

    def __repr__(self):
        return f"PrecursorId(peptide_id: {self.peptide_id}, charge: {self.charge})"


class LfqSettings:
    def __init__(self,
                 peak_scoring_strategy: PeakScoringStrategy = PeakScoringStrategy(),
                 integration_strategy: IntegrationStrategy = IntegrationStrategy(),
                 spectral_angle: float = 0.7,
                 ppm_tolerance: float = 5.0,
                 combine_charge_states: bool = True
                 ):
        self.__lfq_settings_ptr = psc.PyLfqSettings(peak_scoring_strategy.get_py_ptr(),
                                                    integration_strategy.get_py_ptr(),
                                                    spectral_angle, ppm_tolerance, combine_charge_states)

    @classmethod
    def from_py_lfq_settings(cls, lfq_settings: psc.PyLfqSettings):
        instance = cls.__new__(cls)
        instance.__lfq_settings_ptr = lfq_settings
        return instance

    def get_py_ptr(self):
        return self.__lfq_settings_ptr

    @property
    def peak_scoring_strategy(self) -> PeakScoringStrategy:
        return PeakScoringStrategy.from_py_ptr(self.__lfq_settings_ptr.peak_scoring_strategy)

    @property
    def integration_strategy(self) -> IntegrationStrategy:
        return IntegrationStrategy.from_py_integration_strategy(self.__lfq_settings_ptr.integration_strategy)

    @property
    def spectral_angle(self) -> float:
        return self.__lfq_settings_ptr.spectral_angle

    @property
    def ppm_tolerance(self) -> float:
        return self.__lfq_settings_ptr.ppm_tolerance

    @property
    def combine_charge_states(self) -> bool:
        return self.__lfq_settings_ptr.combine_charge_states

    def __repr__(self):
        return (f"LfqSettings(peak_scoring_strategy: {self.peak_scoring_strategy}, "
                f"integration_strategy: {self.integration_strategy}, "
                f"spectral_angle: {self.spectral_angle}, ppm_tolerance: {self.ppm_tolerance}, "
                f"combine_charge_states: {self.combine_charge_states})")


class PrecursorRange:
    def __init__(self, rt: float, mass_lo: float, mass_hi: float, charge: int,
                 isotope: int, peptide: PeptideIx, file_id: int, decoy: bool):
        self.__precursor_range_ptr = psc.PyPrecursorRange(rt, mass_lo, mass_hi, charge, isotope,
                                                          peptide.get_py_ptr(), file_id, decoy)

    @classmethod
    def from_py_precursor_range(cls, precursor_range: psc.PyPrecursorRange):
        instance = cls.__new__(cls)
        instance.__precursor_range_ptr = precursor_range
        return instance

    def get_py_ptr(self):
        return self.__precursor_range_ptr

    @property
    def rt(self) -> float:
        return self.__precursor_range_ptr.rt

    @property
    def mass_lo(self) -> float:
        return self.__precursor_range_ptr.mass_lo

    @property
    def mass_hi(self) -> float:
        return self.__precursor_range_ptr.mass_hi

    @property
    def charge(self) -> int:
        return self.__precursor_range_ptr.charge

    @property
    def isotope(self) -> int:
        return self.__precursor_range_ptr.isotope

    @property
    def peptide(self) -> PeptideIx:
        return PeptideIx.from_py_peptide_ix(self.__precursor_range_ptr.peptide)

    @property
    def file_id(self) -> int:
        return self.__precursor_range_ptr.file_id

    @property
    def decoy(self) -> bool:
        return self.__precursor_range_ptr.decoy

    def __repr__(self):
        return (f"PrecursorRange(rt: {self.rt}, mass_lo: {self.mass_lo}, mass_hi: {self.mass_hi}, "
                f"charge: {self.charge}, isotope: {self.isotope}, peptide: {self.peptide}, "
                f"file_id: {self.file_id}, decoy: {self.decoy})")


class FeatureMap:
    def __init__(self, ranges: List[PrecursorRange], min_rts: List[float], bin_size: int, settings: LfqSettings):
        self.__feature_map_ptr = psc.PyFeatureMap([
            r.get_py_ptr() for r in ranges
        ], min_rts, bin_size, settings.get_py_ptr())

    @classmethod
    def from_py_feature_map(cls, feature_map: psc.PyFeatureMap):
        instance = cls.__new__(cls)
        instance.__feature_map_ptr = feature_map
        return instance

    def get_py_ptr(self):
        return self.__feature_map_ptr

    @property
    def _ranges(self) -> List[PrecursorRange]:
        return [PrecursorRange.from_py_precursor_range(r) for r in self.__feature_map_ptr.ranges]

    @property
    def min_rts(self) -> List[float]:
        return self.__feature_map_ptr.min_rts

    @property
    def bin_size(self) -> int:
        return self.__feature_map_ptr.bin_size

    @property
    def settings(self) -> LfqSettings:
        return LfqSettings.from_py_lfq_settings(self.__feature_map_ptr.settings)

    def get_num_ranges(self) -> int:
        return self.__feature_map_ptr.get_num_ranges()

    def quantify(
            self,
            indexed_db: 'IndexedDatabase',
            ms1: List['ProcessedSpectrum'],
            alignments: List['Alignment'],
    ) -> Dict:
        """Quantify the feature map.

        Args:
            indexed_db: Indexed database
            ms1: List of processed MS-1 (Precursor) spectra
            alignments: List of alignments

        Returns:
            Dict: Dictionary of quantified features
        """
        ret_tmp = self.__feature_map_ptr.quantify(
            indexed_db.get_py_ptr(), [m.get_py_ptr() for m in ms1],
            [a.get_py_ptr() for a in alignments]
        )

        ret_dict = {}

        for key, value in ret_tmp.items():
            key = (PrecursorId.from_py_precursor_id(key[0]), key[1])
            value = (Peak.from_py_ptr(value[0]), value[1])
            ret_dict[key] = value

        return ret_dict

    def quantify_with_mobility(
            self,
            indexed_db: 'IndexedDatabase',
            ms1: List['ProcessedIMSpectrum'],
            alignments: List['Alignment'],
    ) -> Dict:
        """Quantify the feature map using MS1 spectra with ion mobility.

        Args:
            indexed_db: Indexed database
            ms1: List of processed MS-1 spectra with mobility
            alignments: List of alignments

        Returns:
            Dict: Dictionary of quantified features
        """
        ret_tmp = self.__feature_map_ptr.quantify_with_mobility(
            indexed_db.get_py_ptr(), [m.get_py_ptr() for m in ms1],
            [a.get_py_ptr() for a in alignments]
        )

        ret_dict = {}

        for key, value in ret_tmp.items():
            key = (PrecursorId.from_py_precursor_id(key[0]), key[1])
            value = (Peak.from_py_ptr(value[0]), value[1])
            ret_dict[key] = value

        return ret_dict

    def __repr__(self):
        return f"FeatureMap(num_ranges: {self.get_num_ranges()}, bin_size: {self.bin_size}, settings: {self.settings})"


class Query:
    def __init__(self, ranges: List[PrecursorRange], page_lo: int,
                 page_hi: int, bin_size: int, min_rt: float, max_rt: float):
        self.__query_ptr = psc.PyQuery([
            r.get_py_ptr() for r in ranges
        ], page_lo, page_hi, bin_size, min_rt, max_rt)

    @classmethod
    def from_py_query(cls, query: psc.PyQuery):
        instance = cls.__new__(cls)
        instance.__query_ptr = query
        return instance

    def get_py_ptr(self):
        return self.__query_ptr

    @property
    def _ranges(self) -> List[PrecursorRange]:
        return [PrecursorRange.from_py_precursor_range(r) for r in self.__query_ptr.ranges]

    @property
    def page_lo(self) -> int:
        return self.__query_ptr.page_lo

    @property
    def page_hi(self) -> int:
        return self.__query_ptr.page_hi

    @property
    def bin_size(self) -> int:
        return self.__query_ptr.bin_size

    @property
    def min_rt(self) -> float:
        return self.__query_ptr.min_rt

    @property
    def max_rt(self) -> float:
        return self.__query_ptr.max_rt

    def get_num_ranges(self) -> int:
        return self.__query_ptr.get_num_ranges()

    def __repr__(self):
        return f"Query(num_ranges: {self.get_num_ranges()}, page_lo: {self.page_lo}, page_hi: {self.page_hi}, " \
                f"bin_size: {self.bin_size}, min_rt: {self.min_rt}, max_rt: {self.max_rt})"


def build_feature_map(
        features: List[Feature],
        precursor_charge: Tuple[int, int] = (2, 5),
        lfq_settings: LfqSettings = LfqSettings(),
) -> FeatureMap:
    """Build a feature map.

    Args:
        lfq_settings: LFQ settings
        precursor_charge: Precursor charge
        features: List of features

    Returns:
        FeatureMap: Feature map
    """
    py_feature_map = psc.py_build_feature_map(lfq_settings.get_py_ptr(), precursor_charge, [f.get_py_ptr() for f in features])
    return FeatureMap.from_py_feature_map(py_feature_map)

def build_feature_map_psm(
        psms: List[Psm],
        precursor_charge: Tuple[int, int] = (2, 5),
        lfq_settings: LfqSettings = LfqSettings(),
) -> FeatureMap:
    """Build a feature map from PSMs.

    Args:
        lfq_settings: LFQ settings
        precursor_charge: Precursor charge
        psms: List of PSMs

    Returns:
        FeatureMap: Feature map
    """
    py_feature_map = psc.py_build_feature_map_psm(lfq_settings.get_py_ptr(), precursor_charge, [p.sage_feature.get_py_ptr() for p in psms])
    return FeatureMap.from_py_feature_map(py_feature_map)
