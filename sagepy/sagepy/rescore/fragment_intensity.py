from typing import List, Tuple, Dict

import numpy as np
import sagepy_connector
from numpy.typing import NDArray

from sagepy.core import PeptideSpectrumMatch

psc = sagepy_connector.py_intensity

class FragmentIntensity:
    def __init__(self, intensities_observed: List[float], mz_observed: List[float], mz_calculated: List[float],
                 charges: List[int], ordinals: List[int], ion_types: List[bool], prosit_intensity_predicted: List[float]):
        self.__py_ptr = psc.PyFragmentIntensityPrediction(
            intensities_observed, mz_observed, mz_calculated, charges, ordinals, ion_types, prosit_intensity_predicted)

    def get_py_ptr(self):
        return self.__py_ptr

    def from_py_ptr(cls, py_ptr):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    @property
    def intensities_observed(self) -> List[float]:
        return self.__py_ptr.intensities_observed

    @intensities_observed.setter
    def intensities_observed(self, intensities_observed: List[float]):
        self.__py_ptr.intensities_observed = intensities_observed

    @property
    def mz_observed(self) -> List[float]:
        return self.__py_ptr.mz_observed

    @mz_observed.setter
    def mz_observed(self, mz_observed: List[float]):
        self.__py_ptr.mz_observed = mz_observed

    @property
    def mz_calculated(self) -> List[float]:
        return self.__py_ptr.mz_calculated

    @mz_calculated.setter
    def mz_calculated(self, mz_calculated: List[float]):
        self.__py_ptr.mz_calculated = mz_calculated

    @property
    def charges(self) -> List[int]:
        return self.__py_ptr.charges

    @charges.setter
    def charges(self, charges: List[int]):
        self.__py_ptr.charges = charges

    @property
    def ordinals(self) -> List[int]:
        return self.__py_ptr.ordinals

    @ordinals.setter
    def ordinals(self, ordinals: List[int]):
        self.__py_ptr.ordinals = ordinals

    @property
    def ion_types(self) -> List[bool]:
        return self.__py_ptr.ion_types

    @ion_types.setter
    def ion_types(self, ion_types: List[bool]):
        self.__py_ptr.ion_types = ion_types

    @property
    def prosit_intensity_predicted(self) -> List[float]:
        return self.__py_ptr.prosit_intensity_predicted

    @intensities_observed.setter
    def intensities_observed(self, intensities_observed: List[float]):
        self.__py_ptr.intensities_observed = intensities_observed

    def cosine_similarity(self, epsilon: float = 1e-7, reduce_matched: bool = False) -> float:
        return self.__py_ptr.cosine_similarity(epsilon, reduce_matched)

    def spectral_angle_similarity(self, epsilon: float = 1e-7, reduce_matched: bool = False) -> float:
        return self.__py_ptr.spectral_angle_similarity(epsilon, reduce_matched)

    def pearson_correlation(self, epsilon: float = 1e-7, reduce_matched: bool = False) -> float:
        return self.__py_ptr.pearson_correlation(epsilon, reduce_matched)

    def spearman_correlation(self, epsilon: float = 1e-7, reduce_matched: bool = False) -> float:
        return self.__py_ptr.spearman_correlation(epsilon, reduce_matched)

    def spectral_entropy_similarity(self, epsilon: float = 1e-7, reduce_matched: bool = False) -> float:
        return self.__py_ptr.spectral_entropy_similarity(epsilon, reduce_matched)

    def __repr__(self):
        return (f"FragmentIntensity(intensities_observed={self.intensities_observed}, "
                f"mz_observed={self.mz_observed}, mz_calculated={self.mz_calculated}, "
                f"charges={self.charges}, ordinals={self.ordinals}, ion_types={self.ion_types}, "
                f"prosit_intensity_predicted={self.prosit_intensity_predicted})")

    def observed_intensity_map(self) -> Dict[Tuple[int, int, int], float]:
        return self.__py_ptr.observed_intensity_map()

    def predicted_intensity_map(self) -> Dict[Tuple[int, int, int], float]:
        return self.__py_ptr.predicted_intensity_map()

def peptide_spectrum_match_list_to_intensity_feature_matrix(
        psm_list: List[PeptideSpectrumMatch],
        epsilon: float = 1e-7,
        reduce_matched: bool = False,
        num_threads: int = 16,
) -> NDArray:
    features = psc.peptide_spectrum_match_list_to_intensity_feature_matrix_parallel(
        [p.get_py_ptr() for p in psm_list], epsilon, reduce_matched, num_threads
    )
    return np.array(features)
