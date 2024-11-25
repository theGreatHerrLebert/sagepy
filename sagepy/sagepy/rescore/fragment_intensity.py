from typing import List, Tuple, Dict
import sagepy_connector
from sagepy.core import Fragments

psc = sagepy_connector.py_intensity

class FragmentIntensity:
    def __init__(self,
                 fragments_observed: Fragments,
                 prosit_intensity_predicted: List[float]):
        self.__py_ptr = psc.PyFragmentIntensityPrediction(
            fragments_observed, prosit_intensity_predicted
        )

    def get_py_ptr(self):
        return self.__py_ptr

    @classmethod
    def from_py_ptr(cls, py_ptr):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    @property
    def prosit_intensity_predicted(self) -> List[float]:
        return self.__py_ptr.prosit_intensity_predicted

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

    def observed_intensity_map(self) -> Dict[Tuple[int, int, int], float]:
        return self.__py_ptr.observed_intensity_map()

    def predicted_intensity_map(self) -> Dict[Tuple[int, int, int], float]:
        return self.__py_ptr.predicted_intensity_map()

    def __repr__(self):
        return f"FragmentIntensity(prosit_intensity_predicted={self.prosit_intensity_predicted})"