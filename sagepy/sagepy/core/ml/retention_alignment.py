from typing import List
import sagepy_connector
from imspy.simulation.annotation import RustWrapperObject

from sagepy.core import Psm, Feature

psc = sagepy_connector.py_retention_alignment


class Alignment(RustWrapperObject):
    def __init__(self, file_id: int, max_rt: float, slope: float, intercept: float):
        self.__py_ptr = psc.PyAlignment(file_id, max_rt, slope, intercept)

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyAlignment):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyAlignment:
        return self.__py_ptr

    @property
    def file_id(self) -> int:
        return self.__py_ptr.file_id

    @property
    def max_rt(self) -> float:
        return self.__py_ptr.max_rt

    @property
    def slope(self) -> float:
        return self.__py_ptr.slope

    @property
    def intercept(self) -> float:
        return self.__py_ptr.intercept

    def __repr__(self):
        return f"Alignment(file_id={self.file_id}, max_rt={self.max_rt}, slope={self.slope}, intercept={self.intercept})"


def global_alignment(features: List[Feature], n_files: int) -> List[Alignment]:
    """ Perform global alignment.
    Args:
        features: A list of features
        n_files: Number of files
    Returns:
        List[Alignment]: List of Alignment objects
    """

    py_alignments = psc.py_global_alignment([f.get_py_ptr() for f in features], n_files)
    return [Alignment.from_py_ptr(py_alignment) for py_alignment in py_alignments]

def global_alignment_psm(psms: List[Psm]) -> List[Alignment]:
    """ Perform global alignment on PSMs.
    Args:
        psms: A list of PSMs
    Returns:
        List[Alignment]: List of Alignment objects
    """

    n_files = len(set([p.sage_feature.file_id for p in psms]))

    py_alignments = psc.py_global_alignment_psm([p.get_py_ptr() for p in psms], n_files)
    return [Alignment.from_py_ptr(py_alignment) for py_alignment in py_alignments]
