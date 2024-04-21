from typing import List, Union, Tuple

import sagepy_connector
psc = sagepy_connector.py_qfdr


class PeptideSpectrumMatch:
    def __init__(self, spec_id: str, peptide_id: int, proteins: List[str], decoy: bool, score: float, intensity: float, features: Union[None, List[Tuple[str, float]]] = None):
        self.__py_ptr = psc.PyPeptideSpectrumMatch(spec_id, peptide_id, proteins, decoy, score, intensity, features)

    @property
    def spec_id(self) -> str:
        return self.__py_ptr.spec_id

    @property
    def peptide_id(self) -> int:
        return self.__py_ptr.peptide_id

    @property
    def proteins(self) -> List[str]:
        return self.__py_ptr.proteins

    @property
    def decoy(self) -> bool:
        return self.__py_ptr.decoy

    @property
    def score(self) -> float:
        return self.__py_ptr.score

    @property
    def intensity(self) -> float:
        return self.__py_ptr.intensity

    @property
    def features(self) -> Union[None, List[Tuple[str, float]]]:
        return self.__py_ptr.features

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyPeptideSpectrumMatch):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyPeptideSpectrumMatch:
        return self.__py_ptr

    def __repr__(self):
        return (f"PeptideSpectrumMatch(spec_id: {self.spec_id}, peptide_id: {self.peptide_id}, "
                f"proteins: {self.proteins}, decoy: {self.decoy}, intensity: {self.intensity},"
                f"score: {self.score}, features: {self.features})")


class PsmDataset:
    def __init__(self, spec_ids: List[str], matches: List[List[PeptideSpectrumMatch]]):
        self.__py_ptr = psc.PyPsmDataset(spec_ids, matches)

    @property
    def size(self):
        return self.__py_ptr.size

    @property
    def keys(self):
        return self.__py_ptr.keys

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyPsmDataset):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyPsmDataset:
        return self.__py_ptr

    def get_spec_psms(self, spec_id: str) -> List[PeptideSpectrumMatch]:
        return [PeptideSpectrumMatch.from_py_ptr(psm) for psm in self.__py_ptr.get_spec_psms(spec_id)]

    def __repr__(self):
        return f"PsmDataset(scored spectra: {self.size})"
