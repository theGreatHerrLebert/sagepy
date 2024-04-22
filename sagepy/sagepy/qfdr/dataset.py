from typing import List, Union, Tuple

import sagepy_connector
psc = sagepy_connector.py_qfdr


class TDCMethod:
    def __init__(self, method: int):
        self.methods = {"psm", "peptide_psm_only", "peptide_peptide_only", "peptide_psm_and_peptide"}
        assert method in self.methods, f"Invalid method: {method}, allowed values are: {self.methods}"
        self.__py_ptr = psc.PyTDCMethod(method)

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyTDCMethod):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyTDCMethod:
        return self.__py_ptr

    def __repr__(self):
        return f"TDCMethod({self.__py_ptr.to_str()})"


class PeptideSpectrumMatch:
    def __init__(self, spec_id: str, peptide_id: int,
                 proteins: List[str], decoy: bool, score: float,
                 intensity_ms1: Union[None, float],
                 intensity_ms2: Union[None, float] = None, features: Union[None, List[Tuple[str, float]]] = None,
                 q_value: Union[None, float] = None, confidence: Union[None, float] = None,
                 ):
        self.__py_ptr = psc.PyPeptideSpectrumMatch(spec_id, peptide_id, proteins,
                                                   decoy, score, intensity_ms1, intensity_ms2,
                                                   features, q_value, confidence)

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
    def intensity_ms1(self) -> Union[None, float]:
        return self.__py_ptr.intensity_ms1

    @property
    def intensity_ms2(self) -> Union[None, float]:
        return self.__py_ptr.intensity_ms1

    @property
    def features(self) -> Union[None, List[Tuple[str, float]]]:
        return self.__py_ptr.features

    @property
    def q_value(self) -> Union[None, float]:
        return self.__py_ptr.q_value

    @property
    def confidence(self) -> Union[None, float]:
        return self.__py_ptr.confidence

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyPeptideSpectrumMatch):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyPeptideSpectrumMatch:
        return self.__py_ptr

    def __repr__(self):
        return (f"PeptideSpectrumMatch(spec_id: {self.spec_id}, peptide_id: {self.peptide_id}, "
                f"proteins: {self.proteins}, decoy: {self.decoy}, "
                f"intensity_ms1: {self.intensity_ms1}, intensity_ms2: {self.intensity_ms2}, "
                f"score: {self.score}, features: {self.features}, "
                f"q_value: {self.q_value}, confidence: {self.confidence})")


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

    def target_decoy_competition(self, method: TDCMethod) -> List[PeptideSpectrumMatch]:
        return [PeptideSpectrumMatch.from_py_ptr(psm) for psm in self.__py_ptr.tdc(method.get_py_ptr())]

    def __repr__(self):
        return f"PsmDataset(scored spectra: {self.size})"
