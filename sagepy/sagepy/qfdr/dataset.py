import pandas as pd
from typing import List, Union, Tuple

import sagepy_connector

psc = sagepy_connector.py_qfdr


class TDCMethod:
    def __init__(self, method: int):
        self.methods = {"psm", "peptide_psm_only", "peptide_peptide_only", "peptide_psm_peptide"}
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
    def __init__(self,
                 spec_idx: str, peptide_idx: int, proteins: List[str], decoy: bool, hyper_score: float,
                 mono_mass_observed: Union[None, float], mono_mass_predicted: Union[None, float],
                 charge: Union[None, int], peptide_sequence: Union[None, str],
                 retention_time_observed: Union[None, float], retention_time_predicted: Union[None, float],
                 inverse_mobility_observed: Union[None, float], inverse_mobility_predicted: Union[None, float],
                 intensity_ms1: Union[None, float], intensity_ms2: Union[None, float], q_value: Union[None, float]):

        self.__py_ptr = psc.PyPeptideSpectrumMatch(spec_idx, peptide_idx, proteins, decoy, hyper_score, charge,
                                                   mono_mass_observed, mono_mass_predicted,
                                                   peptide_sequence, retention_time_observed, retention_time_predicted,
                                                   inverse_mobility_observed, inverse_mobility_predicted, intensity_ms1,
                                                   intensity_ms2, q_value)

    @property
    def spec_idx(self):
        return self.__py_ptr.spec_idx

    @property
    def peptide_idx(self):
        return self.__py_ptr.peptide_idx

    @property
    def proteins(self):
        return self.__py_ptr.proteins

    @property
    def decoy(self):
        return self.__py_ptr.decoy

    @property
    def hyper_score(self):
        return self.__py_ptr.hyper_score

    @property
    def charge(self):
        return self.__py_ptr.charge

    @property
    def peptide_sequence(self):
        return self.__py_ptr.peptide_sequence

    @property
    def mono_mz_calculated(self):
        return self.__py_ptr.mono_mz_calculated

    @property
    def mono_mass_observed(self):
        return self.__py_ptr.mono_mass_observed

    @property
    def mono_mass_calculated(self):
        return self.__py_ptr.mono_mass_calculated

    @property
    def retention_time_observed(self):
        return self.__py_ptr.retention_time_observed

    @property
    def retention_time_predicted(self):
        return self.__py_ptr.retention_time_predicted

    @property
    def inverse_mobility_observed(self):
        return self.__py_ptr.inverse_mobility_observed

    @property
    def inverse_mobility_predicted(self):
        return self.__py_ptr.inverse_mobility_predicted

    @property
    def intensity_ms1(self):
        return self.__py_ptr.intensity_ms1

    @property
    def intensity_ms2(self):
        return self.__py_ptr.intensity_ms2

    @property
    def q_value(self):
        return self.__py_ptr.q_value

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyPeptideSpectrumMatch):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyPeptideSpectrumMatch:
        return self.__py_ptr

    def __repr__(self):
        return f"PeptideSpectrumMatch({self.spec_idx}, {self.peptide_idx}, {self.proteins}, {self.decoy}, " \
               f"{self.hyper_score}, {self.charge}, {self.peptide_sequence}, {self.mono_mass_observed}, " \
               f"{self.mono_mass_calculated}, {self.retention_time_observed}, {self.retention_time_predicted}, " \
               f"{self.inverse_mobility_observed}, {self.inverse_mobility_predicted}, {self.intensity_ms1}, " \
               f"{self.intensity_ms2}, {self.q_value})"


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

    def best_target_psm(self, spec_id: str) -> Union[None, PeptideSpectrumMatch]:
        psm = self.__py_ptr.best_target_psm(spec_id)
        if psm is None:
            return None
        return PeptideSpectrumMatch.from_py_ptr(psm)

    def best_decoy_psm(self, spec_id: str) -> Union[None, PeptideSpectrumMatch]:
        psm = self.__py_ptr.best_decoy_psm(spec_id)
        if psm is None:
            return None
        return PeptideSpectrumMatch.from_py_ptr(psm)

    def target_decoy_competition(self, method: TDCMethod) -> List[PeptideSpectrumMatch]:
        return [PeptideSpectrumMatch.from_py_ptr(psm) for psm in self.__py_ptr.tdc(method.get_py_ptr())]

    def target_decoy_competition_pandas(self, method: TDCMethod) -> pd.DataFrame:
        matches = [PeptideSpectrumMatch.from_py_ptr(psm) for psm in self.__py_ptr.tdc(method.get_py_ptr())]

        row_list = []

        for match in matches:
            row_dict = {'spec_idx': match.spec_idx, 'peptide_idx': match.peptide_idx, 'proteins': match.proteins,
                        'decoy': match.decoy, 'hyper_score': match.hyper_score, 'charge': match.charge,
                        'peptide_sequence': match.peptide_sequence, 'retention_time_observed': match.retention_time_observed,
                        'retention_time_predicted': match.retention_time_predicted,
                        'inverse_mobility_observed': match.inverse_mobility_observed,
                        'inverse_mobility_predicted': match.inverse_mobility_predicted,
                        'intensity_ms1': match.intensity_ms1, 'intensity_ms2': match.intensity_ms2, 'q_value': match.q_value}
            row_list.append(row_dict)

        return pd.DataFrame(row_list)

    def __repr__(self):
        return f"PsmDataset(scored spectra: {self.size})"
