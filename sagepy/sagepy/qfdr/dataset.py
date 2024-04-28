from typing import List, Union, Tuple

import pandas as pd
import sagepy_connector

psc = sagepy_connector.py_qfdr


class TDCMethod:
    def __init__(self, method: str):
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
                 spec_idx: str,
                 peptide_idx: int,
                 proteins: List[str],
                 decoy: bool,
                 hyper_score: float,
                 rank: int,
                 mono_mass_observed: Union[None, float],
                 sequence: Union[None, str],
                 charge: Union[None, int],
                 retention_time_observed: Union[None, float],
                 retention_time_predicted: Union[None, float],
                 inverse_mobility_observed: Union[None, float],
                 inverse_mobility_predicted: Union[None, float],
                 intensity_ms1: Union[None, float],
                 intensity_ms2: Union[None, float],
                 q_value: Union[None, float],
                 re_score: Union[None, float]):
        self.__py_ptr = psc.PyPeptideSpectrumMatch(
            spec_idx, peptide_idx, proteins, decoy, hyper_score, rank, mono_mass_observed, sequence, charge,
            retention_time_observed, retention_time_predicted, inverse_mobility_observed, inverse_mobility_predicted,
            intensity_ms1, intensity_ms2, q_value, re_score
        )

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

    @hyper_score.setter
    def hyper_score(self, value):
        self.__py_ptr.hyper_score = value

    @property
    def rank(self):
        return self.__py_ptr.rank

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

    @retention_time_predicted.setter
    def retention_time_predicted(self, value):
        self.__py_ptr.retention_time_predicted = value

    @property
    def inverse_mobility_observed(self):
        return self.__py_ptr.inverse_mobility_observed

    @property
    def inverse_mobility_predicted(self):
        return self.__py_ptr.inverse_mobility_predicted

    @inverse_mobility_predicted.setter
    def inverse_mobility_predicted(self, value):
        self.__py_ptr.inverse_mobility_predicted = value

    @property
    def intensity_ms1(self):
        return self.__py_ptr.intensity_ms1

    @property
    def intensity_ms2(self):
        return self.__py_ptr.intensity_ms2

    @property
    def q_value(self):
        return self.__py_ptr.q_value

    @property
    def re_score(self):
        return self.__py_ptr.re_score

    @re_score.setter
    def re_score(self, value):
        self.__py_ptr.re_score = value

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyPeptideSpectrumMatch):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyPeptideSpectrumMatch:
        return self.__py_ptr

    def __repr__(self):
        return f"PeptideSpectrumMatch({self.spec_idx}, {self.peptide_idx}, {self.proteins}, {self.decoy}, " \
               f"{self.hyper_score}, {self.rank} {self.charge}, {self.peptide_sequence}, {self.mono_mass_observed}, " \
               f"{self.mono_mass_calculated}, {self.retention_time_observed}, {self.retention_time_predicted}, " \
               f"{self.inverse_mobility_observed}, {self.inverse_mobility_predicted}, {self.intensity_ms1}, " \
               f"{self.intensity_ms2}, {self.q_value}, {self.re_score})"


def target_decoy_competition(method: str, spectra_idx: List[str], match_idx: List[int], decoy: List[bool],
                             scores: List[float]) -> Tuple[List[str], List[int], List[bool], List[float], List[float]]:
    tdc_method = TDCMethod(method)
    spec_idx, match_idx, decoy, scores, q_values = psc.target_decoy_competition(
        tdc_method.get_py_ptr(), spectra_idx,
        match_idx, decoy, scores)
    return spec_idx, match_idx, decoy, scores, q_values


def target_decoy_competition_pandas(method: str, df: pd.DataFrame) -> pd.DataFrame:
    assert 'spec_idx' in df.columns, "spec_idx column not found"
    assert 'match_idx' in df.columns, "match_idx column not found"
    assert 'decoy' in df.columns, "decoy column not found"
    assert 'score' in df.columns, "score column not found"

    spec_idx, match_idx, target, scores = (df['spec_idx'].tolist(),
                                           df['match_idx'].tolist(), df['decoy'].tolist(), df['score'].tolist())

    spec_idx, match_idx, target, scores, q_values = target_decoy_competition(method, spec_idx,
                                                                             match_idx, target, scores)

    return pd.DataFrame({
        'spec_idx': spec_idx,
        'match_idx': match_idx,
        'decoy': target,
        'score': scores,
        'q_value': q_values
    }).sort_values('q_value')

