from typing import List, Tuple

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

