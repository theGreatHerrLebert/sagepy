from typing import List, Tuple, Optional

import pandas as pd
import sagepy_connector

from sagepy.core import Psm

psc = sagepy_connector.py_qfdr


class TDCMethod:
    def __init__(self, method: str):
        self.methods = {"psm", "peptide_psm_only", "peptide_peptide_only", "peptide_psm_peptide", "picked_peptide", "picked_protein"}
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


def target_decoy_competition(
        spectra_idx: List[str],
        match_idx: List[str],
        match_identity_candidates: List[List[str]],
        decoy: List[bool],
        scores: List[float],
        method: str = "peptide_psm_peptide") -> Tuple[List[str], List[str], List[List[str]], List[bool], List[float], List[float]]:
    """ Perform target-decoy competition.

    Args:
        spectra_idx: a list of spectrum indices
        match_idx: a list of match indices
        match_identity_candidates: a list of match identity candidates
        decoy: a list of decoy flags
        scores: a list of scores
        method: the method to use, allowed values are: psm, peptide_psm_only, peptide_peptide_only, peptide_psm_peptide

    Returns:
        a tuple of spectrum indices, match indices, decoy flags, scores, and q-values
    """
    tdc_method = TDCMethod(method)
    spec_idx, match_idx, match_identity_candidates, decoy, scores, q_values = psc.target_decoy_competition(
        tdc_method.get_py_ptr(), spectra_idx,
        match_idx, decoy, scores, match_identity_candidates)
    return spec_idx, match_idx, match_identity_candidates, decoy, scores, q_values


def target_decoy_competition_pandas(
        df: pd.DataFrame,
        method: str = "peptide_psm_peptide",
        score: Optional[str] = None,
) -> pd.DataFrame:
    """ Perform target-decoy competition on a pandas DataFrame.

    Args:
        df: a pandas DataFrame
        method: the method to use, allowed values are: psm, peptide_psm_only, peptide_peptide_only, peptide_psm_peptide
        score: the target column name (optional)

    Returns:
        a pandas DataFrame with q-values
    """

    # Ensure necessary columns are present
    required_columns = ['spec_idx', 'match_idx', 'match_identity_candidates', 'decoy']
    for col in required_columns:
        assert col in df.columns, f"{col} column not found"

    # Ensure score column is present
    score_col = score if score else 'hyperscore'
    assert score_col in df.columns, f"{score_col} column not found"

    target_score = df[score_col]
    spec_idx, match_idx, match_identity_candidates, target, scores = (
        df['spec_idx'].tolist(),
        df['match_idx'].tolist(),
        df['match_identity_candidates'].tolist(),
        df['decoy'].tolist(),
        target_score.tolist()
    )

    spec_idx, match_idx, match_identity_candidates, target, scores, q_values = target_decoy_competition(
        spec_idx,
        match_idx,
        match_identity_candidates,
        target,
        scores,
        method
    )

    # Create df with TDC results
    df_tdc = pd.DataFrame({
        'spec_idx': spec_idx,
        'match_idx': match_idx,
        'match_identity_candidates': match_identity_candidates,
        'decoy': target,
        f'{score_col}': scores,
        'q_value': q_values
    }).sort_values(by=['q_value'], ascending=True)

    return df_tdc

def assign_sage_spectrum_q(psm_list: List[Psm], use_hyper_score: bool = True):
    """ Assign SAGE spectrum q-values to PSMs.
    Args:
        psm_list: a list of PeptideSpectrumMatch objects
        use_hyper_score: whether to use hyper score or discriminant score for q-value calculation
    """
    # Perform SAGE FDR
    psc.assign_spectrum_q([psm.get_py_ptr() for psm in psm_list], use_hyper_score)

def assign_sage_peptide_q(psm_list: List[Psm], use_hyper_score: bool = True):
    """ Assign SAGE peptide q-values to PSMs.
    Args:
        psm_list: a list of PeptideSpectrumMatch objects
        use_hyper_score: whether to use hyper score or discriminant score for q-value calculation
    """
    # Perform SAGE FDR
    psc.assign_peptide_q([psm.get_py_ptr() for psm in psm_list], use_hyper_score)

def assign_sage_protein_q(psm_list: List[Psm], use_hyper_score: bool = True):
    """ Assign SAGE protein q-values to PSMs.
    Args:
        psm_list: a list of PeptideSpectrumMatch objects
        use_hyper_score: whether to use hyper score or discriminant score for q-value calculation
    """
    # Perform SAGE FDR
    psc.assign_protein_q([psm.get_py_ptr() for psm in psm_list], use_hyper_score)