from typing import Optional, List, Union, Dict
from sagepy.core.scoring import Psm
from sagepy.core import IndexedDatabase, Feature
from sagepy.core.database import PeptideIx
import sagepy_connector
psc = sagepy_connector.py_fdr


class CompetitionPeptideIx:
    def __init__(self, forward: float, reverse: float,
                 forward_ix: Optional[PeptideIx] = None, reverse_ix: Optional[PeptideIx] = None):

        self.__ptr = psc.PyCompetitionPeptideIx(forward, reverse, forward_ix.get_py_ptr() if forward_ix else None,
                                                reverse_ix.get_py_ptr() if reverse_ix else None)

    @classmethod
    def from_py_competition_peptide_ix(cls, competition_peptide_ix: psc.PyCompetitionPeptideIx):
        instance = cls.__new__(cls)
        instance.__ptr = competition_peptide_ix
        return instance

    def get_py_ptr(self):
        return self.__ptr

    @property
    def forward(self) -> float:
        return self.__ptr.forward

    @property
    def reverse(self) -> float:
        return self.__ptr.reverse

    @property
    def forward_ix(self) -> Optional[PeptideIx]:
        maybe_forward = self.__ptr.forward_ix
        if maybe_forward is None:
            return None
        return PeptideIx.from_py_peptide_ix(maybe_forward)

    @property
    def reverse_ix(self) -> Optional[PeptideIx]:
        maybe_reverse = self.__ptr.reverse_ix
        if maybe_reverse is None:
            return None
        return PeptideIx.from_py_peptide_ix(maybe_reverse)

    def __repr__(self):
        return (f"CompetitionPeptideIx(forward={self.forward}, reverse={self.reverse}, "
                f"forward_ix={self.forward_ix}, reverse_ix={self.reverse_ix})")

def sage_fdr(feature_collection: List[Feature], indexed_db: IndexedDatabase, use_hyper_score: bool = True):
    """ Perform SAGE FDR on all levels (spectrum, peptide, protein), calculates q-values and PEPs for a given feature collection.
    Args:
        feature_collection: a list of features
        indexed_db: an indexed database
        use_hyper_score: whether to use hyper score or discriminant score for q-value calculation
    """
    psc.py_sage_fdr(
        [feature.get_py_ptr() for feature in feature_collection],
        indexed_db.get_py_ptr(),
        use_hyper_score
    )

def sage_fdr_psm(psm_collection: Union[List[Psm], Dict[str, List[Psm]]], indexed_db: IndexedDatabase, use_hyper_score: bool = True):
    """ Perform SAGE FDR on all levels (spectrum, peptide, protein), calculates q-values and PEPs for a given feature collection.
    Args:
        psm_collection: a list of features
        indexed_db: an indexed database
        use_hyper_score: whether to use hyper score or discriminant score for q-value calculation
    """

    f_collection = []

    if isinstance(psm_collection, dict):
        for _, values in psm_collection.items():
            f_collection.extend(values)

    else:
        f_collection = psm_collection

    psc.py_sage_fdr_psm(
        [feature.get_py_ptr() for feature in f_collection],
        indexed_db.get_py_ptr(),
        use_hyper_score
    )
