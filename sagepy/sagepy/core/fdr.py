from typing import Optional
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
