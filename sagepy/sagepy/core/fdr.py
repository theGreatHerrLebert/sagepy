from typing import Optional, List, Union, Dict
from sagepy.core.scoring import Psm
from sagepy.core import IndexedDatabase, Feature
from sagepy.core.database import PeptideIx
from sagepy.core.mass import Tolerance
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


def assign_initial_spectrum_q_psms(psm_collection: Union[List[Psm], Dict[str, List[Psm]]]) -> None:
    """Assign the initial spectrum q-values sage uses before RT/IM predictor fitting.

    This is the pre-LDA pass that trains on the best-scoring PSMs before the
    retention-time and ion-mobility regressors are fit.
    """
    if isinstance(psm_collection, dict):
        f_collection = [p for psms in psm_collection.values() for p in psms]
    else:
        f_collection = psm_collection

    psc.py_assign_initial_spectrum_q_psms(
        [p.get_py_ptr() for p in f_collection],
    )


def lda_rescore_psms(
    psm_collection: Union[List[Psm], Dict[str, List[Psm]]],
    precursor_tolerance: Tolerance,
) -> bool:
    """Fit sage's linear discriminant model on the PSMs and write the resulting
    discriminant_score and posterior_error onto each PSM (in-place).

    This mirrors the step `sage` CLI runs when `predict_rt: true`. After this,
    call ``sage_fdr_psm(..., use_hyper_score=False)`` to get q-values from the
    discriminant. RT and IM features are pulled from the PSM directly, so call
    ``predict_sage_rt`` and ``predict_sage_im`` first if you want the model to
    weight those features.

    Args:
        psm_collection: PSMs (list or {spec_id: [Psm]} dict).
        precursor_tolerance: precursor tolerance from the search config; the LDA
            uses it to choose KDE binning over delta_mass.

    Returns:
        True if the LDA fit succeeded, False if sage fell back to the heuristic
        ``(-poisson).ln_1p() + longest_y_pct/3.0`` discriminant.
    """
    if isinstance(psm_collection, dict):
        f_collection = [p for psms in psm_collection.values() for p in psms]
    else:
        f_collection = psm_collection

    return psc.py_lda_score_psms(
        [p.get_py_ptr() for p in f_collection],
        precursor_tolerance.get_py_ptr(),
    )
