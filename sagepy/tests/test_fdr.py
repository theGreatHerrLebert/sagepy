"""Tests for sagepy.core.fdr."""

from sagepy.core.database import PeptideIx
import sagepy.core.fdr as fdr


class _DummyPtrWrapper:
    def __init__(self, value):
        self._value = value

    def get_py_ptr(self):
        return self._value


def test_competition_peptide_ix_properties():
    competition = fdr.CompetitionPeptideIx(1.5, 0.5, PeptideIx(7), PeptideIx(9))

    assert competition.forward == 1.5
    assert competition.reverse == 0.5
    assert competition.forward_ix.idx == 7
    assert competition.reverse_ix.idx == 9
    assert "CompetitionPeptideIx" in repr(competition)


def test_sage_fdr_forwards_feature_ptrs(monkeypatch):
    calls = []

    class _DummyBackend:
        @staticmethod
        def py_sage_fdr(features, indexed_db, use_hyper_score):
            calls.append((features, indexed_db, use_hyper_score))

    monkeypatch.setattr(fdr, "psc", _DummyBackend)

    features = [_DummyPtrWrapper("f1"), _DummyPtrWrapper("f2")]
    indexed_db = _DummyPtrWrapper("db")

    fdr.sage_fdr(features, indexed_db, use_hyper_score=False)

    assert calls == [(["f1", "f2"], "db", False)]


def test_sage_fdr_psm_flattens_dict(monkeypatch):
    calls = []

    class _DummyBackend:
        @staticmethod
        def py_sage_fdr_psm(features, indexed_db, use_hyper_score):
            calls.append((features, indexed_db, use_hyper_score))

    monkeypatch.setattr(fdr, "psc", _DummyBackend)

    psm_collection = {
        "spec1": [_DummyPtrWrapper("p1"), _DummyPtrWrapper("p2")],
        "spec2": [_DummyPtrWrapper("p3")],
    }
    indexed_db = _DummyPtrWrapper("db")

    fdr.sage_fdr_psm(psm_collection, indexed_db)

    assert calls == [(["p1", "p2", "p3"], "db", True)]
