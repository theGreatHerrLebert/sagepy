"""Tests for sagepy.core.peptide module."""
import pytest
from sagepy.core.peptide import Peptide
from sagepy.core.enzyme import Position, EnzymeParameters


class TestPeptide:
    """Tests for the Peptide class."""

    def test_creation(self):
        """Test creating a peptide."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST_PROTEIN']
        )
        assert peptide.sequence == 'PEPTIDE'
        assert peptide.decoy is False
        assert len(peptide.modifications) == 7

    def test_with_modifications(self):
        """Test creating a peptide with modifications."""
        mods = [0.0, 0.0, 0.0, 15.9949, 0.0, 0.0, 0.0]  # Oxidation on 4th residue
        peptide = Peptide(
            decoy=False,
            sequence='PEPTMDE',
            modifications=mods,
            mono_isotopic=815.35,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST_PROTEIN']
        )
        assert peptide.modifications[3] == pytest.approx(15.9949, rel=1e-3)

    def test_with_terminal_mods(self):
        """Test creating a peptide with terminal modifications."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST_PROTEIN'],
            n_term=42.0106,  # Acetylation
            c_term=None
        )
        assert peptide.n_term == pytest.approx(42.0106, rel=1e-3)
        assert peptide.c_term is None

    def test_multiple_proteins(self):
        """Test peptide with multiple proteins."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['PROT1', 'PROT2', 'PROT3']
        )
        assert len(peptide.proteins) == 3
        assert 'PROT1' in peptide.proteins

    def test_missed_cleavages(self):
        """Test missed cleavages property."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDEKR',
            modifications=[0.0] * 9,
            mono_isotopic=1000.0,
            missed_cleavages=2,
            position=Position('internal'),
            proteins=['TEST']
        )
        assert peptide.missed_cleavages == 2

    def test_decoy_flag(self):
        """Test decoy flag."""
        target = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST']
        )
        decoy = Peptide(
            decoy=True,
            sequence='EDITPEP',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST']
        )
        assert target.decoy is False
        assert decoy.decoy is True

    def test_semi_enzymatic(self):
        """Test semi-enzymatic flag."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST'],
            semi_enzymatic=True
        )
        assert peptide.semi_enzymatic is True

    def test_reverse(self):
        """Test reversing a peptide."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST']
        )
        reversed_pep = peptide.reverse(keep_ends=False)
        assert reversed_pep.sequence == 'EDITPEP'

    def test_reverse_keep_ends(self):
        """Test reversing a peptide keeping terminal residues."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST']
        )
        reversed_pep = peptide.reverse(keep_ends=True)
        # First and last residue should be the same
        assert reversed_pep.sequence[0] == peptide.sequence[0]
        assert reversed_pep.sequence[-1] == peptide.sequence[-1]

    def test_from_enzyme_parameters_digest(self):
        """Test creating peptide from EnzymeParameters digest method."""
        params = EnzymeParameters(missed_cleavages=0, min_len=3, max_len=20)
        digests = params.digest('PEPTIDEKPEPTIDER', 'TEST_PROTEIN')
        assert len(digests) > 0
        # Get first digest and convert to peptide
        peptide = Peptide.from_digest(digests[0])
        assert 'TEST_PROTEIN' in peptide.proteins

    def test_repr(self):
        """Test string representation."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('internal'),
            proteins=['TEST']
        )
        repr_str = repr(peptide)
        assert 'PEPTIDE' in repr_str
        assert 'Peptide' in repr_str

    def test_position_property(self):
        """Test position property."""
        peptide = Peptide(
            decoy=False,
            sequence='PEPTIDE',
            modifications=[0.0] * 7,
            mono_isotopic=799.36,
            missed_cleavages=0,
            position=Position('n_term'),
            proteins=['TEST']
        )
        assert peptide.position.position == 'Nterm'
