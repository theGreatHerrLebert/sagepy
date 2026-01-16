"""Tests for the IndexedDatabase save/load functionality."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from sagepy.core.database import (
    IndexedDatabase,
    SageSearchConfiguration,
    EnzymeBuilder,
)
from sagepy.core.ion_series import IonType


# Simple test FASTA with a few proteins
TEST_FASTA = """>sp|P00001|TEST1 Test protein 1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR
>sp|P00002|TEST2 Test protein 2
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVITKDEAEKLFNQDVDAAVRGILR
"""


class TestIndexedDatabaseSaveLoad:
    """Tests for IndexedDatabase persistence."""

    def test_save_and_load_roundtrip(self):
        """Test that saving and loading preserves the database."""
        # Create a database
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder(
                missed_cleavages=1,
                min_len=7,
                max_len=25,
                cleave_at="KR",
                restrict="P",
                c_terminal=True,
            ),
            generate_decoys=True,
            static_mods={"C": "[UNIMOD:4]"},
            variable_mods={},
        )
        original_db = config.generate_indexed_database()

        with tempfile.NamedTemporaryFile(suffix=".sagdb", delete=False) as f:
            path = f.name

        try:
            # Save and reload
            original_db.save(path)
            loaded_db = IndexedDatabase.load(path)

            # Verify properties match
            assert loaded_db.num_peptides == original_db.num_peptides
            assert loaded_db.num_fragments == original_db.num_fragments
            assert loaded_db.bucket_size == original_db.bucket_size
            assert loaded_db.generate_decoys == original_db.generate_decoys
            assert loaded_db.decoy_tag == original_db.decoy_tag

            # Verify peptide sequences match exactly
            original_seqs = original_db.peptides_as_string()
            loaded_seqs = loaded_db.peptides_as_string()
            assert len(original_seqs) == len(loaded_seqs)
            assert all(o == l for o, l in zip(original_seqs, loaded_seqs))

            # Verify masses match
            original_masses = original_db.mono_masses()
            loaded_masses = loaded_db.mono_masses()
            assert np.allclose(original_masses, loaded_masses)

            # Verify fragment indices match (critical for intensity predictions)
            assert np.array_equal(original_db.fragment_indices, loaded_db.fragment_indices)
            assert np.allclose(original_db.fragment_mzs, loaded_db.fragment_mzs)

        finally:
            Path(path).unlink()

    def test_save_and_load_preserves_peptide_order(self):
        """Test that peptide indices are stable across save/load."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder(
                missed_cleavages=2,
                min_len=6,
                max_len=30,
                cleave_at="KR",
                restrict="P",
                c_terminal=True,
            ),
            generate_decoys=True,
        )
        original_db = config.generate_indexed_database()

        with tempfile.NamedTemporaryFile(suffix=".sagdb", delete=False) as f:
            path = f.name

        try:
            original_db.save(path)
            loaded_db = IndexedDatabase.load(path)

            # Check that indexing by peptide index returns the same peptide
            for i in range(min(10, original_db.num_peptides)):
                orig_pep = original_db[i]
                loaded_pep = loaded_db[i]
                assert orig_pep.sequence == loaded_pep.sequence
                assert orig_pep.decoy == loaded_pep.decoy
                assert pytest.approx(orig_pep.mono_isotopic, rel=1e-5) == loaded_pep.mono_isotopic

        finally:
            Path(path).unlink()

    def test_load_nonexistent_file_raises(self):
        """Test that loading a non-existent file raises an error."""
        with pytest.raises(IOError):
            IndexedDatabase.load("/nonexistent/path/database.sagdb")

    def test_load_invalid_file_raises(self):
        """Test that loading an invalid file raises an error."""
        with tempfile.NamedTemporaryFile(suffix=".sagdb", delete=False) as f:
            f.write(b"INVALID_CONTENT")
            path = f.name

        try:
            with pytest.raises(IOError):
                IndexedDatabase.load(path)
        finally:
            Path(path).unlink()

    def test_save_creates_file(self):
        """Test that save creates the file at the specified path."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder.default_trypsin(),
            generate_decoys=False,
        )
        db = config.generate_indexed_database()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test_database.sagdb"
            assert not path.exists()

            db.save(str(path))
            assert path.exists()
            assert path.stat().st_size > 0

    def test_roundtrip_with_modifications(self):
        """Test save/load with static and variable modifications."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder.default_trypsin(),
            generate_decoys=True,
            static_mods={"C": "[UNIMOD:4]"},
            variable_mods={"M": ["[UNIMOD:35]"]},
            max_variable_mods=2,
        )
        original_db = config.generate_indexed_database()

        with tempfile.NamedTemporaryFile(suffix=".sagdb", delete=False) as f:
            path = f.name

        try:
            original_db.save(path)
            loaded_db = IndexedDatabase.load(path)

            # Check modifications are preserved
            orig_mods = original_db.modifications()
            loaded_mods = loaded_db.modifications()
            assert len(orig_mods) == len(loaded_mods)

            for (orig_idx, orig_mod), (loaded_idx, loaded_mod) in zip(orig_mods, loaded_mods):
                assert orig_idx == loaded_idx
                assert np.allclose(orig_mod, loaded_mod)

        finally:
            Path(path).unlink()

    def test_roundtrip_preserves_ion_kinds(self):
        """Test that ion kinds are preserved across save/load."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder.default_trypsin(),
            generate_decoys=False,
            ion_kinds=[IonType.b(), IonType.y()],
        )
        original_db = config.generate_indexed_database()

        with tempfile.NamedTemporaryFile(suffix=".sagdb", delete=False) as f:
            path = f.name

        try:
            original_db.save(path)
            loaded_db = IndexedDatabase.load(path)

            orig_kinds = [k.to_string() for k in original_db.ion_kinds]
            loaded_kinds = [k.to_string() for k in loaded_db.ion_kinds]
            assert orig_kinds == loaded_kinds

        finally:
            Path(path).unlink()


class TestIndexedDatabaseBasics:
    """Basic tests for IndexedDatabase functionality."""

    def test_database_creation(self):
        """Test that a database can be created from configuration."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder.default_trypsin(),
            generate_decoys=True,
        )
        db = config.generate_indexed_database()

        assert db.num_peptides > 0
        assert db.num_fragments > 0
        assert db.generate_decoys is True

    def test_peptide_indexing(self):
        """Test that peptides can be indexed by integer."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder.default_trypsin(),
            generate_decoys=False,
        )
        db = config.generate_indexed_database()

        peptide = db[0]
        assert peptide is not None
        assert len(peptide.sequence) > 0

    def test_peptides_as_string(self):
        """Test peptides_as_string returns numpy array."""
        config = SageSearchConfiguration(
            fasta=TEST_FASTA,
            enzyme_builder=EnzymeBuilder.default_trypsin(),
            generate_decoys=False,
        )
        db = config.generate_indexed_database()

        sequences = db.peptides_as_string()
        assert isinstance(sequences, np.ndarray)
        assert len(sequences) == db.num_peptides
