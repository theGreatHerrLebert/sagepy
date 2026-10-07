"""Tests for per-modification variable mod limits and max_combinations."""
import pytest
from sagepy.core import SageSearchConfiguration, EnzymeBuilder
from sagepy.core.modification import split_variable_mod_entry
from sagepy.core.unimod import variable_unimod_mods_to_set

FASTA = ">sp|P1|X\nMSSTSPKMSSTSPKAAMMMMSSSKGGG\n"


def _config(**kwargs):
    return SageSearchConfiguration(
        fasta=FASTA,
        generate_decoys=False,
        peptide_min_mass=100,
        max_variable_mods=3,
        enzyme_builder=EnzymeBuilder(missed_cleavages=0, min_len=3, max_len=50,
                                     cleave_at="K", restrict="P", c_terminal=True),
        **kwargs,
    )


def _phospho_count(peptide):
    return sum(1 for m in peptide.modifications if abs(m - 79.966) < 0.01)


UNLIMITED = {"M": ["[UNIMOD:35]"], "S": ["[UNIMOD:21]"]}


class TestPerModLimit:
    def test_tuple_and_dict_entries_are_equivalent(self):
        as_tuple = _config(variable_mods={"M": ["[UNIMOD:35]"], "S": [("[UNIMOD:21]", 1)]})
        as_dict = _config(variable_mods={"M": ["[UNIMOD:35]"], "S": [{"mod": "[UNIMOD:21]", "max_count": 1}]})
        assert len(as_tuple._digest()) == len(as_dict._digest())

    def test_limit_is_enforced(self):
        unlimited = _config(variable_mods=UNLIMITED)._digest()
        limited = _config(variable_mods={"M": ["[UNIMOD:35]"], "S": [("[UNIMOD:21]", 1)]})._digest()
        assert max(_phospho_count(p) for p in unlimited) > 1
        assert max(_phospho_count(p) for p in limited) == 1
        assert len(limited) < len(unlimited)

    def test_limit_none_matches_bare_entry(self):
        bare = _config(variable_mods=UNLIMITED)._digest()
        none = _config(variable_mods={"M": ["[UNIMOD:35]"], "S": [("[UNIMOD:21]", None)]})._digest()
        assert len(bare) == len(none)

    def test_numeric_unimod_ids(self):
        named = _config(variable_mods={"S": [("[UNIMOD:21]", 1)]})._digest()
        numeric = _config(variable_mods={"S": [(21, 1)]})._digest()
        assert len(named) == len(numeric)

    def test_variable_mods_getter_roundtrip(self):
        config = _config(variable_mods={"M": ["[UNIMOD:35]"], "S": [("[UNIMOD:21]", 2)]})
        by_residue = {str(k): v for k, v in config.variable_mods.items()}
        assert isinstance(by_residue["ModificationSpecificity(M)"][0], float)
        mass, max_count = by_residue["ModificationSpecificity(S)"][0]
        assert max_count == 2 and abs(mass - 79.966) < 0.01


class TestMaxCombinations:
    def test_default_is_unlimited(self):
        assert _config(variable_mods=UNLIMITED).max_combinations is None

    def test_cap_reduces_variants(self):
        unlimited = len(_config(variable_mods=UNLIMITED)._digest())
        capped = _config(variable_mods=UNLIMITED, max_combinations=2)
        assert capped.max_combinations == 2
        assert len(capped._digest()) < unlimited

    def test_cap_of_one_keeps_only_unmodified(self):
        peptides = _config(variable_mods=UNLIMITED, max_combinations=1)._digest()
        assert all(all(m == 0 for m in p.modifications) for p in peptides)

    def test_zero_is_normalised_to_one(self):
        assert _config(variable_mods=UNLIMITED, max_combinations=0).max_combinations == 1


class TestSplitEntry:
    @pytest.mark.parametrize("entry,expected", [
        ("[UNIMOD:35]", ("[UNIMOD:35]", None)),
        (35, (35, None)),
        (("[UNIMOD:21]", 2), ("[UNIMOD:21]", 2)),
        ({"mod": "[UNIMOD:21]"}, ("[UNIMOD:21]", None)),
        ({"mod": 21, "max_count": 0}, (21, 0)),
    ])
    def test_valid(self, entry, expected):
        assert split_variable_mod_entry(entry) == expected

    @pytest.mark.parametrize("entry", [
        ("[UNIMOD:21]", -1),
        ("[UNIMOD:21]", 1.5),
        ("[UNIMOD:21]", True),
        ("[UNIMOD:21]", 1, 2),
        {"max_count": 1},
        {"mod": "[UNIMOD:21]", "limit": 1},
    ])
    def test_invalid(self, entry):
        with pytest.raises(ValueError):
            split_variable_mod_entry(entry)

    def test_mods_to_set_unwraps_entries(self):
        mods = {"S": [("[UNIMOD:21]", 1)], "M": [35]}
        assert variable_unimod_mods_to_set(mods) == {"[UNIMOD:21]", "[UNIMOD:35]"}
