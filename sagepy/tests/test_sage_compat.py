"""Compatibility tests copied from upstream sage behavior."""

from sagepy.core.database import EnzymeBuilder, SageSearchConfiguration
from sagepy.core.enzyme import Enzyme, EnzymeParameters
from sagepy.core.fasta import Fasta


def test_enzyme_trypsin_matches_sage():
    sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN"
    expected = [
        ("MADEEK", "Nterm"),
        ("LPPGWEK", "Internal"),
        ("MSR", "Internal"),
        ("SSGR", "Internal"),
        ("VYYFNHITNASQWERPSGN", "Cterm"),
    ]

    params = EnzymeParameters(
        missed_cleavages=0,
        min_len=2,
        max_len=50,
        enzyme=Enzyme("KR", True, "P", False),
    )

    observed = [(digest.sequence, digest.position) for digest in params.digest(sequence, "")]
    assert observed == expected


def test_enzyme_trypsin_missed_cleavage_matches_sage():
    sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN"
    expected = [
        "MADEEK",
        "LPPGWEK",
        "R",
        "MSR",
        "SSGR",
        "VYYFNHITNASQWERPSGN",
        "MADEEKLPPGWEK",
        "LPPGWEKR",
        "RMSR",
        "MSRSSGR",
        "SSGRVYYFNHITNASQWERPSGN",
    ]

    params = EnzymeParameters(
        missed_cleavages=1,
        min_len=0,
        max_len=50,
        enzyme=Enzyme("KR", True, "P", False),
    )

    observed = [digest.sequence for digest in params.digest(sequence, "")]
    assert observed == expected


def test_enzyme_nonspecific_windows_match_sage():
    sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW"
    expected = [
        sequence[start : start + window]
        for window in range(5, 8)
        for start in range(len(sequence) - window + 1)
    ]

    params = EnzymeParameters(
        missed_cleavages=0,
        min_len=5,
        max_len=7,
        enzyme=None,
    )

    observed = [digest.sequence for digest in params.digest(sequence, "")]
    assert observed == expected


def test_fasta_digest_matches_sage_database():
    fasta = Fasta(
        """
>sp|AAAAA
MEWKLEQSMREQALLKAQLTQLK
>sp|BBBBB
RMEWKLEQSMREQALLKAQLTQLK
""",
        decoy_tag="rev_",
        generate_decoys=False,
    )

    params = EnzymeParameters(
        missed_cleavages=1,
        min_len=6,
        max_len=10,
        enzyme=Enzyme("KR", True, "P", False),
    )

    observed = [
        (digest.sequence, digest.protein, digest.missed_cleavages, digest.position)
        for digest in fasta._digest(params)
    ]
    expected = [
        ("LEQSMR", "sp|AAAAA", 0, "Internal"),
        ("EQALLK", "sp|AAAAA", 0, "Internal"),
        ("AQLTQLK", "sp|AAAAA", 0, "Cterm"),
        ("MEWKLEQSMR", "sp|AAAAA", 1, "Nterm"),
        ("LEQSMR", "sp|BBBBB", 0, "Internal"),
        ("EQALLK", "sp|BBBBB", 0, "Internal"),
        ("AQLTQLK", "sp|BBBBB", 0, "Cterm"),
        ("MEWKLEQSMR", "sp|BBBBB", 1, "Internal"),
    ]
    assert observed == expected


def test_database_digest_matches_sage_database():
    fasta = """
>sp|AAAAA
MEWKLEQSMREQALLKAQLTQLK
>sp|BBBBB
RMEWKLEQSMREQALLKAQLTQLK
"""

    config = SageSearchConfiguration(
        fasta=fasta,
        bucket_size=128,
        enzyme_builder=EnzymeBuilder(
            missed_cleavages=1,
            min_len=6,
            max_len=10,
            cleave_at="KR",
            restrict="P",
            c_terminal=True,
        ),
        peptide_min_mass=150.0,
        peptide_max_mass=5000.0,
        min_ion_index=2,
        static_mods={},
        variable_mods={},
        max_variable_mods=2,
        decoy_tag="rev_",
        generate_decoys=False,
    )

    peptides = config._digest()
    observed = [peptide.sequence for peptide in peptides]
    expected = ["EQALLK", "LEQSMR", "AQLTQLK", "MEWKLEQSMR"]

    assert observed == expected
    for peptide in peptides:
        assert len(peptide.proteins) == 2
