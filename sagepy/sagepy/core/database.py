import numpy as np

from typing import List, Dict, Tuple, Union

import pandas as pd

from sagepy.core.mass import Tolerance
from sagepy.core.peptide import Peptide

from sagepy.core.enzyme import EnzymeParameters
import sagepy_connector

from sagepy.core.ion_series import IonType
from sagepy.core.modification import ModificationSpecificity
from sagepy.core.unimod import unimod_variable_mods_to_sage_variable_mods, unimod_static_mods_to_sage_static_mods
from .modification import process_variable_start_end_mods

psc = sagepy_connector.py_database


class PeptideIx:
    def __init__(self, idx: int):
        """ PeptideIx class

        Args:
            idx (int): idx in the database
        """
        self.__peptide_ix_ptr = psc.PyPeptideIx(idx)

    @classmethod
    def from_py_peptide_ix(cls, ix: psc.PyPeptideIx) -> 'PeptideIx':
        instance = cls.__new__(cls)
        instance.__peptide_ix_ptr = ix
        return instance

    @property
    def idx(self) -> int:
        return self.__peptide_ix_ptr.idx

    def __repr__(self) -> str:
        return f"PeptideIx({self.idx})"

    def get_py_ptr(self):
        return self.__peptide_ix_ptr


class Theoretical:
    def __init__(self, idx: int, fragment_mz: float):
        """Theoretical class

        Args:
            idx (int): idx in the database
            fragment_mz (float): fragment m/z
        """
        self.__theoretical_ptr = psc.PyTheoretical(idx, fragment_mz)

    @classmethod
    def from_py_theoretical(cls, theoretical: psc.PyTheoretical) -> 'Theoretical':
        instance = cls.__new__(cls)
        instance.__theoretical_ptr = theoretical
        return instance

    @property
    def idx(self):
        return PeptideIx.from_py_peptide_ix(self.__theoretical_ptr.idx)

    @property
    def fragment_mz(self) -> float:
        return self.__theoretical_ptr.fragment_mz

    def __repr__(self) -> str:
        return f"Theoretical(idx={self.idx}, fragment_mz={self.fragment_mz})"

    def get_py_ptr(self):
        return self.__theoretical_ptr


class EnzymeBuilder:
    def __init__(self, missed_cleavages: int = None, min_len: int = None, max_len: int = None, cleave_at: str = None,
                 restrict: str = None, c_terminal: bool = None, semi_enzymatic: bool = None):
        """EnzymeBuilder class

        Args:
            missed_cleavages (int, optional): Number of missed cleavages. Defaults to None.
            min_len (int, optional): Minimum peptide length. Defaults to None.
            max_len (int, optional): Maximum peptide length. Defaults to None.
            cleave_at (str, optional): Cleavage pattern. Defaults to None.
            restrict (str, optional): Restriction pattern. Defaults to None.
            c_terminal (bool, optional): Cleavage at the C-terminal. Defaults to None.
            semi_enzymatic (bool, optional): Allow semi-enzymatic peptides. Defaults to None.
        """
        self.__enzyme_builder_ptr = psc.PyEnzymeBuilder(missed_cleavages, min_len, max_len, cleave_at, restrict,
                                                        c_terminal, semi_enzymatic)

    @staticmethod
    def default_trypsin() -> 'EnzymeBuilder':
        return EnzymeBuilder.from_py_enzyme_builder(psc.PyEnzymeBuilder.from_default_trypsin())

    @classmethod
    def from_py_enzyme_builder(cls, enzyme_builder: psc.PyEnzymeBuilder) -> 'EnzymeBuilder':
        instance = cls.__new__(cls)
        instance.__enzyme_builder_ptr = enzyme_builder
        return instance

    def get_enzyme_parameters(self) -> EnzymeParameters:
        return EnzymeParameters.from_py_enzyme_parameters(self.__enzyme_builder_ptr.get_enzyme_parameters())

    def get_py_ptr(self):
        return self.__enzyme_builder_ptr

    @property
    def missed_cleavages(self) -> int:
        return self.__enzyme_builder_ptr.missed_cleavages

    @property
    def min_len(self) -> int:
        return self.__enzyme_builder_ptr.min_len

    @property
    def max_len(self) -> int:
        return self.__enzyme_builder_ptr.max_len

    @property
    def cleave_at(self) -> str:
        return self.__enzyme_builder_ptr.cleave_at

    @property
    def restrict(self) -> str:
        return self.__enzyme_builder_ptr.restrict

    @property
    def c_terminal(self) -> bool:
        return self.__enzyme_builder_ptr.c_terminal

    @property
    def semi_enzymatic(self) -> bool:
        return self.__enzyme_builder_ptr.semi_enzymatic

    def __repr__(self):
        return f"EnzymeBuilder(missed_cleavages: {self.missed_cleavages}, min_len: {self.min_len}, " \
               f"max_len: {self.max_len}, cleave_at: {self.cleave_at}, restrict: {self.restrict}, " \
               f"c_terminal: {self.c_terminal}, semi_enzymatic: {self.semi_enzymatic})"


class SageSearchConfiguration:
    def __init__(self,
                 fasta: str,
                 bucket_size: int = 8192,
                 enzyme_builder: 'EnzymeBuilder' = EnzymeBuilder.default_trypsin(),
                 peptide_min_mass: float = 500,
                 peptide_max_mass: float = 5_000,
                 ion_kinds: List[IonType] = None,
                 min_ion_index: int = 2,
                 static_mods: Union[Dict[str, str], Dict[str, int]] = None,
                 variable_mods: Union[Dict[str, List[str]], Dict[str, List[int]]] = None,
                 max_variable_mods: int = 2,
                 decoy_tag: str = 'rev_',
                 generate_decoys: bool = True,
                 shuffle_decoys: Union[bool, None] = None,
                 keep_ends: Union[bool, None] = None,
                 prefilter: Union[bool, None] = None,
                 prefilter_chunk_size: Union[int, None] = None,
                 prefilter_low_memory: Union[bool, None] = None,
                 ):
        """SageSearchConfiguration class

        Args:
            fasta (str): The fasta file
            bucket_size (int, optional): The bucket size. Defaults to 8192.
            enzyme_builder (EnzymeBuilder, optional): The enzyme builder. Defaults to EnzymeBuilder.default_trypsin().
            peptide_min_mass (float, optional): The minimum peptide mass. Defaults to 500.
            peptide_max_mass (float, optional): The maximum peptide mass. Defaults to 5000.
            ion_kinds (List[IonType], optional): The ion types. Defaults to None.
            min_ion_index (int, optional): The minimum ion index. Defaults to 2.
            static_mods (Dict[str, str] | Dict[str, int], optional): The static modifications given in UNIMOD notation. Defaults to None.
            variable_mods (Dict[str, List[str]] | Dict[str, List[int]], optional): The variable modifications given in UNIMOD notation. Defaults to None.
            max_variable_mods (int, optional): The maximum number of variable modifications. Defaults to 2.
            decoy_tag (str, optional): The decoy tag. Defaults to 'rev_'.
            generate_decoys (bool, optional): Whether to generate decoys. Defaults to True.
            shuffle_decoys (Union[bool, None], optional): Whether to shuffle decoys. Defaults to None.
            keep_ends (Union[bool, None], optional): Whether to include start and end amino acid for permutation strategy. Defaults to None.
            prefilter (Union[bool, None], optional): Whether to prefilter. Defaults to None.
            prefilter_chunk_size (Union[int, None], optional): The prefilter chunk size. Defaults to None.
            prefilter_low_memory (Union[bool, None], optional): Whether to use low memory for prefilter. Defaults to None.
        """

        if variable_mods is not None:
            # Process variable mods, expanding wildcard start and end modifications to all possible amino acids
            variable_mods = process_variable_start_end_mods(variable_mods)
            variable_mods = unimod_variable_mods_to_sage_variable_mods(variable_mods)

        if static_mods is not None:
            static_mods = unimod_static_mods_to_sage_static_mods(static_mods)

        self.__py_parameter_ptr = psc.PyParameters(
            find_next_power_of_2(bucket_size),
            enzyme_builder.get_py_ptr(),
            peptide_min_mass,
            peptide_max_mass,
            min_ion_index,
            {k.get_py_ptr(): v for k, v in static_mods.items()} if static_mods is not None else {},
            {k.get_py_ptr(): v for k, v in variable_mods.items()} if variable_mods is not None else {},
            max_variable_mods,
            decoy_tag,
            generate_decoys,
            fasta,
            prefilter,
            prefilter_chunk_size,
            prefilter_low_memory,
            ion_kinds,
            shuffle_decoys,
            keep_ends
        )

    @classmethod
    def from_py_parameters(cls, parameters: psc.PyParameters) -> 'SageSearchConfiguration':
        instance = cls.__new__(cls)
        instance.__py_parameter_ptr = parameters
        return instance

    def get_py_ptr(self):
        return self.__py_parameter_ptr

    def _digest(self):
        return [Peptide.from_py_peptide(p) for p in self.__py_parameter_ptr.digest()]

    def generate_indexed_database(self) -> 'IndexedDatabase':
        """Generate the indexed database

        Returns:
            IndexedDatabase: The indexed database
        """
        return IndexedDatabase.from_py_indexed_database(self.__py_parameter_ptr.build_indexed_database())

    @property
    def bucket_size(self):
        return self.__py_parameter_ptr.bucket_size

    @property
    def enzyme_builder(self):
        return EnzymeBuilder.from_py_enzyme_builder(self.__py_parameter_ptr.enzyme_builder)

    @property
    def peptide_min_mass(self):
        return self.__py_parameter_ptr.peptide_min_mass

    @property
    def peptide_max_mass(self):
        return self.__py_parameter_ptr.peptide_max_mass

    @property
    def ion_kinds(self):
        return [IonType.from_py_kind(k) for k in self.__py_parameter_ptr.ion_kinds]

    @property
    def min_ion_index(self):
        return self.__py_parameter_ptr.min_ion_index

    @property
    def static_mods(self):
        return {ModificationSpecificity.from_py_modification_specificity(k): v for k, v in self.__py_parameter_ptr.static_mods.items()}

    @property
    def variable_mods(self):
        return {ModificationSpecificity.from_py_modification_specificity(k): v for k, v in self.__py_parameter_ptr.variable_mods.items()}

    @property
    def max_variable_mods(self):
        return self.__py_parameter_ptr.max_variable_mods

    @property
    def decoy_tag(self):
        return self.__py_parameter_ptr.decoy_tag

    @property
    def generate_decoys(self):
        return self.__py_parameter_ptr.generate_decoys

    @property
    def fasta(self):
        return self.__py_parameter_ptr.fasta

    def __repr__(self):
        return f"SageSearchConfiguration(bucket_size: {self.bucket_size}, enzyme_builder: {self.enzyme_builder}, " \
               f"peptide_min_mass: {self.peptide_min_mass}, peptide_max_mass: {self.peptide_max_mass}, " \
               f"ion_kinds: {self.ion_kinds}, min_ion_index: {self.min_ion_index}, " \
               f"max_variable_mods: {self.max_variable_mods}, decoy_tag: {self.decoy_tag}, " \
               f"generate_decoys: {self.generate_decoys})"


class IndexedDatabase:
    def __init__(self, peptides: List[Peptide], fragments: List[Theoretical], ion_kinds: List[IonType],
                 min_value: List[float],
                 potential_mods: List[Tuple[(ModificationSpecificity, float)]], bucket_size: int,
                 generate_decoys: bool, decoy_tag: str):
        """IndexedDatabase class

        Args:
            peptides (List[Peptide]): The peptides in the database
            fragments (List[Theoretical]): The theoretical fragments in the database
            ion_kinds (List[IonType]): The ion types in the database
            min_value (List[float]): The minimum values in the database
            potential_mods (List[Tuple[(ModificationSpecificity, float)]]): The potential modifications in the database
            bucket_size (int): The bucket size
            generate_decoys (bool): Whether to generate decoys
            decoy_tag (str): The decoy tag
        """
        self.__indexed_database_ptr = psc.PyIndexedDatabase(
            [p.get_py_ptr() for p in peptides],
            [f.get_py_ptr() for f in fragments],
            [k.get_py_ptr() for k in ion_kinds],
            min_value,
            [(k.get_py_ptr(), v) for (k, v) in potential_mods],
            bucket_size,
            generate_decoys,
            decoy_tag,
        )

    @classmethod
    def from_py_indexed_database(cls, indexed_database: psc.PyIndexedDatabase):
        instance = cls.__new__(cls)
        instance.__indexed_database_ptr = indexed_database
        return instance

    def get_py_ptr(self):
        return self.__indexed_database_ptr

    def __getitem__(self, item: Union[int, PeptideIx]) -> Peptide:
        if isinstance(item, int):
            pep_ix = PeptideIx(item)
            return Peptide.from_py_peptide(self.__indexed_database_ptr[pep_ix.get_py_ptr()])
        elif isinstance(item, PeptideIx):
            return Peptide.from_py_peptide(self.__indexed_database_ptr[item.get_py_ptr()])
        else:
            raise ValueError(f"Invalid item type: {type(item)}")

    def query(self, precursor_mass: float, precursor_tolerance: Tolerance, fragment_tolerance: Tolerance):
        return IndexedQuery.from_py_indexed_query(self.__indexed_database_ptr.query(precursor_mass,
                                                                                    precursor_tolerance.get_py_ptr(),
                                                                                    fragment_tolerance.get_py_ptr()))

    @property
    def _peptides(self):
        """Peptides in the database
        CAUTION: This method is for debugging purposes only. Do not use in production, RAM usage will be very high.
        Returns:
            List[Peptide]: The peptides in the database
        """
        return [Peptide.from_py_peptide(p) for p in self.__indexed_database_ptr.peptides]

    def peptides_as_string(self):
        return np.array(self.__indexed_database_ptr.peptides_as_string())

    def mono_masses(self):
        return np.array(self.__indexed_database_ptr.mono_masses())

    def modifications(self):
        return self.__indexed_database_ptr.modifications()

    @property
    def df(self) -> pd.DataFrame:
        """Peptides in the database as a pandas DataFrame

        Returns:
            pd.DataFrame: The peptides in the database as a pandas DataFrame
        """

        data = {
            'sequence': self.peptides_as_string(),
            'monoisotopic_mass': self.mono_masses()
        }

        return pd.DataFrame(data)

    @property
    def num_peptides(self) -> int:
        """Number of peptides in the database

        Returns:
            int: Number of peptides in the database
        """
        return self.__indexed_database_ptr.num_peptides

    @property
    def _fragments(self):
        """Fragments in the database
        CAUTION: This method is for debugging purposes only. Do not use in production, RAM usage will be very high.
        Returns:
            List[Theoretical]: The fragments in the database
        """
        return [Theoretical.from_py_theoretical(f) for f in self.__indexed_database_ptr.fragments]

    @property
    def fragment_indices(self):
        return self.__indexed_database_ptr.fragment_indices

    @property
    def fragment_mzs(self):
        return self.__indexed_database_ptr.fragment_mzs

    def fragment_dict(self):
        return self.__indexed_database_ptr.fragment_dict()

    @property
    def num_fragments(self) -> int:
        """Number of fragments in the database

        Returns:
            int: Number of fragments in the database
        """
        return self.__indexed_database_ptr.num_fragments

    @property
    def ion_kinds(self):
        return [IonType.from_py_kind(k) for k in self.__indexed_database_ptr.ion_kinds]

    @property
    def min_value(self):
        return self.__indexed_database_ptr.min_value

    @property
    def potential_mods(self):
        return [(ModificationSpecificity.from_py_modification_specificity(k), v) for (k, v) in self.__indexed_database_ptr.potential_mods]

    @property
    def bucket_size(self):
        return self.__indexed_database_ptr.bucket_size

    @property
    def generate_decoys(self):
        return self.__indexed_database_ptr.generate_decoys

    @property
    def decoy_tag(self):
        return self.__indexed_database_ptr.decoy_tag

    def __repr__(self):
        return f"IndexedDatabase(peptides: {self.num_peptides}, fragments: {self.num_fragments}, ion_kinds: {self.ion_kinds}, " \
               f"num_buckets: {len(self.min_value)}, " \
               f"bucket_size: {self.bucket_size}, generate_decoys: {self.generate_decoys}, " \
               f"decoy_tag: {self.decoy_tag})"


class IndexedQuery:
    def __init__(self, precursor_mass: float, precursor_tolerance: Tolerance, fragment_tolerance: Tolerance,
                 pre_idx_lo: int, pre_idx_hi: int):
        """IndexedQuery class

        Args:
            precursor_mass (float): The precursor mass
            precursor_tolerance (Tolerance): The precursor tolerance
            fragment_tolerance (Tolerance): The fragment tolerance
            pre_idx_lo (int): The lower index
            pre_idx_hi (int): The upper index
        """
        self.__indexed_query_ptr = psc.PyIndexedQuery(precursor_mass, precursor_tolerance.get_py_ptr(),
                                                      fragment_tolerance.get_py_ptr(), pre_idx_lo, pre_idx_hi)

    @classmethod
    def from_py_indexed_query(cls, indexed_query: psc.PyIndexedQuery):
        instance = cls.__new__(cls)
        instance.__indexed_query_ptr = indexed_query
        return instance

    def get_py_ptr(self):
        return self.__indexed_query_ptr

    @property
    def precursor_mass(self):
        return self.__indexed_query_ptr.precursor_mass

    @property
    def precursor_tolerance(self):
        return Tolerance.from_py_tolerance(self.__indexed_query_ptr.precursor_tolerance)

    @property
    def fragment_tolerance(self):
        return Tolerance.from_py_tolerance(self.__indexed_query_ptr.fragment_tolerance)

    @property
    def pre_idx_lo(self):
        return self.__indexed_query_ptr.pre_idx_lo

    @property
    def pre_idx_hi(self):
        return self.__indexed_query_ptr.pre_idx_hi

    def __repr__(self):
        return f"IndexedQuery(precursor_mass: {self.precursor_mass}, precursor_tolerance: {self.precursor_tolerance}, " \
               f"fragment_tolerance: {self.fragment_tolerance}, pre_idx_lo: {self.pre_idx_lo}, " \
               f"pre_idx_hi: {self.pre_idx_hi})"


def find_next_power_of_2(n):
    """
    Find the next power of two greater than or equal to n.

    Args:
    n (int): A positive integer.

    Returns:
    int: The smallest power of two greater than or equal to n.
    """
    return 1 << (n-1).bit_length()
