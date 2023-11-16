from typing import Union

import numpy as np
import sagepy_connector
psc = sagepy_connector.py_enzyme


class Position:
    def __init__(self, position: Union[str, None] = None):
        """Position class representing location of peptide in protein

        Args:
            position (Union[str, None], optional): The position of the peptide in the protein. Defaults to None.
        """
        if position is not None:
            try:
                self.__position_ptr = psc.PyPosition.from_string(position)
            except ValueError:
                raise ValueError("Invalid position string, allowed values are: c_term, n_term, full, internal")
        else:
            self.__position_ptr = psc.PyPosition.from_string('internal')

    @classmethod
    def from_py_position(cls, position: psc.PyPosition):
        instance = cls.__new__(cls)
        instance.__position_ptr = position
        return instance

    @property
    def position(self):
        return self.__position_ptr.to_string

    def __repr__(self):
        return f"Position({self.position})"

    def get_py_ptr(self):
        return self.__position_ptr


class Enzyme:
    def __init__(self, cleave_pattern: str = 'KR', c_terminal: bool = True, skip_suffix: str = 'P', semi_enzymatic: bool = False):
        """Enzyme class, default enzyme is Trypsin

        Args:
            cleave_pattern (str, optional): The cleavage pattern of the enzyme. Defaults to 'KR'.
            c_terminal (bool, optional): Cleave from the C-terminal. Defaults to True.
            skip_suffix (str, optional): The skip suffix of the enzyme. Defaults to 'P'.
            semi_enzymatic (bool, optional): Is the enzyme semi enzymatic. Defaults to False.
        """
        self.__enzyme_ptr = psc.PyEnzyme(cleave_pattern, c_terminal, skip_suffix, semi_enzymatic)

    @classmethod
    def from_py_enzyme(cls, enzyme: psc.PyEnzyme):
        instance = cls.__new__(cls)
        instance.__enzyme_ptr = enzyme
        return instance

    @property
    def c_terminal(self):
        if self.__enzyme_ptr is None:
            return None
        return self.__enzyme_ptr.c_terminal

    @property
    def skip_suffix(self):
        if self.__enzyme_ptr is None:
            return None
        return self.__enzyme_ptr.skip_suffix

    def cleavage_sites(self, sequence: str):
        if self.__enzyme_ptr is None:
            return None
        return self.__enzyme_ptr.cleavage_sites(sequence)

    def cleave(self, sequence: str, min_length: int = 1, max_length: int = np.inf):
        if self.__enzyme_ptr is None:
            raise ValueError("Enzyme is not defined")
        cleave_sites = self.__enzyme_ptr.cleavage_sites(sequence)

        cleaved_peptides = [sequence[start:end] for start, end in cleave_sites]
        filtered_peptides = [peptide for peptide in cleaved_peptides if min_length <= len(peptide) <= max_length]

        return filtered_peptides

    def get_py_ptr(self):
        return self.__enzyme_ptr

    def __repr__(self):
        return f"Enzyme(c_terminal: {self.c_terminal}, skip_suffix: {self.skip_suffix})"


class Digest:
    def __init__(self, decoy: bool, sequence: str, protein: str, missed_cleavages: int, position: Position):
        """Digest class representing a peptide digest

        Args:
            decoy (bool): Is the digest peptide a decoy
            sequence (str): The sequence of the digest peptide
            protein (str): The protein that the digest peptide is found in
            missed_cleavages (int): The number of missed cleavages
            position (Position): The position of the digest peptide in the protein
        """
        self.__digest_ptr = psc.PyDigest(decoy, sequence, protein, missed_cleavages, position.get_py_ptr())

    @classmethod
    def from_py_digest(cls, digest: psc.PyDigest):
        instance = cls.__new__(cls)
        instance.__digest_ptr = digest
        return instance

    @property
    def decoy(self):
        return self.__digest_ptr.decoy

    @property
    def sequence(self):
        return self.__digest_ptr.sequence

    @property
    def protein(self):
        return self.__digest_ptr.protein

    @property
    def missed_cleavages(self):
        return self.__digest_ptr.missed_cleavages

    @property
    def position(self):
        return self.__digest_ptr.position

    def reverse(self):
        return self.__digest_ptr.reverse()

    def __eq__(self, other):
        return self.__digest_ptr == other.__digest_ptr

    def __hash__(self):
        return hash(self.__digest_ptr)

    def __repr__(self):
        return f"Digest(decoy: {self.decoy}, protein: {self.protein}, " \
               f"missed_cleavages: {self.missed_cleavages}, position: {self.position}, sequence: {self.sequence})"

    def get_py_ptr(self):
        return self.__digest_ptr


class EnzymeParameters:
    def __init__(self, missed_cleavages: int = 0, min_len: int = 5, max_len: int = 50, enzyme: Enzyme = None):
        """EnzymeParameters class representing the parameters of an enzyme

        Args:
            missed_cleavages (int, optional): The number of missed cleavages. Defaults to 0.
            min_len (int, optional): The minimum length of a peptide. Defaults to 5.
            max_len (int, optional): The maximum length of a peptide. Defaults to 50.
            enzyme (Enzyme, optional): The enzyme. Defaults to None.
        """
        if enzyme is None:
            self.__enzyme_parameters_ptr = psc.PyEnzymeParameters(missed_cleavages, min_len, max_len, None)
        else:
            self.__enzyme_parameters_ptr = psc.PyEnzymeParameters(missed_cleavages, min_len, max_len,
                                                                  enzyme.get_py_ptr())

    @classmethod
    def from_py_enzyme_parameters(cls, enzyme_parameters: psc.PyEnzymeParameters):
        instance = cls.__new__(cls)
        instance.__enzyme_parameters_ptr = enzyme_parameters
        return instance

    @property
    def missed_cleavages(self):
        return self.__enzyme_parameters_ptr.missed_cleavages

    @property
    def min_len(self):
        return self.__enzyme_parameters_ptr.min_len

    @property
    def max_len(self):
        return self.__enzyme_parameters_ptr.max_len

    @property
    def enzyme(self):
        return Enzyme.from_py_enzyme(self.__enzyme_parameters_ptr.enzyme)

    def __repr__(self):
        return f"EnzymeParameters(missed_cleavages: {self.missed_cleavages}, min_len: {self.min_len}, " \
               f"max_len: {self.max_len}, enzyme: {self.enzyme})"

    def cleavage_sites(self, sequence: str):
        return self.__enzyme_parameters_ptr.cleavage_sites(sequence)

    def digest(self, sequence: str, protein: str):
        return [Digest.from_py_digest(s) for s in self.__enzyme_parameters_ptr.digest(sequence, protein)]

    def get_py_ptr(self):
        return self.__enzyme_parameters_ptr