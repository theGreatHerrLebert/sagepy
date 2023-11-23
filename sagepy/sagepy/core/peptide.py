from typing import List, Optional

import sagepy_connector

from sagepy.core.enzyme import Position, Digest
from sagepy.utility import mass_to_mod

psc = sagepy_connector.py_peptide


class Peptide:
    def __init__(self,
                 decoy: bool,
                 sequence: str,
                 modifications: List[float],
                 mono_isotopic: float,
                 missed_cleavages: int,
                 position: Position,
                 proteins: List[str],
                 semi_enzymatic: bool = False,
                 n_term: Optional[float] = None,
                 c_term: Optional[float] = None
                 ):
        """Peptide class

        Args:
            decoy (bool): Is the peptide a decoy
            sequence (str): The peptide sequence
            modifications (List[float]): The modifications of the peptide
            mono_isotopic (float): The monoisotopic mass of the peptide
            missed_cleavages (int): The number of missed cleavages
            position (Position): The position of the peptide
            proteins (List[str]): The proteins that the peptide is found in
            n_term (Optional[float], optional): Potential modifications on the N-terminal of the peptide. Defaults to None.
            c_term (Optional[float], optional): Potential modifications on the C-terminal of the peptide. Defaults to None.
        """
        self.__peptide_ptr = psc.PyPeptide(decoy, sequence, modifications,
                                           mono_isotopic, missed_cleavages, position.get_py_ptr(),
                                           proteins, semi_enzymatic, n_term, c_term)

    @classmethod
    def from_digest(cls, digest: Digest) -> 'Peptide':
        instance = cls.__new__(cls)
        instance.__peptide_ptr = psc.PyPeptide.try_new_from_digest(digest.get_py_ptr())
        return instance

    @classmethod
    def from_py_peptide(cls, peptide: psc.PyPeptide):
        instance = cls.__new__(cls)
        instance.__peptide_ptr = peptide
        return instance

    @property
    def decoy(self):
        return self.__peptide_ptr.decoy

    @property
    def sequence(self):
        return self.__peptide_ptr.sequence

    @property
    def modifications(self):
        return self.__peptide_ptr.modifications

    @property
    def mono_isotopic(self):
        return self.__peptide_ptr.monoisotopic

    @property
    def missed_cleavages(self):
        return self.__peptide_ptr.missed_cleavages

    @property
    def position(self):
        return Position.from_py_position(self.__peptide_ptr.position)

    @property
    def proteins(self):
        return self.__peptide_ptr.proteins

    @property
    def n_term(self):
        return self.__peptide_ptr.n_term

    @property
    def c_term(self):
        return self.__peptide_ptr.c_term

    @property
    def semi_enzymatic(self):
        return self.__peptide_ptr.semi_enzymatic

    def get_py_ptr(self):
        return self.__peptide_ptr

    def __repr__(self):
        return f"Peptide(decoy: {self.decoy}, sequence: {self.sequence}, " \
               f"modifications: {self.modifications}, mono_isotopic: {self.mono_isotopic}, " \
               f"missed_cleavages: {self.missed_cleavages}, position: {self.position}, " \
               f"proteins: {self.proteins}, semi_enzymatic: {self.semi_enzymatic}, n_term: {self.n_term}, " \
               f"c_term: {self.c_term})"

    def to_unimod_sequence(self) -> str:
        """ Get Peptide sequence with UNIMOD modification annotations.

        Returns:
            str: Peptide sequence with UNIMOD modification annotations.
        """

        mods = self.modifications
        sequence = self.sequence

        seq = ''

        for i, (s, m) in enumerate(zip(sequence, mods)):
            if m != 0:
                # TODO: check if this is the correct way to handle N- and C-terminal mods
                if i == 0:
                    if mass_to_mod(m) == '[UNIMOD:1]':
                        seq += f'{mass_to_mod(m)}{s}'
                    else:
                        seq += f'{s}{mass_to_mod(m)}'
                else:
                    seq += f'{s}{mass_to_mod(m)}'
            else:
                seq += s

        return seq
