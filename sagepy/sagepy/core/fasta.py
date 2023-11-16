from .enzyme import EnzymeParameters, Digest
import sagepy_connector
psc = sagepy_connector.py_fasta


class Fasta:
    def __init__(self, fasta_str: str, decoy_tag: str = 'decoy_', generate_decoys: bool = False):
        """Fasta class

        Args:
            fasta_str (str): The fasta string
            decoy_tag (str, optional): The decoy tag. Defaults to 'decoy_'.
            generate_decoys (bool, optional): Should decoys be generated. Defaults to False.
        """
        self.__fasta_ptr = psc.PyFasta.parse(fasta_str, decoy_tag, generate_decoys)

    @classmethod
    def from_py_fasta(cls, fasta: psc.PyFasta):
        instance = cls.__new__(cls)
        instance.__fasta_ptr = fasta
        return instance

    def _digest(self, enzyme_parameters: EnzymeParameters):
        return [Digest.from_py_digest(s) for s in self.__fasta_ptr.digest(enzyme_parameters.get_py_ptr())]
