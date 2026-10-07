from typing import Dict, List, Optional, Tuple, Union

import sagepy_connector

psc = sagepy_connector.py_modification


VarModEntry = Union[float, Tuple[float, Optional[int]]]


def split_variable_mod_entry(entry) -> Tuple[Union[str, int, float], Optional[int]]:
    """Split a variable modification entry into the modification and its per-peptide limit.

    An entry is either a bare modification (UNIMOD id/name or mass), a ``(mod, max_count)`` tuple,
    or a ``{"mod": mod, "max_count": max_count}`` dict. ``max_count`` limits how often this
    modification may occur on a single peptide; ``None`` means unrestricted (except by
    ``max_variable_mods``).

    Args:
        entry: The variable modification entry

    Returns:
        Tuple[Union[str, int, float], Optional[int]]: The modification and its max_count
    """
    if isinstance(entry, dict):
        unknown = set(entry) - {"mod", "max_count"}
        if "mod" not in entry or unknown:
            raise ValueError(f"Variable modification dict needs a 'mod' and optional 'max_count' key, got {entry}")
        mod, max_count = entry["mod"], entry.get("max_count")
    elif isinstance(entry, (tuple, list)):
        if len(entry) != 2:
            raise ValueError(f"Variable modification tuple must be (mod, max_count), got {entry}")
        mod, max_count = entry
    else:
        return entry, None

    if max_count is not None and (isinstance(max_count, bool) or not isinstance(max_count, int) or max_count < 0):
        raise ValueError(f"max_count must be a non-negative int or None, got {max_count!r}")
    return mod, max_count


def process_variable_start_end_mods(variable_modifications):
    """Helper function to process variable modifications for start and end of peptides/proteins
    For some reason, the variable modification wildcards are not processed correctly when passed to SAGE
    This function processes the variable modifications and adds the start and end wildcards for amino acids
    Args:
        variable_modifications: The variable modifications

    Returns:
        Dict: The processed variable modifications
    """

    # peptide C, peptide N, protein C, protein N
    targets = ["^", "$", "[", "]"]

    # combine targets with amino acids
    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

    ret_dict = {}

    for key, values in variable_modifications.items():
        if key in targets:
            for amino_acid in AMINO_ACIDS:
                ret_dict[key + amino_acid] = values

    return { **variable_modifications, **ret_dict }


class SAGE_KNOWN_MODS:
    @staticmethod
    def n_terminal_static():
        return "^", 304.207

    @staticmethod
    def lysine_static():
        return "K", 304.207

    @staticmethod
    def cysteine_static():
        return "C", 57.0215

    @staticmethod
    def phospho_serine_static():
        return "S", 79.9663

    @staticmethod
    def phospho_threonine_static():
        return "T", 79.9663

    @staticmethod
    def phospho_tyrosine_static():
        return "Y", 79.9663

    @staticmethod
    def phospho_serine_variable():
        return "S", [79.9663]

    @staticmethod
    def phospho_threonine_variable():
        return "T", [79.9663]

    @staticmethod
    def phospho_tyrosine_variable():
        return "Y", [79.9663]

    @staticmethod
    def methionine_variable():
        return "M", [15.9949]

    @staticmethod
    def q_variable():
        return "^Q", [-17.026549]

    @staticmethod
    def glutamic_acid_n_terminal_variable():
        return "^E", [-18.010565]

    @staticmethod
    def peptide_c_terminal_variable():
        return "$", [49.2, 22.9]

    @staticmethod
    def protein_n_terminus_variable():
        return "[", [42.0]

    @staticmethod
    def protein_c_terminal_variable():
        return "]", [111.0]

    def __repr__(self):
        return (f"SAGE_KNOWN_MODS({self.n_terminal_static()}, "
                f"{self.lysine_static()}, "
                f"{self.cysteine_static()}, "
                f"{self.methionine_variable()}, "
                f"{self.q_variable()}, "
                f"{self.glutamic_acid_n_terminal_variable()}, "
                f"{self.peptide_c_terminal_variable()}, "
                f"{self.protein_n_terminus_variable()}, "
                f"{self.protein_c_terminal_variable()})")


# TODO: need to re-implement based on constant modification list
class ModificationSpecificity:
    def __init__(self, s: str):
        self.__modification_specificity_ptr = psc.PyModificationSpecificity(s)

    @classmethod
    def from_py_modification_specificity(cls, specificity: psc.PyModificationSpecificity):
        instance = cls.__new__(cls)
        instance.__modification_specificity_ptr = specificity
        return instance

    def __repr__(self):
        return f"ModificationSpecificity({self.__modification_specificity_ptr.as_string})"

    def get_py_ptr(self):
        return self.__modification_specificity_ptr


def validate_mods(mods: Dict[str, float]) -> Dict[ModificationSpecificity, float]:

    py_validate_dict = psc.py_validate_mods(mods)

    return {ModificationSpecificity.from_py_modification_specificity(k): v for k, v in py_validate_dict.items()}


def validate_var_mods(mods: Dict[str, List[VarModEntry]]) -> Dict[ModificationSpecificity, List[VarModEntry]]:
    """Validate variable modifications. Each entry is a mass, or a ``(mass, max_count)`` tuple."""

    py_validate_dict = psc.py_validate_var_mods(mods)

    return {ModificationSpecificity.from_py_modification_specificity(k): v for k, v in py_validate_dict.items()}


if __name__ == "__main__":
    static_mods = {k: v for k, v in [SAGE_KNOWN_MODS.cysteine_static()]}
    variable_mods = {k: v for k, v in [SAGE_KNOWN_MODS.methionine_variable()]}

    static = validate_mods(static_mods)
    variab = validate_var_mods(variable_mods)
