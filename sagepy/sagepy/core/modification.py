from typing import Dict, List

import sagepy_connector

psc = sagepy_connector.py_modification


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
        return "[", 42.0,

    @staticmethod
    def protein_c_terminal_variable():
        return "]", 111.0


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


def validate_var_mods(mods: Dict[str, List[float]]) -> Dict[ModificationSpecificity, List[float]]:

    py_validate_dict = psc.py_validate_var_mods(mods)

    return {ModificationSpecificity.from_py_modification_specificity(k): v for k, v in py_validate_dict.items()}


if __name__ == "__main__":
    static_mods = {k: v for k, v in [SAGE_KNOWN_MODS.cysteine_static()]}
    variable_mods = {k: v for k, v in [SAGE_KNOWN_MODS.methionine_variable()]}

    static = validate_mods(static_mods)
    variab = validate_var_mods(variable_mods)
