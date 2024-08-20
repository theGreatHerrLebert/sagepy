from typing import Dict, Union, List, Tuple, Set
from sagepy.core.modification import ModificationSpecificity

from sagepy_connector import py_unimod as unimod
from .modification import validate_mods, validate_var_mods


def unimod_static_mods_to_sage_static_mods(
        unimod_static_mods: Union[Dict[str, str], Dict[str, int]]
) -> Tuple[Dict[ModificationSpecificity, float], Set[str]]:
    """ Translate a dict that maps modification names to Unimod IDs
     to a dict that maps ModificationSpecificity objects and a set of modification names.
     Args:
         unimod_static_mods: A dict that maps modification names to Unimod IDs.
   Returns:
       A tuple containing a dict that maps ModificationSpecificity objects
       to mass values and a set of modification names.
    """
    mods_numeric = type(list(unimod_static_mods.values())[0]) is int
    if mods_numeric:
        unimod_to_mass = unimod.unimod_modification_to_mass_numerical()
    else:
        unimod_to_mass = unimod.unimod_modification_to_mass()

    sage_raw_dict = {}
    static_modifications = []

    for key, value in unimod_static_mods.items():
        mass = unimod_to_mass[value]
        sage_raw_dict[key] = mass

        if mods_numeric:
            static_modifications.append("[UNIMOD:" + f"{str(value)}" + "]")
        else:
            static_modifications.append(value)

    return validate_mods(sage_raw_dict), set(static_modifications)


def unimod_variable_mods_to_sage_variable_mods(
        unimod_variable_mods: Union[Dict[str, str], Dict[str, int]]
) -> Tuple[Dict[ModificationSpecificity, List[float]], Set[str]]:
    """ Translate a dict that maps modification names to Unimod IDs
    to a dict that maps ModificationSpecificity objects to lists of mass values and a set of modification names.

    Args:
        unimod_variable_mods: A dict that maps modification names to Unimod IDs.

    Returns:
        A tuple containing a dict that maps ModificationSpecificity objects
        to lists of mass values and a set of modification names.
    """
    mods_numeric = type(list(unimod_variable_mods.values())[0]) is int

    if mods_numeric:
        unimod_to_mass = unimod.unimod_modification_to_mass_numerical()
    else:
        unimod_to_mass = unimod.unimod_modification_to_mass()

    sage_raw_dict: Dict[str, float] = {}
    variable_modifications = []

    for key, value in unimod_variable_mods.items():
        mass = unimod_to_mass[value]

        if key in sage_raw_dict:
            sage_raw_dict[key].append(mass)
        else:
            sage_raw_dict[key] = [mass]

        if mods_numeric:
            variable_modifications.append("[UNIMOD:" + f"{str(value)}" + "]")
        else:
            variable_modifications.append(value)

    return validate_var_mods(sage_raw_dict), set(variable_modifications)
