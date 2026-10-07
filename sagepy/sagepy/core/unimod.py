from typing import Dict, Union, List
from sagepy.core.modification import ModificationSpecificity

from sagepy_connector import py_unimod as unimod
from .modification import validate_mods, validate_var_mods, split_variable_mod_entry


def modification_title_to_unimod_id() -> Dict[str, str]:
    """ Get a dict that maps modification names to Unimod IDs.

    Returns:
        A dict that maps modification names to Unimod IDs.
    """
    return unimod.title_to_unimod_ids()


def modification_atomic_composition() -> Dict[str, Dict[str, int]]:
    """ Get a dict that maps modification names to atomic compositions.

    Returns:
        A dict that maps modification names to atomic compositions.
    """
    return unimod.modification_atomic_compositions()


def unimod_to_mass() -> Dict[str, float]:
    """ Get a dict that maps Unimod IDs to mass values.

    Returns:
        A dict that maps Unimod IDs to mass values.
    """
    return unimod.unimod_modification_to_mass()


def unimod_to_mass_numerical() -> Dict[int, float]:
    """ Get a dict that maps Unimod IDs given as integer to mass values.

    Returns:
        A dict that maps Unimod IDs to mass values.
    """
    return unimod.unimod_modification_to_mass_numerical()


def unimod_static_mods_to_sage_static_mods(
        unimod_static_mods: Union[Dict[str, str], Dict[str, int]]
) -> Dict[ModificationSpecificity, float]:
    """ Translate a dict that maps modification names to Unimod IDs
     to a dict that maps ModificationSpecificity objects and a set of modification names.
     Args:
         unimod_static_mods: A dict that maps modification names to Unimod IDs.
   Returns:
       A tuple containing a dict that maps ModificationSpecificity objects
       to mass values and a set of modification names.
    """

    if len(unimod_static_mods) == 0:
        return {}

    mods_numeric = type(list(unimod_static_mods.values())[0]) is int
    if mods_numeric:
        mod_to_mass = unimod.unimod_modification_to_mass_numerical()
    else:
        mod_to_mass = unimod.unimod_modification_to_mass()

    sage_raw_dict = {}

    for key, value in unimod_static_mods.items():
        try:
            mass = mod_to_mass[value]
        except KeyError:
            # check if the value can be parsed as a float
            try:
                mass = float(value)
                print(f"Unimod ID {value} for static modification {key} not found. "
                              f"Interpreting as mass shift: {mass}. If this was intentional, be sure to manually map "
                                f"search results with a custom modification mapping to avoid issues during rescoring.")
            except ValueError:
                raise KeyError(f"Unimod ID {value} for modification {key} not found.")

        sage_raw_dict[key] = mass

    return validate_mods(sage_raw_dict)


def unimod_variable_mods_to_sage_variable_mods(
        unimod_variable_mods: Union[Dict[str, List[str]], Dict[str, List[int]]]
) -> Dict[ModificationSpecificity, List[Union[float, tuple]]]:
    """ Translate a dict that maps modification names to Unimod IDs
    to a dict that maps ModificationSpecificity objects to lists of mass values.

    Entries may also carry a per-peptide limit, as ``(unimod_id, max_count)`` or
    ``{"mod": unimod_id, "max_count": max_count}``; those become ``(mass, max_count)``.

    Args:
        unimod_variable_mods: A dict that maps modification names to Unimod IDs.

    Returns:
        A dict that maps ModificationSpecificity objects to lists of masses or (mass, max_count) tuples.
    """

    if len(unimod_variable_mods) == 0:
        return {}

    mod_to_mass_named = unimod.unimod_modification_to_mass()
    mod_to_mass_numeric = None

    sage_raw_dict: Dict[str, List[Union[float, tuple]]] = {}

    for key, values in unimod_variable_mods.items():
        for entry in values:
            value, max_count = split_variable_mod_entry(entry)
            if isinstance(value, int):
                if mod_to_mass_numeric is None:
                    mod_to_mass_numeric = unimod.unimod_modification_to_mass_numerical()
                mod_to_mass = mod_to_mass_numeric
            else:
                mod_to_mass = mod_to_mass_named
            try:
                mass = mod_to_mass[value]
            except KeyError:
                # check if the value can be parsed as a float
                try:
                    mass = float(value)
                    print(f"Unimod ID {value} for variable modification {key} not found. "
                                  f"Interpreting as mass shift: {mass}. If this was intentional, be sure to manually map "
                                  f"search results with a custom modification mapping to avoid issues during rescoring.")
                except ValueError:
                    raise KeyError(f"Unimod ID {value} for modification {key} not found.")

            sage_raw_dict.setdefault(key, []).append(mass if max_count is None else (mass, max_count))

    return validate_var_mods(sage_raw_dict)


def static_unimod_mods_to_set(
        unimod_mods: Union[Dict[str, str], Dict[str, int]]
) -> set:
    """ Translate a dict that maps modification names to Unimod IDs to a set of modification names.

    Args:
        unimod_mods: A dict that maps modification names to Unimod IDs.

    Returns:
        A set of modification names.
    """

    if len(unimod_mods) == 0:
        return set()

    if isinstance(next(iter(unimod_mods.values())), int):
        return {f"[UNIMOD:{value}]" for value in unimod_mods.values()}
    else:
        return set(unimod_mods.values())

def variable_unimod_mods_to_set(
        unimod_mods: Union[Dict[str, List[str]], Dict[str, List[int]]
    ]) -> set:
    """ Translate a dict that maps modification names to Unimod IDs to a set of modification names.

    Args:
        unimod_mods: A dict that maps modification names to Unimod IDs.

    Returns:
        A set of modification names.
    """

    mods = (split_variable_mod_entry(entry)[0] for values in unimod_mods.values() for entry in values)
    return {f"[UNIMOD:{mod}]" if isinstance(mod, int) else mod for mod in mods}
