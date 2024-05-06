from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray

import sagepy_connector
psc = sagepy_connector.py_utility


def mass_to_mod(mass: float) -> str:
    """ Convert a mass to a UNIMOD modification annotation.

    Args:
        mass: a mass in Da

    Returns:
        a UNIMOD modification annotation
    """
    maybe_key = int(np.round(mass))
    # TODO: find a better way to do the map-back
    mod_dict = {
        42: '[UNIMOD:1]',
        57: '[UNIMOD:4]',
        80: '[UNIMOD:21]',
        16: '[UNIMOD:35]',
        119: '[UNIMOD:312]',
    }
    # try to translate to UNIMOD annotation
    try:
        return mod_dict[maybe_key]
    except KeyError:
        raise KeyError(f"Rounded mass not in dict: {maybe_key}")


def prosit_intensities_to_fragments_map(intensities: NDArray) -> Dict[Tuple[int, int, int], float]:
    """ Convert a Prosit intensity array to a fragment map (ion_type, charge, ordinal) -> intensity.

    Args:
        intensities: a Prosit intensity array

    Returns:
        a fragment map (ion_type, charge, ordinal) -> intensity
    """
    fragment_map = psc.flat_prosit_array_to_fragments_map(intensities)
    return fragment_map


def py_fragments_to_fragments_map(fragments, normalize: bool = True) -> Dict[Tuple[int, int, int], float]:
    """ Convert a Fragments object to a fragment map (ion_type, charge, ordinal) -> intensity.

    Args:
        normalize: whether to normalize the intensities
        fragments: a Fragments object

    Returns:
        a fragment map (ion_type, charge, ordinal) -> intensity
    """
    return psc.py_fragments_to_fragments_map(fragments.get_py_ptr(), normalize)
