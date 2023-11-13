import numpy as np
from numpy.typing import NDArray


def one_over_reduced_mobility_to_ccs(one_over_k0, mz, charge, mass_gas=28.013, temp=31.85, t_diff=273.15):
    """
    convert reduced ion mobility (1/k0) to CCS
    :param one_over_k0: reduced ion mobility
    :param charge: charge state of the ion
    :param mz: mass-over-charge of the ion
    :param mass_gas: mass of drift gas
    :param temp: temperature of the drift gas in C°
    :param t_diff: factor to translate from C° to K
    """
    SUMMARY_CONSTANT = 18509.8632163405
    reduced_mass = (mz * charge * mass_gas) / (mz * charge + mass_gas)
    return (SUMMARY_CONSTANT * charge) / (np.sqrt(reduced_mass * (temp + t_diff)) * 1 / one_over_k0)


def ccs_to_one_over_reduced_mobility(ccs, mz, charge, mass_gas=28.013, temp=31.85, t_diff=273.15):
    """
    convert CCS to 1 over reduced ion mobility (1/k0)
    :param ccs: collision cross-section
    :param charge: charge state of the ion
    :param mz: mass-over-charge of the ion
    :param mass_gas: mass of drift gas (N2)
    :param temp: temperature of the drift gas in C°
    :param t_diff: factor to translate from C° to K
    """
    SUMMARY_CONSTANT = 18509.8632163405
    reduced_mass = (mz * charge * mass_gas) / (mz * charge + mass_gas)
    return ((np.sqrt(reduced_mass * (temp + t_diff))) * ccs) / (SUMMARY_CONSTANT * charge)


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


def tokenize_sage_sequence(sequence: str, modifications: NDArray) -> list[str]:
    """ Tokenize a sequence with modifications into a list of tokens that are ionmob compatible.

    Args:
        sequence: a string of amino acids
        modifications: a numpy array of modifications, where 0 is no modification and mods are given as masses

    Returns:
        a list of tokens that are ionmob compatible
    """
    assert len(sequence) == len(modifications), "Sequence and modifications list need to be same length."
    seq_list = []
    start, end = ['<START>'], ['<END>']
    for char, mod in zip(sequence, modifications):
        if mod == 0:
            seq_list.append(char)
        else:
            seq_list.append(char + mass_to_mod(mod))

    return start + seq_list + end
