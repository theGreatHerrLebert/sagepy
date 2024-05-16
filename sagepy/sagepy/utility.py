from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
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


def peptide_spectrum_match_list_to_pandas(
        psms,
        re_score: bool = False,
        use_sequence_as_match_idx: bool = True) -> pd.DataFrame:
    """Convert a list of peptide spectrum matches to a pandas dataframe

    Args:
        psms (List[PeptideSpectrumMatch]): The peptide spectrum matches
        re_score (bool, optional): Should re-score be used. Defaults to False.
        use_sequence_as_match_idx (bool, optional): Should the sequence be used as the match index. Defaults to True.

    Returns:
        pd.DataFrame: The pandas dataframe
    """
    row_list = []
    for match in psms:
        if match.retention_time_predicted is not None:
            delta_rt = match.retention_time_predicted - match.retention_time_observed
        else:
            delta_rt = None

        if re_score:
            score = match.re_score
        else:
            score = match.hyper_score

        if use_sequence_as_match_idx:
            match_idx = match.sequence
        else:
            match_idx = str(match.peptide_idx)

        row_list.append({
            "spec_idx": match.spec_idx,
            "match_idx": match_idx,
            "proteins": match.proteins,
            "decoy": match.decoy,
            "score": score,
            "re_score": match.re_score,
            "hyper_score": match.hyper_score,
            "rank": match.rank,
            "mono_mz_calculated": match.mono_mz_calculated,
            "mono_mass_observed": match.mono_mass_observed,
            "mono_mass_calculated": match.mono_mass_calculated,
            "delta_mass": match.mono_mass_calculated - match.mono_mass_observed,
            "isotope_error": match.isotope_error,
            "average_ppm": match.average_ppm,
            "delta_next": match.delta_next,
            "delta_best": match.delta_best,
            "matched_peaks": match.matched_peaks,
            "longest_b": match.longest_b,
            "longest_y": match.longest_y,
            "longest_y_pct": match.longest_y_pct,
            "missed_cleavages": match.missed_cleavages,
            "matched_intensity_pct": match.matched_intensity_pct,
            "scored_candidates": match.scored_candidates,
            "poisson": match.poisson,
            "sequence": match.sequence,
            "charge": match.charge,
            "retention_time_observed": match.retention_time_observed,
            "retention_time_predicted": match.retention_time_predicted,
            "delta_rt": delta_rt,
            "inverse_mobility_observed": match.inverse_mobility_observed,
            "inverse_mobility_predicted": match.inverse_mobility_predicted,
            "delta_ims": match.inverse_mobility_predicted - match.inverse_mobility_observed,
            "intensity_ms1": match.intensity_ms1,
            "intensity_ms2": match.intensity_ms2,
            "q_value": match.q_value,
            "collision_energy": match.collision_energy,
            "cosine_similarity": match.cosine_similarity,
        })

    return pd.DataFrame(row_list)


def get_features(ds: pd.DataFrame) -> (NDArray, NDArray):
    features = [
        "score", "delta_rt", "delta_ims", "cosine_similarity", "delta_mass",
        "rank", "isotope_error", "average_ppm", "delta_next", "delta_best",
        "matched_peaks", "longest_b", "longest_y", "longest_y_pct", "missed_cleavages",
        "matched_intensity_pct", "scored_candidates", "poisson", "charge",
        "intensity_ms1", "intensity_ms2", "collision_energy"
    ]

    X = ds[features].to_numpy().astype(np.float32)
    Y = ds["decoy"].to_numpy()
    Y = np.array([0 if x else 1 for x in Y]).astype(np.float32)

    return X, Y
