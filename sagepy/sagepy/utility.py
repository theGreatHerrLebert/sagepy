from typing import Dict, Tuple, List, Optional

from numba import jit
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from sagepy.qfdr.tdc import target_decoy_competition_pandas

import sagepy_connector

psc = sagepy_connector.py_utility


@jit(nopython=True)
def calculate_ppm_error(measured_value, reference_value):
    ppm_error = ((measured_value - reference_value) / reference_value) * 1_000_000
    return ppm_error


@jit(nopython=True)
def calculate_ppms(measured_values, reference_values):
    n = len(measured_values)
    ppm_errors = np.empty(n, dtype=np.float64)
    for i in range(n):
        ppm_errors[i] = calculate_ppm_error(measured_values[i], reference_values[i])
    return ppm_errors


@jit(nopython=True)
def mean_ppm(mz_observed, mz_calculated) -> float:
    return np.mean(calculate_ppms(mz_observed, mz_calculated))


@jit(nopython=True)
def median_ppm(mz_observed, mz_calculated) -> float:
    return np.median(calculate_ppms(mz_observed, mz_calculated))


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
        if match.retention_time_predicted is not None and match.projected_rt is not None:
            delta_rt = match.retention_time_predicted - match.projected_rt
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

        if match.inverse_mobility_predicted is not None:
            delta_ims = match.inverse_mobility_predicted - match.inverse_mobility_observed
        else:
            delta_ims = None

        if match.beta_score is not None:
            beta_score = match.beta_score
        else:
            beta_score = None

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
            "delta_ims": delta_ims,
            "intensity_ms1": match.intensity_ms1,
            "intensity_ms2": match.intensity_ms2,
            "q_value": match.q_value,
            "collision_energy": match.collision_energy,
            "cosine_similarity": match.cosine_similarity,
            "mean_ppm": mean_ppm(match.fragments_observed.mz_experimental, match.fragments_observed.mz_calculated),
            "median_ppm": median_ppm(match.fragments_observed.mz_experimental, match.fragments_observed.mz_calculated),
            "fragments_observed": match.fragments_observed,
            "fragments_predicted": match.fragments_predicted,
            "projected_rt": match.projected_rt,
            "beta_score": beta_score
        })

    return pd.DataFrame(row_list)


def get_features(ds: pd.DataFrame, score: Optional[str] = None) -> (NDArray, NDArray):

    score = score if score is not None else "score"

    features = [
        f"{score}", "delta_rt", "delta_ims", "cosine_similarity", "delta_mass",
        "rank", "isotope_error", "average_ppm", "delta_next", "delta_best",
        "matched_peaks", "longest_b", "longest_y", "longest_y_pct", "missed_cleavages",
        "matched_intensity_pct", "poisson", "charge",
        "intensity_ms1", "intensity_ms2", "collision_energy"
    ]
    ds = ds.copy()

    ds["intensity_ms1"] = ds["intensity_ms1"].apply(lambda x: np.log(x + 1))
    ds["intensity_ms2"] = ds["intensity_ms2"].apply(lambda x: np.log(x + 1))

    X = ds[features].to_numpy().astype(np.float32)
    Y = ds["decoy"].to_numpy()
    Y = np.array([0 if x else 1 for x in Y]).astype(np.float32)

    return X, Y


def apply_mz_calibration(psm, fragments: pd.DataFrame, use_median: bool = True,
                         target_q: float = 0.01) -> float:
    """Apply a mass calibration to a set of peptide spectrum matches

    Args:
        psm: The peptide spectrum matches
        fragments: The fragments
        use_median: Whether to use the median ppm error
        target_q: The target q-value

    Returns:
        float: The ppm error
    """

    psms = []
    for _, item in psm.items():
        psms.extend(item)

    P = peptide_spectrum_match_list_to_pandas(psms)
    TDC = target_decoy_competition_pandas(P)
    TDC = TDC[TDC.q_value <= target_q]

    B = pd.merge(P, TDC, on=["spec_idx", "match_idx"])

    ppm_error = 0.0

    if use_median:
        ppm_error = np.median(B.median_ppm)
    else:
        ppm_error = np.mean(B.mean_ppm)

    fragments.apply(lambda row: row.processed_spec.calibrate_mz_ppm(ppm_error), axis=1)

    return ppm_error
