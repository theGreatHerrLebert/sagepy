import numpy as np
import pandas as pd

import scipy.stats

from numpy.typing import NDArray
from typing import Optional, Tuple, List

from sagepy.core import PeptideSpectrumMatch
from sagepy.qfdr.tdc import target_decoy_competition_pandas
from sagepy.utility import peptide_spectrum_match_collection_to_pandas


def dict_to_dense_array(peak_dict, array_length=174):
    """
    Convert a dictionary of peaks to a fixed-length array.
    Args:
        peak_dict: a dictionary of peaks (ion_type, charge, ordinal) -> intensity
        array_length: the length of the fixed-length array

    Returns:
        A fixed-length array of intensities.
    """
    # initialize a fixed-length array of zeros
    intensities = np.zeros(array_length)

    half_length = array_length // 2  # first half for b ions, second half for y ions
    block_size = 29  # number of ordinals per charge state

    for (ion_type, charge, ordinal), intensity in peak_dict.items():
        # map (b=0, y=1) ions to the correct index
        index = ion_type * half_length + (charge - 1) * block_size + (ordinal - 1)
        intensities[index] = intensity

    return intensities


def spectral_entropy_similarity(observed_intensities, predicted_intensities,
                                epsilon: float = 1e-7) -> float:
    """
    Calculate the spectral entropy similarity between observed and predicted intensities
    Args:
        observed_intensities: the observed intensities
        predicted_intensities: the predicted intensities
        epsilon: small value to avoid division by zero

    Returns:
        The spectral entropy similarity between observed and predicted intensities.
    """

    valid_ion_mask = predicted_intensities > epsilon

    observed_filtered = observed_intensities[valid_ion_mask]
    predicted_filtered = predicted_intensities[valid_ion_mask]

    entropy_merged = scipy.stats.entropy(observed_filtered + predicted_filtered)
    entropy_obs = scipy.stats.entropy(observed_filtered)
    entropy_pred = scipy.stats.entropy(predicted_filtered)

    # calculate the spectral entropy similarity
    entropy = 1 - (2 * entropy_merged - entropy_obs - entropy_pred) / np.log(4)

    # handle cases where the entropy is NaN (set them to 0)
    if np.isnan(entropy):
        entropy = 0

    return entropy


def spectral_correlation(
        observed_intensities,
        predicted_intensities,
        method: str = "pearson", epsilon: float = 1e-7,
) -> float:
    """
    Calculate the spectral correlation between observed and predicted intensities.
    Args:
        observed_intensities: intensities observed in the spectrum
        predicted_intensities: intensities predicted by the model
        method: correlation method (pearson or spearman)
        epsilon: small value to avoid division by zero

    Returns:
        The spectral correlation between observed and predicted intensities.
    """

    if method not in ["pearson", "spearman"]:
        raise ValueError(f"Invalid correlation method: {method}. Choose 'pearson' or 'spearman'.")

    valid_ion_mask = predicted_intensities > epsilon

    observed_filtered = observed_intensities[valid_ion_mask]
    predicted_filtered = predicted_intensities[valid_ion_mask]

    observed_filtered = observed_filtered[~np.isnan(observed_filtered)]
    predicted_filtered = predicted_filtered[~np.isnan(predicted_filtered)]

    if len(observed_filtered) <= 2 or len(predicted_filtered) <= 2:
        return 0

    if method == "pearson":
        corr, _ = scipy.stats.pearsonr(observed_filtered, predicted_filtered)
    else:
        corr, _ = scipy.stats.spearmanr(observed_filtered, predicted_filtered)

    if np.isnan(corr):
        corr = 0

    return corr


def spectral_coverage(observed_intensities, predicted_intensities):
    """
    Calculate the spectral coverage between observed and predicted intensities
    Args:
        observed_intensities: the observed intensities
        predicted_intensities: the predicted intensities

    Returns:
        The spectral coverage between observed and predicted intensities.
    """
    predicted = predicted_intensities
    observed = observed_intensities

    intensity_covered = 0.0
    total_intensity = 0.0

    for key, value in predicted.items():
        total_intensity += value

        if key in observed:
            intensity_covered += value

    return intensity_covered / total_intensity


def cosim_to_spectral_angle_sim(cosim: float) -> float:
    return 1 - ((np.arccos(cosim) * (180 / np.pi)) / 180)


def get_features(
        ds: pd.DataFrame,
        score: Optional[str] = None,
        replace_nan: bool = True,
) -> Tuple[NDArray, NDArray]:
    """
    Get features and labels from a dataset.
    Args:
        ds: a pandas DataFrame containing the dataset.
        score: the name of the target score column.
        replace_nan: if True, replace NaN values with 0.

    Returns:
        A tuple containing the features and labels.
    """

    score = score if score is not None else "score"

    # The currently used features for the model fit
    # TODO: extend this list with additional features
    features = [
        f"{score}",
        "delta_rt",
        "delta_ims",
        "cosine_similarity",
        "delta_mass",
        "rank",
        "isotope_error",
        "average_ppm",
        "delta_next",
        "delta_best",
        "matched_peaks",
        "longest_b",
        "longest_y",
        "longest_y_pct",
        "missed_cleavages",
        "matched_intensity_pct",
        "poisson",
        "charge",
        "intensity_ms1",
        "intensity_ms2",
        "collision_energy"
    ]
    ds = ds.copy()

    # Log-transform the intensity columns
    ds["intensity_ms1"] = ds["intensity_ms1"].apply(lambda x: np.log1p(x))
    ds["intensity_ms2"] = ds["intensity_ms2"].apply(lambda x: np.log1p(x))

    # avoid none values for cosine similarity
    ds["cosine_similarity"] = ds["cosine_similarity"].apply(lambda x: 0.0 if x is None else x)

    # transform cosine similarity to spectral angle similarity
    ds["cosine_similarity"] = ds["cosine_similarity"].apply(cosim_to_spectral_angle_sim)

    X = ds[features].to_numpy().astype(np.float32)

    if replace_nan:
        X = np.nan_to_num(X)

    Y = ds["decoy"].to_numpy()
    Y = np.array([0 if x else 1 for x in Y]).astype(np.float32)

    return X, Y


def generate_training_data(
        psm_list: List[PeptideSpectrumMatch],
        method: str = "psm",
        q_max: float = 0.01,
        balance: bool = True,
        replace_nan: bool = True,
) -> Tuple[NDArray, NDArray]:
    """ Generate training data.
    Args:
        psm_list: List of PeptideSpectrumMatch objects
        method: Method to use for training data generation
        q_max: Maximum q-value allowed for positive examples
        balance: Whether to balance the dataset
        replace_nan: Whether to replace NaN values with 0

    Returns:
        Tuple[NDArray, NDArray]: X_train and Y_train
    """
    # create pandas table from psms
    PSM_pandas = peptide_spectrum_match_collection_to_pandas(psm_list)

    # calculate q-values to get inital "good" hits
    PSM_q = target_decoy_competition_pandas(PSM_pandas, method=method)
    PSM_pandas = PSM_pandas.drop(columns=["q_value", "score"])

    # merge data with q-values
    TDC = pd.merge(
        PSM_q, PSM_pandas,
        left_on=["spec_idx", "match_idx", "decoy"],
        right_on=["spec_idx", "match_idx", "decoy"]
    )

    # select best positive examples
    TARGET = TDC[(TDC.decoy == False) & (TDC.q_value <= q_max)]
    X_target, Y_target = get_features(TARGET, replace_nan=replace_nan)

    # select all decoys
    DECOY = TDC[TDC.decoy]
    X_decoy, Y_decoy = get_features(DECOY, replace_nan=replace_nan)

    # balance the dataset such that the number of target and decoy examples are equal
    if balance:
        num_target = np.min((len(DECOY), len(TARGET)))
        target_indices = np.random.choice(np.arange(len(X_target)), size=num_target)
        X_target = X_target[target_indices, :]
        Y_target = Y_target[target_indices]

    # combine target and decoy examples
    X_train = np.vstack((X_target, X_decoy))
    Y_train = np.hstack((Y_target, Y_decoy))

    return X_train, Y_train


def split_psm_list(psm_list: List[PeptideSpectrumMatch], num_splits: int = 5) -> List[List[PeptideSpectrumMatch]]:
    """ Split PSMs into multiple splits.

    Args:
        psm_list: List of PeptideSpectrumMatch objects
        num_splits: Number of splits

    Returns:
        List[List[PeptideSpectrumMatch]]: List of splits
    """

    # floor division by num_splits
    split_size = len(psm_list) // num_splits

    # remainder for last split
    remainder = len(psm_list) % num_splits

    splits = []

    start_index = 0

    for i in range(num_splits):
        end_index = start_index + split_size + (1 if i < remainder else 0)
        splits.append(psm_list[start_index:end_index])
        start_index = end_index

    return splits
