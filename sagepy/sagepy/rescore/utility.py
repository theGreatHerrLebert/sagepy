import numpy as np
import pandas as pd
from numpy.typing import NDArray
from typing import Optional, Tuple, List

from sagepy.core import PeptideSpectrumMatch
from sagepy.qfdr.tdc import target_decoy_competition_pandas
from sagepy.utility import peptide_spectrum_match_collection_to_pandas


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
