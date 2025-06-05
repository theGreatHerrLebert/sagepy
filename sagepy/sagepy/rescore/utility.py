import numpy as np
import random
import pandas as pd

from numpy.typing import NDArray
from typing import Optional, Tuple, List

from sagepy.core.scoring import Psm
from sagepy.utility import psm_collection_to_pandas

from collections import defaultdict
from typing import Callable

peptide_key_maker = lambda psm: (psm.sequence,)
ion_key_maker = lambda psm: (psm.sequence, psm.charge)


def assign_random_groups(group_count, number_of_assignments, seed=None):
    """
    Assign random groups to a number of assignments.
    Args:
        group_count: int, number of groups to assign to
        number_of_assignments: int, number of assignments to make
        seed: None | int, random seed for reproducibility

    Returns:
        list: a list of group assignments, where each assignment is an integer in the range [0, group_count - 1]
    """
    if seed is not None:
        random.seed(seed)
    return [random.randint(0, group_count - 1) for _ in range(number_of_assignments)]


def split_into_chunks(
        psms: list,
        splits_count: int = 3,
        key_maker: Callable = peptide_key_maker,
        seed: None | int = None,
) -> list[list]:
    """
    Split a list of PSMs into multiple chunks based on a key maker function.
    Args:
        psms: list of Psm objects
        splits_count: int, number of splits to create
        key_maker: Callable, a function that takes a Psm object and returns a key for grouping
        seed: None | int, random seed for reproducibility

    Returns:
        list[list]: a list of lists, where each inner list contains Psm objects that share the same key
    """
    grouped_psms = defaultdict(list)
    for psm in psms:
        grouped_psms[key_maker(psm)].append(psm)

    split_assignments = assign_random_groups(splits_count, len(grouped_psms), seed=seed)

    splits = [[] for _ in range(splits_count)]
    for split_assignment, (group, grouped_psms) in zip(split_assignments, grouped_psms.items()):
        splits[split_assignment].extend(grouped_psms)

    return splits


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

    score = score if score is not None else "hyperscore"

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
        "collision_energy",
        "cosine_similarity",
        "spectral_angle_similarity",
        "pearson_correlation",
        "spearman_correlation",
        "spectral_entropy_similarity",
    ]
    ds = ds.copy()

    # Log-transform the intensity columns
    ds["intensity_ms1"] = ds["intensity_ms1"].apply(lambda x: np.log1p(x))
    ds["intensity_ms2"] = ds["intensity_ms2"].apply(lambda x: np.log1p(x))

    # avoid none values for cosine similarity
    ds["cosine_similarity"] = ds["cosine_similarity"].apply(lambda x: 0.0 if x is None else x)

    X = ds[features].to_numpy().astype(np.float32)

    if replace_nan:
        X = np.nan_to_num(X)

    Y = ds["decoy"].to_numpy()
    Y = np.array([0 if x else 1 for x in Y]).astype(np.float32)

    return X, Y


def generate_training_data(
        psm_list: List[Psm],
        method: str = "peptide_q",
        q_max: float = 0.01,
        balance: bool = True,
        replace_nan: bool = True,
        num_threads: int = 16,
        **kwargs
) -> Tuple[NDArray, NDArray]:
    """ Generate training data.
    Args:
        psm_list: List of PeptideSpectrumMatch objects
        method: Method to use for training data generation
        q_max: Maximum q-value allowed for positive examples
        balance: Whether to balance the dataset
        replace_nan: Whether to replace NaN values with 0
        num_threads: Number of threads to use for feature extraction

    Returns:
        Tuple[NDArray, NDArray]: X_train and Y_train
    """
    # create pandas table from psms
    PSM_pandas = psm_collection_to_pandas(psm_list, num_threads=num_threads)

    if method == "spectrum_q":
        TARGET = PSM_pandas[(PSM_pandas.decoy == False) & (PSM_pandas.spectrum_q <= q_max) & (PSM_pandas["rank"] == 1)]
    elif method == "peptide_q":
        TARGET = PSM_pandas[(PSM_pandas.decoy == False) & (PSM_pandas.peptide_q <= q_max) & (PSM_pandas["rank"] == 1)]

    elif method == "decoy_quantile":
        cutoff = PSM_pandas[PSM_pandas.decoy].hyperscore.quantile(1 - q_max)
        TARGET = PSM_pandas[(PSM_pandas.decoy == False) & (PSM_pandas["rank"] == 1) & (PSM_pandas.hyperscore >= cutoff)]
    else:
        raise ValueError(f"Unknown method: {method}. Use 'spectrum_q' or 'peptide_q'.")

    X_target, Y_target = get_features(TARGET, replace_nan=replace_nan)

    # select all decoys
    DECOY = PSM_pandas[PSM_pandas.decoy & (PSM_pandas["rank"] == 1)]
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


def get_list_index_by_sequence(psms, num_splits: int = 5, seed: int = 35):
    # Set random seed for reproducibility
    np.random.seed(seed)

    # Extract unique sequences
    unique_sequences = list({psm.sequence for psm in psms})

    # Shuffle and split
    shuffled = np.random.permutation(unique_sequences)
    split = np.array_split(shuffled, num_splits)

    # Create mapping from sequence -> split index
    index_dict = {seq: i for i, group in enumerate(split) for seq in group}
    return index_dict


def split_psm_list_broken(psm_list: List[Psm], num_splits: int = 5) -> List[List]:
    # Get sequence-to-split mapping
    seq_to_split = get_list_index_by_sequence(psm_list, num_splits)

    # Preallocate split containers
    splits = [[] for _ in range(num_splits)]

    # Assign PSMs to their respective splits
    for psm in psm_list:
        split_idx = seq_to_split[psm.sequence]
        splits[split_idx].append(psm)

    return splits

def split_psm_list(psm_list: List[Psm], num_splits: int = 5, seed: None| int = None, key_maker: Callable = peptide_key_maker, **kwargs) -> List[List[Psm]]:
    """
    Split PSMs into multiple splits.

    Args:
        psm_list: List of PeptideSpectrumMatch objects
        num_splits: Number of splits
        seed: Optional seed for reproducibility
        key_maker: Callable function to create keys for grouping PSMs

    Returns:
        List[List[PeptideSpectrumMatch]]: List of splits

    """
    return split_into_chunks(psm_list, num_splits, seed=seed, key_maker=key_maker)

"""
def split_psm_list(psm_list: List[Psm], num_splits: int = 5) -> List[List[Psm]]:
     Split PSMs into multiple splits.

    Args:
        psm_list: List of PeptideSpectrumMatch objects
        num_splits: Number of splits

    Returns:
        List[List[PeptideSpectrumMatch]]: List of splits

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
"""

def transform_psm_to_mokapot_pin(psm_df, seq_modified: bool = False):
    """ Transform a PSM DataFrame to a mokapot PIN DataFrame.
    Args:
        psm_df: a DataFrame containing PSMs
        seq_modified: whether the sequences are modified

    Returns:
        A DataFrame containing the PSMs in mokapot PIN format.
    """

    columns_map = {
        # target columns mapping for mokapot
        'spec_idx': 'SpecId',
        'decoy': 'Label',
        'charge': 'Charge',
        'sequence_modified': 'Peptide',
        'proteins': 'Proteins',

        # feature mapping for re-scoring
        'hyperscore': 'Feature1',
        'isotope_error': 'Feature2',
        'delta_mass': 'Feature3',
        'delta_rt': 'Feature4',
        'delta_ims': 'Feature5',
        'matched_peaks': 'Feature6',
        'matched_intensity_pct': 'Feature7',
        'intensity_ms1': 'Feature8',
        'intensity_ms2': 'Feature9',
        'average_ppm': 'Feature10',
        'poisson': 'Feature11',
        'spectral_entropy_similarity': 'Feature12',
        'pearson_correlation': 'Feature13',
        'spearman_correlation': 'Feature14',
        'spectral_angle_similarity': 'Feature15',
        'collision_energy': 'Feature16',
        'delta_next': 'Feature17',
        'delta_best': 'Feature18',
        'longest_b': 'Feature19',
        'longest_y': 'Feature20',
        'longest_y_pct': 'Feature21',
        'cosine_similarity': 'Feature22',
        'rank': 'Feature23',
        'missed_cleavages': 'Feature24',
    }

    if not seq_modified:
        columns_map['sequence'] = 'Peptide'
        columns_map.pop('sequence_modified')

    psm_df = psm_df[list(columns_map.keys())]
    df_pin = psm_df.rename(columns=columns_map)
    df_pin_clean = df_pin.dropna(axis=1, how='all')
    df_pin_clean = df_pin_clean.dropna()

    df_pin_clean['Label'] = df_pin_clean['Label'].apply(lambda x: -1 if x else 1)
    df_pin_clean['ScanNr'] = range(1, len(df_pin_clean) + 1)

    return df_pin_clean
