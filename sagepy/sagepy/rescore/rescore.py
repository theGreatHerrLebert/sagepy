from sklearn.preprocessing import StandardScaler, MinMaxScaler

from tqdm import tqdm
from typing import Union, List, Dict

from sagepy.core import Psm
from sagepy.rescore.utility import get_features, generate_training_data, split_psm_list
from sagepy.utility import psm_collection_to_pandas


def rescore_psms(
        psm_collection: Union[List[Psm], Dict[str, List[Psm]]],
        model,
        use_min_max_scaler: bool = False,
        num_splits: int = 3,
        verbose: bool = True,
        balance: bool = True,
        replace_nan: bool = True,
        score: str = "hyperscore",
        num_threads: int = 16,
) -> List[Psm]:
    """ Re-score PSMs using a model (e.g. Random Forest, Gradient Boosting, etc.).
    Args:
        psm_collection: A collection of PSMs
        model: A model to use for re-scoring, needs to comply to the sklearn API
        use_min_max_scaler: Whether to use MinMaxScaler instead of StandardScaler
        num_splits: Number of splits (folds) to use for cross-validation
        verbose: Whether to print progress
        balance: Whether to balance the dataset (equal number of target and decoy examples)
        replace_nan: Whether to replace NaN values with 0
        score: Score to use for re-scoring
        num_threads: Number of threads to use for feature extraction

    Returns:
        List[PeptideSpectrumMatch]: List of PeptideSpectrumMatch objects
    """

    psm_list = []

    if isinstance(psm_collection, dict):
        for spec_id, psm_candidates in psm_collection.items():
            psm_list.extend(psm_candidates)
    else:
        psm_list = psm_collection

    # get features for all PSMs, which will be a matrix of shape (n_samples, n_features)
    X_all, _ = get_features(psm_collection_to_pandas(psm_list, num_threads=num_threads), score=score, replace_nan=replace_nan)

    # use a scaler to scale the features
    if use_min_max_scaler:
        scaler = MinMaxScaler()
    else:
        scaler = StandardScaler()
    scaler.fit(X_all)

    # split the PSMs into num_splits folds to perform cross-validation
    splits = split_psm_list(psm_list=psm_list, num_splits=num_splits)

    predictions = []

    for i in tqdm(range(num_splits), disable=not verbose, desc='Re-scoring PSMs', ncols=100):

        target = splits[i]
        features = []

        for j in range(num_splits):
            if j != i:
                features.extend(splits[j])

        # generate training data
        X_train, Y_train = generate_training_data(features, balance=balance, replace_nan=replace_nan, num_threads=num_threads)

        # get features for target that we want to re-score
        X, _ = get_features(psm_collection_to_pandas(target), replace_nan=replace_nan)
        model.fit(scaler.transform(X_train), Y_train)

        # try to use decision function, otherwise use predict_proba
        try:
            Y_pred = model.decision_function(scaler.transform(X))
            predictions.extend(Y_pred)  # Use decision scores directly
        except AttributeError:
            Y_pred = model.predict_proba(scaler.transform(X))
            predictions.extend(Y_pred[:, 1])  # Use class probabilities (second column for binary classification)

    # assign the re-scored values to the PSMs
    for score, match in zip(predictions, psm_list):
        match.re_score = score

    return psm_list
