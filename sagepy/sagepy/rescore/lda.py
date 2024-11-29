import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler

from tqdm import tqdm
from typing import Union, List, Dict

from sagepy.core import Psm
from sagepy.rescore.utility import get_features, generate_training_data, split_psm_list
from sagepy.utility import psm_collection_to_pandas


def rescore_lda(
        psm_collection: Union[List[Psm], Dict[str, List[Psm]]],
        num_splits: int = 5,
        verbose: bool = True,
        balance: bool = True,
        replace_nan: bool = True,
        score: str = "hyperscore",
        num_threads: int = 16,
) -> List[Psm]:
    """ Re-score PSMs using Linear Discriminant Analysis (LDA).
    Args:
        psm_collection: A collection of PSMs
        num_splits: Number of splits (folds) to use for cross-validation
        verbose: Whether to print progress
        balance: Whether to balance the dataset (equal number of target and decoy examples)
        replace_nan: Whether to replace NaN values with 0
        score: Score to use for rescoring
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


    X_all, _ = get_features(psm_collection_to_pandas(psm_list, num_threads=num_threads), score=score, replace_nan=replace_nan)
    scaler = StandardScaler()
    scaler.fit(X_all)

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
        X, _ = get_features(psm_collection_to_pandas(target, num_threads=num_threads), replace_nan=replace_nan)

        # experimenting with different settings for LDA showed that shrinkage should be used, which tries to
        # keep model weights small and helps to prevent overfitting
        lda = LinearDiscriminantAnalysis(solver="eigen", shrinkage="auto")
        lda.fit(scaler.transform(X_train), Y_train)

        try:
            # check for flip sign of LDA classification return to be compatible with good score ascending
            score_flip = 1.0 if Y_train[np.argmax(np.squeeze(lda.transform(scaler.transform(X_train))))] == 1.0 else -1.0
        except:
            score_flip = 1.0

        Y_pred = np.squeeze(lda.transform(scaler.transform(X))) * score_flip
        predictions.extend(Y_pred)

    for score, match in zip(predictions, psm_list):
        match.re_score = score

    return psm_list
