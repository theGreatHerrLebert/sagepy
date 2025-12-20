from sklearn.preprocessing import StandardScaler, MinMaxScaler
from tqdm import tqdm
from typing import Union, List, Dict

import numpy as np

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

        # NEW: quantile training defaults
        train_method: str = "decoy_quantile",
        q_max: float = 0.01,

        **kwargs,
) -> List[Psm]:
    """
    Re-score PSMs using cross-validated training.

    This version defaults to the "decoy_quantile" trick for selecting positive
    training examples (targets with score >= decoy quantile cutoff), which is
    much less likely to produce empty TARGET sets in some folds than peptide_q/spectrum_q.

    Args:
        psm_collection: list or dict of PSMs
        model: sklearn-compatible model (fit/predict_proba or decision_function)
        use_min_max_scaler: use MinMaxScaler instead of StandardScaler
        num_splits: number of CV folds
        verbose: show progress bar
        balance: balance TARGET/DECOY counts during training
        replace_nan: replace NaN with 0
        score: which score column is used as feature (Feature1)
        num_threads: threads used for psm_collection_to_pandas
        train_method: training positive selection method (default "decoy_quantile")
        q_max: for decoy_quantile: tail size (e.g. 0.01 keeps top 1% decoys as cutoff)

    Returns:
        List[Psm] with match.re_score filled
    """

    # Flatten collection
    if isinstance(psm_collection, dict):
        psm_list: List[Psm] = []
        for _, psm_candidates in psm_collection.items():
            psm_list.extend(psm_candidates)
    else:
        psm_list = psm_collection

    # Fit scaler on *all* PSMs (unsupervised step)
    all_df = psm_collection_to_pandas(psm_list, num_threads=num_threads)
    X_all, _ = get_features(all_df, score=score, replace_nan=replace_nan)

    scaler = MinMaxScaler() if use_min_max_scaler else StandardScaler()
    scaler.fit(X_all)

    # Split for CV (usually grouped by peptide sequence via your split_psm_list)
    splits = split_psm_list(psm_list=psm_list, num_splits=num_splits, **kwargs)

    predictions: List[float] = []
    final_psms: List[Psm] = []

    for i in tqdm(range(num_splits), disable=not verbose, desc='Re-scoring PSMs', ncols=100):
        target_fold = splits[i]
        final_psms.extend(target_fold)

        train_psms: List[Psm] = []
        for j in range(num_splits):
            if j != i:
                train_psms.extend(splits[j])

        # Generate training set using decoy-quantile positives by default
        X_train, Y_train = generate_training_data(
            train_psms,
            method=train_method,
            q_max=q_max,
            balance=balance,
            replace_nan=replace_nan,
            num_threads=num_threads,
            **kwargs,
        )

        # Guard: XGBoost (and others) can’t train on a single class
        uniq = np.unique(Y_train)
        if uniq.size < 2:
            raise ValueError(
                f"Fold {i}: training labels have only one class ({uniq}). "
                f"Try fewer splits, different key_maker, or adjust q_max/train_method."
            )

        # Features for held-out fold (IMPORTANT: pass score=score)
        X_fold, _ = get_features(
            psm_collection_to_pandas(target_fold, num_threads=num_threads),
            score=score,
            replace_nan=replace_nan,
        )

        model.fit(scaler.transform(X_train), Y_train)

        # Use decision_function if available, else probability
        if hasattr(model, "decision_function"):
            Y_pred = model.decision_function(scaler.transform(X_fold))
            predictions.extend([float(v) for v in Y_pred])
        else:
            Y_pred = model.predict_proba(scaler.transform(X_fold))
            predictions.extend([float(v) for v in Y_pred[:, 1]])

    # Assign rescored values
    for s, match in zip(predictions, final_psms):
        match.re_score = s

    return final_psms