# sagepy/rescore/rescore.py

from __future__ import annotations

from collections import defaultdict
from typing import Union, List, Dict, Callable, Optional

import numpy as np
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, MinMaxScaler

from sagepy.core import Psm
from sagepy.rescore.utility import get_features, generate_training_data, split_psm_list
from sagepy.utility import psm_collection_to_pandas


# ----------------------------
# Group-stratified splitter
# ----------------------------

def _default_key_maker(psm: Psm):
    # peptide-level grouping (prevents leakage across same peptide)
    return (psm.sequence,)


def _split_psm_list_stratified_groups(
    psm_list: List[Psm],
    *,
    num_splits: int,
    seed: int = 35,
    key_maker: Callable[[Psm], object] = _default_key_maker,
    pos_ok: Optional[Callable[[Psm], bool]] = None,
) -> List[List[Psm]]:
    """
    Split PSMs into folds by grouping (key_maker), but *stratify by group label*
    (whether the group contains at least one "positive-capable" PSM).

    This reduces the risk that some training folds have zero positives after filtering.
    """
    rng = np.random.default_rng(seed)

    if pos_ok is None:
        # Safe default: rank-1 targets count as "positive-capable"
        pos_ok = lambda p: (not getattr(p, "decoy", False)) and (getattr(p, "rank", 1) == 1)

    # Group PSMs
    groups: Dict[object, List[Psm]] = defaultdict(list)
    for psm in psm_list:
        groups[key_maker(psm)].append(psm)

    tagged = []
    for gk, psms in groups.items():
        has_pos = any(pos_ok(p) for p in psms)
        tagged.append((gk, psms, has_pos))

    # Shuffle then assign: place pos-groups first so they get spread out
    rng.shuffle(tagged)
    tagged.sort(key=lambda x: (not x[2], -len(x[1])))  # pos first, then larger groups

    splits: List[List[Psm]] = [[] for _ in range(num_splits)]
    pos_counts = np.zeros(num_splits, dtype=int)
    size_counts = np.zeros(num_splits, dtype=int)

    for _, psms, has_pos in tagged:
        if has_pos:
            # pick fold with lowest pos count; tie-break by current size
            candidates = np.where(pos_counts == pos_counts.min())[0]
            fold = int(candidates[np.argmin(size_counts[candidates])])
            pos_counts[fold] += 1
        else:
            # keep sizes balanced for non-pos groups
            fold = int(np.argmin(size_counts))

        splits[fold].extend(psms)
        size_counts[fold] += len(psms)

    return splits


# ----------------------------
# Main rescoring function
# ----------------------------

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

    # Training-positive selection (defaults to the quantile trick)
    train_method: str = "decoy_quantile",
    q_max: float = 0.01,

    # NEW: reduce single-class folds by stratifying peptide-groups across folds
    stratify_groups: bool = True,
    stratify_seed: int = 35,
    key_maker: Callable[[Psm], object] = _default_key_maker,

    **kwargs,
) -> List[Psm]:
    """
    Re-score PSMs using cross-validated training.

    Defaults:
      - train_method="decoy_quantile" (robust positives)
      - stratify_groups=True (group-stratified fold assignment by peptide)

    Returns:
        List[Psm] with match.re_score set.
    """

    # Flatten collection
    if isinstance(psm_collection, dict):
        psm_list: List[Psm] = []
        for _, psm_candidates in psm_collection.items():
            psm_list.extend(psm_candidates)
    else:
        psm_list = psm_collection

    if num_splits < 2:
        raise ValueError("num_splits must be >= 2")

    # Fit scaler on all PSMs (unsupervised)
    all_df = psm_collection_to_pandas(psm_list, num_threads=num_threads)
    X_all, _ = get_features(all_df, score=score, replace_nan=replace_nan)

    scaler = MinMaxScaler() if use_min_max_scaler else StandardScaler()
    scaler.fit(X_all)

    # Build a fold splitter
    if stratify_groups:
        # Define "positive-capable" at the group level.
        # For decoy_quantile, we just need rank-1 targets to exist at all.
        # (The quantile cutoff is computed per training fold anyway.)
        def pos_ok(p: Psm) -> bool:
            return (not getattr(p, "decoy", False)) and (getattr(p, "rank", 1) == 1)

        splits = _split_psm_list_stratified_groups(
            psm_list,
            num_splits=num_splits,
            seed=stratify_seed,
            key_maker=key_maker,
            pos_ok=pos_ok,
        )
    else:
        # Fall back to original random group assignment splitter
        splits = split_psm_list(psm_list=psm_list, num_splits=num_splits, key_maker=key_maker, **kwargs)

    predictions: List[float] = []
    final_psms: List[Psm] = []

    for i in tqdm(range(num_splits), disable=not verbose, desc="Re-scoring PSMs", ncols=100):
        target_fold = splits[i]
        final_psms.extend(target_fold)

        train_psms: List[Psm] = []
        for j in range(num_splits):
            if j != i:
                train_psms.extend(splits[j])

        # Generate training set (quantile trick by default)
        X_train, Y_train = generate_training_data(
            train_psms,
            method=train_method,
            q_max=q_max,
            balance=balance,
            replace_nan=replace_nan,
            num_threads=num_threads,
            **kwargs,
        )

        # Guard against single-class training labels (XGBoost logistic will crash otherwise)
        uniq = np.unique(Y_train)
        if uniq.size < 2:
            raise ValueError(
                f"Fold {i}: training labels have only one class ({uniq}). "
                f"Try num_splits=2/3, or confirm you have both targets+decoys at rank==1, "
                f"or loosen filters / change key_maker."
            )

        # Features for held-out fold (IMPORTANT: pass score=score)
        fold_df = psm_collection_to_pandas(target_fold, num_threads=num_threads)
        X_fold, _ = get_features(fold_df, score=score, replace_nan=replace_nan)

        model.fit(scaler.transform(X_train), Y_train)

        # decision_function preferred; otherwise predict_proba
        if hasattr(model, "decision_function"):
            y_pred = model.decision_function(scaler.transform(X_fold))
            predictions.extend([float(v) for v in y_pred])
        else:
            y_pred = model.predict_proba(scaler.transform(X_fold))
            predictions.extend([float(v) for v in y_pred[:, 1]])

    # Assign rescored values
    for s, match in zip(predictions, final_psms):
        match.re_score = s

    return final_psms