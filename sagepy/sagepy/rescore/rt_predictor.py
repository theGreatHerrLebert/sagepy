import re
from typing import List

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import Ridge, Lasso

from sagepy.core import Psm
from sagepy.qfdr.tdc import target_decoy_competition_pandas
from sagepy.utility import psm_collection_to_pandas


def tokenize_peptide(sequence: str) -> List[str]:
    """
    Tokenize a peptide sequence into amino acid tokens.
    Args:
        sequence: A peptide sequence string

    Returns:
        A list of amino acid tokens, including modifications
    """
    tokens = re.findall(r'[A-Z](?:\[UNIMOD:\d+\])?', sequence)
    return tokens


def sequence_to_vector(sequence: str, token_alphabet: List[str]) -> np.ndarray:
    """
    Convert a peptide sequence into a vector representation based on a token alphabet
    Args:
        sequence: A peptide sequence string
        token_alphabet: A list of amino acid tokens

    Returns:
        A vector representation of the peptide sequence
    """
    tokens = tokenize_peptide(sequence)
    vector = np.zeros(len(token_alphabet))
    for token in tokens:
        if token in token_alphabet:
            idx = token_alphabet.index(token)
            vector[idx] += 1
    return vector


def create_token_alphabet(sequences):
    """
    Create a token alphabet from a set of peptide sequences.
    Args:
        sequences: A list of peptide sequences

    Returns:
        A list of unique amino acid tokens
    """
    unique_tokens = set()
    for seq in sequences:
        unique_tokens.update(tokenize_peptide(seq))
    token_alphabet = sorted(unique_tokens)  # Sort to maintain consistent ordering
    return token_alphabet


def prepare_data(sequences, retention_times, token_alphabet):
    """
    Prepare the dataset for training a retention time predictor
    Args:
        sequences: A list of peptide sequences
        retention_times: A list of retention times
        token_alphabet: A list of amino acid tokens

    Returns:
        X: A matrix of input features
        y: A vector of target values
    """
    # Convert sequences to vectors based on the provided token_alphabet
    X = np.array([sequence_to_vector(seq, token_alphabet) for seq in sequences])
    y = np.array(retention_times)

    return X, y


def train_ridge_regression_model(X, y, alpha=1.0, verbose=False):
    """
    Train a ridge regression model with L2 regularization.
    Args:
        X: tokenized peptide sequences
        y: retention times
        alpha: regularization strength
        verbose: whether to print the test MSE

    Returns:
        Trained model and train/test data
    """

    model = Ridge(alpha=alpha)
    model.fit(X, y)

    if verbose:
        y_pred = model.predict(X)
        mse = mean_squared_error(y, y_pred)
        print(f"Test MSE (Ridge, alpha={alpha}): {mse}")

    return model


def train_lasso_regression_model(X, y, alpha=1.0, verbose=False):
    """
    Train a Lasso regression model with L1 regularization.
    Args:
        X: Peptide sequences
        y: Retention times
        alpha: Regularization strength
        verbose: Whether to print the test MSE

    Returns:
        Trained model and train/test data
    """

    model = Lasso(alpha=alpha)
    model.fit(X, y)

    if verbose:
        y_pred = model.predict(X)
        mse = mean_squared_error(y, y_pred)
        print(f"Test MSE (Lasso, alpha={alpha}): {mse}")

    return model


def transform_sequences(sequences, token_alphabet):
    """
    Transform a list of peptide sequences into a matrix of token counts based on a token alphabet.
    Args:
        sequences: A list of peptide sequences
        token_alphabet: A list of amino acid tokens

    Returns:
        A matrix of token counts
    """
    X_new = np.array([sequence_to_vector(seq, token_alphabet) for seq in sequences])
    return X_new

def predict_retention_times_psm(psm_collection: List[Psm], fdr_threshold: float = 0.01, alpha: float = 0.2):
    """
    Predict retention times for peptide spectrum matches using a ridge regression model
    Args:
        psm_collection: A list of PeptideSpectrumMatch objects
        fdr_threshold: The false discovery rate threshold for selecting target hits
        alpha: The regularization strength for the ridge regression model

    Returns:
        None, the retention times are updated in place in the PeptideSpectrumMatch objects
    """

    # Convert the peptide spectrum matches to a pandas DataFrame
    PSM_pandas = psm_collection_to_pandas(psm_collection)

    # Create a token alphabet from the peptide sequences
    token_alphabet = create_token_alphabet(PSM_pandas["sequence"])

    # Prepare the dataset for training
    TDC = target_decoy_competition_pandas(PSM_pandas, method="psm")
    TDC = pd.merge(TDC.drop(columns=["score"]),
                                   PSM_pandas.drop(columns=["q_value"]), on=["spec_idx", "match_idx", "decoy"])

    # we can select target hits with a q-value of 0.01, translating to 1 percent FDR
    FDR_controlled = TDC[(TDC.q_value <= fdr_threshold) & (TDC.decoy == False)]
    X, y = prepare_data(FDR_controlled["sequence"], FDR_controlled["retention_time_observed"], token_alphabet)

    # Train a ridge regression model
    model = train_ridge_regression_model(X, y, alpha=alpha)

    # Predict retention times for all PSMs
    X_new = transform_sequences(PSM_pandas["sequence"], token_alphabet)
    predicted_times = model.predict(X_new)

    for rt_pred, psm in zip(predicted_times, psm_collection):
        psm.retention_time_predicted = rt_pred
