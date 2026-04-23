from typing import Tuple
from numpy.typing import NDArray
import numpy as np
import sagepy_connector

psc = sagepy_connector.py_utility


def calculate_pep_single(
        scores: NDArray,
        decoys: NDArray,
        bins: int = 1000,
        bw_adjust: float = 1.0,
        monotonic: bool = True
) -> Tuple[NDArray, float, float]:
    scores = np.asarray(scores, dtype=np.float64)
    decoys = np.asarray(decoys, dtype=bool)
    pep_bins, min_score, score_step = psc.calculate_pep_single(
        scores.tolist(),
        decoys.tolist(),
        bins,
        bw_adjust,
        monotonic,
    )
    return np.asarray(pep_bins, dtype=np.float64), min_score, score_step


def posterior_error(pep_bins: NDArray,
                    min_score: float,
                    score_step: float,
                    score: float) -> float:
    pep_bins = np.asarray(pep_bins, dtype=np.float64)
    return psc.posterior_error(pep_bins.tolist(), min_score, score_step, score)


def calculate_pep(scores: NDArray,
                  decoys: NDArray,
                  bins: int = 1000,
                  bw_adjust: float = 1.0,
                  monotonic: bool = True) -> NDArray:
    scores = np.asarray(scores, dtype=np.float64)
    decoys = np.asarray(decoys, dtype=bool)
    pep = psc.calculate_pep(
        scores.tolist(),
        decoys.tolist(),
        bins,
        bw_adjust,
        monotonic,
    )
    return np.asarray(pep, dtype=np.float64)

if __name__ == "__main__":
    # create 1000 radom scores between 0 and 50
    scores = np.random.uniform(0, 50, 50000)
    decoys = np.random.choice([True, False], 50000)

    # sort scores ascending
    scores = np.sort(scores)

    # sort decoys where true comes first
    decoys = np.sort(decoys)[::-1]

    pep_bins, min_score, score_step = calculate_pep_single(scores, decoys)
    pep = posterior_error(pep_bins, min_score, score_step, scores[0])

    peps = calculate_pep(scores, decoys)

    from matplotlib import pyplot as plt
    plt.plot(scores, peps)
    plt.show()
