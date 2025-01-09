from typing import Tuple
from numpy.typing import NDArray
import numpy as np
from numba import njit, prange


@njit
def std(sample: NDArray) -> float:
    """Calculate the standard deviation of the sample.

    Args:
        sample: numpy array of values

    Returns:
        float: standard deviation of the sample
    """
    mean = np.mean(sample)
    variance = np.mean((sample - mean) ** 2)
    return np.sqrt(variance)


@njit(parallel=True)
def kde_pdf(sample: NDArray,
            bandwidth: float,
            x: float) -> float:
    """Calculate the KDE PDF for a given x.

    Args:
        sample: numpy array of values
        bandwidth: bandwidth parameter
        x: value for which to calculate the PDF

    Returns:
        float: KDE PDF for the given
    """
    n = len(sample)
    constant = (2.0 * np.pi) ** 0.5 * bandwidth * n
    sum_pdf = 0.0

    for i in prange(n):
        sum_pdf += np.exp(-0.5 * ((x - sample[i]) / bandwidth) ** 2)

    return sum_pdf / constant


@njit
def calculate_pep_single(
        scores: NDArray,
        decoys: NDArray,
        bins: int = 1000,
        bw_adjust: float = 1.0,
        monotonic: bool = True
) -> Tuple[NDArray, float, float]:
    """Calculate the PEP using KDE and binning with linear interpolation.

    Args:
        scores: numpy array of scores
        decoys: numpy array of boolean values indicating decoys
        bins: number of bins for the PEP calculation
        bw_adjust: bandwidth adjustment factor
        monotonic: whether to enforce monotonicity

    Returns:
        Tuple[NDArray, float, float]: PEP values, minimum score, score step
    """
    d = scores[decoys]
    t = scores[~decoys]

    pi = len(d) / len(scores)
    sigma_d = std(d)
    sigma_t = std(t)

    bandwidth_d = bw_adjust * sigma_d * (4. / 3. / len(d)) ** 0.2
    bandwidth_t = bw_adjust * sigma_t * (4. / 3. / len(t)) ** 0.2

    min_score = np.min(scores)
    max_score = np.max(scores)
    score_step = (max_score - min_score) / (bins - 1)

    pep_bins = np.zeros(bins)

    for bin_idx in range(bins):
        score = min_score + bin_idx * score_step
        decoy_pdf = kde_pdf(d, bandwidth_d, score) * pi
        target_pdf = kde_pdf(t, bandwidth_t, score) * (1.0 - pi)
        pep_bins[bin_idx] = decoy_pdf / (decoy_pdf + target_pdf)

    if monotonic:
        for i in range(len(pep_bins) - 2, -1, -1):
            pep_bins[i] = max(pep_bins[i], pep_bins[i + 1])

    return pep_bins, min_score, score_step


@njit
def posterior_error(pep_bins: NDArray,
                    min_score: float,
                    score_step: float,
                    score: float) -> float:
    """Interpolate PEP for a given score.

    Args:
        pep_bins: numpy array of PEP values
        min_score: minimum score
        score_step: score step
        score: score for which to calculate the PEP

    Returns:
        float: interpolated PEP value
    """
    bin_lo = int((score - min_score) / score_step)
    bin_hi = min(bin_lo + 1, len(pep_bins) - 1)

    lower = pep_bins[bin_lo]
    upper = pep_bins[bin_hi]

    bin_lo_score = bin_lo * score_step + min_score
    linear = (score - bin_lo_score) / score_step

    delta = upper - lower
    return lower + (delta * linear)

# caclulate pep for all scores
@njit
def calculate_pep(scores: NDArray,
                  decoys: NDArray,
                  bins: int = 1000,
                  bw_adjust: float = 1.0,
                  monotonic: bool = True) -> NDArray:
    """Calculate PEP for all scores.

    Args:
        scores: numpy array of scores
        decoys: numpy array of boolean values indicating decoys
        bins: number of bins for the PEP calculation
        bw_adjust: bandwidth adjustment factor
        monotonic: whether to enforce monotonicity

    Returns:
        numpy array: PEP values for all scores
    """
    pep_bins, min_score, score_step = calculate_pep_single(scores, decoys, bins, bw_adjust, monotonic)
    pep = np.zeros(len(scores))
    for i in range(len(scores)):
        pep[i] = posterior_error(pep_bins, min_score, score_step, scores[i])
    return pep

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