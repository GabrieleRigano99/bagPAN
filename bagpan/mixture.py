"""A small pure-Python 1D two-component Gaussian mixture (EM), used as an
opt-in alternative (`--classification-method mixture`) to bagpan's fixed
`--percent` threshold for deciding whether an orthogroup "has" a functional
category. In the spirit of PPanGGOLiN's frequency-based partitioning
(persistent/shell/cloud via a statistical model rather than an arbitrary
cutoff) - a deliberately simplified analogue: a plain 1D Gaussian mixture
over per-orthogroup hit-fractions, not PPanGGOLiN's full graph/Markov-Random-
Field machinery.
"""

from __future__ import annotations

import math
import random
from typing import List, Tuple

_EPS = 1e-6


def _gaussian_pdf(x: float, mean: float, stdev: float) -> float:
    stdev = stdev if stdev > _EPS else _EPS
    coeff = 1.0 / (stdev * math.sqrt(2 * math.pi))
    exponent = -((x - mean) ** 2) / (2 * stdev * stdev)
    return coeff * math.exp(exponent)


def fit_gaussian_mixture_1d(
    values: List[float], n_iter: int = 100, seed: int = 0
) -> Tuple[List[float], List[float], List[float]]:
    """Fit a 2-component 1D Gaussian mixture via EM. Returns
    (means, stdevs, weights), each a length-2 list, sorted by ascending mean.
    """
    vals = list(values)
    n = len(vals)
    if n < 2:
        raise ValueError("need at least 2 values to fit a 2-component mixture")

    lo, hi = min(vals), max(vals)
    if hi == lo:
        return [lo, lo], [_EPS, _EPS], [0.5, 0.5]

    means = [lo + 0.25 * (hi - lo), lo + 0.75 * (hi - lo)]
    stdevs = [(hi - lo) / 4, (hi - lo) / 4]
    weights = [0.5, 0.5]
    random.Random(seed)  # reserved for future randomized restarts

    for _ in range(n_iter):
        responsibilities = []
        for x in vals:
            densities = [weights[k] * _gaussian_pdf(x, means[k], stdevs[k]) for k in range(2)]
            total = sum(densities) or _EPS
            responsibilities.append([d / total for d in densities])

        for k in range(2):
            rk = sum(r[k] for r in responsibilities)
            if rk < _EPS:
                continue
            mean_k = sum(r[k] * x for r, x in zip(responsibilities, vals)) / rk
            var_k = sum(r[k] * (x - mean_k) ** 2 for r, x in zip(responsibilities, vals)) / rk
            means[k] = mean_k
            stdevs[k] = math.sqrt(var_k) if var_k > _EPS else _EPS
            weights[k] = rk / n

    order = sorted(range(2), key=lambda k: means[k])
    return [means[k] for k in order], [stdevs[k] for k in order], [weights[k] for k in order]


def classify_by_mixture(values: List[float], n_iter: int = 100, seed: int = 0) -> List[bool]:
    """Fit the mixture and assign each value to the higher-mean component
    (posterior probability >= 0.5) as the "hit" call. Degenerates to a plain
    >= 0.5 cutoff when there's no real bimodal signal to fit (too few points,
    or every value identical).
    """
    vals = list(values)
    if len(vals) < 2 or min(vals) == max(vals):
        return [v >= 0.5 for v in vals]

    means, stdevs, weights = fit_gaussian_mixture_1d(vals, n_iter=n_iter, seed=seed)
    result = []
    for x in vals:
        densities = [weights[k] * _gaussian_pdf(x, means[k], stdevs[k]) for k in range(2)]
        total = sum(densities) or _EPS
        posterior_high = densities[1] / total
        result.append(posterior_high >= 0.5)
    return result