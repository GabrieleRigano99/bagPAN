"""Standard-library-only replacements for the scipy/statsmodels calls
FunFinder_Pangenome.py used (stats.fisher_exact, multipletests(method='fdr_bh')),
so bagpan doesn't need scipy/statsmodels installed.
"""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple

Table2x2 = Tuple[Tuple[int, int], Tuple[int, int]]


def _hypergeom_pmf(a: int, row1: int, row2: int, col1: int) -> float:
    n = row1 + row2
    return (math.comb(row1, a) * math.comb(row2, col1 - a)) / math.comb(n, col1)


def fisher_exact(table: Table2x2) -> float:
    """Two-sided Fisher's exact test p-value for a 2x2 contingency table,
    via direct enumeration of the hypergeometric distribution (fine for the
    small orthogroup counts this tool works with).
    """
    (a, b), (c, d) = table
    row1, row2 = a + b, c + d
    col1 = a + c
    n = row1 + row2
    if n == 0 or row1 == 0 or row2 == 0 or col1 == 0 or col1 == n:
        return 1.0
    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    observed_p = _hypergeom_pmf(a, row1, row2, col1)
    # small relative+absolute tolerance so floating point noise doesn't drop
    # a table equally-as-extreme as the observed one from the two-sided sum
    tolerance = observed_p * 1e-7 + 1e-12
    p_value = 0.0
    for x in range(lo, hi + 1):
        p_x = _hypergeom_pmf(x, row1, row2, col1)
        if p_x <= observed_p + tolerance:
            p_value += p_x
    return min(p_value, 1.0)


def benjamini_hochberg(p_values: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg FDR-adjusted q-values, in the same order as the
    input p-values.
    """
    n = len(p_values)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: p_values[i])
    sorted_p = [p_values[i] for i in order]

    q_sorted = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        i = rank - 1
        q = sorted_p[i] * n / rank
        running_min = min(running_min, q)
        q_sorted[i] = running_min

    result = [0.0] * n
    for sorted_pos, orig_idx in enumerate(order):
        result[orig_idx] = min(q_sorted[sorted_pos], 1.0)
    return result