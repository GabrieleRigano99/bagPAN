"""Pangenome accumulation (rarefaction) curve and a Heaps'-law-style
power-law fit of pangenome "openness" - pure Python (stdlib `statistics`,
`itertools`, `random`, `math` only), no scipy/numpy needed.

For k = 1..N species, averages over many random species orderings (all
N! orderings when that's small enough, otherwise a random sample - mirroring
FunFinder_Pangenome.py's own --max_subsamples pattern) the core-genome size
(orthogroups present in every one of the first k species in that ordering)
and pan-genome size (orthogroups present in at least one of them).
"""

from __future__ import annotations

import itertools
import math
import random
import statistics
from typing import Dict, List, Optional, Tuple

from bagpan.orthogroups import OrthogroupSet


def _orthogroups_by_species(orthogroup_set: OrthogroupSet) -> Dict[str, set]:
    by_species: Dict[str, set] = {name: set() for name in orthogroup_set.species_names}
    for cluster, present in orthogroup_set.species_present.items():
        for species in present:
            by_species[species].add(cluster)
    return by_species


def accumulation_curve(
    orthogroup_set: OrthogroupSet, n_permutations: int = 200, seed: int = 0
) -> Dict[int, Dict[str, float]]:
    """Returns {k: {"core_mean", "core_std", "pan_mean", "pan_std"}} for
    k = 1..N, N = number of species in orthogroup_set.
    """
    species_names = list(orthogroup_set.species_names)
    n = len(species_names)
    by_species = _orthogroups_by_species(orthogroup_set)

    rng = random.Random(seed)
    if math.factorial(n) <= n_permutations:
        orderings: List[Tuple[str, ...]] = list(itertools.permutations(species_names))
    else:
        seen = set()
        orderings = []
        while len(orderings) < n_permutations:
            order = tuple(rng.sample(species_names, n))
            if order not in seen:
                seen.add(order)
                orderings.append(order)

    core_by_k: Dict[int, List[int]] = {k: [] for k in range(1, n + 1)}
    pan_by_k: Dict[int, List[int]] = {k: [] for k in range(1, n + 1)}

    for order in orderings:
        core_set: Optional[set] = None
        pan_set: Optional[set] = None
        for k, species in enumerate(order, start=1):
            sp_orthogroups = by_species[species]
            if core_set is None:
                core_set = set(sp_orthogroups)
                pan_set = set(sp_orthogroups)
            else:
                core_set &= sp_orthogroups
                pan_set |= sp_orthogroups
            core_by_k[k].append(len(core_set))
            pan_by_k[k].append(len(pan_set))

    return {
        k: {
            "core_mean": statistics.mean(core_by_k[k]),
            "core_std": statistics.pstdev(core_by_k[k]),
            "pan_mean": statistics.mean(pan_by_k[k]),
            "pan_std": statistics.pstdev(pan_by_k[k]),
        }
        for k in range(1, n + 1)
    }


def fit_power_law(xs: List[float], ys: List[float]) -> Optional[Tuple[float, float]]:
    """Closed-form log-log least-squares fit of y = kappa * x^gamma.
    Returns (kappa, gamma), or None if there isn't enough data to fit.
    """
    pairs = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(pairs) < 2:
        return None

    n = len(pairs)
    sum_x = sum(p[0] for p in pairs)
    sum_y = sum(p[1] for p in pairs)
    sum_xx = sum(p[0] * p[0] for p in pairs)
    sum_xy = sum(p[0] * p[1] for p in pairs)

    denom = n * sum_xx - sum_x * sum_x
    if denom == 0:
        return None

    gamma = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - gamma * sum_x) / n
    return math.exp(intercept), gamma