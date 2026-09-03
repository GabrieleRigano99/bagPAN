import math

from bagpan.curves import accumulation_curve, fit_power_law
from bagpan.orthogroups import classify_orthogroups

LOCUS_PREFIX_TO_SPECIES = {"X": "sx", "Y": "sy", "Z": "sz"}
SPECIES = ["sx", "sy", "sz"]

CLUSTERS = {
    "OG1": ["X_1", "Y_1", "Z_1"],  # core
    "OG2": ["X_2", "Y_2"],  # accessory: sx, sy
    "OG3": ["X_3"],  # singleton: sx
    "OG4": ["Z_2"],  # singleton: sz
}


def _build_set():
    return classify_orthogroups(CLUSTERS, LOCUS_PREFIX_TO_SPECIES, SPECIES)


def test_accumulation_curve_hand_computed():
    # 3 species -> 3! = 6 permutations, all enumerated exactly (deterministic,
    # no sampling), so these means can be hand-computed exactly:
    #   sx orthogroups = {OG1,OG2,OG3} (3)   sy = {OG1,OG2} (2)   sz = {OG1,OG4} (2)
    # k=1: each species is "first" in 2/6 permutations -> mean = (3+2+2)/3 = 7/3
    # k=2: each of the 3 unordered pairs is "the first two" in 2/6 perms:
    #   {sx,sy}: core=|{OG1,OG2}|=2 pan=|{OG1,OG2,OG3}|=3
    #   {sx,sz}: core=|{OG1}|=1     pan=|{OG1,OG2,OG3,OG4}|=4
    #   {sy,sz}: core=|{OG1}|=1     pan=|{OG1,OG2,OG4}|=3
    #   mean core = (2+1+1)/3 = 4/3   mean pan = (3+4+3)/3 = 10/3
    # k=3 (all species, same every permutation): core=|{OG1}|=1  pan=|all 4|=4
    result = accumulation_curve(_build_set(), n_permutations=200, seed=0)

    assert set(result) == {1, 2, 3}
    assert abs(result[1]["core_mean"] - 7 / 3) < 1e-9
    assert abs(result[1]["pan_mean"] - 7 / 3) < 1e-9  # k=1: core == pan
    assert abs(result[2]["core_mean"] - 4 / 3) < 1e-9
    assert abs(result[2]["pan_mean"] - 10 / 3) < 1e-9
    assert abs(result[3]["core_mean"] - 1) < 1e-9
    assert abs(result[3]["pan_mean"] - 4) < 1e-9


def test_accumulation_curve_pan_is_monotonic_nondecreasing_and_core_nonincreasing():
    result = accumulation_curve(_build_set(), n_permutations=200, seed=1)
    ks = sorted(result)
    for a, b in zip(ks, ks[1:]):
        assert result[b]["pan_mean"] >= result[a]["pan_mean"] - 1e-9
        assert result[b]["core_mean"] <= result[a]["core_mean"] + 1e-9


def test_fit_power_law_recovers_known_parameters():
    kappa_true, gamma_true = 2.0, 0.5
    xs = [1, 2, 3, 4, 5, 6, 7, 8]
    ys = [kappa_true * (x**gamma_true) for x in xs]
    fitted = fit_power_law(xs, ys)
    assert fitted is not None
    kappa_fit, gamma_fit = fitted
    assert abs(kappa_fit - kappa_true) < 1e-6
    assert abs(gamma_fit - gamma_true) < 1e-6


def test_fit_power_law_insufficient_data_returns_none():
    assert fit_power_law([1], [1]) is None
    assert fit_power_law([], []) is None


def test_fit_power_law_ignores_nonpositive_points():
    # a zero/negative y (degenerate edge case) shouldn't crash log()
    fitted = fit_power_law([1, 2, 3], [0, 4, 9])
    assert fitted is not None
    assert all(math.isfinite(v) for v in fitted)