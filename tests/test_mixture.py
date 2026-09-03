from bagpan.mixture import classify_by_mixture, fit_gaussian_mixture_1d


def _bimodal_data():
    low = [0.01, 0.02, 0.03, 0.04, 0.05, 0.02, 0.03]
    high = [0.90, 0.93, 0.95, 0.97, 0.92, 0.96, 0.91]
    return low, high


def test_fit_gaussian_mixture_separates_bimodal_clusters():
    low, high = _bimodal_data()
    means, stdevs, weights = fit_gaussian_mixture_1d(low + high)
    assert len(means) == 2 and means[0] < means[1]
    assert abs(means[0] - sum(low) / len(low)) < 0.05
    assert abs(means[1] - sum(high) / len(high)) < 0.05
    assert abs(sum(weights) - 1.0) < 1e-6


def test_classify_by_mixture_matches_known_groups():
    low, high = _bimodal_data()
    values = low + high
    labels = classify_by_mixture(values)
    assert labels[: len(low)] == [False] * len(low)
    assert labels[len(low):] == [True] * len(high)


def test_classify_by_mixture_degenerate_identical_values():
    assert classify_by_mixture([0.3, 0.3, 0.3]) == [False, False, False]
    assert classify_by_mixture([0.7, 0.7]) == [True, True]


def test_classify_by_mixture_too_few_values_falls_back_to_cutoff():
    assert classify_by_mixture([]) == []
    assert classify_by_mixture([0.9]) == [True]
    assert classify_by_mixture([0.1]) == [False]


def test_fit_gaussian_mixture_requires_at_least_two_values():
    import pytest

    with pytest.raises(ValueError):
        fit_gaussian_mixture_1d([0.5])