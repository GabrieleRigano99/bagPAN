from bagpan.stats import benjamini_hochberg, fisher_exact


def test_fisher_exact_matches_hand_computed_hypergeometric():
    # table [[3,1],[1,3]]: row1=row2=col1=4, n=8, C(8,4)=70.
    # P(a) = C(4,a)*C(4,4-a)/70 for a in 0..4 -> [1,16,36,16,1]/70.
    # observed a=3 -> p=16/70; two-sided sum of p(x)<=16/70 is (1+16+16+1)/70
    # = 34/70 = 0.4857142857142857 (this also matches scipy.stats.fisher_exact
    # on the same table - a standard textbook example).
    p = fisher_exact(((3, 1), (1, 3)))
    assert abs(p - 34 / 70) < 1e-9


def test_fisher_exact_symmetric_and_bounded():
    table = ((5, 0), (0, 5))
    p = fisher_exact(table)
    assert 0.0 <= p <= 1.0
    # extremely skewed table -> very small p-value
    assert p < 0.01


def test_fisher_exact_identical_rows_gives_p_one():
    p = fisher_exact(((4, 4), (4, 4)))
    assert abs(p - 1.0) < 1e-9


def test_fisher_exact_degenerate_row_returns_one():
    assert fisher_exact(((0, 0), (3, 5))) == 1.0
    assert fisher_exact(((3, 0), (0, 0))) == 1.0


def test_benjamini_hochberg_hand_computed_example():
    # p = [0.01, 0.02, 0.03, 0.5], n=4, ranks 1..4 (already ascending).
    # q(4)=0.5*4/4=0.5
    # q(3)=0.03*4/3=0.04
    # q(2)=0.02*4/2=0.04
    # q(1)=0.01*4/1=0.04 -> running min stays 0.04
    q = benjamini_hochberg([0.01, 0.02, 0.03, 0.5])
    expected = [0.04, 0.04, 0.04, 0.5]
    for got, want in zip(q, expected):
        assert abs(got - want) < 1e-9


def test_benjamini_hochberg_is_monotonic_and_order_preserving():
    p_values = [0.2, 0.001, 0.9, 0.04]
    q_values = benjamini_hochberg(p_values)
    assert len(q_values) == len(p_values)
    # smallest p must not get a larger q than any p it's compared against
    assert q_values[1] <= q_values[0]
    for q in q_values:
        assert 0.0 <= q <= 1.0


def test_benjamini_hochberg_empty():
    assert benjamini_hochberg([]) == []