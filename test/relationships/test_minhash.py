import pytest
import searchsifter.relationships.minhash as mh
import searchsifter.relationships.jaccard as jc


@pytest.mark.parametrize("a, b, result", [
    ({1, 2}, {2}, 0.5),
    ({1, 2}, {2, 3}, 1/3),
    ({1}, {2}, 0),
    ({1}, {1}, 1)
])
def test_jaccard(a, b, result):
    assert jc.jaccard(a, b) == result
    assert jc.jaccard(b, a) == result


@pytest.mark.parametrize("a, b, result", [
    ({1, 2}, {2}, 0.5),
    ({1}, {1, 2}, 1),
    ({1, 2}, {2, 3}, 0.5),
    ({1}, {2}, 0),
    ({1}, {1}, 1),
])
def test_jaccard_containment(a, b, result):
    assert jc.jaccard_containment(a, b) == result


@pytest.mark.parametrize("a, b, result", [
    ({1, 2}, {2}, 0.5),
    ({1, 2}, {2, 3}, 1/3),
    ({1}, {2}, 0),
    ({1}, {1}, 1)
])
def test_minhash(a, b, result):
    s, t = mh.signature(a, 5), mh.signature(b, 5)
    assert mh.minhash(s, t, 5) == result
    assert mh.minhash(t, s, 5) == result


@pytest.mark.parametrize("a, b, result", [
    ({1, 2}, {2}, 0.5),
    ({1}, {1, 2}, 1),
    ({1, 2}, {2, 3}, 0.5),
    ({1}, {2}, 0),
    ({1}, {1}, 1),
])
def test_minhash_containment(a, b, result):
    s, t = mh.signature(a, 5), mh.signature(b, 5)
    assert mh.minhash_containment(s, t) == result


@pytest.fixture
def a():
    return set(range(1, 100))


@pytest.fixture
def b():
    return set(range(50, 100))


@pytest.fixture
def c():
    return set(range(75, 100))


def test_intersection(a, b, c):
    assert mh.intersection_signature(a, b) == set(range(50, 100))
    assert mh.intersection_signature(a, b, c) == set(range(75, 100))


def test_union(a, b):
    assert mh.union_signature(a, b, 100) == a
    assert len(mh.union_signature(a, b, 20)) == 20
