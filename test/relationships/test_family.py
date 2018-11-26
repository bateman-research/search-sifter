import pytest
from searchsifter.Family import Family, _merge_ranges


@pytest.fixture
def f1():
    f = Family()
    f.add_region('a1', 10, 20)
    f.add_region('a2', 10, 20)
    return f


@pytest.fixture
def f2():
    f = Family()
    f.add_region('a1', 15, 20)
    f.add_region('a1', 30, 40)
    return f


@pytest.fixture
def f3():
    f = Family()
    f.add_region('a1', 21, 50)
    f.add_region('a1', 5, 9)
    return f


def test_union_1(f1, f2):
    result = {'a1' : [range(10, 21), range(30, 41)],
              'a2' : [range(10, 21)]}
    assert f1.union(f2) == f2.union(f1) == result


def test_union_2(f1, f3):
    result = {'a1' : [range(5, 51)], 'a2' : [range(10, 21)]}
    assert f1.union(f3) == f3.union(f1) == result


def test_union_3(f1, f2, f3):
    result = {'a1' : [range(5, 51)], 'a2': [range(10, 21)]}
    assert (f1.union(f2, f3) ==
            f1.union(f3, f2) ==
            f2.union(f1, f3) ==
            f2.union(f3, f1) ==
            f3.union(f1, f2) ==
            f3.union(f2, f1) == result )


def test_union_4(f1):
    result = {'a1' : [range(10, 21)], 'a2' : [range(10, 21)]}
    assert f1.union() == result


@pytest.mark.parametrize("rs, result", [
    ([range(0, 5), range(3, 10)], [range(0, 10)]),
    ([range(0, 10), range(2, 5)], [range(0, 10)]),
    ([range(3, 10), range(0, 10)], [range(0, 10)]),
    ([range(8, 10), range(1, 2)], [range(1, 2), range(8, 10)]),
    ([range(0, 5), range(5, 10)], [range(0, 10)]),
    ([range(0, 2), range(2, 8), range(8, 10)], [range(0, 10)]),
    ([range(0, 0), range(0, 0)], [range(0, 0)]),
])
def test_overlap(rs, result):
    assert _merge_ranges(rs) == result
