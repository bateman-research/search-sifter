import pytest
import searchsifter.relationships.pfam_db as pf


def test_same_family():
    f1 = pf.PfamFamily.from_accession("PF00042")
    f2 = pf.PfamFamily.from_accession("PF00042")
    assert f1 is f2
    assert f1.jaccard_index(f2) == 1
    assert f1.jaccard_containment(f2) == 1

def test_different_family():
    f1 = pf.PfamFamily.from_accession("PF00042")
    f2 = pf.PfamFamily.from_accession("PF01383")
    assert f1.jaccard_index(f2) == 0
    assert f1.jaccard_containment(f2) == 0

def test_same_clan():
    f1 = pf.PfamFamily.from_accession("PF00042")
    f2 = pf.PfamFamily.from_accession("PF11563")
    assert 0 < f1.jaccard_index(f2) < 1
    assert 0 < f1.jaccard_containment(f2) < 1
    assert f1.jaccard_index(f2) == f2.jaccard_index(f1)
    assert f1.jaccard_containment(f2) != f2.jaccard_containment(f1)
