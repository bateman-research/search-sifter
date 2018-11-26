from collections import defaultdict
import operator
from functools import reduce
# NB .relationships.minhash is imported at the END of this module to solve a
# circular dependancy issue, using the statement below.

# from .relationships import minhash as mh
# from .relationships import jaccard as jc


class Family(object):
    """
    A protein family.

    This class allows a protein family, defined by contiguous regions on
    one or more proteins, to be compared to other protein families, and to
    Pfam, via MinHash.
    """
    def __init__(self, _regions=None):
        if _regions == None:
            self._regions = defaultdict(list)
        else:
            self._regions = _regions
        self._signatures = {}
        self._fhashes = {}
        self._finalised = False

    def add_region(self, acc, start, end):
        """
        Add a protein sequence region to the protein family.

        This method should only be called before the Family has been finalised.

        Parameters
        ----------
        acc : str
            The accession of the protein to be added to the family.
        start, end : int
            The start and end of the region of the protein sequence which is
            homologous with the protein family.

        Raises
        ------
        RuntimeError
            If `finalise` has been called on the Family.
        """
        if not self._finalised:
            self._regions[acc].append((start, end))
            self._signatures = {}
        else:
            raise RuntimeError()

    def finalise(self):
        """
        Declare the Family immutable.

        Once called, regions can no longer be added to the Family.
        """
        self._finalised = True

    def proteins(self):
        """
        Get the accessions of proteins which have regions in the Family.

        Returns
        -------
        iterable of str
        """
        return self._regions.keys()

    def regions(self):
        """
        Iterate over the Family's regions.

        Yields
        ------
        protein : str
            The protein accession.
        start, end : int
            The start and end coordinates of the region.
        """
        for protein, regions in self._regions.items():
            for start, end in regions:
                yield protein, start, end

    def jaccard_index(self, other):
        """
        Compute the Jaccard index between the proteins which are members of
        this protein family, and another protein family.

        Parameters
        ----------
        other : Family

        Returns
        -------
        float
        """
        return jc.jaccard(self.proteins(), other.proteins())

    def jaccard_containment(self, other):
        """
        Compute the Jaccard index between the proteins which are members of
        this protein family, and another protein family.

        Parameters
        ----------
        other : Family

        Returns
        -------
        float
        """
        return jc.jaccard_containment(self.proteins(), other.proteins())

    def signature(self, hashes):
        """
        Compute the MinHash signature for the members of this protein family.

        Parameters
        ----------
        hashes : relationships.pfam.Hashes
            The Hashes object with which the signature will be used.

        Returns
        -------
        list of tuple
            The first element of the tuples is a hash, the second is the
            object from which that hash was generated.
        """
        try:
            sig = self._signatures[hashes]
        except KeyError:
            sig = hashes.signature(self)
            self._signatures[hashes] = sig
        return sig

    def full_hash(self, hashes):
        try:
            fhash = self._fhashes[hashes]
        except KeyError:
            fhash = hashes.full_hash(self)
            self._fhashes[hashes] = fhash
        return fhash

    def estimate_jaccard_with_pfam(self, hashes):
        """
        Estimate the Jaccard index between the members of this protein family
        and all protein families in Pfam.

        Parameters
        ----------
        hashes : relationships.pfam.Hashes
            The Hashes object containing the signatures with which this family
            will be compared.

        Returns
        -------
        float
        """
        return hashes.estimate_jaccard(self)

    def estimate_containment_with_pfam(self, hashes):
        """
        Estimate the Jaccard containment between the members of this protein
        family and all protein families in Pfam.

        Parameters
        ----------
        hashes : relationships.pfam.Hashes
            The Hashes object containing the signatures with which this family
            will be compared.

        Returns
        -------
        float
        """
        return hashes.estimate_containment(self)

    def overlap(self, other):
        """Compute the residue-level overlap with another family.

        Return the accessions and ranges which both of these families cover.

        Parameters
        ----------
        other : Family

        Returns
        -------
        dict of list of range
            A dictionary keyed by protein accessions for which there are
            overlapping regions. The values are lists of `range`, representing
            the regions for which the key protein has overlaps.
        """
        overlaps = defaultdict(set)
        if len(self.proteins()) < len(other.proteins()):
            f1, f2 = self, other
        else:
            f2, f1 = self, other

        intersection = f1.proteins() & f2.proteins()
        for p1, p2 in zip(_dict_slice(f1._regions, intersection),
                          _dict_slice(f2._regions, intersection)):
            acc1, r1 = p1
            acc2, r2 = p2
            assert acc1 == acc2
            sorted_r1, sorted_r2 = map(sorted, [r1, r2])
            for s1, e1 in sorted_r1:
                for s2, e2 in sorted_r2:
                    if s2 > e1:
                        break
                    o = overlap(s1, e1, s2, e2)
                    if len(o):
                        overlaps[acc1].add(o)
            overlaps[acc1] = _merge_ranges(overlaps[acc1])
        return overlaps

    def union(self, *others):
        """Compute the residue-level union with another family.

        Return the accessions and ranges which either of these families cover.

        Parameters
        ----------
        other : searchsifter.Family

        Returns
        -------
        dict of list of range
            A dictionary keyed by protein accessions for which there are
            regions covered by either family. The values are lists of
            `range`, representing the regions for which the key protein is
            covered.
        """
        families = [self] + list(others)
        proteins = set(reduce(operator.or_, map(lambda o: o.proteins(), families)))
        out = {}
        for protein in proteins:
            rs = map(lambda r: r._regions.get(protein), families)
            rs = map(lambda x: set(x) if x is not None else set(), rs)
            rs = {range(s, e + 1) for s, e in reduce(operator.or_, rs)}
            u = _merge_ranges(rs)
            out[protein] = u
        return out

    def proteins_covered(self):
        return len(self.proteins())

    def residues_covered(self):
        """Compute the number of residues covered by this family.

        Returns
        -------
        int
            The number of residues covered by this family.
        """
        return sum(len(r) for _, rs in self._regions.items()
                          for r in _merge_ranges(range(s, e + 1)
                                                 for s, e in rs))


def _dict_slice(d, s):
    return sorted((k, v) for k, v in d.items() if k in s)


def _merge_ranges(ranges):
    # Given a list of range, produce a list which covers the same ranges,
    # with the fewest range objects possible.
    u = []
    # Iterate through the ranges and combine those which overlap.
    for r in sorted(ranges, key=lambda r: (r.start, r.stop)):
        if len(u) == 0:
            u.append(r)
        else:
            r1 = u[-1]
            r2 = r
            if (len(overlap(r1.start, r1.stop - 1, r2.start, r2.stop - 1)) or
                r1.stop == r2.start):

                u.pop()
                u.append(range(min(r1.start, r2.start), max(r1.stop, r2.stop)))
            else:
                u.append(r)
    return u


def overlap(start_1, end_1, start_2, end_2):
    """Return the `range` covered by two sets of coordinates.

    The coordinates should be supplied inclusive, that is, the end coordinates
    are included in the region. The `range` returned will be exclusive, in
    keeping with the correct usage of that type.

    Parameters
    ----------
    start_1, start_2 : int
        The start coordinates of each region.
    end_1, end_2 : int
        The end coordinates of each region.

    Returns
    -------
    range
        The `range` covered by both regions. If there is no overlap then start
        of the range will be equal to or greater than the end, thus having zero
        or negative length.
    """
    return range(max(start_1, start_2),
                 min(end_1, end_2) + 1)

from .relationships import minhash as mh
from .relationships import jaccard as jc
