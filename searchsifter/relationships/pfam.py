"""Generate MinHash signatures from Pfam, or load from disk.
"""
import json
from searchsifter.relationships import minhash as mh
from ..Family import Family
from collections import defaultdict
from . import pfam_db

pfam_db = None


class Hashes(object):
    """Generate or load Pfam MinHash signatures from disk.

    For each Pfam family, the signature is calculated using the accessions
    of proteins which have regions matching the family.

    If signatures are not loaded from disk, they will be generated from the
    MySQL database, or from a Pfam flat file.

    Parameters
    ----------
    n : int
        The signature length. If `from_file` is specified, signatures will be
        truncated to length `n`.
    from_file : file_like
        A JSON encoded file produced by `save_to_file` from which signatures
        should be loaded.
    pfam : module
        The module which should be used for loading Pfam data, i.e., pfam_db
        or pfam_file.
    pfam_args : list
        The args which should be supplied to the various constructors in the
        module providing Pfam data. If using the MySQL database, this should
        be empty, if using the flat file it should contain the filename of the
        file.
    """
    def __init__(self, n=100, from_file=None, pfam=pfam_db, pfam_args=None, _hashes=None):
        self.n = n
        if pfam_args is None:
            pfam_args = []
        self.pfam_args = pfam_args
        self.pfam = pfam
        if from_file is not None:
            self.load_from_file(from_file)
        elif _hashes is not None:
            self._hashes = _hashes
        else:
            self._hashes = None

    def __iter__(self):
        """Iterate over Pfam family accessions and family signatures.

        Yields
        ------
        str
            A Pfam accession.
        (int, object)
            The first element of the tuple is the signature corresponding to
            the Pfam family. The second element is the object which was hashed
            to produce the signature. For example, a protein accession.
        """
        return iter(self.hashes.items())

    @property
    def hashes(self):
        """Get the dictionary of Pfam family accessions and signatures.

        If the signatures have already been generated, they are cached.

        Returns
        -------
        dict
            key-value pairs described in `__iter__`.
        """
        if self._hashes is None:
            self._hashes = {}
            for family, proteins in self.pfam.Families(*self.pfam_args):
                self._hashes[family] = mh.signature(proteins, self.n)
        return self._hashes

    def save_to_file(self, filename):
        """Save the signatures to a path.

        There must not already be a file located at `filename`

        Parameters
        ----------
        filename : str"""
        with open(filename, 'x') as save_file:
            json.dump({k: list(v) for k, v in self.hashes.items()}, save_file)

    def load_from_file(self, hash_file, n=None):
        """Load the signatures from a file.

        Parameters
        ----------
        hash_file : file_like"""
        if n is None:
            n = self.n
        self._hashes = {f: set(map(tuple, sorted(p)[:n])) for f, p
                        in json.load(hash_file).items()}

    def estimate_jaccard(self, family):
        """Estimate the Jaccard index between a family and Pfam.

        For each family in Pfam, the an estimate for the index is computed
        using MinHash.

        Parameters
        ----------
        family : searchsifter.Family

        Returns
        -------
        dict
            A dictionary with `str` Pfam family accession keys and `float`
            estimated Jaccard Index values."""
        A = family.signature(self)
        jaccards = {}

        for acc, B in self:
            jaccards[acc] = mh.minhash(A, B, self.n)
        return jaccards

    def estimate_containment(self, family):
        """Estimate the Jaccard containment between a family and Pfam.

        For each family in Pfam, the an estimate for the containment is
        computed using MinHash.

        Parameters
        ----------
        family : searchsifter.Family

        Returns
        -------
        dict
            A dictionary with `str` Pfam family accession keys and `float`
            estimated Jaccard Containment values."""
        B = family.full_hash(self)
        containments = {}

        for acc, A in self:
            containments[acc] = mh.minhash_containment(A, B)
        return containments

    def signature(self, family):
        """Get the MinHash signature for a Family.

        Parameters
        ----------
        family : searchsifter.Family

        Returns
        -------
        (int, object)
            The first element of the tuple is the signature corresponding to
            the Pfam family. The second element is the object which was hashed
            to produce the signature. For example, a protein accession.
        """
        return mh.signature(family.proteins(), self.n)

    def full_hash(self, family):
        return mh.set_hashes(family.proteins())

    def clan_hashes(self):
        """Get a `Hashes` object containing signatures for clans.

        This is calculated as the union signature of families which are
        members of the clan.

        Returns
        -------
        Hashes
        """
        chs = self._clan_hashes()
        return type(self)(self.n, _hashes=chs)

    def _clan_hashes(self):
        chs = {}
        for clan, members in self.pfam.Clans(*self.pfam_args):
            hs = [self.hashes[f] for f in members if f in self.hashes]
            if len(hs) == 1:
                chs[clan] = hs[0]
            elif len(hs) > 1:
                chs[clan] = mh.union_signature(*(hs + [self.n]))
        return chs


class ResidueHashes(Hashes):
    """Generate or load Pfam MinHash signatures from disk.

    For each Pfam family, the signature is calculated using the accessions of
    proteins which have regions matching the family, and a numbered 'chunk' of
    that protein, corresponding to the coordinates in which the protein matches
    the family.

    For example, if PF99999 is matched by protein P12345 from residue 52 to 80
    and the window size is 25, P12345 matches PF99999 in chunks 3 and 4. Hence
    the hashes of (P12345, 3) and (P12345, 4) will possibly be included in the
    hash of PF99999.

    Parameters
    ----------
    w : int
        The window size.
    n : int
        The signature length. If `from_file` is specified, signatures will be
        truncated to length `n`.
    from_file : file_like
        A JSON encoded file produced by `save_to_file` from which signatures
        should be loaded.
    """
    def __init__(self, w, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.w = w

    @property
    def hashes(self):
        """See `Hashes.hashes`."""
        if self._hashes is None:
            self._hashes = {}
            for family, regions in self.pfam.FamiliesRegions(*self.pfam_args):
                ci = _chunk_iterator(regions, self.w)
                self._hashes[family] = mh.signature(set(ci), self.n)
        return self._hashes

    def signature(self, family):
        """See `Hashes.signature`."""
        return mh.signature(set(_chunk_iterator(family.regions(), self.w)),
                            self.n)

    def full_hash(self, family):
        return mh.set_hashes(set(_chunk_iterator(family.regions(), self.w)))

    def load_from_file(self, hash_file, n=None):
        """See `Hashes.load_from_file`."""
        if n is None:
            n = self.n
        self._hashes = {f: {(h, tuple(c)) for h, c in sorted(p)[:n]} for f, p
                        in json.load(hash_file).items()}

    def clan_hashes(self):
        return type(self)(self.w, self.n, _hashes=self._clan_hashes())

    @classmethod
    def hashes_with_windows(cls, ws, n, pfam=pfam_db, pfam_args=None):
        """Create ResidueHashes objects with different values of `w`."""
        if pfam_args is None:
            pfam_args = []
        hs = defaultdict(dict)
        for family, regions in pfam.FamiliesRegions(*pfam_args):
            for w in ws:
                ci = _chunk_iterator(regions, w)
                hs[w][family] = mh.signature(set(ci), n)
        return {w: cls(w, n, _hashes=h) for w, h in hs.items()}


class Sizes(object):
    def __init__(self, from_file=None, pfam=pfam_db, pfam_args=None):
        if from_file is not None:
            self.load_from_file(from_file)
        else:
            self._sizes = None
            if pfam_args is None:
                pfam_args = []
            self.pfam_args = pfam_args
            self.pfam = pfam

    @property
    def sizes(self):
        if self._sizes is None:
            self._sizes = {}
            for family, proteins in self.pfam.Families(*self.pfam_args):
                self._sizes[family] = len(proteins)
        return self._sizes

    # This can't be an instance variable. Otherwise it binds to the instance
    # somehow!
    @property
    def size_method(self):
        return Family.proteins_covered

    def save_to_file(self, filename):
        with open(filename, 'x') as size_file:
            json.dump(self.sizes, size_file)

    def load_from_file(self, size_file):
        self._sizes = json.load(size_file)


class ResidueSizes(Sizes):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def sizes(self):
        if self._sizes is None:
            self._sizes = {}
            for family, regions in self.pfam.FamiliesRegions(*self.pfam_args):
                self._sizes[family] = sum(e - s + 1 for _, s, e in regions)
        return self._sizes

    @property
    def size_method(self):
        return Family.residues_covered


def chunks_from_coordinates(start, end, w):
    """Given start and end coordinates and window size, calculate chunks.

    Parameters
    ----------
    start, end : int
        The start and end coordinates.
    w : int
        The window size.
    """
    return(list(range(start // w, end // w + 1)))


def _decode(s):
    try:
        return s.decode("utf-8")
    except AttributeError:
        return s


def _chunk_iterator(regions, w):
    for acc, start, end in regions:
        for chunk in chunks_from_coordinates(start, end, w):
            yield acc, chunk
