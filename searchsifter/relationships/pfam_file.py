"""Load Pfam data from a flat file."""
import gzip
from ..Family import Family
import sys

PFAM_FILETYPE_STOCKHOLM = "stockholm"
PFAM_FILETYPE_REGIONS = "regions"


class FamiliesRegions(object):
    """Iterate over Pfam families and their members.

    For each Pfam family, yields the accessions and coordinates of each region
    matching the family.
    """
    def __init__(self, filename, filetype=PFAM_FILETYPE_REGIONS):
        self.filename = filename
        self.filetype = filetype

    def __iter__(self):
        """
        Yields
        ------
        str
            The family's Pfam accesion.
        set of (str, int, int)
            A tuple with the first element, a protein's accession, and the
            second and third, the start and end coordinates of the alignment
            of the protein to the Pfam family.
        """
        yield from pfam_file_iter(self.filename, self.filetype)


class Families(object):
    """Iterate over Pfam families and their members.

    For each Pfam family, yield the family accession and the accessions
    of all proteins with a matching region for the family.
    """
    def __init__(self, filename):
        self.families_regions = FamiliesRegions(filename=filename)

    def __iter__(self):
        """
        Yields
        ------
        str
            The family's Pfam accession.
        set of str
            The accessions of Pfam family members.
        """
        for pfam_acc, regions in self.families_regions:
            yield pfam_acc, {acc for acc, _, _ in regions}


class PfamFamilies(object):
    """Load Pfam families from a flat file.

    Parameters
    ----------
    filename : str
        Path to the Pfam file.
    filetype : stockholm (default) or regions
    accs : list of str, optional
        A list of families to load. If not given, all families will be loaded.
    """
    def __init__(self, filename, filetype=PFAM_FILETYPE_STOCKHOLM, accs=None):
        self.filename = filename
        self.filetype = filetype
        self.accs = accs
        self.families = None
        self._load()

    def _load(self):
        self.families = {}
        for current_fam, current_members in pfam_file_iter(self.filename, self.filetype):
            if self.accs is None or current_fam in self.accs:
                new_fam = Family()
                for acc, start, end in current_members:
                    new_fam.add_region(acc, start, end)
                new_fam.finalise()
                try:
                    assert current_fam not in self.families
                except AssertionError:
                    print("Pfam file misformatted at {}, continuing to next family".format(current_fam), file=sys.stderr)
                    continue
                self.families[current_fam] = new_fam

    def __getitem__(self, key):
        """Returns the Family object for the given accession.

        If the family is not present, returns an empty family."""
        try:
            return self.families[key]
        except KeyError:
            print("Family missing, returning empty family", file=sys.stderr)
            empty_fam = Family()
            empty_fam.finalise()
            self.families[key] = empty_fam
            return empty_fam


def pfam_file_iter(filename, filetype):
    if filetype == PFAM_FILETYPE_REGIONS:
        current_fam = None
        current_members = set()
        with gzip.open(filename, 'rt', encoding="latin_1") as pfam_file:
            next(pfam_file)
            for line in pfam_file:
                components = line.split('\t')
                acc = components[0]
                family = components[4]
                start, end = map(int, components[5:7])
                if current_fam is None:
                    current_fam = family
                if family == current_fam:
                    current_members.add((acc, start, end))
                else:
                    yield current_fam, current_members
                    current_fam = family
                    current_members = set()
                    current_members.add((acc, start, end))
            yield current_fam, current_members
    elif filetype == PFAM_FILETYPE_STOCKHOLM:
        acc_by_id = {}

        current_fam = None
        current_members = set()
        with gzip.open(filename, 'rt', encoding="latin_1") as pfam_file:
            for line in pfam_file:
                components = line.split()
                if len(line) and line[0] == "#":
                    if components[0] == "#=GF":
                        if components[1] == "AC":
                            current_fam = components[2].split('.')[0]
                    elif components[0] == "#=GS":
                        if components[2] == 'AC':
                            acc = components[3].split('.')[0]
                            id_ = components[1].split('/')[0]
                            acc_by_id[id_] = acc
                elif components[0] == "//":
                    yield current_fam, current_members
                    acc_by_id = {}
                    current_members = set()
                else:
                    id_, range_ = components[0].split('/')
                    start, end = map(int, range_.split('-'))
                    maybe_acc = id_.split('.')
                    if len(maybe_acc) > 1:
                        acc = maybe_acc[0]
                    else:
                        acc = acc_by_id[id_]
                    current_members.add((acc, start, end))
    else:
        raise RuntimeError
