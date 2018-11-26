"""Load Pfam data from the MySQL database."""
from pkg_resources import Requirement, resource_filename
import json
import pymysql as mc
import os
from ..Family import Family


def load_config(file_=None):
    """Load database configuration from a file-like object.

    By default, the configuration is loaded from the path stored in the
    environment variable SEARCH_SIFTER_DB_CONFIG.

    Parameters
    ----------
    file_ : file-like object, optional
        A file like object containing a JSON dictionary which will be unpacked
        and passed directly to MySQL as the connection parameters. If not
        given, SEARCH_SIFTER_DB_CONFIG is used to find the configuration.
    """
    if file_ is None:
        try:
            file_ = open(os.environ["SEARCH_SIFTER_DB_CONFIG"])
        except KeyError:
            try:
                file_ = open(resource_filename(
                    Requirement.parse("SearchSifter"),
                    "config/db.json"))
            except FileNotFoundError:
                raise RuntimeError("Couln't load a Pfam DB config")
    with file_ as f:
        return json.load(f)


class Cnx(object):
    def __enter__(self, config=None):
        if config is None:
            config = load_config()
        self.cnx = mc.connect(**config)
        return self.cnx

    def __exit__(self, type, valu7e, traceback):
        self.cnx.close()


_families_query = ("select distinct pfamA_acc, pfamseq_acc "
                   "from pfamA_reg_full_significant "
                   "order by pfamA_acc; ")

_families_regions_query = ("select pfamA_acc, pfamseq_acc, ali_start, ali_end "
                           "from pfamA_reg_full_significant "
                           "order by pfamA_acc; ")

_family_query = ("select pfamseq_acc, ali_start, ali_end "
                 "from pfamA_reg_full_significant "
                 "where pfamA_acc=%s; ")

_clans_query = ("select clan_acc, pfamA_acc "
                "from clan_membership "
                "order by clan_acc; ")


class Classification(object):
    """Iterate over sets of MySQL row values grouped by some key value.

    Parameters
    ----------
    _query : str
        A MySQL SELECT query having the following properties:
            - The leftmost selected column is the key value which will be
              used to group the remaining columns.
            - The query uses the ORDER BY keyword to sort on the key value.
    """

    def __init__(self, _query):
        self._query = _query

    def __iter__(self):
        """Iterate over each distinct key in the query.

        For each distinct key, yield a set of values which have this key.

        Yields
        ------
        key : str
        values : set of tuple
            A tuple composed of each of the value columns (ie, all but the
            leftmost column) in the SELECT query.
        """
        current_group = None
        current_members = None
        with Cnx() as cnx:
            cursor = cnx.cursor()
            try:
                cursor.execute(self._query)
                for row in cursor:
                    group, *memberl = map(_decode, row)
                    member = tuple(memberl)
                    if group != current_group:
                        if current_group is not None:
                            yield current_group, current_members
                        current_group = group
                        current_members = set()
                    current_members.add(member)
                if current_group is not None:
                    yield current_group, current_members
            except:
                raise
            finally:
                cursor.close()


class Families(Classification):
    """Iterate over Pfam families and their members.

    For each Pfam family, yield the family accession and the accessions
    of all proteins with a matching region for the family.
    """
    def __init__(self, *args, _query=_families_query, **kwargs):
        super().__init__(_query, *args, **kwargs)

    def __iter__(self):
        """
        Yields
        ------
        str
            The family's Pfam accession.
        set of str
            The accessions of Pfam family members.
        """
        for family, member in super().__iter__():
            yield family, {m[0] for m in member}


class FamiliesRegions(Classification):
    """Iterate over Pfam families and their members.

    For each Pfam family, yields the accessions and coordinates of each region
    matching the family.
    """
    def __init__(self, *args, _query=_families_regions_query, **kwargs):
        super().__init__(_query, *args, **kwargs)

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
        return super().__iter__()


def _cache_decorator(m):
    # Fetch the clan-family associations, and cache for future queries.
    def new_m(self, *args, **kwargs):
        if self.clans is None:
            self.clans = {k: v for k, v in self}
            self.reversed_clans = {f: c for c, fs in self.clans.items() for f in fs}
        return m(self, *args, **kwargs)
    return new_m


class Clans(Classification):
    """Find the members of a clan.
    """
    def __init__(self, *args, _query=_clans_query, **kwargs):
        super().__init__(_query, *args, **kwargs)
        self.clans = None

    def __iter__(self):
        """Iterate over Pfam clans and their members.

        For each clan, yields the set of Pfam family accessions which comprise
        it.

        Yields
        ------
        str
            The clans's Pfam accesion.
        set of str
            The accessions of the Pfam families in the clan.
        """
        for clan, members in super().__iter__():
            yield(clan, {m[0] for m in members})

    @_cache_decorator
    def families_in_clan(self, clan):
        """Get the members of a clan.

        Parameters
        ----------
        clan : str

        Returns
        -------
        set of str
        """
        return self.clans[clan]

    @_cache_decorator
    def clan_for_family(self, family):
        """Get the clan of which a family is a member, if any.

        Parameters
        ----------
        family : str

        Returns
        -------
        str or None
            If the family is a member of a clan, returns the clan's accession,
            otherwise None.
        """
        try:
            return self.reversed_clans[family]
        except KeyError:
            return None


class PfamFamily(Family):
    """Fetch a single Pfam family from the MySQL database.

    The family is characterised by its accession and by the accessions and
    coordinates of the regions of proteins which match the family.

    For each accession, there should be at most one instance of this class.

    Methods
    -------
    from_accession(acc)
        Get the instance for a particular accession.
    """
    _instances = {}

    def __init__(self, accession, _query=_family_query):
        """.. warning:: This class should be instantiated using
        `from_accession`"""
        super().__init__()
        self.accession = accession
        self._query = _query

    @classmethod
    def from_accession(cls, acc):
        """Instantiate `PfamFamily` if necessary, and return it.

        Parameters
        ----------
        acc : str
            A Pfam family accession.

        Returns
        -------
        PfamFamily
        """
        try:
            return cls._instances[acc]
        except KeyError:
            new_instance = cls(acc)
            new_instance._load()
            cls._instances[acc] = new_instance
            new_instance.finalise()
            return new_instance

    def _load(self):
        with Cnx() as cnx:
            cursor = cnx.cursor()
            try:
                cursor.execute(self._query, (self.accession,))
                for row in cursor:
                    proteinb, start, end = row
                    protein = _decode(proteinb)
                    self.add_region(protein, start, end)
            except:
                raise
            finally:
                cursor.close()


class PfamClan(Family):
    _instances = {}

    def __init__(self, accession):
        super().__init__()
        self.accession = accession

    @classmethod
    def from_accession(cls, acc):
        try:
            return cls._instances[acc]
        except KeyError:
            new_instance = cls(acc)
            new_instance._load()
            cls._instances[acc] = new_instance
            new_instance.finalise()
            return new_instance

    def _load(self):
        member_accs = Clans().families_in_clan(self.accession)
        members = [PfamFamily.from_accession(a) for a in member_accs]
        if len(members) > 1:
            clan_regions = members[0].union(*members[1:])
            self._regions = Family(_regions={k: [(v.start, v.stop) for v in vs]
                                             for k, vs
                                             in clan_regions.items()})._regions


def _decode(s):
    try:
        return s.decode("utf-8")
    except AttributeError:
        return s
