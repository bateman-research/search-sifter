from searchsifter.relationships.pfam import Hashes, PfamFamily, Clans
from searchsifter import Family
import operator
from functools import reduce
import abc
import copy


SIFTS = "sifts"
SEARCH = "search"
ID = "id"
FAMILY_HASH_SCORES = "family_hash_scores"
NORMALISED_FAMILY_HASH_SCORES = "normalised_family_hash_scores"
CLAN_HASH_SCORES = "clan_hash_scores"
CLANS = "clans"
ABSOLUTE_SIZES = "absolute_sizes"
RELATIVE_SIZES = "relative_sizes"
FAMILY_RELATIVE_SIZES = "family_relative_sizes"
SEARCH_SIZE = "search_size"


class Sifter(abc.ABC):
    def __init__(self, **sinks):
        self.sinks = sinks

    def _distribute_results(self, **sinks):
        for sink, data in sinks.items():
            self.sinks[sink].sift(data)

    @abc.abstractmethod
    def sift(self, data):
        out_data = copy.copy(data)
        out_data[SIFTS].append(self)
        return out_data


def package(search, id_):
    return {SEARCH: search,
            ID: id_,
            SIFTS: []}


class HashSifter(Sifter):
    def __init__(self, hashes, low_threshold, high_threshold,
                 low_sink, high_sink, other_sink, sift_method=None,
                 sizes=None):
        if sift_method is None:
            raise RuntimeError
        super().__init__(low_sink=low_sink,
                         high_sink=high_sink,
                         other_sink=other_sink)
        self.hashes = hashes
        self.sift_method = sift_method
        self.low_threshold = low_threshold
        self.high_threshold = high_threshold
        self.sizes = sizes

    def sift(self, data):
        out_data = super().sift(data)
        scores = self.sift_method(data[SEARCH], self.hashes)
        out_data[FAMILY_HASH_SCORES] = {f: s for f, s in scores.items() if s > 0}

        # If *all* families are below the low threshold, send the result to
        # the low sink. Else, if *any* of the families are above the high
        # threshold, send to the high sink. Else, send to the other sink.
        if reduce(operator.and_, (v <= self.low_threshold for v in scores.values())):
            self._distribute_results(low_sink=out_data)
        elif (reduce(operator.or_, (v >= self.high_threshold for v in scores.values()))
            and sum(v >= self.low_threshold for v in scores.values()) == 1):
            self._distribute_results(high_sink=out_data)
        else:
            self._distribute_results(other_sink=out_data)


class Terminator(Sifter):
    def __init__(self):
        self.results = {}

    def sift(self, data):
        out_data = super().sift(data)
        self.results[data[ID]] = out_data


def make_hash_sifter(sift_method):
    class S(HashSifter):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, sift_method=sift_method, **kwargs)
    return S


EstimateJaccardSifter = make_hash_sifter(
    Family.estimate_jaccard_with_pfam)

EstimateContainmentSifter = make_hash_sifter(
    Family.estimate_containment_with_pfam)


class ClanSifter(Sifter):
    def __init__(self, threshold, same_sink, different_sink):
        super().__init__(same_sink=same_sink, different_sink=different_sink)
        self.threshold = threshold

    def sift(self, data):
        out_data = super().sift(data)
        threshold_families = {f: s for f, s in data[FAMILY_HASH_SCORES].items() if s > self.threshold}
        if len(threshold_families) <= 1:
            self._distribute_results(same_sink=out_data)

        clans = {f: Clans().clan_for_family(f)
                    for f, s in threshold_families.items()}
        out_data[CLANS] = clans
        if len(set(clans.values())) == 1 and None not in clans:
            self._distribute_results(same_sink=out_data)
        else:
            self._distribute_results(different_sink=out_data)


class SizeSifter(Sifter):
    def __init__(self, sizes, bigger_sink, smaller_sink,
                 abs_threshold=None, rel_threshold=None):
        super().__init__(bigger_sink=bigger_sink, smaller_sink=smaller_sink)
        if abs_threshold is None and rel_threshold is None:
            raise ValueError
        self.abs_threshold = abs_threshold
        self.rel_threshold = rel_threshold
        self.sizes = sizes

    def sift(self, data):
        out_data = super().sift(data)

        abs_sizes = {f: self.sizes.sizes[f] for f, _ in data[FAMILY_HASH_SCORES].items()}
        out_data[ABSOLUTE_SIZES] = abs_sizes
        search_size = self.sizes.size_method(data[SEARCH])
        out_data[SEARCH_SIZE] = search_size
        rel_sizes = {f: search_size / s for f, s in abs_sizes.items()}
        out_data[RELATIVE_SIZES] = rel_sizes
        normalised_scores = {f: (s * self.sizes.sizes[f]) / search_size
                             for f, s in data[FAMILY_HASH_SCORES].items()
                             if s > 0}
        out_data[NORMALISED_FAMILY_HASH_SCORES] = normalised_scores
        bigger = True
        if self.abs_threshold is not None:
            if sum(s < self.abs_threshold for s in abs_sizes.values()) > 0:
                bigger |= False
        if self.rel_threshold is not None:
            if sum(r < self.rel_threshold for r in rel_sizes.values()) > 0:
                bigger |= False
        if bigger:
            self._distribute_results(bigger_sink=out_data)
        else:
            self._distribute_results(smaller_sink=out_data)


class NewFamilySifter(Sifter):
    def __init__(self, new_sink, other_sink, abs_threshold, rel_threshold):
        super().__init__(new_sink=new_sink, other_sink=other_sink)
        self.abs_threshold = abs_threshold
        self.rel_threshold = rel_threshold

    def sift(self, data):
        out_data = super().sift(data)
        search_size = data[SEARCH_SIZE]
        normalised_scores = data[NORMALISED_FAMILY_HASH_SCORES]

        if search_size < self.abs_threshold:
            self._distribute_results(other_sink=out_data)
        elif len([s for s in normalised_scores.values() if s > self.rel_threshold]) > 0:
            self._distribute_results(other_sink=out_data)
        else:
            self._distribute_results(new_sink=out_data)

class BiggerBetterSifter(Sifter):
    def __init__(self, bigger_threshold, better_threshold,
                 bigger_better_sink, smaller_worse_sink):
        super().__init__(bigger_better_sink=bigger_better_sink,
                         smaller_worse_sink=smaller_worse_sink)
        self.bigger_threshold = bigger_threshold
        self.better_threshold = better_threshold

    def sift(self, data):
        out_data = super().sift(data)
