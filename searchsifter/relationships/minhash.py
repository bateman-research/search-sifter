"""Generate and compare MinHash signatures."""
import binascii
from functools import reduce
import operator
from .jaccard import _zero_on_divide_by_zero


def _crc32_hash(e):
    return binascii.crc32(str(e).encode())


def signature(s, n, hash_function=_crc32_hash):
    """
    Calculate the set signature of an iterable.

    Parameters
    ----------
    s : iterable
    n : int
        Signature length.
    hash_function : function
        The hash function to use. By default, CRC32.

    Returns
    -------
    set of int
        The signature as a set.

    """
    return set(sorted((hash_function(e), e) for e in s)[:n])


def set_hashes(s, hash_function=_crc32_hash):
    """
    Compute and return hashes for every element of an iterable.

    Parameters
    ----------
    s : iterable
    hash_function : function
        The hash function to use. By default, CRC32.

    Returns
    -------
    set of int
    """
    return {(hash_function(e), e) for e in s}


def union_signature(*args):
    """
    Calculate the signature of the union of several sets from their signatures.

    Parameters
    ----------
    s, t... : iterable of int
    n : int
        Signature length.

    Returns
    -------
    set of int
        The signature of the union.
    """
    *sets, n = args
    if len(sets) == 2:
        return set(sorted(sets[0] | sets[1])[:n])
    else:
        return set(sorted(reduce(operator.or_, list(sets)))[:n])


def intersection_signature(*sets):
    """
    Calculate the signature of the intersection of two sets from their
    signatures.

    Paramaters
    ----------
    s, t,... : iterable of int

    Returns
    -------
    set of int
        The signature of the union.
    """
    return reduce(operator.and_, sets)


@_zero_on_divide_by_zero
def minhash(s, t, n):
    """
    Calculate the MinHash estimate of the Jaccard Index.

    Parameters
    ----------
    s, t : iterable of int
        Set signatures, as returned by `signature`.
    n : int
        The length of signature.

    Returns
    -------
    float
    """
    X = union_signature(s, t, n)
    Y = intersection_signature(X, s, t)
    return len(Y) / len(X)


@_zero_on_divide_by_zero
def minhash_containment(s, t):
    """
    Calculate the MinHash estimate of the Jaccard Containment.

    Parameters
    ----------
    s, t : iterable of int
        Set signatures, as returned by `signature`.
    n : int
        The length of signature.

    Returns
    -------
    float
    """
    return len(s & t) / len(s)
