"""Calculate Jaccard Index and Containment between sets."""


def _zero_on_divide_by_zero(func):
    def wrapped_function(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ZeroDivisionError:
            return 0
    return wrapped_function


@_zero_on_divide_by_zero
def jaccard(s, t):
    """
    Calculate the Jaccard Index of two sets.

    J(S, T) = |S ^ T| / |S v T|

    Parameters
    ----------
    s, t : iterable

    Returns
    -------
    float
    """
    s_, t_ = set(s), set(t)
    return len(s_ & t_) / len(s_ | t_)


@_zero_on_divide_by_zero
def jaccard_containment(s, t):
    """
    Calculate the Jaccard Containment of two sets.

    C(S, T) = |S ^ T| / |T|

    Parameters
    ----------
    s, t : iterable

    Returns
    -------
    float
    """
    s_, t_ = set(s), set(t)
    return len(s_ & t_) / len(s_)
