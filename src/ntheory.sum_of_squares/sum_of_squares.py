from functools import lru_cache

from sympy.external.gmpy import bit_scan1, is_square, sqrtrem, jacobi, remove, sqrt as isqrt
from sympy.utilities.misc import as_int
from .factor_ import divisors, divisor_sigma, factorint
from .generate import nextprime
from .primetest import isprime

# def _comp(a,b,c,d):
#     return a*c + b*d, a*d - c*b

def _prime_as_sum_of_two_squares(p) -> tuple[int, int]:
    """
    Represent a prime `p` as a unique sum of two squares; this can
    only be done if the prime is congruent to 1 mod 4.

    Parameters
    ==========

    p : Integer
        A prime that is congruent to 1 mod 4

    Returns
    =======

    (int, int) : Pair of positive integers ``(x, y)`` satisfying ``x**2 + y**2 = p``.

    Examples
    ========

    >>> from sympy.solvers.diophantine.diophantine import prime_as_sum_of_two_squares
    >>> prime_as_sum_of_two_squares(5)
    (1, 2)

    Reference
    =========

    .. [1] Representing a number as a sum of four squares, [online],
           Available: https://schorn.ch/lagrange.html

    See Also
    ========

    sum_of_two_squares

    """
    if p % 8 == 5:
        # Legendre symbol (2/p) == -1 if p % 8 in [3, 5]
        b = 2
    elif p % 12 == 5:
        # Legendre symbol (3/p) == -1 if p % 12 in [5, 7]
        b = 3
    elif p % 5 in [2, 3]:
        # Legendre symbol (5/p) == -1 if p % 5 in [2, 3]
        b = 5
    else:
        b = 7
        while jacobi(b, p) == 1:
            b = nextprime(b)

    b = pow(b, p >> 2, p)
    a = p
    while b**2 > p:
        a, b = b, a % b
    return (int(a % b), int(b))  # convert from long


def sum_of_two_squares(n) -> tuple[int, int] | None:
    """

    Parameters
    ==========

    p : Integer
        non-negative integer

    Returns
    =======

    (int, int) | None : Pair of positive integers ``(x, y)`` satisfying ``x**2 + y**2 = n``.
                        a,b are sorted in ascending order. ``None`` if no such ``(a,b,c)``.

    Raises
    ======

    ValueError
        If ``n`` is a negative integer

    Examples
    ========

    Reference
    =========

    See Also
    ========

    sum_of_squares

    """
    n = as_int(n)
    if n < 0:
        raise ValueError()
    if n == 0:
        return (0, 0)
    n, r = remove(n, 2)
    t = 1 << (r>>1)
    s = t if r % 2 else 0
    for p, e in factorint(n).items():
        if e % 2 == 0:
            pe = p**(e>>1)
            s, t = t*pe, s*pe
        elif p % 4 == 1:
            a, b = _prime_as_sum_of_two_squares(p)
            for _ in range(e):
                s, t = s*a + t*b, s*b - t*a
        else:
            return None
    return tuple(sorted([abs(s), abs(t)]))


def sum_of_three_squares(n):
    r"""
    Returns a 3-tuple $(a, b, c)$ such that $a^2 + b^2 + c^2 = n$ and
    $a, b, c \geq 0$.

    Returns None if $n = 4^a(8m + 7)$ for some `a, m \in \mathbb{Z}`. See
    [1]_ for more details.

    Parameters
    ==========

    n : Integer
        non-negative integer

    Returns
    =======

    (int, int, int) | None : 3-tuple non-negative integers ``(a, b, c)`` satisfying ``a**2 + b**2 + c**2 = n``.
                             a,b,c are sorted in ascending order. ``None`` if no such ``(a,b,c)``.

    Raises
    ======

    ValueError
        If ``n`` is a negative integer

    Examples
    ========

    >>> from sympy.solvers.diophantine.diophantine import sum_of_three_squares
    >>> sum_of_three_squares(44542)
    (18, 37, 207)

    References
    ==========

    .. [1] Representing a number as a sum of three squares, [online],
        Available: https://schorn.ch/lagrange.html

    See Also
    ========

    sum_of_squares :
        ``sum_of_three_squares(n)`` is one of the solutions output by ``sum_of_squares(n, 3, zeros=True)``

    """
    # https://math.stackexchange.com/questions/483101/rabin-and-shallit-algorithm/651425#651425
    # discusses these numbers (except for 1, 2, 3) as the exceptions of H&L's conjecture that
    # Every sufficiently large number n is either a square or the sum of a prime and a square.
    special = {1: (0, 0, 1), 2: (0, 1, 1), 3: (1, 1, 1), 10: (0, 1, 3), 34: (3, 3, 4),
               58: (0, 3, 7), 85: (0, 6, 7), 130: (0, 3, 11), 214: (3, 6, 13), 226: (8, 9, 9),
               370: (8, 9, 15), 526: (6, 7, 21), 706: (15, 15, 16), 730: (0, 1, 27),
               1414: (6, 17, 33), 1906: (13, 21, 36), 2986: (21, 32, 39), 9634: (56, 57, 57)}
    n = as_int(n)
    if n < 0:
        raise ValueError("n should be a non-negative integer")
    if n == 0:
        return (0, 0, 0)
    n, v = remove(n, 4)
    v = 1 << v
    if n % 8 == 7:
        return
    if n in special:
        return tuple([v*i for i in special[n]])

    s, rem = sqrtrem(n)
    if not rem:
        return (0, 0, v*s)
    if n % 8 == 3:
        if not s % 2:
            s -= 1
        for x in range(s, -1, -2):
            N = (n - x**2) // 2
            if isprime(N):
                # n % 8 == 3 and x % 2 == 1 => N % 4 == 1
                y, z = _prime_as_sum_of_two_squares(N)
                return tuple(sorted([v*x, v*(y + z), v*abs(y - z)]))
        # We will never reach this point because there must be a solution.
        assert False

    # assert n % 4 in [1, 2]
    if not((n % 2) ^ (s % 2)):
        s -= 1
    for x in range(s, -1, -2):
        N = n - x**2
        if isprime(N):
            # assert N % 4 == 1
            y, z = _prime_as_sum_of_two_squares(N)
            return tuple(sorted([v*x, v*y, v*z]))
    # We will never reach this point because there must be a solution.
    assert False


def sum_of_four_squares(n):
    r"""
    Returns a 4-tuple `(a, b, c, d)` such that `a^2 + b^2 + c^2 + d^2 = n`.
    Here `a, b, c, d \geq 0`.

    Parameters
    ==========

    n : Integer
        non-negative integer

    Returns
    =======

    (int, int, int, int) : 4-tuple non-negative integers ``(a, b, c, d)`` satisfying ``a**2 + b**2 + c**2 + d**2 = n``.
                           a,b,c,d are sorted in ascending order.

    Raises
    ======

    ValueError
        If ``n`` is a negative integer

    Examples
    ========

    >>> from sympy.solvers.diophantine.diophantine import sum_of_four_squares
    >>> sum_of_four_squares(3456)
    (8, 8, 32, 48)
    >>> sum_of_four_squares(1294585930293)
    (0, 1234, 2161, 1137796)

    References
    ==========

    .. [1] Representing a number as a sum of four squares, [online],
        Available: https://schorn.ch/lagrange.html

    See Also
    ========

    sum_of_squares :
        ``sum_of_four_squares(n)`` is one of the solutions output by ``sum_of_squares(n, 4, zeros=True)``

    """
    n = as_int(n)
    if n < 0:
        raise ValueError("n should be a non-negative integer")
    if n == 0:
        return (0, 0, 0, 0)
    # remove factors of 4 since a solution in terms of 3 squares is
    # going to be returned; this is also done in sum_of_three_squares,
    # but it needs to be done here to select d
    n, v = remove(n, 4)
    v = 1 << v
    if n % 8 == 7:
        d = 2
        n = n - 4
    elif n % 8 in (2, 6):
        d = 1
        n = n - 1
    else:
        d = 0
    x, y, z = sum_of_three_squares(n)  # sorted
    return tuple(sorted([v*d, v*x, v*y, v*z]))


def sum_of_squares(n, k, zeros=False):
    pass


def sum_of_powers(n, p, k, zeros=False):
    pass


@lru_cache
def sum_of_squares_count(n, k):
    """
    https://en.wikipedia.org/wiki/Sum_of_squares_function
    https://mathworld.wolfram.com/SumofSquaresFunction.html
    A122141
    """
    n = as_int(n)
    k = as_int(k)
    if k <= 0:
        raise ValueError()
    if n < 0:
        raise ValueError()
    if n == 0:
        return 1
    if k == 1:
        return 2 if is_square(n) else 0
    if k == 2:
        n >>= bit_scan1(n)
        res = 1
        for p, e in factorint(n).items():
            if p % 4 == 1:
                res *= e + 1
            elif e % 2:
                return 0
        return res << 2
    if k == 4:
        n, r = remove(n, 2)
        return int((3 if r else 1)*divisor_sigma(n) << 3)
    if k == 6:
        n, r = remove(n, 2)
        s = (n % 4) - 2
        t = (1 << (2*r + 2)) + s
        return int(t*sum(s*d**2*((d % 4)-2) for d in divisors(n))) << 2
    if k == 8:
        if n % 2:
            return int(divisor_sigma(n, 3) << 4)
        return sum((-1 if d % 2 else 1)*d**3 for d in divisors(n)) << 4
    # otherwise
    z = sum_of_squares_count(n, k - 1)
    s = sum(sum_of_squares_count(n - i**2, k - 1) for i in range(1, isqrt(n) + 1))
    return int(z + (s << 1))


@lru_cache
def _rec_nontrivial_sum_of_squares_count(n, i, t):
    if n == 0:
        return 1 if t == 0 else 0
    if i < 1 or t < 1:
        return 0
    return _rec_nontrivial_sum_of_squares_count(n, min(i - 1, isqrt(n)), t) +\
           _rec_nontrivial_sum_of_squares_count(n - i**2, i, t - 1)


def nontrivial_sum_of_squares_count(n, k):
    """ A243148
    """
    n = as_int(n)
    k = as_int(k)
    if k <= 0:
        raise ValueError()
    if n < 0:
        raise ValueError()
    if n == 0:
        return 0
    if k == 1:
        return 1 if is_square(n) else 0
    if k == 2:
        n, r = remove(n, 2)
        B = 1
        for p, e in factorint(n).items():
            if p % 4 == 1:
                B *= e + 1
            elif e % 2:
                return 0
        if B % 2:
            B -= -1 if r % 2 else 1
        return B >> 1
    # otherwise
    return int(_rec_nontrivial_sum_of_squares_count(n, isqrt(n), k))

