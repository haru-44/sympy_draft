import math

from sympy.external.gmpy import kronecker
from sympy.solvers.diophantine.diophantine import diop_DN
from sympy.utilities.misc import as_int


def _class_number_neg_naive(disc):
    r""" class number of quadratic field `\mathbb{Q}(\sqrt{d})`.
    Let ``disc`` be a discriminant of `\mathbb{Q}(\sqrt{d})` and assume it is negative.

    .. math ::
        h(d) = \frac{1}{w} \sum{n=1}^{|D|/2} \left( \frac{D}{n} \right)

    where `D` = ``disc``.

    References
    ==========

    .. [1] ****

    """
    h = sum(kronecker(disc, n) for n in range(1, abs(disc) // 2 + 1))
    if disc % 8 == 5:
        h //= 3
    elif disc % 4 == 0:
        h //= 2
    return int(h)


def class_number_of_quadratic_field(*, d=None, disc=None):
    r""" Calculate class number of quadratic field `\mathbb{Q}(\sqrt{d})`.

    The discriminant ``disc`` of ``mathbb{Q}(\sqrt{d})`` is defined as follows [1]_.

    .. math ::
        disc(d) =
        \begin{cases}
            d, & \mbox{if } d \equiv 1 \pmod{4}\\
            4d, & \mbox{if } d \equiv 2,3 \pmod{4}
        \end{cases}

    This function will compute the class number if either ``d`` or ``disc`` is given.
    Let ``d`` be an integer that is neither 0 nor 1 and square-free integer.
    Therefore, ``disc`` must satisfy either of the following:

    * ``disc % 4 == 1`` and ``disc`` is a square-free integer.
    * ``disc % 16 in [8, 12]`` and ``disc//4`` is a square-free integer.

    However, inputting ``d`` or ``disc`` that does not satisfy the assumption will not result in an error.

    According to the Dirichlet class number formula [2]_, the class number `h(d)` can be calculated as follows:

    .. math ::
        h(d) =
        \begin{cases}
            \frac{w \sqrt{|d|}}{2 \pi} L(1, \chi), & \mbox{if } d < 0\\
            \frac{\sqrt{d}}{\ln \varepsilon} L(1, \chi), & \mbox{if } d > 0
        \end{cases}

    Parameters
    ==========

    d : Integer
        ``d`` be an integer that is neither 0 nor 1 and square-free integer.
    disc : Integer
        discriminant of `\mathbb{Q}(\sqrt{d})`

    Returns
    =======

    int : class number of quadratic field `\mathbb{Q}(\sqrt{d})`

    Examples
    ========

    There are only a finite number of negative numbers ``d`` whose class number is 1, as listed below.
    These are called Gauss numbers or Heegner numbers.

    >>> from ****
    >>> ds = [1, 2, 3, 7, 11, 19, 43, 67, 163] # A003173
    >>> all(class_number_of_quadratic_field(d=-d) == 1 for d in ds)
    True

    The same is true for ``disc`` as follows:

    >>> discs = [3, 4, 7, 8, 11, 19, 43, 67, 163] # A014602
    >>> all(class_number_of_quadratic_field(disc=-disc) == 1 for disc in discs)
    True

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Quadratic_field
    .. [2] https://en.wikipedia.org/wiki/Class_number_formula
    .. [3] https://oeis.org/A003173
    .. [4] https://oeis.org/A014602

    """
    if not((d is None) ^ (disc is None)):
        raise ValueError("Please specify only one of d or disc")
    if disc is None:
        d = as_int(d)
        disc = d if d % 4 == 1 else 4*d
    else:
        disc = as_int(disc)
    if disc == 0:
        raise ValueError("disc must be non-zero")
    if disc < 0:
        # imaginary quadratic field
        if -4 <= disc:
            return 1
        return _class_number_neg_naive(disc)
    else:
        # real quadratic field
        raise ValueError("Not yet implemented")

#######################


def _class_number_pos_naive(disc):
    r""" class number of quadratic field `\mathbb{Q}(\sqrt{d})`.
    ``disc`` は `\mathbb{Q}(\sqrt{d})` の判別式とし、正であることを仮定する。

    References
    ==========

    .. [1] ****

    """
    x, y = sorted((x, y) for x, y in diop_DN(disc, 4) + diop_DN(disc, -4) if 0 < y and 0 < x)[0]
    # 基本単数
    eps_0 = (x + y * math.sqrt(disc)) / 2
    h = 1
    for a in range(1, disc):
        k = kronecker(disc, a)
        if k == 1:
            h /= math.sin(a * math.pi / disc)
        elif k == -1:
            h *= math.sin(a * math.pi / disc)
    return round(math.log(h, eps_0) / 2)

class QuadraticFormClassGroup:
    def __init__(self, disc):
        self.disc = disc

    def reduction(self, a, b, c):
        """ 簡約形式を返す
        """
        while True:
            if c < a:
                a, b, c = c, -b, a
            if a < b or b <= -a:
                b = b % (2 * a)
                if a < b:
                    b -= 2 * a
                c = (b**2 - self.disc) // (4 * a)
            elif a == c and -a < b < 0:
                return a, -b, c
            else:
                return a, b, c

    def composition(self, a1, b1, c1, a2, b2, c2):
        """ 合成
        """
        h, _, v = gcdext(a1, a2)
        g, u, w = gcdext(h, (b1 + b2) // 2)
        a3 = (a1 * a2) // g**2
        b3 = b2 + 2 * a2 * ((b1 - b2) // 2 * v * u - c2 * w) // g
        c3 = (b3**2 - self.disc) // (4 * a3)
        return self.reduction(a3, b3, c3)

    def composition_double(self, a, b, c):
        g, _, w = gcdext(a, b)
        k = a // g
        a_ = k**2
        b_ = b - 2 * k * c * w
        c_ = (b_**2 - self.disc) // (4 * a_)
        return self.reduction(a_, b_, c_)

    def n_times(self, a, b, c, n):
        ra, rb, rc = a, b, c
        n -= 1
        while 0 < n:
            if n % 2 == 1:
                ra, rb, rc = self.composition(ra, rb, rc, a, b, c)
            a, b, c = self.composition_double(a, b, c)
            n >>= 1
        return ra, rb, rc

    def enumerate_elements(self):
        """ 枚挙
        """
        for a in range(1, isqrt(-self.disc // 3) + 1):
            for b in filter(lambda x: x <= a, sqrt_mod(self.disc, 4 * a, all_roots=True)):
                c = (b**2 - self.disc) // (4 * a)
                if c < a or gcd(a, b, c) != 1:
                    continue
                yield a, b, c
                if a != b and a != c and b != 0:
                    yield a, -b, c

    def is_identity(self, a, b, c):
        if a != 1:
            return False
        odd = self.disc % 2
        return b == odd and 4 * c == odd - self.disc

    def counting_naive(self):
        return sum(1 for _ in self.enumerate_elements())


def test_class_number_of_quadratic_field():
    # imaginary
    ds = [1, 2, 3, 7, 11, 19, 43, 67, 163] # A003173
    assert all(class_number_of_quadratic_field(d=-d) == 1 for d in ds)

    discs_list = [(1, [3, 4, 7, 8, 11, 19, 43, 67, 163]), # A014602
                  (3, [23, 31, 59, 83, 107, 139, 211, 283,
                       307, 331, 379, 499, 547, 643, 883,
                       907]), # A006203
                  (23, [647, 1039, 1103, 1279, 1447, 1471,
                        1811, 1979, 2411, 2671, 3491, 3539,
                        3847, 3923, 4211, 4783, 5387, 5507,
                        5531, 6563, 6659, 6703, 7043, 9587,
                        9931, 10867, 10883, 12203, 12739,
                        13099, 13187, 15307, 15451, 16267,
                        17203, 17851, 18379, 20323]), # A046020
                 ]
    for h, discs in discs_list:
        assert all(class_number_of_quadratic_field(disc=-disc) == h for disc in discs)

    # real
    ds_list = [(1, [2, 3, 5, 6, 7, 11, 13, 14, 17, 19,
                    21, 22, 23, 29, 31, 33, 37]), # A003172
               (2, [10, 15, 26, 30, 34, 35, 39, 42, 51]), # A029702
               (3, [79, 142, 223, 229, 254, 257, 321, 326,
                    359, 443, 469, 473, 659, 733, 761, 839]), # A029703
               (4, [82, 130, 145, 170, 195, 210, 219, 231,
                     255, 274, 290, 291, 322, 323, 330, 370,
                     390, 410, 434, 435]), # A029704
               (5, [401, 439, 499, 727, 817, 982, 1093, 1126,
                    1327, 1393, 1429, 1486, 1641, 1766, 1897,
                    2027, 2081, 2153]), # A029705
               (6, [235, 346, 427, 506, 574, 697, 785, 786,
                    842, 874, 894, 895, 898, 899, 906, 985,
                    1086, 1191, 1211]), # A218038
               (7, [577, 1009, 1087, 1294, 1601, 1761, 1934,
                    2029, 2251, 2302, 2467, 2913, 4139, 4229,
                    4702, 5039, 5273]), # A218039
               (8, [226, 399, 442, 646, 799, 870, 910, 994,
                    1023, 1122, 1155, 1239, 1290, 1299, 1351,
                    1443, 1446, 1590]), # A218040
               (9, [1129, 1654, 3137, 3719, 4409, 4534, 5521,
                    5623, 5878, 6809, 7573, 7873]), # A218041
               (10, [1111, 1226, 2031, 2335, 2362, 2602, 2986,
                     3129, 3246, 3379, 3585, 3598, 3599, 3722,
                     3782, 3966, 4097]) # A218042
               ]
    for h, ds in ds_list:
        assert all(class_number_of_quadratic_field(d=d) == h for d in ds)

