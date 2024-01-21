import math

##
from gmpy2 import kronecker, gcdext
from math import isqrt, gcd
##
from sympy.solvers.diophantine.diophantine import diop_DN
from sympy.utilities.misc import as_int
from sympy.ntheory.residue_ntheory import sqrt_mod


def _class_number_neg_naive(disc):
    r""" class number of quadratic field `\mathbb{Q}(\sqrt{d})`.
    ``disc`` は `\mathbb{Q}(\sqrt{d})` の判別式とし、負であることを仮定する。
    類数公式

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

    # def get_order(self, a, b, c):
    #     odd = self.disc % 2
    #     dic = {(1, odd, odd - self.disc): 0, (a, abs(b), c): 1}
    #     ra, rb, rc = self.composition_double(a, b, c)
    #     for k in range(2, 100):
    #         rb_abs = abs(rb)
    #         if (ra, rb_abs, rc) in dic:
    #             pass
    #         dic[(ra, rb_abs, rc)] = k
    #         ra, rb, rc = self.composition(ra, rb, rc, a, b, c)

    # def counting_order(self):
    #     cnt = 0
    #     odd = self.disc % 2
    #     identity = (1, odd, odd - self.disc)
    #     upper = 2 * abs(self.disc) // 3
    #     m = 1
    #     for a, b, c in self.enumerate_elements():
    #         cnt += 1
    #         if b < 0 or self.is_identity(a, b, c):
    #             continue
    #         if m > 1:
    #             a, b, c = self.n_times(a, b, c, m)
    #             if self.is_identity(a, b, c):
    #                 continue
    #         self.get_order(a, b, c)
    #     return cnt



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


def class_number_of_quadratic_field(*, d=None, disc=None):
    r""" class number of quadratic field `\mathbb{Q}(\sqrt{d})`.

    `\mathbb{Q}(\sqrt{d})`の判別式``disc``は次のように定義される。

    .. math ::
        disc(d) =
        \begin{cases}
            d, & \mbox{if } d \equiv 1 \pmod{4}\\
            4d, & \mbox{otherwise}
        \end{cases}

    この関数は、``d``, ``disc``のどちらか一方を与えれば、類数を計算できる。

    dは、0でも1でもない整数で、平方因子を持たないとする。
    よって、``disc``は次のどちらかを満たさなければならない。

    * ``disc % 4 == 1``かつ``disc``は平方因子を持たない。
    * ``disc % 16 in [8, 12]``かつ``disc//4``は平方因子を持たない。

    ただし、前提を満たさない``d``あるいは``disc``を入力しても、エラーは発生しない。

    Parameters
    ==========

    d : 0でも1でもない整数で、平方因子を持たない。
    disc : `\mathbb{Q}(\sqrt{d})` の 判別式。

    Returns
    =======

    int

    Raises
    ======

    ValueError

    Examples
    ========

    類数が1となるnegative number dは、次に挙げる有限個しかない。
    これらは、Gauss numbers あるいは Heegner numbers と呼ばれる。

    >>> from ****
    >>> ds = [1, 2, 3, 7, 11, 19, 43, 67, 163] # A003173
    >>> all(class_number_of_quadratic_field(d=-d) == 1 for d in ds)
    True

    同じことだが、``disc``で指定する場合は次のようになる

    >>> discs = [3, 4, 7, 8, 11, 19, 43, 67, 163] # A014602
    >>> all(class_number_of_quadratic_field(disc=-disc) == 1 for disc in discs)
    True

    References
    ==========

    .. [1] ****

    """
    if not((d is None) ^ (disc is None)):
        raise ValueError("dとdiscどちから一方のみを指定してください")
    if disc is None:
        d = as_int(d)
        disc = d if d % 4 == 1 else 4*d
    else:
        disc = as_int(disc)
    if disc == 0:
        raise ValueError()
    if disc < 0:
        # 虚2次体
        if -4 <= disc:
            return 1
        if -1000 <= disc:
            return _class_number_neg_naive(disc)
        return QuadraticFormClassGroup(disc).counting_naive()
    else:
        # 実2次体
        if 30_000 < disc:
            raise ValueError("大きすぎなので計算できません")
        return _class_number_pos_naive(disc)



# import time
# time_sta = time.time()
# for _ in range(100):
#     len(QuadraticFormClassGroup(-800))
# time_end = time.time()

# print(time_end- time_sta)


# time_sta = time.time()
# for _ in range(100):
#     _class_number_neg_naive(-800)
# time_end = time.time()

# print(time_end- time_sta)

# qf = QuadraticFormClassGroup(-800)
# for a,b,c in qf.enumerate_elements():
#     if qf.is_identity(a,b,c):
#         continue
#     print(a,b,c)
#     for k in range(1, 12):
#         print(k, qf.n_times(a,b,c,k))

# qf = QuadraticFormClassGroup(-63)
# for k in range(1, 5):
#     print(k, qf.n_times(2, 1, 8, k))

print(_class_number_neg_naive(-27))