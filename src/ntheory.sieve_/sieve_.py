from array import array as _array
from sympy.utilities.misc import as_int
from .factor_ import trailing, multiplicity


class SieveBase:
    """ Basic class for sieve
    """

    def __init__(self, initial_list, n):
        self._sieve_list = initial_list
        self._initial_length = len(self._sieve_list)
        if n is not None:
            self.extend(n)

    def reset(self):
        """ reset the sieve
        """
        self._sieve_list = self._sieve_list[:self._initial_length]


class FactorSieve(SieveBase):
    """ Prime factor sieve

    Examples
    ========

    >>> from sympy.ntheory.sieve_ import FactorSieve
    >>> fsieve = FactorSieve(100) # Creates the sieve up to 100
    >>> fsieve[35] # Returns the smallest prime factor
    5
    >>> fsieve[7] # Returns itself if prime
    7
    >>> fsieve[143] # The sieve is automatically expanded
    11

    """

    def __init__(self, n=None):
        super().__init__(_array('L', [1, 3, 5]), n)

    def extend(self, n):
        """ Extend the sieve to n

        Parameters
        ==========

        n : positive integer

        Raises
        ======

        ValueError
            If ``2**32 < n`` or ``n`` is not an integer.

        """
        n = as_int(n)
        begin = len(self._sieve_list)
        n >>= 1
        if n < begin:
            return
        self._sieve_list += _array('L', range(2 * begin + 1, 2 * n, 2))
        for i in range(1, begin):
            val = self._sieve_list[i]
            i2 = 2 * i + 1
            k = (begin + i) // i2
            for j in range(k * i2 + i, n, i2):
                self._sieve_list[j] = min(val, self._sieve_list[j])
        for i in range(begin, (n + 1) // 3):
            val = self._sieve_list[i]
            i2 = 2 * i + 1
            for j in range(i2 + i, n, i2):
                self._sieve_list[j] = min(val, self._sieve_list[j])

    def __getitem__(self, n):
        n = as_int(n)
        if n < 1:
            raise ValueError("n must be a positive integer")
        if n % 2 == 0:
            return 2
        if len(self._sieve_list) <= n >> 1:
            self.extend(n + 1)
        return self._sieve_list[n >> 1]

    def factorint(self, n):
        """ Returns the prime factorization of ``n``

        Examples
        ========

        >>> from sympy.ntheory.sieve_ import FactorSieve
        >>> fsieve = FactorSieve()
        >>> fsieve.factorint(100)
        {2: 2, 5: 2}

        If ``n`` is negative, the result is
        the same as when ``n`` is positive.

        >>> fsieve.factorint(64) == fsieve.factorint(-64)
        True

        If ``n`` is ``0`` or ``1``, an empty set is returned.

        >>> fsieve.factorint(0) == fsieve.factorint(1) == {}
        True

        See Also
        ========

        sympy.ntheory.factor_.factorint

        """
        factors = {}
        n = as_int(n)
        t = trailing(n)
        if t:
            n >>= t
            factors[2] = t
        while n > 1:
            p = self[n]
            t = multiplicity(p, n)
            factors[p] = t
            n //= p**t
        return factors


class TotientSieve(SieveBase):
    """ Totient function sieve

    Examples
    ========

    >>> from sympy.ntheory.sieve_ import TotientSieve
    >>> tsieve = TotientSieve(100) # Creates the sieve up to 100
    >>> tsieve[35]
    24
    >>> tsieve[143] # The sieve is automatically expanded
    120

    See Also
    ========

    sympy.ntheory.factor_.totient

    """

    def __init__(self, n=None):
        super().__init__(_array('L', [1, 2, 4]), n)

    def extend(self, n):
        """ Extend the sieve to n

        Parameters
        ==========

        n : positive integer

        Raises
        ======

        ValueError
            If ``2**32 < n`` or ``n`` is not an integer.

        """
        n = as_int(n)
        begin = len(self._sieve_list)
        n >>= 1
        if n < begin:
            return
        self._sieve_list += _array('L', range(2 * begin + 1, 2 * n, 2))
        for i in range(begin):
            val = self._sieve_list[i]
            i2 = 2 * i + 1
            k = (begin + i) // i2
            for j in range(k * i2 + i, n, i2):
                self._sieve_list[j] -= val
        for i in range(begin, (n + 1) // 3):
            val = self._sieve_list[i]
            i2 = 2 * i + 1
            for j in range(i2 + i, n, i2):
                self._sieve_list[j] -= val

    def __getitem__(self, n):
        n = as_int(n)
        if n < 1:
            raise ValueError("n must be a positive integer")
        if len(self._sieve_list) <= n >> 1:
            self.extend(n + 1)
        t = trailing(n)
        if t:
            return self._sieve_list[n >> (t + 1)] << (t - 1)
        return self._sieve_list[n >> 1]


class MobiusSieve(SieveBase):
    """ Mobius function sieve

    Examples
    ========

    >>> from sympy.ntheory.sieve_ import MobiusSieve
    >>> msieve = MobiusSieve(100) # Creates the sieve up to 100
    >>> msieve[35]
    1
    >>> msieve[105] # The sieve is automatically expanded
    -1

    See Also
    ========

    sympy.ntheory.residue_ntheory.mobius

    """
    def __init__(self, n=None):
        super().__init__(_array('i', [1, -1, -1]), n)

    def extend(self, n):
        """ Extend the sieve to n

        Parameters
        ==========

        n : positive integer

        Raises
        ======

        ValueError
            If ``n`` is not an integer.

        """
        n = as_int(n)
        begin = len(self._sieve_list)
        n >>= 1
        if n < begin:
            return
        self._sieve_list += _array('i', [0] * (n - begin))
        for i in range(begin):
            val = self._sieve_list[i]
            i2 = 2 * i + 1
            k = (begin + i) // i2
            for j in range(k * i2 + i, n, i2):
                self._sieve_list[j] -= val
        for i in range(begin, (n + 1) // 3):
            val = self._sieve_list[i]
            i2 = 2 * i + 1
            for j in range(i2 + i, n, i2):
                self._sieve_list[j] -= val

    def __getitem__(self, n):
        n = as_int(n)
        if n < 1:
            raise ValueError("n must be a positive integer")
        if len(self._sieve_list) <= n >> 1:
            self.extend(n + 1)
        t = trailing(n)
        if t > 1:
            return 0
        return pow(-1, t % 2) * self._sieve_list[n >> (t + 1)]
