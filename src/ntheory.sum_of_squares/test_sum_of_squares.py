from sympy.core.random import randint

from sympy.ntheory.sum_of_squares import (_prime_as_sum_of_two_squares,
                                          sum_of_two_squares, sum_of_three_squares, sum_of_four_squares,
                                          sum_of_squares, sum_of_powers,
                                          sum_of_squares_count,
                                          nontrivial_sum_of_squares_count)

from sympy.testing.pytest import raises


def test_prime_as_sum_of_two_squares():
    for i in [5, 13, 17, 29, 37, 41, 2341, 3557, 34841, 64601]:
        a, b = _prime_as_sum_of_two_squares(i)
        assert a**2 + b**2 == i
    ans = _prime_as_sum_of_two_squares(800029)
    assert ans == (450, 773) and type(ans[0]) is int


def test_sum_of_two_squares():
    raises(ValueError, lambda: sum_of_two_squares(-1))

    A001481 = [0, 1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20, 25,
               26, 29, 32, 34, 36, 37, 40, 41, 45, 49, 50, 52,
               53, 58, 61, 64, 65, 68, 72, 73, 74, 80, 81, 82,
               85, 89, 90, 97, 98]
    for n in range(100):
        ret = sum_of_two_squares(n)
        if ret is None:
            assert n not in A001481
        else:
            x, y = ret
            assert x**2 + y**2 == n


def test_sum_of_three_squares():
    for i in [0, 1, 2, 34, 123, 34304595905, 34304595905394941, 343045959052344,
              800, 801, 802, 803, 804, 805, 806]:
        a, b, c = sum_of_three_squares(i)
        assert a**2 + b**2 + c**2 == i
        assert a >= 0

    # error
    raises(ValueError, lambda: sum_of_three_squares(-1))

    assert sum_of_three_squares(7) is None
    assert sum_of_three_squares((4**5)*15) is None
    # if there are two zeros, there might be a solution
    # with only one zero, e.g. 25 => (0, 3, 4) or
    # with no zeros, e.g. 49 => (2, 3, 6)
    assert sum_of_three_squares(25) == (0, 0, 5)
    assert sum_of_three_squares(4) == (0, 0, 2)


def test_sum_of_four_squares():
    # this should never fail
    n = randint(1, 100000000000000)
    assert sum(i**2 for i in sum_of_four_squares(n)) == n

    # error
    raises(ValueError, lambda: sum_of_four_squares(-1))

    for n in range(1000):
        result = sum_of_four_squares(n)
        assert len(result) == 4
        assert all(r >= 0 for r in result)
        assert sum(r**2 for r in result) == n
        assert list(result) == sorted(result)


def test_sum_of_squares():
    pass


def test_sum_of_powers():
    pass


def test_sum_of_squares_count():
    raises(ValueError, lambda: sum_of_squares_count(-1, 4))
    raises(ValueError, lambda: sum_of_squares_count(0, 0))

    A000122 = [1, 2, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
               0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for n in range(len(A000122)):
        assert sum_of_squares_count(n, 1) == A000122[n]

    A004018 = [1, 4, 4, 0, 4, 8, 0, 0, 4, 4, 8, 0, 0, 8, 0,
               0, 4, 8, 4, 0, 8, 0, 0, 0, 0, 12, 8, 0, 0, 8,
               0, 0, 4, 0, 8, 0, 4, 8, 0, 0, 8, 8, 0, 0, 0,
               8, 0, 0, 0, 4, 12, 0, 8, 8, 0, 0, 0, 0, 8, 0]
    for n in range(len(A004018)):
        assert sum_of_squares_count(n, 2) == A004018[n]

    A000118 = [1, 8, 24, 32, 24, 48, 96, 64, 24, 104, 144, 96,
               96, 112, 192, 192, 24, 144, 312, 160, 144, 256,
               288, 192, 96, 248, 336, 320, 192, 240, 576, 256,
               24, 384, 432, 384, 312, 304, 480, 448, 144, 336]
    for n in range(len(A000118)):
        assert sum_of_squares_count(n, 4) == A000118[n]

    A000141 = [1, 12, 60, 160, 252, 312, 544, 960, 1020, 876, 1560,
               2400, 2080, 2040, 3264, 4160, 4092, 3480, 4380, 7200,
               6552, 4608, 8160, 10560, 8224, 7812, 10200, 13120,
               12480, 10104, 14144, 19200, 16380, 11520, 17400, 24960]
    for n in range(len(A000141)):
        assert sum_of_squares_count(n, 6) == A000141[n]

    A000143 = [1, 16, 112, 448, 1136, 2016, 3136, 5504, 9328, 12112,
               14112, 21312, 31808, 35168, 38528, 56448, 74864, 78624,
               84784, 109760, 143136, 154112, 149184, 194688, 261184,
               252016, 246176, 327040, 390784, 390240, 395136, 476672]
    for n in range(len(A000143)):
        assert sum_of_squares_count(n, 8) == A000143[n]

    A066535 = [1, 2, 4, 8, 24, 112, 544, 2368, 9328, 34802, 129064,
               491768, 1938336, 7801744, 31553344, 127083328, 509145568,
               2035437440, 8148505828, 32728127192, 131880275664,
               532597541344, 2153312518240, 8710505815360, 35250721087168,
               142743029326162, 578472382307304]
    for n in range(1, len(A066535)):
        assert sum_of_squares_count(n, n) == A066535[n]


def test_nontrivial_sum_of_squares_count():
    raises(ValueError, lambda: nontrivial_sum_of_squares_count(-1, 4))
    raises(ValueError, lambda: nontrivial_sum_of_squares_count(0, 0))

    A010052 = [1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
               0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    for n in range(1, len(A010052)):
        assert nontrivial_sum_of_squares_count(n, 1) == A010052[n]

    A025426 = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1,
               0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1,
               0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 1, 1, 0, 0, 0,
               0, 1, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0]
    for n in range(len(A025426)):
        assert nontrivial_sum_of_squares_count(n, 2) == A025426[n]

    A025427 = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1,
               1, 0, 1, 1, 0, 1, 0, 1, 2, 0, 1, 1, 0, 0, 2, 1, 1, 1, 0,
               2, 0, 0, 2, 1, 1, 1, 1, 1, 0, 1, 1, 1, 2, 0, 1, 3, 0, 1,
               2, 0, 2, 0, 1, 2, 0, 0, 1, 3, 1, 1, 2, 1, 0, 1, 1, 2, 2]
    for n in range(len(A025427)):
        assert nontrivial_sum_of_squares_count(n, 3) == A025427[n]

    A025428 = [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1,
               1, 1, 1, 1, 1, 0, 1, 1, 1, 3, 0, 1, 2, 0, 1, 2, 1, 2, 2,
               1, 2, 1, 0, 3, 2, 1, 2, 1, 2, 1, 2, 2, 1, 4, 1, 2, 3, 0,
               2, 4, 1, 3, 2, 1, 4, 1, 1, 3, 3, 2, 2, 4, 2, 1, 3, 2, 3]
    for n in range(len(A025428)):
        assert nontrivial_sum_of_squares_count(n, 4) == A025428[n]

