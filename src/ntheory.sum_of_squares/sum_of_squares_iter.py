def _rec(n, k, bottom):
    if k == 1:
        m, rem = sqrtrem(n)
        if rem == 0 and bottom <= m:
            yield [m]
        return
    for i in range(bottom, isqrt(n) + 1):
        for j in _rec(n - i**2, k - 1, i):
            yield [i] + j
