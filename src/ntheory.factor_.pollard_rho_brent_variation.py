# factorintとの接続をどうするか

def pollard_rho_brent_variation(n, s=2, a=1, max_steps=None):
    r""" Use Pollard's rho method (Brent Variation) to find a nontrivial factor of ``n``.

    Explanation
    ===========

    The original Pollard's rho method uses Floyd's cycle-finding algorithm,
    which is replaced by Brent's cycle-finding algorithm to speed up the process.
    In fact, the number of times the pseudo-random function is computed is less than in the original rho method.
    Moreover, by computing the gcd every ``INTERVAL_GCD`` times instead of every time, further speed-up is expected.

    Parameters
    ==========

    n : Integer
        Positive integer to be factored
    s : int
        initial value
    a : int
        Using ``(x**2+a) % n`` as a pseudo-random function
    max_steps : int | None
        Maximum number of steps to search. ``n`` is set if ``None``.
        The number of steps is expected to be `O(n^{1/4})`, so it should be found by ``n``.

    Returns
    =======

    int | None : A nontrivial divisor of ``n``. ``None`` if not found.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import pollard_rho_brent_variation
    >>> pollard_rho_brent_variation(1099564581221)
    524309
    >>> pollard_rho_brent_variation(1099564581221, max_steps=50) is None
    True

    References
    ==========

    .. [1] Richard P. Brent, An improved Monte Carlo factorization algorithm (1980),
           BIT, Volume 20, Issue 2, pp. 176-184,
           https://doi.org/10.1007/BF01933190
           https://maths-people.anu.edu.au/~brent/pub/pub051.html

    """
    n = int(n)
    if n < 5:
        raise ValueError('pollard_rho_brent_variation should receive n > 4')
    s %= n
    a %= n
    reach = q = 1
    INTERVAL_GCD = 32
    if max_steps is None:
        max_steps = n
    while True:
        x = s
        # Skip to `reach`
        for _ in range(reach):
            s = pow(s, 2, n) + a
        for k in range(0, reach, INTERVAL_GCD):
            store_s = s
            for _ in range(min(INTERVAL_GCD, reach - k)):
                s = pow(s, 2, n) + a
                q = (q * (x - s)) % n
            g = gcd(q, n)
            if 1 < g:
                if g < n:
                    # Found a nontrivial divisor
                    return g
                # Start over from the saved s
                while True:
                    store_s = pow(store_s, 2, n) + a
                    g = gcd(x - store_s, n)
                    if 1 < g:
                        return g if g < n else None
            if max_steps <= k:
                return None
        reach <<= 1

