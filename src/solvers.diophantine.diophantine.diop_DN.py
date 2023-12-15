# cornacchia の挙動を変える

def diop_DN(D, N, t=symbols("t", integer=True)):
    """
    Solves the equation `x^2 - Dy^2 = N`.

    Explanation
    ===========

    Mainly concerned with the case `D > 0, D` is not a perfect square,
    which is the same as the generalized Pell equation. The LMM
    algorithm [1]_ is used to solve this equation.

    Returns one solution tuple, (`x, y)` for each class of the solutions.
    Other solutions of the class can be constructed according to the
    values of ``D`` and ``N``.

    Usage
    =====

    ``diop_DN(D, N, t)``: D and N are integers as in `x^2 - Dy^2 = N` and
    ``t`` is the parameter to be used in the solutions.

    Details
    =======

    ``D`` and ``N`` correspond to D and N in the equation.
    ``t`` is the parameter to be used in the solutions.

    Examples
    ========

    >>> from sympy.solvers.diophantine.diophantine import diop_DN
    >>> diop_DN(13, -4) # Solves equation x**2 - 13*y**2 = -4
    [(3, 1), (393, 109), (36, 10)]

    The output can be interpreted as follows: There are three fundamental
    solutions to the equation `x^2 - 13y^2 = -4` given by (3, 1), (393, 109)
    and (36, 10). Each tuple is in the form (x, y), i.e. solution (3, 1) means
    that `x = 3` and `y = 1`.

    >>> diop_DN(986, 1) # Solves equation x**2 - 986*y**2 = 1
    [(49299, 1570)]

    See Also
    ========

    find_DN(), diop_bf_DN()

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
        Robertson, July 31, 2004, Pages 16 - 17. [online], Available:
        https://web.archive.org/web/20160323033128/http://www.jpr2718.org/pell.pdf
    """
    sol = []
    if D < 0:
        if N == 0:
            return [(0, 0)]
        if N < 0:
            return sol
        # N > 0:
        for d in divisors(square_factor(N)):
            for x, y in cornacchia(1, -D, N // d**2):
                sol.append((d*x, d*y))
                if D == -1:
                    sol.append((d*y, d*x))
        return sol

    if D == 0:
        if N == 0:
            return [(0, t)]
        if N < 0:
            return sol
        sN, _exact = integer_nthroot(N, 2)
        if _exact:
            return [(sN, t)]
        return sol

    # D > 0
    sD, _exact = integer_nthroot(D, 2)
    if _exact:
        if N == 0:
            return [(sD*t, t)]
        for y in range(floor(sign(N)*(N - 1)/(2*sD)) + 1):
            if (_N := D*y**2 + N) > 0:
                sq, _exact = integer_nthroot(_N, 2)
                if _exact:
                    sol.append((sq, y))
        return sol

    if 1 < N**2 < D:
        # It is much faster to call `_special_diop_DN`.
        return _special_diop_DN(D, N)

    if N == 0:
        return [(0, 0)]

    if abs(N) == 1:
        pqa = PQa(0, 1, D)
        *_, prev_B, prev_G = next(pqa)
        for j, (*_, a, _, _B, _G) in enumerate(pqa):
            if a == 2*sD:
                break
            prev_B, prev_G = _B, _G
        if j % 2:
            if N == 1:
                sol.append((prev_G, prev_B))
            return sol
        if N == -1:
            return [(prev_G, prev_B)]
        for _ in range(j):
            *_, _B, _G = next(pqa)
        return [(_G, _B)]

    for f in divisors(N):
        m, r = divmod(N, f**2)
        if r:
            continue
        am = abs(m)
        if am == 2:
            continue
        for sqm in sqrt_mod(D, am, all_roots=True):
            z = symmetric_residue(sqm, am)
            pqa = PQa(z, am, D)
            *_, prev_B, prev_G = next(pqa)
            for _ in range(length(z, am, D) - 1):
                _, q, *_, _B, _G = next(pqa)
                if abs(q) == 1:
                    if prev_G**2 - D*prev_B**2 == m:
                        sol.append((f*prev_G, f*prev_B))
                    elif a := diop_DN(D, -1):
                        sol.append((f*(prev_G*a[0][0] + prev_B*D*a[0][1]),
                                    f*(prev_G*a[0][1] + prev_B*a[0][0])))
                    break
                prev_B, prev_G = _B, _G
    return sol


def cornacchia(a, b, m) -> set[tuple[int, int]]:
    r"""
    Solves `ax^2 + by^2 = m` where `\gcd(a, b) = 1 = gcd(a, m)` and `a, b > 0`.

    Explanation
    ===========

    Uses the algorithm due to Cornacchia. The method only finds primitive
    solutions, i.e. ones with `\gcd(x, y) = 1`. So this method cannot be used to
    find the solutions of `x^2 + y^2 = 20` since the only solution to former is
    `(x, y) = (4, 2)` and it is not primitive. When `a = b`, only the
    solutions with `x \leq y` are found. For more details, see the References.

    Examples
    ========

    >>> from sympy.solvers.diophantine.diophantine import cornacchia
    >>> cornacchia(2, 3, 35) # equation 2x**2 + 3y**2 = 35
    {(2, 3), (4, 1)}
    >>> cornacchia(1, 1, 25) # equation x**2 + y**2 = 25
    {(4, 3)}

    References
    ===========

    .. [1] A. Nitaj, "L'algorithme de Cornacchia"
    .. [2] Solving the diophantine equation ax**2 + by**2 = m by Cornacchia's
        method, [online], Available:
        http://www.numbertheory.org/php/cornacchia.html

    See Also
    ========

    sympy.utilities.iterables.signed_permutations
    """
    a = as_int(a)
    b = as_int(b)
    m = as_int(m)

    if a <= 0 or b <= 0:
        raise ValueError("")
    if gmpy_gcd(a, b) != 1 or gmpy_gcd(a, m) != 1:
        raise ValueError("")

    sols = set()
    for t in sqrt_mod_iter(-b*invert(a, m), m):
        if t < m // 2:
            continue
        t, r = m, t % m
        while (m1 := m - a*r**2) <= 0:
            t, r = r, t % r
        m1, _r = divmod(m1, b)
        if _r:
            continue
        s, _exact = integer_nthroot(m1, 2)
        if _exact:
            if a == b and r < s:
                r, s = s, r
            sols.add((int(r), int(s)))
    return sols
