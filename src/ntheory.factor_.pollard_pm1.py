def pollard_pm1(n, B=10, a=2, retries=0, seed=1234, B2=None):
    n = int(n)
    if n < 4 or B < 3:
        raise ValueError('pollard_pm1 should receive n > 3 and B > 2')
    B2 = B2 or B
    sieve.extend(isqrt(B2))
    randint = _randint(seed + B)

    # computing a**lcm(1,2,3,..B) % n for B > 2
    # it looks weird, but it's right: primes run [2, B]
    # and the answer's not right until the loop is done.
    for _ in range(retries + 1):
        aM = a
        for p in sieve.primerange(2, B + 1):
            e = int(math.log(B, p))
            aM = pow(aM, pow(p, e), n)
        g = gcd(aM - 1, n)
        if g == n:
            continue
        if g != 1:
            return int(g)
        if B2 <= B:
            continue
        q = nextprime(p)
        aMq = pow(aM, q, n)
        M = aMq - 1
        dic = {}
        from itertools import pairwise
        for p1, p2 in pairwise(sieve.primerange(q, B2 + 1)):
            if p2 - p1 not in dic:
                dic[p2 - p1] = pow(aM, p2 - p1, n)
            aMq = aMq * dic[p2 - p1] % n
            M = M * (aMq - 1) % n
        g = gcd(M, n)
        if 1 < g < n:
            return int(g)
        # get a new a:
        # since the exponent, lcm(1..B), is even, if we allow 'a' to be 'n-1'
        # then (n - 1)**even % n will be 1 which will give a g of 0 and 1 will
        # give a zero, too, so we set the range as [2, n-2]. Some references
        # say 'a' should be coprime to n, but either will detect factors.
        a = randint(2, n - 2)
