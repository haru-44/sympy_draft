from sympy import Function, S, Sum, divisor_sigma, Mul, factorint, Pow
from sympy.ntheory.factor_ import perfect_power
from sympy.ntheory import isprime

class ramanujan_tau(Function):
    """ https://en.wikipedia.org/wiki/Ramanujan_tau_function
    """
    is_integer = True
    
    @classmethod
    def eval(cls, n):
        if n.is_integer is False:
            raise TypeError("n should be an integer")
        if n.is_positive is False:
            raise ValueError("n should be a positive integer")
        if n is S.One:
            return S.One

    def _eval_rewrite(self, rule, args, **hints):
        n = args[0]
        if rule == Mul:
            if n.is_Integer is True:
                return Mul(*[ramanujan_tau(p**e) for p, e in factorint(n).items()])
        if rule == ramanujan_tau:
            if n.is_Integer is True:
                pe = perfect_power(n)
                if not pe:
                    return
                p, e = pe
                if isprime(p) and 1 < e:
                    return ramanujan_tau(p)*ramanujan_tau(p**(e-1)) - p**11*ramanujan_tau(p**(e-2))
                return
            p, e = n.as_base_exp()
            if p.is_prime is True and e.is_Integer and 1 < e:
                return ramanujan_tau(p)*ramanujan_tau(p**(e-1)) - p**11*ramanujan_tau(p**(e-2))
        if rule == Sum:
            from sympy import Symbol
            i = Symbol('i', integer=True, positive=True)
            return n**4 * divisor_sigma(n) -\
                24*Sum(i**2 * (35*i**2 - 52*i*n + 18*n**2) * divisor_sigma(i) * divisor_sigma(n - i), (i, 1, n - 1))
    
    def doit(self, deep=False, **hints):
        n = self.args[0]
        if deep:
            n = n.doit(deep=deep, **hints)
        if n.is_Integer is not True:
            return self
        factors = factorint(n)
        sigma_sieve = [1] * (max(factors.keys()) + 1)
        for i in range(2, len(sigma_sieve)):
            for j in range(i, len(sigma_sieve), i):
                sigma_sieve[j] += i
        # sigma_sieve[i] = divisor_sigma(i)
        result = 1
        for p, e in factors.items():
            tau_p = p**4 * sigma_sieve[p] - 24*sum(i**2*(35*i**2 - 52*i*p + 18*p**2)*sigma_sieve[i]*sigma_sieve[p - i]
                           for i in range(1, p))
            t1, t2 = tau_p, 1
            for _ in range(e - 1):
                t1, t2 = tau_p*t1 - p**11*t2, t1
            result *= t1
        return result
    
    def _latex(self, printer, exp=None):
        if exp is None:
            return r'\tau\left( %s \right)' % printer._print(self.args[0])
        return r'\tau^{%s}\left( %s \right)' % (exp, printer._print(self.args[0]))
