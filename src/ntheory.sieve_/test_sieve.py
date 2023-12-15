from sympy.ntheory.sieve_ import FactorSieve, TotientSieve, MobiusSieve
from sympy.ntheory.factor_ import factorint, totient
from sympy.ntheory.residue_ntheory import mobius


def test_factorsieve():
    fs = FactorSieve(51)
    for n in range(2, 100):
        f = factorint(n)
        ff = sorted(f.keys())
        assert fs[n] == ff[0]
        assert fs.factorint(n) == f
    fs.reset()
    fs.extend(10)



def test_totientsieve():
    ts = TotientSieve(51)
    for n in range(1, 100):
        assert ts[n] == totient(n)


def test_mobiussieve():
    ms = MobiusSieve(51)
    for n in range(1, 100):
        assert ms[n] == mobius(n)
