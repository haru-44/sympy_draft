# わざわざ関数化する必要あるか？

from sympy.ntheory.factor_ import factorint
from sympy.functions.combinatorial.numbers import partition

def abelian_groups_count(n):
    prod = 1
    for e in factorint(n).values():
        prod *= partition(e)
    return int(prod)

def nonabelian_groups_count(n):
    return groups_count(n) - abelian_groups_count(n)
