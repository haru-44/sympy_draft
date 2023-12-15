# オイラー擬素数を使った部分、オイラー・ヤコビ擬素数を使った方が良かったと思い至る
# 場合分けすればもっと大きくできるが、どこまでするか

if n < 2809:
    return True
# if n < 65077:
#     # There are only five Euler pseudoprimes with a least prime factor greater than 47
#     return pow(2, n >> 1, n) in [1, n - 1] and n not in [8321, 31621, 42799, 49141, 49981]

if n % 6 == 5 and n < 1551941:
    pp = pow(3, n >> 1, n)
    j = 1 if n % 12 in [1, 11] else n - 1
    return pp == j and n not in [432821, 973241]
if n % 6 == 1 and n < 80581:
    pp = pow(2, n >> 1, n)
    j = 1 if n % 8 in [1, 7] else n - 1
    return pp == j and n not in [42799, 49141, 65281]

