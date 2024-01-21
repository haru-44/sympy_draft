from itertools import count

@recurrence_memo([1, 1])
def _partition_rec(n, prev):
    v = 0
    penta = 0 # pentagonal number: 1, 5, 12, ...
    for i in count():
        penta += 3*i + 1
        np = n - penta
        if np < 0:
            break
        s = prev[np]
        np -= i + 1
        # np = n - gp where gp = generalized pentagonal: 2, 7, 15, ...
        if 0 <= np:
            s += prev[np]
        v += -s if i % 2 else s
    return v
