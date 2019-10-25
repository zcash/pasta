#!/usr/bin/python

from itertools import product

def gamma_delta():
    yield (-1,  0)
    yield ( 1,  0)
    yield ( 0, -1)
    yield ( 0,  1)

def sums_mod4(gd):
    return (sum([gd[j][0] * (2**j) for j in range(k)]) % 4, sum([gd[j][1] * (2**j) for j in range(k)]) % 4)

for k in (1, 2):
    M_k = [list(s) for s in product(gamma_delta(), repeat=k)]
    assert(len(M_k) == 4**k)

    for (gd, gd_dash) in product(M_k, repeat=2):
        if gd[0] != gd_dash[0]:
            assert(sums_mod4(gd) != sums_mod4(gd_dash))

print("QED")
