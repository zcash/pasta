#!/usr/bin/python

from itertools import product

def c_d():
    yield (-1,  0)
    yield ( 1,  0)
    yield ( 0, -1)
    yield ( 0,  1)

def sums_mod4(cd):
    return ((2**(k+1) + sum([cd[j][0] * (2**j) for j in range(k)])) % 4,
            (2**(k+1) + sum([cd[j][1] * (2**j) for j in range(k)])) % 4)

for k in (1, 2):
    M_k = [list(s) for s in product(c_d(), repeat=k)]
    assert(len(M_k) == 4**k)

    for (cd, cd_dash) in product(M_k, repeat=2):
        if cd[0] != cd_dash[0]:
            assert(sums_mod4(cd) != sums_mod4(cd_dash))

print("QED")
