#!/usr/bin/python3
# -*- coding: utf-8 -*-

# This checks the cases k = 1 and k = 2 needed to complete the proof of a
# hexary version of the Halo injectivity lemma.
#
# Inputs: r ∈ {0, 1}^λ, s ∈ [0, 2)^λ, P ∈ E \ {O}
#
# (a : Fq, b : Fq, c : Fq) := (2, 2, 2)
#
# for i from λ - 1 down to 0:
#     let (d_i, e_i, f_i) = { (0, 0, 2r_i - 1), if s_i = 0
#                           { (0, 2r_i - 1, 0), if s_i = 1
#                           { (2r_i - 1, 0, 0), if s_i = 2
#     (a, b, c) := (2a + d_i, 2b + e_i, 2c + f_i)
#     Output [a.ζ^2 + b.ζ + c] P.
#
# For each i ∈ [0, λ), the mapping (r_i, s_i) ↦ (d_i, e_i, f_i) is injective,
# and exactly one of d_i, e_i, f_i ∈ {-1, 0, +1} is 0.
# Let M_k = {(d, e, f : {-1, 0, +1}^k) for all i, exactly one of d_i, e_i, f_i is 0}.
# So (r, s) ↦ (d, e, f) : M_λ is also injective.
#
# Lemma: For k ≥ 0, (d, e, f) ∈ M_k ↦ (∑_{j ∈ [0, k-1)} d_j 2^j,
#                                      ∑_{j ∈ [0, k-1)} e_j 2^j,
#                                      ∑_{j ∈ [0, k-1)} f_j 2^j) is injective.
#
# Proof (sketch). If (d, e, f) and (d', e', f') coincide on a prefix of length m,
# then the statement reduces to a smaller instance of the lemma with that prefix
# deleted and k reduced by m. So we need only consider the case
# (d_0, e_0, f_0) ≠ (d'_0, e'_0, f'_0) and show that the resulting sums always
# differ. In fact they always differ modulo 4:
#
# (d_0, e_0, f_0) ≠ (d'_0, e'_0, f'_0) ⇒ (∑_{j ∈ [0, k-1)}  d_j 2^j,  ∑_{j ∈ [0, k-1)}  e_j 2^j,  ∑_{j ∈ [0, k-1)}  f_j 2^j)
#                                      ≠ (∑_{j ∈ [0, k-1)} d'_j 2^j,  ∑_{j ∈ [0, k-1)} e'_j 2^j,  ∑_{j ∈ [0, k-1)} f'_j 2^j)
#
# Therefore, it is sufficient to verify this property exhaustively for k = 1 and
# k = 2, since terms with j ≥ 2 do not affect the sums modulo 4.

from itertools import product
from collections import namedtuple

def_triple = namedtuple('def_triple', ['d', 'e', 'f'])

def d_e_f():
    yield def_triple(-1,  0,  0)
    yield def_triple( 1,  0,  0)
    yield def_triple( 0, -1,  0)
    yield def_triple( 0,  1,  0)
    yield def_triple( 0,  0,  1)
    yield def_triple( 0,  0, -1)

def sums_mod4(def_):
    return ((2**(k+1) + sum([def_[j].d * (2**j) for j in range(k)])) % 4,
            (2**(k+1) + sum([def_[j].e * (2**j) for j in range(k)])) % 4,
            (2**(k+1) + sum([def_[j].f * (2**j) for j in range(k)])) % 4)

for k in (1, 2):
    M_k = [list(s) for s in product(d_e_f(), repeat=k)]
    assert(len(M_k) == 6**k)

    for (def_, def_dash) in product(M_k, repeat=2):
        if def_[0] != def_dash[0]:
            assert(sums_mod4(def_) != sums_mod4(def_dash))

print("QED")
