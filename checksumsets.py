#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Dependency: <https://pypi.org/project/bintrees/> (pip install bintrees)


# From the Halo paper:

# Let A = [0, 2^{λ/2 + 1} + 2^{λ/2} - 1]. It is straightforward to verify that a, b ∈ A
# at the end of Algorithm 3 for any input \mathbf{r}.
# {In fact a, b ∈ [2^{λ/2} + 1, 2^{λ/2 + 1} + 2^{λ/2} - 1], but it is convenient to
# define A to start at 0.}
#
# Next we need to show that the mapping (a ⦂ A, b ⦂ A) ↦ (a ζ_q + b) mod q is injective.
# This will depend on the specific values of ζ_q and q, and can be cast as a sumset problem.
#
# We use the notation v·A + A for { (av + b) mod q : a, b ∈ A }, and A - A for -1·A + A.
# {We take sumsets using this notation to implicitly be subsets of F_q.}
# The question is then whether |ζ_q·A + A| = |A|^2.
#
# For intuition, note that if av + b = a'v + b' (mod q), with a ≠ a', we would have
# v = (b' - b)/(a - a') (mod q). Thus the number of v ∈ F_q for which |v·A + A| < |A|^2
# is at most (|A - A| - 1)^2. {We thank Robert Israel for this observation. [HI2019]}
# Since in our case (|A - A| - 1)^2 ≈ 9·2^130 is small compared to q ≈ 2^254, we would
# heuristically expect that |ζ_q·A + A| = |A|^2 unless there is some reason why ζ_q does
# not "behave like a random element of F_q".
#
# Of course ζ_q is *not* a random element of F_q, and so the above argument can only be
# used for intuition. Even when (|A - A| - 1)^2 is small compared to q, there are clearly
# values of ζ_q and q for which it would not hold. To prove that it holds in the needed
# cases for the Tweedledum and Tweedledee curves used in our implementation, we take a
# different tack.
#
# Define a distance metric δ_q on F_q so that δ_q(x, y) is the minimum distance between
# x and y around the ring of integers modulo q in either direction, i.e.
#
#    δ_q(x, y) = min(z, q - z) where z = (x - y) mod q
#
# Now let D_{q,ζ_q}(m) be the minimum δ_q-distance between any two elements of ζ_q·[0, m],
# i.e.
#
#    D_{q,ζ_q}(m) = min{ δ_q(a ζ_q, a' ζ_q ) : a, a' ∈ [0, m] }
#
# An algorithm to compute D_{q,ζ_q}(m) is implemented by checksumsets.py in [Hopw2019]
# [i.e. this file]; it works by iteratively finding each m at which D_{q,ζ_q}(m)
# decreases. [...]
#
# Now if D_{q,ζ_q}(2^{λ/2 + 1} + 2^{λ/2} - 1) ≥ 2^{λ/2 + 1} + 2^{λ/2}, then copies of
# A will "fit within the gaps" in ζ_q·A. That is, ζ_q·A + A will have |A|^2 elements,
# because all of the sets { ζ_q·{a} + A : a ∈ A } will be disjoint.
#
# The algorithm is based on the observation that the problem of deciding when
# D_{q,ζ_q}(m) next decreases is self-similar to deciding when it first decreases.
# It computes the exact min-distance at each decrease (not just a lower bound),
# which facilitates detecting any bugs in the algorithm. Also, we check correctness
# of the partial results up to a given bound on m, against a naive algorithm.

BRUTEFORCE_THRESHOLD = 100000

def D(q, zeta, mm):
    if DEBUG: print("(q, zeta, mm) =", (q, zeta, mm))
    Dcheck = bruteforce_D(q, zeta, min(mm, BRUTEFORCE_THRESHOLD))

    (u, m) = (0, 1)  # (u + am) : a ∈ Nat is the current arithmetic progression
    n = q            # the previous min-distance
    d = zeta         # the current min-distance

    while True:
        # Consider values of x where D_{q,ζ_q}(x) decreases, i.e. where
        # D_{q,ζ_q}(x) < D_{q,ζ_q})(x-1).
        #
        # We keep track of an arithmetic progression (u + am) such that the next
        # value at which D_{q,ζ_q}(x) decreases will be for x in this progression,
        # at the point at which xζ gets close to (but not equal to) 0.
        #
        # TODO: explain why the target is always 0.
        #
        # If we set s = floor(n/d), then D_{q,ζ_q}(x) can decrease at a = s
        # and potentially also at a = s+1.
        assert (m*zeta) % q in (d, q-d)
        s = n // d
        x0 = u + s*m
        d0 = n % d
        if DEBUG: print("(x0, d0, u, m, n, d, s) =", (x0, d0, u, m, n, d, s))
        assert dist(0, x0*zeta, q) == d0
        if x0-1 < len(Dcheck): assert Dcheck[x0-1] == d
        if x0 > mm: return d
        if x0 < len(Dcheck): assert Dcheck[x0] == d0

        x1 = u + (s+1)*m
        d1 = (s+1)*d - n
        if d1 < d0:
            if DEBUG: print("(x1, d1, u, m, n, d, s+1) =", (x1, d1, u, m, n, d, s+1))
            assert dist(0, x1*zeta, q) == d1
            if x1-1 < len(Dcheck): assert Dcheck[x1-1] == d0
            if x1 > mm: return d0
            if x1 < len(Dcheck): assert Dcheck[x1] == d1

            # TODO: This is painfully non-obvious! Need to draw some diagrams to explain it.
            (u, m, n, d) = (x0, x1, d0, d1)
        else:
            (u, m, n, d) = (x1, x0, d-d0, d0)


def bruteforce_D(q, zeta, mm):
    # Can't use sortedcontainers because its data structures are backed by
    # lists-of-lists, not trees. We must have O(log n) insert, prev, and succ.
    from bintrees import RBTree as sortedset
    from collections import deque

    resD = deque([zeta])
    lastd = zeta
    S = sortedset()
    S.insert(0, None)
    S.insert(q, None)
    for x in range(1, mm+1):
        v = (x*zeta) % q
        S.insert(v, None)
        vp = S.prev_key(v)
        vs = S.succ_key(v)
        d = min(v-vp, vs-v)
        resD.append(d)
        if DEBUG and d < lastd: print((x, d))
        lastd = d

    return list(resD)

def dist(x, y, q):
    z = (x-y+q) % q
    return min(z, q-z)

def check_sumset(name, q, zeta, limit):
    print("===== %s =====" % (name,))
    Dq = D(q, zeta, limit-1)
    print("D_%s = %s" % (name, Dq), end=' ')
    assert Dq >= limit
    print(">=", limit)


halflambda = 64
limit = 3<<halflambda

# Tweedledum and Tweedledee
p = (1<<254) + 4707489545178046908921067385359695873
q = (1<<254) + 4707489544292117082687961190295928833
zeta_p = 9513155655832138286304767221959569637168364952810827555227185832555034233288
zeta_q = 24775483399512474214391554062650059912556682109176536098332128018848638018813

# Tests
DEBUG = False
assert(D(65537, 6123, 10000) == 3)
assert(D(1299721, 538936, 10000) == 41)
assert(D(179424691, 134938504, 100000) == 121)

DEBUG = False
check_sumset("p", p, zeta_p, limit)
check_sumset("q", q, zeta_q, limit)
