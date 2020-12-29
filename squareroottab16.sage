#!/usr/bin/env sage

# This implements a prototype of Palash Sarkar's square root algorithm
# from <https://eprint.iacr.org/2020/1407>, for the Pasta fields.

import sys
from copy import copy

if sys.version_info[0] == 2:
    range = xrange

DEBUG = False
VERBOSE = False
EXPENSIVE = False

SUBGROUP_TEST = True
OP_COUNT = True

class Cost:
    def __init__(self, sqrs, muls):
        self.sqrs = sqrs
        self.muls = muls

    def __repr__(self):
        return repr((self.sqrs, self.muls))

    def __add__(self, other):
        return Cost(self.sqrs + other.sqrs, self.muls + other.muls)

    def divide(self, divisor):
        return Cost((self.sqrs / divisor).numerical_approx(), (self.muls / divisor).numerical_approx())


class SqrtField:
    def __init__(self, p, z, base_cost, hash_xor=None, hash_mod=None):
        n = 32
        m = p >> n
        assert p == 1 + m * 2^n
        if EXPENSIVE: assert Mod(z, p).multiplicative_order() == p-1
        g = Mod(z, p)^m
        if EXPENSIVE: assert g.multiplicative_order() == 2^n

        gtab = [[0]*16 for i in range(8)]
        gi = g
        for i in range(8):
            if DEBUG: assert gi == g^(16^i), (i, gi)
            acc = Mod(1, p)
            for j in range(16):
                 if DEBUG: assert acc == g^(16^i * j), (i, j, acc)
                 gtab[i][j] = acc
                 acc *= gi
            gi = acc

        if hash_xor is None:
            (hash_xor, hash_mod) = self.find_perfect_hash(gtab[7])
        (self.hash_xor, self.hash_mod) = (hash_xor, hash_mod)

        # Now invert gtab[7].
        invtab = [1]*hash_mod
        for j in range(16):
            h = self.hash(gtab[7][j])
            # 1 is the last value to be assigned, so this ensures there are no collisions.
            assert invtab[h] == 1
            invtab[h] = (16-j) % 16

        gtab[7] = gtab[7][:8]

        (self.p, self.n, self.m, self.g, self.gtab, self.invtab, self.base_cost) = (
              p,      n,      m,      g,      gtab,      invtab,      base_cost)

    def hash(self, x):
        return ((int(x) & 0xFFFFFFFF) ^^ self.hash_xor) % self.hash_mod

    def find_perfect_hash(self, gt):
        gt = [int(x) & 0xFFFFFFFF for x in gt]
        assert len(set(gt)) == len(gt)

        def is_ok(c_invtab, c_xor, c_mod):
            for j in range(16):
                hash = (gt[j] ^^ c_xor) % c_mod
                if c_invtab[hash] == c_mod:
                    return False
                c_invtab[hash] = c_mod

            return True

        hash_xor = None
        hash_mod = 10000
        for c_xor in range(0, 0x200000):
            c_invtab = [0]*hash_mod
            for c_mod in range(16, hash_mod):
                if is_ok(c_invtab, c_xor, c_mod):
                    (hash_xor, hash_mod) = (c_xor, c_mod)
                    print("0x%X: %d" % (hash_xor, hash_mod))
                    break

        print("best is hash_xor=0x%X, hash_mod=%d" % (hash_xor, hash_mod))
        return (hash_xor, hash_mod)

    def sarkar_sqrt(self, u):
        if VERBOSE: print("u = %r" % (u,))

        # This would actually be done using the addition chain.
        v = u^((self.m-1)/2)
        cost = copy(self.base_cost)

        uv = u * v
        x7 = uv * v
        cost.muls += 2
        if DEBUG: assert x7 == u^self.m
        if EXPENSIVE:
            x7_order = x7.multiplicative_order()
            if VERBOSE: print("x7_order = %r" % (x7_order,))
            # x7_order is 2^n iff u is nonsquare, otherwise it divides 2^(n-1).
            assert x7.divides(2^self.n)

        x6 = x7^(1<<4)
        x5 = x6^(1<<4)
        x4 = x5^(1<<4)
        x3 = x4^(1<<4)
        x2 = x3^(1<<4)
        x1 = x2^(1<<4)
        x0 = x1^(1<<4)
        if DEBUG:
            assert x0 == x7^(1<<(self.n-1-3))
            assert x1 == x7^(1<<(self.n-1-7))
            assert x2 == x7^(1<<(self.n-1-11))
            assert x3 == x7^(1<<(self.n-1-15))
            assert x4 == x7^(1<<(self.n-1-19))
            assert x5 == x7^(1<<(self.n-1-23))
            assert x6 == x7^(1<<(self.n-1-27))

        cost.sqrs += 4*7

        # i = 0, 1
        t_ = self.invtab[self.hash(x0)]  # = t >> 24
        if DEBUG: assert 1 == x0 * self.g^(t_ << 28), (x0, t_)
        assert t_ < 0x10, t_
        alpha = x1 * self.gtab[6][t_]
        cost.muls += 1

        # i = 2
        t_ += self.invtab[self.hash(alpha)] << 4  # = t >> 20
        if DEBUG: assert 1 == x1 * self.g^(t_ << 24), (x1, t_)
        assert t_ < 0x100, t_
        alpha = x2 * self.gtab[5][t_ % 16] * self.gtab[6][t_ >> 4]
        cost.muls += 2

        # i = 3
        t_ += self.invtab[self.hash(alpha)] << 8  # = t >> 16
        if DEBUG: assert 1 == x2 * self.g^(t_ << 20), (x2, t_)
        assert t_ < 0x1000, t_
        alpha = x3 * self.gtab[4][t_ % 16] * self.gtab[5][(t_ >> 4) % 16] * self.gtab[6][t_ >> 8]
        cost.muls += 2

        # i = 4
        t_ += self.invtab[self.hash(alpha)] << 12  # = t >> 12
        if DEBUG: assert 1 == x3 * self.g^(t_ << 16), (x3, t_)
        assert t_ < 0x10000, t_
        alpha = x4 * self.gtab[3][t_ % 16] * self.gtab[4][(t_ >> 4) % 16] * self.gtab[5][(t_ >> 8) % 16] * self.gtab[6][t_ >> 12]
        cost.muls += 4

        # i = 5
        t_ += self.invtab[self.hash(alpha)] << 16  # = t >> 8
        if DEBUG: assert 1 == x4 * self.g^(t_ << 12), (x4, t_)
        assert t_ < 0x100000, t_
        alpha = x5 * self.gtab[2][t_ % 16] * self.gtab[3][(t_ >> 4) % 16] * self.gtab[4][(t_ >> 8) % 16] * self.gtab[5][(t_ >> 12) % 16] * self.gtab[6][t_ >> 16]
        cost.muls += 5

        # i = 6
        t_ += self.invtab[self.hash(alpha)] << 20  # = t >> 4
        if DEBUG: assert 1 == x5 * self.g^(t_ << 8), (x5, t_)
        assert t_ < 0x1000000, t_
        alpha = x6 * self.gtab[1][t_ % 16] * self.gtab[2][(t_ >> 4) % 16] * self.gtab[3][(t_ >> 8) % 16] * self.gtab[4][(t_ >> 12) % 16] * self.gtab[5][(t_ >> 16) % 16] * self.gtab[6][t_ >> 20]
        cost.muls += 6

        # i = 7
        t_ += self.invtab[self.hash(alpha)] << 24  # = t
        if DEBUG: assert 1 == x6 * self.g^(t_ << 4), (x6, t_)
        assert t_ < 0x10000000, t_
        alpha = x7 * self.gtab[0][t_ % 16] * self.gtab[1][(t_ >> 4) % 16] * self.gtab[2][(t_ >> 8) % 16] * self.gtab[3][(t_ >> 12) % 16] * self.gtab[4][(t_ >> 16) % 16] * self.gtab[5][(t_ >> 20) % 16] * self.gtab[6][t_ >> 24]
        cost.muls += 7

        t_ += self.invtab[self.hash(alpha)] << 28  # = t << 1
        if DEBUG: assert 1 == x7 * self.g^t_, (x7, t_)
        t_ >>= 1
        assert t_ < 0x80000000, t_
        res = uv * self.gtab[0][t_ % 16] * self.gtab[1][(t_ >> 4) % 16] * self.gtab[2][(t_ >> 8) % 16] * self.gtab[3][(t_ >> 12) % 16] * self.gtab[4][(t_ >> 16) % 16] * self.gtab[5][(t_ >> 20) % 16] * self.gtab[6][(t_ >> 24) % 16] * self.gtab[7][t_ >> 28]
        cost.muls += 8

        if res^2 != u:
            res = None
        cost.sqrs += 1
        if DEBUG:
            issq = u.is_square()
            assert issq == (res is not None)
            if EXPENSIVE: assert issq == (x7_order != 2^self.n), (issq, x7_order)
        return (res, cost)


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

# see addchain_sqrt.py for base costs of u^{(m-1)/2}
F_p = SqrtField(p, 5, Cost(223, 23), hash_xor=0xC847, hash_mod=17)
F_q = SqrtField(q, 5, Cost(223, 24), hash_xor=0xFA68, hash_mod=17)


print("p = %r" % (p,))

x = Mod(0x1234567890123456789012345678901234567890123456789012345678901234, p)
print(F_p.sarkar_sqrt(x))

x = Mod(0x2345678901234567890123456789012345678901234567890123456789012345, p)
print(F_p.sarkar_sqrt(x))

# nonsquare
x = Mod(0x3456789012345678901234567890123456789012345678901234567890123456, p)
print(F_p.sarkar_sqrt(x))

if SUBGROUP_TEST:
    for i in range(33):
        x = F_p.g^(2^i)
        print(F_p.sarkar_sqrt(x))

if OP_COUNT:
    total_cost = Cost(0, 0)
    iters = 1000
    for i in range(iters):
        x = GF(p).random_element()
        (_, cost) = F_p.sarkar_sqrt(x)
        total_cost += cost

    print total_cost.divide(iters)
