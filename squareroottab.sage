#!/usr/bin/env sage

# This implements a prototype of Palash Sarkar's square root algorithm
# from <https://eprint.iacr.org/2020/1407>, for the Pasta fields.

from copy import copy

DEBUG = True
VERBOSE = False
EXPENSIVE = False


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
    def __init__(self, p, z, base_cost, hash_bits=None, hash_mod=None):
        n = 32
        m = p >> n
        assert p == 1 + m * 2^n
        if EXPENSIVE: assert Mod(z, p).multiplicative_order() == p-1
        g = Mod(z, p)^m
        if EXPENSIVE: assert g.multiplicative_order() == 2^n

        gtab = [[0]*256 for i in range(4)]
        gi = g
        for i in range(4):
            if DEBUG: assert gi == g^(256^i), (i, gi)
            acc = Mod(1, p)
            for j in range(256):
                 if DEBUG: assert acc == g^(256^i * j), (i, j, acc)
                 gtab[i][j] = acc
                 acc *= gi
            gi = acc

        (self.hash_bits, self.hash_mod) = self.find_perfect_hash(gtab[3]) if hash_bits is None else (hash_bits, hash_mod)

        # Now invert gtab[3].
        invtab = [0]*self.hash_mod
        for j in range(256):
            invtab[self.hash(gtab[3][j])] = (256-j) % 256

        minus1 = Mod(-1, p)

        (self.p, self.n, self.m, self.g, self.gtab, self.invtab, self.minus1, self.base_cost) = (
              p,      n,      m,      g,      gtab,      invtab,      minus1,      base_cost)

    def hash(self, x):
        return (int(x) % (1 << self.hash_bits)) % self.hash_mod

    def find_perfect_hash(self, gt):
        def is_ok(c_bits, c_mod):
            c_invtab = [False]*c_mod
            for j in range(256):
                hash = (int(gt[j]) % (1 << c_bits)) % c_mod
                if c_invtab[hash]:
                    return False
                c_invtab[hash] = True

            return True

        hash_bits = None
        hash_mod = None
        for c_bits in range(8, 32):
            for c_mod in range(256, 10000):
                if is_ok(c_bits, c_mod):
                    if VERBOSE: print(c_bits, c_mod)
                    if hash_mod is None or c_mod < hash_mod:
                        (hash_bits, hash_mod) = (c_bits, c_mod)
                    break

        print("best is hash_bits=%d, hash_mod=%d" % (hash_bits, hash_mod))
        return (hash_bits, hash_mod)

    def mul_by_g_to(self, acc, t, j, k, cost):
        if DEBUG: expected = acc * self.g^t

        t >>= 8*j
        for i in range(j, k):
            acc *= self.gtab[i][t % 256]
            t >>= 8
            cost.muls += 1

        if DEBUG: assert acc == expected, (t, acc, expected)
        return acc

    def eval(self, alpha):
        s = self.invtab[self.hash(alpha)] << 24
        return s

    def sarkar_sqrt(self, u):
        if VERBOSE: print("u = %r" % (u,))

        # This would actually be done using the addition chain.
        v = u^((self.m-1)/2)
        cost = copy(self.base_cost)

        x3 = u * v^2
        cost.sqrs += 1
        cost.muls += 1
        if DEBUG: assert x3 == u^self.m
        if EXPENSIVE: assert x3.multiplicative_order().divides(2^self.n)

        x2 = x3^(1<<8)
        x1 = x2^(1<<8)
        x0 = x1^(1<<8)
        if DEBUG:
            assert x0 == x3^(1<<(self.n-1-7))
            assert x1 == x3^(1<<(self.n-1-15))
            assert x2 == x3^(1<<(self.n-1-23))

        cost.sqrs += 8+8+8

        # i = 0
        s = self.eval(x0)

        # i = 1
        t = s >> 8
        assert t & 0xFF00FFFF == 0, "0x%x" % (t,)
        alpha = self.mul_by_g_to(x1, t, 2, 3, cost)
        s = self.eval(alpha)

        # i = 2
        t = (s+t) >> 8
        assert t & 0xFF0000FF == 0, "0x%x" % (t,)
        alpha = self.mul_by_g_to(x2, t, 1, 3, cost)
        s = self.eval(alpha)

        # i = 3
        t = (s+t) >> 8
        alpha = self.mul_by_g_to(x3, t, 0, 3, cost)
        s = self.eval(alpha)

        t = (s+t) >> 1
        res = self.mul_by_g_to(u * v, t, 0, 4, cost)
        cost.muls += 1

        if res^2 != u:
            res = None
        if DEBUG: assert u.is_square() == (res is not None)
        return (res, cost)


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

# see addchain.py for base costs of u^{(m-1)/2}
F_p = SqrtField(p, 5, Cost(223, 23), hash_bits=24, hash_mod=1415)
F_q = SqrtField(q, 5, Cost(222, 25), hash_bits=30, hash_mod=2384)


print("p = %r" % (p,))

x = Mod(0x1234567890123456789012345678901234567890123456789012345678901234, p)
print(F_p.sarkar_sqrt(x))

x = Mod(0x2345678901234567890123456789012345678901234567890123456789012345, p)
print(F_p.sarkar_sqrt(x))

# nonsquare
x = Mod(0x3456789012345678901234567890123456789012345678901234567890123456, p)
print(F_p.sarkar_sqrt(x))

if True:
    total_cost = Cost(0, 0)
    iters = 10000
    for i in range(iters):
        x = GF(p).random_element()
        (_, cost) = F_p.sarkar_sqrt(x)
        total_cost += cost

    print total_cost.divide(iters)
