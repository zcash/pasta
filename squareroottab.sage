#!/usr/bin/env sage

# This implements a prototype of Palash Sarkar's square root algorithm
# from <https://eprint.iacr.org/2020/1407>, for the Pasta fields.

from copy import copy

if sys.version_info[0] == 2:
    range = xrange

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
    def __init__(self, p, z, base_cost, hash_xor=None, hash_mod=None):
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

        if hash_xor is None:
            (hash_xor, hash_mod) = self.find_perfect_hash(gtab[3])
        (self.hash_xor, self.hash_mod) = (hash_xor, hash_mod)

        # Now invert gtab[3].
        invtab = [0]*hash_mod
        for j in range(256):
            invtab[self.hash(gtab[3][j])] = (256-j) % 256

        gtab[3] = gtab[3][:128]

        (self.p, self.n, self.m, self.gtab, self.invtab, self.base_cost) = (
              p,      n,      m,      gtab,      invtab,      base_cost)

    def hash(self, x):
        return ((int(x) & 0xFFFFFFFF) ^^ self.hash_xor) % self.hash_mod

    def find_perfect_hash(self, gt):
        gt = [int(x) & 0xFFFFFFFF for x in gt]
        assert len(set(gt)) == len(gt)

        def is_ok(c_invtab, c_xor, c_mod):
            for j in range(256):
                hash = (gt[j] ^^ c_xor) % c_mod
                if c_invtab[hash] == c_mod:
                    return False
                c_invtab[hash] = c_mod

            return True

        hash_xor = None
        hash_mod = 10000
        for c_xor in range(1, 0x200000):
            c_invtab = [0]*hash_mod
            for c_mod in range(256, hash_mod):
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
        x3 = uv * v
        cost.muls += 2
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

        # i = 0, 1
        t_ = self.invtab[self.hash(x0)]  # = t >> 16
        assert t_ < 0x100, t_
        alpha = x1 * self.gtab[2][t_]
        cost.muls += 1

        # i = 2
        t_ += self.invtab[self.hash(alpha)] << 8  # = t >> 8
        assert t_ < 0x10000, t_
        alpha = x2 * self.gtab[1][t_ % 256] * self.gtab[2][t_ >> 8]
        cost.muls += 2

        # i = 3
        t_ += self.invtab[self.hash(alpha)] << 16  # = t
        assert t_ < 0x1000000, t_
        alpha = x3 * self.gtab[0][t_ % 256] * self.gtab[1][(t_ >> 8) % 256] * self.gtab[2][t_ >> 16]
        cost.muls += 3

        t_ = ((self.invtab[self.hash(alpha)] << 24) + t_) >> 1  # = t
        assert t_ < 0x80000000, t_
        res = uv * self.gtab[0][t_ % 256] * self.gtab[1][(t_ >> 8) % 256] * self.gtab[2][(t_ >> 16) % 256] * self.gtab[3][t_ >> 24]
        cost.muls += 4

        if res^2 != u:
            res = None
        cost.sqrs += 1
        if DEBUG: assert u.is_square() == (res is not None)
        return (res, cost)


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

# see addchain.py for base costs of u^{(m-1)/2}
F_p = SqrtField(p, 5, Cost(223, 23), hash_xor=0x11BE,   hash_mod=1098)
F_q = SqrtField(q, 5, Cost(223, 24), hash_xor=0x116A9E, hash_mod=1206)


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
    iters = 1000
    for i in range(iters):
        x = GF(p).random_element()
        (_, cost) = F_p.sarkar_sqrt(x)
        total_cost += cost

    print total_cost.divide(iters)
