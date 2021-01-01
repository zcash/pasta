#!/usr/bin/env sage

# This implements a prototype of Palash Sarkar's square root algorithm
# from <https://eprint.iacr.org/2020/1407>, for the Pasta fields.

import sys

if sys.version_info[0] == 2:
    range = xrange

DEBUG = True
VERBOSE = False
EXPENSIVE = False

SUBGROUP_TEST = True
OP_COUNT = True


class Cost:
    def __init__(self, sqrs=0, muls=0, invs=0):
        self.sqrs = sqrs
        self.muls = muls
        self.invs = invs

    def sqr(self, x):
        self.sqrs += 1
        return x^2

    def mul(self, x, y):
        self.muls += 1
        return x * y

    def div(self, x, y):
        self.invs += 1
        self.muls += 1
        return x / y

    def batch_inv0(self, xs):
        self.invs += 1
        self.muls += 3*(len(xs)-1)
        # This should use Montgomery's trick (with constant-time substitutions to handle zeros).
        return [0 if x == 0 else x^-1 for x in xs]

    def __repr__(self):
        return "%dS + %dM + %dI" % (self.sqrs, self.muls, self.invs)

    def __add__(self, other):
        return Cost(self.sqrs + other.sqrs, self.muls + other.muls, self.invs + other.invs)

    def include(self, other):
        self.sqrs += other.sqrs
        self.muls += other.muls

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
        invtab = [1]*hash_mod
        for j in range(256):
            h = self.hash(gtab[3][j])
            # 1 is the last value to be assigned, so this ensures there are no collisions.
            assert invtab[h] == 1
            invtab[h] = (256-j) % 256

        gtab[3] = gtab[3][:129]

        (self.p, self.n, self.m, self.g, self.gtab, self.invtab, self.base_cost) = (
              p,      n,      m,      g,      gtab,      invtab,      base_cost)

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

    """
    Return (sqrt(u),   True ), if u is square in the field.
           (sqrt(g*u), False), otherwise.
    """
    def sarkar_sqrt(self, u, c):
        if VERBOSE: print("u = %r" % (u,))

        # This would actually be done using the addition chain.
        v = u^((self.m-1)/2)
        c.include(self.base_cost)

        uv = c.mul(u, v)
        (res, zero_if_square) = self.sarkar_sqrt_common(u, 1, uv, v, c)
        return (res, zero_if_square)

    """
    Return (sqrt(N/D),   True,  c), if N/D is square in the field.
           (sqrt(g*N/D), False, c), otherwise.

    This avoids the full cost of computing N/D.
    """
    def sarkar_divsqrt(self, N, D, c):
        if DEBUG:
            u = N/D
            if VERBOSE: print("N/D = %r/%r\n    = %r" % (N, D, u))

        # We need to calculate uv and v, where v = u^((m-1)/2), u = N/D, and p-1 = m * 2^n.
        # We can rewrite as follows:
        #
        #      v = (N/D)^((m-1)/2)
        #        = N^((m-1)/2) * D^(p-1 - (m-1)/2)    [Fermat's Little Theorem]
        #        =      "      * D^(m * 2^n - (m-1)/2)
        #        =      "      * D^((2^(n+1) - 1)*(m-1)/2 + 2^n)
        #        = (N * D^(2^(n+1) - 1))^((m-1)/2) * D^(2^n)
        #
        # Let  w = (N * D^(2^(n+1) - 1))^((m-1)/2) * D^(2^n - 1).
        # Then v = w * D, and uv = N * v/D = N * w.
        #
        # We calculate:
        #
        #      s = D^(2^n - 1) using an addition chain
        #      t = D^(2^(n+1) - 1) = s^2 * D
        #      w = (N * t)^((m-1)/2) * s using another addition chain
        #
        # then u and uv as above. The addition chains are given in addchain_sqrt.py .
        # The overall cost of this part is similar to a single full-width exponentiation,
        # regardless of n.

        s = D^(2^self.n - 1)
        c.sqrs += 31
        c.muls += 5
        t = c.mul(c.sqr(s), D)
        if DEBUG: assert t == D^(2^(self.n+1) - 1)
        w = c.mul(c.mul(N, t)^((self.m-1)/2), s)
        c.include(self.base_cost)
        v = c.mul(w, D)
        uv = c.mul(N, w)

        if DEBUG:
            assert v == u^((self.m-1)/2)
            assert uv == u * v

        (res, zero_if_square) = self.sarkar_sqrt_common(N, D, uv, v, c)

        if DEBUG:
            (res_ref, zero_if_square_ref) = self.sarkar_sqrt(u, Cost())
            assert res == res_ref
            assert (zero_if_square == 0) == (zero_if_square_ref == 0)

        return (res, zero_if_square)

    def sarkar_sqrt_common(self, N, D, uv, v, c):
        x3 = uv * v
        c.muls += 2
        if DEBUG:
            u = N/D
            assert x3 == u^self.m
        if EXPENSIVE:
            x3_order = x3.multiplicative_order()
            if VERBOSE: print("x3_order = %r" % (x3_order,))
            # x3_order is 2^n iff u is nonsquare, otherwise it divides 2^(n-1).
            assert x3.divides(2^self.n)

        x2 = x3^(1<<8)
        x1 = x2^(1<<8)
        x0 = x1^(1<<8)
        if DEBUG:
            assert x0 == x3^(1<<(self.n-1-7))
            assert x1 == x3^(1<<(self.n-1-15))
            assert x2 == x3^(1<<(self.n-1-23))

        c.sqrs += 8+8+8

        # i = 0, 1
        t_ = self.invtab[self.hash(x0)]  # = t >> 16
        if DEBUG: assert 1 == x0 * self.g^(t_ << 24), (x0, t_)
        assert t_ < 0x100, t_
        alpha = x1 * self.gtab[2][t_]
        c.muls += 1

        # i = 2
        t_ += self.invtab[self.hash(alpha)] << 8  # = t >> 8
        if DEBUG: assert 1 == x1 * self.g^(t_ << 16), (x1, t_)
        assert t_ < 0x10000, t_
        alpha = x2 * self.gtab[1][t_ % 256] * self.gtab[2][t_ >> 8]
        c.muls += 2

        # i = 3
        t_ += self.invtab[self.hash(alpha)] << 16  # = t
        if DEBUG: assert 1 == x2 * self.g^(t_ << 8), (x2, t_)
        assert t_ < 0x1000000, t_
        alpha = x3 * self.gtab[0][t_ % 256] * self.gtab[1][(t_ >> 8) % 256] * self.gtab[2][t_ >> 16]
        c.muls += 3

        t_ += self.invtab[self.hash(alpha)] << 24  # = t << 1
        if DEBUG: assert 1 == x3 * self.g^t_, (x3, t_)
        t_ = (t_ + 1) >> 1
        assert t_ <= 0x80000000, t_
        res = uv * self.gtab[0][t_ % 256] * self.gtab[1][(t_ >> 8) % 256] * self.gtab[2][(t_ >> 16) % 256] * self.gtab[3][t_ >> 24]
        c.muls += 4

        zero_if_square = c.mul(c.sqr(res), D) - N
        if DEBUG:
            assert (zero_if_square == 0) == u.is_square()
            if EXPENSIVE: assert (zero_if_square == 0) == (x3_order != 2^self.n), (zero_if_square, x3_order)
            if zero_if_square != 0:
                assert(res^2 == u * self.g)

        return (res, zero_if_square)


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

# see addchain_sqrt.py for base costs of u^{(m-1)/2}
F_p = SqrtField(p, 5, Cost(223, 23), hash_xor=0x11BE,   hash_mod=1098)
F_q = SqrtField(q, 5, Cost(223, 24), hash_xor=0x116A9E, hash_mod=1206)

print("p = %r" % (p,))

x = Mod(0x1234567890123456789012345678901234567890123456789012345678901234, p)
print(F_p.sarkar_sqrt(x, Cost()))
Dx = Mod(0x123456, p)
print(F_p.sarkar_divsqrt(x*Dx, Dx, Cost()))

x = Mod(0x2345678901234567890123456789012345678901234567890123456789012345, p)
print(F_p.sarkar_sqrt(x, Cost()))

# nonsquare
x = Mod(0x3456789012345678901234567890123456789012345678901234567890123456, p)
print(F_p.sarkar_sqrt(x, Cost()))

if SUBGROUP_TEST:
    for i in range(33):
        x = F_p.g^(2^i)
        print(F_p.sarkar_sqrt(x, Cost()))

if OP_COUNT:
    cost = Cost()
    iters = 50
    for i in range(iters):
        x = GF(p).random_element()
        y = GF(p).random_element()
        (_, _) = F_p.sarkar_divsqrt(x, y, cost)

    print cost.divide(iters)
