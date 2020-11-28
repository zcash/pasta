#!/usr/bin/env sage

# This implements a prototype of Palash Sarkar's square root algorithm
# from <https://eprint.iacr.org/2020/1407>, for the Pasta fields.

from copy import copy

DEBUG = False

def count_bits(x):
    return len(format(x, 'b'))

def count_ones(x):
    return sum([int(b) for b in format(x, 'b')])


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
    def __init__(self, p, z, base_cost):
        n = 32
        m = p >> n
        assert p == 1 + m * 2^n
        assert Mod(z, p).multiplicative_order() == p-1
        g = Mod(z, p)^m

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

        minus1 = Mod(-1, p)

        (self.p, self.n, self.m, self.g, self.gtab, self.minus1, self.base_cost) = (
              p,      n,      m,      g,      gtab,      minus1,      base_cost)

        for k in range(32):
            self.g_to_power_of_2(k)

    def g_to_power_of_2(self, k):
        expected = self.g^(2^k)
        res = self.gtab[k // 8][1<<(k % 8)]
        assert res == expected
        return res

    def mul_by_g_to(self, acc, t, cost):
        if DEBUG: print(t, count_bits(t), count_ones(t))
        expected = acc * self.g^t

        for i in range(4):
            acc *= self.gtab[i][t % 256]
            t >>= 8
            cost.muls += 1

        assert acc == expected
        return acc

    def eval(self, alpha, cost):
        delta = alpha
        s = 0
        while delta != 1:
            # find(delta)
            mu = delta
            i = 0
            while mu != self.minus1:
                mu *= mu
                cost.sqrs += 1
                i += 1
            # end find

            k = self.n-1-i
            s += 1<<k
            if i > 0:
                delta *= self.g_to_power_of_2(k)
                cost.muls += 1
            else:
                delta = -delta

        return s

    def sarkar_sqrt(self, u):
        if DEBUG: print("u = %r" % (u,))

        # This would actually be done using the addition chain.
        v = u^((self.m-1)/2)
        cost = copy(self.base_cost)

        x = u * v^2
        cost.sqrs += 1
        cost.muls += 1
        assert x == u^self.m

        x3 = x
        x2 = x3^(1<<7)
        x1 = x2^(1<<8)
        x0 = x1^(1<<8)
        assert x0 == x^(1<<(self.n-1-8))
        assert x1 == x^(1<<(self.n-1-16))
        assert x2 == x^(1<<(self.n-1-24))
        assert x3 == x^(1<<(self.n-1-31))

        cost.sqrs += 7+8+8

        # i = 0
        s = self.eval(x0, cost)

        # i = 1
        t = s >> 8
        alpha = self.mul_by_g_to(x1, t, cost)
        s = self.eval(alpha, cost)

        # i = 2
        t = (s+t) >> 8
        alpha = self.mul_by_g_to(x2, t, cost)
        s = self.eval(alpha, cost)

        # i = 3
        t = (s+t) >> 7
        alpha = self.mul_by_g_to(x3, t, cost)
        s = self.eval(alpha, cost)

        t = (s+t) >> 1
        res = self.mul_by_g_to(u * v, t, cost)
        cost.muls += 1

        if res^2 != u:
            res = None
        assert u.is_square() == (res is not None)
        return (res, cost)


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

# see addchain.py for base costs of u^{(m-1)/2}
F_p = SqrtField(p, 5, Cost(223, 23))
F_q = SqrtField(q, 5, Cost(222, 25))


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
        #print(x)
        (_, cost) = F_p.sarkar_sqrt(x)
        total_cost += cost

    print total_cost.divide(iters)
