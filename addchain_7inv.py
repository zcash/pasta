#!/usr/bin/env python3

class Chain(object):
    SQR_COST = 0.8
    MUL_COST = 1

    def __init__(self):
        self.muls = 0
        self.sqrs = 0
        self.computed = set((1,))

    def sqr(self, x, n=1):
        for i in range(n):
            assert(x in self.computed)
            self.sqrs += 1
            x = 2*x
            self.computed.add(x)

        return x

    def mul(self, x, y):
        assert(x in self.computed)
        assert(y in self.computed)
        assert(x != y)
        self.muls += 1
        r = x+y
        self.computed.add(r)
        return r

    def sqrmul(self, x, n, y):
        return self.mul(self.sqr(x, n), y)

    def cost(self):
        return self.sqrs*self.SQR_COST + self.muls*self.MUL_COST

    def __repr__(self):
        return "%dS + %dM (%.1f)" % (self.sqrs, self.muls, self.cost())


pr48 = 0b101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101
x_p  = 0b10110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110111001111010101101111111110001111011101000101101110001101010111001101101100100000010001111000010010110110110110110110110110110111
#                                                                                                       |                                  a       b   cx   d  e      f   g h      i   j      k    l   m     n   o      p      q    r                         s     t u
#                                                                                                                                         11    1111 111

pc = Chain()
p1    = 1
p10   = pc.sqr(p1)
p11   = pc.mul(p10, p1)
p101  = pc.mul(p11, p10)
p111  = pc.mul(p101, p10)
p1001 = pc.mul(p111, p10)
p1011 = pc.mul(p1001, p10)
p1101 = pc.mul(p1011, p10)
p1111 = pc.mul(p1101, p10)
pr2   = pc.sqrmul(p1011, 2, p1)
pr4   = pc.sqrmul(pr2, 6, pr2)
pr8   = pc.sqrmul(pr4, 12, pr4)
pr16  = pc.sqrmul(pr8, 24, pr8)
pr32  = pc.sqrmul(pr16, 48, pr16)
pr42a = pc.sqrmul(pr32, 35, p11)
pr42b = pc.sqrmul(pr42a, 8, p1111)
pr42c = pc.sqrmul(pr42b, 4, p111)
pr42x = pc.sqrmul(pr42c, 1, pr16)
pr42d = pc.sqrmul(pr42x, 4, p1111)
pr42e = pc.sqrmul(pr42d, 3, p111)
pr42f = pc.sqrmul(pr42e, 7, p1111)
pr42g = pc.sqrmul(pr42f, 4, p111)
pr42h = pc.sqrmul(pr42g, 2, p1)
pr42i = pc.sqrmul(pr42h, 7, p1011)
pr42j = pc.sqrmul(pr42i, 4, p111)
pr42k = pc.sqrmul(pr42j, 7, p1101)
pr42l = pc.sqrmul(pr42k, 5, p1011)
pr42m = pc.sqrmul(pr42l, 4, p1001)
pr42n = pc.sqrmul(pr42m, 6, pr2)
pr42o = pc.sqrmul(pr42n, 4, p1001)
pr42p = pc.sqrmul(pr42o, 7, p1)
pr42q = pc.sqrmul(pr42p, 7, p1111)
pr42r = pc.sqrmul(pr42q, 5, p1)
pr42s = pc.sqrmul(pr42r, 26, pr8)
pr42t = pc.sqrmul(pr42s, 6, pr2)
pr42u = pc.sqrmul(pr42t, 2, p11)
assert pr42u == x_p, format(pr42u, 'b')
print(pc)


qr48 = 0b101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101101
x_q  = 0b10110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110111001111010101101111111110001111011101001000111011000001110000101101000111101001100000110110000010110110110110110110110110110111
#                                                                                                       |                                  a       b   cx   d  e      f   g    h     i  j       k         l      m n   o      p  q                            r     s t
#                                                                                                                                         11    1111 111

qc = Chain()
q1    = 1
q10   = qc.sqr(q1)
q11   = qc.mul(q10, q1)
q101  = qc.mul(q11, q10)
q111  = qc.mul(q101, q10)
q1001 = qc.mul(q111, q10)
q1111 = qc.sqrmul(q111, 1, q1)
qr2   = qc.sqrmul(q1001, 2, q1001)
qr4   = qc.sqrmul(qr2, 6, qr2)
qr8   = qc.sqrmul(qr4, 12, qr4)
qr16  = qc.sqrmul(qr8, 24, qr8)
qr32  = qc.sqrmul(qr16, 48, qr16)
qr42a = qc.sqrmul(qr32, 35, q11)
qr42b = qc.sqrmul(qr42a, 8, q1111)
qr42c = qc.sqrmul(qr42b, 4, q111)
qr42x = qc.sqrmul(qr42c, 1, qr16)
qr42d = qc.sqrmul(qr42x, 4, q1111)
qr42e = qc.sqrmul(qr42d, 3, q111)
qr42f = qc.sqrmul(qr42e, 7, q1111)
qr42g = qc.sqrmul(qr42f, 4, q111)
# diverges here
qr42h = qc.sqrmul(qr42g, 5, q1001)
qr42i = qc.sqrmul(qr42h, 6, q111)
qr42j = qc.sqrmul(qr42i, 3, q11)
qr42k = qc.sqrmul(qr42j, 8, q111)
qr42l = qc.sqrmul(qr42k, 10, qr2)
qr42m = qc.sqrmul(qr42l, 7, q1111)
qr42n = qc.sqrmul(qr42m, 2, q1)
qr42o = qc.sqrmul(qr42n, 4, q11)
qr42p = qc.sqrmul(qr42o, 7, q11)
qr42q = qc.sqrmul(qr42p, 3, q11)
qr42r = qc.sqrmul(qr42q, 29, qr8)
qr42s = qc.sqrmul(qr42r, 6, qr2)
qr42t = qc.sqrmul(qr42s, 2, q11)
assert qr42t == x_q, format(qr42t, 'b')
print(qc)

