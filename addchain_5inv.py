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


x_p = 0b11001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001101001110100111101110000011001001101000010000101001100000111000101110000011110000111100111111000011001100110011001100110011001101
#                                                                                                                                         a       b   c d      e     f  g    h      i   j       k   l   m        n       o     p q                                 r s
#                                                                                                                                          11001100110011
#                                                                                                                                               111   1

pc = Chain()
p1    = 1
p10   = pc.sqr(p1)
p11   = pc.mul(p10, p1)
p101  = pc.mul(p10, p11)
p110  = pc.sqr(p11)
p111  = pc.mul(p110, p1)
p1001 = pc.mul(p111, p10)
p1111 = pc.mul(p1001, p110)
pr2   = pc.sqrmul(p110, 3, p11)
pr4   = pc.sqrmul(pr2, 8, pr2)
pr8   = pc.sqrmul(pr4, 16, pr4)
pr16  = pc.sqrmul(pr8, 32, pr8)
pr32  = pc.sqrmul(pr16, 64, pr16)
pr32a = pc.sqrmul(pr32, 5, p1001)
pr32b = pc.sqrmul(pr32a, 8, p111)
pr32c = pc.sqrmul(pr32b, 4, p1)
pr32d = pc.sqrmul(pr32c, 2, pr4)
pr32e = pc.sqrmul(pr32d, 7, p11)
pr32f = pc.sqrmul(pr32e, 6, p1001)
pr32g = pc.sqrmul(pr32f, 3, p101)
pr32h = pc.sqrmul(pr32g, 5, p1)
pr32i = pc.sqrmul(pr32h, 7, p101)
pr32j = pc.sqrmul(pr32i, 4, p11)
pr32k = pc.sqrmul(pr32j, 8, p111)
pr32l = pc.sqrmul(pr32k, 4, p1)
pr32m = pc.sqrmul(pr32l, 4, p111)
pr32n = pc.sqrmul(pr32m, 9, p1111)
pr32o = pc.sqrmul(pr32n, 8, p1111)
pr32p = pc.sqrmul(pr32o, 6, p1111)
pr32q = pc.sqrmul(pr32p, 2, p11)
pr32r = pc.sqrmul(pr32q, 34, pr8)
pr32s = pc.sqrmul(pr32r, 2, p1)
assert pr32s == x_p
print(pc)

x_q = 0b11001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001101001110100111101110000011001001101000010100001110111010010010101101011010011111001000101000000011001100110011001100110011001101
#                                                                                                                                         a       b   c d      e     f  g      h      i   j    k    l  m   n  o     p   q     r                                    s t
#                                                                                                                                          11001100110011
#                                                                                                                                               111   1

qc = Chain()
q1    = 1
q10   = qc.sqr(q1)
q11   = qc.mul(q10, q1)
q101  = qc.mul(q10, q11)
q110  = qc.sqr(q11)
q111  = qc.mul(q110, q1)
q1001 = qc.mul(q111, q10)
q1111 = qc.mul(q1001, q110)
qr2   = qc.sqrmul(q110, 3, q11)
qr4   = qc.sqrmul(qr2, 8, qr2)
qr8   = qc.sqrmul(qr4, 16, qr4)
qr16  = qc.sqrmul(qr8, 32, qr8)
qr32  = qc.sqrmul(qr16, 64, qr16)
qr32a = qc.sqrmul(qr32, 5, q1001)
qr32b = qc.sqrmul(qr32a, 8, q111)
qr32c = qc.sqrmul(qr32b, 4, q1)
qr32d = qc.sqrmul(qr32c, 2, qr4)
qr32e = qc.sqrmul(qr32d, 7, q11)
qr32f = qc.sqrmul(qr32e, 6, q1001)
qr32g = qc.sqrmul(qr32f, 3, q101)
# diverges here
qr32h = qc.sqrmul(qr32g, 7, q101)
qr32i = qc.sqrmul(qr32h, 7, q111)
qr32j = qc.sqrmul(qr32i, 4, q111)
qr32k = qc.sqrmul(qr32j, 5, q1001)
qr32l = qc.sqrmul(qr32k, 5, q101)
qr32m = qc.sqrmul(qr32l, 3, q11)
qr32n = qc.sqrmul(qr32m, 4, q101)
qr32o = qc.sqrmul(qr32n, 3, q101)
qr32p = qc.sqrmul(qr32o, 6, q1111)
qr32q = qc.sqrmul(qr32p, 4, q1001)
qr32r = qc.sqrmul(qr32q, 6, q101)
qr32s = qc.sqrmul(qr32r, 37, qr8)
qr32t = qc.sqrmul(qr32s, 2, q1)
assert qr32t == x_q
print(qc)

