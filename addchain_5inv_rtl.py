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


#                                                                                                                                 ---> up to here is a multiple of 0b110011 = 51 :-)
x_p = 0b11001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001101001110100111101110000011001001101000010000101001100000111000101110000011110000111100111111000011001100110011001100110011001101

pchain = Chain()
pi = pa = 1
for i in range(1, 128):
    pi = pchain.sqr(pi)
    if '01001110100111101110000011001001101000010000101001100000111000101110000011110000111100111111000011001100110011001100110011001101'[127-i] == '1':
        pa = pchain.mul(pa, pi)

b                                = '1000000010000000100000001000000010000000100000001000000010000000100000001000000010000000100000001000000010000000100000001'
assert int(b, 2)*0b110011 == 0b110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011
pj = pchain.sqrmul(pi, 1, pi)
pk = pchain.sqrmul(pj, 4, pj)
for k in range(1, 122):
    pk = pchain.sqr(pk)
    if b[121-k] == '1':
        pa = pchain.mul(pa, pk)

assert pa == x_p, "\n" + format(pa, '0255b') + "\n" + format(x_p, '0255b')
print(pchain)

#                                                                                                                                 ---> up to here is a multiple of 0b110011 = 51 :-)
x_q = 0b11001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001101001110100111101110000011001001101000010100001110111010010010101101011010011111001000101000000011001100110011001100110011001101

qchain = Chain()
qi = qa = 1
for i in range(1, 128):
    qi = qchain.sqr(qi)
    if '01001110100111101110000011001001101000010100001110111010010010101101011010011111001000101000000011001100110011001100110011001101'[127-i] == '1':
        qa = qchain.mul(qa, qi)

b                                = '1000000010000000100000001000000010000000100000001000000010000000100000001000000010000000100000001000000010000000100000001'
assert int(b, 2)*0b110011 == 0b110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011
qj = qchain.sqrmul(qi, 1, qi)
qk = qchain.sqrmul(qj, 4, qj)
for k in range(1, 122):
    qk = qchain.sqr(qk)
    if b[121-k] == '1':
        qa = qchain.mul(qa, qk)

assert qa == x_q, "\n" + format(qa, '0255b') + "\n" + format(x_q, '0255b')
print(qchain)
