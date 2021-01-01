#!/usr/bin/env python3

# Addition chains needed for calculation of square roots on the Pasta fields:
#
#   r = (m_p - 1)/2 for p = 1 + m_p * 2^32
#   s = (m_q - 1)/2 for q = 1 + m_q * 2^32

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


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001
n = 32
assert(p % (1<<n) == 1)
assert(q % (1<<n) == 1)

# print(format(p >> (n+1), 'b'))

r = (1<<221) + 0b100010010001101001100011111100000010010100110011111001000110111001100100101101001100001110110
#                a      b      c   d     e  f         g    h   i  j   k    l   m   n     o    p   q      r  st

assert(r == p >> (n+1))
rch = Chain()
r1    = 1
r10   = rch.sqr(r1)
r11   = rch.mul(r10, r1)
r110  = rch.sqr(r11)
r111  = rch.mul(r110, r1)
r1001 = rch.mul(r111, r10)
r1101 = rch.mul(r111, r110)
ra    = rch.sqrmul(r1, 129, r1)
rb    = rch.sqrmul(ra,  7, r1001)
rc    = rch.sqrmul(rb,  7, r1101)
rd    = rch.sqrmul(rc,  4, r11)
re    = rch.sqrmul(rd,  6, r111)
rf    = rch.sqrmul(re,  3, r111)
rg    = rch.sqrmul(rf, 10, r1001)
rh    = rch.sqrmul(rg,  5, r1001)
ri    = rch.sqrmul(rh,  4, r1001)
rj    = rch.sqrmul(ri,  3, r111)
rk    = rch.sqrmul(rj,  4, r1001)
rl    = rch.sqrmul(rk,  5, r11)
rm    = rch.sqrmul(rl,  4, r111)
rn    = rch.sqrmul(rm,  4, r11)
ro    = rch.sqrmul(rn,  6, r1001)
rp    = rch.sqrmul(ro,  5, r1101)
rq    = rch.sqrmul(rp,  4, r11)
rr    = rch.sqrmul(rq,  7, r111)
rs    = rch.sqrmul(rr,  3, r11)
rt    = rch.sqr(rs)
assert rt == r, format(rt, 'b')
print(rch)


# print(format(q >> (n+1), 'b'))

s = (1<<221) + 0b100010010001101001100011111100000010011001010010101000110111011000110001000110111010110010000
#                a      b      c   d     e  f         g   h    i    j  k   l   m    n   o    p   q    r  s   t
#                                                                1001
#                                                                   1001

assert(s == q >> (n+1))
sch = Chain()
s1    = 1
s10   = sch.sqr(s1)
s11   = sch.mul(s10, s1)
s111  = sch.sqrmul(s11, 1, s1)
s1001 = sch.mul(s111, s10)
s1011 = sch.mul(s1001, s10)
s1101 = sch.mul(s1011, s10)
sa    = sch.sqrmul(s1, 129, s1)
sb    = sch.sqrmul(sa,  7, s1001)
sc    = sch.sqrmul(sb,  7, s1101)
sd    = sch.sqrmul(sc,  4, s11)
se    = sch.sqrmul(sd,  6, s111)
sf    = sch.sqrmul(se,  3, s111)
sg    = sch.sqrmul(sf, 10, s1001)
sh    = sch.sqrmul(sg,  4, s1001)
si    = sch.sqrmul(sh,  5, s1001)
sj    = sch.sqrmul(si,  5, s1001)
sk    = sch.sqrmul(sj,  3, s1001)
sl    = sch.sqrmul(sk,  4, s1011)
sm    = sch.sqrmul(sl,  4, s1011)
sn    = sch.sqrmul(sm,  5, s11)
so    = sch.sqrmul(sn,  4, s1)
sp    = sch.sqrmul(so,  5, s11)
sq    = sch.sqrmul(sp,  4, s111)
sr    = sch.sqrmul(sq,  5, s1011)
ss    = sch.sqrmul(sr,  3, s1)
st    = sch.sqr(ss, 4)
assert st == s, format(st, 'b')
print(sch)


t = (1<<32) - 1

assert(s == q >> (n+1))
tch = Chain()
t1  = 1
t2  = tch.sqrmul(t1, 1, t1)
t4  = tch.sqrmul(t2, 2, t2)
t8  = tch.sqrmul(t4, 4, t4)
t16 = tch.sqrmul(t8, 8, t8)
t32 = tch.sqrmul(t16, 16, t16)
assert t32 == t, format(t32, 'b')
print(tch)
