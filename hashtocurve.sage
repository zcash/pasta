#!/usr/bin/env sage

# Simplified SWU for a = 0 as described in [WB19] <https://eprint.iacr.org/2019/403> and
# <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html#name-simplified-swu-for-ab-0-2>.

import sys
from math import ceil, log
from struct import pack

import hashlib
if sys.version_info < (3, 6):
    try:
        import sha3
    except ImportError:
        print('Please run:\n`sage -c "import sys; print(sys.executable)"` -m pip install pysha3\n')
        raise
from hashlib import shake_128

if sys.version_info[0] == 2:
    range = xrange
    as_byte = ord
else:
    as_byte = lambda x: x


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

    def sqrt(self, x):
        self.sqrs += 247
        self.muls += 35
        return x.sqrt() if x.is_square() else 0

    def __add__(self, other):
        return Cost(self.sqrs + other.sqrs, self.muls + other.muls, self.invs + other.invs)

    def __repr__(self):
        return "%dS + %dM + %dI" % (self.sqrs, self.muls, self.invs)

# E: a short Weierstrass elliptic curve
def find_z_sswu(E):
    (0, 0, 0, A, B) = E.a_invariants()
    F = E.base_field()

    R.<x> = F[]                       # Polynomial ring over F
    g = x^3 + F(A) * x + F(B)         # y^2 = g(x) = x^3 + A * x + B
    ctr = F.gen()
    while True:
        for Z_cand in (F(ctr), F(-ctr)):
            if Z_cand.is_square():
                # Criterion 1: Z is non-square in F.
                continue
            if Z_cand == F(-1):
                # Criterion 2: Z != -1 in F.
                continue
            if not (g - Z_cand).is_irreducible():
                # Criterion 3: g(x) - Z is irreducible over F.
                continue
            if g(B / (Z_cand * A)).is_square():
                # Criterion 4: g(B / (Z * A)) is square in F.
                return Z_cand
        ctr += 1


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
E_isop = EllipticCurve(GF(p), [10949663248450308183708987909873589833737836120165333298109615750520499732811, 1265])
E_p    = EllipticCurve(GF(p), [0, 5])
Z_isop = find_z_sswu(E_isop)
assert Z_isop == Mod(-13, p)

k = 128
L = (len(format(p, 'b')) + k + 7) // 8
assert L == 48

CONSTANT_TIME = True

def select_z_nz(s, ifz, ifnz):
    # This should be constant-time in a real implementation.
    return ifz if (s == 0) else ifnz

def map_to_curve_simple_swu(E, Z, us, c):
    # would be precomputed
    (0, 0, 0, A, B) = E.a_invariants()
    mBdivA = -B / A
    BdivZA = B / (Z * A)
    Z2 = Z^2

    # 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
    #        = inv0((Z^2 * u^2 + Z) * u^2)
    u2s = [c.sqr(u) for u in us]
    tas = [c.mul((Z2*u2 + Z), u2) for u2 in u2s]
    tv1s = c.batch_inv0(tas)

    Qs = []
    for i in range(len(us)):
        (u, u2, tv1) = (us[i], u2s[i], tv1s[i])

        # 2. x1 = (-B / A) * (1 + tv1)
        # 3. If tv1 == 0, set x1 = B / (Z * A)
        x1 = select_z_nz(tv1, BdivZA, mBdivA * (1 + tv1))

        # 4. gx1 = x1^3 + A * x1 + B
        #        = x1*(x1^2 + A) + B
        x1_2 = c.sqr(x1)
        gx1 = c.mul(x1, x1_2 + A) + B

        # 5. x2 = Z * u^2 * x1
        tb = c.mul(Z, u2)
        x2 = c.mul(tb, x1)

        # 6. gx2 = x2^3 + A * x2 + B
        #        = x2*(x2^2 + A) + B
        x2_2 = c.sqr(x2)
        gx2 = c.mul(x2, x2_2 + A) + B

        # 7. If is_square(gx1), set x = x1 and y = sqrt(gx1)
        # 8. Else set x = x2 and y = sqrt(gx2)
        y1 = c.sqrt(gx1)
        y1_2 = c.sqr(y1)
        if CONSTANT_TIME or y1_2 != gx1:
            y2 = c.sqrt(gx2)
            x = select_z_nz(y1_2 - gx1, x1, x2)
            y = select_z_nz(y1_2 - gx1, y1, y2)
        else:
            (x, y) = (x1, y1)

        # 9. If sgn0(u) != sgn0(y), set y = -y
        y = select_z_nz((int(u) % 2) - (int(y) % 2), y, -y)
        Qs.append(E((x, y)))

    return Qs


# iso_Ep = Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 10949663248450308183708987909873589833737836120165333298109615750520499732811*x + 1265 over Fp
def iso_map(x, y, c):
    c.muls += 2+1+1 + 2+1+1+2
    # batch inversion
    c.muls += 3
    c.invs += 1
    return (((( 6432893846517566412420610278260439325191790329320346825767705947633326140075 *x +
               23989696149150192365340222745168215001509815558210986772351135915822265203574)*x +
               10492611921771203378452795982353351666191589197598957448093274638589204800759)*x +
               12865787693035132824841220556520878650383580658640693651535411895266652280192) /
             ((                                                                               x +
               13271109177048389296812780941310096270046944650307955939477485891950613419807)*x +
               22768321103861051515190775253992702316905399997697804654926324362758820947460),
            (((11793638718615538422771118843477472096184948937087302513907460903994431256804 *x +
               11994848074575096182670111372584107500754907779105493386175567957911132601787)*x +
               28823569610051396102362669851238297121581474897215657071023781420043761726004)*x +
                1072148974419594402070101713043406554198631721553391137627950991272221023311) * y /
            (((                                                                               x +
                5432652610908059517272798285879155923388888734491153551238890455750936314542)*x +
               10408918692925056833786833257634153023990087029210292532869619559576527581706)*x +
               28948022309329048855892746252171976963363056481941560715954676764349967629797))


def expand_message_xof(msg, DST, len_in_bytes):
    assert len(DST) < 256
    len_in_bytes = int(len_in_bytes)

    # This is horrible but matches the reference code.
    xof = shake_128()
    xof.update(msg)
    xof.update(pack(">H", len_in_bytes))
    xof.update(pack("B", len(DST)))
    xof.update(DST)
    return xof.digest(len_in_bytes)

def hash_to_field(msg, DST, count):
    uniform_bytes = expand_message_xof(msg, DST, L*count)
    return [Mod(OS2IP(uniform_bytes[L*i : L*(i+1)]), p) for i in range(count)]

def OS2IP(bs):
    acc = 0
    for b in bs:
        acc = (acc<<8) + as_byte(b)
    return acc

def hash_to_curve(msg, DST, uniform=True):
    c = Cost()
    us = hash_to_field(msg, DST, 2 if uniform else 1)
    #print("u = ", u)
    Qs = map_to_curve_simple_swu(E_isop, Z_isop, us, c)

    if uniform:
        # Complete addition using affine coordinates: I + 2M + 2S
        # (S for x1^2; compute numerator and denominator of the division for the correct case;
        # I + M to divide; S + M to compute x and y of the result.)
        R = Qs[0] + Qs[1]
        #print("R = ", R)
        c.invs += 1
        c.sqrs += 2
        c.muls += 2
    else:
        R = Qs[0]

    # no cofactor clearing needed since Pallas and Vesta are prime-order
    (x, y) = R.xy()
    P = E_p(iso_map(x, y, c))
    return (P, c)


#print(hash_to_curve("hello", "blah"))

iters = 100
for i in range(iters):
    (res, cost) = hash_to_curve(pack(">I", i), "blah", uniform=True)
    print(res, cost)
