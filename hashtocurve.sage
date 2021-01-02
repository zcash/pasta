#!/usr/bin/env sage

# Simplified SWU for a = 0 as described in [WB2019] <https://eprint.iacr.org/2019/403> and
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

load('squareroottab.sage')

DEBUG = False

# E: a short Weierstrass elliptic curve
def find_z_sswu(E):
    (0, 0, 0, A, B) = E.a_invariants()
    F = E.base_field()

    R.<x> = F[]                       # Polynomial ring over F
    g = x^3 + F(A) * x + F(B)         # y^2 = g(x) = x^3 + A * x + B
    ctr = F.gen()
    while True:
        for Z_cand in (F(ctr), F(-ctr)):
            if is_good_Z(F, g, A, B, Z_cand):
                return Z_cand
        ctr += 1

def is_good_Z(F, g, A, B, Z):
    # Criterion 1: Z is non-square in F.
    if Z.is_square():
        return False

    # Criterion 2: Z != -1 in F.
    if Z == F(-1):
        return False

    # Criterion 3: g(x) - Z is irreducible over F.
    if not (g - Z).is_irreducible():
        return False

    # Criterion 4: g(B / (Z * A)) is square in F.
    if not g(F(B) / (Z * F(A))).is_square():
        return False

    return True


# <https://core.ac.uk/download/pdf/82012348.pdf>
class ChudnovskyPoint:
    def __init__(self, x, y, z, z2, z3):
        (self.x, self.y, self.z, self.z2, self.z3) = (x, y, z, z2, z3)

    def add(self, other, E, c):
        # Addition on y^2 = x^3 + Ax + B with Chudnovsky input and output.
        # FIXME: should be unified addition.
        # <https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-1986-cc-2>
        (X1, Y1, Z1, Z1_2, Z1_3) = ( self.x,  self.y,  self.z,  self.z2,  self.z3)
        (X2, Y2, Z2, Z2_2, Z2_3) = (other.x, other.y, other.z, other.z2, other.z3)

        U1 = c.mul(X1, Z2_2)
        U2 = c.mul(X2, Z1_2)
        S1 = c.mul(Y1, Z2_3)
        S2 = c.mul(Y2, Z1_3)
        P = U2 - U1
        R = S2 - S1
        P_2 = c.sqr(P)
        P_3 = c.mul(P_2, P)
        R_2 = c.sqr(R)
        T = U1 + U2
        X3 = R_2 - c.mul(T, P_2)
        Y3 = (c.mul(R, -2*R_2 + c.mul(3*P_2, T)) - c.mul(P_3, S1 + S2))/2
        Z3 = c.mul(c.mul(Z1, Z2), P)
        Z3_2 = c.sqr(Z3)
        Z3_3 = c.mul(Z3_2, Z3)
        R = ChudnovskyPoint(X3, Y3, Z3, Z3_2, Z3_3)

        if DEBUG: assert R.to_sage(E) == self.to_sage(E) + other.to_sage(E)
        return R

    def to_sage(self, E):
        return E((self.x / self.z2, self.y / self.z3))

    def to_jacobian(self):
        return (self.x, self.y, self.z)

    def __repr__(self):
        return "%r : %r : %r : %r : %r" % (self.x, self.y, self.z, self.z2, self.z3)


assert p == 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
assert q == 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001
Fp = GF(p)
Fq = GF(q)

E_isop_A = 10949663248450308183708987909873589833737836120165333298109615750520499732811
E_isoq_A = 17413348858408915339762682399132325137863850198379221683097628341577494210225
E_isop_B = 1265
E_isoq_B = 1265
E_isop   = EllipticCurve(Fp, [E_isop_A, E_isop_B])
E_isoq   = EllipticCurve(Fq, [E_isoq_A, E_isoq_B])
E_p      = EllipticCurve(Fp, [0, 5])
E_q      = EllipticCurve(Fq, [0, 5])

k = 128
Lp = (len(format(p, 'b')) + k + 7) // 8
Lq = (len(format(q, 'b')) + k + 7) // 8
assert Lp == 48 and Lq == 48
L = Lp

Z_isop = find_z_sswu(E_isop)
Z_isoq = find_z_sswu(E_isoq)
assert Z_isop == Mod(-13, p)
assert Z_isoq == Mod(-13, q)

h_p = F_p.g
h_q = F_q.g


def select_z_nz(s, ifz, ifnz):
    # This should be constant-time in a real implementation.
    return ifz if (s == 0) else ifnz

def map_to_curve_simple_swu(F, E, Z, u, c):
    # would be precomputed
    h = F.g
    (0, 0, 0, A, B) = E.a_invariants()
    mBdivA = -B / A
    BdivZA = B / (Z * A)
    Z2 = Z^2
    assert (Z/h).is_square()
    theta = sqrt(Z/h)

    # 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
    # 2. x1 = (-B / A) * (1 + tv1)
    # 3. If tv1 == 0, set x1 = B / (Z * A)
    # 4. gx1 = x1^3 + A * x1 + B
    #
    # We use the "Avoiding inversions" optimization in [WB2019, section 4.2]
    # (not to be confused with section 4.3):
    #
    #   here       [WB2019]
    #   -------    ---------------------------------
    #   Z          \xi
    #   u          t
    #   Z * u^2    \xi * t^2 (called u, confusingly)
    #   x1         X_0(t)
    #   x2         X_1(t)
    #   gx1        g(X_0(t))
    #   gx2        g(X_1(t))
    #
    # Using the "here" names:
    #    x1 = N_x1/D = [B*(Z^2 * u^4 + Z * u^2 + 1)] / [-A*(Z^2 * u^4 + Z * u^2]
    #   gx1 = U/V    = [N_x1^3 + A * N_x1 * D^2 + B * D^3] / D^3

    # Z and B are small so we don't count multiplication by them as a mul; A is large.
    Zu2 = Z * c.sqr(u)
    ta = c.sqr(Zu2) + Zu2
    N_x1 = B * (ta + 1)
    D = c.mul(-A, ta)
    N2_x1 = c.sqr(N_x1)
    D2 = c.sqr(D)
    D3 = c.mul(D2, D)
    U = select_z_nz(ta, BdivZA, c.mul(N2_x1 + c.mul(A, D2), N_x1) + B*D3)
    V = select_z_nz(ta, 1, D3)

    if DEBUG:
        x1 = N_x1/D
        gx1 = U/V
        tv1 = (0 if ta == 0 else 1/ta)
        assert x1 == (BdivZA if tv1 == 0 else mBdivA * (1 + tv1))
        assert gx1 == x1^3 + A * x1 + B

    # 5. x2 = Z * u^2 * x1
    N_x2 = c.mul(Zu2, N_x1)  # same D

    # 6. gx2 = x2^3 + A * x2 + B  [optimized out; see below]
    # 7. If is_square(gx1), set x = x1 and y = sqrt(gx1)
    # 8. Else set x = x2 and y = sqrt(gx2)
    (y1, zero_if_gx1_square) = F.sarkar_divsqrt(U, V, c)

    # This magic also comes from a generalization of [WB2019, section 4.2].
    #
    # The Sarkar square root algorithm with input s gives us a square root of
    # h * s for free when s is not square, where h is a fixed nonsquare.
    # We know that Z/h is a square since both Z and h are nonsquares.
    # Precompute \theta as a square root of Z/h, or choose Z = h so that \theta = 1.
    #
    # We have gx2 = g(Z * u^2 * x1) = Z^3 * u^6 * gx1
    #                               = (Z * u^3)^2 * (Z/h * h * gx1)
    #                               = (Z * \theta * u^3)^2 * (h * gx1)
    #
    # When gx1 is not square, y1 is a square root of h * gx1, and so Z * \theta * u^3 * y1
    # is a square root of gx2. Note that we don't actually need to compute gx2.

    y2 = c.mul(theta, c.mul(Zu2, c.mul(u, y1)))
    if DEBUG and zero_if_gx1_square != 0:
        assert y1^2 == h * gx1, (y1_2, Z, gx1)
        assert y2^2 == x2^3 + A * x2 + B, (y2, x2, A, B)

    N_x = select_z_nz(zero_if_gx1_square, N_x1, N_x2)
    y = select_z_nz(zero_if_gx1_square, y1, y2)

    # 9. If sgn0(u) != sgn0(y), set y = -y
    y = select_z_nz((int(u) % 2) - (int(y) % 2), y, -y)

    return ChudnovskyPoint(c.mul(N_x, D), c.mul(y, D3), D, D2, D3)


# iso_Ep = Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 10949663248450308183708987909873589833737836120165333298109615750520499732811*x + 1265 over Fp
def isop_map_affine(x, y, c):
    c.muls += 2+1+1 + 2+1+1+2
    # batch inversion
    c.muls += 3
    c.invs += 1
    Nx = ((( 6432893846517566412420610278260439325191790329320346825767705947633326140075 *x +
            23989696149150192365340222745168215001509815558210986772351135915822265203574)*x +
            10492611921771203378452795982353351666191589197598957448093274638589204800759)*x +
            12865787693035132824841220556520878650383580658640693651535411895266652280192)
    Dx =  ((                                                                               x +
            13271109177048389296812780941310096270046944650307955939477485891950613419807)*x +
            22768321103861051515190775253992702316905399997697804654926324362758820947460)

    Ny = (((11793638718615538422771118843477472096184948937087302513907460903994431256804 *x +
            11994848074575096182670111372584107500754907779105493386175567957911132601787)*x +
            28823569610051396102362669851238297121581474897215657071023781420043761726004)*x +
             1072148974419594402070101713043406554198631721553391137627950991272221023311) * y
    Dy = (((                                                                               x +
             5432652610908059517272798285879155923388888734491153551238890455750936314542)*x +
            10408918692925056833786833257634153023990087029210292532869619559576527581706)*x +
            28948022309329048855892746252171976963363056481941560715954676764349967629797)

    return (Nx / Dx, Ny / Dy)

# The same isogeny but with input in Chudnovsky coordinates (Jacobian plus z^2 and z^3)
# and output in Jacobian <https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html>,
# according to "Avoiding inversions" in [WB2019, section 4.3].
def isop_map_jacobian(P, c):
    (x, y, z, z2, z3) = (P.x, P.y, P.z, P.z2, P.z3)
    z4 = c.sqr(z2)
    z6 = c.sqr(z3)

    Nx = ((( 6432893846517566412420610278260439325191790329320346825767705947633326140075    *x +
            23989696149150192365340222745168215001509815558210986772351135915822265203574*z2)*x +
            10492611921771203378452795982353351666191589197598957448093274638589204800759*z4)*x +
            12865787693035132824841220556520878650383580658640693651535411895266652280192*z6)
    c.muls += 6
    Dx =  ((                                                                              z2 *x +
            13271109177048389296812780941310096270046944650307955939477485891950613419807*z4)*x +
            22768321103861051515190775253992702316905399997697804654926324362758820947460*z6)
    c.muls += 4

    Ny = (((11793638718615538422771118843477472096184948937087302513907460903994431256804    *x +
            11994848074575096182670111372584107500754907779105493386175567957911132601787*z2)*x +
            28823569610051396102362669851238297121581474897215657071023781420043761726004*z4)*x +
             1072148974419594402070101713043406554198631721553391137627950991272221023311*z6) * y
    c.muls += 7
    Dy = (((                                                                                  x +
             5432652610908059517272798285879155923388888734491153551238890455750936314542*z2)*x +
            10408918692925056833786833257634153023990087029210292532869619559576527581706*z4)*x +
            28948022309329048855892746252171976963363056481941560715954676764349967629797*z6) * z3
    c.muls += 6

    zo = c.mul(Dx, Dy)
    xo = c.mul(c.mul(Nx, Dy), zo)
    yo = c.mul(c.mul(Ny, Dx), c.sqr(zo))

    assert isop_map_affine(x / z2, y / z3, Cost()) == (xo / zo^2, yo / zo^3)
    return (xo, yo, zo)


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


def hash_to_curve_jacobian(msg, DST):
    c = Cost()
    us = hash_to_field(msg, DST, 2)
    #print("u = ", u)
    Q0 = map_to_curve_simple_swu(F_p, E_isop, Z_isop, us[0], c)
    Q1 = map_to_curve_simple_swu(F_p, E_isop, Z_isop, us[1], c)

    R = Q0.add(Q1, E_isop, c)

    # no cofactor clearing needed since Pallas and Vesta are prime-order
    (Px, Py, Pz) = isop_map_jacobian(R, c)
    P = E_p((Px / Pz^2, Py / Pz^3))
    return (P, c)


print(hash_to_curve_jacobian("hello", "blah"))
print("")

iters = 100
for i in range(iters):
    (R, cost) = hash_to_curve_jacobian(pack(">I", i), "blah")
    print(R, cost)
