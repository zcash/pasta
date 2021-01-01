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

DEBUG = True

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

def map_to_curve_simple_swu(F, E, Z, us, c):
    # would be precomputed
    h = F.g
    (0, 0, 0, A, B) = E.a_invariants()
    mBdivA = -B / A
    BdivZA = B / (Z * A)
    Z2 = Z^2
    assert (Z/h).is_square()
    theta = sqrt(Z/h)

    Qs = []
    for u in us:
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
        #   gx1        g(X_0(t))
        #   gx2        g(X_1(t))
        #
        #   X0(u)  = N/D = [B*(Z^2 * u^4 + Z * u^2 + 1)] / [-A*(Z^2 * u^4 + Z * u^2]
        # g(X0(u)) = U/V = [N^3 + A * N * D^2 + B * D^3] / D^3

        Zu2 = Z * c.sqr(u)  # Z is small
        ta = c.sqr(Zu2) + Zu2
        N = c.mul(B, ta + 1)
        D = c.mul(-A, ta)
        N2 = c.sqr(N)
        D2 = c.sqr(D)
        D3 = c.mul(D2, D)
        U = select_z_nz(ta, BdivZA, c.mul(N2 + A*D2, N) + B*D3)
        V = select_z_nz(ta, 1, D3)

        if DEBUG:
            x1 = N/D
            gx1 = U/V
            tv1 = (0 if ta == 0 else 1/ta)
            assert x1 == (BdivZA if tv1 == 0 else mBdivA * (1 + tv1))
            assert gx1 == x1^3 + A * x1 + B

        # 5. x2 = Z * u^2 * x1
        x2 = c.mul(Zu2, x1)

        # 6. gx2 = x2^3 + A * x2 + B  [optimized out; see below]
        # 7. If is_square(gx1), set x = x1 and y = sqrt(gx1)
        # 8. Else set x = x2 and y = sqrt(gx2)
        (y1, zero_if_gx1_square) = F.sarkar_divsqrt(U, V, c)

        # This magic also comes from a generalization of [WB2019, section 4.2].
        #
        # The Sarkar square root algorithm with input s gives us a square root of
        # h * s for free when s is not square, provided we choose h to be a generator
        # of the order 2^n multiplicative subgroup (where n = 32 for Pallas and Vesta).
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

        x = select_z_nz(zero_if_gx1_square, x1, x2)
        y = select_z_nz(zero_if_gx1_square, y1, y2)

        # 9. If sgn0(u) != sgn0(y), set y = -y
        y = select_z_nz((int(u) % 2) - (int(y) % 2), y, -y)
        Qs.append(E((x, y)))

    return Qs


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

# The same isogeny but in Jacobian coordinates <https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html>,
# according to "Avoiding inversions" in [WB2019, section 4.3].
def isop_map_jacobian(x, y, z, c):
    z2 = c.sqr(z)
    z3 = c.mul(z, z2)
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


# Unified addition on y^2 = x^3 + Ax + B with affine input and Jacobian output.
# The inputs must not be the point at infinity; the output may be.
def unified_mmadd_jacobian(A, Px, Py, Qx, Qy, c):
    # Addition using Jacobian coordinates for general A
    # <https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-mmadd-2007-bl>
    H = Qx - Px
    I = 4*c.sqr(H)
    J = c.mul(H, I)
    r = 2*(Qy - Py)
    V = c.mul(Px, I)

    # Doubling using Jacobian coordinates for general A
    # <https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-mdbl-2007-bl>
    XX = c.sqr(Px)
    YY = c.sqr(Py)
    YYYY = c.sqr(YY)
    S = 2*(c.sqr(Px + YY) - XX - YYYY)
    M = 3*XX + A

    # Common part between doubling and addition. J = 0 for doubling.
    M_or_r = select_z_nz(H, M, r)
    S_or_V = select_z_nz(H, S, V)
    Rx = c.sqr(M_or_r) - J - 2*S_or_V
    Ry = c.mul(M_or_r, S_or_V - Rx) - select_z_nz(H, 8*YYYY, 2*c.mul(Py, J))

    # If Q = -P (i.e. H = 0 and Py + Qy = 0), then the result is the point at infinity, represented by Rz = 0.
    U = select_z_nz(Py + Qy, 0, Qy)
    Rz = 2*select_z_nz(H, U, H)

    return (Rx, Ry, Rz)


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


def hash_to_curve_affine(msg, DST, uniform=True):
    c = Cost()
    us = hash_to_field(msg, DST, 2 if uniform else 1)
    #print("u = ", u)
    Qs = map_to_curve_simple_swu(F_p, E_isop, Z_isop, us, c)

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
    P = E_p(isop_map_affine(x, y, c))
    return (P, c)

def hash_to_curve_jacobian(msg, DST):
    c = Cost()
    us = hash_to_field(msg, DST, 2)
    #print("u = ", u)
    Qs = map_to_curve_simple_swu(F_p, E_isop, Z_isop, us, c)

    R = Qs[0] + Qs[1]
    #print("R = ", R)

    (Q0x, Q0y) = Qs[0].xy()
    (Q1x, Q1y) = Qs[1].xy()
    (Rx, Ry, Rz) = unified_mmadd_jacobian(E_isop_A, Q0x, Q0y, Q1x, Q1y, c)
    assert E_isop((Rx / Rz^2, Ry / Rz^3)) == R

    # no cofactor clearing needed since Pallas and Vesta are prime-order
    (Px, Py, Pz) = isop_map_jacobian(Rx, Ry, Rz, c)
    P = E_p((Px / Pz^2, Py / Pz^3))
    return (P, c)


print(hash_to_curve_affine("hello", "blah", uniform=True))
print(hash_to_curve_jacobian("hello", "blah"))
print("")

iters = 100
for i in range(iters):
    (R_affine, cost_affine) = hash_to_curve_affine(pack(">I", i), "blah", uniform=True)
    (R_jacobian, cost_jacobian) = hash_to_curve_jacobian(pack(">I", i), "blah")
    assert R_affine == R_jacobian  # Sage will normalize them
    print(R_affine, cost_affine, cost_jacobian)
