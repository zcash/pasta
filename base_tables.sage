#!/usr/bin/env sage

# Let's say we want to interpolate between h curve points (x, y) over a
# curve y^2 = x^3 + b in a PLONK circuit, for small h.
# The obvious way to do it involves 2h fixed columns:
# use \sum\limits_{j=0}^{h-1} x_j . l_j to interpolate x, and similarly for y.
#
# We want to use only h+1 columns. Here's how:
# - Interpolate x as above.
# - Witness y and check y^2 = x^3 + b.
# - Witness u such that u^2 = y+z.
#
# where z is some field element that "makes the signs come out right".
# The purpose of this script is to find z.

load('hashtocurve.sage')

if sys.version_info[0] == 2:
    from string import maketrans
else:
    maketrans = str.maketrans


def hash_to_pallas(domain_prefix, msg):
    (P, _, _) = hash_to_pallas_jacobian(msg, domain_prefix + "-pallas_XMD:BLAKE2b_SSWU_RO_")
    return P

def I2LEOSP_32(j):
    return pack("<I", j)

TO_UNDERSCORES = maketrans("-", "_")

def to_identifier(s):
    return s.translate(TO_UNDERSCORES).upper()

def gen_bases():
    b = []

    # ValueCommit
    b.append( ("ORCHARD_CV_V", hash_to_pallas("z.cash:Orchard_cv", b"v")) )
    b.append( ("ORCHARD_CV_R", hash_to_pallas("z.cash:Orchard_cv", b"r")) )

    # Sinsemilla Q
    for D in ( b"z.cash:Orchard-NoteCommit-M",
               b"z.cash:Orchard-NoteCommit-r",
               b"z.cash:Orchard-CommitIvk-M",
               b"z.cash:Orchard-CommitIvk-r" ):
        b.append( ("Q_" + to_identifier(D.split(b':')[-1]), hash_to_pallas("z.cash:SinsemillaQ", D)) )

    # Sinsemilla S
    for j in range(1<<10):
        b.append( ("S_%d" % (j,), hash_to_pallas("z.cash:SinsemillaS", I2LEOSP_32(j))) )

    #print(b)
    return b

bases = gen_bases()

assert p == 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
h = 8

def find_z(ps):
    for z in range(1, 1000*(1<<(2*h))):
        if sum([(y+z).is_square() and not (-y+z).is_square() for (x, y) in ps]) == h:
            return z

    print("There's a glitch in the matrix.")
    return None

def test():
    for (name, p) in bases:
        ps = [(i * p).xy() for i in range(1, h+1)]
        print(name, ps, find_z(ps))

test()
