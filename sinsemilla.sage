#!/usr/bin/env sage

load('hashtocurve.sage')

VERBOSE = False

k = 10
c = len(format((q-1)//2, 'b'))-1
assert c == 253

def grouphash(D, M):
    (P, _, _) = hash_to_pallas_jacobian(M, D + b"-pallas_XMD:BLAKE2b_SSWU_RO_")
    return P

S = [grouphash(b"z.cash:SinsemillaS", pack('<I', j)) for j in range(0, 1<<k)]

def pad(n, M):
    assert n % k == 0
    m = ""
    for b in M:
        m += format(as_byte(b), '08b')[::-1]

    r = [int(m[i:i+k][::-1], 2) for i in range(0, n, k)]
    return r

def sinsemilla_hash(D, m):
    Acc = grouphash(b"z.cash:SinsemillaQ", D)
    if VERBOSE: print(Acc)
    for x in m:
        Acc = (Acc + S[x]) + Acc
        if VERBOSE: print(Acc)

    return Acc

print(sinsemilla_hash(b"z.cash:test-Sinsemilla", pad(40, b"hello")))
