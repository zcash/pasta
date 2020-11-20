#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from itertools import combinations

if sys.version_info[0] == 2:
    range = xrange
    from string import maketrans
else:
    maketrans = str.maketrans


# Let Ep/Fp : y^2 = x^3 + bp
# Let Eq/Fq : y^2 = x^3 + bq

# p and q should each be ~ L bits.

DEFAULT_TWOADICITY = 32
DEFAULT_STRETCH = 0

COEFFICIENT_RANGE = (5,)
#COEFFICIENT_RANGE = range(1, 100)

GCD_PRIMES = (5, 7, 11, 13, 17)

# Set to a prime, or 0 to disable searching for isogenies.
ISOGENY_DEGREE_MAX = 3
#ISOGENY_DEGREE_MAX = 37

DEFAULT_TWIST_SECURITY = 120
REQUIRE_PRIMITIVE = True
REQUIRE_HALFZERO = True  # nearpowerof2 strategy only


# <https://cryptojedi.org/papers/pfcpo.pdf> section 2:
# [...] the order of a curve satisfying the norm equation 3V^2 = 4p - T^2 has one
# of the six forms {p+1 +/- T, p+1 +/- (T +/- 3V)/2} [IEEE Std 1363-2000, section
# A.14.2.3, item 6].
#
# We choose 4p = 3V^2 + T^2, where (V-1)/2 and (T-1)/2 are both multiples of 2^twoadicity.
#
# Then 4p = (3(V-1)^2 + 6(V-1) + 3) + ((T-1)^2 + 2(T-1) + 1)
#         = 3(V-1)^2 + 6(V-1) + (T-1)^2 + 2(T-1) + 4
#       p = 3((V-1)/2)^2 + 3(V-1)/2 + ((T-1)/2)^2 + (T-1)/2 + 1
#
# So p-1 will be a multiple of 2^twoadicity, and so will q-1 for q in
# { p + 1 - T, p + 1 + (T-3V)/2 }.
#
# We'd also like both p and q to be 1 (mod 6), so that we have efficient endomorphisms
# on both curves.

def low_hamming_order(L, twoadicity, wid, processes):
    Vlen = (L-1)//2 + 1
    Vbase = 1 << Vlen
    Tlen = (L-1)//4
    Tbase = 1 << Tlen
    trailing_zeros = twoadicity+1
    for w in range(wid, Tlen-trailing_zeros, processes):
        for Vc in combinations(range(trailing_zeros, Vlen), w):
            V = Vbase + sum([1 << i for i in Vc]) + 1
            assert(((V-1)/2) % (1<<twoadicity) == 0)
            for Tw in range(1, w+1):
                for Tc in combinations(range(trailing_zeros, Tlen), Tw):
                    T = Tbase + sum([1 << i for i in Tc]) + 1
                    assert(((T-1)/2) % (1<<twoadicity) == 0)
                    if T % 6 != 1:
                        continue
                    p4 = 3*V^2 + T^2
                    assert(p4 % 4 == 0)
                    p = p4//4
                    assert(p % (1<<twoadicity) == 1)
                    if p % 6 == 1 and is_pseudoprime(p):
                        yield (p, T, V)


def near_powerof2_order(L, twoadicity, wid, processes):
    trailing_zeros = twoadicity+1
    Vbase = isqrt((1<<(L+1))//3) >> trailing_zeros
    for Voffset in symmetric_range(100000, base=wid, step=processes):
        V = ((Vbase + Voffset) << trailing_zeros) + 1
        assert(((V-1)/2) % (1 << twoadicity) == 0)
        tmp = (1<<(L+1)) - 3*V^2
        if tmp < 0: continue
        Tbase = isqrt(tmp) >> trailing_zeros
        for Toffset in symmetric_range(100000):
            T = ((Tbase + Toffset) << trailing_zeros) + 1
            assert(((T-1)/2) % (1<<twoadicity) == 0)
            if T % 6 != 1:
                continue
            p4 = 3*V^2 + T^2
            assert(p4 % 4 == 0)
            p = p4//4
            assert(p % (1<<twoadicity) == 1)
            if REQUIRE_HALFZERO and p>>(L//2) != 1<<(L - 1 - L//2):
                continue

            if p > 1<<(L-1) and p % 6 == 1 and is_pseudoprime(p):
                yield (p, T, V)

def symmetric_range(n, base=0, step=1):
    for i in range(base, n, step):
        yield -i
        yield i+1

SWAP_SIGNS = maketrans("+-", "-+")

def find_nice_curves(strategy, L, twoadicity, stretch, requireisos, sortpq, twistsec, wid, processes):
    for (p, T, V) in strategy(L, max(0, twoadicity-stretch), wid, processes):
        if p % (1<<twoadicity) != 1:
            continue

        sys.stdout.write('.')
        sys.stdout.flush()

        for (q, qdesc) in ((p + 1 - T,          "p + 1 - T"),
                           (p + 1 + (T-3*V)//2, "p + 1 + (T-3*V)/2")):
            if strategy == near_powerof2_order and REQUIRE_HALFZERO and q>>(L//2) != 1<<(L - 1 - L//2):
                continue

            if q not in (p, p+1, p-1) and q > 1<<(L-1) and q % 6 == 1 and q % (1<<twoadicity) == 1 and is_prime(q) and is_prime(p):
                if sortpq and q < p:
                    (p, q) = (q, p)
                    qdesc = qdesc.translate(SWAP_SIGNS)

                (Ep, bp) = find_curve(p, q)
                if bp == None: continue
                (Eq, bq) = find_curve(q, p, (bp,))
                if bq == None: continue

                sys.stdout.write('*')
                sys.stdout.flush()

                primp = (Mod(bp, p).multiplicative_order() == p-1)
                if REQUIRE_PRIMITIVE and not primp: continue
                primq = (Mod(bq, q).multiplicative_order() == q-1)
                if REQUIRE_PRIMITIVE and not primq: continue

                (twsecp, twembedp) = twist_security(p, q)
                if twsecp < twistsec: continue
                (twsecq, twembedq) = twist_security(q, p)
                if twsecq < twistsec: continue

                (secp, embedp) = curve_security(p, q)
                (secq, embedq) = curve_security(q, p)

                zetap = GF(p).zeta(3)
                zetap = min(zetap, zetap^2)
                assert(zetap**3 == Mod(1, p))

                zetaq = GF(q).zeta(3)
                P = Ep.gens()[0]
                zP = endo(Ep, zetap, P)
                if zP != int(zetaq)*P:
                    zetaq = zetaq^2
                    assert(zP == int(zetaq)*P)
                assert(zetaq**3 == Mod(1, q))

                Q = Eq.gens()[0]
                assert(endo(Eq, zetaq, Q) == int(zetap)*Q)

                embeddivp = (q-1)/embedp
                embeddivq = (p-1)/embedq
                twembeddivp = (2*p + 1 - q)/twembedp
                twembeddivq = (2*q + 1 - p)/twembedq

                iso_Ep = find_iso(Ep)
                iso_Eq = find_iso(Eq)
                if requireisos and (iso_Ep is None or iso_Eq is None):
                    continue

                yield (p, q, bp, bq, zetap, zetaq, qdesc, primp, primq, secp, secq, twsecp, twsecq,
                       embeddivp, embeddivq, twembeddivp, twembeddivq, iso_Ep, iso_Eq)

def endo(E, zeta, P):
   (xp, yp) = P.xy()
   return E(zeta*xp, yp)

def find_curve(p, q, b_range=None):
    for b in (b_range or COEFFICIENT_RANGE):
        E = EllipticCurve(GF(p), [0, b])
        if E.count_points() == q:
            return (E, b)
    return (None, None)

def find_gcd_primes(p):
    return (r for r in GCD_PRIMES if gcd(p-1, r) == 1)

pi_12 = (pi/12).numerical_approx()

def curve_security(p, q):
    sys.stdout.write('!')
    sys.stdout.flush()
    r = factor(q)[-1][0]
    return (log(pi_12 * r, 4), embedding_degree(p, r))

def twist_security(p, q):
    return curve_security(p, 2*(p+1) - q)

def embedding_degree(p, r):
    sys.stdout.write('#')
    sys.stdout.flush()
    assert(gcd(p, r) == 1)
    Z_q = Integers(r)
    u = Z_q(p)
    d = r-1
    V = factor(d)
    for (v, k) in V:
        while d % v == 0:
            if u^(d/v) != 1: break
            d /= v

    return d

def find_iso(E):
    # Based on <https://eprint.iacr.org/2019/403.pdf> Appendix A.
    # Look for isogenous curves having j-invariant not in {0, 1728}.
    for degree in primes(ISOGENY_DEGREE_MAX+1):
        sys.stdout.write('~')
        sys.stdout.flush()
        for iso in E.isogenies_prime_degree(degree):
            if iso.codomain().j_invariant() not in (0, 1728):
                return iso.dual()

    return None


def format_weight(x, detail=True):
    X = format(abs(x), 'b')
    if detail:
        assert(X.endswith('1'))
        detailstr = " (bitlength %d, weight %d, 2-adicity %d)" % (len(X), sum([int(c) for c in X]),
                                                                  len(X) - len(X[:-1].rstrip('0')))
    else:
        detailstr = " (bitlength %d)" % (len(X),)

    return "%s0b%s%s" % ("-" if x < 0 else "", X, detailstr)


def main():
    args = sys.argv[1:]
    strategy = near_powerof2_order if "--nearpowerof2" in args else low_hamming_order
    processes = 1 if "--sequential" in args else cpu_count()
    if processes >= 6:
        processes -= 2
    requireisos = "--requireisos" in args
    sortpq = "--sortpq" in args
    twistsec = 0 if "--ignoretwist" in args else DEFAULT_TWIST_SECURITY
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1:
        print("""
Usage: sage amicable.sage [--sequential] [--requireisos] [--sortpq] [--ignoretwist] [--nearpowerof2] <min-bitlength> [<min-2adicity> [<stretch]]

Arguments:
  --sequential     Use only one thread, avoiding non-determinism in the output order.
  --requireisos    Require isogenies useful for a "simplified SWU" hash to curve.
  --sortpq         Sort p smaller than q.
  --ignoretwist    Ignore twist security.
  --nearpowerof2   Search for primes near a power of 2, rather than with low Hamming weight.
  <min-bitlength>  Both primes should have this minimum bit length.
  <min-2adicity>   Both primes should have this minimum 2-adicity.
  <stretch>        Find more candidates, by filtering from 2-adicity smaller by this many bits.
""")
        return

    L          = int(args[0])
    twoadicity = int(args[1]) if len(args) > 1 else DEFAULT_TWOADICITY
    stretch    = int(args[2]) if len(args) > 2 else DEFAULT_STRETCH

    print("Using %d processes." % (processes,))
    pool = Pool(processes=processes)

    try:
        for wid in range(processes):
            pool.apply_async(worker, (strategy, L, twoadicity, stretch, requireisos, sortpq, twistsec, wid, processes))

        while True:
            sleep(1000)
    except (KeyboardInterrupt, SystemExit):
        pass
    finally:
        pool.terminate()

def worker(*args):
    try:
        real_worker(*args)
    except (KeyboardInterrupt, SystemExit):
        pass
    except:
        print_exc()

def real_worker(*args):
    for (p, q, bp, bq, zetap, zetaq, qdesc, primp, primq, secp, secq, twsecp, twsecq,
         embeddivp, embeddivq, twembeddivp, twembeddivq, iso_Ep, iso_Eq) in find_nice_curves(*args):
        output  = "\n"
        output += "p   = %s\n" % format_weight(p)
        output += "q   = %s\n" % format_weight(q)
        output += "    = %s\n" % qdesc
        output += "ζ_p = %s (mod p)\n" % format_weight(int(zetap), detail=False)
        output += "ζ_q = %s (mod q)\n" % format_weight(int(zetaq), detail=False)

        output += "Ep/Fp : y^2 = x^3 + %d\n" % (bp,)
        output += "Eq/Fq : y^2 = x^3 + %d\n" % (bq,)

        output += "gcd(p-1, α) = 1 for α ∊ {%s}\n" % (", ".join(map(str, find_gcd_primes(p))),)
        output += "gcd(q-1, α) = 1 for α ∊ {%s}\n" % (", ".join(map(str, find_gcd_primes(q))),)

        output += "%d is %ssquare and %sprimitive in Fp\n" % (bp, "" if Mod(bp, p).is_square() else "non", "" if primp else "non")
        output += "%d is %ssquare and %sprimitive in Fq\n" % (bq, "" if Mod(bq, q).is_square() else "non", "" if primq else "non")

        output += "Ep security = %.1f, embedding degree = (q-1)/%d\n" % (secp, embeddivp)
        output += "Eq security = %.1f, embedding degree = (p-1)/%d\n" % (secq, embeddivq)

        output += "Ep twist security = %.1f, embedding degree = (2p + 1 - q)/%d\n" % (twsecp, twembeddivp)
        output += "Eq twist security = %.1f, embedding degree = (2q + 1 - p)/%d\n" % (twsecq, twembeddivq)

        if iso_Ep is not None:
            output += "iso_Ep = %r\n" % (iso_Ep,)
            output += "iso_Ep maps = %r\n" % (iso_Ep.rational_maps(),)
        elif ISOGENY_DEGREE_MAX > 0:
            output += "No Ep isogenies for simplified SWU with degree ≤ %d\n" % (ISOGENY_DEGREE_MAX,)

        if iso_Eq is not None:
            output += "iso_Eq = %r\n" % (iso_Eq,)
            output += "iso_Eq maps = %r\n" % (iso_Eq.rational_maps(),)
        elif ISOGENY_DEGREE_MAX > 0:
            output += "No Eq isogenies for simplified SWU with degree ≤ %d\n" % (ISOGENY_DEGREE_MAX,)

        print(output)  # one syscall to minimize tearing

main()
