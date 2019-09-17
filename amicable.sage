# -*- coding: utf-8 -*-
import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import combinations

# Let Ep/Fp : y^2 = x^3 + bp
# Let Eq/Fq : y^2 = x^3 + bq

# p and q should each be ~ L bits.

DEFAULT_TWOADICITY = 21
DEFAULT_STRETCH = 0

COEFFICIENT_RANGE = (5,)
#COEFFICIENT_RANGE = xrange(1, 10000)

ACCEPTABLE_PRIMES = (5,)
#ACCEPTABLE_PRIMES = Primes()

# <https://eprint.iacr.org/2011/465.pdf>
# It is well known that if g is neither a square nor a cube in Fp, then all
# possible group orders an elliptic curve E : y^2 = x^3 + b can have over Fp
# occur as the order of one of the 6 twists with b \in {1, g, g^2, g^3, g^4, g^5}.

# <https://math.stackexchange.com/questions/127251/when-is-a-not-a-cube-mod-p>:
#   If p = 2 (mod 3) then all elements are cubes.
#   If p = 1 (mod 3) then a is a cube iff a^((p-1)/3) = 1.

# <https://cryptojedi.org/papers/pfcpo.pdf> section 2:
# [...] the order of a curve satisfying the norm equation 3V^2 = 4p - t^2 has one
# of the six forms {p+1 +/- t, p+1 +/- (t +/- 3V)/2} [IEEE Std 1363-2000, section
# A.14.2.3, item 6].
#
# We choose 4p = 3V^2 + t^2, where (V-1)/2 and (t-1)/2 are both multiples of 2^twoadicity.
#
# Then 4p = (3(V-1)^2 + 6(V-1) + 3) + ((t-1)^2 + 2(t-1) + 1)
#         = 3(V-1)^2 + 6(V-1) + (t-1)^2 + 2(t-1) + 4
#       p = 3((V-1)/2)^2 + 3(V-1)/2 + ((t-1)/2)^2 + (t-1)/2 + 1
#
# So p-1 will be a multiple of 2^twoadicity, and so will (p+1-t)-1 = (p-1)-(t-1).
#
# We'd also like both p and q to be 1 (mod 3), so that we have efficient endomorphisms
# on both curves. We explicitly check p = 1 (mod 3), and then if t is chosen to be
# 1 (mod 3) then p+1-t will be 1 (mod 3) (but we must still check q since it
# is not necessarily that order).

def low_hamming_order(L, twoadicity, wid, processes):
    Vlen = (L-1)//2 + 1
    Vbase = 1 << Vlen
    tlen = (L-1)//4
    tbase = 1 << tlen
    trailing_zeros = twoadicity+1
    for w in xrange(wid, tlen-trailing_zeros, processes):
        for Vc in combinations(xrange(trailing_zeros, Vlen), w):
            V = Vbase + sum([1 << i for i in Vc]) + 1
            assert(((V-1)/2) % (1<<twoadicity) == 0)
            for tw in xrange(1, w+1):
                for tc in combinations(xrange(trailing_zeros, tlen), tw):
                    t = tbase + sum([1 << i for i in tc]) + 1
                    assert(((t-1)/2) % (1<<twoadicity) == 0)
                    if t % 6 != 1:
                        continue
                    p4 = 3*V^2 + t^2
                    assert(p4 % 4 == 0)
                    p = p4//4
                    assert(p % (1<<twoadicity) == 1)
                    if p % 6 == 1 and is_pseudoprime(p):
                        yield (p, t, V)

def near_powerof2_order(L, twoadicity, wid, processes):
    trailing_zeros = twoadicity+1
    Vbase = isqrt((1<<(L+1))//3) >> trailing_zeros
    for Voffset in symmetric_range(10000, base=wid, step=processes):
        V = ((Vbase + Voffset) << trailing_zeros) + 1
        assert(((V-1)/2) % (1 << twoadicity) == 0)
        tmp = (1<<(L+1)) - 3*V^2
        if tmp < 0: continue
        tbase = isqrt(tmp) >> trailing_zeros
        for toffset in symmetric_range(10000):
            t = ((tbase + toffset) << trailing_zeros) + 1
            assert(((t-1)/2) % (1<<twoadicity) == 0)
            if t % 6 != 1:
                continue
            p4 = 3*V^2 + t^2
            assert(p4 % 4 == 0)
            p = p4//4
            assert(p % (1<<twoadicity) == 1)
            if p > 1<<(L-1) and p % 6 == 1 and is_pseudoprime(p):
                yield (p, t, V)

def find_nonsquare_noncube(p):
    for g_int in xrange(2, 1000):
        g = Mod(g_int, p)
        if g^((p-1)//3) != 1 and g^((p-1)//2) != 1:
            return g
    return None

def symmetric_range(n, base=0, step=1):
    for i in xrange(base, n, step):
        yield -i
        yield i+1

def find_nice_curves(strategy, L, twoadicity, stretch, wid, processes):
    for (p, t, V) in strategy(L, max(0, twoadicity-stretch), wid, processes):
        sys.stdout.write('.')
        sys.stdout.flush()

        for (q, qdesc) in ((p + 1 - t,          "p + 1 - t"),
                           (p + 1 + (t-3*V)//2, "p + 1 + (t-3*V)/2")):
            if q > 1<<(L-1) and q % 6 == 1 and q % (1<<twoadicity) == 1 and is_prime(q) and is_prime(p):
                bp = find_coefficient(p, q)
                if bp == None: continue
                bq = find_coefficient(q, p)
                if bq == None: continue
                gp = find_nonsquare_noncube(p)
                if gp == None: continue
                gq = find_nonsquare_noncube(q)
                if gq == None: continue

                aq = gq**((q-1)//3)
                assert(aq**3 == Mod(1, q))
                ap = gp**((p-1)//3)
                assert(ap**3 == Mod(1, p))
                yield (p, q, bp, bq, ap, aq, qdesc)

def find_coefficient(p, q):
    for b in COEFFICIENT_RANGE:
        E = EllipticCurve(GF(p), [0, b])
        if E.count_points() == q:
            return b
    return None

def find_lowest_prime(p):
    for r in ACCEPTABLE_PRIMES:
        if gcd(p-1, r) == 1:
            return r

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
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1:
        print("Usage: sage amicable.sage [--sequential] [--nearpowerof2] <min-bitlength> [<min-2adicity> [<stretch]]\n")
        return

    L          = int(args[0])
    twoadicity = int(args[1]) if len(args) > 1 else DEFAULT_TWOADICITY
    stretch    = int(args[2]) if len(args) > 2 else DEFAULT_STRETCH

    print("Using %d processes." % (processes,))
    pool = Pool(processes=processes)

    try:
        for wid in xrange(processes):
            pool.apply_async(worker, (strategy, L, twoadicity, stretch, wid, processes))

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
    for (p, q, bp, bq, ap, aq, qdesc) in find_nice_curves(*args):
        output  = "\n"
        output += "p   = %s\n" % format_weight(p)
        output += "q   = %s\n" % format_weight(q)
        output += "    = %s\n" % qdesc
        output += "α_p = %s (mod p)\n" % format_weight(int(ap), detail=False)
        output += "α_q = %s (mod q)\n" % format_weight(int(aq), detail=False)

        output += "Ep/Fp : y^2 = x^3 + %d (%ssquare)\n" % (bp, "" if Mod(bp, p).is_square() else "non")
        output += "Eq/Fq : y^2 = x^3 + %d (%ssquare)\n" % (bq, "" if Mod(bq, q).is_square() else "non")

        output += "gcd(p-1, %d) = 1\n" % find_lowest_prime(p)
        output += "gcd(q-1, %d) = 1\n" % find_lowest_prime(q)

        print(output)  # one syscall to minimize tearing

main()
