#!/usr/bin/env sage

# Find the smallest element > 1 of { \omega^j : j \in [0, 2^32) }, over the Pasta Fp and Fq.
#
# This is a bit clunky at the moment since the threads work independently on subsets
# of the space, so it requires you to scan the output by eye to get the actual smallest
# element for each field.

import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc

if sys.version_info[0] == 2:
    range = xrange


p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

def check(ps):
    workers = cpu_count()//len(ps)
    pool = Pool(processes=workers*len(ps))

    try:
        for (which, p) in ps.items():
            print("Checking %s = %r" % (which, p))
            t = p >> 32
            omega = GF(p).multiplicative_generator()^t
            assert omega.multiplicative_order() == 1<<32

            for wid in range(1, workers+1):
                pool.apply_async(worker, (which, p, omega, wid, workers))

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

def real_worker(which, p, omega, wid, workers):
    print("Worker %d for %s" % (wid, which))

    lowest = 1<<240
    dot = workers*65536

    x = omega^wid
    m = omega^workers

    for i in range(wid, 1<<32, workers):
        if i % dot == 1:
            sys.stdout.write('.')
            sys.stdout.flush()

        if int(x) < lowest:
            lowest = int(x)
            print("\n%s: i = %r, %r (%d bits)" % (which, i, lowest, len(format(lowest, 'b'))))

        x *= m

    return lowest


check({"p": p, "q": q})
