Tweedledum/Tweedledee supporting evidence
-----------------------------------------

This repository contains supporting evidence that the amicable pair of
prime-order curves:

* Ep : y^2 = x^3 + 5 over GF(p) of order q, called (provisional) Tweedledum;
* Eq : y^2 = x^3 + 5 over GF(q) of order p, called (provisional) Tweedledee;

with

* p = 2^254 + 11429413694214642624661040171709366273
* q = 2^254 + 11429413694209135470422256387130130433

satisfy *some* of the [SafeCurves criteria](https://safecurves.cr.yp.to/index.html).

The criteria that are *not* satisfied are, in summary:

* large CM discriminant (both curves have CM discriminant 3, as a consequence of how
  they were constructed);
* completeness (complete formulae are possible, but not according to the Safe curves
  criterion);
* ladder support (not possible for prime-order curves);
* Elligator 2 support (indistinguishability is possible using
  [Elligator Squared](https://ifca.ai/pub/fc14/paper_25.pdf), but not using Elligator 2);
* twist security.

(Provisional) Tweedledum/Tweedledee is the first cycle output by
``sage amicable.sage --nearpowerof2 255 30``.

**Which cycle we call Tweedledum/Tweedledee is subject to change as we make further
optimizations to Halo.**

Prerequisites:

* apt-get install sagemath
* pip install sortedcontainers

Run ``sage verify.sage Ep`` and ``sage verify.sage Eq``; or ``./run.sh`` to run both
and also print out the results.
