Tweedledum/Tweedledee supporting evidence
-----------------------------------------

This repository contains supporting evidence that the amicable pair of
prime-order curves:

* Ep : y^2 = x^3 + 5 over GF(p) of order q, called Tweedledum;
* Eq : y^2 = x^3 + 5 over GF(q) of order p, called Tweedledee;

with

* p = 57896044618658097711785492504343953925989756877607163761872965584918954377217
* q = 57896044618658097711785492504343953925989756877607147991657089165100807356417

satisfy *some* of the [SafeCurves criteria](https://safecurves.cr.yp.to/index.html).

The criteria that are *not* satisfied are, in summary:

* large CM discriminant (both curves have CM discriminant 3, as a consequence of how
  they were constructed);
* completeness (complete formulae are possible, but not according to the Safe curves
  criterion);
* ladder support (not possible for prime-order curves);
* Elligator 2 support (indistinguishability is possible using
  [Elligator Squared](https://ifca.ai/pub/fc14/paper_25.pdf), but not using Elligator 2);
* twist security;
* rigidity.

Prerequisites:

* apt-get install sagemath
* pip install sortedcontainers

Run ``sage verify.sage Ep`` and ``sage verify.sage Eq``; or ``./run.sh`` to run both
and also print out the results.
