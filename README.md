Tweedledum/Tweedledee supporting evidence
-----------------------------------------

This repository contains supporting evidence that the amicable pair of
prime-order curves:

* Ep : y^2 = x^3 + 5 over GF(p) of order q, called Tweedledum;
* Eq : y^2 = x^3 + 5 over GF(q) of order p, called Tweedledee;

with

* p = 2^254 + 4707489545178046908921067385359695873
* q = 2^254 + 4707489544292117082687961190295928833

satisfy *some* of the [SafeCurves criteria](https://safecurves.cr.yp.to/index.html).

The criteria that are *not* satisfied are, in summary:

* large-magnitude CM discriminant (both curves have CM discriminant of absolute value 3,
  as a consequence of how they were constructed);
* completeness (complete formulae are possible, but not according to the Safe curves
  criterion);
* ladder support (not possible for prime-order curves);
* Elligator 2 support (indistinguishability is possible using
  [Elligator Squared](https://ifca.ai/pub/fc14/paper_25.pdf), but not using Elligator 2).

Tweedledum/Tweedledee is the first cycle output by
``sage amicable.sage --sequential --nearpowerof2 255 32``.

(The `--sequential` option makes the output completely deterministic and so resolves
ambiguity about which result is "first". For exploratory searches it is faster not to
use `--sequential`.)

**The cycle we call Tweedledum/Tweedledee has changed from the initial (September 2019) draft of the paper.**

Prerequisites:

* apt-get install sagemath
* pip install sortedcontainers

Run ``sage verify.sage Ep`` and ``sage verify.sage Eq``; or ``./run.sh`` to run both
and also print out the results.
