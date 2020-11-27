Pallas/Vesta supporting evidence
--------------------------------

This repository contains supporting evidence that the amicable pair of
prime-order curves:

* Ep : y^2 = x^3 + 5 over GF(p) of order q, called Pallas;
* Eq : y^2 = x^3 + 5 over GF(q) of order p, called Vesta;

with

* p = 2^254 + 45560315531419706090280762371685220353
* q = 2^254 + 45560315531506369815346746415080538113

satisfy *some* of the [SafeCurves criteria](https://safecurves.cr.yp.to/index.html).

The criteria that are *not* satisfied are, in summary:

* large-magnitude CM discriminant (both curves have CM discriminant of absolute value 3,
  as a consequence of how they were constructed);
* completeness (complete formulae are possible, but not according to the Safe curves
  criterion);
* ladder support (not possible for prime-order curves);
* Elligator 2 support (indistinguishability is possible using
  [Elligator Squared](https://ifca.ai/pub/fc14/paper_25.pdf), but not using Elligator 2);
* twist security above 100 bits for Pallas.

Pallas/Vesta is the first cycle output by
``sage amicable.sage --sequential --requireisos --sortpq --ignoretwist --nearpowerof2 255 32``.

(The `--sequential` option makes the output completely deterministic and so resolves
ambiguity about which result is "first". For exploratory searches it is faster not to
use `--sequential`.)

Prerequisites:

* ``apt-get install sagemath``

Run ``sage verify.sage Ep`` and ``sage verify.sage Eq``; or ``./run.sh`` to run both
and also print out the results.

The output of ``amicable.sage`` with the above options includes isogenies of degree 3,
suitable for use with the "simplified SWU" method for hashing to an elliptic curve.
This is based on code from Appendix A of [Wahby and Boneh 2019](https://eprint.iacr.org/2019/403.pdf).

To check the correctness of the endomorphism optimization described in the Halo paper, run
``python3 injectivitylemma.py`` and ``python3 checksumsets.py``. To also generate animations
showing the minimum distances between multiples of ζ used in the proof, run ``./animation.sh``.

``animation.sh`` has the following prerequisites:

* ``apt-get install ffmpeg ffcvt``
* ``pip3 install bintrees Pillow``

``checksumsets.py`` on its own only requires the ``bintrees`` Python package.
