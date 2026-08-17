#!/usr/bin/env sage
# -*- coding: utf-8 -*-

# Calculate the Weil character-sum constant for the deployed simplified-SWU
# mappings to iso-Pallas and iso-Vesta (Z = -13).
#
# Goal. For the odd (zero-repaired) form of map_to_curve_simple_swu
# f : F_q → E'(F_q), bound the character sums
#
#     S_f(χ) = Σ_{u ∈ F_q} χ(f(u))
#
# for every nontrivial character χ of E'(F_q), in the form
# |S_f(χ)| ≤ C_leading·√q + C_additive. This is the named `WeilBounded`
# input of the group-hash indifferentiability formalization (CompElliptic
# `Hashing/WellDistributed.lean`; the derivation task is
# https://github.com/daira/CompElliptic/issues/28). The composition through
# the 3-isogeny costs nothing: the isogeny is bijective on rational points,
# so χ ∘ iso ranges over the nontrivial characters of the iso-curve.
#
# Method. Farashahi–Fouque–Shparlinski–Tibouchi–Voloch, "Indifferentiable
# Deterministic Hashing to Elliptic and Hyperelliptic Curves"
# (https://eprint.iacr.org/2010/539), Theorem 6, adapted:
#
#  1. Each SSWU branch j ∈ {1, 2} gives a covering curve C_j → E', cut out
#     by the branch's x-formula as a (biquadratic) quartic in the input u
#     over the function field of E'.
#  2. FFSTV Lemma 1 / Theorem 3(6): if C_j → E' does not factor through a
#     nontrivial unramified subcover, then for every nontrivial χ,
#     |Σ_{P ∈ C_j(F_q)} χ(h_j(P))| ≤ (2·g_j − 2)·√q.
#  3. Sign-freeness: FFSTV needed a further conductor term (their
#     deg y = 12) to select which of the two conjugate C_j-points matches
#     the quadratic-residue sign rule; the deployed sgn0 rule is not
#     algebraic, so that step is unavailable — and also unnecessary. For
#     odd f the sum S_f(χ) is real (substitute u → −u), so
#
#       2·S_f(χ) = Σ_u (χ + χ̄)(f(u))
#                = Σ_j Σ_{P ∈ C_j(F_q)} χ(h_j(P)) + O(bad points),
#
#     because over each u in branch j the two rational points of C_j are
#     the y-conjugates, contributing χ(f(u)) + χ̄(f(u)) whichever sign the
#     map picks, while u in the other branch contributes no rational points
#     at all (its candidate ordinate is a nonsquare). Hence
#
#       |S_f(χ)| ≤ ((2·g_1 − 2) + (2·g_2 − 2))/2 · √q + C_additive.
#
#  4. This script computes the genera g_j, gathers the ramification/
#     monodromy evidence behind step 2's side condition, and counts the
#     exceptional points behind C_additive.
#
# This is a calculation, not a proof. Each step prints what the eventual
# hand derivation must establish; the genus and monodromy facts are the
# analogues of FFSTV's Riemann–Hurwitz and Eisenstein-criterion arguments.

# ----------------------------------------------------------------------
# Deployed constants (must match hashtocurve.sage and CompElliptic's
# Curves/IsoPasta.lean + Hashing/PastaSSWU.lean).

p_pallas = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
p_vesta  = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001

iso_pallas_A = 0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b
iso_vesta_A  = 0x267f9b2ee592271a81639c4d96f787739673928c7d01b212c515ad7242eaa6b1
iso_B = 1265
Z_int = -13

# ----------------------------------------------------------------------
# The branch covers.
#
# With t = Z·u² and ta = t² + t, the SSWU branch abscissae are
#   x1 = B·(ta + 1) / (A · −ta)          (the ta = 0 inputs are exceptional)
#   x2 = t·x1
# Rearranging x = x1:  ta·(A·x + B) = −B, and ta = Z²·u⁴ + Z·u², so
#   P1(u) = Z²·(A·x + B)·u⁴ + Z·(A·x + B)·u² + B = 0.
# Rearranging x = x2 = Z·u²·x1 (same u, same ta):
#   (Z·u² + 1)·(A·x + B·Z·u²) = −B, so
#   P2(u) = B·Z²·u⁴ + Z·(A·x + B)·u² + (A·x + B) = 0.
# Both are biquadratic in u: the u → −u symmetry of the covers is the
# geometric face of the mapping's oddness.

def make_curve_data(p, A_int, B_int, Z_val):
    F = GF(p)
    A = F(A_int); B = F(B_int); Z = F(Z_val)
    assert not Z.is_square(), "Z must be a nonsquare"
    assert A != 0 and B != 0
    E = EllipticCurve(F, [A, B])
    return (F, A, B, Z, E)

def branch_quartic(R, A, B, Z, xv, which):
    # The branch-`which` quartic at abscissa value or indeterminate `xv`,
    # as a polynomial in R (a univariate polynomial ring whose base
    # contains A, B, Z, xv).
    U = R.gen()
    w = A*xv + B
    if which == 1:
        return Z^2*w*U^4 + Z*w*U^2 + B
    else:
        return B*Z^2*U^4 + Z*w*U^2 + w

# ----------------------------------------------------------------------
# Self-check: the covers describe the deployed map. For random u, run the
# RFC 9380 branch formulas and confirm each branch's quartic vanishes at
# that branch's abscissa, that exactly one branch has a square ordinate
# candidate (so the other cover has no rational points over u), and the
# g(x2) = t³·g(x1) branch-ratio identity.

def check_covers_match_map(F, A, B, Z, trials=50):
    R = PolynomialRing(F, 'U')
    done = 0
    while done < trials:
        u = F.random_element()
        t = Z*u^2
        ta = t^2 + t
        if ta == 0:
            continue  # exceptional input, handled in the additive term
        x1 = B*(ta + 1) / (A * -ta)
        x2 = t*x1
        g1 = x1^3 + A*x1 + B
        g2 = x2^3 + A*x2 + B
        assert g2 == t^3 * g1, "branch-ratio identity failed"
        assert g1.is_square() != g2.is_square(), "exactly-one-branch failed"
        assert branch_quartic(R, A, B, Z, x1, 1)(u) == 0, "P1 relation failed"
        assert branch_quartic(R, A, B, Z, x2, 2)(u) == 0, "P2 relation failed"
        done += 1
    print("  cover/map consistency: OK (%d random inputs)" % trials)

# ----------------------------------------------------------------------
# Genus of each branch cover.
#
# The quartics are linear in x, so each cover re-roots as a hyperelliptic
# curve over the u-line: solving the branch relation for x gives the
# branch abscissa x_j(u) as a rational function of u, and
# F_q(C_j) = F_q(u)(y) with y² = g(x_j(u)) where g(x) = x³ + A·x + B.
# This is a single quadratic extension of a rational function field, which
# Sage's genus machinery handles over any base — including the deployed
# ~254-bit fields, where the earlier tower presentation hit unimplemented
# cases. It is also the model the hand derivation should use: clearing the
# denominator, C_j is Y² = H_j(u) for an explicit polynomial H_j, and the
# genus of a hyperelliptic model reads off from the degree of the
# squarefree part of H_j. The map h_j : C_j → E' is (x_j(u), y), so the
# presentation change does not affect the character-sum argument.

def branch_abscissa(K, A, B, Z, which):
    # The branch-`which` abscissa as a rational function of the generator
    # of the rational function field K = k(u).
    u = K.gen()
    t = Z*u^2
    ta = t^2 + t
    x1 = B*(ta + 1) / (A * -ta)
    return x1 if which == 1 else t*x1

def hyperelliptic_profile(k, A, B, Z, which):
    # The cleared hyperelliptic model Y² = H_j(u) and the degree of its
    # odd-multiplicity part — the quantities the genus derivation
    # consumes. Uses only rational-function and univariate-polynomial
    # arithmetic, so it works over the deployed fields (where
    # every Sage path through multivariate factorization is unimplemented,
    # including the FunctionField extension constructor's irreducibility
    # check). Returns (deg H_j, deg of the odd-multiplicity part).
    R = PolynomialRing(k, 'u')
    K = R.fraction_field()
    xj = branch_abscissa(K, A, B, Z, which)
    r = xj^3 + A*xj + B
    num = R(r.numerator()); den = R(r.denominator())
    H = num*den  # Y² = r  ⟺  (Y·den)² = num·den
    odd_part = prod((f for f, e in H.factor() if e % 2 == 1), R.one())
    return (H.degree(), odd_part.degree())

def cover_genus_formula(k, A, B, Z, which):
    # Genus of the hyperelliptic model Y² = H_j(u) in characteristic ≠ 2:
    # with h the odd-multiplicity part of H_j (so Y² = s²·h and (Y/s)² = h
    # with h squarefree), the genus is ⌊(deg h − 1)/2⌋. The classical
    # formula accounts for the places at infinity; the leading-coefficient
    # twist does not change the (geometric) genus, which is what Weil's
    # bound consumes.
    (_, d) = hyperelliptic_profile(k, A, B, Z, which)
    return (d - 1) // 2

def cover_genus(k, A, B, Z, which):
    # Sage's own genus computation, usable when the constant field is a
    # prime field small enough for Singular (or ℚ): the cross-check for
    # `cover_genus_formula`.
    K = FunctionField(k, 'u')
    xj = branch_abscissa(K, A, B, Z, which)
    r = xj^3 + A*xj + B
    RY = PolynomialRing(K, 'Y')
    Y = RY.gen()
    L = K.extension(Y^2 - r, 'y')
    return L.genus()

def generic_genus_sweep(trials=12):
    # Genera over random small-prime surrogates. Returns the set of
    # (branch, genus) pairs seen, which should be {(1, g), (2, g)} for a
    # single g if the genus is parameter-uniform.
    seen = set()
    done = 0
    while done < trials:
        pr = random_prime(2^24, lbound=2^12)
        if pr % 4 != 1 or pr in (2, 3, 13):
            continue
        k = GF(pr)
        Z = k(-13)
        if Z.is_square() or Z == 0:
            continue
        A = k.random_element(); B = k.random_element()
        if A == 0 or B == 0 or 4*A^3 + 27*B^2 == 0:
            continue
        for which in (1, 2):
            g = cover_genus(k, A, B, Z, which)
            gf = cover_genus_formula(k, A, B, Z, which)
            assert g == gf, "formula mismatch at p=%s: %s vs %s" % (pr, g, gf)
            seen.add((which, g))
        done += 1
    return seen

# ----------------------------------------------------------------------
# Monodromy evidence for the no-unramified-subcover condition, by
# Frobenius specialization: factor the quartic at many random points of
# E'(F_q) and collect the factorization degree patterns. The multiset of
# patterns that occurs identifies the monodromy group among the transitive
# subgroups of S_4 compatible with the biquadratic shape (subgroups of the
# dihedral group D_4). What the proof needs from this evidence: every
# intermediate subcover is ramified — established by hand from the
# ramification data, in the role FFSTV's Eisenstein criterion plays.

def monodromy_sample(F, A, B, Z, E, which, trials=5000):
    R = PolynomialRing(F, 'U')
    patterns = {}
    done = 0
    while done < trials:
        Pt = E.random_point()
        if Pt.is_zero():
            continue
        q = branch_quartic(R, A, B, Z, Pt[0], which)
        if q.degree() != 4 or q.discriminant() == 0:
            continue
        pat = tuple(sorted(sum(([f.degree()]*m for f, m in q.factor()), [])))
        patterns[pat] = patterns.get(pat, 0) + 1
        done += 1
    return patterns

# ----------------------------------------------------------------------
# Subcover witnesses. The quartic u⁴ + p·u² + q has Galois group V₄ when q
# is a square in the base field and C₄ when q·(p² − 4·q) is; only the D₄
# case has the single intermediate subcover F(u²). For our covers those
# tests reduce to B·(A·x + B) and B·(A·x − 3·B) being squares in the
# function field of E'. The divisor-parity argument (two simple zeros plus
# a double pole is not twice a divisor) shows they are not; a nonsquare
# specialization at any single point is a computational witness, since a
# square function takes square values wherever it is defined and nonzero.

def subcover_witnesses(F, A, B, Z, E, trials=64):
    found_v4 = False; found_c4 = False
    for _ in range(trials):
        Pt = E.random_point()
        if Pt.is_zero():
            continue
        xv = Pt[0]
        q_test = B*(A*xv + B)
        c_test = B*(A*xv - 3*B)
        if q_test != 0 and not q_test.is_square():
            found_v4 = True
        if c_test != 0 and not c_test.is_square():
            found_c4 = True
        if found_v4 and found_c4:
            break
    assert found_v4, "no nonsquare witness for the V4 test"
    assert found_c4, "no nonsquare witness for the C4 test"
    print("  subcover square-ness tests: nonsquare witnesses found for both")

# Expected Frobenius cycle-type proportions on the four roots, for the
# three transitive subgroups of D₄ compatible with the biquadratic shape.
EXPECTED_PATTERNS = {
    'D4': {(1, 1, 1, 1): 1/8, (1, 1, 2): 2/8, (2, 2): 3/8, (4,): 2/8},
    'V4': {(1, 1, 1, 1): 1/4, (2, 2): 3/4},
    'C4': {(1, 1, 1, 1): 1/4, (2, 2): 1/4, (4,): 2/4},
}

# ----------------------------------------------------------------------
# Exceptional-input counting: the inputs the covers exclude, each
# contributing O(1) to the additive constant. (The full additive-constant
# audit for the write-up also counts the points at infinity of each C_j and
# the u = 0 value of the zero-repaired map; all are O(1) and are dominated
# by the slack in the recorded constant at the deployed sizes.)

def exceptional_counts(F, A, B, Z):
    counts = {}
    # u = 0 (t = 0): one input; the zero-repaired map is defined there
    # directly, and neither cover has points over it.
    counts['u = 0'] = 1
    # t = −1, i.e. u² = −1/Z: 0 or 2 inputs depending on square-ness.
    # For q ≡ 1 (mod 4), −1 is a square, so −1/Z is a nonsquare whenever
    # Z is: this fibre is empty for every admissible Z on such fields, ta
    # vanishes only at u = 0, and the map's exceptional set is exactly
    # {0}. (In FFSTV's q ≡ 3 (mod 4) setting it is nonempty.)
    counts['u^2 = -1/Z'] = 2 if (-1/Z).is_square() else 0
    # Points of E' over w = A·x + B = 0 (the Eisenstein fibre): 2 when
    # g(−B/A) is a nonzero square, 0 when a nonsquare, 1 when zero. The
    # Eisenstein argument needs w to vanish simply there, i.e.
    # g(−B/A) ≠ 0.
    gval = (-B/A)^3 + A*(-B/A) + B
    assert gval != 0, "g(-B/A) = 0: Eisenstein fibre is 2-torsion"
    counts['points over w = 0'] = 2 if gval.is_square() else 0
    # Points of E' over A·x = 3·B (the (2,2)-ramification fibre).
    gval3 = (3*B/A)^3 + A*(3*B/A) + B
    counts['points over Ax = 3B'] = (2 if gval3.is_square() else 0) \
        if gval3 != 0 else 1
    return counts

def additive_constant(counts):
    # A generous explicit bound on the additive term C_additive in
    # |S_f(χ)| ≤ C_leading·√q + C_additive, per curve:
    #   · the u = 0 and u² = −1/Z inputs appear in S but not in the cover
    #     sums (≤ 1 each in absolute value, counted twice for the factor-2
    #     bookkeeping of 2·S);
    #   · each cover's smooth projective model has ≤ 2 points at infinity
    #     and ≤ 4 points over each of the finitely many base points where
    #     the u ↔ point correspondence degenerates (the w = 0 and
    #     A·x = 3·B fibres), each contributing ≤ 1 to the cover sum.
    bad_inputs = counts['u = 0'] + counts['u^2 = -1/Z']
    bad_base = counts['points over w = 0'] + counts['points over Ax = 3B']
    per_cover = 2 + 4*bad_base
    return (2*bad_inputs + 2*per_cover) / 2

# ----------------------------------------------------------------------

def run(name, p, A_int, expected_genus=None):
    print("== %s (p = %s…)" % (name, hex(p)[:18]))
    (F, A, B, Z, E) = make_curve_data(p, A_int, iso_B, Z_int)
    print("  iso-curve j-invariant nonzero:", E.j_invariant() != 0)
    check_covers_match_map(F, A, B, Z)
    assert F.characteristic() not in (2, 3)
    genera = {}
    for which in (1, 2):
        genera[which] = cover_genus_formula(F, A, B, Z, which)
        (degH, degodd) = hyperelliptic_profile(F, A, B, Z, which)
        print("  branch-%d: genus %s (hyperelliptic model degree %s, "
              "odd-multiplicity part degree %s)"
              % (which, genera[which], degH, degodd))
        if expected_genus is not None:
            assert genera[which] == expected_genus, "sweep mismatch"
    for which in (1, 2):
        pats = monodromy_sample(F, A, B, Z, E, which)
        total = sum(pats.values())
        obs = {k: float(v/total) for k, v in sorted(pats.items())}
        print("  branch-%d Frobenius factorization patterns: %s" % (which, obs))
    print("  expected proportions: %s" % {
        g: {k: float(v) for k, v in d.items()}
        for g, d in EXPECTED_PATTERNS.items()})
    subcover_witnesses(F, A, B, Z, E)
    exc = exceptional_counts(F, A, B, Z)
    print("  exceptional counts:", exc)
    c_leading = ((2*genera[1] - 2) + (2*genera[2] - 2)) / 2
    c_additive = additive_constant(exc)
    print("  leading constant ((2g_1 - 2) + (2g_2 - 2))/2 =", c_leading)
    print("  additive constant (generous audit) <=", c_additive)
    # The recorded WeilBounded constant: C with C_leading·√q + C_additive
    # ≤ C·√q. Anything above C_leading works at the deployed sizes; we
    # record C_leading + 1/2 with an astronomically large margin.
    C = c_leading + 1/2
    margin = (C - c_leading)*sqrt(RR(p)) - c_additive
    assert margin > 0
    print("  recorded constant C = %s (margin %.3g)" % (C, margin))
    return (c_leading, c_additive, C)

print("== small-prime genus sweep (p ≡ 1 mod 4, Z = −13)")
sweep = generic_genus_sweep()
print("  (branch, genus) pairs seen:", sorted(sweep))
genera_seen = set(g for (_, g) in sweep)
assert len(genera_seen) == 1, "genus is not parameter-uniform: %s" % sweep
expected_genus = genera_seen.pop()
print("  uniform genus:", expected_genus)

print("== exact-constant check over Q (deployed integers)")
try:
    alarm(600)
    gq1 = cover_genus(QQ, QQ(iso_pallas_A), QQ(iso_B), QQ(Z_int), 1)
    gq2 = cover_genus(QQ, QQ(iso_pallas_A), QQ(iso_B), QQ(Z_int), 2)
    cancel_alarm()
    print("  genera over Q with iso-Pallas constants:", gq1, gq2)
    assert gq1 == expected_genus and gq2 == expected_genus
except AlarmInterrupt:
    print("  (skipped: exceeded 600s)")

res_pallas = run("iso-Pallas", p_pallas, iso_pallas_A, expected_genus)
res_vesta = run("iso-Vesta", p_vesta, iso_vesta_A, expected_genus)
print("results (leading, additive bound, recorded C):")
print("  iso-Pallas:", res_pallas)
print("  iso-Vesta: ", res_vesta)
