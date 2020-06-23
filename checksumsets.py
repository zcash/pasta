#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Dependencies:
#   <https://pypi.org/project/bintrees/> (pip install bintrees)
#   <https://pypi.org/project/Pillow/> (pip install Pillow), if --animate is used

import sys
from dataclasses import dataclass
from typing import Optional
from collections import deque
from math import log


# From the Halo paper:

# Let A = [0, 2^{λ/2 + 1} + 2^{λ/2} - 1]. It is straightforward to verify that a, b ∈ A
# at the end of Algorithm 3 for any input \mathbf{r}.
# {In fact a, b ∈ [2^{λ/2} + 1, 2^{λ/2 + 1} + 2^{λ/2} - 1], but it is convenient to
# define A to start at 0.}
#
# Next we need to show that the mapping (a ⦂ A, b ⦂ A) ↦ (a ζ_q + b) mod q is injective.
# This will depend on the specific values of ζ_q and q, and can be cast as a sumset problem.
#
# We use the notation v·A + A for { (av + b) mod q : a, b ∈ A }, and A - A for -1·A + A.
# {We take sumsets using this notation to implicitly be subsets of F_q.}
# The question is then whether |ζ_q·A + A| = |A|^2.
#
# For intuition, note that if av + b = a'v + b' (mod q), with a ≠ a', we would have
# v = (b' - b)/(a - a') (mod q). Thus the number of v ∈ F_q for which |v·A + A| < |A|^2
# is at most (|A - A| - 1)^2. {We thank Robert Israel for this observation. [HI2019]}
# Since in our case (|A - A| - 1)^2 ≈ 9·2^130 is small compared to q ≈ 2^254, we would
# heuristically expect that |ζ_q·A + A| = |A|^2 unless there is some reason why ζ_q does
# not "behave like a random element of F_q".
#
# Of course ζ_q is *not* a random element of F_q, and so the above argument can only be
# used for intuition. Even when (|A - A| - 1)^2 is small compared to q, there are clearly
# values of ζ_q and q for which it would not hold. To prove that it holds in the needed
# cases for the Tweedledum and Tweedledee curves used in our implementation, we take a
# different tack.
#
# Define a distance metric δ_q on F_q so that δ_q(x, y) is the minimum distance between
# x and y around the ring of integers modulo q in either direction, i.e.
#
#    δ_q(x, y) = min(z, q - z) where z = (x - y) mod q
#
# Now let D_{q,ζ_q}(m) be the minimum δ_q-distance between any two elements of ζ_q·[0, m],
# i.e.
#
#    D_{q,ζ_q}(m) = min{ δ_q(a ζ_q, a' ζ_q ) : a, a' ∈ [0, m] }
#
# An algorithm to compute D_{q,ζ_q}(m) is implemented by checksumsets.py in [Hopw2019]
# [i.e. this file]; it works by iteratively finding each m at which D_{q,ζ_q}(m)
# decreases. [...]
#
# Now if D_{q,ζ_q}(2^{λ/2 + 1} + 2^{λ/2} - 1) ≥ 2^{λ/2 + 1} + 2^{λ/2}, then copies of
# A will "fit within the gaps" in ζ_q·A. That is, ζ_q·A + A will have |A|^2 elements,
# because all of the sets { ζ_q·{a} + A : a ∈ A } will be disjoint.
#
# The algorithm is based on the observation that the problem of deciding when
# D_{q,ζ_q}(m) next decreases is self-similar to deciding when it first decreases.
# It computes the exact min-distance at each decrease (not just a lower bound),
# which facilitates detecting any bugs in the algorithm. Also, we check correctness
# of the partial results up to a given bound on m, against a naive algorithm.

BRUTEFORCE_THRESHOLD = 100000

DEBUG = False

@dataclass
class State:
    u: Optional[int]
    m: Optional[int]
    n: Optional[int]
    d: Optional[int]


def D(q, zeta, mm, animator=None):
    if DEBUG: print("(q, zeta, mm) =", (q, zeta, mm))
    Dcheck = [] if BRUTEFORCE_THRESHOLD == 0 else bruteforce_D(q, zeta, min(mm, BRUTEFORCE_THRESHOLD))

    # (u + am) : a ∈ Nat is the current arithmetic progression
    # n is the previous min-distance
    # d is the current min-distance
    cur = State(u=0, m=1, n=q, d=zeta)
    old = None

    while True:
        # Consider values of x where D_{q,ζ_q}(x) decreases, i.e. where
        # D_{q,ζ_q}(x) < D_{q,ζ_q})(x-1).
        #
        # We keep track of an arithmetic progression (u + am) such that the next
        # value at which D_{q,ζ_q}(x) decreases will be for x in this progression,
        # at the point at which xζ gets close to (but not equal to) 0.
        #
        # TODO: explain why the target is always 0.
        #
        # If we set s = floor(n/d), then D_{q,ζ_q}(x) can decrease at a = s
        # and potentially also at a = s+1.
        assert (cur.m*zeta) % q in (cur.d, q - cur.d)
        s = cur.n // cur.d
        x0 = cur.u + s*cur.m
        d0 = cur.n % cur.d
        if DEBUG: print("\n(x0, d0, cur, s) =", (x0, d0, cur, s))
        assert dist(0, x0*zeta, q) == d0
        if x0-1 < len(Dcheck): assert Dcheck[x0-1] == cur.d
        if x0 > mm:
            if animator is not None:
                animator.render(q, zeta, old, cur, None, s)
            return cur.d
        if x0 < len(Dcheck): assert Dcheck[x0] == d0

        x1 = cur.u + (s+1)*cur.m
        d1 = (s+1)*cur.d - cur.n
        if d1 < d0:
            if DEBUG: print("(x1, d1, cur, s+1) =", (x1, d1, cur, s+1))
            assert dist(0, x1*zeta, q) == d1
            if x1-1 < len(Dcheck): assert Dcheck[x1-1] == d0
            if x1 > mm:
                if animator is not None:
                    animator.render(q, zeta, old, cur, None, s+1)
                return d0
            if x1 < len(Dcheck): assert Dcheck[x1] == d1

            # This is the case where the smaller new distance is past zero.
            # The next iteration should consider the region of size d0 starting at x = x0
            # (i.e. just before we went past zero) and increasing by x1, i.e. dividing
            # that region by intervals of d1.
            new = State(u=x0, m=x1, n=d0, d=d1)
        else:
            # This is the case where the smaller new distance is short of zero.
            # The next iteration should check the region of size cur.d - d0 starting at x = x1
            # (i.e. the wraparound past zero) and increasing by x0, i.e. dividing that
            # region by intervals of d0.
            new = State(u=x1, m=x0, n=cur.d - d0, d=d0)

        assert dist(0, new.u*zeta, q) in (new.n, q - new.n)
        #if dist(0, new.u*zeta, q) != new.n: print("hmm")

        if animator is not None:
           animator.render(q, zeta, old, cur, new, s)

        (old, cur) = (cur, new)


def bruteforce_D(q, zeta, mm):
    # Can't use sortedcontainers because its data structures are backed by
    # lists-of-lists, not trees. We must have O(log n) insert, prev, and succ.
    from bintrees import RBTree as sortedset

    resD = deque([zeta])
    lastd = zeta
    S = sortedset()
    S.insert(0, None)
    S.insert(q, None)
    for x in range(1, mm+1):
        v = (x*zeta) % q
        S.insert(v, None)
        vp = S.prev_key(v)
        vs = S.succ_key(v)
        d = min(v-vp, vs-v)
        resD.append(d)
        #if DEBUG and d < lastd: print((x, d))
        lastd = d

    return list(resD)

def dist(x, y, q):
    z = (x-y+q) % q
    return min(z, q-z)

def signed_mod(x, q):
    r = x % q
    return r if r <= q//2 else r-q


class Animator:
    fontfile = '/usr/share/texlive/texmf-dist/fonts/truetype/google/roboto/Roboto-Regular.ttf'

    frame_duration = 20 # ms
    pause_frames   = 35
    zoom_frames    = 45

    width          = 800 # pixels
    height         = 400 # pixels
    oversample     = 3
    line_halfwidth = 1 # subpixels

    ground_colour  = '#ffffff' # white
    scale_colour   = '#0000a0' # blue
    midline_colour = '#c00000' # red
    old_colour     = '#a0a0a0' # grey
    cur_colour     = '#000000' # black
    new_colour     = '#008000' # green
    final_colour   = '#c00000' # red

    def __init__(self, name):
        # We don't want to depend on PIL unless an Animator is instantiated.
        from PIL import Image, ImageDraw, ImageColor, ImageFont
        self.Image = Image
        self.ImageDraw = ImageDraw
        self.ImageColor = ImageColor

        self.font = ImageFont.truetype(self.fontfile, 20*self.oversample, index=0, encoding='unic')
        self.font_super = ImageFont.truetype(self.fontfile, 12*self.oversample, index=0, encoding='unic')
        self.images = deque()
        self.name = name

    def render(self, q, zeta, old, cur, new, s):
        n = min(cur.n, q//2)
        for aa in range(1, s+1):
            self.render_zoomed(q, zeta, old, cur, None, aa,  n)

        if new is None:
            self.render_zoomed(q, zeta, old, cur, new,  s+1, n, final=True)
            return

        self.render_zoomed(    q, zeta, old, cur, new,  s+1, n, frames=self.pause_frames)

        step = (1.0*n/new.n - 1.0)/self.zoom_frames
        for zoom in range(1, self.zoom_frames):
            n_scale = int(n/(1.0 + zoom*step))
            self.render_zoomed(q, zeta, old, cur, new,  s+1, n_scale)

        self.render_zoomed(    q, zeta, old, cur, new,  s+1, new.n, frames=self.pause_frames)

    def render_zoomed(self, q, zeta, old, cur, new, aa, n_scale, frames=1, final=False):
        px = self.oversample
        lx = self.line_halfwidth
        (w, h) = (self.width * px, self.height * px)
        scale = (w/2)/n_scale
        xmid = w//2
        ymid = (40*px + h)//2

        image = self.Image.new('RGB', (w, h), color=self.ground_colour)
        image.convert('P')
        draw = self.ImageDraw.Draw(image)

        bits = int(log(n_scale, 2))
        for tick in range(bits-3, bits+1):
            xoff = int(scale*(1<<tick))
            draw.text((xmid-xoff-21*px, 7*px), '−2', self.scale_colour, font=self.font)
            draw.text((xmid-xoff, px), str(tick), self.scale_colour, font=self.font_super)
            draw.rectangle((xmid-xoff-lx, 30*px, xmid-xoff+lx, 40*px), fill=self.scale_colour)

            draw.text((xmid+xoff-21*px, 7*px), '+2', self.scale_colour, font=self.font)
            draw.text((xmid+xoff, px), str(tick), self.scale_colour, font=self.font_super)
            draw.rectangle((xmid+xoff-lx, 30*px, xmid+xoff+lx, 40*px), fill=self.scale_colour)

        draw.rectangle((xmid-lx, 0, xmid+lx, h), fill=self.midline_colour)
        draw.rectangle((0, 40*px-lx, w, 40*px+lx), fill=self.scale_colour)

        if old is not None:
            old_aa = (old.n // old.d)+1
            for a in range(old_aa+1):
                x = signed_mod(zeta*(old.u + a*old.m), q)
                xpos = w//2 + int(scale*x)
                draw.rectangle((xpos-lx, 40*px, xpos+lx, h), fill=self.old_colour)

        for a in range(aa+1):
            x = signed_mod(zeta*(cur.u + a*cur.m), q)
            xpos = w//2 + int(scale*x)
            draw.rectangle((xpos-lx, 40*px, xpos+lx, h), fill=self.cur_colour)

        if new is not None:
            x = signed_mod(zeta*new.u, q)
            xpos = w//2 + int(scale*x)
            draw.rectangle((xpos, ymid-lx, xmid, ymid+lx), fill=self.new_colour)
            draw.rectangle((xpos-lx, 40*px, xpos+lx, h), fill=self.new_colour)

        if final:
            x = signed_mod(zeta*(cur.u + aa*cur.m), q)
            xpos = w//2 + int(scale*x)
            draw.rectangle((xpos, ymid-lx, xmid, ymid+lx), fill=self.final_colour)
            draw.rectangle((xpos-lx, 40*px, xpos+lx, h), fill=self.final_colour)

        image = image.resize((self.width, self.height), self.Image.ANTIALIAS)
        for f in range(frames):
            self.images.append(image)
        sys.stdout.write('.')
        sys.stdout.flush()

    def save(self):
        filename = 'animation-%s.gif' % (self.name,)
        print("Saving %s..." % (filename,))

        first, *rest = list(self.images)
        # Save as animated GIF. We can convert to a more space-efficient format separately.
        first.save(fp=filename, format='GIF', append_images=rest,
                   save_all=True, duration=self.frame_duration, loop=1)
        del self.images


def check_sumset(name, q, zeta, limit, animator=None):
    print("===== %s =====" % (name,))

    Dq = D(q, zeta, limit-1, animator)
    print("\nD_%s = %s" % (name, Dq))
    print("    %s" % ('≥' if Dq >= limit else '<'), limit)

    if animator is not None:
        animator.save()

    assert Dq >= limit


def main():
    args = sys.argv[1:]
    if "--help" in args:
        print("Usage: checksumsets.py [--animate]")
        return

    halflambda = 64
    limit = 3<<halflambda

    # Tweedledum and Tweedledee
    p = (1<<254) + 4707489545178046908921067385359695873
    q = (1<<254) + 4707489544292117082687961190295928833
    zeta_p = 9513155655832138286304767221959569637168364952810827555227185832555034233288
    zeta_q = 24775483399512474214391554062650059912556682109176536098332128018848638018813

    # Tests
    global DEBUG
    DEBUG = False
    assert(D(65537, 6123, 10000, None) == 3)
    assert(D(1299721, 538936, 10000, None) == 41)
    assert(D(179424691, 134938504, 100000, None) == 121)

    p_params = ("p", p, zeta_p, limit)
    q_params = ("q", q, zeta_q, limit)

    DEBUG = False
    for (name, prime, zeta, limit) in (p_params, q_params):
        animator = None
        if "--animate" in args:
            animator = Animator(name)

        check_sumset(name, prime, zeta, limit, animator=animator)

main()
