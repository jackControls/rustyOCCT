#!/usr/bin/env python3
"""Fixtures for S9c.2 of REVIEW_NOTES.md: Booleans of prisms whose
cylindrical walls meet in S7b.1's procedural curves (two cylinders of
different radii or skew axes: quartics), before any of its kernel code.

The cases, the frames, the reference (`curved_boolean_reference.py`), the
checks' machinery (inclusion and exclusion with each operation sliced
apart, every face's classes against its closed-form area, the area
identity, a second slicing direction for parallel axes, Monte Carlo) and
the near-coincidence scan are `generate_curved_boolean_fixtures.py`'s,
imported; only its frames are used (their stored axes the kernel's bit for
bit: `boolean-procedural-frames.tsv` records them). S9c.2a's cases lie in
exact frames (`XY`, `TURN`, `DOWN` against `SIDE`: the world's axes
permuted or reversed, crossing axes perpendicular), S9c.2b's in turned ones
(`TILT`, `LEAN`, `R125` against `XY`).

`boolean-procedural-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); `boolean-procedural-expected.tsv`
gives per case `expect KIND STEP` (`solid`, `empty` or `degenerate`, then
`reason TEXT`; `STEP` S9c.2a or S9c.2b) and `result N volume area cx cy cz`
or `empty`, as the curved fixtures.

In the frame of the two axes `a` (the thick cylinder's) and `b`, with `e`
along `a x b`, a point of both walls has `eta = p.e` in both intervals
`[eA - rA, eA + rA]` and `[eB - rB, eB + rB]`; the classes are the
decisions': one strictly inside the other, two rings (`rings`,
`rings_crossing`, `hole_wall`, `oblique`); overlapping in part, one loop
(`bite`, `equal_offset`, `blind`, `skew`); an end on the other's end, a
tangency (`inside_tangent`: a node; `outside_tangent`: a point), declared
`degenerate` for every operation (S9c.1's amendment (c): a tangency between
the inputs is refused whether or not the result involves it). A pin
parallel to a round hole in `R125` (`parallel_hole`) is S9c.2b's too: the
two cylinders are not circular in a common measure.

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* the common of perpendicular cylinders by one quadrature in `eta`
  (mpmath's tanh-sinh between the integrand's kinks):
  `int len([xA - wA, xA + wA] & [bLo, bHi]) len([zB - wB, zB + wB] & [aLo,
  aHi]) d eta` with `wA = sqrt(rA^2 - (eta - eA)^2)` and `wB` alike (for
  cylinders crossing each other's whole section `4 int wA wB d eta`), its
  first moments the same way, and its area from the walls' own angles (the
  length of each generatrix inside the other solid) and the caps' parts
  inside the other, so the three operations' volumes, areas and centres
  (`fuse = A + B - common`, `cut = A - common`, boundaries likewise), for
  every pair of perpendicular cylinders, a box's round hole taken as the box
  less its cylinder;
* Legendre's closed form for crossing axes, `8 rA / 3 ((rA^2 + rB^2) E(k) -
  (rA^2 - rB^2) K(k))`, `k = rB / rA` (mpmath's `ellipe(m)` and `ellipk(m)`
  of the parameter `m = k^2`), against the quadrature and the reference, the
  formula checked at `rA = rB` against `16 r^3 / 3`, the common's centre at
  the axes' crossing by symmetry;
* in turned frames the common of crossing cylinders at angle `phi` as the
  perpendicular one over `sin phi` (each slice a parallelogram of the two
  strips), centred on the axes' common perpendicular, and parallel
  cylinders by their lens (`generate_curved_boolean_fixtures.lens`);

within 1e-30 in exact frames and 1e-15 in turned ones (the stored axes not
exactly orthonormal); and the curved generator's checks per pair, within
1e-30 and 5 standard errors. `generate_curved_boolean_fixtures.py --check`
must still pass (its own fixtures, S9a's and S9b's). The scan for near
coincidences (a class of a face or a result's volume positive but below
1e-9 of the case's size, slicing breakpoints closer than that) must find
none. Pairs run in worker processes (`--workers`).
"""
import argparse
import os
import struct

import mpmath as mp

from curve_surface_reference import stored_axes
import curved_boolean_reference as ref
import generate_curved_boolean_fixtures as curved
from generate_curved_boolean_fixtures import FRAMES, at, disc, group, lens, run, square

ROOT = curved.ROOT
PREFIX = 'boolean-procedural'

SOLID3 = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
NODE = 'two cylinders tangent inside: their meeting has a node'
TOUCH = 'two cylinders tangent outside at a point'
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}


def cases():
    """S9c.2a in exact frames (a thick cylinder along `z`, `XY`, `TURN` or
    `DOWN`, against a pipe along `x`, `SIDE`, offset along `y`): a pipe
    through with offset axes and with crossing axes (two rings; the pipe the
    object in the second), a partial bite and equal radii offset (one loop),
    a pipe ending inside (its cap inside), a box's round hole crossed by a
    pipe (rings on the hole's wall), internal and external tangency. S9c.2b
    in turned frames: oblique crossing axes (`TILT`), skew unequal cylinders
    (`LEAN`), a pin in `R125` parallel to a box's round hole cutting its
    wall."""
    out = []
    thick = ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -4.0, 4.0)
    out += group('rings', thick, ([disc(0.5, 0.0, 1.0)], at('SIDE', (-5, 0, 0)), 0.0, 10.0), SOLID3)
    out += group('rings_crossing', ([disc(0.0, 0.3, 1.5)], at('SIDE', (-6, 0, 0)), 0.0, 12.0),
                 ([disc(0.0, 0.0, 2.5)], at('TURN', (0, 0, 0)), -3.0, 5.0), SOLID3)
    out += group('bite', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -3.0, 3.0),
                 ([disc(1.2, 0.25, 1.5)], at('SIDE', (-4, 0, 0)), 0.0, 8.0), SOLID3)
    out += group('equal_offset', ([disc(0.0, 0.0, 2.0)], at('DOWN', (0, 0, 0)), -3.0, 3.0),
                 ([disc(0.7, -0.2, 2.0)], at('SIDE', (-5, 0, 0)), 0.0, 10.0),
                 {'fuse': 'solid', 'common': 'solid'})
    out += group('blind', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -3.0, 3.0),
                 ([disc(-0.4, 0.3, 1.0)], at('SIDE', (-6, 0, 0)), 0.0, 6.5),
                 {'fuse': 'solid', 'cut': 'solid'})
    out += group('hole_wall', ([square(0.0, 0.0, 10.0, 10.0), disc(5.0, 5.0, 2.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(5.5, 2.1, 0.8)], at('SIDE', (-3, 0, 0)), 0.0, 16.0), SOLID3)
    out += group('inside_tangent', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -3.0, 3.0),
                 ([disc(1.0, 0.0, 1.0)], at('SIDE', (-5, 0, 0)), 0.0, 10.0),
                 {op: ('degenerate', NODE) for op in ('fuse', 'cut', 'common')})
    out += group('outside_tangent', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -3.0, 3.0),
                 ([disc(3.0, 0.0, 1.0)], at('SIDE', (-5, 0, 0)), 0.0, 10.0),
                 {op: ('degenerate', TOUCH) for op in ('fuse', 'cut', 'common')})
    out += group('oblique', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -7.0, 7.0),
                 ([disc(0.0, 0.0, 1.2)], at('TILT', (0, 0, 0.5)), -9.0, 9.0),
                 {'fuse': 'solid', 'common': 'solid'})
    out += group('skew', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -6.0, 6.0),
                 ([disc(0.9, 0.0, 1.3)], at('LEAN', (0, 0, 0.3)), -8.0, 8.0),
                 {'common': 'solid', 'cut': 'solid'})
    out += group('parallel_hole', ([square(0.0, 0.0, 10.0, 10.0), disc(5.0, 5.0, 2.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(1.3, 0.0, 1.1)], at('R125', (5, 5, -1)), 0.0, 7.0),
                 {'cut': 'solid', 'fuse': 'solid'})
    return out


def step(case):
    return 'S9c.2a' if all(curved.frame_name(f) in EXACT for f in case.frames) else 'S9c.2b'


# ------------------------------------------------------------------ closed forms

def quad(f, points):
    """The integral of `f` over the sorted, merged `points` (tanh-sinh on
    each piece: square-root ends are harmless), and its error estimate."""
    pts = sorted(set(points))
    total, err = mp.mpf(0), mp.mpf(0)
    for a, b in zip(pts, pts[1:]):
        if b-a > mp.mpf(10)**-35:
            v, e = mp.quad(f, [a, b], error=True, maxdegree=10)
            total, err = total+v, err+e
    return total, err


def span(c, w, lo, hi):
    """`[c - w, c + w] & [lo, hi]`: its length and first moment."""
    a, b = max(c-w, lo), min(c+w, hi)
    if b <= a:
        return mp.mpf(0), mp.mpf(0)
    return b-a, (b*b-a*a)/2


def half_width(r, e, eta):
    q = r*r-(eta-e)**2
    return mp.sqrt(q) if q > 0 else mp.mpf(0)


def eta_kinks(r, e, c, w0):
    """The `eta` where `c + s w(eta)` (`w` the half-width of radius `r` about
    `e`) crosses `w0`, both signs."""
    out = []
    q = r*r-(w0-c)**2
    if q > 0:
        out += [e-mp.sqrt(q), e+mp.sqrt(q)]
    return out


class Perpendicular:
    """A solid cylinder `A` (radius `rA` about the line `xi = xA`, `eta =
    eA`, its axis along `zeta` in `[aLo, aHi]`) and `B` (radius `rB` about
    `eta = eB`, `zeta = zB`, along `xi` in `[bLo, bHi]`) in world
    coordinates `origin + xi X + eta E + zeta Z`: their common's volume,
    world first moments, and areas (`A`'s wall and caps inside `B`, `B`'s
    inside `A`), each by mpmath quadrature at the working precision."""

    def __init__(self, rA, xA, eA, aLo, aHi, rB, eB, zB, bLo, bHi, frame):
        mpf = mp.mpf
        self.rA, self.xA, self.eA, self.aLo, self.aHi = (mpf(v) for v in (rA, xA, eA, aLo, aHi))
        self.rB, self.eB, self.zB, self.bLo, self.bHi = (mpf(v) for v in (rB, eB, zB, bLo, bHi))
        self.origin, self.X, self.E, self.Z = (tuple(mpf(c) for c in v) for v in frame)
        self.lo = max(self.eA-self.rA, self.eB-self.rB)
        self.hi = min(self.eA+self.rA, self.eB+self.rB)

    def spans(self, eta):
        wA, wB = half_width(self.rA, self.eA, eta), half_width(self.rB, self.eB, eta)
        return span(self.xA, wA, self.bLo, self.bHi), span(self.zB, wB, self.aLo, self.aHi)

    def eta_points(self):
        pts = [self.lo, self.hi]
        for w0 in (self.bLo, self.bHi):
            pts += eta_kinks(self.rA, self.eA, self.xA, w0)
        for w0 in (self.aLo, self.aHi):
            pts += eta_kinks(self.rB, self.eB, self.zB, w0)
        return [p for p in pts if self.lo <= p <= self.hi]

    def volume(self):
        """Volume and world first moments of the common."""
        if self.hi <= self.lo:
            return mp.mpf(0), (mp.mpf(0),)*3, mp.mpf(0)
        pts = self.eta_points()
        out, err = [], mp.mpf(0)
        for k in range(4):
            def f(eta, k=k):
                (l1, m1), (l2, m2) = self.spans(eta)
                return (l1*l2, m1*l2, eta*l1*l2, l1*m2)[k]
            v, e = quad(f, pts)
            out.append(v)
            err = max(err, e)
        V, Mxi, Meta, Mzeta = out
        M = tuple(self.origin[i]*V+self.X[i]*Mxi+self.E[i]*Meta+self.Z[i]*Mzeta for i in range(3))
        return V, M, err

    def areas(self):
        """(A's boundary inside B, B's boundary inside A): walls and caps."""
        rA, rB = self.rA, self.rB
        # A's wall: xi = xA + rA cos t, eta = eA + rA sin t.
        tp = [mp.mpf(0), 2*mp.pi]
        for w0 in (self.bLo, self.bHi):
            c = (w0-self.xA)/rA
            if abs(c) <= 1:
                tp += [mp.acos(c), 2*mp.pi-mp.acos(c)]
        for eta in [self.eB-rB, self.eB+rB]+[p for w0 in (self.aLo, self.aHi)
                                              for p in eta_kinks(rB, self.eB, self.zB, w0)]:
            s = (eta-self.eA)/rA
            if abs(s) <= 1:
                a = mp.asin(s)
                tp += [a % (2*mp.pi), (mp.pi-a) % (2*mp.pi)]

        def a_wall(t):
            xi, eta = self.xA+rA*mp.cos(t), self.eA+rA*mp.sin(t)
            if not self.bLo < xi < self.bHi:
                return mp.mpf(0)
            return rA*span(self.zB, half_width(rB, self.eB, eta), self.aLo, self.aHi)[0]
        # B's wall: eta = eB + rB cos p, zeta = zB + rB sin p.
        pp = [mp.mpf(0), 2*mp.pi]
        for w0 in (self.aLo, self.aHi):
            s = (w0-self.zB)/rB
            if abs(s) <= 1:
                a = mp.asin(s)
                pp += [a % (2*mp.pi), (mp.pi-a) % (2*mp.pi)]
        for eta in [self.eA-rA, self.eA+rA]+[p for w0 in (self.bLo, self.bHi)
                                              for p in eta_kinks(rA, self.eA, self.xA, w0)]:
            c = (eta-self.eB)/rB
            if abs(c) <= 1:
                pp += [mp.acos(c), 2*mp.pi-mp.acos(c)]

        def b_wall(p):
            eta, zeta = self.eB+rB*mp.cos(p), self.zB+rB*mp.sin(p)
            if not self.aLo < zeta < self.aHi:
                return mp.mpf(0)
            return rB*span(self.xA, half_width(rA, self.eA, eta), self.bLo, self.bHi)[0]
        pts = self.eta_points()
        extra = []
        for w0 in (self.bLo, self.bHi):
            extra += eta_kinks(rA, self.eA, self.xA, w0)
        for w0 in (self.aLo, self.aHi):
            extra += eta_kinks(rB, self.eB, self.zB, w0)
        pts = [p for p in pts+extra if self.lo <= p <= self.hi]

        def b_caps(eta):
            wA = half_width(rA, self.eA, eta)
            n = sum(1 for c in (self.bLo, self.bHi) if abs(c-self.xA) < wA)
            return n*span(self.zB, half_width(rB, self.eB, eta), self.aLo, self.aHi)[0]

        def a_caps(eta):
            wB = half_width(rB, self.eB, eta)
            n = sum(1 for c in (self.aLo, self.aHi) if abs(c-self.zB) < wB)
            return n*span(self.xA, half_width(rA, self.eA, eta), self.bLo, self.bHi)[0]
        in_b = quad(a_wall, tp)[0]+(quad(a_caps, pts)[0] if self.hi > self.lo else 0)
        in_a = quad(b_wall, pp)[0]+(quad(b_caps, pts)[0] if self.hi > self.lo else 0)
        return in_b, in_a


def legendre(rA, rB):
    """The common of two crossing perpendicular cylinders, `rA >= rB` (at
    equal radii `K(1)` is infinite and its factor zero: `E(1) = 1` alone)."""
    rA, rB = mp.mpf(rA), mp.mpf(rB)
    m = (rB/rA)**2
    first = 8*rA/3*(rA*rA+rB*rB)*mp.ellipe(m)
    return first if rA == rB else first-8*rA/3*(rA*rA-rB*rB)*mp.ellipk(m)


def cylinder(r, lo, hi, centre):
    """Volume, world first moments and area of a solid cylinder, `centre`
    the midpoint of its axis."""
    r, h = mp.mpf(r), mp.mpf(hi)-mp.mpf(lo)
    V = mp.pi*r*r*h
    return V, tuple(V*mp.mpf(c) for c in centre), 2*mp.pi*r*r+2*mp.pi*r*h


def box_with_hole(side, height, cx, cy, r):
    """The fixtures' box `[0, side]^2 x [0, height]` less a round hole."""
    side, height, r = mp.mpf(side), mp.mpf(height), mp.mpf(r)
    V = (side*side-mp.pi*r*r)*height
    Mx = side*side*height*side/2-mp.pi*r*r*height*mp.mpf(cx)
    My = side*side*height*side/2-mp.pi*r*r*height*mp.mpf(cy)
    S = 2*(side*side-mp.pi*r*r)+4*side*height+2*mp.pi*r*height
    return V, (Mx, My, V*height/2), S


def ops_from_common(A, B, common, cut_object='A'):
    """{op: (volume, moments, area)} from both inputs' measures and the
    common's (volume, moments, area of A's boundary inside B, area of B's
    boundary inside A)."""
    (va, ma, sa), (vb, mb, sb) = A, B
    vc, mc, a_in, b_in = common
    out = {'common': (vc, mc, None if a_in is None else a_in+b_in),
           'fuse': (va+vb-vc, tuple(x+y-z for x, y, z in zip(ma, mb, mc)),
                    None if a_in is None else sa+sb-a_in-b_in)}
    if cut_object == 'A':
        out['cut'] = (va-vc, tuple(x-z for x, z in zip(ma, mc)), None if a_in is None else sa-a_in+b_in)
    else:
        out['cut'] = (vb-vc, tuple(x-z for x, z in zip(mb, mc)), None if a_in is None else sb-b_in+a_in)
    return out


WORLD = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def perpendicular_forms():
    """{pair: ({op: (volume, moments, area)}, quadrature error, crossing
    (rA, rB) or None)} for S9c.2a's pairs (`xi = x`, `eta = y`, `zeta = z`;
    the pipes' `SIDE` frames put their `(u, v)` on `(y, z)`, their height
    on `x` from the frame's origin)."""
    out = {}

    def pair(name, thick, pipe, cut_object='A', crossing=None):
        (rA, aLo, aHi), (rB, eB, zB, bLo, bHi) = thick, pipe
        P = Perpendicular(rA, 0, 0, aLo, aHi, rB, eB, zB, bLo, bHi, WORLD)
        vc, mc, err = P.volume()
        a_in, b_in = P.areas()
        A = cylinder(rA, aLo, aHi, (0, 0, (mp.mpf(aLo)+aHi)/2))
        B = cylinder(rB, bLo, bHi, ((mp.mpf(bLo)+bHi)/2, eB, zB))
        if cut_object == 'B':
            # The pipe is the object: cut = pipe - common.
            A, B, a_in, b_in = B, A, b_in, a_in
        out[name] = (ops_from_common(A, B, (vc, mc, a_in, b_in)), err, crossing, P)

    pair('rings', (2, -4, 4), (1, 0.5, 0, -5, 5))
    # The pipe is the object; the thick cylinder (TURN) the tool.
    pair('rings_crossing', (2.5, -3, 5), (1.5, 0, 0.3, -6, 6), cut_object='B', crossing=(2.5, 1.5))
    pair('bite', (2, -3, 3), (1.5, 1.2, 0.25, -4, 4))
    # DOWN: the thick cylinder's heights [-3, 3] are z in [-3, 3].
    pair('equal_offset', (2, -3, 3), (2, 0.7, -0.2, -5, 5))
    pair('blind', (2, -3, 3), (1, -0.4, 0.3, -6, 0.5))
    pair('inside_tangent', (2, -3, 3), (1, 1, 0, -5, 5))
    # A box with a round hole crossed by a pipe: the box less the hole's
    # cylinder H (radius 2 about x = y = 5, z in [0, 4]).
    H = Perpendicular(2, 5, 5, 0, 4, 0.8, 5.5, 2.1, -3, 13, WORLD)
    vh, mh, err = H.volume()
    h_in, p_in = H.areas()
    rB, L = mp.mpf(0.8), 10
    vbox = mp.pi*rB*rB*L
    mbox = (vbox*5, vbox*mp.mpf(5.5), vbox*mp.mpf(2.1))
    vc = vbox-vh
    mc = tuple(x-y for x, y in zip(mbox, mh))
    A = box_with_hole(10, 4, 5, 5, 2)
    B = cylinder(rB, -3, 13, (5, mp.mpf(5.5), mp.mpf(2.1)))
    a_in = 2*mp.pi*rB*rB+h_in                      # the box's side faces and the hole's wall inside the pipe
    b_in = 2*mp.pi*rB*L-p_in                       # the pipe's wall in the box, less inside the hole
    out['hole_wall'] = (ops_from_common(A, B, (vc, mc, a_in, b_in)), err, None, H)
    return out


def turned_forms():
    """{pair: {op: (volume, moments, area or None)}} for S9c.2b's pairs, from
    the ideal orthonormal frames: crossing cylinders at `sin phi = 3/5` (the
    perpendicular common over `sin phi`, centred on the axes' crossing
    along their common perpendicular), parallel ones by their lens."""
    out = {}
    s = mp.mpf(3)/5
    # oblique: XY r 2 over z in [-7, 7]; TILT r 1.2 about (0, 0, 0.5) along
    # (0, 3, 4)/5 over [-9, 9]; the common perpendicular along x.
    P = Perpendicular(2, 0, 0, -50, 50, 1.2, 0, 0, -50, 50, WORLD)
    vp, mp_, _ = P.volume()
    vc = vp/s
    eta = mp_[1]/vp                                # y moment is the eta moment here (E = y)
    centre = (eta, mp.mpf(0), mp.mpf(0.5))
    A = cylinder(2, -7, 7, (0, 0, 0))
    B = cylinder(1.2, -9, 9, (0, 0, mp.mpf(0.5)))
    common = (vc, tuple(vc*c for c in centre), None, None)
    out['oblique'] = (ops_from_common(A, B, common), (2, 1.2))
    # skew: XY r 2 over z in [-6, 6]; LEAN r 1.3 about (0, 0.9, 0.3) along
    # (3, 0, 4)/5 over [-8, 8]; the common perpendicular along y.
    P = Perpendicular(2, 0, 0, -50, 50, 1.3, 0.9, 0, -50, 50, WORLD)
    vp, mp_, _ = P.volume()
    vc = vp/s
    centre = (mp.mpf(0), mp_[1]/vp, mp.mpf(0.3))
    A = cylinder(2, -6, 6, (0, 0, 0))
    B = cylinder(1.3, -8, 8, (0, mp.mpf(0.9), mp.mpf(0.3)))
    out['skew'] = (ops_from_common(A, B, (vc, tuple(vc*c for c in centre), None, None)), None)
    # parallel_hole: the box [0, 10]^2 x [0, 4] less a hole r 2 about (5, 5);
    # a pin r 1.1 about (6.2, 5.5) (1.3 along (12, 5)/13), z in [-1, 6].
    d, r1, r2 = mp.mpf(1.3), mp.mpf(2), mp.mpf(1.1)
    la, lg, _ = lens(r1, r2, d)
    xc = (d*d+r1*r1-r2*r2)/(2*d)
    a1, a2 = mp.acos(xc/r1), mp.acos((d-xc)/r2)
    u = (mp.mpf(12)/13, mp.mpf(5)/13)
    pin = (5+d*u[0], 5+d*u[1])
    area_c = mp.pi*r2*r2-la                        # the pin's section in the box's material
    m2 = [mp.pi*r2*r2*pin[i]-la*(5+lg*u[i]) for i in range(2)]
    vc = 4*area_c
    mc = (4*m2[0], 4*m2[1], vc*2)
    A = box_with_hole(10, 4, 5, 5, 2)
    B = cylinder(r2, -1, 6, (pin[0], pin[1], mp.mpf(2.5)))
    a_in = 2*area_c+4*2*r1*a1                      # the box's caps and the hole's wall inside the pin
    b_in = 4*(2*mp.pi*r2-2*r2*a2)                  # the pin's wall outside the hole, within the box's height
    out['parallel_hole'] = (ops_from_common(A, B, (vc, mc, a_in, b_in)), None)
    return out


# ------------------------------------------------------------------ checks

def deviation(res, want):
    n, vol, area, centre = res
    V, Mm, S = want
    size = max(abs(V)**(mp.mpf(1)/3), 1)
    dev = abs(vol-V)/abs(V)
    if S is not None:
        dev = max(dev, abs(area-S)/abs(S))
    return max(dev, max(abs(centre[i]-Mm[i]/V) for i in range(3))/size)


def reference_checks(results):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        assert not mp.isnan(value), key
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    by_pair = {r[0]: r for r in results}
    ops_of = {}
    for c in cases():
        ops_of.setdefault(c.pair_name, []).append(c.operation)
    # Legendre's form, and its value at equal radii.
    for r in (1, 2, mp.mpf(1.5)):
        note('legendre_equal_radii', abs(legendre(r, r)-16*mp.mpf(r)**3/3)/mp.mpf(r)**3)
    for name, (forms, err, crossing, P) in perpendicular_forms().items():
        note('quadrature_estimate', err/max(1, P.rA)**3)
        _, _, _, res, _, _, _ = by_pair[name]
        if crossing is not None:
            L = legendre(*crossing)
            note('legendre_against_quadrature', abs(L-forms['common'][0])/L)
            note('closed_forms_exact_frames', abs(res['common'][1]-L)/L)
            # By symmetry the common's centre is the axes' crossing.
            n, vol, area, centre = res['common']
            note('closed_forms_exact_frames', max(abs(centre[0]), abs(centre[1]), abs(centre[2]-mp.mpf(0.3))))
        for op in ref.OPS:
            if res[op][0] == 0:
                continue
            note('closed_forms_exact_frames', deviation(res[op], forms[op]))
    for name, (forms, crossing) in turned_forms().items():
        _, _, _, res, _, _, _ = by_pair[name]
        if crossing is not None:
            L = legendre(*crossing)/(mp.mpf(3)/5)
            note('closed_forms_turned_frames', abs(res['common'][1]-L)/L)
        for op in ops_of[name]+['common']:
            note('closed_forms_turned_frames', deviation(res[op], forms[op]))
    for name, rows, _, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = ['# case\trow (S9c.2, curved_boolean_reference.py: expect KIND STEP, reason TEXT for a degenerate '
           'case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, res, _, near, _ = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {step(case)}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (obj, tool) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in cases():
        for label, c in (('obj', case.obj), ('tool', case.tool)):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\t{label} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c.obj, c.tool, [], c.opposite])[2].append(c.operation)
    return [(name, obj, tool, ops, opposite, mc_n) for name, (obj, tool, ops, opposite) in pairs.items()]


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
          'legendre_equal_radii': 1e-35, 'legendre_against_quadrature': 1e-32, 'quadrature_estimate': 1e-30,
          'inclusion_exclusion': 1e-30, 'face_classes': 1e-30, 'shared_both_sides': 1e-30,
          'area_identity': 1e-30, 'second_axis': 1e-30, 'monte_carlo_sigma': 5,
          'volume_quadrature_estimate': 1e-30, 'face_quadrature_estimate': 1e-30}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    args = parser.parse_args()
    listed = cases()
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert all(curved.frame_name(f) in FRAMES for f in c.frames)
        assert stored_axes(c.obj.frame)[1:] != stored_axes(c.tool.frame)[1:], f'{c.name}: frames with equal axes'
        assert any(b.circle is not None for b in c.obj.boundaries) and \
            any(b.circle is not None for b in c.tool.boundaries), f'{c.name}: not two cylinders'
    results = run(jobs(args.samples), curved.evaluate, args.workers)
    for name, _, _, _, _, near, _ in results:
        if near:
            raise SystemExit(f'{name}: near coincidences {near}')
    worst, covered = reference_checks(results)
    for key, value in worst.items():
        assert value <= LIMITS[key], (key, mp.nstr(value, 3))
    files = generate(results)
    for name, contents in files.items():
        target = ROOT/'fixtures'/name
        if args.check:
            if target.read_text() != contents:
                raise SystemExit(f'{target} is stale')
        else:
            target.write_text(contents)
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{sum(1 for c in listed if step(c) == s)} {s}' for s in ('S9c.2a', 'S9c.2b')))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    by_pair = {r[0]: r for r in results}
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
