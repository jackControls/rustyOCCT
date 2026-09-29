#!/usr/bin/env python3
"""Fixtures for S9d.4b.2 of REVIEW_NOTES.md: Booleans of a whole torus
(`Solid::torus_with`, the full tube and turn) against a prism with arcs
(cylinders; one stadium), a whole sphere (`Solid::sphere_with`), a cone or
frustum (`Solid::cone_with`) and another whole torus, before any of its
kernel code.

`boolean-torus-curved-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a torus's block is its `frame` and
`torus MAJOR MINOR LOW HIGH ANGLE` rows, a sphere's `frame` and `sphere R LOW
HIGH`, a cone's `frame` and `cone BOTTOM TOP HEIGHT`, either input or both);
`boolean-torus-curved-expected.tsv` gives per case `expect KIND S9d.4b.2`
(`solid`, `empty` or `degenerate`, then `reason TEXT`) and `result N volume
area cx cy cz` or `empty` from `torus_curved_boolean_reference.py`, as
S9d.4b.1's; `boolean-torus-curved-frames.tsv` the stored axes' bits. Frames
are the curved generator's.

`T` is the torus of radii 5/2 and 1 on `XY` at the origin (its hole of
radius 3/2, its outer equator of radius 7/2) unless said. Coaxial pairs,
meeting in circles: a cylinder of radius 7/4 through the hole cutting the
tube's inner side (`pipe_hole`, all three operations; the pipe first,
`hole_pipe`), a cylinder of radius 13/4 around it between the heights -1/2
and 2 (`sleeve`: its wall and its lower cap cut the tube), a sphere of
radius 11/4 about the centre (`ball_coax`) and one of radius 9/4 on the axis
at the height 3/2 (`ball_axis`), a cone of radius 4 and height 4 from below
through the hole (`cone_coax`), a torus of radii 13/4 and 3/4 half a unit
up (`tori_coax`). Other pairs: a rod of radius 1/2 across the whole torus
along `x` (`SIDE`; its two crossings of the tube, the common two solids,
`rod_ends` the rod less the torus three), a vertical cylinder through the
tube (`bore`: a graph over the cylinder's angle), a stadium prism across the
tube cutting it through (`stadium`: flat walls and arcs), a sphere on the
tube's top (`ball_top`, all three), a sphere of radius 5/4 on the core
circle swallowing the tube there (`ball_core`, and `core_ball` the sphere
first), a frustum standing through the tube (`spike`: parallel axes), a cone
in `LEAN` whose apex lies inside the tube (`cone_lean`), a torus about `x`
(`SIDE`) through the hole linked with `T` without meeting it (`tori_link`:
the fuse two solids, the common empty), a small torus ringing the tube and
pressing into it (`tori_ring`, all three), a torus of parallel axis
overlapping the tube twice (`tori_side`: the common two solids), a torus
about `x` linked with `T` whose tube crosses `T`'s once (`tori_cross`);
apart, a cylinder in
the hole (`pin_hole`) and a sphere in the hole (`ball_hole`); turned frames,
the pipe with the torus in `TILT` (`pipe_tilt`) and a sphere on the tube of
a torus in `TILT` (`ball_tilt`). `degenerate` with reasons: a coaxial
cylinder tangent to the outer equator (`pipe_equator`), a sphere touching
the tube from outside at a point (`ball_touch`), a coaxial torus touching
`T` along its outer equator (`tori_kiss`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms by one quadrature (mpmath's tanh-sinh between the kinks)
  along the torus's axis in its ideal frame, each slice the torus's annulus
  (S9d.4a's) against the other input's section: discs or annuli about a
  parallel axis (a cylinder, a sphere anywhere, a cone, a parallel or
  coaxial torus) by lenses of signed discs (areas, first moments, each
  circle's angle inside the other: the torus's wall by `r rho / q` times
  its circles' angles, a cylinder's by its radius, a sphere's by
  Archimedes, a cone's by `rho sqrt(1 + k^2)`, a torus's by its own
  element; caps by their discs' lenses at their levels), a rod across the
  axis by its strip, a rectangle `[x0, x1] x [y0 - eta, y0 + eta]` against
  the annulus (S9d.1's rectangle-in-disc antiderivatives; the rod's wall by
  its lines' chords in the annulus times `a / eta`, its caps by their
  chords), inputs apart by their sums; within 1e-30 in exact frames, 1e-15
  in turned ones;
* every operation's volume, moments and area two ways (each face's two
  families of curves), both inputs' measures from their faces against their
  closed forms, inclusion and exclusion (`fuse = A + B - common`, `cut = A -
  common`), the area identity `area(fuse) + area(common) = area(A) +
  area(B)`, every face's two families' totals equal, within 1e-30; every
  quadrature node's structure its interval's (no event missed); the
  declared degenerate pairs within 1e-18 and their closed forms within
  1e-20 (a tangency along a circle: the curves through it meet the other
  surface in double roots, and the closed forms' lenses round there);
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* the solids by the sweep of the torus's normal slices, the fuse's by the
  rule (one where the inputs overlap);
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, a family's events closer than
  that) must find none but in the declared degenerate pairs; every other
  pair's margins (the least sine between surfaces where they meet, the
  least distance to a surface a face never crosses, the least sine of an
  edge crossing the torus and its least distance where stationary, the
  vertices' distances) are reported.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import math
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, encode_boolean_case
from curve_surface_reference import stored_axes
import torus_curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, stadium
from generate_sphere_boolean_fixtures import arc_in_rect, frame_name, ideal_axes, rect_disc
from generate_spheres_boolean_fixtures import lens
from generate_cones_boolean_fixtures import roots_on

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-torus-curved'
STEP = 'S9d.4b.2'
BOOLEAN_OPERATION = 93
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
TWO_PI = ref.TWO_PI
HP = math.pi/2
Z = mp.mpf(0)


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def circle(cx, cy, r):
    return Boundary(circle=(float(cx), float(cy), float(r)))


def make(spec, op):
    kind = spec[0]
    if kind == 'torus':
        _, R, r, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], torus=(float(R), float(r), 0.0, TWO_PI, TWO_PI))
    if kind == 'sphere':
        _, r, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), -HP, HP))
    if kind == 'cone':
        _, r0, r1, h, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def torus(R, r, frame):
    return ('torus', R, r, frame)


def sphere(r, frame):
    return ('sphere', r, frame)


def cone(r0, r1, h, frame):
    return ('cone', r0, r1, h, frame)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def spec_frame(spec):
    return {'torus': 3, 'sphere': 2, 'cone': 4, 'prism': 2}[spec[0]]


class Boolean:
    def __init__(self, name, operation, obj, tool, form, kind='solid', reason=None):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.form = form
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.specs = (obj, tool)
        self.frames = (obj[spec_frame(obj)], tool[spec_frame(tool)])

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, form):
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, form, kind, reason))
    return out


EQUATOR = 'a coaxial cylinder tangent to the torus along its outer equator'
TOUCH = 'a sphere tangent to the torus at a point'
KISS = 'two coaxial tori tangent along a circle'
ALL = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
CC = {'cut': 'solid', 'common': 'solid'}


def cases():
    """S9d.4b.2's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    XY = at('XY', O)
    T = torus(2.5, 1, XY)
    out = []
    pipe = prism([circle(0, 0, 1.75)], XY, -2, 2)
    out += group('pipe_hole', T, pipe, ALL, 'slices')
    out += group('hole_pipe', pipe, T, {'cut': 'solid'}, 'slices')
    out += group('sleeve', T, prism([circle(0, 0, 3.25)], XY, -0.5, 2), {'common': 'solid'}, 'slices')
    out += group('ball_coax', T, sphere(2.75, XY), CC, 'slices')
    out += group('ball_axis', T, sphere(2.25, at('XY', (0, 0, 1.5))), {'common': 'solid'}, 'slices')
    out += group('cone_coax', T, cone(4, 0, 4, at('XY', (0, 0, -2))), CC, 'slices')
    out += group('tori_coax', T, torus(3.25, 0.75, at('XY', (0, 0, 0.5))), {'fuse': 'solid', 'common': 'solid'},
                 'slices')
    rod = prism([circle(0, 0, 0.5)], at('SIDE', (-4.5, 0, 0)), 0, 9)
    out += group('rod', T, rod, CC, 'strip')
    out += group('rod_ends', rod, T, {'cut': 'solid'}, 'strip')
    out += group('bore', T, prism([circle(2.75, 0.25, 0.5)], XY, -2, 2), CC, 'slices')
    out += group('stadium', T, prism([stadium(0.75, 0.5)], at('XY', (2.5, 0, -2)), 0, 4), CC, None)
    out += group('ball_top', T, sphere(0.75, at('XY', (2.5, 0, 1))), ALL, 'slices')
    core = sphere(1.25, at('XY', (0, 2.5, 0)))
    out += group('ball_core', T, core, {'cut': 'solid'}, 'slices')
    out += group('core_ball', core, T, {'cut': 'solid'}, 'slices')
    out += group('spike', T, cone(0.75, 0.25, 3, at('XY', (2.75, 0.25, -1.5))), CC, 'slices')
    out += group('cone_lean', T, cone(1, 0, 3, at('LEAN', (0.75, 0.25, -1.75))), {'common': 'solid'}, None)
    out += group('tori_link', T, torus(2.5, 0.75, at('SIDE', (0, 2.5, 0))), {'fuse': 'solid', 'common': 'empty'},
                 'apart')
    out += group('tori_ring', T, torus(1.25, 0.5, at('SIDE', (0, 2.5, 0))), ALL, None)
    out += group('tori_side', T, torus(2.5, 0.75, at('XY', (2.5, 0, 0.125))), CC, 'slices')
    out += group('tori_cross', T, torus(2.5, 0.75, at('SIDE', (0, 4.25, 0))), CC, None)
    out += group('pin_hole', T, prism([circle(0, 0, 1.25)], XY, -2, 2), {'fuse': 'solid', 'common': 'empty'},
                 'apart')
    out += group('ball_hole', T, sphere(1.25, XY), {'fuse': 'solid'}, 'apart')
    out += group('pipe_equator', T, prism([circle(0, 0, 3.5)], XY, -2, 2), {'fuse': ('degenerate', EQUATOR)},
                 'slices')
    out += group('ball_touch', T, sphere(0.5, at('XY', (4, 0, 0))), {'fuse': ('degenerate', TOUCH)}, 'slices')
    out += group('tori_kiss', T, torus(4.25, 0.75, XY), {'fuse': ('degenerate', KISS)}, 'slices')
    TILT = at('TILT', O)
    out += group('pipe_tilt', torus(2.5, 1, TILT), prism([circle(0, 0, 1.75)], TILT, -2, 2), {'common': 'solid'},
                 'slices')
    out += group('ball_tilt', torus(2.5, 1, TILT), sphere(0.75, (2.5, 0.6, 0.8)+FRAMES['XY']), {'cut': 'solid'},
                 'slices')
    return out


def exact_pair(case):
    """Whether the closed forms hold exactly: no torus, cone or prism in a
    turned frame (a whole sphere's model does not depend on its frame's
    axes)."""
    return all(frame_name(f) in EXACT for f, s in zip(case.frames, case.specs) if s[0] != 'sphere')


# ------------------------------------------------------------------ closed forms

class Ideal:
    """An input's ideal frame: origin and orthonormal axes (mpf)."""

    def __init__(self, frame):
        self.o = tuple(mp.mpf(v) for v in frame[:3])
        self.x, self.y, self.n = ideal_axes(frame_name(frame))

    def world(self, u, v, w):
        return tuple(self.o[i]+u*self.x[i]+v*self.y[i]+w*self.n[i] for i in range(3))

    def local(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        return tuple(sum(d[i]*a[i] for i in range(3)) for a in (self.x, self.y, self.n))

    def vec_local(self, D):
        return tuple(sum(D[i]*a[i] for i in range(3)) for a in (self.x, self.y, self.n))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def measures(spec):
    """An input's volume, world moments and area (ideal frame)."""
    kind = spec[0]
    if kind == 'torus':
        _, R, r, frame = spec
        R, r = mp.mpf(R), mp.mpf(r)
        fr = Ideal(frame)
        V = 2*mp.pi**2*R*r*r
        return V, tuple(c*V for c in fr.o), 4*mp.pi**2*R*r
    if kind == 'sphere':
        _, r, frame = spec
        r = mp.mpf(r)
        V = 4*mp.pi*r**3/3
        return V, tuple(mp.mpf(c)*V for c in frame[:3]), 4*mp.pi*r*r
    if kind == 'cone':
        _, r0, r1, h, frame = spec
        r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
        fr = Ideal(frame)
        V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3
        zc = h*(r0*r0+2*r0*r1+3*r1*r1)/(4*(r0*r0+r0*r1+r1*r1))
        A = mp.pi*(r0+r1)*mp.sqrt(h*h+(r1-r0)**2)+mp.pi*(r0*r0+r1*r1)
        return V, tuple(c*V for c in fr.world(Z, Z, zc)), A
    _, boundaries, frame, lo, hi = spec
    fr = Ideal(frame)
    cx, cy, a = (mp.mpf(v) for v in boundaries[0].circle)
    h = mp.mpf(hi)-mp.mpf(lo)
    V = mp.pi*a*a*h
    return V, tuple(c*V for c in fr.world(cx, cy, (mp.mpf(lo)+mp.mpf(hi))/2)), 2*mp.pi*a*h+2*mp.pi*a*a


class Sections:
    """The other input's sections by the torus's slices `w = s` (its ideal
    chart), each a signed set of discs about one centre, its wall's circles
    with their elements, and its caps (level, centre, radius)."""

    def __init__(self, spec, K):
        self.kind = spec[0]
        kind = spec[0]
        if kind == 'sphere':
            _, r, frame = spec
            self.c = K.local(tuple(mp.mpf(v) for v in frame[:3]))
            self.r = mp.mpf(r)
            self.range = (self.c[2]-self.r, self.c[2]+self.r)
            self.caps = []
            return
        fr = Ideal(frame := spec[spec_frame(spec)])
        sg = dot(fr.n, K.n)
        assert abs(abs(sg)-1) < mp.mpf(10)**-30, 'the sections are discs about a parallel axis'
        self.sg = mp.mpf(1) if sg > 0 else mp.mpf(-1)
        if kind == 'torus':
            _, R, r, _ = spec
            self.c = K.local(fr.o)
            self.R, self.r = mp.mpf(R), mp.mpf(r)
            self.range = (self.c[2]-self.r, self.c[2]+self.r)
            self.caps = []
        elif kind == 'cone':
            _, r0, r1, h, _ = spec
            self.c = K.local(fr.o)
            self.r0, self.r1, self.h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
            self.k = (self.r1-self.r0)/self.h
            ends = (self.c[2], self.c[2]+self.sg*self.h)
            self.range = (min(ends), max(ends))
            self.caps = [(ends[0], self.c, self.r0), (ends[1], self.c, self.r1)]
        else:
            _, boundaries, _, lo, hi = spec
            cx, cy, a = (mp.mpf(v) for v in boundaries[0].circle)
            self.c = K.local(fr.world(cx, cy, Z))
            self.a = a
            ends = (self.c[2]+self.sg*mp.mpf(lo), self.c[2]+self.sg*mp.mpf(hi))
            self.range = (min(ends), max(ends))
            self.caps = [(e, self.c, a) for e in ends]

    def at(self, s):
        """([(sign, radius)], [(radius, element)]) at level `s`, about `c`."""
        lo, hi = self.range
        if not lo < s < hi:
            return [], []
        if self.kind == 'sphere':
            rho = mp.sqrt(self.r**2-(s-self.c[2])**2)
            return [(1, rho)], [(rho, self.r)]
        if self.kind == 'torus':
            q = mp.sqrt(self.r**2-(s-self.c[2])**2)
            return [(1, self.R+q), (-1, self.R-q)], [(self.R+q, self.r*(self.R+q)/q), (self.R-q, self.r*(self.R-q)/q)]
        if self.kind == 'cone':
            t = self.sg*(s-self.c[2])
            rho = self.r0+self.k*t
            return [(1, rho)], [(rho, rho*mp.sqrt(1+self.k**2))]
        return [(1, self.a)], [(self.a, self.a)]

    def radii_funcs(self):
        """The circles' radii as functions of the level (for the kinks)."""
        c = self.c[2]
        if self.kind == 'sphere':
            return [lambda s: mp.sqrt(max(self.r**2-(s-c)**2, Z))]
        if self.kind == 'torus':
            q = lambda s: mp.sqrt(max(self.r**2-(s-c)**2, Z))
            return [lambda s: self.R+q(s), lambda s: self.R-q(s)]
        if self.kind == 'cone':
            return [lambda s: self.r0+self.k*self.sg*(s-c)]
        return [lambda s: self.a]


def annulus(R, r, s):
    if abs(s) >= r:
        return [], []
    q = mp.sqrt(r*r-s*s)
    return [(1, R+q), (-1, R-q)], [(R+q, r*(R+q)/q), (R-q, r*(R-q)/q)]


def disc_lens(r1, c1, r2, c2):
    """Lens of the discs `(c1, r1)` and `(c2, r2)` (2D centres): area, first
    moments (vector), the first circle's angle inside the second disc and the
    second's inside the first."""
    d = mp.sqrt((c2[0]-c1[0])**2+(c2[1]-c1[1])**2)
    A, m, a1, a2 = lens(r1, r2, d)
    if d > 0:
        u = ((c2[0]-c1[0])/d, (c2[1]-c1[1])/d)
    else:
        u = (Z, Z)
    return A, (A*c1[0]+m*u[0], A*c1[1]+m*u[1]), a1, a2


def slice_levels(K, P, kinks):
    """Kinks of the slices: the sections' ends, and levels where a circle of
    the torus and one of the other are tangent (sign changes and near zeros
    of `d - (r1 + r2)` and `d - |r1 - r2|` on a fine scan, bisected)."""
    R, r = K.R, K.r
    lo = max(-r, P.range[0])
    hi = min(r, P.range[1])
    pts = {-r, r, P.range[0], P.range[1]}
    if hi <= lo:
        return sorted(pts), lo, hi
    d = mp.sqrt(P.c[0]**2+P.c[1]**2)
    q = lambda s: mp.sqrt(max(r*r-s*s, Z))
    kf = [lambda s: R+q(s), lambda s: R-q(s)]
    for f1 in kf:
        for f2 in P.radii_funcs():
            for g in (lambda s, f1=f1, f2=f2: d-(f1(s)+f2(s)), lambda s, f1=f1, f2=f2: d-abs(f1(s)-f2(s))):
                for s in roots_on(g, lo, hi, n=2000):
                    pts.add(s)
    for level, _, _ in P.caps:
        pts.add(level)
    pts |= set(kinks)
    return sorted(pts), lo, hi


def slices_form(case):
    """The torus against an input whose slices are discs about a parallel
    axis, by one quadrature along the torus's axis of exact lenses."""
    kspec, pspec = (case.specs[0], case.specs[1]) if case.specs[0][0] == 'torus' else (case.specs[1],
                                                                                      case.specs[0])
    if case.specs[0][0] == 'torus' and case.specs[1][0] == 'torus':
        kspec, pspec = case.specs
    k_first = case.specs[0] is kspec
    K = Ideal(kspec[3])
    R, r = mp.mpf(kspec[1]), mp.mpf(kspec[2])
    K.R, K.r = R, r
    P = Sections(pspec, K)
    pts, lo, hi = slice_levels(K, P, [])
    origin = (Z, Z)
    pc = (P.c[0], P.c[1])

    def g(s):
        kd, kw = annulus(R, r, s)
        pd, pw = P.at(s)
        A = mu = mv = Z
        for sk, rk in kd:
            for sp, rp in pd:
                a, m, _, _ = disc_lens(rk, origin, rp, pc)
                A += sk*sp*a
                mu += sk*sp*m[0]
                mv += sk*sp*m[1]
        pin = Z
        if P.kind in ('sphere', 'cone', 'prism'):
            for rp, el in pw:
                ang = sum(sk*disc_lens(rp, pc, rk, origin)[2] for sk, rk in kd)
                pin += el*ang
        return A, mu, mv, s*A, pin
    inner = [p for p in pts if lo <= p <= hi]
    vals = []
    for k in range(5):
        if len(inner) < 2:
            vals.append(Z)
        else:
            vals.append(mp.quad(lambda s: g(s)[k], inner))
    V, Mu, Mv_, Mw, pin = vals

    # The walls of tori by their latitudes (`s = r sin phi`, the element `r
    # rho dphi` smooth where `r rho / q` is not).
    def wall(c2, R2, r2, inside, levels):
        phis = {-mp.pi/2, mp.pi/2, 3*mp.pi/2}
        for s in levels:
            x = (s-c2)/r2
            if -1 < x < 1:
                p0 = mp.asin(x)
                phis |= {p0, mp.pi-p0}
        phis = sorted(phis)

        def f(phi):
            s = c2+r2*mp.sin(phi)
            rho = R2+r2*mp.cos(phi)
            return r2*rho*inside(rho, s)
        return mp.quad(f, phis)

    def k_inside(rho, s):
        pd, _ = P.at(s)
        return sum(sp*disc_lens(rho, origin, rp, pc)[2] for sp, rp in pd)
    kin = wall(Z, R, r, k_inside, pts)
    if P.kind == 'torus':
        def p_inside(rho, s):
            kd, _ = annulus(R, r, s)
            return sum(sk*disc_lens(rho, pc, rk, origin)[2] for sk, rk in kd)
        pin = wall(P.c[2], P.R, P.r, p_inside, pts)
    for level, c, rad in P.caps:
        kd, _ = annulus(R, r, level)
        for sk, rk in kd:
            pin += sk*disc_lens(rad, (c[0], c[1]), rk, origin)[0]
    return assemble(kspec, pspec, k_first, K, V, (Mu, Mv_, Mw), kin, pin)


def assemble(kspec, pspec, k_first, K, V, mom, kin, pin):
    """Every operation from the common's volume and chart moments and each
    input's area inside the other."""
    VK, MK, AK = measures(kspec)
    VP, MP, AP = measures(pspec)
    Mc = tuple(K.o[i]*V+sum(mom[j]*(K.x, K.y, K.n)[j][i] for j in range(3)) for i in range(3))
    common = (V, Mc, kin+pin)
    fuse = (VK+VP-V, tuple(a+b-c for a, b, c in zip(MK, MP, Mc)), (AK-kin)+(AP-pin))
    if k_first:
        cut = (VK-V, tuple(a-c for a, c in zip(MK, Mc)), (AK-kin)+pin)
    else:
        cut = (VP-V, tuple(a-c for a, c in zip(MP, Mc)), (AP-pin)+kin)
    out = {}
    for op, (v, m, a) in (('common', common), ('fuse', fuse), ('cut', cut)):
        out[op] = None if v < mp.mpf(10)**-30 else (v, a, tuple(x/v for x in m))
    return out


def strip_form(case):
    """The torus against a rod across its axis: each slice the rectangle
    `[x0, x1] x [y0 - eta, y0 + eta]` (along and across the rod's axis)
    against the annulus."""
    kspec, pspec = (case.specs[0], case.specs[1]) if case.specs[0][0] == 'torus' else (case.specs[1],
                                                                                      case.specs[0])
    k_first = case.specs[0] is kspec
    K = Ideal(kspec[3])
    R, r = mp.mpf(kspec[1]), mp.mpf(kspec[2])
    _, boundaries, frame, lo_, hi_ = pspec
    fr = Ideal(frame)
    cx, cy, a = (mp.mpf(v) for v in boundaries[0].circle)
    e = K.vec_local(fr.n)
    assert abs(e[2]) < mp.mpf(10)**-30, 'a rod across the axis'
    ep = (-e[1], e[0])
    base = K.local(fr.world(cx, cy, Z))
    y0 = base[0]*ep[0]+base[1]*ep[1]
    xb = base[0]*e[0]+base[1]*e[1]
    x0, x1 = xb+mp.mpf(lo_), xb+mp.mpf(hi_)
    z0 = base[2]
    eta = lambda s: mp.sqrt(max(a*a-(s-z0)**2, Z))
    lo, hi = max(-r, z0-a), min(r, z0+a)
    pts = {lo, hi}
    q = lambda s: mp.sqrt(max(r*r-s*s, Z))
    rads = [lambda s: R+q(s), lambda s: R-q(s)]
    for rf in rads:
        for sg in (1, -1):
            for g in (lambda s, rf=rf, sg=sg: (y0+sg*eta(s))**2-rf(s)**2,
                      lambda s, rf=rf, sg=sg: x0**2+(y0+sg*eta(s))**2-rf(s)**2,
                      lambda s, rf=rf, sg=sg: x1**2+(y0+sg*eta(s))**2-rf(s)**2):
                pts |= set(roots_on(g, lo, hi, n=2000))
        for xc in (x0, x1):
            pts |= set(roots_on(lambda s, rf=rf, xc=xc: xc**2-rf(s)**2, lo, hi, n=2000))
    pts = sorted(p for p in pts if lo <= p <= hi)

    def chord_len(c, a0, a1, s):
        """Length of `{t in [a0, a1]: R - q < |(t, c)| < R + q}`."""
        kd, _ = annulus(R, r, s)
        tot = Z
        for sk, rk in kd:
            Q2 = rk*rk-c*c
            if Q2 <= 0:
                continue
            Q = mp.sqrt(Q2)
            tot += sk*max(Z, min(a1, Q)-max(a0, -Q))
        return tot

    def g(s):
        h = eta(s)
        kd, kw = annulus(R, r, s)
        A = mx = my = Z
        for sk, rk in kd:
            ar, m1, m2 = rect_disc(x0, x1, y0-h, y0+h, rk)
            A += sk*ar
            mx += sk*m1
            my += sk*m2
        pin = chord_len(x0, y0-h, y0+h, s)+chord_len(x1, y0-h, y0+h, s)
        mu, mv = mx*e[0]+my*ep[0], mx*e[1]+my*ep[1]
        return A, mu, mv, s*A, pin
    vals = [mp.quad(lambda s: g(s)[k], pts) for k in range(5)]
    V, Mu, Mv_, Mw, pin = vals
    # The torus's wall by its latitudes (`s = r sin phi`, element `r rho`),
    # the rod's by its angle (`s = z0 + a sin psi`, its line across at `y0 +
    # a cos psi`, element `a`).
    phis = {-mp.pi/2, mp.pi/2, 3*mp.pi/2}
    psis = {-mp.pi/2, mp.pi/2, 3*mp.pi/2}
    for s in pts:
        for c2, r2, out in ((Z, r, phis), (z0, a, psis)):
            x = (s-c2)/r2
            if -1 < x < 1:
                p0 = mp.asin(x)
                out |= {p0, mp.pi-p0}
    for rf in rads:
        for g2 in (lambda p, rf=rf: (y0+a*mp.cos(p))**2-rf(z0+a*mp.sin(p))**2,
                   lambda p, rf=rf: x0**2+(y0+a*mp.cos(p))**2-rf(z0+a*mp.sin(p))**2,
                   lambda p, rf=rf: x1**2+(y0+a*mp.cos(p))**2-rf(z0+a*mp.sin(p))**2):
            psis |= set(p for p in roots_on(g2, -mp.pi/2, 3*mp.pi/2, n=2000) if abs(z0+a*mp.sin(p)) < r)

    def kf(phi):
        s = r*mp.sin(phi)
        if not z0-a < s < z0+a:
            return Z
        h = eta(s)
        rho = R+r*mp.cos(phi)
        return r*rho*arc_in_rect(x0, x1, y0-h, y0+h, rho)

    def pf(psi):
        s = z0+a*mp.sin(psi)
        if not -r < s < r:
            return Z
        return a*chord_len(y0+a*mp.cos(psi), x0, x1, s)
    kin = mp.quad(kf, sorted(phis))
    pin += mp.quad(pf, sorted(psis))
    return assemble(kspec, pspec, k_first, K, V, (Mu, Mv_, Mw), kin, pin)


def apart_form(case):
    VA, MA, AA = measures(case.specs[0])
    VB, MB, AB = measures(case.specs[1])
    return {'common': None, 'fuse': (VA+VB, AA+AB, tuple((a+b)/(VA+VB) for a, b in zip(MA, MB))),
            'cut': (VA, AA, tuple(a/VA for a in MA))}


def pair_forms(case):
    if case.form == 'slices':
        return slices_form(case)
    if case.form == 'strip':
        return strip_form(case)
    if case.form == 'apart':
        return apart_form(case)
    return None


# ------------------------------------------------------------------ Monte Carlo

def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fa, fb = ref.FloatModel(pair.A), ref.FloatModel(pair.B)
    ba, bb = pair.A.box(), pair.B.box()
    box0 = [min(ba[0][i], bb[0][i]) for i in range(3)]
    box1 = [max(ba[1][i], bb[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        a, b = fa.contains(X), fb.contains(X)
        for op in ref.OPS:
            if sets[op](a, b):
                e = acc[op]
                e[0] += 1
                for i in range(3):
                    e[1][i] += X[i]
                    e[2][i] += X[i]*X[i]
    out = {}
    for op, (k, s1, s2) in acc.items():
        p = k/n
        vol = vbox*p
        sv = vbox*math.sqrt(max(p*(1-p), 1e-300)/n)
        if k > 1:
            c = [s/k for s in s1]
            sc = [math.sqrt(max(s2[i]/k-c[i]*c[i], 0)/k) for i in range(3)]
        else:
            c, sc = None, None
        out[op] = (vol, sv, c, sc)
    return out


# ------------------------------------------------------------------ the pairs' work

def dev3(a, b, size):
    """Largest relative deviation of (volume, world moments, area)."""
    (va, ma, aa), (vb, mb, ab) = a, b
    return max(abs(va-vb)/size**3, max(abs(x-y) for x, y in zip(ma, mb))/size**4, abs(aa-ab)/size**2)


def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, first, ops, mc_n = job
    obj, tool = first.obj, first.tool
    pair = ref.Pair(obj, tool)
    size = pair.size
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    checks = {}
    closed = {'A': pair.A.closed(), 'B': pair.B.closed()}
    checks['inputs'] = max(dev3(pair.input_measures(role, k), closed[role], size) for role in 'AB' for k in (0, 1))
    checks['second_family'] = max(dev3(pair.measures(op, 0), pair.measures(op, 1), size) for op in ref.OPS)
    (VA, MA, AA), (VB, MB, AB) = closed['A'], closed['B']
    dev = Z
    for k in (0, 1):
        Vc, Mc, Ac = pair.measures('common', k)
        Vf, Mf, Af = pair.measures('fuse', k)
        Vx, Mx, _ = pair.measures('cut', k)
        dev = max(dev, dev3((Vf, Mf, Z), (VA+VB-Vc, tuple(a+b-c for a, b, c in zip(MA, MB, Mc)), Z), size),
                  dev3((Vx, Mx, Z), (VA-Vc, tuple(a-c for a, c in zip(MA, Mc)), Z), size))
        checks['area_identity'] = max(checks.get('area_identity', Z), abs(Af+Ac-AA-AB)/size**2)
    checks['inclusion_exclusion'] = dev
    face = Z
    sw = pair.sweeps()
    for (k, role, fname), s in sw.items():
        if k == 0:
            other = sw[(1, role, fname)]
            face = max(face, abs(sum(s.result[c][0] for c in s.result)-sum(other.result[c][0] for c in other.result)))
    checks['face_classes'] = face/size**2
    mcr = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mcr[op]
        V, mm = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in mm]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = ref.near_coincidences(pair)
    counts = {k: pair.swept_count(k) for k in ('common', 'A-B', 'B-A', 'fuse')}
    overlap = pair.volume('common')[0] > mp.mpf(10)**-20*size**3
    rule = 1 if overlap else 2
    assert counts['fuse'] == rule, (name, 'the fuse swept', counts['fuse'], 'by the rule', rule)
    stats = {'missed': pair.missed(), 'quadrature': pair.quad_error()/size**4, 'margins': ref.margins(pair),
             'solids': counts, 'events': {f'{k}{r}:{f}': len(s.events) for (k, r, f), s in sw.items()}}
    return name, rows, res, checks, near, stats


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, pair_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15, 'closed_forms_degenerate': 1e-20,
          'degenerate_pairs': 1e-18,
          'inputs': 1e-30, 'second_family': 1e-30, 'inclusion_exclusion': 1e-30, 'area_identity': 1e-30,
          'face_classes': 1e-30, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-28}


def reference_checks(results, forms):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, Z), value)
        covered[key] = covered.get(key, 0)+1
    by_pair = {r[0]: r for r in results}
    firsts = {}
    for c in cases():
        firsts.setdefault(c.pair_name, c)
    for name, form in forms.items():
        if form is None:
            continue
        _, _, res, _, _, _ = by_pair[name]
        kind = 'exact_frames' if exact_pair(firsts[name]) else 'turned_frames'
        if any(c.kind == 'degenerate' for c in cases() if c.pair_name == name):
            kind = 'degenerate'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}', dev)
    degenerate = {c.pair_name for c in cases() if c.kind == 'degenerate'}
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            if name in degenerate and key != 'monte_carlo_sigma':
                # A tangency along a circle: the curves through it meet the
                # other surface in double roots, rounding's pieces there.
                note('degenerate_pairs', value)
            else:
                note(key, value)
        note('volume_quadrature_estimate', stats['quadrature'])
        if name not in degenerate:
            assert stats['missed'] == 0, (name, 'quadrature nodes off their interval\'s structure', stats['missed'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, torus_curved_boolean_reference.py: expect KIND {STEP}, reason TEXT for a '
           'degenerate case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, _, _, _ = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP}')
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


def jobs(mc_n, only=None):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items() if only is None or name in only]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert c.obj.torus is not None or c.tool.torus is not None, f'{c.name}: a torus'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', nargs='*', help='evaluate these pairs only (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    if args.only:
        results = run(jobs(args.samples, set(args.only)), evaluate, args.workers)
        firsts = {c.pair_name: c for c in listed}
        for name, rows, res, checks, near, stats in results:
            print(name, rows, {k: mp.nstr(v, 3) for k, v in checks.items()}, near,
                  {k: (mp.nstr(v, 3) if not isinstance(v, (dict, tuple, int)) else v) for k, v in stats.items()})
            form = pair_forms(firsts[name])
            if form:
                for op, want in form.items():
                    n, vol, area, centre = res[op]
                    if want is None:
                        print('  form', op, 'empty; reference', n)
                    else:
                        print('  form', op, mp.nstr(vol-want[0], 3), mp.nstr(area-want[1], 3),
                              [mp.nstr(centre[i]-want[2][i], 3) for i in range(3)])
        return
    results = run(jobs(args.samples), evaluate, args.workers)
    degenerate = {c.pair_name for c in listed if c.kind == 'degenerate'}
    near = [f'{name}: near coincidences {n}' for name, _, _, _, n, _ in results if n and name not in degenerate]
    if near:
        raise SystemExit('\n'.join(near))
    firsts = {}
    for c in listed:
        firsts.setdefault(c.pair_name, c)
    forms = dict(run(list(firsts.values()), closed_job, args.workers))
    worst, covered = reference_checks(results, forms)
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
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames')
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('closed forms:', sum(1 for f in forms.values() if f is not None), 'of', len(forms), 'pairs')
    by_pair = {r[0]: r for r in results}
    least = {}
    for name, _, _, _, n, stats in results:
        if n:
            print('declared near coincidences:', name, '; '.join(n[:4]))
        if name not in degenerate:
            for k, v in stats['margins'].items():
                if v < least.get(k, (mp.inf, None))[0]:
                    least[k] = (v, name)
    print('least margins (non-degenerate pairs):', ', '.join(f'{k} {mp.nstr(v, 3)} ({n})'
                                                             for k, (v, n) in sorted(least.items())))
    for name, _, _, _, _, stats in results:
        print('solids', name, stats['solids'], 'events', sum(stats['events'].values()))
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
