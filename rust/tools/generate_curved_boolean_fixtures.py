#!/usr/bin/env python3
"""Fixtures for S9c.1 of REVIEW_NOTES.md: Booleans of prisms whose profiles
hold arcs and circles (cylindrical walls) as well as lines, in any relative
position, every pair of faces meeting in lines, circles or ellipses (plane
and plane, plane and cylinder, two cylinders parallel, coaxial or of equal
radii with crossing axes).

`boolean-curved-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); the two frames' axes always
differ (S9b's and S9c's domain) and at least one profile holds an arc or a
circle. Only frames whose stored axes the kernel's `Frame3::new` gives bit
for bit as `stored_axes` are used (`boolean-curved-frames.tsv` records
them, as `split-frames.tsv`): exact ones (`XY`, `SIDE`, `DOWN`, `TURN`: the
world's axes permuted or reversed) and turned ones (`TILT`, `TILT2` and
`TILTX`, normals `(0, 3, 4)` and `(0, -4, 3)`; `LEAN`, normal `(3, 0, 4)`;
`R125`, the `XY` normal with `x` along `(12, 5, 0)`). The split fixtures'
`ROT` (`x` along `(3, 4, 0)`) is not used: on this platform the kernel's
`x` differs from `stored_axes` in its last bit (the platform `hypot`).

`boolean-curved-expected.tsv` gives per case, from `curved_boolean_reference.py`:

* `expect KIND S9c.1`: the declared outcome (`solid`: one or more solids;
  `empty`: none; `degenerate`: a tangency of a plane and a cylinder along a
  generatrix, of two cylinders, or a result touching itself, `Degenerate`
  in the decisions), then for a degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

Before writing, `reference_checks` compares the reference with
independent results (limits relative to the case's size): closed forms
(Steinmetz solids of perpendicular and oblique axes, `16 r^3 / (3 sin phi)`
and `16 r^2 / sin phi`; an oblique cylinder between two planes by
Cavalieri; a quarter cylinder; a box corner in a cylinder; coaxial and
parallel cylinders, lenses; a cylinder cut by an oblique plane) within
1e-30 in exact frames and 1e-15 in turned ones (whose stored axes are not
exactly orthonormal); `fuse = A + B - common` and `cut = A - common`
(volumes and moments, `A` and `B` in closed form, each operation's slices
computed apart), `area(fuse) + area(common) = area(A) + area(B)` (pairs
without faces shared with opposite orientations), every face's classes
summing to its closed-form area, both sides' shared areas equal, and a
second slicing direction for parallel axes, within 1e-30; Monte-Carlo
estimates (200,000 uniform points per pair) of every volume and centre
within 5 standard errors; S9a's and S9b's own references on their 45 and 45
fixtures (same frames; polyhedra in any position): equal solid counts and
the 25 printed digits within 1e-15 (S9a's frame coordinates take the stored
axes as orthonormal) and 1e-24. A scan for near coincidences (a class of a
face or a result's volume positive but below 1e-9 of the case's size,
slicing breakpoints closer than that) must find none: two fixtures were
moved off such coincidences of the slicing (two events at one height).
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
import curved_boolean_reference as ref

ROOT = Path(__file__).resolve().parents[1]
BOOLEAN_OPERATION = 93

FRAMES = {
    'XY': (0.0, 0.0, 1.0, 1.0, 0.0, 0.0),
    'SIDE': (1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
    'DOWN': (0.0, 0.0, -1.0, 1.0, 0.0, 0.0),
    'TURN': (0.0, 0.0, 1.0, 0.0, 1.0, 0.0),
    'TILT': (0.0, 3.0, 4.0, 1.0, 0.0, 0.0),
    'TILT2': (0.0, 3.0, 4.0, 0.0, -4.0, 3.0),
    'TILTX': (0.0, -4.0, 3.0, 1.0, 0.0, 0.0),
    'LEAN': (3.0, 0.0, 4.0, 0.0, 1.0, 0.0),
    'R125': (0.0, 0.0, 1.0, 12.0, 5.0, 0.0),
}


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def disc(cx, cy, r):
    return Boundary(circle=(cx, cy, r))


def stadium(half, r):
    """A stadium of straight length `2 half` and radius `r` about the origin."""
    return Boundary(points=[(-half, -r), (half, -r), (half, r), (-half, r)],
                    segments=[None, (half, 0.0, r, True), None, (-half, 0.0, r, True)])


def prism(boundaries, frame, start, end, op):
    return Case('', 1e-7, op, frame, start, end, boundaries)


class Boolean:
    def __init__(self, name, operation, obj, tool, kind='solid', reason=None, opposite=False, near=None):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.opposite, self.near = opposite, near
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = prism(*obj, op=91)
        self.obj.name = name
        self.tool = prism(*tool, op=92)
        self.tool.name = name
        self.frames = (obj[1], tool[1])

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, **kw):
    """Cases of one pair: `ops` maps each operation to its declared kind
    (or `(kind, reason)`)."""
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, kind, reason, **kw))
    return out


SOLID3 = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
TANGENT_PLANE = 'a plane tangent to a cylinder along a generatrix'
TANGENT_CYLINDERS = 'two cylinders tangent along a generatrix'
STEINMETZ_TOUCH = 'the result touches itself at the two points where the equal cylinders are tangent'
# S9c.2b's decisions (REVIEW_NOTES.md): equal cylinders whose axes meet in
# stored turned frames have exactly equal extents across the common
# perpendicular, a double tangency, refused.
STORED_NODE = ('equal cylinders with meeting axes in turned stored frames: their models touch '
               'at the two ends of their common extent (S9c.2b)')


def cases():
    """Every class of the S9c.1 decisions, in exact and turned frames: a
    tilted cylinder through a box (ellipse edges), a box corner in a
    cylinder (exact and leaning), Steinmetz solids of perpendicular and
    oblique equal cylinders (exact, tilted and both tilted), parallel
    cylinders sharing their caps' planes, coaxial cylinders, a tilted pin
    through a box's hole without touching, a coaxial pin filling a round
    hole (coincident cylinders of opposite orientations), a tilted cylinder
    across a hole's wall, a stadium in a turned and in a tilted frame, a
    quarter cylinder (planes through the axis), and the degenerate classes:
    a plane tangent along a generatrix, cylinders tangent outside and
    inside, a Steinmetz cut touching itself at two points."""
    box = [square(0.0, 0.0, 10.0, 10.0)]
    out = []
    out += group('tilted_cylinder_box', (box, at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 1.5)], at('TILT', (5, 3, 0)), -2.0, 8.0), SOLID3)
    out += group('corner_in_cylinder', ([disc(0.0, 0.0, 4.0)], at('XY', (0, 0, 0)), 0.0, 6.0),
                 ([square(1.5, 2.0, 8.0, 10.0)], at('SIDE', (1, 0, 0)), 0.0, 8.0), SOLID3)
    out += group('corner_in_cylinder_leaning', ([disc(0.0, 0.0, 4.0)], at('XY', (0, 0, 0)), 0.0, 6.0),
                 ([square(0.0, 0.0, 6.0, 6.0)], at('LEAN', (2.5, 1, 3)), 0.0, 6.0), SOLID3)
    out += group('steinmetz', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -5.0, 5.0),
                 ([disc(0.0, 0.0, 2.0)], at('SIDE', (-5, 0, 0)), 0.0, 10.0),
                 {'fuse': 'solid', 'cut': ('degenerate', STEINMETZ_TOUCH), 'common': 'solid'})
    out += group('steinmetz_oblique', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -7.0, 7.0),
                 ([disc(0.0, 0.0, 2.0)], at('TILT', (0, 0, 0)), -9.0, 9.0),
                 {'fuse': ('degenerate', STORED_NODE), 'common': ('degenerate', STORED_NODE)})
    out += group('steinmetz_tilted', ([disc(0.0, 0.0, 1.5)], at('TILT', (0, 0, 0)), -6.0, 6.0),
                 ([disc(0.0, 0.0, 1.5)], at('TILTX', (0, 0, 0)), -6.0, 6.0),
                 {'common': ('degenerate', STORED_NODE)})
    out += group('parallel_cylinders', ([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 2.0)], at('R125', (4, 0.5, 0)), 0.0, 5.0), SOLID3)
    out += group('coaxial', ([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 1.5)], at('R125', (0, 0, -1)), 0.0, 7.0), SOLID3)
    out += group('pin_through_hole', ([square(0.0, 0.0, 10.0, 10.0), square(3.0, 3.0, 7.0, 7.0)],
                                      at('XY', (0, 0, 0)), 0.0, 2.0),
                 ([disc(0.0, 0.0, 0.8)], at('TILT', (5, 5, 1)), -6.0, 6.0),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'empty'})
    out += group('pin_fills_hole', ([square(0.0, 0.0, 10.0, 10.0), disc(5.0, 5.0, 2.0)],
                                    at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(5.0, -5.0, 2.0)], at('DOWN', (0, 0, 7)), 0.0, 9.0),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'empty'}, opposite=True)
    out += group('across_hole', ([square(0.0, 0.0, 10.0, 10.0), square(3.0, 3.0, 7.0, 7.0)],
                                 at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(0.0, 0.0, 1.0)], at('TILT', (5, 7, 2)), -5.0, 5.0), SOLID3)
    out += group('stadium_turned', ([square(0.0, 0.0, 10.0, 6.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([stadium(3.0, 1.0)], at('R125', (8, 3, -1)), 0.0, 3.0), SOLID3)
    out += group('stadium_tilted', (box, at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([stadium(3.0, 1.0)], at('TILT2', (5, 5, 2.5)), -1.0, 1.0), SOLID3)
    out += group('quarter', ([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(0.0, -1.0, 5.0, 6.0)], at('SIDE', (0, 0, 0)), 0.0, 5.0), {'common': 'solid'})
    out += group('quarter_turned', ([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(0.0, 0.0, 5.0, 5.0)], at('R125', (0, 0, -1)), 0.0, 7.0), {'common': 'solid'})
    out += group('tangent_plane', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(2.0, 1.0, 6.0, 4.0)], at('SIDE', (-3, 0, 0)), 0.0, 6.0),
                 {op: ('degenerate', TANGENT_PLANE) for op in ('fuse', 'cut', 'common')})
    out += group('tangent_cylinders', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 1.0)], at('TURN', (3, 0, -1)), 0.0, 7.0),
                 {op: ('degenerate', TANGENT_CYLINDERS) for op in ('fuse', 'common')})
    out += group('tangent_inside', ([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(2.0, 0.0, 1.0)], at('DOWN', (0, 0, 6)), 0.0, 7.0),
                 {'cut': ('degenerate', TANGENT_CYLINDERS)})
    return out


# ------------------------------------------------------------------ closed forms

def corner_region(a, b, r):
    """Area, first moments and boundary lengths of `{x >= a, y >= b, x^2 +
    y^2 <= r^2}`."""
    a, b, r = mp.mpf(a), mp.mpf(b), mp.mpf(r)
    x1, y1 = mp.sqrt(r*r-b*b), mp.sqrt(r*r-a*a)
    Fx = lambda x: (x*mp.sqrt(r*r-x*x)+r*r*mp.asin(x/r))/2
    area = Fx(x1)-Fx(a)-b*(x1-a)
    mx = (-(r*r-x1*x1)**mp.mpf(1.5)+(r*r-a*a)**mp.mpf(1.5))/3-b*(x1*x1-a*a)/2
    my = (((r*r-b*b)*x1-x1**3/3)-((r*r-b*b)*a-a**3/3))/2
    arc = r*(mp.atan2(y1, a)-mp.atan2(b, x1))
    return area, mx, my, (y1-b)+(x1-a)+arc


def lens(r1, r2, d):
    """Area, centroid distance from the first centre and perimeter of two
    circles' lens."""
    r1, r2, d = mp.mpf(r1), mp.mpf(r2), mp.mpf(d)
    xc = (d*d+r1*r1-r2*r2)/(2*d)
    a1, a2 = mp.acos(xc/r1), mp.acos((d-xc)/r2)
    s1, s2 = r1*r1*(a1-mp.sin(a1)*mp.cos(a1)), r2*r2*(a2-mp.sin(a2)*mp.cos(a2))
    g1 = 2*r1*mp.sin(a1)**3/(3*(a1-mp.sin(a1)*mp.cos(a1)))
    g2 = 2*r2*mp.sin(a2)**3/(3*(a2-mp.sin(a2)*mp.cos(a2)))
    return s1+s2, (s1*g1+s2*(d-g2))/(s1+s2), 2*r1*a1+2*r2*a2


def closed_forms():
    """{(pair, operation): (volume, area, centre)} in exact arithmetic
    (mpmath); turned frames' values are the ideal orthonormal ones."""
    pi = mp.pi
    out = {}
    # Steinmetz, perpendicular: r = 2, lengths 10.
    r = mp.mpf(2)
    vc, ac = 16*r**3/3, 16*r**2
    va = vb = pi*r*r*10
    aa = ab = 2*pi*r*r+2*pi*r*10
    out[('steinmetz', 'common')] = (vc, ac, (0, 0, 0))
    out[('steinmetz', 'fuse')] = (va+vb-vc, aa+ab-ac, (0, 0, 0))
    # Oblique: sin(phi) = 0.6, lengths 14 and 18.
    sphi = mp.mpf(3)/5
    vc, ac = 16*r**3/(3*sphi), 16*r**2/sphi
    va, vb = pi*r*r*14, pi*r*r*18
    aa, ab = 2*pi*r*r+2*pi*r*14, 2*pi*r*r+2*pi*r*18
    out[('steinmetz_oblique', 'common')] = (vc, ac, (0, 0, 0))
    out[('steinmetz_oblique', 'fuse')] = (va+vb-vc, aa+ab-ac, (0, 0, 0))
    r = mp.mpf(1.5)
    out[('steinmetz_tilted', 'common')] = (16*r**3/3, 16*r**2, (0, 0, 0))
    # An oblique cylinder between two planes (Cavalieri): r 1.5, height 5, cos 0.8.
    r, h, c = mp.mpf(1.5), mp.mpf(5), mp.mpf(4)/5
    out[('tilted_cylinder_box', 'common')] = (pi*r*r*h/c, 2*pi*r*r/c+2*pi*r*h/c, (5, 3+mp.mpf(2.5)*3/4, 2.5))
    # A box corner in a cylinder: {x >= 1, y >= 1.5, x^2 + y^2 <= 16}, z in [2, 6].
    area, mx, my, length = corner_region(1, 1.5, 4)
    out[('corner_in_cylinder', 'common')] = (4*area, 2*area+4*length, (mx/area, my/area, 4))
    # A quarter cylinder, r 3, height 5 (the second turned by atan2(5, 12)).
    r, h = mp.mpf(3), mp.mpf(5)
    g = 4*r/(3*pi)
    q = (pi*r*r*h/4, pi*r*r/2+pi*r*h/2+2*r*h)
    out[('quarter', 'common')] = q+((g, g, 2.5),)
    phi = mp.atan2(5, 12)+pi/4
    out[('quarter_turned', 'common')] = q+((g*mp.sqrt(2)*mp.cos(phi), g*mp.sqrt(2)*mp.sin(phi), 2.5),)
    # Coaxial: r 3 over [0, 5], r 1.5 over [-1, 6].
    out[('coaxial', 'common')] = (pi*mp.mpf(2.25)*5, 2*pi*mp.mpf(2.25)+2*pi*mp.mpf(1.5)*5, (0, 0, 2.5))
    out[('coaxial', 'cut')] = (pi*mp.mpf(6.75)*5, 2*pi*mp.mpf(6.75)+2*pi*3*5+2*pi*mp.mpf(1.5)*5, (0, 0, 2.5))
    out[('coaxial', 'fuse')] = (pi*9*5+pi*mp.mpf(2.25)*2, 54*pi, (0, 0, 2.5))
    # Parallel cylinders: r 3 at the origin, r 2 at (4, 0.5), both over [0, 5].
    d = mp.sqrt(mp.mpf('16.25'))
    la, lg, lp = lens(3, 2, d)
    u = (4/d, mp.mpf(0.5)/d)
    out[('parallel_cylinders', 'common')] = (5*la, 2*la+5*lp, (lg*u[0], lg*u[1], 2.5))
    return out


def plane_cut_check():
    """A cylinder (r 2, z in [0, 6]) above the plane `z = 3 - 0.75 y` (a
    large tilted box): volume 12 pi, area 21 pi, centre (0, 1/4, 4.40625)."""
    A = prism([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), 0.0, 6.0, 91)
    B = prism([square(-20.0, -20.0, 20.0, 20.0)], at('TILT', (0, 0, 3)), 0.0, 20.0, 92)
    return ('sliced_cylinder', A, B, 'common', (12*mp.pi, 21*mp.pi, (0, mp.mpf(1)/4, mp.mpf('4.40625'))))


# ------------------------------------------------------------------ Monte Carlo

class FloatPrism:
    def __init__(self, p):
        self.o = [float(c) for c in p.o]
        self.inv = [[float(c) for c in r] for r in p.inv]
        self.lo, self.hi = float(p.lo), float(p.hi)
        self.elements = []
        for el in p.profile.elements:
            if el.kind == 'seg':
                self.elements.append(('seg', float(el.p[0]), float(el.p[1]), float(el.e[0]), float(el.e[1])))
            else:
                self.elements.append(('round', float(el.c[0]), float(el.c[1]), float(el.r), float(el.start),
                                      float(el.sweep), el.kind == 'circle'))

    def contains(self, X):
        d = (X[0]-self.o[0], X[1]-self.o[1], X[2]-self.o[2])
        w = self.inv[2][0]*d[0]+self.inv[2][1]*d[1]+self.inv[2][2]*d[2]
        if not self.lo < w < self.hi:
            return False
        u = self.inv[0][0]*d[0]+self.inv[0][1]*d[1]+self.inv[0][2]*d[2]
        v = self.inv[1][0]*d[0]+self.inv[1][1]*d[1]+self.inv[1][2]*d[2]
        ru, rv = 0.8150192, 0.5794340
        count = 0
        for el in self.elements:
            if el[0] == 'seg':
                _, px, py, ex, ey = el
                den = ru*ey-rv*ex
                if den == 0:
                    continue
                wx, wy = px-u, py-v
                t = (wx*ey-wy*ex)/den
                s = (wx*rv-wy*ru)/den
                if t > 0 and 0 <= s < 1:
                    count += 1
            else:
                _, cx, cy, r, start, sweep, full = el
                Wx, Wy = u-cx, v-cy
                b = ru*Wx+rv*Wy
                disc_ = b*b-(Wx*Wx+Wy*Wy-r*r)
                if disc_ <= 0:
                    continue
                sq = math.sqrt(disc_)
                for t in (-b-sq, -b+sq):
                    if t > 0:
                        if full:
                            count += 1
                        else:
                            phi = math.atan2(v+t*rv-cy, u+t*ru-cx)
                            if (phi-start) % (2*math.pi) <= sweep:
                                count += 1
        return count % 2 == 1


def monte_carlo(obj, tool, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    A, B = ref.Prism(obj), ref.Prism(tool)
    fa, fb = FloatPrism(A), FloatPrism(B)
    lo, hi = [], []
    for p in (A, B):
        b0, b1 = p.bounds()
        lo.append([float(c) for c in b0])
        hi.append([float(c) for c in b1])
    box0 = [min(lo[0][i], lo[1][i]) for i in range(3)]
    box1 = [max(hi[0][i], hi[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        a, b = fa.contains(X), fb.contains(X)
        for op in ref.OPS:
            if ref.SET[op](a, b):
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

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, obj, tool, ops, opposite, mc_n = job
    pair = ref.Pair(obj, tool)
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    A, B = pair.A, pair.B
    size = pair.size
    va, ma, aa = A.measures()
    vb, mb, ab = B.measures()
    vols = pair.volumes()
    checks = {}
    # Inclusion and exclusion (operations computed apart).
    dev = mp.mpf(0)
    vf, mf = vols['fuse']
    vc, mc = vols['common']
    vt, mt = vols['cut']
    for x, y in ((vf, va+vb-vc), (vt, va-vc)):
        dev = max(dev, abs(x-y)/size**3)
    for i in range(3):
        dev = max(dev, abs(mf[i]-(ma[i]+mb[i]-mc[i]))/size**4, abs(mt[i]-(ma[i]-mc[i]))/size**4)
    checks['inclusion_exclusion'] = dev
    # Faces: classes against closed-form areas; the area identity.
    fdev = mp.mpf(0)
    for tag, f, cls, _ in pair.face_areas():
        fdev = max(fdev, abs(sum(cls.values())-f.closed_area())/size**2)
    checks['face_classes'] = fdev
    same_a = sum((cls['same'] for tag, f, cls, _ in pair.face_areas() if tag == 'A'), mp.mpf(0))
    same_b = sum((cls['same'] for tag, f, cls, _ in pair.face_areas() if tag == 'B'), mp.mpf(0))
    opp_a = sum((cls['opp'] for tag, f, cls, _ in pair.face_areas() if tag == 'A'), mp.mpf(0))
    opp_b = sum((cls['opp'] for tag, f, cls, _ in pair.face_areas() if tag == 'B'), mp.mpf(0))
    checks['shared_both_sides'] = max(abs(same_a-same_b), abs(opp_a-opp_b))/size**2
    if not opposite:
        checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')-aa-ab)/size**2
    # A second slicing direction (parallel axes).
    if pair.slicing.parallel:
        axes = [(ref.F(1), ref.F(0), ref.F(0)), (ref.F(0), ref.F(1), ref.F(0)), (ref.F(0), ref.F(0), ref.F(1))]
        used = min(axes, key=lambda e: abs(ref.dot(A.n, e)))
        other = [e for e in axes if e != used and ref.cross(A.n, e) != (0, 0, 0)][0]
        second = ref.Slicing(A, B, other).measure()
        dev = mp.mpf(0)
        for op in ref.OPS:
            dev = max(dev, abs(second[op][0]-vols[op][0])/size**3)
            dev = max(dev, max(abs(x-y) for x, y in zip(second[op][1], vols[op][1]))/size**4)
        checks['second_axis'] = dev
    # Monte Carlo.
    mc = monte_carlo(obj, tool, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mc[op]
        v = float(vols[op][0])
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/vols[op][0]) for x in vols[op][1]]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(pair)
    stats = {'volume_quadrature': pair.slicing.quad_error/size**4,
             'face_quadrature': max(sw.quad_error for _, _, _, sw in pair.face_areas())/size**2,
             'breaks': len(pair.slicing.breaks)}
    return name, rows, {op: res[op][0] for op in ref.OPS}, res, checks, near, stats


def near_coincidences(pair):
    """Signs of accidental near coincidences: a face class, or a result's
    volume, positive but thinner than 1e-9 of the case's size; slicing
    breakpoints closer than that (and farther than 1e-30)."""
    size = pair.size
    out = []
    for tag, f, cls, _ in pair.face_areas():
        for c, v in cls.items():
            if 0 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{tag} {ref.describe(f)} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volumes()[op][0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    b = pair.slicing.breaks
    for p, q in zip(b, b[1:]):
        if q-p < mp.mpf(10)**-9*size:
            out.append(f'breakpoints {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def check_job(job):
    """A closed-form-only pair (the sliced cylinder)."""
    name, obj, tool, op, want = job
    pair = ref.Pair(obj, tool)
    return name, op, pair.result(op), want


def s9a_job(case):
    rows, _ = ref.rows(case.obj, case.operation, case.tool)
    return case.name, rows[0]


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def expected_file(path):
    out = {}
    for line in path.read_text().splitlines()[1:]:
        name, row = line.split('\t')
        if row.split()[0] in ('result', 'empty'):
            out[name] = row
    return out


def reference_checks(results, workers):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    by_pair = {r[0]: r for r in results}
    forms = closed_forms()
    exact_frames = {'XY', 'SIDE', 'DOWN', 'TURN'}
    frames_of = {c.pair_name: c.frames for c in cases()}
    for (name, op), (V, S, C) in forms.items():
        _, _, _, res, _, _, _ = by_pair[name]
        n, vol, area, centre = res[op]
        size = max(abs(V)**(mp.mpf(1)/3), 1)
        dev = max(abs(vol-V)/abs(V), abs(area-S)/abs(S), max(abs(centre[i]-C[i]) for i in range(3))/size)
        kind = 'exact' if all(frame_name(f) in exact_frames for f in frames_of[name]) else 'turned'
        note(f'closed_forms_{kind}_frames', dev)
    name, op, (n, vol, area, centre), (V, S, C) = run([plane_cut_check()], check_job, 1)[0]
    note('closed_forms_turned_frames', max(abs(vol-V)/V, abs(area-S)/S, max(abs(centre[i]-C[i]) for i in range(3))))
    for name, rows, _, res, checks, near, stats in results:
        for key, value in checks.items():
            if key == 'monte_carlo_sigma':
                worst[key] = max(worst.get(key, mp.mpf(0)), value)
                covered[key] = covered.get(key, 0)+1
            else:
                note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    # S9a's and S9b's references on their fixtures.
    import generate_boolean_fixtures as s9a
    import generate_polyhedral_fixtures as s9b
    for label, module, path in (('s9a_reference', s9a, 'boolean-expected.tsv'),
                                ('s9b_reference', s9b, 'boolean-polyhedra-expected.tsv')):
        want = expected_file(ROOT/'fixtures'/path)
        for case_name, row in run(module.cases(), s9a_job, workers):
            a, b = row.split(), want[case_name].split()
            assert a[0] == b[0] and a[1:2] == b[1:2], (label, case_name, a, b)
            dev = mp.mpf(0)
            for x, y in zip(a[2:], b[2:]):
                dev = max(dev, abs(mp.mpf(x)-mp.mpf(y))/max(1, abs(mp.mpf(y))))
            note(label, dev)
    return worst, covered


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = ['# case\trow (S9c.1, curved_boolean_reference.py: expect KIND S9c.1, reason TEXT for a degenerate '
           'case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, res, _, near, _ = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9c.1')
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
    return {'boolean-curved-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-curved-expected.tsv': '\n'.join(out)+'\n',
            'boolean-curved-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c.obj, c.tool, [], c.opposite])[2].append(c.operation)
    return [(name, obj, tool, ops, opposite, mc_n) for name, (obj, tool, ops, opposite) in pairs.items()]


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
        assert stored_axes(c.obj.frame)[1:] != stored_axes(c.tool.frame)[1:], f'{c.name}: frames with equal axes'
        assert any(b.circle is not None or b.segments is not None for b in c.obj.boundaries+c.tool.boundaries), \
            f'{c.name}: no arc'
    results = run(jobs(args.samples), evaluate, args.workers)
    by_pair = {r[0]: r for r in results}
    declared = {c.pair_name: c.near for c in listed}
    for name, _, _, _, _, near, _ in results:
        if near and not declared[name]:
            raise SystemExit(f'{name}: near coincidences {near}')
    worst, covered = reference_checks(results, args.workers)
    limits = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
              'inclusion_exclusion': 1e-30, 'face_classes': 1e-30, 'shared_both_sides': 1e-30,
              'area_identity': 1e-30, 'second_axis': 1e-30, 'monte_carlo_sigma': 5,
              's9a_reference': 1e-15, 's9b_reference': 1e-24, 'volume_quadrature_estimate': 1e-30,
              'face_quadrature_estimate': 1e-30}
    for key, value in worst.items():
        assert value <= limits[key], (key, mp.nstr(value, 3))
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    for name, _, _, _, _, near, _ in results:
        if near:
            print('declared near coincidences:', name, '; '.join(near[:4]))


if __name__ == '__main__':
    main()
