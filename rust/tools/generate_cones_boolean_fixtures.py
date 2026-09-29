#!/usr/bin/env python3
"""Fixtures for S9d.3b of REVIEW_NOTES.md: Booleans of a cone or frustum
(`Solid::cone_with`) against a prism with arcs and circles (cylindrical
walls), a sphere or cap (`Solid::sphere_with`), and another cone or frustum,
before any of its kernel code.

`boolean-cones-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a cone's block is its `frame` and
`cone BOTTOM TOP HEIGHT` rows, a sphere's its `frame` and `sphere R LOW
HIGH`, either input or both); `boolean-cones-expected.tsv` gives per case
`expect KIND STEP` (`solid`, `empty` or `degenerate`, then `reason TEXT`)
and `result N volume area cx cy cz` or `empty` from
`cones_boolean_reference.py`, as the cone fixtures; `boolean-cones-
frames.tsv` the stored axes' bits. Frames are the curved generator's
(`FRAMES`, their stored axes the kernel's bit for bit). A cap or zone only
in an exact frame centred at the world's origin, its axis a cone's (the
reference slices along it).

`STEP` is the sub-step that would take the case if S9d.3b were split as
S9d.2 was: `S9d.3b.1` where some input's rulings (a cone's, or a
cylinder's generatrices) meet every quadric of the other transversally
wherever their curve runs, so every curve is a circle or a graph over that
input's angle (`Curve3::Meet` with a cone carrier, S9d.2a's rings), and
`S9d.3b.2` where both inputs' rulings are tangent to the other's quadric on
their faces (loops with turning points on both, S9d.2b's graphs over the
height), decided by the exact discriminants of the rulings against the
other's quadric (`ruling_tangencies`).

Pairs, with `F` the frustum of radii 2 and 1 and height 2 and `A` the cone
of radius 1.5 and height 3 (its apex on top), both on `XY` at the origin
unless said. A cone and a cylinder: a coaxial pipe through `F` (`pipe_coax`,
circles, all three operations; the pipe less `F` two solids, `pipe_ends`),
a vertical cylinder through `F`'s wall and base off the axis (`bite`: an arc
of the quartic over the cylinder's angle, ending on the base rim), a rod
across `F` through its wall (`rod`: two rings; the rod less `F` its two
ends, `rod_ends`), a rod across `A`'s tip (`tip_rod`: two rings over the
cone's angle, `A` less the rod its tip and base), `A`'s apex inside a rod
(`apex_rod`, one ring), a box with a coaxial round hole whose wall cuts
`F` (`hole_box`, circles), a stadium through `F` (`stadium`: hyperbolas
and quartics), a coaxial pipe with `F` in `TILT` (`pipe_tilt`). A cone and
a sphere: coaxial through `F`'s top disc (`ball_coax`, circles, all three),
on `A`'s apex region (`ball_apex`), off the axis through `F`'s wall
(`ball_side`, a loop), a dome (a hemisphere) through which a cone passes
(`dome_cone`, the cap's disc cut in a circle), a sphere against `F` in
`TILT` (`ball_tilt`). Two cones: coaxial tips interpenetrating
(`cones_apex`, a circle), coaxial frusta of opposite slopes (`frusta_cross`),
a frustum standing on `F`'s top disc (`frusta_stack`: coplanar discs of
opposite orientation, the fuse one solid, the common empty), a cone across
`F` through its wall (`cones_cross`, two rings; the crossing cone less `F`
its two ends, `cones_ends`), `A` in `LEAN` against `F` (`cones_lean`), two
cones of parallel axes (`cones_parallel`, a loop). And `degenerate` with
reasons: a cylinder tangent to a cone along a ruling (`rod_ruling`), a
sphere inscribed in a cone, tangent along a circle (`ball_inscribed`), two
cones apex to apex (`apex_to_apex`), two equal coaxial cones
(`cones_equal`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms by one quadrature (mpmath's tanh-sinh between the kinks)
  along the cone's axis in its ideal frame: an input whose sections are
  discs about a parallel axis (a cylinder, a sphere or cap, a cone; coaxial
  or not) by the lens of two discs (areas, first moments by circular
  segments, each circle's angle inside the other: the cone's wall by
  `sqrt(1 + k^2) r` times its angle, a cylinder's by its radius, a
  sphere's by Archimedes, a cone's by its own element; end discs by the
  lens, coplanar ones apart); a cylinder across the axis by the strip `|y -
  y0| <= sqrt(a^2 - (z - z0)^2)` against the disc (two circular segments,
  the wall's angle inside, the cylinder's wall by its generatrices' chords
  inside the disc over its angle);
  a box with a coaxial hole as S9d.3a's box less the coaxial core; equal
  cones as the cone; within 1e-30 in exact frames and 1e-15 in turned ones
  (their stored axes not exactly orthonormal);
* `fuse = A + B - common` and `cut = A - common` (each operation sliced
  apart), both inputs' slicings against their closed forms, the area
  identity `area(fuse) + area(common) + 2 opp = area(A) + area(B)`, every
  face's classes (swept) summing to its area, both inputs' areas, both
  sides' shared areas equal, the cone's wall inside the other from the
  slices where they are normal to its axis, a second slicing direction
  (another valid chart: pairs without a cap, coplanar faces or a
  declared degeneracy), within 1e-30;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, slicing breakpoints closer
  than that), and every edge of one input and every vertex at least 1e-3
  from tangency with, or incidence on, each surface of the other (the
  critical values of the surface's equation along the edge, relative to
  the surface's radius squared or the case's size; an edge in a coplanar
  face of the other apart): none but in the declared degenerate pairs.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import itertools
import math
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, encode_boolean_case
from curve_surface_reference import stored_axes
import cones_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, stadium
from generate_sphere_boolean_fixtures import frame_name, ideal_axes
from generate_cone_boolean_fixtures import FloatCone, box_form, segment
from generate_spheres_boolean_fixtures import FloatArcPrism, FloatBall, lens

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-cones'
STEP = 'S9d.3b'
BOOLEAN_OPERATION = 93
HP = math.pi/2
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
Z = mp.mpf(0)


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def circle(cx, cy, r):
    return Boundary(circle=(float(cx), float(cy), float(r)))


def make(spec, op):
    if spec[0] == 'cone':
        _, r0, r1, h, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))
    if spec[0] == 'sphere':
        _, r, frame, low, high = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), low, high))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def cone(r0, r1, h, frame):
    return ('cone', r0, r1, h, frame)


def sphere(r, frame, low=-HP, high=HP):
    return ('sphere', r, frame, low, high)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def spec_frame(spec):
    return spec[4] if spec[0] == 'cone' else spec[2]


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
        self.frames = (spec_frame(obj), spec_frame(tool))

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


RULING = 'a cylinder tangent to the cone along a ruling'
INSCRIBED = 'a sphere tangent to the cone along a circle'
# Declared after S9d.3b.2's kernel: the sphere touches the frustum's base
# plane at a point inside its disc.
TOUCHING = 'a sphere tangent to the base plane inside its disc'
APEXES = 'two cones touching at their apexes'
EQUAL = 'two equal cones'
ALL = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}


def cases():
    """S9d.3b's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    F_ = cone(2, 1, 2, at('XY', O))
    A_ = cone(1.5, 0, 3, at('XY', O))
    out = []
    pipe = prism([circle(0, 0, 1.5)], at('XY', O), -1, 3)
    out += group('pipe_coax', F_, pipe, ALL, ('axial',))
    out += group('pipe_ends', pipe, F_, {'cut': 'solid'}, ('axial',))
    out += group('bite', F_, prism([circle(2, 0.5, 0.75)], at('XY', O), -1, 3), {'cut': 'solid', 'common': 'solid'},
                 ('axial',))
    rod = prism([circle(0.25, 1, 0.5)], at('SIDE', O), -4, 4)
    out += group('rod', F_, rod, {'cut': 'solid', 'common': 'solid'}, ('strip',))
    out += group('rod_ends', rod, F_, {'cut': 'solid'}, ('strip',))
    out += group('tip_rod', A_, prism([circle(0, 2.25, 0.5)], at('SIDE', O), -3, 3), {'cut': 'solid'}, ('strip',))
    out += group('apex_rod', A_, prism([circle(0.125, 2.75, 0.5)], at('SIDE', O), -3, 3), {'common': 'solid'},
                 ('strip',))
    out += group('hole_box', F_, prism([square(-2.5, -2.5, 2.5, 2.5), circle(0, 0, 1.5)], at('XY', O), 0.5, 1.5),
                 {'cut': 'solid', 'common': 'solid'}, ('holed',))
    out += group('stadium', F_, prism([stadium(1.0, 0.75)], at('XY', (0.25, 0.5, 0)), -1, 3),
                 {'cut': 'solid', 'common': 'solid'}, None)
    out += group('pipe_tilt', cone(2, 1, 2, at('TILT', O)), prism([circle(0, 0, 1.5)], at('TILT', O), -1, 3),
                 {'common': 'solid'}, ('axial',))
    out += group('ball_coax', F_, sphere(1.5, at('XY', (0, 0, 2))), ALL, ('axial',))
    out += group('ball_apex', A_, sphere(1, at('XY', (0, 0, 3))), {'common': 'solid'}, ('axial',))
    out += group('ball_side', F_, sphere(1, at('XY', (1.75, 0.5, 1))),
                 {'cut': ('degenerate', TOUCHING), 'common': ('degenerate', TOUCHING)},
                 ('axial',))
    out += group('dome_cone', sphere(2, at('XY', O), 0.0, HP), cone(1.25, 0, 3.5, at('XY', (0, 0, -0.5))),
                 {'cut': 'solid', 'common': 'solid'}, ('axial',))
    out += group('ball_tilt', cone(2, 1, 2, at('TILT', O)), sphere(1, at('XY', (1.5, -0.25, 1.25))),
                 {'common': 'solid'}, ('axial',))
    out += group('cones_apex', A_, cone(1.5, 0, 3, at('DOWN', (0, 0, 4))), {'fuse': 'solid', 'common': 'solid'},
                 ('axial',))
    out += group('frusta_cross', F_, cone(1, 2, 2, at('XY', O)), {'common': 'solid'}, ('axial',))
    out += group('frusta_stack', F_, cone(1.5, 0.5, 1.5, at('XY', (0, 0, 2))), {'fuse': 'solid', 'common': 'empty'},
                 ('axial',))
    crossing = cone(0.75, 0.25, 4, at('SIDE', (-2.25, 0.25, 1)))
    out += group('cones_cross', F_, crossing, ALL, None)
    out += group('cones_ends', crossing, F_, {'cut': 'solid'}, None)
    out += group('cones_lean', F_, cone(1.5, 0, 3, at('LEAN', (0.5, 0.25, -0.5))), {'common': 'solid'}, None)
    out += group('cones_parallel', F_, cone(1.5, 0, 3, at('XY', (1.75, 0.5, -0.5))), {'common': 'solid'},
                 ('axial',))
    out += group('rod_ruling', cone(3, 0, 4, at('XY', O)), prism([circle(0, 1.25, 1.25)], at('LEAN', (-3, 0, 0)), -1, 6),
                 {'common': ('degenerate', RULING)}, None)
    out += group('ball_inscribed', cone(3, 0, 4, at('XY', O)), sphere(0.75, at('XY', (0, 0, 2.75))),
                 {'cut': ('degenerate', INSCRIBED)}, ('axial',))
    out += group('apex_to_apex', A_, cone(1.5, 0, 3, at('DOWN', (0, 0, 6))), {'fuse': ('degenerate', APEXES)},
                 ('axial',))
    out += group('cones_equal', F_, F_, {'fuse': ('degenerate', EQUAL)}, ('equal',))
    return out


def exact_pair(case):
    """Whether the closed forms hold exactly: no cone or prism in a turned
    frame (a whole sphere's model does not depend on its frame's axes)."""
    return all(frame_name(f) in EXACT for f, s in zip(case.frames, case.specs) if s[0] != 'sphere')


# ------------------------------------------------------------------ closed forms

def local_frame(frame):
    x, y, n = ideal_axes(frame_name(frame))
    return tuple(mp.mpf(v) for v in frame[:3]), (x, y, n)


def to_local(o, axes, X):
    return tuple(sum((X[j]-o[j])*axes[i][j] for j in range(3)) for i in range(3))


def to_world(o, axes, V, m):
    return tuple(o[j]*V+sum(m[i]*axes[i][j] for i in range(3)) for j in range(3))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cone_measures(r0, r1, h):
    r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
    V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3
    Mw = mp.pi*h*h*(r0*r0+2*r0*r1+3*r1*r1)/12
    wall = mp.pi*(r0+r1)*mp.sqrt(h*h+(r1-r0)**2)
    return V, Mw, wall+mp.pi*(r0*r0+r1*r1)


def roots_on(f, lo, hi, n=4000):
    """Sign changes of `f` on `[lo, hi]` refined by bisection (kinks)."""
    out = []
    if hi <= lo:
        return out
    xs = [lo+(hi-lo)*k/n for k in range(n+1)]
    vals = [f(x) for x in xs]
    for k in range(n):
        a, b = xs[k], xs[k+1]
        fa, fb = vals[k], vals[k+1]
        if fa == 0:
            out.append(a)
            continue
        if fa*fb < 0:
            for _ in range(200):
                m = (a+b)/2
                fm = f(m)
                if fa*fm <= 0:
                    b = m
                else:
                    a, fa = m, fm
            out.append((a+b)/2)
    for k in range(1, n):
        if abs(vals[k]) <= abs(vals[k-1]) and abs(vals[k]) <= abs(vals[k+1]) and vals[k-1]*vals[k+1] > 0:
            a, b = xs[k-1], xs[k+1]
            for _ in range(200):
                m1, m2 = a+(b-a)/3, b-(b-a)/3
                if abs(f(m1)) < abs(f(m2)):
                    b = m2
                else:
                    a = m1
            m = (a+b)/2
            if abs(f(m)) < mp.mpf(10)**-25:
                out.append(m)
    return out


class Axial:
    """An input whose sections by planes normal to the cone's axis are discs
    of radius `rho(z)` about a point at distance `e` from the axis, for `z`
    in `[z0, z1]`: its lateral element per angle `lateral(z)`, its end discs
    `[(z, outward sign, radius)]`, its measures (volume, world moments,
    area)."""

    def __init__(self, rho, z0, z1, e, u, lateral, discs, measures, turning=()):
        self.rho, self.z0, self.z1, self.e, self.u = rho, z0, z1, e, u
        self.lateral, self.discs, self.measures, self.turning = lateral, discs, measures, turning


def axial_of(spec, o, axes):
    """An input as `Axial` in the cone's ideal frame `o, axes` (its axis
    parallel to the cone's), or None."""
    X, Y, N = axes
    if spec[0] == 'sphere':
        _, R, frame, low, high = spec
        R = mp.mpf(R)
        c = to_local(o, axes, tuple(mp.mpf(v) for v in frame[:3]))
        lo = -R if low == -HP else R*mp.sin(low)
        hi = R if high == HP else R*mp.sin(high)
        if low != -HP or high != HP:
            sn = ideal_axes(frame_name(frame))[2]
            assert abs(dot(sn, N)-1) < mp.mpf(10)**-30
        e = mp.sqrt(c[0]**2+c[1]**2)
        u = (c[0]/e, c[1]/e) if e > 0 else (mp.mpf(1), Z)
        rho = lambda z: mp.sqrt(max(R*R-(z-c[2])**2, Z))
        discs = [(c[2]+lo, -1, rho(c[2]+lo))] if low != -HP else []
        discs += [(c[2]+hi, 1, rho(c[2]+hi))] if high != HP else []
        V = mp.pi*((R*R*hi-hi**3/3)-(R*R*lo-lo**3/3))
        Mz = mp.pi*((R*R*hi**2/2-hi**4/4)-(R*R*lo**2/2-lo**4/4))
        area = 2*mp.pi*R*(hi-lo)+sum(mp.pi*r*r for _, _, r in discs)
        cw = tuple(mp.mpf(v) for v in frame[:3])
        sn = ideal_axes(frame_name(frame))[2]
        mom = tuple(V*cw[i]+Mz*sn[i] for i in range(3))
        return Axial(rho, c[2]+lo, c[2]+hi, e, u, lambda z: R, discs, (V, mom, area), turning=(c[2]-R, c[2]+R))
    if spec[0] == 'cone':
        _, r0, r1, h, frame = spec
        po, (px, py, pn) = local_frame(frame)
        s = dot(pn, N)
        if abs(abs(s)-1) > mp.mpf(10)**-30:
            return None
        base = to_local(o, axes, po)
        r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
        k = (r1-r0)/h
        e = mp.sqrt(base[0]**2+base[1]**2)
        u = (base[0]/e, base[1]/e) if e > 0 else (mp.mpf(1), Z)
        if s > 0:
            z0, z1 = base[2], base[2]+h
            rho = lambda z: r0+k*(z-base[2])
            discs = [(z0, -1, r0), (z1, 1, r1)]
        else:
            z0, z1 = base[2]-h, base[2]
            rho = lambda z: r0+k*(base[2]-z)
            discs = [(z1, 1, r0), (z0, -1, r1)]
        discs = [d for d in discs if d[2] > 0]
        V, Mw, A = cone_measures(r0, r1, h)
        mom = tuple(V*po[i]+Mw*pn[i] for i in range(3))
        sk = mp.sqrt(1+k*k)
        return Axial(rho, z0, z1, e, u, lambda z: sk*rho(z), discs, (V, mom, A))
    _, bounds, frame, w0, w1 = spec
    po, (px, py, pn) = local_frame(frame)
    s = dot(pn, N)
    if len(bounds) != 1 or bounds[0].circle is None or abs(abs(s)-1) > mp.mpf(10)**-30:
        return None
    cx, cy, a = (mp.mpf(v) for v in bounds[0].circle)
    w0, w1 = mp.mpf(w0), mp.mpf(w1)
    cw0 = tuple(po[i]+cx*px[i]+cy*py[i]+w0*pn[i] for i in range(3))
    cw1 = tuple(po[i]+cx*px[i]+cy*py[i]+w1*pn[i] for i in range(3))
    l0, l1 = to_local(o, axes, cw0), to_local(o, axes, cw1)
    z0, z1 = min(l0[2], l1[2]), max(l0[2], l1[2])
    e = mp.sqrt(l0[0]**2+l0[1]**2)
    u = (l0[0]/e, l0[1]/e) if e > 0 else (mp.mpf(1), Z)
    V = mp.pi*a*a*(w1-w0)
    cen = tuple((cw0[i]+cw1[i])/2 for i in range(3))
    area = 2*mp.pi*a*(w1-w0)+2*mp.pi*a*a
    return Axial(lambda z: a, z0, z1, e, u, lambda z: a, [(z0, -1, a), (z1, 1, a)],
                 (V, tuple(V*c for c in cen), area))


def axial_form(r0, r1, h, B):
    """The cone (ideal frame, base at 0) and an `Axial` input: the common's
    volume and moments, the cone's surface inside, the other's inside, the
    shared discs of the same and opposite orientation."""
    r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
    k = (r1-r0)/h
    rK = lambda z: r0+k*(z-0)
    sk = mp.sqrt(1+k*k)
    lo, hi = max(Z, B.z0), min(h, B.z1)
    L = lambda z: lens(max(rK(z), Z), B.rho(z), B.e)
    out = {'volume': Z, 'moments': (Z, Z, Z), 'K_in': Z, 'B_in': Z, 'same': Z, 'opp': Z}
    if hi > lo:
        kinks = {lo, hi}
        for sg1, sg2 in ((1, 1), (1, -1), (-1, 1)):
            f = lambda z, a=sg1, b=sg2: a*rK(z)+b*B.rho(z)-B.e
            kinks.update(roots_on(f, lo, hi))
        kinks.update(t for t in B.turning if lo < t < hi)
        pts = sorted(kinks)
        V = mp.quad(lambda z: L(z)[0], pts)
        Mu = mp.quad(lambda z: L(z)[1], pts)
        Mz = mp.quad(lambda z: z*L(z)[0], pts)
        out['volume'] = V
        out['moments'] = (B.u[0]*Mu, B.u[1]*Mu, Mz)
        out['K_in'] = mp.quad(lambda z: sk*rK(z)*L(z)[2], pts)
        out['B_in'] = mp.quad(lambda z: B.lateral(z)*L(z)[3], pts)
    for z, sign, r in ((Z, -1, r0), (h, 1, r1)):
        if r <= 0:
            continue
        co = [d for d in B.discs if d[0] == z]
        if co:
            zb, sb, rb = co[0]
            shared = lens(r, rb, B.e)[0]
            out['same' if sb == sign else 'opp'] += shared
        elif B.z0 < z < B.z1:
            out['K_in'] += lens(r, B.rho(z), B.e)[0]
    for zb, sb, rb in B.discs:
        if zb in (Z, h) and ((zb == 0 and r0 > 0) or (zb == h and r1 > 0)):
            continue
        if 0 < zb < h:
            out['B_in'] += lens(rK(zb), rb, B.e)[0]
    return out


def strip_form(r0, r1, h, y0, z0, a):
    """The cone (ideal frame) and the cylinder of radius `a` whose axis is
    the line `y = y0, z = z0` along `x` (its caps outside the cone)."""
    r0, r1, h, y0, z0, a = (mp.mpf(v) for v in (r0, r1, h, y0, z0, a))
    k = (r1-r0)/h
    rK = lambda z: max(r0+k*z, Z)
    sk = mp.sqrt(1+k*k)
    q = lambda z: mp.sqrt(max(a*a-(z-z0)**2, Z))
    lo, hi = max(Z, z0-a), min(h, z0+a)
    out = {'volume': Z, 'moments': (Z, Z, Z), 'K_in': Z, 'B_in': Z, 'same': Z, 'opp': Z}
    if hi <= lo:
        return out

    def band(z):
        r, qq = rK(z), q(z)
        s1, s2 = segment(r, y0-qq), segment(r, y0+qq)
        return [s1[i]-s2[i] for i in range(3)]
    kinks = {lo, hi}
    for sg in (1, -1):
        for sr in (1, -1):
            f = lambda z, sg=sg, sr=sr: (y0+sg*q(z))-sr*rK(z)
            kinks.update(roots_on(f, lo, hi))
    pts = sorted(kinks)
    out['volume'] = mp.quad(lambda z: band(z)[0], pts)
    out['moments'] = (Z, mp.quad(lambda z: band(z)[1], pts), mp.quad(lambda z: z*band(z)[0], pts))
    out['K_in'] = mp.quad(lambda z: sk*rK(z)*band(z)[2], pts)
    # The cylinder's wall by its generatrices: at angle phi the line y = y0
    # + a cos(phi), z = z0 + a sin(phi), its chord inside the disc.

    def chord(ph):
        z = z0+a*mp.sin(ph)
        if not 0 < z < h:
            return Z
        y, r = y0+a*mp.cos(ph), rK(z)
        return 2*mp.sqrt(r*r-y*y) if y*y < r*r else Z
    turns = {Z, 2*mp.pi}
    for zz in (Z, h):
        if abs(zz-z0) < a:
            v = mp.asin((zz-z0)/a)
            turns.update(t % (2*mp.pi) for t in (v, mp.pi-v))
    for sr in (1, -1):
        turns.update(roots_on(lambda ph, sr=sr: y0+a*mp.cos(ph)-sr*rK(z0+a*mp.sin(ph)), Z, 2*mp.pi))
    out['B_in'] = a*mp.quad(chord, sorted(turns))
    for z, r in ((Z, r0), (h, r1)):
        if r > 0 and abs(z-z0) < a:
            out['K_in'] += segment(r, y0-q(z))[0]-segment(r, y0+q(z))[0]
    return out


def assemble(A, B, common, A_in, B_in, same=Z, opp=Z):
    """Every operation by inclusion and exclusion (S9d.3a's, coplanar faces
    apart)."""
    VA, MA, AA = A
    VB, MB, AB = B
    Vc, Mc = common
    A_out, B_out = AA-A_in-same-opp, AB-B_in-same-opp
    ops = {'common': (Vc, Mc, A_in+B_in+same),
           'fuse': (VA+VB-Vc, tuple(a+b-c for a, b, c in zip(MA, MB, Mc)), A_out+B_out+same),
           'cut': (VA-Vc, tuple(a-c for a, c in zip(MA, Mc)), A_out+opp+B_in)}
    out = {}
    for op, (V, m, area) in ops.items():
        out[op] = None if V <= mp.mpf(10)**-30 else (V, area, tuple(x/V for x in m))
    return out


def closed_forms(first):
    """{op: (volume, area, centre) or None}, or None where no closed form is
    taken."""
    if first.form is None:
        return None
    obj, tool = first.specs
    ks = [s for s in (obj, tool) if s[0] == 'cone']
    kspec = ks[0]
    other = tool if kspec is obj else obj
    _, r0, r1, h, kframe = kspec
    o, axes = local_frame(kframe)
    V, Mw, AK = cone_measures(r0, r1, h)
    K = (V, tuple(V*o[i]+Mw*axes[2][i] for i in range(3)), AK)
    kind = first.form[0]
    if kind == 'equal':
        forms = assemble(K, K, (K[0], K[1]), Z, Z, same=AK)
        return forms
    if kind == 'axial':
        B = axial_of(other, o, axes)
        assert B is not None, first.name
        f = axial_form(r0, r1, h, B)
        Bm = B.measures
    elif kind == 'strip':
        _, bounds, frame, w0, w1 = other
        po, (px, py, pn) = local_frame(frame)
        assert abs(dot(pn, axes[2])) < mp.mpf(10)**-30 and abs(dot(pn, axes[1])) < mp.mpf(10)**-30
        cx, cy, a = (mp.mpf(v) for v in bounds[0].circle)
        c = to_local(o, axes, tuple(po[i]+cx*px[i]+cy*py[i] for i in range(3)))
        ends = [to_local(o, axes, tuple(po[i]+cx*px[i]+cy*py[i]+mp.mpf(w)*pn[i] for i in range(3)))[0]
                for w in (w0, w1)]
        assert all(abs(x) > max(r0, r1)+a for x in ends), 'the caps outside the cone'
        f = strip_form(r0, r1, h, c[1], c[2], a)
        L = mp.mpf(w1)-mp.mpf(w0)
        VB = mp.pi*a*a*L
        cen = tuple(po[i]+cx*px[i]+cy*py[i]+(mp.mpf(w0)+mp.mpf(w1))/2*pn[i] for i in range(3))
        Bm = (VB, tuple(VB*v for v in cen), 2*mp.pi*a*L+2*mp.pi*a*a)
    elif kind == 'holed':
        _, bounds, frame, w0, w1 = other
        sq, hole = bounds
        (u0, v0), (u1, v1) = sq.points[0], sq.points[2]
        hx, hy, a = (mp.mpf(v) for v in hole.circle)
        assert frame == kframe and hx == 0 and hy == 0
        box = [(u0, u1), (v0, v1), (w0, w1)]
        fb = box_form(r0, r1, h, box)
        core = axial_of(prism([circle(0, 0, a)], frame, w0, w1), o, axes)
        fc = axial_form(r0, r1, h, core)
        f = {'volume': fb['volume']-fc['volume'],
             'moments': tuple(x-y for x, y in zip(fb['moments'], fc['moments'])),
             'K_in': fb['cone_in']-fc['K_in'], 'same': Z, 'opp': Z}
        caps = sum((lens(max(mp.mpf(r0)+(mp.mpf(r1)-mp.mpf(r0))/mp.mpf(h)*z, Z), a, Z)[0]
                    for z, _, _ in core.discs if 0 < z < h), Z)
        wall = fc['B_in']-caps
        f['B_in'] = fb['box_in']-caps+wall
        su, sv, L = mp.mpf(u1-u0), mp.mpf(v1-v0), mp.mpf(w1)-mp.mpf(w0)
        Ap = su*sv-mp.pi*a*a
        VB = Ap*L
        mu = mp.mpf(u0+u1)/2*su*sv
        mv = mp.mpf(v0+v1)/2*su*sv
        Bm = (VB, to_world(o, axes, VB, (mu*L, mv*L, (mp.mpf(w0)+mp.mpf(w1))/2*VB)),
              2*Ap+2*(su+sv)*L+2*mp.pi*a*L)
    Vc = f['volume']
    Mc = to_world(o, axes, Vc, f['moments'])
    same, opp = f['same'], f['opp']
    if kspec is obj:
        return assemble(K, Bm, (Vc, Mc), f['K_in'], f['B_in'], same, opp)
    return assemble(Bm, K, (Vc, Mc), f['B_in'], f['K_in'], same, opp)


# ------------------------------------------------------------------ Monte Carlo

def float_input(x):
    if x.kind == 'cone':
        return FloatCone(x)
    if x.kind == 'sphere':
        return FloatBall(x.S)
    return FloatArcPrism(x.P)


def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fa, fb = float_input(pair.A), float_input(pair.B)
    box0 = [min(fa.box[0][i], fb.box[0][i]) for i in range(3)]
    box1 = [max(fa.box[1][i], fb.box[1][i]) for i in range(3)]
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


# ------------------------------------------------------------------ clearances and sub-steps

def surface_scale(x, s, size):
    """A surface's squared radius (quadrics) or its normal's length times
    the case's size (planes), to make its equation's values lengths
    squared."""
    if s.plane:
        a = s.a()
        return mp.sqrt(ref.M(ref.dot(a, a)))*size
    if x.kind == 'cone':
        return ref.M(max(x.r0, x.r1))**2
    if x.kind == 'sphere':
        return ref.M(x.r2)
    for e in x.P.elements:
        if e[0] == 'arc' and x.surfaces[x.el_surf[e[-2]]] is s:
            return ref.M(e[2])**2
    return size**2


def clearance(pair):
    """The smallest critical value of each surface's equation along every
    edge of the other input (and at its vertices), relative to
    `surface_scale`: incidences and tangencies of edges with surfaces."""
    best = None
    size = pair.size

    def note(v):
        nonlocal best
        best = v if best is None else min(best, v)
    for x, y in ((pair.A, pair.B), (pair.B, pair.A)):
        for s in y.surfaces:
            sc = surface_scale(y, s, size)
            for v in x.vertices():
                note(abs(s.value(ref.Mv(v)))/sc)
            for e in x.edges():
                if on_plane(e, s):
                    # An edge in a coplanar face of the other (a stack).
                    continue
                if isinstance(e, ref.LineEdge):
                    c0, c1, c2 = (ref.M(c) for c in s.line(e.P0, e.e))
                    ts = [Z, mp.mpf(1)]
                    if c2 != 0:
                        t = -c1/(2*c2)
                        if 0 < t < 1:
                            ts.append(t)
                    for t in ts:
                        note(abs(c0+c1*t+c2*t*t)/sc)
                else:
                    for v in conic_minima(e, s):
                        note(v/sc)
    return best


def on_plane(e, s):
    """Whether an edge lies in a plane surface."""
    if not s.plane:
        return False
    a = s.a()
    if isinstance(e, ref.LineEdge):
        return s.value(e.P0) == 0 and ref.dot(a, e.e) == 0
    return s.value(e.P) == 0 and ref.dot(a, e.E1) == 0 and ref.dot(a, e.E2) == 0


def conic_minima(e, s, n=720):
    """`|Q|` of a surface at the local extrema of `Q` along a conic edge (a
    circle about the plane chart's origin), sampled and refined by golden
    sections; a prism's arc only within its span and at its ends."""
    Pm, E1, E2 = ref.Mv(e.P), ref.Mv(e.E1), ref.Mv(e.E2)
    A = tuple(ref.Pol([e.P[i], e.E1[i]]) for i in range(3))
    B = tuple(ref.Pol([v]) for v in e.E2)
    f0, f1, f2 = (x.value(Z) for x in e.surf.family(A, B, 'pol'))

    def point(th):
        c, s_ = mp.cos(th), mp.sin(th)
        d = tuple(c*E1[i]+s_*E2[i] for i in range(3))
        cc0, cc1, cc2 = e.surf.num_line(Pm, d)
        rho = max(ref.solve_quadratic(cc0, cc1, cc2))
        return tuple(Pm[i]+rho*d[i] for i in range(3))
    g = lambda th: s.value(point(th))
    if e.arc is not None:
        _, t0, t1 = e.arc
        lo, hi = min(t0, t1), max(t0, t1)
    else:
        lo, hi = Z, 2*mp.pi
    xs = [lo+(hi-lo)*k/n for k in range(n+1)]
    vals = [g(x) for x in xs]
    out = [abs(vals[0]), abs(vals[-1])] if e.arc is not None else []
    cyclic = e.arc is None
    for k in range(0 if cyclic else 1, n):
        prev, nxt = vals[(k-1) % n], vals[k+1]
        for sg in (1, -1):
            if sg*vals[k] <= sg*prev and sg*vals[k] <= sg*nxt:
                a, b = xs[k]-(hi-lo)/n, xs[k+1]
                for _ in range(80):
                    m1, m2 = a+(b-a)*0.382, b-(b-a)*0.382
                    if sg*g(m1) < sg*g(m2):
                        b = m2
                    else:
                        a = m1
                out.append(abs(g((a+b)/2)))
    return out


def ruling_tangencies(x, y, size):
    """Whether some ruling of `x`'s curved faces is tangent to a quadric of
    `y` at a point on both faces."""
    for face in x.faces():
        if not isinstance(face, ref.Ruled) or face.ring != 'trig':
            continue
        for j, s in y.quad_surfaces():
            c0, c1, c2 = s.family(face.A, face.B, 'trig')
            disc = c1*c1+ref.neg((c0*c2).scale(ref.F(4)))
            for th in ref.ring_roots(disc, face.t0, face.t1, 'trig'):
                A, B = face.num(th)
                a0, a1, a2 = (v.value(th) for v in (c0, c1, c2))
                if a2 == 0:
                    continue
                w = -a1/(2*a2)
                if not ref.M(face.wr[0]) <= w <= ref.M(face.wr[1]):
                    continue
                X = tuple(A[i]+w*B[i] for i in range(3))
                if on_face(y, j, X):
                    return True
    return False


def on_face(y, j, X):
    """Whether a point of `y`'s quadric `j` lies on its face (within its
    other bounds)."""
    if y.kind == 'cone':
        return 0 <= y.uvw(X)[2] <= ref.M(y.h)
    if y.kind == 'sphere':
        return all(s.value(X) >= 0 for s in y.surfaces[1:])
    u, v, w = y.uvw(X)
    if not ref.M(y.P.lo) <= w <= ref.M(y.P.hi):
        return False
    for e in y.P.elements:
        if e[0] == 'arc' and y.el_surf[e[-2]] == j:
            _, C, r, t0, t1, full, _, _ = e
            if full:
                return True
            ang = mp.atan2(v-ref.M(C[1]), u-ref.M(C[0]))
            return ref.s2.span_param(t0, t1, ang) is not None
    return False


def substep(pair, kind):
    """`S9d.3b.1` or `S9d.3b.2` (see the module's docstring)."""
    if pair.identical:
        return STEP+'.1'
    size = pair.size
    a = ruling_tangencies(pair.A, pair.B, size)
    b = ruling_tangencies(pair.B, pair.A, size)
    has_a = any(isinstance(f, ref.Ruled) and f.ring == 'trig' for f in pair.A.faces())
    has_b = any(isinstance(f, ref.Ruled) and f.ring == 'trig' for f in pair.B.faces())
    graph = (has_a and not a) or (has_b and not b)
    return STEP+('.1' if graph else '.2')


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, first, ops, mc_n = job
    obj, tool = first.obj, first.tool
    pair = ref.Pair(obj, tool)
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    size = pair.size
    r = pair.sliced()
    VA, MA, AA = pair.A.closed()
    VB, MB, AB = pair.B.closed()
    checks = {}
    dev = mp.mpf(0)
    for (v, m), (V, Mm) in ((r['A'], (VA, MA)), (r['B'], (VB, MB))):
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inputs_sliced'] = dev
    vc, mc = r['common']
    dev = mp.mpf(0)
    for key, (V, Mm) in (('fuse', (VA+VB-vc, [a+b-c for a, b, c in zip(MA, MB, mc)])),
                         ('A-B', (VA-vc, [a-c for a, c in zip(MA, mc)])),
                         ('B-A', (VB-vc, [a-c for a, c in zip(MB, mc)]))):
        v, m = r[key]
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inclusion_exclusion'] = dev
    faces = pair.faces()
    checks['face_classes'] = max(abs(sum(cls.values())-area) for _, _, cls, area in faces)/size**2
    ca, cb = pair.classes('A'), pair.classes('B')
    checks['shared_both_sides'] = max(abs(ca['same']-cb['same']), abs(ca['opp']-cb['opp']))/size**2
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')+2*ca['opp']-AA-AB)/size**2
    checks['input_areas'] = max(abs(sum(sum(cls.values()) for w, _, cls, _ in faces if w == k)-A)
                                for k, A in (('A', AA), ('B', AB)))/size**2
    if not pair.identical:
        dev = None
        for (role, tag, cls, _) in faces:
            key = ('curved', role)
            x = pair.A if role == 'A' else pair.B
            if tag == ('wall',) and x.kind == 'cone' and key in r:
                d = max(abs(r[key]['in']-cls['in']), abs(r[key]['out']-cls['out']))/size**2
                dev = d if dev is None else max(dev, d)
        if dev is not None:
            checks['wall_two_ways'] = dev
        caps = any(x.kind == 'sphere' and x.ball.S.planes for x in pair.inputs)
        # Coplanar faces make a tilted slice's boundaries share a line.
        caps = caps or any(ref.same_plane(p, q)[0] for _, p in pair.A.plane_surfaces()
                           for _, q in pair.B.plane_surfaces())
        # A tangency along a curve leaves a tilted slice's sections tangent
        # at every level (a degenerate pencil).
        caps = caps or any(c.kind == 'degenerate' for c in cases() if c.pair_name == name)
        second = None
        if not caps:
            try:
                second = ref.Pair(obj, tool, avoid=pair.chart.d)
            except AssertionError:
                second = None
        if second is not None:
            s2r = second.sliced()
            dev = mp.mpf(0)
            for key in ('A', 'B', 'common', 'fuse', 'A-B', 'B-A'):
                dev = max(dev, abs(s2r[key][0]-r[key][0])/size**3,
                          max(abs(x-y) for x, y in zip(s2r[key][1], r[key][1]))/size**4)
            for x in pair.inputs:
                key = ('curved', x.role)
                if x.kind == 'sphere':
                    for k in ('in', 'out'):
                        dev = max(dev, abs(s2r[key][k]-r[key][k])/size**2)
            checks['second_direction'] = dev
    mcr = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mcr[op]
        V, m = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in m]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(pair)
    clear = clearance(pair)
    if clear is not None and clear < mp.mpf(10)**-3:
        near.append(f'an edge or vertex within {mp.nstr(clear, 3)} of a surface')
    step = substep(pair, first.kind)
    stats = {'volume_quadrature': pair.quad_error/size**4, 'breaks': len(getattr(pair, 'breaks', [])),
             'sweep_breaks': pair.sweep_breaks, 'clearance': clear, 'step': step,
             'chart': None if pair.chart is None else tuple(str(c) for c in pair.chart.d)}
    return name, rows, res, checks, near, stats


def evaluate_safe(job):
    """`evaluate`, a failure reported with its pair's name."""
    import traceback
    try:
        return evaluate(job)
    except Exception:
        return ('error', job[0], traceback.format_exc())


def near_coincidences(pair):
    """A face class, or a result's volume, positive but thinner than 1e-9 of
    the case's size; slicing breakpoints closer than that (and farther than
    1e-30)."""
    size = pair.size
    out = []
    for w, tag, cls, _ in pair.faces():
        for c, v in cls.items():
            if 0 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{w} {tag} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volume(op)[0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    b = getattr(pair, 'raw_breaks', [])
    for p, q in zip(b, b[1:]):
        if mp.mpf(10)**-30*size < q-p < mp.mpf(10)**-9*size:
            out.append(f'breakpoints {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, closed_forms(first)


PER_PAIR = {}
LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
          'inputs_sliced': 1e-30, 'inclusion_exclusion': 1e-30, 'face_classes': 1e-30,
          'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30, 'wall_two_ways': 1e-30,
          'second_direction': 1e-30, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-30}


def reference_checks(results, forms):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1
    by_pair = {r[0]: r for r in results}
    firsts = {}
    for c in cases():
        firsts.setdefault(c.pair_name, c)
    for name, form in forms.items():
        if form is None:
            continue
        _, _, res, _, _, _ = by_pair[name]
        kind = 'exact' if exact_pair(firsts[name]) else 'turned'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            assert n > 0, (name, op, 'the reference is empty')
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}_frames', dev)
            PER_PAIR[name] = max(PER_PAIR.get(name, mp.mpf(0)), dev)
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, cones_boolean_reference.py: expect KIND STEP (S9d.3b.1 circles and rings, '
           'S9d.3b.2 loops), reason TEXT for a degenerate case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, _, _, stats = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {stats["step"]}')
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
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items()]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert c.obj.cone is not None or c.tool.cone is not None, f'{c.name}: a cone'
        for s in (c.obj, c.tool):
            if s.sphere is not None and (s.sphere[1] != -HP or s.sphere[2] != HP):
                assert frame_name(s.frame) in EXACT and s.frame[:3] == (0.0, 0.0, 0.0), \
                    f'{c.name}: a cap or zone in an exact frame at the origin'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', help='evaluate only these pairs, comma separated (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        wanted = set(args.only.split(','))
        todo = [j for j in todo if j[0] in wanted]
    results = run(todo, evaluate_safe, args.workers)
    errors = [r for r in results if r[0] == 'error']
    for _, name, tb in errors:
        print('FAILED', name, tb[-1500:])
    if errors:
        raise SystemExit(f'{len(errors)} pairs failed')
    degenerate = {c.pair_name for c in listed if c.kind == 'degenerate'}
    near = [f'{name}: near coincidences {n}' for name, _, _, _, n, _ in results if n and name not in degenerate]
    firsts = {}
    for c in listed:
        firsts.setdefault(c.pair_name, c)
    done = {r[0] for r in results}
    forms = dict(run([f for k, f in firsts.items() if k in done], closed_job, args.workers))
    worst, covered = reference_checks(results, forms)
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    for name, rows, _, checks, n, stats in results:
        print(name, stats['step'], 'closed', mp.nstr(PER_PAIR[name], 3) if name in PER_PAIR else None, 'chart', stats['chart'], 'breaks', stats['breaks'], 'sweep', stats['sweep_breaks'],
              'clearance', None if stats['clearance'] is None else mp.nstr(stats['clearance'], 3),
              '; '.join(f'{op} {r[0]}' for op, r in rows.items()))
        if n:
            print('  near coincidences:', '; '.join(n[:4]))
    if near:
        raise SystemExit('\n'.join(near))
    for key, value in worst.items():
        assert value <= LIMITS[key], (key, mp.nstr(value, 3))
    if args.only:
        return
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
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames (whole spheres in any)')


if __name__ == '__main__':
    main()
