#!/usr/bin/env python3
"""Fixtures for S9d.1 of REVIEW_NOTES.md: Booleans of a sphere (whole, a cap
or a zone, `Solid::sphere_with`) against a polyhedral prism in any relative
position, before any of its kernel code.

`boolean-sphere-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a sphere's block is its `frame`
and `sphere R LOW HIGH` rows); `boolean-sphere-expected.tsv` gives per case
`expect KIND S9d.1` (`solid`, `empty` or `degenerate`, then `reason TEXT`)
and `result N volume area cx cy cz` or `empty` from
`sphere_boolean_reference.py`, as the curved fixtures;
`boolean-sphere-frames.tsv` the stored axes' bits. Frames are the curved
generator's (`FRAMES`, their stored axes the kernel's bit for bit). Caps
and zones are only in exact frames centred at the world's origin, where
every reading of their end planes (the affine `w = h`, the plane through
`o + h n` normal to the stored `n`) is the same plane exactly; whole
spheres also in turned ones.

Pairs: a sphere through a box's face (`face_cap`, the common a cap); a
box's corner in a sphere (`corner_box`; `corner_tilt`, the box in `TILT`;
`corner_lean`, the box in `LEAN` and the sphere in `TILTX`); a box's edge
through a sphere (`edge_wedge`); a square post through a zone's flat end
(`zone_post`); a box through a zone's flat end and its band, the zone
along `x` (`zone_side`); a hemisphere (a cap to the equator) crossed by a
box through its disc and dome (`hemisphere_box`); a dome on a box's top
face, their discs coplanar with opposite orientations (`dome_on_box`, the
common empty) and inside a box on its bottom face, the same orientation
(`hemisphere_in_box`); a sphere inside a box and a box inside a sphere,
cavities (`sphere_in_box`, `box_in_sphere`); a bar in `TILT` through a
sphere (`bar_tilt`, the common; `bar_through`, the bar less the sphere,
two solids); a thin slab cutting a sphere in two (`slab`); an L-shaped
prism severed at its corner (`lshape`, its profile ear-clipped into four
triangles, two solids); an octant, the box's corner at the centre
(`octant`); and `degenerate` with reasons: a face tangent to the sphere
(`face_tangent`), the box's eight vertices on it (`inscribed`), an edge
tangent to it (`edge_tangent`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms of every pair, by one quadrature (mpmath's tanh-sinh between
  the kinks) of exact 2D forms in the box's own frame: a rectangle's part of
  a disc (closed-form antiderivatives in `x` between the breakpoints where
  the circle crosses the rectangle's lines), the angle of a circle inside a
  rectangle, each box face's part of the sphere's section; so the common's
  volume, moments and both inputs' surfaces inside the other, and every
  operation by inclusion and exclusion (coplanar discs by hand); within
  1e-30 in exact frames and 1e-15 in turned ones (their stored axes not
  exactly orthonormal);
* `fuse = A + B - common` and `cut = A - common` (each operation sliced
  apart, `A` and `B` in closed form), the common a second way (each piece
  clipped by the disc), both inputs' slicings against their closed forms,
  the area identity `area(fuse) + area(common) + 2 opp = area(A) +
  area(B)`, every face's classes summing to its area, both sides' shared
  areas equal, a second slicing direction `(2, -3, 5)` (pairs without
  coplanar faces), within 1e-30;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, slicing breakpoints closer
  than that) must find none but in the declared degenerate pairs.

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
import sphere_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-sphere'
STEP = 'S9d.1'
BOOLEAN_OPERATION = 93
HP = math.pi/2
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
SECOND = (2, -3, 5)


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def make(spec, op):
    if spec[0] == 'sphere':
        _, r, frame, low, high = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), low, high))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def sphere(r, frame, low=-HP, high=HP):
    return ('sphere', r, frame, low, high)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


class Boolean:
    def __init__(self, name, operation, obj, tool, kind='solid', reason=None, opposite=False, coplanar=False):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.opposite, self.coplanar = opposite, coplanar
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.frames = (obj[2], tool[2])

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, **kw):
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, kind, reason, **kw))
    return out


TANGENT_FACE = 'a face of the prism tangent to the sphere'
VERTICES_ON = 'vertices of the prism on the sphere'
TANGENT_EDGE = 'an edge of the prism tangent to the sphere'
L_PROFILE = Boundary(points=[(0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)])


def cases():
    """S9d.1's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    out = []
    out += group('face_cap', sphere(2, at('XY', O)), prism([square(-3, -3, 3, 3)], at('XY', O), 1, 4),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    corner = [square(1, 0.5, 6, 6)]
    out += group('corner_box', prism(corner, at('XY', O), 1.5, 6), sphere(3, at('XY', O)),
                 {'cut': 'solid', 'common': 'solid'})
    out += group('corner_tilt', prism(corner, at('TILT', O), 1.5, 6), sphere(3, at('XY', O)),
                 {'common': 'solid'})
    out += group('corner_lean', prism(corner, at('LEAN', (1, 2, -1)), 1.5, 6), sphere(3, at('TILTX', (1, 2, -1))),
                 {'cut': 'solid', 'common': 'solid'})
    out += group('edge_wedge', sphere(2.5, at('XY', O)), prism([square(0.5, 0.75, 5, 5)], at('XY', O), -5, 5),
                 {'common': 'solid'})
    out += group('zone_post', sphere(2, at('XY', O), -0.5, 1.0),
                 prism([square(-0.5, -0.5, 0.5, 0.5)], at('XY', O), 1, 4),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('zone_side', sphere(2, at('SIDE', O), -0.5, 1.0),
                 prism([square(1, 0.5, 4, 3)], at('XY', O), -0.5, 0.5), {'common': 'solid'})
    out += group('hemisphere_box', sphere(2, at('XY', O), 0.0, HP),
                 prism([square(-1, -1, 3, 3)], at('XY', O), -1, 1), {'cut': 'solid', 'common': 'solid'})
    out += group('dome_on_box', sphere(1.5, at('XY', O), 0.0, HP),
                 prism([square(-2, -2, 2, 2)], at('XY', O), -1, 0), {'fuse': 'solid', 'common': 'empty'},
                 opposite=True, coplanar=True)
    out += group('hemisphere_in_box', sphere(1.5, at('XY', O), 0.0, HP),
                 prism([square(-2, -2, 2, 2)], at('XY', O), 0, 1), {'cut': 'solid', 'common': 'solid'},
                 coplanar=True)
    out += group('sphere_in_box', prism([square(-2, -2, 2, 3)], at('XY', O), -2, 2),
                 sphere(1, at('XY', (0, 0.25, 0))), {'cut': 'solid'})
    out += group('box_in_sphere', sphere(2, at('XY', O)), prism([square(-1, -1, 1, 1)], at('XY', O), -0.5, 1),
                 {'cut': 'solid'})
    centre = (0.25, -0.5, 0.75)
    bar = prism([square(-0.5, -0.5, 0.5, 0.5)], at('TILT', centre), -4, 4)
    out += group('bar_tilt', sphere(2, at('XY', centre)), bar, {'common': 'solid'})
    out += group('bar_through', bar, sphere(2, at('XY', centre)), {'cut': 'solid'})
    out += group('slab', sphere(2, at('XY', O)), prism([square(-3, -3, 3, 3)], at('XY', O), -0.25, 0.25),
                 {'cut': 'solid'})
    out += group('lshape', prism([L_PROFILE], at('XY', O), -1, 1), sphere(1.5, at('XY', (0.5, 0.5, 0))),
                 {'cut': 'solid'})
    out += group('octant', sphere(2, at('XY', O)), prism([square(0, 0, 5, 5)], at('XY', O), 0, 5),
                 {'common': 'solid'})
    out += group('face_tangent', sphere(2, at('XY', O)), prism([square(-3, -3, 3, 3)], at('XY', O), 2, 4),
                 {'fuse': ('degenerate', TANGENT_FACE), 'common': ('degenerate', TANGENT_FACE)})
    out += group('inscribed', sphere(3, at('XY', O)), prism([square(-1, -2, 1, 2)], at('XY', O), -2, 2),
                 {'cut': ('degenerate', VERTICES_ON)})
    out += group('edge_tangent', sphere(5, at('XY', O)), prism([square(-6, 3, 6, 8)], at('XY', O), 4, 9),
                 {'common': ('degenerate', TANGENT_EDGE)})
    return out


def exact_pair(case):
    return all(frame_name(f) in EXACT for f in case.frames)


# ------------------------------------------------------------------ closed forms

def rect_disc(x0, x1, y0, y1, rho):
    """Area and first moments of `[x0, x1] x [y0, y1]` inside the disc of
    radius `rho` about the origin: closed-form antiderivatives in `x`
    between the breakpoints where the circle crosses the rectangle's lines."""
    zero = (mp.mpf(0), mp.mpf(0), mp.mpf(0))
    if rho <= 0 or x1 <= x0 or y1 <= y0:
        return zero
    r2 = rho*rho
    a, b = max(x0, -rho), min(x1, rho)
    if b <= a:
        return zero
    pts = [a, b]
    for y in (y0, y1):
        if y*y < r2:
            xq = mp.sqrt(r2-y*y)
            pts += [-xq, xq]
    pts = sorted(p for p in pts if a <= p <= b)
    q = lambda x: mp.sqrt(max(r2-x*x, mp.mpf(0)))

    def Q1(x):
        return (x*q(x)+r2*mp.asin(max(min(x/rho, mp.mpf(1)), mp.mpf(-1))))/2
    area = mx = my = mp.mpf(0)
    for p0, p1 in zip(pts, pts[1:]):
        if p1 <= p0:
            continue
        xm = (p0+p1)/2
        qm = q(xm)
        hc, lc = y1 < qm, y0 > -qm
        hi = y1 if hc else qm
        lo = y0 if lc else -qm
        if hi <= lo:
            continue

        def I(x):
            L = (y1*x if hc else Q1(x))-(y0*x if lc else -Q1(x))
            X = (y1*x*x/2 if hc else -q(x)**3/3)-(y0*x*x/2 if lc else q(x)**3/3)
            sq = lambda c: c*c*x if c is not None else r2*x-x**3/3
            Y = (sq(y1 if hc else None)-sq(y0 if lc else None))/2
            return L, X, Y
        e1, e0 = I(p1), I(p0)
        area += e1[0]-e0[0]
        mx += e1[1]-e0[1]
        my += e1[2]-e0[2]
    return area, mx, my


def arc_in_rect(x0, x1, y0, y1, rho):
    """The angle of the circle of radius `rho` about the origin inside the
    rectangle."""
    if rho <= 0:
        return mp.mpf(0)
    T = 2*mp.pi
    angs = []
    for x in (x0, x1):
        if abs(x) < rho:
            t = mp.acos(x/rho)
            angs += [t, -t]
    for y in (y0, y1):
        if abs(y) < rho:
            t = mp.asin(y/rho)
            angs += [t, mp.pi-t]
    inside = lambda t: x0 < rho*mp.cos(t) < x1 and y0 < rho*mp.sin(t) < y1
    if not angs:
        return T if inside(mp.mpf(0)) else mp.mpf(0)
    angs = sorted(a % T for a in angs)
    total = mp.mpf(0)
    for i in range(len(angs)):
        a = angs[i]
        b = angs[i+1] if i+1 < len(angs) else angs[0]+T
        if b > a and inside((a+b)/2):
            total += b-a
    return total


def box_sphere(box, c, r, zone=None):
    """A box `[u0, u1] x [v0, v1] x [w0, w1]` (orthonormal coordinates)
    against the ball of centre `c` and radius `r`, cut to `zone = (lo, hi)`
    in `w` (None: whole): the common's volume and first moments, the
    sphere's surface inside the box (its spherical part and its end discs
    strictly inside the box's `w` range), and the box's faces inside the
    ball and zone (faces in an end plane left out: coplanar, by hand)."""
    (u0, u1), (v0, v1), (w0, w1) = [tuple(mp.mpf(x) for x in rng) for rng in box]
    cu, cv, cw = (mp.mpf(x) for x in c)
    r = mp.mpf(r)
    zlo, zhi = (cw-r, cw+r) if zone is None else (mp.mpf(zone[0]), mp.mpf(zone[1]))
    X0, X1, Y0, Y1 = u0-cu, u1-cu, v0-cv, v1-cv
    rho = lambda w: mp.sqrt(max(r*r-(w-cw)**2, mp.mpf(0)))
    wa, wb = max(w0, zlo, cw-r), min(w1, zhi, cw+r)
    out = {'volume': mp.mpf(0), 'moments': (mp.mpf(0),)*3, 'sphere_in': mp.mpf(0), 'box_in': mp.mpf(0),
           'ends_in': mp.mpf(0)}
    if wb > wa:
        kinks = set()
        for t2 in [X0**2, X1**2, Y0**2, Y1**2]+[a*a+b*b for a in (X0, X1) for b in (Y0, Y1)]:
            if t2 < r*r:
                d = mp.sqrt(r*r-t2)
                kinks.update([cw-d, cw+d])
        pts = [wa]+sorted(k for k in kinks if wa < k < wb)+[wb]
        rd = lambda w: rect_disc(X0, X1, Y0, Y1, rho(w))
        V = mp.quad(lambda w: rd(w)[0], pts)
        Mu = mp.quad(lambda w: rd(w)[1], pts)+cu*V
        Mv = mp.quad(lambda w: rd(w)[2], pts)+cv*V
        Mw = mp.quad(lambda w: w*rd(w)[0], pts)
        out['volume'], out['moments'] = V, (Mu, Mv, Mw)
        out['sphere_in'] = r*mp.quad(lambda w: arc_in_rect(X0, X1, Y0, Y1, rho(w)), pts)
    # The box's faces.
    wlo, whi = max(w0, zlo), min(w1, zhi)
    box_in = mp.mpf(0)
    for u in (u0, u1):
        s2 = r*r-(u-cu)**2
        if s2 > 0 and whi > wlo:
            box_in += rect_disc(Y0, Y1, wlo-cw, whi-cw, mp.sqrt(s2))[0]
    for v in (v0, v1):
        s2 = r*r-(v-cv)**2
        if s2 > 0 and whi > wlo:
            box_in += rect_disc(X0, X1, wlo-cw, whi-cw, mp.sqrt(s2))[0]
    for w in (w0, w1):
        if zlo < w < zhi:
            box_in += rect_disc(X0, X1, Y0, Y1, rho(w))[0]
    out['box_in'] = box_in
    if zone is not None:
        for z in (zlo, zhi):
            if w0 < z < w1:
                out['ends_in'] += rect_disc(X0, X1, Y0, Y1, rho(z))[0]
    return out


def ideal_axes(name):
    """The frame's ideal orthonormal axes (x, y, n), as Frame3::new's rule
    in exact arithmetic."""
    f = FRAMES[name]
    nrm = lambda v: tuple(c/mp.sqrt(sum(x*x for x in v)) for c in v)
    n = nrm(tuple(mp.mpf(c) for c in f[:3]))
    hint = nrm(tuple(mp.mpf(c) for c in f[3:]))
    y = nrm(ref.cross(n, hint))
    x = ref.cross(y, n)
    return x, y, n


def sphere_measures(r, zone=None):
    """Volume, axial moment (about the centre) and area of a sphere or zone
    `[lo, hi]` along its axis."""
    r = mp.mpf(r)
    lo, hi = (-r, r) if zone is None else (mp.mpf(zone[0]), mp.mpf(zone[1]))
    V = mp.pi*((r*r*hi-hi**3/3)-(r*r*lo-lo**3/3))
    Mz = mp.pi*((r*r*hi**2/2-hi**4/4)-(r*r*lo**2/2-lo**4/4))
    A = 2*mp.pi*r*(hi-lo)
    for z in (lo, hi):
        if abs(z) < r:
            A += mp.pi*(r*r-z*z)
    return V, Mz, A


def box_measures(box):
    (u0, u1), (v0, v1), (w0, w1) = [tuple(mp.mpf(x) for x in rng) for rng in box]
    a, b, c = u1-u0, v1-v0, w1-w0
    V = a*b*c
    return V, (V*(u0+u1)/2, V*(v0+v1)/2, V*(w0+w1)/2), 2*(a*b+b*c+c*a)


def zone_of(s):
    """A sphere's zone along its axis from the centre (its stored heights)."""
    _, r, _, low, high = s
    lo = None if low == -HP else r*math.sin(low)
    hi = None if high == HP else r*math.sin(high)
    if lo is None and hi is None:
        return None
    return (-r if lo is None else lo, r if hi is None else hi)


def pair_forms(case):
    """{op: (volume, area, centre)} in closed form for a pair of a sphere
    and a box (or the L of two boxes), in the box's ideal frame: the
    sphere's axis (for a zone) one of the box's axes, permuted to `w`."""
    specs = {}
    for c, role in ((case.obj, 'A'), (case.tool, 'B')):
        specs[role] = c
    S_role = 'A' if case.obj.sphere is not None else 'B'
    P_role = 'B' if S_role == 'A' else 'A'
    sc, pc = specs[S_role], specs[P_role]
    x, y, n = ideal_axes(frame_name(pc.frame))
    po = tuple(mp.mpf(v) for v in pc.frame[:3])
    so = tuple(mp.mpf(v) for v in sc.frame[:3])
    r, low, high = sc.sphere
    zone = zone_of(('sphere', r, sc.frame, low, high))
    axes = [x, y, n]
    perm = [0, 1, 2]
    if zone is not None:
        sn = ideal_axes(frame_name(sc.frame))[2]
        k = max(range(3), key=lambda i: abs(ref.dot(axes[i], sn)))
        assert abs(ref.dot(axes[k], sn)-1) < mp.mpf(10)**-30, 'a zone along one of the box axes'
        perm = [i for i in range(3) if i != k]+[k]
    local = lambda X: tuple(ref.dot(ref.sub(X, po), axes[i]) for i in range(3))
    world_m = lambda V, m: tuple(po[j]*V+sum(m[i]*axes[perm[i]][j] for i in range(3)) for j in range(3))
    cl = local(so)
    c = tuple(cl[i] for i in perm)
    zl = None if zone is None else (c[2]+zone[0], c[2]+zone[1])
    boxes = []
    pts = pc.boundaries[0].points
    if len(pts) == 4:
        boxes.append([(pts[0][0], pts[2][0]), (pts[0][1], pts[2][1]), (pc.start, pc.end)])
    else:
        # The L: [0, 4] x [0, 1] and [0, 1] x [1, 4].
        boxes += [[(0, 4), (0, 1), (pc.start, pc.end)], [(0, 1), (1, 4), (pc.start, pc.end)]]
    Vc, Mc, s_in, p_in, s_tot_extra = mp.mpf(0), (mp.mpf(0),)*3, mp.mpf(0), mp.mpf(0), mp.mpf(0)
    VP, MP, AP = mp.mpf(0), (mp.mpf(0),)*3, mp.mpf(0)
    for b in boxes:
        bp = [b[i] for i in perm]
        f = box_sphere(bp, c, r, zl)
        Vc += f['volume']
        Mc = tuple(p+q for p, q in zip(Mc, f['moments']))
        s_in += f['sphere_in']+f['ends_in']
        p_in += f['box_in']
        bv, bm, ba = box_measures(bp)
        VP += bv
        MP = tuple(p+q for p, q in zip(MP, bm))
        AP += ba
    if len(boxes) == 2:
        # The shared face y = 1, x in [0, 1]: both boxes count it.
        AP -= 2*2
        s2 = mp.mpf(r)**2-(1-c[1])**2
        if s2 > 0:
            p_in -= 2*rect_disc(0-c[0], 1-c[0], mp.mpf(pc.start)-c[2], mp.mpf(pc.end)-c[2], mp.sqrt(s2))[0]
    VS, MzS, AS = sphere_measures(r, zone)
    MS = tuple(c[i]*VS for i in range(2))+(c[2]*VS+MzS,)
    # Coplanar discs (by hand): a hemisphere's disc on the box's face.
    on_same = on_opp = mp.mpf(0)
    if case.coplanar:
        disc = mp.pi*mp.mpf(r)**2
        if case.opposite:
            on_opp = disc
        else:
            on_same = disc
    parts = {'S': (VS, MS, AS, s_in), 'P': (VP, MP, AP, p_in)}
    first, second = ('S', 'P') if S_role == 'A' else ('P', 'S')
    VA, MA, AA, A_in = parts[first]
    VB, MB, AB, B_in = parts[second]
    A_out = AA-A_in-on_same-on_opp
    B_out = AB-B_in-on_same-on_opp
    ops = {'common': (Vc, Mc, A_in+B_in+on_same),
           'fuse': (VA+VB-Vc, tuple(a+b-c_ for a, b, c_ in zip(MA, MB, Mc)), A_out+B_out+on_same),
           'cut': (VA-Vc, tuple(a-c_ for a, c_ in zip(MA, Mc)), A_out+on_opp+B_in)}
    out = {}
    for op, (V, m, A) in ops.items():
        if V <= mp.mpf(10)**-30:
            out[op] = None
            continue
        wm = world_m(V, m)
        out[op] = (V, A, tuple(v/V for v in wm))
    return out


# ------------------------------------------------------------------ Monte Carlo

class FloatSphere:
    def __init__(self, S):
        self.c = [float(v) for v in S.c]
        self.r2 = float(S.r2)
        self.planes = [([float(v) for v in a], float(b)) for a, b, _ in S.planes]
        r = float(S.r)
        self.box = ([v-r for v in self.c], [v+r for v in self.c])

    def contains(self, X):
        if sum((X[i]-self.c[i])**2 for i in range(3)) >= self.r2:
            return False
        return all(a[0]*X[0]+a[1]*X[1]+a[2]*X[2] > b for a, b in self.planes)


class FloatPrism:
    def __init__(self, P):
        self.o = [float(v) for v in P.o]
        inv = [ref.cross(P.y, P.n), ref.cross(P.n, P.x), ref.cross(P.x, P.y)]
        self.inv = [[float(c/P.det) for c in row] for row in inv]
        self.lo, self.hi = float(P.lo), float(P.hi)
        self.ring = [(float(u), float(v)) for u, v in P.ring]
        pts = [[float(c) for c in p] for p in P.P.values()]
        self.box = ([min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)])

    def contains(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        w = sum(self.inv[2][i]*d[i] for i in range(3))
        if not self.lo < w < self.hi:
            return False
        u = sum(self.inv[0][i]*d[i] for i in range(3))
        v = sum(self.inv[1][i]*d[i] for i in range(3))
        inside = False
        n = len(self.ring)
        for k in range(n):
            (x0, y0), (x1, y1) = self.ring[k], self.ring[(k+1) % n]
            if (y0 > v) != (y1 > v):
                xc = x0+(v-y0)*(x1-x0)/(y1-y0)
                if xc > u:
                    inside = not inside
        return inside


def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fs, fp = FloatSphere(pair.S), FloatPrism(pair.P)
    box0 = [min(fs.box[0][i], fp.box[0][i]) for i in range(3)]
    box1 = [max(fs.box[1][i], fp.box[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        s, p = fs.contains(X), fp.contains(X)
        a, b = (s, p) if pair.sphere_first else (p, s)
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
    VS, MS, AS = pair.S.closed()
    VP, MP, AP = pair.P.closed()
    checks = {}
    # The inputs' slicings against their closed forms, and inclusion and
    # exclusion with each operation sliced apart.
    dev = mp.mpf(0)
    for (v, m), (V, Mm) in ((r['S'], (VS, MS)), (r['P'], (VP, MP))):
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inputs_sliced'] = dev
    vc, mc = r['common2']
    dev = abs(r['common'][0]-vc)/size**3
    dev = max(dev, max(abs(x-y) for x, y in zip(r['common'][1], mc))/size**4)
    checks['common_two_ways'] = dev
    dev = mp.mpf(0)
    for key, (V, Mm) in (('fuse', (VS+VP-vc, [a+b-c for a, b, c in zip(MS, MP, mc)])),
                         ('S-P', (VS-vc, [a-c for a, c in zip(MS, mc)])),
                         ('P-S', (VP-vc, [a-c for a, c in zip(MP, mc)]))):
        v, m = r[key]
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inclusion_exclusion'] = dev
    # Faces: classes against exact areas; the area identity; shared areas.
    faces = pair.faces()
    checks['face_classes'] = max(abs(sum(cls.values())-area) for _, _, cls, area in faces)/size**2
    cs, cp = pair.classes('S'), pair.classes('P')
    checks['shared_both_sides'] = max(abs(cs['same']-cp['same']), abs(cs['opp']-cp['opp']))/size**2
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')+2*cs['opp']-AS-AP)/size**2
    checks['input_areas'] = max(abs(sum(sum(cls.values()) for w, _, cls, _ in faces if w == k)-A)
                                for k, A in (('S', AS), ('P', AP)))/size**2
    if not first.coplanar:
        second = ref.Slicing(pair.S, pair.P, SECOND, size).results()
        dev = mp.mpf(0)
        for key in ('S', 'P', 'common', 'fuse', 'S-P', 'P-S'):
            dev = max(dev, abs(second[key][0]-r[key][0])/size**3,
                      max(abs(x-y) for x, y in zip(second[key][1], r[key][1]))/size**4)
        for k in ('in', 'out'):
            dev = max(dev, abs(second['sphere_face'][k]-r['sphere_face'][k])/size**2)
        checks['second_direction'] = dev
    mc = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mc[op]
        V, m = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in m]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(pair)
    stats = {'volume_quadrature': pair.slicing.quad_error/size**4, 'breaks': len(pair.slicing.breaks)}
    return name, rows, res, checks, near, stats


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
    b = pair.slicing.raw_breaks
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
    return first.pair_name, pair_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
          'inputs_sliced': 1e-30, 'common_two_ways': 1e-30, 'inclusion_exclusion': 1e-30,
          'face_classes': 1e-30, 'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30,
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
        _, _, res, _, _, _ = by_pair[name]
        kind = 'exact' if exact_pair(firsts[name]) else 'turned'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}_frames', dev)
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, sphere_boolean_reference.py: expect KIND {STEP}, reason TEXT for a degenerate '
           'case, then result N volume area cx cy cz or empty)']
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


def jobs(mc_n):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items()]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert (c.obj.sphere is None) != (c.tool.sphere is None), f'{c.name}: a sphere and a prism'
        s = c.obj if c.obj.sphere is not None else c.tool
        if s.sphere[1] != -HP or s.sphere[2] != HP:
            assert frame_name(s.frame) in EXACT and s.frame[:3] == (0.0, 0.0, 0.0), \
                f'{c.name}: a cap or zone in an exact frame at the origin'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
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
    by_pair = {r[0]: r for r in results}
    for name, _, _, _, n, stats in results:
        if n:
            print('declared near coincidences:', name, '; '.join(n[:4]))
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
