#!/usr/bin/env python3
"""Fixtures for S9f.3a of REVIEW_NOTES.md: Booleans of spline prisms against
spheres and hemispheres in any relative position, either the object: a
spline wall meets the sphere in a curve that is, along the wall's ruling at
the spline's parameter `t`, a root of `A w^2 + 2 B(t) w + C(t)` (`A = n .
n`, degrees `p` and `2 p`), a graph over `t` between its turning points
(the roots of `B^2 - A C`) or, about a turning point inside both faces, a
graph over the height (S9f.2b.2's loops).

`boolean-spline-sphere-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`; a sphere's block its `frame`
and `sphere R LOW HIGH` rows). The profiles are S9a.2's and S9f.1's
(`generate_spline_any_boolean_fixtures.profiles`: the bulge, dome, blob,
wave, lens and `knot`, R4's C1 knot of multiplicity two), in `XY` and in
`TILT` (`generate_curved_boolean_fixtures.FRAMES`); the spheres whole in
`XY` about binary64 centres (a whole sphere's model does not depend on its
frame's axes but through the kernel's split, which is its own), and
hemispheres (latitudes 0 and pi/2, whose heights are exact) in `XY` and on
their `SIDE`. `boolean-spline-sphere-frames.tsv` records the stored axes.

`boolean-spline-sphere-expected.tsv` gives per case, from
`spline_sphere_boolean_reference.py`:

* `expect KIND S9f.3a`: the declared outcome (`solid`, `empty`, or
  `degenerate`: the sphere tangent to the wall, a turning point at a knot
  or on a face's boundary, `Degenerate` in the decisions); then for a
  degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

Pairs: a sphere straddling the bulge's wall (`bulge_ball`, a loop with two
turning points inside both faces), the lens's lower cubic (`lens_ball`),
the blob's wall and its top cap (`blob_cap_ball`, a loop cut by the cap)
and `knot`'s second span in `TILT` (`knot_tilt_ball`); a sphere as the
object across the wave's knot (`ball_wave`); a large sphere over the dome
(`dome_ball`, its meeting a graph over `t` from edge to edge, the top cap
inside it); a sphere below the blob in `TILT` (`tilt_blob_low`, arches
from the bottom cap's edge back to it); a hemisphere through the bulge's
bottom cap on its axis (`hemi_bulge`: its split great circle's crossing of
the wall in a tower field) and one on its side whose rim's plane holds
the wall's axis (`side_hemi_bulge`: the rim's crossings in a tower field,
one turning point inside the hemisphere and one beyond its rim); and
`degenerate` with reasons: a sphere touching the dome's arch at its apex
(`dome_touch`), a turning point at `knot`'s knot (`knot_turn`), a loop
turning back on the bulge's top cap's edge (`cap_turn`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* each region's (the prism's section, the sphere's, their common) volume
  and first moments by two slicings, along the caps' normal and along
  `(2, -3, 5)`, with their own breakpoints and sections, within 1e-30;
* both inputs' slicings against their closed forms; every face's classes
  summing to its area (the caps' profile and the walls' closed forms, the
  sphere's `4 pi r^2` or `2 pi r^2` by Archimedes, the disc's `pi r^2`);
  the area identity `area(fuse) + area(common) = area(A) + area(B)` against
  the closed forms; within 1e-30;
* Monte Carlo estimates (100,000 points per pair, a fixed seed, the
  profile as a polygon of 512 points a span) of every operation's volume
  and centre within 5 standard errors;
* solid counts: each case's declared count, the reference's by its
  slices' union-find;
* margins that flag near coincidences: classes or volumes positive but
  below 1e-9 of the size, breakpoints closer than that; and the
  meeting's own: every turning point on a spline wall at least 1e-3
  outside a face (`turn`) or inside both by at least that (`loop`), at
  least 1e-3 from a knot or a vertex of the profile (`knot`), the
  discriminant's critical values at least 1e-3 from zero in the gap they
  measure (`touch`), the meeting's points at a cap or the rim crossing it
  at a sine of at least 1e-3 and the surfaces there at a sine of at least
  1e-3 (`cross`), the rim's crossings of the walls and caps at a sine of
  at least 1e-3 (`rim`), a hemisphere's plane at least 1e-3 from the caps
  where parallel and its trace crossing the spans at a sine of at least
  1e-3 where it holds the axis (`trace`), the prism's vertical edges
  crossing the sphere at a sine of at least 1e-3 (`pierce`), every vertex
  of either input (the kernel's poles of a whole sphere for each of its
  four splits) at least 1e-3 from the other's faces (`vertex`). The
  declared degenerate pairs must fail their margin (below 1e-12) and are
  exempt from the rest.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import math
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Case, encode_boolean_case
from curve_surface_reference import stored_axes
import spline_sphere_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES
import generate_spline_any_boolean_fixtures as s9f1

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-spline-sphere'
STEP = 'S9f.3a'
BOOLEAN_OPERATION = 93
HP = ref.HALF_PI

TOUCH = 'a sphere tangent to a spline wall'
KNOT_TURN = "a spline wall's meeting with a sphere turning back at a knot"
TURN_EDGE = "a spline wall's meeting with a sphere turning back on a face's boundary"

dot, cross, sub, add, scale, M, Mv = ref.dot, ref.cross, ref.sub, ref.add, ref.scale, ref.M, ref.Mv


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def world(name, origin, local):
    """A frame's point at local coordinates, in binary64 (the stored axes:
    a sphere's centre on a prism in a turned frame)."""
    o, x, y, n = stored_axes(at(name, origin))
    return tuple(float(o[i])+local[0]*x[i]+local[1]*y[i]+local[2]*n[i] for i in range(3))


def profiles():
    return s9f1.profiles()


def ball(r, frame, low=-HP, high=HP):
    return ('sphere', r, frame, low, high)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def make(spec, op):
    if spec[0] == 'sphere':
        _, r, frame, low, high = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), low, high))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


class Boolean:
    def __init__(self, name, operation, obj, tool, kind='solid', solids=1, reason=None):
        self.name, self.operation, self.kind, self.solids, self.reason = name, operation, kind, solids, reason
        self.step = STEP
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops):
    """Cases of one pair: `ops` maps each operation to its solid count, or to
    `('empty',)`, or to `('degenerate', reason)`."""
    out = []
    for op, want in ops.items():
        if isinstance(want, int):
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'solid', want))
        elif want[0] == 'empty':
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'empty', 0))
        else:
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'degenerate', None, want[1]))
    return out


ONE3 = {'fuse': 1, 'cut': 1, 'common': 1}


def cases():
    p = profiles()
    O = (0, 0, 0)
    out = []
    # Loops: spheres straddling a spline wall inside the heights.
    out += group('bulge_ball', prism(p['bulge'], at('XY', O), 0.0, 3.0), ball(1.0, at('XY', (11.125, 3.0, 1.5))),
                 ONE3)
    out += group('lens_ball', prism(p['lens'], at('XY', O), 0.0, 4.0), ball(0.8, at('XY', (5.0, 2.875, 2.0))),
                 ONE3)
    out += group('blob_cap_ball', prism(p['blob'], at('XY', O), 0.0, 3.0), ball(1.0, at('XY', (8.125, 4.0, 2.625))),
                 ONE3)
    out += group('knot_tilt_ball', prism(p['knot'], at('TILT', O), 0.0, 3.0),
                 ball(0.75, at('XY', world('TILT', O, (2.75, 5.125, 1.5)))), ONE3)
    out += group('ball_wave', ball(1.0, at('XY', (6.5, 6.625, 1.5))), prism(p['wave'], at('XY', O), 0.0, 3.0), ONE3)
    # Branches over the spline's parameter: turning points outside the faces.
    out += group('dome_ball', prism(p['dome'], at('XY', O), 0.0, 3.0), ball(3.0, at('XY', (2.0, 0.5, 4.25))),
                 ONE3)
    out += group('tilt_blob_low', prism(p['blob'], at('TILT', O), 0.0, 3.0),
                 ball(2.0, at('XY', world('TILT', O, (7.875, 4.25, -1.0)))), ONE3)
    # Hemispheres: the split great circle's tower on the prism's axis, the
    # rim's on its side.
    out += group('hemi_bulge', prism(p['bulge'], at('XY', O), 0.0, 3.0),
                 ball(1.5, at('XY', (11.25, 3.0, -0.5)), 0.0, HP), ONE3)
    out += group('side_hemi_bulge', prism(p['bulge'], at('XY', O), 0.0, 3.0),
                 ball(0.9, at('SIDE', (10.6, 4.6, 1.5)), 0.0, HP), ONE3)
    # Declared degenerate.
    out += group('dome_touch', prism(p['dome'], at('XY', O), 0.0, 3.0), ball(1.0, at('XY', (2.0, 3.0, 1.5))),
                 {op: ('degenerate', TOUCH) for op in ('cut', 'common')})
    out += group('knot_turn', prism(p['knot'], at('XY', O), 0.0, 3.0), ball(2.5, at('XY', (6.5, 6.0, 1.5))),
                 {op: ('degenerate', KNOT_TURN) for op in ('cut', 'common')})
    out += group('cap_turn', prism(p['bulge'], at('XY', O), 0.0, 3.0), ball(1.0, at('XY', (11.0, 3.0, 3.0))),
                 {op: ('degenerate', TURN_EDGE) for op in ('cut', 'common')})
    return out


# ------------------------------------------------------------------ Monte Carlo

class FloatProfile:
    """The profile as a polygon (512 points a span), its edges bucketed by
    `v` for a crossing count."""

    def __init__(self, S, per=512):
        edges = []
        for el in S.elements:
            n = per if el.spline else 1
            pts = [tuple(float(c) for c in el.point(mp.mpf(k)/n)) for k in range(n+1)]
            edges += list(zip(pts, pts[1:]))
        self.v0 = min(min(a[1], b[1]) for a, b in edges)
        self.v1 = max(max(a[1], b[1]) for a, b in edges)
        self.n = 1024
        self.h = (self.v1-self.v0)/self.n
        self.buckets = [[] for _ in range(self.n)]
        for a, b in edges:
            lo, hi = sorted((a[1], b[1]))
            for k in range(max(0, int((lo-self.v0)/self.h)), min(self.n, int((hi-self.v0)/self.h)+1)):
                self.buckets[k].append((a, b))

    def inside(self, u, v):
        if not self.v0 < v < self.v1:
            return False
        k = min(self.n-1, int((v-self.v0)/self.h))
        count = 0
        for (x0, y0), (x1, y1) in self.buckets[k]:
            if (y0 <= v < y1) or (y1 <= v < y0):
                if u < x0+(v-y0)*(x1-x0)/(y1-y0):
                    count += 1
        return count % 2 == 1


def monte_carlo(pair, name, n=100000):
    """Each operation's volume and centre by `n` uniform points in the
    pair's box: (volume, centre, their standard errors)."""
    S, B = pair.S, pair.ball
    fp = FloatProfile(S)
    inv = [[float(c) for c in r] for r in S.inv]
    o = [float(c) for c in S.o]
    lo, hi = float(S.lo), float(S.hi)
    c = [float(x) for x in B.c]
    r2 = float(B.r)**2
    ns = [float(x) for x in B.ns] if B.half else None
    box = ([float(x) for x in pair.lo], [float(x) for x in pair.hi])
    vol = 1.0
    for i in range(3):
        vol *= box[1][i]-box[0][i]
    rng = random.Random(zlib.crc32(name.encode()))
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box[0][i]+(box[1][i]-box[0][i])*rng.random() for i in range(3)]
        d = [X[i]-o[i] for i in range(3)]
        u, v, w = (sum(inv[k][i]*d[i] for i in range(3)) for k in range(3))
        inP = lo < w < hi and fp.inside(u, v)
        e = [X[i]-c[i] for i in range(3)]
        inB = e[0]*e[0]+e[1]*e[1]+e[2]*e[2] < r2 and (ns is None or sum(a*b for a, b in zip(e, ns)) > 0)
        obj, tool = (inP, inB) if pair.obj_prism else (inB, inP)
        for op, inside in (('fuse', obj or tool), ('cut', obj and not tool), ('common', obj and tool)):
            if inside:
                a = acc[op]
                a[0] += 1
                for i in range(3):
                    a[1][i] += X[i]
                    a[2][i] += X[i]*X[i]
    out = {}
    for op, (k, s1, s2) in acc.items():
        p = k/n
        v = vol*p
        se = vol*math.sqrt(max(p*(1-p), 1e-300)/n)
        if k > 1:
            mean = [x/k for x in s1]
            sd = [math.sqrt(max(s2[i]/k-mean[i]**2, 0.0)/k) for i in range(3)]
        else:
            mean, sd = None, None
        out[op] = (v, se, mean, sd)
    return out


# ------------------------------------------------------------------ margins

INF = mp.mpf('inf')
ROTATIONS = [((2, 3, 6), 7), ((1, 4, 8), 9), ((2, 6, 9), 11), ((1, 2, 2), 3)]


def sine(a, b):
    c = cross(a, b)
    return mp.sqrt(dot(c, c)/(dot(a, a)*dot(b, b)))


def wall_normal(S, el, t):
    du, dv = el.deriv(t)
    T = tuple(du*S.xm[i]+dv*S.ym[i] for i in range(3))
    return cross(T, S.nm)


def hemi_depth(ball, X):
    """How far a point of the sphere lies on the hemisphere's side (signed)."""
    if not ball.half:
        return INF
    return dot(Mv(ball.ns), sub(X, ball.cm))


def prism_distance(S, X):
    """The distance from a point to the prism's boundary (binary64-grade)."""
    u, v, w = S.local(X)
    lo, hi = M(S.lo), M(S.hi)
    bd = S.profile.boundary_distance(u, v)
    inside = S.profile.inside(u, v)
    if lo <= w <= hi:
        return min(w-lo, hi-w, bd) if inside else bd
    dh = lo-w if w < lo else w-hi
    return dh if inside else mp.sqrt(dh*dh+bd*bd)


def ball_distance(ball, X):
    d = sub(X, ball.cm)
    rad = mp.sqrt(dot(d, d))
    if not ball.half:
        return abs(rad-ball.rm)
    h = dot(Mv(ball.ns), d)
    if h >= 0:
        return min(abs(rad-ball.rm), h if rad <= ball.rm else INF)
    flat = mp.sqrt(max(rad*rad-h*h, 0))
    return -h if flat <= ball.rm else mp.sqrt(h*h+(flat-ball.rm)**2)


def margins(pair):
    S, ball = pair.S, pair.ball
    out = {k: INF for k in ('turn', 'loop', 'knot', 'touch', 'cross', 'rim', 'trace', 'pierce', 'vertex')}
    lo, hi = M(S.lo), M(S.hi)
    A = M(dot(S.n, S.n))
    m = Mv(S.m)
    for el in S.elements:
        P, dP = S.wall(el)
        Aq, B, C = ref.wall_quadratic(S, el, ball)
        D = ref.psub(ref.pmul(B, B), ref.pc(C, Aq))
        Bm, Cm, Dm = ref.pm(B), ref.pm(C), ref.pm(D)
        if el.spline:
            for t in ref.exact_roots(D, slack=mp.mpf(0)):
                w = -ref.peval(Bm, t)/A
                u, v = el.point(t)
                X = S.world(u, v, w)
                hd = hemi_depth(ball, X)
                outside = max(0, lo-w, w-hi, -hd if hd != INF else 0)
                if outside == 0:
                    out['loop'] = min(out['loop'], w-lo, hi-w, hd)
                else:
                    out['turn'] = min(out['turn'], outside)
                out['knot'] = min(out['knot'], min(mp.sqrt((u-k[0])**2+(v-k[1])**2) for k in S.profile.vm))
        for t in ref.exact_roots(ref.pder(D), slack=mp.mpf(0)):
            out['touch'] = min(out['touch'], 2*mp.sqrt(abs(ref.peval(Dm, t)))/A)
        # The meeting at the caps.
        for wq in (S.lo, S.hi):
            for t in ref.exact_roots(ref.padd(ref.padd([Aq*wq*wq], ref.pc(B, 2*wq)), C), slack=mp.mpf(0)):
                u, v = el.point(t)
                X = S.world(u, v, M(wq))
                if hemi_depth(ball, X) < 0:
                    continue
                Nw, Ns = wall_normal(S, el, t), sub(X, ball.cm)
                out['cross'] = min(out['cross'], sine(Nw, Ns), 1-sine(cross(Nw, Ns), m))
        if ball.half:
            ns = Mv(ball.ns)
            nsn = dot(ball.ns, S.n)
            hp = ref.vpoly_dot(ball.ns, tuple(ref.psub(P[i], [ball.c[i]]) for i in range(3)))
            pts = []
            if nsn != 0:
                wh = ref.pc(hp, -1/nsn)
                for t in ref.exact_roots(ref.padd(ref.padd(ref.pc(ref.pmul(wh, wh), Aq), ref.pc(ref.pmul(B, wh), 2)),
                                                  C), slack=mp.mpf(0)):
                    pts.append((t, ref.peval(ref.pm(wh), t)))
                whm = ref.pm(wh)
                out['trace'] = min(out['trace'], min(abs(ref.peval(whm, t)-k) for t in (mp.mpf(0), mp.mpf(1))
                                                     for k in (lo, hi)))
            else:
                for t in ref.exact_roots(hp, slack=mp.mpf(0)):
                    du, dv = el.deriv(t)
                    T = tuple(du*S.xm[i]+dv*S.ym[i] for i in range(3))
                    out['trace'] = min(out['trace'], sine(T, cross(ns, S.nm)))
                    b, cc = ref.peval(Bm, t), ref.peval(Cm, t)
                    disc = b*b-A*cc
                    if disc > 0:
                        pts += [(t, (-b+sg*mp.sqrt(disc))/A) for sg in (-1, 1)]
            for t, w in pts:
                if not lo <= w <= hi:
                    continue
                u, v = el.point(t)
                X = S.world(u, v, w)
                Nw, Ns = wall_normal(S, el, t), sub(X, ball.cm)
                Tr = cross(ns, Ns)
                out['rim'] = min(out['rim'], abs(dot(Tr, Nw))/mp.sqrt(dot(Tr, Tr)*dot(Nw, Nw)))
                out['cross'] = min(out['cross'], sine(Nw, Ns), 1-sine(cross(Nw, Ns), ns))
    if ball.half:
        ns = Mv(ball.ns)
        for wq in (S.lo, S.hi):
            for X in ref.plane_circle(S.m, add(S.o, scale(S.n, wq)), ball):
                u, v, _ = S.local(X)
                if S.profile.inside(u, v):
                    Tr = cross(ns, sub(X, ball.cm))
                    out['rim'] = min(out['rim'], abs(dot(Tr, m))/mp.sqrt(dot(Tr, Tr)*dot(m, m)))
    # Vertical edges against the sphere.
    for vtx in S.profile.vm:
        P0 = S.world(vtx[0], vtx[1], mp.mpf(0))
        Q = sub(P0, ball.cm)
        b, cc = dot(S.nm, Q), dot(Q, Q)-ball.rm**2
        disc = b*b-A*cc
        if disc > 0:
            for sg in (-1, 1):
                w = (-b+sg*mp.sqrt(disc))/A
                X = add(P0, scale(S.nm, w))
                if lo <= w <= hi and hemi_depth(ball, X) >= 0:
                    out['pierce'] = min(out['pierce'], 1-sine(S.nm, sub(X, ball.cm)))
    # Vertices against the other's faces.
    for vtx in S.profile.vm:
        for w in (lo, hi):
            out['vertex'] = min(out['vertex'], ball_distance(ball, S.world(vtx[0], vtx[1], w)))
    poles = []
    if ball.half:
        poles.append(add(ball.cm, scale(Mv(ball.ns), ball.rm)))
    else:
        for row, length in ROTATIONS:
            for sg in (-1, 1):
                poles.append(add(ball.cm, scale(Mv(row), sg*ball.rm/length)))
    for X in poles:
        out['vertex'] = min(out['vertex'], prism_distance(S, X))
    return out


MARGIN = mp.mpf(10)**-3
DEGENERATE_MARGIN = {TOUCH: 'touch', KNOT_TURN: 'knot', TURN_EDGE: 'loop'}


def near_coincidences(pair):
    size = pair.size
    out = []
    for tag, name, cls, _ in pair.faces():
        for c, v in cls.items():
            if 0 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{tag} {name} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volumes()[op][0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    for way in (pair.first, pair.second):
        b = pair.sliced(way)['breaks']
        for p, q in zip(b, b[1:]):
            if q-p < mp.mpf(10)**-9*size:
                out.append(f'breakpoints {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, obj, tool, ops = job
    pair = ref.Pair(obj, tool)
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    size = pair.size
    checks = {}
    one, two = pair.sliced(pair.first), pair.sliced(pair.second)
    dev = mp.mpf(0)
    for k in 'PBC':
        dev = max(dev, abs(one[k][0]-two[k][0])/size**3,
                  max(abs(x-y) for x, y in zip(one[k][1], two[k][1]))/size**4)
    checks['second_slicing'] = dev
    vs, ms, as_ = pair.S.closed()
    vb, mb, ab = pair.ball.closed()
    checks['closed_prism'] = max(abs(one['P'][0]-vs)/size**3,
                                 max(abs(x-y) for x, y in zip(one['P'][1], ms))/size**4)
    checks['closed_sphere'] = max(abs(one['B'][0]-vb)/size**3,
                                  max(abs(x-y) for x, y in zip(one['B'][1], mb))/size**4)
    checks['face_classes'] = max(abs(cls['in']+cls['out']-closed) for _, _, cls, closed in pair.faces())/size**2
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')-as_-ab)/size**2
    mc = monte_carlo(pair, name)
    z = 0.0
    vols = pair.volumes()
    for op in ref.OPS:
        v, se, mean, sd = mc[op]
        V, mom = vols[op]
        z = max(z, abs(v-float(V))/se)
        if mean is not None and V > 0:
            for i in range(3):
                if sd[i] > 0:
                    z = max(z, abs(mean[i]-float(mom[i]/V))/sd[i])
    near = near_coincidences(pair)
    stats = {'volume_quadrature': max(one['error'], two['error'])/size**4,
             'breaks': (len(one['breaks']), len(two['breaks']))}
    return name, rows, checks, z, near, margins(pair), stats


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


def generate(results, listed):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, spline_sphere_boolean_reference.py: expect KIND {STEP}, reason TEXT for a '
           'degenerate case, then result N volume area cx cy cz or empty)']
    for case in listed:
        blocks.append(case.encode())
        rows = by_pair[case.pair_name][1]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        if case.kind == 'solid':
            assert n == case.solids, f'{case.name}: declared {case.solids} solids, the reference gives {row[0]}'
        elif case.kind == 'empty':
            assert n == 0, f'{case.name}: declared empty, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (obj, tool) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in listed:
        for label, c in (('obj', case.obj), ('tool', case.tool)):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\t{label} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(listed):
    pairs = {}
    for c in listed:
        pairs.setdefault(c.pair_name, [c.obj, c.tool, []])[2].append(c.operation)
    return [(name, obj, tool, ops) for name, (obj, tool, ops) in pairs.items()]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        spheres = [x for x in (c.obj, c.tool) if x.sphere is not None]
        assert len(spheres) == 1, f'{c.name}: one sphere'
        prism_case = c.tool if c.obj.sphere is not None else c.obj
        assert s9f1_has_spline(prism_case), f'{c.name}: a spline prism'
        frame_name(c.obj.frame), frame_name(c.tool.frame)


def s9f1_has_spline(case):
    return any(b.segments and any(isinstance(s, ref.cref.Spline) for s in b.segments) for b in case.boundaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--only', help='evaluate the pairs whose name holds this text, write nothing')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    if args.only:
        listed = [c for c in listed if args.only in c.name]
    results = run(jobs(listed), evaluate, args.workers)
    declared = {c.pair_name: c for c in listed}
    worst_margin, degenerate_margin = {}, {}
    failed = []
    worst_z = 0.0
    for name, rows, checks, z, near, margin, stats in results:
        first = declared[name]
        kinds = {c.kind for c in listed if c.pair_name == name}
        if args.only:
            print(name, {op: r[0] for op, r in rows.items()})
            print('  checks', {k: mp.nstr(v, 3) for k, v in checks.items()}, 'monte carlo z', round(z, 2))
            print('  margins', {k: mp.nstr(v, 3) for k, v in margin.items()}, 'stats',
                  {k: (mp.nstr(v, 3) if isinstance(v, mp.mpf) else v) for k, v in stats.items()})
            if near:
                print('  near', near)
        if 'degenerate' in kinds:
            key = DEGENERATE_MARGIN[first.reason]
            if not margin[key] < mp.mpf(10)**-12:
                failed.append(f'{name}: declared degenerate, {key} margin {margin[key]}')
            degenerate_margin[key] = max(degenerate_margin.get(key, mp.mpf(0)), margin[key])
            continue
        worst_z = max(worst_z, z)
        if z > 5:
            failed.append(f'{name}: Monte Carlo {z:.2f} standard errors')
        if near:
            failed.append(f'{name}: near coincidences {near}')
        for key, value in margin.items():
            if value < MARGIN:
                failed.append(f'{name}: {key} margin {mp.nstr(value, 3)}')
            worst_margin[key] = min(worst_margin.get(key, INF), value)
    if failed:
        raise SystemExit('\n'.join(failed))
    if args.only:
        return
    worst = {}
    for name, rows, checks, z, near, margin, stats in results:
        for key, value in checks.items():
            worst[key] = max(worst.get(key, mp.mpf(0)), value)
        worst['volume_quadrature_estimate'] = max(worst.get('volume_quadrature_estimate', mp.mpf(0)),
                                                  stats['volume_quadrature'])
    for key, value in worst.items():
        assert value <= 1e-30, (key, mp.nstr(value, 3))
    files = generate(results, cases())
    for fname, contents in files.items():
        target = ROOT/'fixtures'/fname
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
    print(f'Monte Carlo: largest deviation {worst_z:.2f} standard errors')
    print('smallest margins (non-degenerate pairs):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst_margin.items())))
    print('declared degenerate margins (largest):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(degenerate_margin.items())))


if __name__ == '__main__':
    main()
