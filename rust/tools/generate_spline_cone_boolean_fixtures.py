#!/usr/bin/env python3
"""Fixtures for S9f.3b of REVIEW_NOTES.md: Booleans of spline prisms against
cones and frustums in any relative position, either the object: a spline
wall meets the cone in a curve that is, along the wall's ruling at the
spline's parameter `t`, a root of `A w^2 + 2 B(t) w + C(t)` (`A = q_u^2 +
q_v^2 - k^2 q_w^2` of either sign, degrees `p` and `2 p`): for `A > 0` a
graph over `t` between its turning points (the roots of `B^2 - A C`) or,
about a turning point inside both faces, a graph over the height (S9f.2b.2's
loops); for `A < 0` two branches over the whole run, one on each nappe.

`boolean-spline-cone-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`; a cone's block its `frame` and
`cone BOTTOM TOP HEIGHT` rows). The profiles are S9a.2's and S9f.1's
(`generate_spline_any_boolean_fixtures.profiles`), in `XY` and in `TILT`
(`generate_curved_boolean_fixtures.FRAMES`); the cones in exact frames
(axes along the world's: `XY`, `DOWN`, `SIDE` and `YAX`, the last's axis
along `y`) about binary64 points. `boolean-spline-cone-frames.tsv` records
the stored axes.

`boolean-spline-cone-expected.tsv` gives per case, from
`spline_cone_boolean_reference.py`:

* `expect KIND S9f.3b`: the declared outcome (`solid`, `empty`, or
  `degenerate`: the apex on a spline wall, a wall along the cone's ruling,
  the cone tangent to a wall, `Degenerate` in the decisions); then for a
  degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

Pairs: `A < 0`: a frustum on the bulge's axis across its wall and both caps
(`bulge_frustum`, the other nappe far above), a cone pointing down into the
lens prism with its apex inside it (`lens_apex`: the other nappe above the
apex meets the lens's walls inside the prism's heights, outside the cone's
face), a cone hanging from above the blob with its apex below the prism
(`blob_down`), a frustum as the object across the wave's middle span
(`cone_wave`); `A > 0`: a thin frustum across the dome's axis piercing its
arch (`dome_pierce`, a loop with two turning points inside both faces), a
frustum across the dome's axis whose bottom rim's plane holds the axis and
crosses the arch's generatrices (`dome_side_cone`: the rim's points in a
tower field, branches over the run from them over the arch, no turning
point on the quadric, the prism's top cap cutting the cone in a
hyperbola), and a frustum on `z` against `knot` in `TILT`
(`knot_tilt_cone`, a loop); and
`degenerate` with reasons: a cone's apex on the dome's arch
(`apex_wall`), the bulge in `TILT` against a cone whose slope is the
stored axis's (`ruling_tilt`, `A = 0`), and a cone across the dome's axis
touching its arch (`dome_touch`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* each region's (the prism's section, the cone's, their common) volume and
  first moments by two slicings with their own breakpoints and sections,
  within 1e-30;
* both inputs' slicings against their closed forms; every face's classes
  summing to its area (the caps' profile by chords against its closed form,
  the walls', the cone's wall by its area element against `pi (b + t)`
  times the slant, the discs' `pi r^2`); the area identity `area(fuse) +
  area(common) = area(A) + area(B)` against the closed forms; within 1e-30;
* Monte Carlo estimates (100,000 points per pair, a fixed seed, the
  profile as a polygon of 512 points a span) of every operation's volume
  and centre within 5 standard errors;
* solid counts: each case's declared count, the reference's by its
  slices' union-find;
* margins that flag near coincidences: classes or volumes positive but
  below 1e-9 of the size, breakpoints closer than that; and the
  meeting's own: `|A|` at least 1e-3 (`axis`), the apex (real or virtual)
  at least 1e-3 from every spline wall's surface (`apex`), every turning
  point on a spline wall at least 1e-3 outside a face (`turn`) or inside
  both by at least that (`loop`), at least 1e-3 from a knot or a vertex of
  the profile (`knot`), the discriminant's critical values at least 1e-3
  from zero in the gap they measure (`touch`), the meeting's points at a
  cap crossing it at a sine of at least 1e-3 and the surfaces there at a
  sine of at least 1e-3 (`cross`), the rims' crossings of the walls and
  caps at a sine of at least 1e-3 (`rim`), an end plane at least 1e-3 from
  the caps where parallel and its trace crossing the spans at a sine of at
  least 1e-3 where it holds the axis (`trace`), the prism's vertical edges
  crossing the cone at a sine of at least 1e-3 (`pierce`), every vertex of
  either input (the cone's apex and its rims' seam points) at least 1e-3
  from the other's faces (`vertex`). The declared degenerate pairs must
  fail their margin (below 1e-12) and are exempt from the rest.

Pairs run in worker processes (`--workers`).
"""
import argparse
from fractions import Fraction as F
import math
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Case
from curve_surface_reference import stored_axes
import spline_cone_boolean_reference as ref
import generate_spline_sphere_boolean_fixtures as s9f3a
from generate_curved_boolean_fixtures import FRAMES
import generate_spline_any_boolean_fixtures as s9f1

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-spline-cone'
STEP = 'S9f.3b'
BOOLEAN_OPERATION = 93

APEX = "a cone's apex on the other input's surface"
RULING = "a spline wall along a cone's ruling"
TOUCH = 'a cone tangent to a spline wall'

# Exact frames for the cones beside the shared ones: `YAX` along `y`.
CONE_FRAMES = {'XY': FRAMES['XY'], 'DOWN': FRAMES['DOWN'], 'SIDE': FRAMES['SIDE'],
               'YAX': (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)}
ALL_FRAMES = dict(FRAMES, **CONE_FRAMES)
# The kernel's first seams (`curved::SEAMS`): the object's rims at the
# first, the tool's at the second.
SEAMS = (F(2, 7), F(3, 11))

dot, cross, sub, add, scale, M, Mv = ref.dot, ref.cross, ref.sub, ref.add, ref.scale, ref.M, ref.Mv


def at(name, origin):
    return tuple(float(c) for c in origin)+ALL_FRAMES[name]


def world(name, origin, local):
    o, x, y, n = stored_axes(at(name, origin))
    return tuple(float(o[i])+local[0]*x[i]+local[1]*y[i]+local[2]*n[i] for i in range(3))


def profiles():
    return s9f1.profiles()


def cone(frame, bottom, top, height):
    return ('cone', frame, bottom, top, height)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def make(spec, op):
    if spec[0] == 'cone':
        _, frame, bottom, top, height = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(bottom), float(top), float(height)))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


class Boolean(s9f3a.Boolean):
    def __init__(self, name, operation, obj, tool, kind='solid', solids=1, reason=None):
        self.name, self.operation, self.kind, self.solids, self.reason = name, operation, kind, solids, reason
        self.step = STEP
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name


def group(name, obj, tool, ops):
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


def ruling_cone():
    """A cone on `z` whose slope is the `TILT` frame's stored axis's ratio
    (`k = n_y / n_z`, exactly: `A = 0`), its apex at its base."""
    _, _, _, n = stored_axes(at('TILT', (0, 0, 0)))
    return float(2*n[1]), float(2*n[2])


def cases():
    p = profiles()
    O = (0, 0, 0)
    out = []
    # A < 0: every ruling meets both nappes.
    out += group('bulge_frustum', prism(p['bulge'], at('XY', O), 0.0, 3.0),
                 cone(at('XY', (10.75, 3.0, -0.5)), 1.0, 0.5, 4.0), ONE3)
    out += group('lens_apex', prism(p['lens'], at('XY', O), 0.0, 4.0),
                 cone(at('XY', (5.0, 4.0, -1.0)), 3.0, 0.0, 3.0), ONE3)
    out += group('blob_down', prism(p['blob'], at('XY', O), 0.0, 3.0),
                 cone(at('DOWN', (7.0, 4.0, 4.0)), 2.5, 0.0, 5.5), ONE3)
    out += group('cone_wave', cone(at('XY', (5.0, 5.25, 0.5)), 1.5, 0.75, 2.0),
                 prism(p['wave'], at('XY', O), 0.0, 3.0), ONE3)
    # A > 0: turning points.
    out += group('dome_pierce', prism(p['dome'], at('XY', O), 0.0, 3.0),
                 cone(at('YAX', (3.0, 0.5, 1.5)), 0.5, 0.25, 3.5), ONE3)
    out += group('dome_side_cone', prism(p['dome'], at('XY', O), 0.0, 3.0),
                 cone(at('YAX', (2.125, 1.25, 2.75)), 1.5, 0.375, 1.25), ONE3)
    out += group('knot_tilt_cone', prism(p['knot'], at('TILT', O), 0.0, 3.0),
                 cone(at('XY', world('TILT', O, (5.0, 4.0, 1.5))[:2]+(-2.5,)), 0.5, 1.5, 4.0), ONE3)
    # Declared degenerate.
    out += group('apex_wall', prism(p['dome'], at('XY', O), 0.0, 3.0), cone(at('XY', (2.0, 2.0, 1.5)), 0.0, 1.0, 1.0),
                 {op: ('degenerate', APEX) for op in ('cut', 'common')})
    top, height = ruling_cone()
    out += group('ruling_tilt', prism(p['bulge'], at('TILT', O), 0.0, 3.0),
                 cone(at('XY', world('TILT', O, (10.5, 3.0, -1.0))), 0.0, top, height),
                 {op: ('degenerate', RULING) for op in ('cut', 'common')})
    out += group('dome_touch', prism(p['dome'], at('XY', O), 0.0, 3.0),
                 cone(at('SIDE', (-0.5, 0.875, 1.5)), 0.0, 1.25, 2.5),
                 {op: ('degenerate', TOUCH) for op in ('cut', 'common')})
    return out


# ------------------------------------------------------------------ Monte Carlo

def monte_carlo(pair, name, n=100000):
    """Each operation's volume and centre by `n` uniform points in the
    pair's box: (volume, centre, their standard errors)."""
    S, K = pair.S, pair.cone
    fp = s9f3a.FloatProfile(S)
    inv = [[float(c) for c in r] for r in S.inv]
    o = [float(c) for c in S.o]
    lo, hi = float(S.lo), float(S.hi)
    co = [float(c) for c in K.o]
    cx, cy, cn = ([float(c) for c in v] for v in (K.x, K.y, K.n))
    b, k, h = float(K.b), float(K.k), float(K.h)
    box = ([float(x) for x in pair.lo], [float(x) for x in pair.hi])
    vol = 1.0
    for i in range(3):
        vol *= box[1][i]-box[0][i]
    rng = random.Random(zlib.crc32(name.encode()))
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box[0][i]+(box[1][i]-box[0][i])*rng.random() for i in range(3)]
        d = [X[i]-o[i] for i in range(3)]
        u, v, w = (sum(inv[k_][i]*d[i] for i in range(3)) for k_ in range(3))
        inP = lo < w < hi and fp.inside(u, v)
        e = [X[i]-co[i] for i in range(3)]
        cu, cv, cw = (sum(a*c for a, c in zip(e, r)) for r in (cx, cy, cn))
        inK = 0 < cw < h and cu*cu+cv*cv < (b+k*cw)**2
        obj, tool = (inP, inK) if pair.obj_prism else (inK, inP)
        for op, inside in (('fuse', obj or tool), ('cut', obj and not tool), ('common', obj and tool)):
            if inside:
                a = acc[op]
                a[0] += 1
                for i in range(3):
                    a[1][i] += X[i]
                    a[2][i] += X[i]*X[i]
    out = {}
    for op, (kk, s1, s2) in acc.items():
        p = kk/n
        v = vol*p
        se = vol*math.sqrt(max(p*(1-p), 1e-300)/n)
        if kk > 1:
            mean = [x/kk for x in s1]
            sd = [math.sqrt(max(s2[i]/kk-mean[i]**2, 0.0)/kk) for i in range(3)]
        else:
            mean, sd = None, None
        out[op] = (v, se, mean, sd)
    return out


# ------------------------------------------------------------------ margins

INF = mp.mpf('inf')
sine = s9f3a.sine
wall_normal = s9f3a.wall_normal
prism_distance = s9f3a.prism_distance


def slab_depth(K, X):
    """How far a point lies inside the cone's end planes' slab (signed)."""
    w = K.local(X)[2]
    return min(w, K.hm-w)


def cone_distance(K, X):
    """The distance from a point to the cone's boundary (binary64-grade)."""
    u, v, w = K.local(X)
    rad = mp.sqrt(u*u+v*v)
    # The lateral surface: in the meridian plane, the segment from (b, 0) to
    # (t, h).
    p0, p1 = (K.bm, mp.mpf(0)), (K.tm, K.hm)
    dx, dy = p1[0]-p0[0], p1[1]-p0[1]
    s = ((rad-p0[0])*dx+(w-p0[1])*dy)/(dx*dx+dy*dy)
    s = min(max(s, 0), 1)
    side = mp.sqrt((rad-p0[0]-s*dx)**2+(w-p0[1]-s*dy)**2)
    best = side
    for e, r in ((mp.mpf(0), K.bm), (K.hm, K.tm)):
        if r > 0:
            dw = abs(w-e)
            best = min(best, dw if rad <= r else mp.sqrt(dw*dw+(rad-r)**2))
    return best


def margins(pair):
    S, K = pair.S, pair.cone
    out = {k: INF for k in ('axis', 'apex', 'turn', 'loop', 'knot', 'touch', 'cross', 'rim', 'trace', 'pierce',
                            'vertex')}
    lo, hi = M(S.lo), M(S.hi)
    A = M(pair.A)
    out['axis'] = abs(A)/M(dot(S.n, S.n))
    m = Mv(S.m)
    # The apex against every spline wall's surface (in the profile).
    au, av, _ = S.local(Mv(K.apex))
    for el in S.profile.elements:
        if el.kind == 'spline':
            out['apex'] = min(out['apex'], el.distance(au, av))

    def grad(X):
        return K.gradient(X)
    for el in S.elements:
        P, dP = S.wall(el)
        Aq, B, C, rows, q = ref.cone_quadratic(S, el, K)
        D = ref.psub(ref.pmul(B, B), ref.pc(C, Aq))
        Bm, Dm = ref.pm(B), ref.pm(D)
        if el.spline and Aq != 0:
            for t in ref.exact_roots(D, slack=mp.mpf(0)):
                w = -ref.peval(Bm, t)/A
                u, v = el.point(t)
                X = S.world(u, v, w)
                sd = slab_depth(K, X)
                outside = max(0, lo-w, w-hi, -sd)
                if outside == 0:
                    out['loop'] = min(out['loop'], w-lo, hi-w, sd)
                else:
                    out['turn'] = min(out['turn'], outside)
                out['knot'] = min(out['knot'], min(mp.sqrt((u-kn[0])**2+(v-kn[1])**2) for kn in S.profile.vm))
        if Aq != 0:
            for t in ref.exact_roots(ref.pder(D), slack=mp.mpf(0)):
                out['touch'] = min(out['touch'], 2*mp.sqrt(abs(ref.peval(Dm, t)))/abs(A))
        # The meeting at the caps.
        for wq in (S.lo, S.hi):
            for t in ref.exact_roots(ref.padd(ref.padd([Aq*wq*wq], ref.pc(B, 2*wq)), C), slack=mp.mpf(0)):
                u, v = el.point(t)
                X = S.world(u, v, M(wq))
                if slab_depth(K, X) < 0:
                    continue
                Nw, Nc = wall_normal(S, el, t), grad(X)
                out['cross'] = min(out['cross'], sine(Nw, Nc), 1-sine(cross(Nw, Nc), m))
        # The rims against the wall, and the end planes' traces.
        for e, r in K.ends:
            if r == 0:
                continue
            if q[2] != 0:
                wh = ref.pc(ref.psub([e], rows[2]), 1/q[2])
                whm = ref.pm(wh)
                out['trace'] = min(out['trace'], min(abs(ref.peval(whm, t)-kk) for t in (mp.mpf(0), mp.mpf(1))
                                                     for kk in (lo, hi)))
            else:
                for t in ref.exact_roots(ref.psub(rows[2], [e]), slack=mp.mpf(0)):
                    du, dv = el.deriv(t)
                    T = tuple(du*S.xm[i]+dv*S.ym[i] for i in range(3))
                    out['trace'] = min(out['trace'], sine(T, cross(K.nm, S.nm)))
            for t, w in ref.rim_wall_points(S, el, K, e, r, rows, q):
                if not lo <= w <= hi:
                    continue
                u, v = el.point(t)
                X = S.world(u, v, w)
                c = add(K.om, scale(K.nm, M(e)))
                Nw = wall_normal(S, el, t)
                Tr = cross(K.nm, sub(X, c))
                out['rim'] = min(out['rim'], abs(dot(Tr, Nw))/mp.sqrt(dot(Tr, Tr)*dot(Nw, Nw)))
                Nc = grad(X)
                out['cross'] = min(out['cross'], sine(Nw, Nc), 1-sine(cross(Nw, Nc), K.nm))
    # The rims against the caps; end planes parallel to the caps.
    for e, r in K.ends:
        if r == 0:
            continue
        c = add(K.om, scale(K.nm, M(e)))
        for wq in (S.lo, S.hi):
            p0 = Mv(ref.add(S.o, ref.scale(S.n, wq)))
            for X in ref.circle_plane(c, K.nm, M(r), m, p0):
                u, v, _ = S.local(X)
                if S.profile.inside(u, v):
                    Tr = cross(K.nm, sub(X, c))
                    out['rim'] = min(out['rim'], abs(dot(Tr, m))/mp.sqrt(dot(Tr, Tr)*dot(m, m)))
            if cross(K.n, S.m) == (0, 0, 0):
                out['trace'] = min(out['trace'], abs(dot(K.nm, sub(p0, c))))
    # Vertical edges against the cone.
    for vtx in S.profile.vm:
        P0 = S.world(vtx[0], vtx[1], mp.mpf(0))
        for X in ref.line_points(P0, S.nm, K):
            w = S.local(X)[2]
            if lo <= w <= hi and slab_depth(K, X) >= 0:
                Nc = grad(X)
                out['pierce'] = min(out['pierce'], abs(dot(S.nm, Nc))/mp.sqrt(dot(S.nm, S.nm)*dot(Nc, Nc)))
    # Vertices against the other's faces.
    for vtx in S.profile.vm:
        for w in (lo, hi):
            out['vertex'] = min(out['vertex'], cone_distance(K, S.world(vtx[0], vtx[1], w)))
    corners = []
    if K.real_apex:
        corners.append(Mv(K.apex))
    s = SEAMS[0] if not pair.obj_prism else SEAMS[1]
    cs = (M((1-s*s)/(1+s*s)), M(2*s/(1+s*s)))
    for e, r in K.ends:
        if r > 0:
            corners.append(add(K.om, add(scale(K.xm, cs[0]*M(r)), add(scale(K.ym, cs[1]*M(r)), scale(K.nm, M(e))))))
    for X in corners:
        out['vertex'] = min(out['vertex'], prism_distance(S, X))
    return out


MARGIN = mp.mpf(10)**-3
DEGENERATE_MARGIN = {APEX: 'apex', RULING: 'axis', TOUCH: 'touch'}


def near_coincidences(pair):
    """S9f.3a's, a class within rounding of zero (below 1e-30 of the size,
    a disc wholly on one side) taken as zero."""
    size = pair.size
    out = []
    for tag, name, cls, _ in pair.faces():
        for c, v in cls.items():
            if mp.mpf(10)**-30*size**2 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{tag} {name} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volumes()[op][0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    # Events of distinct slices closer than that, but two of the prism's
    # own vertices or cap edges' extremes (a turned prism's corners on one
    # slice but for rounding, which the integrand does not feel).
    own = ('vertex', 'cap_extreme')
    for way in (pair.first, pair.second):
        ev = sorted(way.events(), key=lambda e: e[0])
        for (p, tp, _), (q, tq, _) in zip(ev, ev[1:]):
            if mp.mpf(10)**-30*size < q-p < mp.mpf(10)**-9*size and not (tp in own and tq in own):
                out.append(f'breakpoints {mp.nstr(p, 12)} ({tp}) and {mp.nstr(q, 12)} ({tq})')
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
    vk, mk, ak = pair.cone.closed()
    checks['closed_prism'] = max(abs(one['P'][0]-vs)/size**3,
                                 max(abs(x-y) for x, y in zip(one['P'][1], ms))/size**4)
    checks['closed_cone'] = max(abs(one['B'][0]-vk)/size**3,
                                max(abs(x-y) for x, y in zip(one['B'][1], mk))/size**4)
    checks['face_classes'] = max(abs(cls['in']+cls['out']-closed) for _, _, cls, closed in pair.faces())/size**2
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')-as_-ak)/size**2
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
             'breaks': (len(one['breaks']), len(two['breaks'])), 'A': float(pair.A)}
    return name, rows, checks, z, near, margins(pair), stats


run = s9f3a.run


def frame_name(frame):
    for k, v in ALL_FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


def generate(results, listed):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, spline_cone_boolean_reference.py: expect KIND {STEP}, reason TEXT for a '
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
        cones = [x for x in (c.obj, c.tool) if x.cone is not None]
        assert len(cones) == 1, f'{c.name}: one cone'
        prism_case = c.tool if c.obj.cone is not None else c.obj
        assert s9f3a.s9f1_has_spline(prism_case), f'{c.name}: a spline prism'
        frame_name(c.obj.frame), frame_name(c.tool.frame)
        assert frame_name(cones[0].frame) in CONE_FRAMES, f'{c.name}: a cone in an exact frame'


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
