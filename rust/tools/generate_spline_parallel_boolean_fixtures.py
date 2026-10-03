#!/usr/bin/env python3
"""Fixtures for S9f.2a of REVIEW_NOTES.md: Booleans of spline prisms against
prisms with arc, circle or spline walls whose axes are exactly parallel, in
frames that differ (turned about the axis by a quarter or half turn
exactly, by `R125`'s rounded rotation, `TILT2` against `TILT` about the
tilted axis, or equal axes whose origins' offset is not binary64 in the
object's frame: S9a.2 takes one exact frame).

`boolean-spline-parallel-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); in every case one profile holds
a spline and the other an arc, a circle or a spline, the two normals are
exactly parallel and the frames differ. The profiles are S9a.2's and S9f.1's
(`generate_spline_any_boolean_fixtures.profiles`: the bulge, dome, blob,
wave, capsule, lens and lens hole), a quartic `hump`, discs, stadiums and a
square with a round hole. Frames are S9c.1's
(`generate_curved_boolean_fixtures.FRAMES`) and `FLIP` (the `XY` normal,
`x` along `-x`: a half turn), stored bit for bit by the kernel's
`Frame3::new`; `boolean-spline-parallel-frames.tsv` records them.

`boolean-spline-parallel-expected.tsv` gives per case, from
`curved_boolean_reference.py` with its parallel walls:

* `expect KIND S9f.2a`: the declared outcome (`solid`, `empty`, or
  `degenerate`: a spline wall tangent to a cylinder or to another spline
  wall along a generatrix, `Degenerate` in the decisions), then for a
  degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* every operation two ways: the slicing's volume and first moments against
  the divergence theorem over the face sweeps' kept pieces; `fuse = A + B -
  common` and `cut = A - common`; `area(fuse) + area(common) = area(A) +
  area(B)`; every face's classes summing to its closed-form area; both
  sides' shared areas equal; a second slicing direction;
* S9a.2's `SplinePair` (atoms by Green's theorem over classified pieces,
  meetings by root finding) on every pair whose map between the frames is
  exact (`XY` against `TURN` or `FLIP`: the tool's profile turned exactly
  into the object's frame; equal axes: the offset taken exactly), within
  1e-24 in exact frames and 1e-15 in `TILT` (S9a's frame coordinates take
  the stored axes as orthonormal);
* solid counts: each case's declared count, the reference's by its
  slicing's union-find;
* margins that flag near coincidences: S9c.1's scan (a face's class or a
  result's volume positive but below 1e-9 of the size, slicing breakpoints
  closer than that), and geometric ones: every 2D crossing of a spline with
  the other's arc or spline at a sine of at least 1e-3 (`cross`), every
  local minimum of their distance that is no crossing at least 1e-3
  (`touch`), S9f.1's margins of each spline prism against the other's
  planes (`tangent`, `axis`, `crease`), line edges crossing spline walls
  (`edge`) and every vertex of either prism (spline knots included) at
  least 1e-3 from the other's faces (`vertex`). The declared degenerate
  pairs must fail their margin (below 1e-12) and are exempt from the rest.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import os
from pathlib import Path
import struct

import mpmath as mp

from identity_reference import Boundary, Case, Spline, encode_boolean_case
from curve_surface_reference import stored_axes
import curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES as CURVED_FRAMES, disc, square
from generate_boolean_fixtures import path, spline
import generate_spline_any_boolean_fixtures as s9f1

ROOT = Path(__file__).resolve().parents[1]
BOOLEAN_OPERATION = 93
STEP = 'S9f.2a'

FRAMES = dict(CURVED_FRAMES, FLIP=(0.0, 0.0, 1.0, -1.0, 0.0, 0.0))

TOUCH_CYLINDER = 'a spline wall tangent to a cylinder along a generatrix'
TOUCH_SPLINE = 'two spline walls tangent along a generatrix'


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def profiles():
    """S9f.1's profiles, a quartic hump (a graph over [0, 8], poles (8, 0),
    (6, 6), (4, 1), (2, 5), (0, 0), closed by its base), discs, a stadium
    along `v` and a square with a round hole."""
    p = s9f1.profiles()
    p['hump'] = [path([(0.0, 0.0), (8.0, 0.0)], [
        None, spline(4, [(8.0, 0.0), (6.0, 6.0), (4.0, 1.0), (2.0, 5.0), (0.0, 0.0)])])]
    return p


def stadium_v(half, r):
    """A stadium of straight length `2 half` along `v` and radius `r` about
    the origin."""
    return Boundary(points=[(r, -half), (r, half), (-r, half), (-r, -half)],
                    segments=[None, (0.0, half, r, True), None, (0.0, -half, r, True)])


def ring(x0, y0, x1, y1, cx, cy, r):
    """A square with a round hole."""
    return [square(x0, y0, x1, y1), Boundary(circle=(cx, cy, r))]


prism = s9f1.prism


class Boolean(s9f1.Boolean):
    pass


def group(name, obj, tool, ops, **kw):
    """Cases of one pair: `ops` maps each operation to its solid count, or to
    `('empty',)`, or to `('degenerate', reason)`."""
    out = []
    for op, want in ops.items():
        if isinstance(want, int):
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'solid', want, **kw))
        elif want[0] == 'empty':
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'empty', 0, **kw))
        else:
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'degenerate', None, want[1], **kw))
    return out


ONE3 = {'fuse': 1, 'cut': 1, 'common': 1}


def cases():
    """Every class of S9f.2a's decisions: spline walls against a disc's
    cylinder (turned by `R125`, a quarter and a half turn), a stadium's arc
    and line walls about the tilted axis, a square's lines and its round
    hole, a capsule (its spline and arc) against a disc (arc against arc
    too); spline walls against spline walls (quarter and half turns, `R125`,
    `TILT2` against `TILT`, cubics against cubics, a quartic against a
    cubic, two capsules); rounding offsets in `TILT` (a spline against a
    disc and against a spline); and the declared degenerate classes: a
    dome touching a disc at its apex, and touching a half-turned dome
    there."""
    p = profiles()
    out = []
    # A spline against arcs and circles.
    out += group('bulge_disc_r125', (p['bulge'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 1.5)], at('R125', (11, 2, -1)), 0.0, 7.0), ONE3)
    out += group('dome_stadium_tilt2', (p['dome'], at('TILT', (0, 0, 0)), 0.0, 5.0),
                 ([stadium_v(1.2, 0.6)], at('TILT2', (2.0, 0.0, -1.5)), 1.0, 4.0), ONE3)
    out += group('blob_ring_turn', (p['blob'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 (ring(-3.0, -3.0, 3.0, 3.0, 0.0, 0.0, 1.75), at('TURN', (7.5, 3.25, 1.0)), 0.0, 3.0), ONE3)
    out += group('lens_disc_flip', (p['lens_hole'], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(0.0, 0.0, 1.25)], at('FLIP', (6.5, 6.0, -1.0)), 0.0, 6.0), ONE3)
    out += group('capsule_disc_turn', (p['capsule'], at('XY', (0, 0, 0)), 0.0, 2.0),
                 ([disc(0.0, 0.0, 0.75)], at('TURN', (-0.75, 0.5, 0.5)), 0.0, 1.0), ONE3)
    # A spline against a spline.
    out += group('bulge_dome_turn', (p['bulge'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 (p['dome'], at('TURN', (12.5, 1.0, -1.0)), 0.0, 7.0), ONE3)
    out += group('blob_lens_r125', (p['blob'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 (p['lens'], at('R125', (5.25, -2.5, 1.0)), 0.0, 3.0), ONE3)
    out += group('wave_wave_flip', (p['wave'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 (p['wave'], at('FLIP', (10.0, 12.5, 1.0)), 0.0, 3.0), {'fuse': 1, 'cut': 1, 'common': 2})
    out += group('dome_dome_tilt2', (p['dome'], at('TILT', (0, 0, 0)), 0.0, 5.0),
                 (p['dome'], at('TILT2', (3.0, 0.4, -0.8)), -1.0, 6.0), ONE3)
    out += group('hump_lens_turn', (p['hump'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 (p['lens'], at('TURN', (9.0, -2.5, -1.0)), 0.0, 5.0), ONE3)
    out += group('capsule_capsule_turn', (p['capsule'], at('XY', (0, 0, 0)), 0.0, 2.0),
                 (p['capsule'], at('TURN', (0.25, -0.5, 0.5)), 0.0, 1.0), ONE3)
    # Rounding offsets in one turned frame.
    out += group('bulge_disc_offset', (p['bulge'], at('TILT', (1, -2, 0.5)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 1.5)], at('TILT', (12.1, 0.3, 0.1)), -1.0, 3.0), ONE3)
    out += group('bulge_dome_offset', (p['bulge'], at('TILT', (1, -2, 0.5)), 0.0, 5.0),
                 (p['dome'], at('TILT', (10.6, 0.1, -0.7)), 0.0, 4.0), ONE3)
    # Declared degenerate: tangent walls.
    out += group('dome_disc_touch', (p['dome'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([disc(0.0, 0.0, 2.0)], at('TURN', (2.0, 4.0, -1.0)), 0.0, 7.0),
                 {op: ('degenerate', TOUCH_CYLINDER) for op in ('fuse', 'common')})
    out += group('dome_flip_touch', (p['dome'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 (p['dome'], at('FLIP', (4.0, 4.0, 1.0)), 0.0, 3.0),
                 {op: ('degenerate', TOUCH_SPLINE) for op in ('fuse', 'cut')})
    return out


# ------------------------------------------------------------------ S9a.2's SplinePair

def turn_boundaries(boundaries, L):
    """The boundaries under the exact linear map `L` (entries 0 and +-1, a
    rotation): points, arcs' and circles' centres, spline poles."""
    f = lambda p: (L[0][0]*p[0]+L[0][1]*p[1], L[1][0]*p[0]+L[1][1]*p[1])
    out = []
    for b in boundaries:
        if b.circle is not None:
            cx, cy, r = b.circle
            c = f((cx, cy))
            out.append(Boundary(circle=(c[0], c[1], r)))
            continue
        segs = None
        if b.segments is not None:
            segs = []
            for s in b.segments:
                if s is None:
                    segs.append(None)
                elif isinstance(s, Spline):
                    segs.append(Spline(s.degree, tuple(f(q) for q in s.poles), s.knots, s.mults))
                else:
                    cx, cy, r, ccw = s
                    c = f((cx, cy))
                    segs.append((c[0], c[1], r, ccw))
        out.append(Boundary(points=[f(q) for q in b.points], segments=segs))
    return out


def exact_map(obj, tool):
    """The tool's frame coordinates in the object's: the linear part `L`
    (entries 0 and +-1) and the offset `(a, b, c)`, exactly, or None when
    the map is no such turn (or the normals differ)."""
    A, B = ref.Prism(obj), ref.Prism(tool)
    if A.n != B.n:
        return None
    m = ref.parallel_map(A, B)
    L = ((m[0][1], m[0][2]), (m[1][1], m[1][2]))
    if any(v not in (0, 1, -1) for row in L for v in row):
        return None
    c = ref.apply(A.inv, ref.sub(B.o, A.o))[2]
    return L, (m[0][0], m[1][0], c)


def s9a2_result(obj, operation, tool):
    """S9a.2's reference with the tool's profile turned exactly into the
    object's frame and offset exactly: (solids, volume, area, centre)."""
    import boolean_reference as s9a
    L, (a, b, c) = exact_map(obj, tool)

    class TurnedPair(s9a.SplinePair):
        def __init__(self, obj, tool):
            self.axes = tuple(tuple(F(v) for v in w) for w in stored_axes(obj.frame))
            self.offset = (a, b, c)
            self.A = s9a.profile_elements(obj.boundaries, obj.tolerance)
            self.B = s9a.profile_elements(turn_boundaries(tool.boundaries, L), tool.tolerance, a, b)
            self.heights = {'A': tuple(sorted((F(obj.start), F(obj.end)))),
                            'B': tuple(sorted((F(tool.start)+c, F(tool.end)+c)))}
            coords = [abs(s9a.M(v)) for e in self.A+self.B for v in (*e.p, *e.q)]
            coords += [abs(s9a.M(e.c[i]))+s9a.M(e.r) for e in self.A+self.B if e.kind == 'A' for i in range(2)]
            coords += [abs(s9a.M(v)) for e in self.A+self.B if e.kind == 'S' for ctrl in e.ctrls for p in ctrl
                       for v in p]
            self.scale = max([s9a.M(1)]+coords)
            self.eps = s9a.M(10)**-25*self.scale
            self._meetings = {}
            self._bands()
            self.slice_atoms = self.atoms
            self.classes, self.atoms = self._green()

    solids, V, S, centre, _ = TurnedPair(obj, tool).result(operation)
    return solids, V, S, centre


# ------------------------------------------------------------------ margins

INF = mp.mpf('inf')


def curved_margins(pair):
    """The 2D crossings' sines (`cross`) and the near misses' distances
    (`touch`) of every spline with the other's arcs, circles and splines."""
    out = {'cross': INF, 'touch': INF}
    A, B = pair.A, pair.B
    for c in pair.slicing.cross2d:
        out['cross'] = min(out['cross'], c.sine)
    m = ref.parallel_map(A, B)
    inv = ref.invert_map(m)
    size = pair.size
    curved = lambda p: [el for el in p.profile.elements if el.kind != 'seg']
    reach = F(1, 100)
    for ea in curved(A):
        for eb in curved(B):
            if ea.kind == 'spline' and eb.kind == 'spline':
                P = ref.Curve2(ea.ctrl)
                Q = ref.Curve2(tuple(ref.map_point(m, c) for c in eb.ctrl))
                for _, _, dist in ref.span_near_misses(P, Q, reach):
                    out['touch'] = min(out['touch'], dist)
            elif ea.kind == 'spline':
                ctrl = tuple(ref.map_point(inv, c) for c in ea.ctrl)
                for _, dist in ref.span_circle_near(ctrl, eb, size):
                    out['touch'] = min(out['touch'], dist)
            elif eb.kind == 'spline':
                ctrl = tuple(ref.map_point(m, c) for c in eb.ctrl)
                for _, dist in ref.span_circle_near(ctrl, ea, size):
                    out['touch'] = min(out['touch'], dist)
    return out


def plane_margins(pair):
    """S9f.1's margins of each spline prism against the other's planes, line
    edges and vertices (`tangent`, `axis`, `crease`, `edge`, `vertex`)."""
    out = {k: INF for k in ('tangent', 'axis', 'crease', 'edge', 'vertex')}
    for S, Q in ((pair.A, pair.B), (pair.B, pair.A)):
        spans = [el for el in S.profile.elements if el.kind == 'spline']
        nS = S.n
        nlen = mp.sqrt(ref.M(ref.dot(nS, nS)))
        lo, hi = ref.M(S.lo), ref.M(S.hi)
        for f in Q.faces:
            if f.kind not in ('cap', 'wall'):
                continue
            N, P0 = f.normal, f.point
            a0, ax, ay, an = ref.dot(N, ref.sub(S.o, P0)), ref.dot(N, S.x), ref.dot(N, S.y), ref.dot(N, nS)
            Nlen = mp.sqrt(ref.M(ref.dot(N, N)))
            for el in spans:
                g = ref.padd(ref.padd(ref.pscale(el.X, ax), ref.pscale(el.Y, ay)), [a0])
                dg = ref.pder(g)
                if ref.pzero(dg):
                    continue
                crit = [t for t in ref.real_roots(dg) if -mp.mpf(10)**-30 <= t <= 1+mp.mpf(10)**-30]
                for tau in crit:
                    value = ref.peval([ref.M(c) for c in g], tau)
                    if an == 0:
                        out['tangent'] = min(out['tangent'], abs(value)/mp.sqrt(ref.M(ax*ax+ay*ay)))
                    else:
                        wstar = -value/ref.M(an)
                        out['crease'] = min(out['crease'], min(abs(wstar-lo), abs(wstar-hi))*nlen)
            if an != 0 and spans:
                out['axis'] = min(out['axis'], abs(ref.M(an))/(Nlen*nlen))
                if ax == 0 and ay == 0:
                    w0 = -a0/an
                    if w0 not in (S.lo, S.hi):
                        out['crease'] = min(out['crease'], min(abs(ref.M(w0)-lo), abs(ref.M(w0)-hi))*nlen)
        # Line cap edges of Q crossing spline walls (vertical edges run along
        # them).
        edges = []
        for el in Q.profile.elements:
            if el.kind != 'seg':
                continue
            for h in (Q.lo, Q.hi):
                edges.append((Q.world_exact(el.p[0], el.p[1], h), Q.world_exact(el.q[0], el.q[1], h)))
        for P, R in edges:
            lp, lr = ref.apply(S.inv, ref.sub(P, S.o)), ref.apply(S.inv, ref.sub(R, S.o))
            D = ref.Mv(ref.sub(R, P))
            if (lr[0]-lp[0], lr[1]-lp[1]) == (0, 0):
                continue
            for el in spans:
                for tau, t in el.point_events(lp[:2], (lr[0]-lp[0], lr[1]-lp[1])):
                    w = ref.M(lp[2])+t*ref.M(lr[2]-lp[2])
                    if not (0 <= tau <= 1 and 0 <= t <= 1 and lo <= w <= hi):
                        continue
                    dx, dy = el.tangent(tau)
                    T = tuple(dx*S.xm[i]+dy*S.ym[i] for i in range(3))
                    Nw = ref.cross(T, S.nm)
                    sine = abs(ref.dot(D, Nw))/mp.sqrt(ref.dot(D, D)*ref.dot(Nw, Nw))
                    out['edge'] = min(out['edge'], sine)
        # Vertices (a spline prism's knots included) against the other's
        # surface.
        for vtx in Q.profile.vertices:
            for h in (Q.lo, Q.hi):
                out['vertex'] = min(out['vertex'], s9f1.surface_distance(S, Q.world_exact(vtx[0], vtx[1], h)))
    return out


def margins(pair):
    return {**plane_margins(pair), **curved_margins(pair)}


MARGIN = mp.mpf(10)**-3
DEGENERATE_MARGIN = {TOUCH_CYLINDER: 'cross', TOUCH_SPLINE: 'cross'}


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, obj, tool, ops, opposite = job
    pair = ref.Pair(obj, tool)
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    A, B = pair.A, pair.B
    size = pair.size
    va, ma, aa = A.measures()
    vb, mb, ab = B.measures()
    vols = pair.volumes()
    checks = {}
    dev = mp.mpf(0)
    vf, mf = vols['fuse']
    vc, mc = vols['common']
    vt, mt = vols['cut']
    for x, y in ((vf, va+vb-vc), (vt, va-vc)):
        dev = max(dev, abs(x-y)/size**3)
    for i in range(3):
        dev = max(dev, abs(mf[i]-(ma[i]+mb[i]-mc[i]))/size**4, abs(mt[i]-(ma[i]-mc[i]))/size**4)
    checks['inclusion_exclusion'] = dev
    div = pair.divergence_volumes()
    dev = mp.mpf(0)
    for op in ref.OPS:
        dev = max(dev, abs(div[op][0]-vols[op][0])/size**3,
                  max(abs(x-y) for x, y in zip(div[op][1], vols[op][1]))/size**4)
    checks['divergence'] = dev
    fdev = mp.mpf(0)
    for tag, f, cls, _ in pair.face_areas():
        fdev = max(fdev, abs(sum(cls.values())-f.closed_area())/size**2)
    checks['face_classes'] = fdev
    total = lambda t, c: sum((cls[c] for tag, f, cls, _ in pair.face_areas() if tag == t), mp.mpf(0))
    checks['shared_both_sides'] = max(abs(total('A', 'same')-total('B', 'same')),
                                      abs(total('A', 'opp')-total('B', 'opp')))/size**2
    if not opposite:
        checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')-aa-ab)/size**2
    axes = [(F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))]
    used = min(axes, key=lambda e: abs(ref.dot(A.n, e)))
    other = [e for e in axes if e != used and ref.cross(A.n, e) != (0, 0, 0)][0]
    second = ref.Slicing(A, B, other).measure()
    dev = mp.mpf(0)
    for op in ref.OPS:
        dev = max(dev, abs(second[op][0]-vols[op][0])/size**3)
        dev = max(dev, max(abs(x-y) for x, y in zip(second[op][1], vols[op][1]))/size**4)
    checks['second_axis'] = dev
    near = s9f1.near_coincidences(pair)
    stats = {'volume_quadrature': pair.slicing.quad_error/size**4,
             'face_quadrature': max(sw.quad_error for _, _, _, sw in pair.face_areas())/size**2,
             'crossings': len(pair.slicing.cross2d)}
    return name, rows, {op: res[op][0] for op in ref.OPS}, res, checks, near, margins(pair), stats


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


def s9a2_job(c):
    n, V, S, C = s9a2_result(c.obj, c.operation, c.tool)
    return c.name, n, V, S, C


def reference_checks(results, workers):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, res, checks, near, margin, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    # S9a.2's SplinePair where the map between the frames is exact.
    by_pair = {r[0]: r for r in results}
    listed = [c for c in cases() if c.kind != 'degenerate' and exact_map(c.obj, c.tool) is not None]
    for name, n, V, S, C in run(listed, s9a2_job, workers):
        c = next(x for x in listed if x.name == name)
        got = by_pair[c.pair_name][3][c.operation]
        assert n == got[0], (c.name, n, got[0])
        if n == 0:
            continue
        size = max(abs(V)**(mp.mpf(1)/3), 1)
        key = 's9a2_exact_frames' if frame_name(c.obj.frame) == 'XY' else 's9a2_turned_frames'
        note(key, max(abs(got[1]-V)/abs(V), abs(got[2]-S)/abs(S),
                      max(abs(got[3][i]-C[i]) for i in range(3))/size))
    return worst, covered


def generate(results, listed, prefix):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, curved_boolean_reference.py with parallel walls: expect KIND {STEP}, reason '
           'TEXT for a degenerate case, then result N volume area cx cy cz or empty)']
    for case in listed:
        blocks.append(case.encode())
        _, rows, _, res, _, near, _, _ = by_pair[case.pair_name]
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
    return {f'{prefix}-cases.txt': '\n'.join(blocks)+'\n',
            f'{prefix}-expected.tsv': '\n'.join(out)+'\n',
            f'{prefix}-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(listed):
    pairs = {}
    for c in listed:
        pairs.setdefault(c.pair_name, [c.obj, c.tool, [], c.opposite])[2].append(c.operation)
    return [(name, obj, tool, ops, opposite) for name, (obj, tool, ops, opposite) in pairs.items()]


def has_spline(case):
    return any(b.segments is not None and any(isinstance(s, Spline) for s in b.segments) for b in case.boundaries)


def has_curve(case):
    return has_spline(case) or any(b.circle is not None or (b.segments and any(isinstance(s, tuple)
                                                                                for s in b.segments))
                                   for b in case.boundaries)


def validate(listed):
    """The fixture list's own rules (see the module's docstring)."""
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert has_spline(c.obj) or has_spline(c.tool), f'{c.name}: a spline'
        S, Q = (c.obj, c.tool) if has_spline(c.obj) else (c.tool, c.obj)
        assert has_curve(Q), f'{c.name}: the other profile holds an arc, a circle or a spline'
        A, B = ref.Prism(c.obj), ref.Prism(c.tool)
        assert ref.cross(A.n, B.n) == (0, 0, 0), f'{c.name}: exactly parallel axes'
        if stored_axes(c.obj.frame)[1:] == stored_axes(c.tool.frame)[1:]:
            a, b, cc = s9f1.exact_offset(c.obj, c.tool)
            assert any(F(float(v)) != v for v in (a, b, cc)), f'{c.name}: an S9a pair (a binary64 offset)'
        frame_name(c.obj.frame), frame_name(c.tool.frame)


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
    for name, rows, _, _, checks, near, margin, stats in results:
        first = declared[name]
        kinds = {c.kind for c in listed if c.pair_name == name}
        if args.only:
            print(name, {op: r[0] for op, r in rows.items()})
            print('  checks', {k: mp.nstr(v, 3) for k, v in checks.items()})
            print('  margins', {k: mp.nstr(v, 3) for k, v in margin.items()}, 'stats', stats)
            if near:
                print('  near', near)
        if 'degenerate' in kinds:
            key = DEGENERATE_MARGIN[first.reason]
            if not margin[key] < mp.mpf(10)**-12:
                failed.append(f'{name}: declared degenerate, {key} margin {margin[key]}')
            degenerate_margin[key] = max(degenerate_margin.get(key, mp.mpf(0)), margin[key])
            continue
        if near and not first.near:
            failed.append(f'{name}: near coincidences {near}')
        for key, value in margin.items():
            if value < MARGIN:
                failed.append(f'{name}: {key} margin {mp.nstr(value, 3)}')
            worst_margin[key] = min(worst_margin.get(key, INF), value)
    if failed:
        raise SystemExit('\n'.join(failed))
    if args.only:
        return
    worst, covered = reference_checks(results, args.workers)
    limits = {'inclusion_exclusion': 1e-30, 'divergence': 1e-30, 'face_classes': 1e-30,
              'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'second_axis': 1e-30,
              'volume_quadrature_estimate': 1e-30, 'face_quadrature_estimate': 1e-30,
              's9a2_exact_frames': 1e-24, 's9a2_turned_frames': 1e-15}
    for key, value in worst.items():
        assert value <= limits[key], (key, mp.nstr(value, 3))
    files = generate(results, cases(), 'boolean-spline-parallel')
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
    print('smallest margins (non-degenerate pairs):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst_margin.items())))
    print('declared degenerate margins (largest):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(degenerate_margin.items())))
    print('crossings:', sum(r[7]['crossings'] for r in results))


if __name__ == '__main__':
    main()
