#!/usr/bin/env python3
"""Fixtures for S9f.2b of REVIEW_NOTES.md: Booleans of spline prisms against
prisms with arc or circle walls whose axes cross (are not parallel): a
spline wall meets a cylinder in a curve that is, along the wall's ruling at
the spline's parameter `tau`, a root of `A w^2 + 2 B(tau) w + C(tau)` (degrees
0, `p` and `2 p`), a graph over `tau` between its turning points (the roots
of `B^2 - A C`).

`boolean-spline-crossing-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); in every case one profile holds
a spline and the other an arc or a circle and no spline, and the two
normals cross. The profiles are S9a.2's and S9f.1's
(`generate_spline_any_boolean_fixtures.profiles`: the bulge, dome, blob,
wave, capsule, lens and `knot`, R4's C1 knot of multiplicity two), discs,
a stadium and a square with a round hole. Frames are S9c.1's
(`generate_curved_boolean_fixtures.FRAMES`) and `STEEP` (the normal (0, 5,
12) / 13, leaning 22.6 degrees, `x` along the world's), stored bit for bit
by the kernel's `Frame3::new`; `boolean-spline-crossing-frames.tsv`
records them.

`boolean-spline-crossing-expected.tsv` gives per case, from
`curved_boolean_reference.py` with its crossing walls:

* `expect KIND STEP`: the declared outcome (`solid`, `empty`, or
  `degenerate`: the cylinder tangent to the wall, or a turning point at a
  knot, `Degenerate` in the decisions) and the sub-step that decides it:
  `S9f.2b.1` (every meeting a graph over the spline's parameter inside the
  faces) or `S9f.2b.2` (a meeting turning back inside the faces: loops,
  refused by S9f.2b.1's kernel); then for a degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* every operation two ways: the slicing's volume and first moments against
  the divergence theorem over the face sweeps' kept pieces; `fuse = A + B -
  common` and `cut = A - common`; `area(fuse) + area(common) = area(A) +
  area(B)`; every face's classes summing to its closed-form area; both
  sides' shared areas equal;
* perpendicular pairs in exact frames (the spline prism in `XY`, the
  cylinder on its `SIDE`): the common's volume and first moments as `int
  L(y) W(y) dy` (`L` the profile's chord along `x` at `y` within the
  cylinder's length, `W` the disc's height chord there within the spline
  prism's heights), a product of two chords independent of the slicing's
  polygons;
* solid counts: each case's declared count, the reference's by its
  slicing's union-find;
* margins that flag near coincidences: S9c.1's scan (a face's class or a
  result's volume positive but below 1e-9 of the size, slicing breakpoints
  closer than that), S9f.1's and S9f.2a's margins of each spline prism
  against the other's planes and edges (`tangent`, `axis`, `crease`, `edge`,
  `vertex`), and the crossing walls' own: every turning point of a meeting
  at least 1e-3 outside a face (`turn`; a declared loop's inside both faces
  by at least that, `loop`), the meeting's points at a cap crossing it at a
  sine of at least 1e-3 and the walls there at a sine of at least 1e-3
  (`cross`), the discriminant's critical values at least 1e-3 from zero in
  the gap they measure (`touch`), every turning point at least 1e-3 from a
  knot or vertex of the profile (`knot`), every vertical edge of either
  prism crossing the other's curved wall at a sine of at least 1e-3
  (`pierce`). The declared degenerate pairs must fail their margin (below
  1e-12) and are exempt from the rest.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import os
from pathlib import Path
import struct

import mpmath as mp

from identity_reference import Boundary, Spline
from curve_surface_reference import stored_axes
import curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES as CURVED_FRAMES, disc, square, stadium
import generate_spline_any_boolean_fixtures as s9f1
import generate_spline_parallel_boolean_fixtures as s9f2a

ROOT = Path(__file__).resolve().parents[1]
STEP = 'S9f.2b.1'
LOOPS = 'S9f.2b.2'

FRAMES = dict(CURVED_FRAMES, STEEP=(0.0, 5.0, 12.0, 1.0, 0.0, 0.0))

TOUCH = 'a cylinder tangent to a spline wall'
KNOT_TURN = "a spline wall's meeting with a cylinder turning back at a knot"


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def profiles():
    return s9f1.profiles()


prism = s9f1.prism


class Boolean(s9f1.Boolean):
    def __init__(self, name, operation, obj, tool, kind='solid', solids=1, reason=None, step=STEP, **kw):
        super().__init__(name, operation, obj, tool, kind, solids, reason, **kw)
        self.step = step


def group(name, obj, tool, ops, step=STEP, **kw):
    """Cases of one pair: `ops` maps each operation to its solid count, or to
    `('empty',)`, or to `('degenerate', reason)`."""
    out = []
    for op, want in ops.items():
        if isinstance(want, int):
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'solid', want, step=step, **kw))
        elif want[0] == 'empty':
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'empty', 0, step=step, **kw))
        else:
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'degenerate', None, want[1], step=step, **kw))
    return out


ONE3 = {'fuse': 1, 'cut': 1, 'common': 1}


def cases():
    """Every class of S9f.2b's decisions: meetings that run cap to cap as
    graphs over the spline's parameter (tilted rods along the bulge's,
    capsule's and wave's walls, the wave's and capsule's crossing their
    knots, a perpendicular rod covering the dome, a steep cylinder holding
    most of the blob as the object, a rod ending inside the blob, a
    stadium's arc and edges against the bulge, a tilted ring's hole around
    the dome); loops (S9f.2b.2: a perpendicular rod through the bulge's
    wall, a tilted rod through the lens's); and the declared degenerate
    classes: a rod touching the dome's apex, and a rod whose meeting with
    `knot`'s wall turns back at its knot."""
    p = profiles()
    out = []
    out += group('bulge_tilt', (p['bulge'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([disc(0.0, 0.0, 0.4)], at('TILT', (11.1, 0.1, -2.5)), 0.0, 10.0), ONE3)
    out += group('dome_side', (p['dome'], at('XY', (0, 0, 0)), 0.0, 2.0),
                 ([disc(1.0, -1.0, 2.0)], at('SIDE', (-1.0, 0.0, 0.0)), 0.0, 6.0), ONE3)
    out += group('blob_tilt_end', (p['blob'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([disc(0.0, 0.0, 0.4)], at('TILT', (7.8, 1.6, -1.7)), 0.0, 4.5), ONE3)
    out += group('wave_lean', (p['wave'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([disc(0.0, 0.0, 0.75)], at('LEAN', (6.725, 7.0, 0.0)), 0.0, 6.0), ONE3)
    out += group('capsule_tilt', (p['capsule'], at('XY', (0, 0, 0)), 0.0, 2.0),
                 ([disc(0.0, 0.0, 0.95)], at('TILT', (-0.75, -2.9, -3.0)), 0.0, 10.0), ONE3)
    out += group('steep_blob', ([disc(0.0, 0.0, 3.25)], at('STEEP', (5.25, 3.375, -1.5)), 0.0, 6.0),
                 (p['blob'], at('XY', (0, 0, 0)), 0.0, 3.0), ONE3)
    out += group('bulge_stadium_tilt', (p['bulge'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([stadium(1.0, 0.4)], at('TILT', (9.85, 0.75, -1.0)), 0.0, 7.0), ONE3)
    out += group('ring_dome', ([square(-4.0, -4.0, 4.0, 4.0), disc(0.0, 0.0, 1.6)], at('TILT', (2.0, 1.7, 1.1)),
                               0.0, 1.0),
                 (p['dome'], at('XY', (0, 0, 0)), 0.0, 3.0), ONE3)
    # S9f.2b.2: loops.
    out += group('bulge_side_loop', (p['bulge'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([disc(3.0, 1.5, 1.0)], at('SIDE', (-1.0, 0.0, 0.0)), 0.0, 14.0), ONE3, step=LOOPS)
    out += group('lens_tilt_loop', (p['lens'], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(0.0, 0.0, 0.5)], at('TILT', (5.0, 1.0, -0.5)), 0.0, 6.0), ONE3, step=LOOPS)
    # Declared degenerate.
    out += group('dome_touch', (p['dome'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([disc(3.0, 1.5, 1.0)], at('SIDE', (-1.0, 0.0, 0.0)), 0.0, 6.0),
                 {op: ('degenerate', TOUCH) for op in ('fuse', 'common')})
    out += group('knot_turn', (p['knot'], at('XY', (0, 0, 0)), 0.0, 3.0),
                 ([disc(5.5, 2.0, 1.5)], at('SIDE', (-1.0, 0.0, 0.0)), 0.0, 12.0),
                 {op: ('degenerate', KNOT_TURN) for op in ('cut', 'common')})
    return out


# ------------------------------------------------------------------ perpendicular pairs

def perpendicular(obj, tool):
    """The spline prism in `XY` and a disc on its `SIDE` (in either order):
    (spline case, disc case), else None."""
    for S, C in ((obj, tool), (tool, obj)):
        if frame_name(S.frame) == 'XY' and frame_name(C.frame) == 'SIDE' and \
                len(C.boundaries) == 1 and C.boundaries[0].circle is not None and s9f2a.has_spline(S):
            return S, C
    return None


def product_common(S, C):
    """The common of a spline prism in `XY` (heights `[a, b]`) and a disc
    prism on its `SIDE` (its axis along `x`, the disc `(y - cy)^2 + (z -
    cz)^2 < r^2` over `x` in `[x0, x1]`): volume and first moments as `int
    L(y) W(y) dy` and its moments, `L` the profile's chord along `x` at `y`
    clipped to `[x0, x1]`, `W` the disc's chord in `z` at `y` clipped to
    `[a, b]`."""
    A = ref.Prism(S)
    o = [F(c) for c in C.frame[:3]]
    cy, cz, r = (F(c) for c in C.boundaries[0].circle)
    cy, cz = cy+o[1], cz+o[2]
    x0, x1 = o[0]+F(C.start), o[0]+F(C.end)
    a, b = A.lo, A.hi
    prof = A.profile
    ys = [v[1] for v in prof.vertices]+[cy-r, cy+r]
    for el in prof.elements:
        if el.kind == 'spline':
            ys += [el.point(t)[1] for t in el.tangent_params((F(1), F(0)))]
    for h in (a, b):
        q = r*r-(h-cz)**2
        if q > 0:
            ys += [M for M in (ref.M(cy)-mp.sqrt(ref.M(q)), ref.M(cy)+mp.sqrt(ref.M(q)))]
    # Where a chord's end crosses the cylinder's ends.
    for el in prof.elements:
        if el.kind == 'spline':
            for x in (x0, x1):
                ys += [y for y in (el.point(t)[1] for t in ref.real_roots(ref.psub(el.X, [x])))]
        elif el.kind == 'seg' and el.e[0] != 0:
            for x in (x0, x1):
                s = (x-el.p[0])/el.e[0]
                if 0 <= s <= 1:
                    ys.append(el.p[1]+s*el.e[1])
    lo, hi = ref.M(cy-r), ref.M(cy+r)
    breaks = ref.merge_breaks([ref.M(y) for y in ys], lo, hi, mp.mpf(10)**-30)

    def f(y):
        q = ref.M(r)**2-(y-ref.M(cy))**2
        if q <= 0:
            return [mp.mpf(0)]*4
        s = mp.sqrt(q)
        z0, z1 = max(ref.M(a), ref.M(cz)-s), min(ref.M(b), ref.M(cz)+s)
        if z1 <= z0:
            return [mp.mpf(0)]*4
        W, Wz = z1-z0, (z1*z1-z0*z0)/2
        L = Lx = mp.mpf(0)
        for (t0, _), (t1, _) in prof.chords((mp.mpf(0), y), (mp.mpf(1), mp.mpf(0))):
            u0, u1 = max(t0, ref.M(x0)), min(t1, ref.M(x1))
            if u1 > u0:
                L += u1-u0
                Lx += (u1*u1-u0*u0)/2
        return [L*W, Lx*W, y*L*W, L*Wz]

    total = [mp.mpf(0)]*4
    for p, q in zip(breaks, breaks[1:]):
        est, _ = ref.integrate(f, p, q, mp.mpf(10)**-36)
        if est is not None:
            total = [x+y for x, y in zip(total, est)]
    return total[0], tuple(total[1:])


# ------------------------------------------------------------------ margins

INF = mp.mpf('inf')


def spline_and_cylinders(pair):
    A, B = pair.A, pair.B
    S, C = (A, B) if A.profile.splines else (B, A)
    spans = [el for el in S.profile.elements if el.kind == 'spline']
    circs = [el for el in C.profile.elements if el.kind in ('arc', 'circle')]
    return S, C, spans, circs


def wall_normal(S, el, tau):
    dx, dy = el.tangent(tau)
    T = tuple(dx*S.xm[i]+dy*S.ym[i] for i in range(3))
    return ref.cross(T, S.nm)


def cylinder_gradient(C, circ, X):
    """The gradient of the cylinder's function `(u - cu)^2 + (v - cv)^2 -
    r^2` in the world (halved)."""
    u, v, _ = C.local(X)
    du, dv = u-circ.cm[0], v-circ.cm[1]
    return tuple(du*C.invm[0][i]+dv*C.invm[1][i] for i in range(3))


def sine(a, b):
    c = ref.cross(a, b)
    return mp.sqrt(ref.dot(c, c)/(ref.dot(a, a)*ref.dot(b, b)))


def in_faces(S, C, circ, X):
    """How far `X` lies outside the spline wall's heights and the cylinder
    face (its heights and its arc's angles), 0 inside both."""
    ws = S.local(X)[2]
    u, v, wc = C.local(X)
    out_s = max(0, ref.M(S.lo)-ws, ws-ref.M(S.hi))
    phi = mp.atan2(v-circ.cm[1], u-circ.cm[0])
    out_c = max(0, ref.M(C.lo)-wc, wc-ref.M(C.hi), circ.range_bad(phi))
    depth = min(ws-ref.M(S.lo), ref.M(S.hi)-ws, wc-ref.M(C.lo), ref.M(C.hi)-wc)
    return max(out_s, out_c), depth


def crossing_margins(pair):
    """The crossing walls' margins (see the module's docstring): {'turn',
    'loop', 'cross', 'touch', 'knot', 'pierce'}."""
    out = {k: INF for k in ('turn', 'loop', 'cross', 'touch', 'knot', 'pierce')}
    S, C, spans, circs = spline_and_cylinders(pair)
    knots = [ref.Mv(v) for v in S.profile.vertices]
    for el in spans:
        for circ in circs:
            turns, D = ref.wall_cylinder_turns(S, el, C, circ)
            for tau, X in turns:
                if not 0 <= tau <= 1:
                    continue
                outside, depth = in_faces(S, C, circ, X)
                if outside == 0:
                    out['loop'] = min(out['loop'], depth)
                else:
                    out['turn'] = min(out['turn'], outside)
                x, y = el.point(tau)
                out['knot'] = min(out['knot'], min(mp.sqrt((x-k[0])**2+(y-k[1])**2) for k in knots))
            # The discriminant's critical values, as the gap between the
            # ruling's two meetings they measure.
            (_, q0), (_, q1), _ = ref._wall_local(S, el, C)
            A = ref.M(q0*q0+q1*q1)
            Dm = [ref.M(c) for c in D]
            for tau in ref.real_roots(ref.pder(D)):
                if 0 <= tau <= 1:
                    out['touch'] = min(out['touch'], 2*mp.sqrt(abs(ref.peval(Dm, tau)))/A)
            # The meeting's points at the caps of either prism.
            for H in (S, C):
                m = ref.Mv(ref.cross(H.x, H.y))
                for h in (H.lo, H.hi):
                    for tau, X in ref.wall_cylinder_events(S, el, C, circ, H, h):
                        if not 0 <= tau <= 1 or in_faces(S, C, circ, X)[0] > 0:
                            continue
                        Nw, Nc = wall_normal(S, el, tau), cylinder_gradient(C, circ, X)
                        T = ref.cross(Nw, Nc)
                        out['cross'] = min(out['cross'], sine(Nw, Nc), 1-sine(T, m))
    # Vertical edges: the cylinder prism's against the spline walls, the
    # spline prism's against the cylinders.
    for vtx in C.profile.vertices:
        P = C.world_exact(vtx[0], vtx[1], C.lo)
        lp = ref.apply(S.inv, ref.sub(P, S.o))
        d = ref.apply(S.inv, C.n)
        for el in spans:
            for tau, t in el.point_events(lp[:2], d[:2]):
                w = ref.M(lp[2])+t*ref.M(d[2])
                wc = ref.M(C.lo)+t
                if 0 <= tau <= 1 and ref.M(S.lo) <= w <= ref.M(S.hi) and ref.M(C.lo) <= wc <= ref.M(C.hi):
                    out['pierce'] = min(out['pierce'], 1-sine(ref.Mv(C.n), wall_normal(S, el, tau)))
    for vtx in S.profile.vertices:
        P = S.world_exact(vtx[0], vtx[1], 0)
        lp = C.local(ref.Mv(P))
        d = ref.apply(C.inv, S.n)
        for circ in circs:
            a = ref.M(d[0])**2+ref.M(d[1])**2
            u0, v0 = lp[0]-circ.cm[0], lp[1]-circ.cm[1]
            b = u0*ref.M(d[0])+v0*ref.M(d[1])
            cc = u0*u0+v0*v0-circ.rm**2
            disc = b*b-a*cc
            if disc < 0:
                continue
            for sgn in (-1, 1):
                w = (-b+sgn*mp.sqrt(disc))/a
                X = tuple(ref.M(P[i])+w*S.nm[i] for i in range(3))
                if in_faces(S, C, circ, X)[0] == 0:
                    out['pierce'] = min(out['pierce'], 1-sine(S.nm, cylinder_gradient(C, circ, X)))
    return out


def margins(pair):
    return {**s9f2a.plane_margins(pair), **crossing_margins(pair)}


MARGIN = mp.mpf(10)**-3
DEGENERATE_MARGIN = {TOUCH: 'touch', KNOT_TURN: 'knot'}


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
    perp = perpendicular(obj, tool)
    if perp is not None:
        V, Mom = product_common(*perp)
        checks['perpendicular_product'] = max([abs(V-vc)/size**3]+[abs(x-y)/size**4 for x, y in zip(Mom, mc)])
    near = s9f1.near_coincidences(pair)
    stats = {'volume_quadrature': pair.slicing.quad_error/size**4,
             'face_quadrature': max(sw.quad_error for _, _, _, sw in pair.face_areas())/size**2,
             'breaks': len(pair.slicing.breaks)}
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


def reference_checks(results):
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
    return worst, covered


def generate(results, listed, prefix):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = ['# case\trow (S9f.2b, curved_boolean_reference.py with crossing walls: expect KIND STEP, reason '
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
        out.append(f'{case.name}\texpect {case.kind} {case.step}')
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


has_spline = s9f2a.has_spline


def has_round(case):
    return any(b.circle is not None or (b.segments and any(isinstance(s, tuple) for s in b.segments))
               for b in case.boundaries)


def validate(listed):
    """The fixture list's own rules (see the module's docstring)."""
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert has_spline(c.obj) != has_spline(c.tool), f'{c.name}: one spline prism'
        S, Q = (c.obj, c.tool) if has_spline(c.obj) else (c.tool, c.obj)
        assert has_round(Q), f'{c.name}: the other profile holds an arc or a circle'
        A, B = ref.Prism(c.obj), ref.Prism(c.tool)
        assert ref.cross(A.n, B.n) != (0, 0, 0), f'{c.name}: crossing axes'
        assert c.step in (STEP, LOOPS)
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
        loops = first.step == LOOPS
        if loops != (margin['loop'] < INF):
            failed.append(f'{name}: declared {first.step}, turning points inside the faces: {margin["loop"]}')
        for key, value in margin.items():
            if value < MARGIN:
                failed.append(f'{name}: {key} margin {mp.nstr(value, 3)}')
            worst_margin[key] = min(worst_margin.get(key, INF), value)
    if failed:
        raise SystemExit('\n'.join(failed))
    if args.only:
        return
    worst, covered = reference_checks(results)
    limits = {'inclusion_exclusion': 1e-30, 'divergence': 1e-30, 'face_classes': 1e-30,
              'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'perpendicular_product': 1e-30,
              'volume_quadrature_estimate': 1e-30, 'face_quadrature_estimate': 1e-30}
    for key, value in worst.items():
        assert value <= limits[key], (key, mp.nstr(value, 3))
    files = generate(results, cases(), 'boolean-spline-crossing')
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
    steps = {}
    for c in listed:
        steps[c.step] = steps.get(c.step, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(steps.items())))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('smallest margins (non-degenerate pairs):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst_margin.items())))
    print('declared degenerate margins (largest):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(degenerate_margin.items())))


if __name__ == '__main__':
    main()
