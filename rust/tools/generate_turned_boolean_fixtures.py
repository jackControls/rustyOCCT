#!/usr/bin/env python3
"""Fixtures for S9c.2b.1 of REVIEW_NOTES.md: Booleans of two cylinders in
turned frames (affine models on the stored axes, elliptic in the world)
whose axes cross, every vertex still a quadratic surd (a generatrix of one
input meeting the other cylinder) and the quartic section crossing no cap's
circle, before any of its kernel code.

The frames, the reference (`curved_boolean_reference.py`), the per-pair
checks (inclusion and exclusion with each operation sliced apart, every
face's classes against its closed-form area, the area identity, Monte
Carlo), the near-coincidence scan and the closed-form machinery
(`Perpendicular`, `legendre`, `ops_from_common`) are
`generate_curved_boolean_fixtures.py`'s and
`generate_procedural_boolean_fixtures.py`'s, imported; only the curved
generator's frames are used (their stored axes the kernel's bit for bit:
`boolean-turned-frames.tsv` records them). Every pair holds two cylinders
with axes not parallel, at least one in a turned frame (`TILT`, `TILT2`,
`TILTX`, `LEAN`, `R125`), so the two are not circular in a common measure.

`boolean-turned-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); `boolean-turned-expected.tsv`
gives per case `expect KIND S9c.2b.1` (`solid` or `degenerate`, then
`reason TEXT`) and `result N volume area cx cy cz`, as the curved and
procedural fixtures.

Classes (the S9c.2 decisions: in the frame of the two axes, with `e` along
`a x b`, the intervals `[eA - rA, eA + rA]` and `[eB - rB, eB + rB]`): a
partial bite, one loop (`bite_tiltx`, `XY` against `TILTX`; `skew_lean`,
unequal skew cylinders in `LEAN` against `TILTX`); a thin pipe through a
thick cylinder, two rings, the pipe the object (`rings_pipe`, `TILT2`
against `XY`; its cut two stubs); a pipe ending inside the other, its cap's
disc inside and its circle clear of the other's wall, one ring
(`blind_tilt`, `TILT` against `SIDE`, perpendicular); a box with a round
hole crossed by a tilted pipe, two rings on the hole's wall (`hole_tiltx`,
its common two pieces); and equal radii with axes meeting (`node_lean`,
`LEAN` against `TILTX`; `node_r125`, `R125` against `SIDE`), S9c.1's two
ellipses in the ideal frames but in the stored frames' models two rings
whose intervals' ends lie 2.4e-17 and 8.3e-17 apart, so two loops within
the resolution of each other: declared `degenerate`. (Where one frame's
stored `x` or `y` lies along `e` exactly, as `XY` against `TILT` or `LEAN`
and `TILT` against `TILTX`, the models' intervals are equal exactly: the
cylinders touch at both ends, a double tangency, and meet in two conics.)

Before writing, besides the curved generator's per-pair checks (1e-30 and
5 standard errors):

* closed forms in the ideal orthonormal frames, within 1e-15 (the stored
  axes are not exactly orthonormal): the common of cylinders crossing each
  other whole at angle `phi` as the perpendicular one (one quadrature in
  `eta`) over `sin phi`, centred on the axes' common perpendicular, and so
  the three operations' volumes and centres; the perpendicular pair
  (`blind_tilt`) by the quadrature with its caps, areas included; a box's
  round hole as the box less its cylinder, the pipe between the box's two
  parallel walls by Cavalieri; equal radii by Legendre's form at `k = 1`,
  `16 r^3 / (3 sin phi)`, and its area `16 r^2 / sin phi`;
* the curve's class from the ideal intervals (as declared, their ends at
  least 0.1 apart unless declared a node) and, for a node, the stored
  models' intervals (their ends apart, but within 1e-15 of each other: the
  loops' separation about the square root of that times the radius);
* every cap's circle against the other input's cylinders: its crossings
  with the cylinder's model (sign changes of the model's function along the
  circle, refined by bisection) must lie outside the cylinder's face by at
  least 0.1 along its axis, and a circle without crossings must keep off
  the cylinder by 1e-3 of the radius squared (no tangency between samples),
  so the quartic meets no cap's circle.

The near-coincidence scan must find none outside the declared nodes. Pairs
run in worker processes (`--workers`).
"""
import argparse
import math
import os
import struct

import mpmath as mp

from curve_surface_reference import stored_axes
import curved_boolean_reference as ref
import generate_curved_boolean_fixtures as curved
import generate_procedural_boolean_fixtures as procedural
from generate_curved_boolean_fixtures import FRAMES, at, disc, group, run, square
from generate_procedural_boolean_fixtures import cylinder, box_with_hole, legendre, ops_from_common

ROOT = curved.ROOT
PREFIX = 'boolean-turned'
STEP = 'S9c.2b.1'

SOLID3 = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
NODE = "near node: the stored frames' models meet in two loops within the resolution of each other"
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
TURNED = {'TILT', 'TILT2', 'TILTX', 'LEAN', 'R125'}
# Each pair's curve, from the decisions' intervals.
CURVES = {'bite_tiltx': 'loop', 'rings_pipe': 'rings', 'blind_tilt': 'rings', 'hole_tiltx': 'rings',
          'skew_lean': 'loop', 'node_lean': 'node', 'node_r125': 'node'}
# The pipe of `blind_tilt` (radius 1 along x) about the world point (0, -0.14, 0.48):
# eta = -0.4, zeta = 0.3 in TILT's (y, n).
BLIND = (-0.14, 0.48)


def cases():
    """S9c.2b.1's classes in turned frames (see the module's docstring)."""
    out = []
    out += group('bite_tiltx', ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -6.0, 6.0),
                 ([disc(1.2, 0.0, 1.5)], at('TILTX', (0, 0, 0)), -8.0, 8.0), SOLID3)
    out += group('rings_pipe', ([disc(0.0, 0.6, 1.0)], at('TILT2', (0, 0, 0)), -8.0, 8.0),
                 ([disc(0.0, 0.0, 2.5)], at('XY', (0, 0, 0)), -7.0, 7.0), SOLID3)
    out += group('blind_tilt', ([disc(0.0, 0.0, 2.0)], at('TILT', (0, 0, 0)), -3.0, 3.0),
                 ([disc(BLIND[0], BLIND[1], 1.0)], at('SIDE', (-6, 0, 0)), 0.0, 6.5),
                 {'fuse': 'solid', 'cut': 'solid'})
    out += group('hole_tiltx', ([square(0.0, 0.0, 10.0, 10.0), disc(5.0, 5.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 10.0),
                 ([disc(0.8, 0.0, 0.6)], at('TILTX', (5, 5, 5)), -9.0, 9.0), SOLID3)
    out += group('skew_lean', ([disc(0.0, 0.0, 2.0)], at('LEAN', (0, 0, 0)), -6.0, 6.0),
                 ([disc(1.2, 0.0, 1.3)], at('TILTX', (0, 0, 0)), -8.0, 8.0),
                 {'cut': 'solid', 'common': 'solid'})
    out += group('node_lean', ([disc(0.0, 0.0, 1.5)], at('LEAN', (0, 0, 0)), -6.0, 6.0),
                 ([disc(0.0, 0.0, 1.5)], at('TILTX', (0, 0, 0)), -8.0, 8.0),
                 {'common': ('degenerate', NODE)})
    out += group('node_r125', ([disc(0.0, 0.0, 1.5)], at('R125', (0, 0, 0)), -5.0, 5.0),
                 ([disc(0.0, 0.0, 1.5)], at('SIDE', (-6, 0, 0)), 0.0, 12.0),
                 {'fuse': ('degenerate', NODE)})
    return out


# ------------------------------------------------------------------ ideal geometry

def unit(v):
    n = mp.sqrt(sum(c*c for c in v))
    return tuple(c/n for c in v)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def ideal_axes(frame):
    """The ideal orthonormal (x, y, n) of a frame (`stored_axes`'s steps in
    mpmath from the frame's integer directions)."""
    name = curved.frame_name(frame)
    v = FRAMES[name]
    n = unit(tuple(mp.mpf(c) for c in v[:3]))
    y = unit(cross(n, tuple(mp.mpf(c) for c in v[3:])))
    return cross(y, n), y, n


def axis_of(prism, circle):
    """A cylinder's ideal axis: a point and a unit direction."""
    x, y, n = ideal_axes(prism.frame)
    o = tuple(mp.mpf(c) for c in prism.frame[:3])
    cx, cy = mp.mpf(circle[0]), mp.mpf(circle[1])
    return tuple(o[i]+cx*x[i]+cy*y[i] for i in range(3)), n


def along(prism, circle, h):
    p, n = axis_of(prism, circle)
    return tuple(p[i]+mp.mpf(h)*n[i] for i in range(3))


def crossing_common(pa, a, rA, pb, b, rB):
    """The common of two infinite cylinders with crossing axes (points `pa`,
    `pb`, unit directions `a`, `b`): the perpendicular one (one quadrature
    in `eta`, `procedural.Perpendicular`) over `sin phi`, its centre on the
    axes' common perpendicular at the mean `eta`. Returns (volume, world
    first moments)."""
    w = cross(a, b)
    s = mp.sqrt(dot(w, w))
    e = tuple(c/s for c in w)
    eA, eB = dot(pa, e), dot(pb, e)
    P = procedural.Perpendicular(rA, 0, 0, -50, 50, rB, eB-eA, 0, -50, 50, procedural.WORLD)
    vp, m, _ = P.volume()
    eta = m[1]/vp                                  # WORLD: eta along y, from A's axis
    # The foot of the common perpendicular on A's axis.
    d = tuple(pb[i]-pa[i] for i in range(3))
    c = dot(a, b)
    t = (dot(d, a)-c*dot(d, b))/(1-c*c)
    foot = tuple(pa[i]+t*a[i] for i in range(3))
    centre = tuple(foot[i]+eta*e[i] for i in range(3))
    vc = vp/s
    return vc, tuple(vc*x for x in centre)


def whole_cylinder(case_prism):
    """Volume, world first moments and area of a disc prism (ideal)."""
    circle = case_prism.boundaries[0].circle
    lo, hi = case_prism.start, case_prism.end
    return cylinder(circle[2], lo, hi, along(case_prism, circle, (mp.mpf(lo)+hi)/2))


def closed_forms():
    """{pair: {op: (volume, moments, area or None)}} in the ideal frames."""
    by_pair = {}
    for c in cases():
        by_pair.setdefault(c.pair_name, c)
    out = {}

    def crossing(name, cut_object='A'):
        c = by_pair[name]
        ca, cb = c.obj.boundaries[-1].circle, c.tool.boundaries[-1].circle
        pa, a = axis_of(c.obj, ca)
        pb, b = axis_of(c.tool, cb)
        vc, mc = crossing_common(pa, a, ca[2], pb, b, cb[2])
        return c, vc, mc

    for name in ('bite_tiltx', 'rings_pipe', 'skew_lean'):
        c, vc, mc = crossing(name)
        A, B = whole_cylinder(c.obj), whole_cylinder(c.tool)
        out[name] = ops_from_common(A, B, (vc, mc, None, None))
    # A box with a round hole (radius 3 about (5, 5), z in [0, 10]) and a
    # pipe (radius 0.6 about (5.8, 5, 5) along (0, -4, 3) / 5): the pipe
    # between the box's walls y = 0 and y = 10 (Cavalieri: its axis 12.5
    # long between them, centred at (5.8, 5, 5)) less its common with the
    # hole's cylinder.
    c, vh, mh = crossing('hole_tiltx')
    r = mp.mpf(0.6)
    vbox = mp.pi*r*r*10/(mp.mpf(4)/5)
    centre = (mp.mpf(5.8), mp.mpf(5), mp.mpf(5))
    vc = vbox-vh
    mc = tuple(vbox*centre[i]-mh[i] for i in range(3))
    A = box_with_hole(10, 10, 5, 5, 3)
    B = whole_cylinder(c.tool)
    out['hole_tiltx'] = ops_from_common(A, B, (vc, mc, None, None))
    # Perpendicular (TILT against SIDE): the quadrature with the pipe's cap,
    # in the frame X = x, E = TILT's y, Z = TILT's n.
    c = by_pair['blind_tilt']
    x, y, n = ideal_axes(c.obj.frame)
    world = ((0, 0, 0), (1, 0, 0), y, n)
    py, pz = mp.mpf(BLIND[0]), mp.mpf(BLIND[1])
    eB, zB = py*y[1]+pz*y[2], py*n[1]+pz*n[2]
    P = procedural.Perpendicular(2, 0, 0, -3, 3, 1, eB, zB, -6, 0.5, world)
    vc, mc, _ = P.volume()
    a_in, b_in = P.areas()
    A, B = whole_cylinder(c.obj), whole_cylinder(c.tool)
    out['blind_tilt'] = ops_from_common(A, B, (vc, mc, a_in, b_in))
    # Equal radii with meeting axes: Legendre at k = 1 over sin phi, centred
    # at the axes' crossing; each wall's part inside the other 8 r^2 / sin phi.
    for name in ('node_lean', 'node_r125'):
        c = by_pair[name]
        ca, cb = c.obj.boundaries[0].circle, c.tool.boundaries[0].circle
        pa, a = axis_of(c.obj, ca)
        pb, b = axis_of(c.tool, cb)
        s = mp.sqrt(dot(cross(a, b), cross(a, b)))
        r = mp.mpf(ca[2])
        vc = legendre(r, r)/s
        half = 8*r*r/s
        A, B = whole_cylinder(c.obj), whole_cylinder(c.tool)
        out[name] = ops_from_common(A, B, (vc, (0, 0, 0), half, half))
    return out


# ------------------------------------------------------------------ classes and caps

def cylinders(case_prism):
    return [b.circle for b in case_prism.boundaries if b.circle is not None]


def intervals(case):
    """The two cylinders' `eta` intervals (centre, half width), ideal and in
    the stored models (`e` from the stored normals)."""
    ca, cb = cylinders(case.obj)[-1], cylinders(case.tool)[-1]
    pa, a = axis_of(case.obj, ca)
    pb, b = axis_of(case.tool, cb)
    e = unit(cross(a, b))
    ideal = [(dot(pa, e), mp.mpf(ca[2])), (dot(pb, e), mp.mpf(cb[2]))]
    stored = []
    for c, circle in ((case.obj, ca), (case.tool, cb)):
        o, x, y, n = (tuple(mp.mpf(v) for v in w) for w in stored_axes(c.frame))
        stored.append((o, x, y, n, circle))
    es = unit(cross(stored[0][3], stored[1][3]))
    model = []
    for o, x, y, n, circle in stored:
        centre = tuple(o[i]+mp.mpf(circle[0])*x[i]+mp.mpf(circle[1])*y[i] for i in range(3))
        model.append((dot(centre, es), mp.mpf(circle[2])*mp.sqrt(dot(x, es)**2+dot(y, es)**2)))
    return ideal, model


def classify(ivs):
    """('rings' | 'loop' | 'node', the smallest distance between the
    intervals' ends)."""
    (c1, h1), (c2, h2) = ivs
    ends = [c1-h1, c1+h1, c2-h2, c2+h2]
    gap = min(abs(ends[0]-ends[2]), abs(ends[1]-ends[3]), abs(ends[0]-ends[3]), abs(ends[1]-ends[2]))
    lo1, hi1, lo2, hi2 = ends
    if (lo1 < lo2 and hi2 < hi1) or (lo2 < lo1 and hi1 < hi2):
        return 'rings', gap
    if hi1 <= lo2 or hi2 <= lo1:
        return 'apart', gap
    return 'loop', gap


def cap_circles(case):
    """Every cap's circle of one input against every cylinder of the other,
    in the stored models (floats): (the least distance along the cylinder's
    axis from its face of a crossing, the least |f| / r^2 of a circle
    without crossings, crossings counted)."""
    A, B = ref.Prism(case.obj), ref.Prism(case.tool)
    least_out, least_clear, count = math.inf, math.inf, 0
    N = 8192
    for P, Q in ((A, B), (B, A)):
        fq = curved.FloatPrism(Q)
        o = [float(c) for c in P.o]
        x, y, n = ([float(c) for c in v] for v in (P.x, P.y, P.n))
        for el in P.profile.elements:
            if el.kind != 'circle':
                continue
            cx, cy, r = float(el.c[0]), float(el.c[1]), float(el.r)
            for h in (float(P.lo), float(P.hi)):
                def point(t):
                    u, v = cx+r*math.cos(t), cy+r*math.sin(t)
                    return [o[i]+u*x[i]+v*y[i]+h*n[i] for i in range(3)]

                for el2 in Q.profile.elements:
                    if el2.kind != 'circle':
                        continue
                    qx, qy, qr = float(el2.c[0]), float(el2.c[1]), float(el2.r)

                    def local(X):
                        d = [X[i]-fq.o[i] for i in range(3)]
                        return [sum(fq.inv[k][i]*d[i] for i in range(3)) for k in range(3)]

                    def f(t):
                        u, v, _ = local(point(t))
                        return ((u-qx)**2+(v-qy)**2-qr*qr)/(qr*qr)

                    ts = [2*math.pi*k/N for k in range(N+1)]
                    fs = [f(t) for t in ts]
                    roots = []
                    for k in range(N):
                        if (fs[k] < 0) != (fs[k+1] < 0):
                            a, b, fa = ts[k], ts[k+1], fs[k]
                            for _ in range(80):
                                m = (a+b)/2
                                fm = f(m)
                                if (fm < 0) == (fa < 0):
                                    a, fa = m, fm
                                else:
                                    b = m
                            roots.append((a+b)/2)
                    if not roots:
                        least_clear = min(least_clear, min(abs(v) for v in fs))
                    for t in roots:
                        w = local(point(t))[2]
                        count += 1
                        least_out = min(least_out, max(fq.lo-w, w-fq.hi))
    return least_out, least_clear, count


def geometry_checks(case):
    """The pair's class and its caps' circles: (ideal class and gap, model
    class and gap, cap-circle margins)."""
    ideal, model = intervals(case)
    return classify(ideal), classify(model), cap_circles(case)


def geometry_job(case):
    return case.pair_name, geometry_checks(case)


# ------------------------------------------------------------------ checks

def reference_checks(results, geometry):
    """Largest deviations of each check (relative to the case's size) and
    the geometry's margins."""
    worst, covered = {}, {}

    def note(key, value):
        assert not mp.isnan(value), key
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    by_pair = {r[0]: r for r in results}
    ops_of = {}
    for c in cases():
        ops_of.setdefault(c.pair_name, []).append(c.operation)
    for name, forms in closed_forms().items():
        _, _, _, res, _, _, _ = by_pair[name]
        for op in ref.OPS:
            if res[op][0] == 0:
                continue
            note('closed_forms_turned_frames', procedural.deviation(res[op], forms[op]))
    for name, rows, _, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    margins = {'class_gap': mp.inf, 'node_model_gap': mp.mpf(0), 'cap_crossing_outside': math.inf,
               'cap_circle_clear': math.inf, 'cap_crossings': 0}
    for name, ((kind, gap), (mkind, mgap), (out, clear, count)) in geometry.items():
        want = CURVES[name]
        if want == 'node':
            assert gap < mp.mpf(10)**-30, (name, 'ideal intervals not equal', gap)
            assert 0 < mgap < mp.mpf(10)**-15, (name, 'model intervals not within the resolution', mgap)
            margins['node_model_gap'] = max(margins['node_model_gap'], mgap)
        else:
            assert kind == want and mkind == want, (name, kind, mkind, want)
            assert gap > mp.mpf(1)/10, (name, 'interval ends close', gap)
            margins['class_gap'] = min(margins['class_gap'], gap)
        if count:
            assert out > 0.1, (name, 'a cap circle meets the other cylinder on its face', out)
            margins['cap_crossing_outside'] = min(margins['cap_crossing_outside'], out)
        if clear < math.inf:
            assert clear > 1e-3, (name, 'a cap circle nearly tangent to the other cylinder', clear)
            margins['cap_circle_clear'] = min(margins['cap_circle_clear'], clear)
        margins['cap_crossings'] += count
    return worst, covered, margins


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, curved_boolean_reference.py: expect KIND STEP, reason TEXT for a degenerate '
           'case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, res, _, near, _ = by_pair[case.pair_name]
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
        pairs.setdefault(c.pair_name, [c.obj, c.tool, [], c.opposite])[2].append(c.operation)
    return [(name, obj, tool, ops, opposite, mc_n) for name, (obj, tool, ops, opposite) in pairs.items()]


LIMITS = dict(procedural.LIMITS)


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    assert set(CURVES) == {c.pair_name for c in listed}
    for c in listed:
        frames = [curved.frame_name(f) for f in c.frames]
        assert all(f in FRAMES for f in frames)
        assert any(f in TURNED for f in frames), f'{c.name}: no turned frame'
        _, _, _, na = stored_axes(c.obj.frame)
        _, _, _, nb = stored_axes(c.tool.frame)
        assert cross(na, nb) != (0.0, 0.0, 0.0), f'{c.name}: parallel axes'
        assert cylinders(c.obj) and cylinders(c.tool), f'{c.name}: not two cylinders'
        assert (c.kind == 'degenerate') == (CURVES[c.pair_name] == 'node'), c.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    firsts = {}
    for c in listed:
        firsts.setdefault(c.pair_name, c)
    geometry = dict(run(list(firsts.values()), geometry_job, args.workers))
    results = run(jobs(args.samples), curved.evaluate, args.workers)
    for name, _, _, _, _, near, _ in results:
        if near and CURVES[name] != 'node':
            raise SystemExit(f'{name}: near coincidences {near}')
    worst, covered, margins = reference_checks(results, geometry)
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('geometry:', ', '.join(f'{k} {v}' if isinstance(v, int) else f'{k} {float(v):.3g}'
                                 for k, v in margins.items()))
    for name, ((kind, gap), (mkind, mgap), (out, clear, count)) in geometry.items():
        print(f'  {name}: {kind} (ends {mp.nstr(gap, 3)} apart; model {mkind}, {mp.nstr(mgap, 3)}), '
              f'cap circles: {count} crossings, outside by {out:.3g}, clear by {clear:.3g}')
    for name, _, _, _, _, near, _ in results:
        if near:
            print('declared near coincidences (a node):', name, '; '.join(near[:4]))
    by_pair = {r[0]: r for r in results}
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
