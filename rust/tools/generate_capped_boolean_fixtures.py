#!/usr/bin/env python3
"""Fixtures for S9c.2b.2 of REVIEW_NOTES.md: Booleans of two cylinders
whose quartic section crosses a cap's circle within both faces (the vertex
there a nested surd in exact frames, a root of the circle's quartic in
turned ones), before any of its kernel code.

The frames, the reference (`curved_boolean_reference.py`), the per-pair
checks (inclusion and exclusion with each operation sliced apart, every
face's classes against its closed-form area, the area identity, Monte
Carlo), the near-coincidence scan and the closed-form machinery
(`Perpendicular`, `ops_from_common`, the ideal frames and the classes from
the intervals) are `generate_curved_boolean_fixtures.py`'s,
`generate_procedural_boolean_fixtures.py`'s and
`generate_turned_boolean_fixtures.py`'s, imported; only the curved
generator's frames are used (their stored axes the kernel's bit for bit:
`boolean-capped-frames.tsv` records them).

`boolean-capped-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); `boolean-capped-expected.tsv`
gives per case `expect KIND S9c.2b.2` (`solid`, `empty` or `degenerate`,
then `reason TEXT`) and `result N volume area cx cy cz` or `empty`, as the
curved, procedural and turned fixtures.

Pairs (the classes the S9c.2 decisions' intervals give). Exact frames (`XY`
against `SIDE`, perpendicular): a pipe across the rim of a thicker
cylinder's top cap, its axis just below the cap's plane, two rings each
crossing the cap's circle twice (`rim_pipe`); a partial bite at the rim,
one loop crossing it twice (`rim_bite`); a pipe ending partway through the
wall, its end cap's circle crossing the wall (`blind_wall`); a box's round
hole crossed at its top rim by a pipe along the box's top face, two rings
on the hole's wall (`hole_rim`, its common two pieces). Turned frames,
against an `XY` cylinder or box: a `TILTX` pipe across the top cap's rim,
two rings, one crossing the cap's circle (`rim_tiltx`); a `LEAN` pipe
entering through the rim and ending inside, a blind hole drilled obliquely
at the rim (`enter_lean`); a `TILT` pipe ending partway through the wall,
its end cap's circle crossing it (`blind_tilt`); a `TILTX` pipe crossing a
box's round hole at its top rim and leaving through the box's bottom
(`hole_rim_tiltx`). No turned pair is perpendicular: in stored turned
frames a perpendicular pair's cap plane is nearly parallel to the other's
generatrices (its section two lines 1e-17 from parallel), which the scan
finds as breakpoints about 1e-17 apart. Every case is declared `solid`.

Before writing, besides the curved generator's per-pair checks (1e-30 and
5 standard errors):

* closed forms: every perpendicular pair of discs by `Perpendicular` (one
  quadrature in `eta` of the product of the two strips, clipped by both
  inputs' caps; its moments, and its areas from the walls' own angles and
  the caps' parts inside the other) within 1e-30, so the three operations'
  volumes, areas and centres; `hole_rim` as the pipe's part below the
  box's top face (a circular segment times the box's width) less its
  common with the hole's cylinder, areas included (the top face's strip
  less the hole's disc between two chords); the turned pairs have none
  (their crossings are not whole, their axes not perpendicular: the
  reference's own checks and Monte Carlo stand for them);
* the curve's class from the ideal and the models' intervals (as
  declared, ends at least 0.1 apart);
* every cap's circle against the other input's cylinders (the opposite of
  the turned generator's check): crossings of the model found by sign
  changes along the circle (8192 samples) refined by bisection, each at
  least 0.1 from the other face's ends along its axis and crossing at a
  slope `|df/dt| >= 0.05` (`f` the model's function over its radius
  squared: no tangency), at least one per pair inside the other face (the
  section's vertex on the cap's circle); a circle without crossings keeps
  off the cylinder by 1e-3 of its radius squared.

The near-coincidence scan must find none. Pairs run in worker processes
(`--workers`).
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
import generate_turned_boolean_fixtures as turned
from generate_curved_boolean_fixtures import FRAMES, at, disc, group, run, square
from generate_procedural_boolean_fixtures import Perpendicular, box_with_hole, cylinder, ops_from_common
from generate_turned_boolean_fixtures import cross, dot, ideal_axes, whole_cylinder

ROOT = curved.ROOT
PREFIX = 'boolean-capped'
STEP = 'S9c.2b.2'

EXACT = turned.EXACT
TURNED = turned.TURNED
# Each pair's curve, from the decisions' intervals.
CURVES = {'rim_pipe': 'rings', 'rim_bite': 'loop', 'blind_wall': 'rings', 'hole_rim': 'rings',
          'rim_tiltx': 'rings', 'enter_lean': 'rings', 'blind_tilt': 'rings', 'hole_rim_tiltx': 'rings'}
# Pairs of two discs with perpendicular axes: closed forms by `Perpendicular`.
PERPENDICULAR = ('rim_pipe', 'rim_bite', 'blind_wall')
# The margins the geometry checks require.
FACE_MARGIN = 0.1
SLOPE = 0.05
CLEAR = 1e-3
SAMPLES = 8192


def cases():
    """S9c.2b.2's pairs (see the module's docstring)."""
    thick = ([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -3.0, 3.0)
    box = [square(0.0, 0.0, 10.0, 10.0), disc(5.0, 5.0, 2.0)]
    out = []
    out += group('rim_pipe', thick, ([disc(0.3, 2.9, 1.0)], at('SIDE', (-5, 0, 0)), 0.0, 10.0),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('rim_bite', thick, ([disc(1.2, 2.6, 1.5)], at('SIDE', (-4, 0, 0)), 0.0, 8.0),
                 {'fuse': 'solid', 'common': 'solid'})
    out += group('blind_wall', thick, ([disc(0.4, 0.3, 1.0)], at('SIDE', (-6, 0, 0)), 0.0, 4.5),
                 {'cut': 'solid', 'common': 'solid'})
    out += group('hole_rim', (box, at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(5.5, 3.9, 0.8)], at('SIDE', (-3, 0, 0)), 0.0, 16.0),
                 {'cut': 'solid', 'common': 'solid'})
    out += group('rim_tiltx', thick, ([disc(0.5, 0.0, 0.8)], at('TILTX', (0, -1.93, 3)), -8.0, 4.0),
                 {'fuse': 'solid', 'cut': 'solid'})
    out += group('enter_lean', thick, ([disc(0.0, 0.0, 0.7)], at('LEAN', (1.95, 0.3, 3)), -2.5, 3.0),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('blind_tilt', thick, ([disc(0.3, 0.0, 0.9)], at('TILT', (0, 2, 0)), 0.4, 5.0),
                 {'cut': 'solid', 'common': 'solid'})
    out += group('hole_rim_tiltx', (box, at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([disc(0.6, 0.0, 0.6)], at('TILTX', (5, 3.1, 4)), -9.0, 3.0),
                 {'fuse': 'solid', 'common': 'solid'})
    return out


def exact_pair(case):
    return all(curved.frame_name(f) in EXACT for f in case.frames)


# ------------------------------------------------------------------ closed forms

def perpendicular(case):
    """The common of a perpendicular pair of discs as `Perpendicular`, in the
    ideal frames: `Z` along the object's axis, `X` along the tool's, `E =
    Z x X` (the common perpendicular), from the object's frame origin."""
    (ca,), (cb,) = turned.cylinders(case.obj), turned.cylinders(case.tool)
    _, _, Z = ideal_axes(case.obj.frame)
    _, _, X = ideal_axes(case.tool.frame)
    assert abs(dot(Z, X)) < mp.mpf(10)**-30, (case.name, 'axes not perpendicular')
    E = cross(Z, X)
    O = tuple(mp.mpf(c) for c in case.obj.frame[:3])
    pa, _ = turned.axis_of(case.obj, ca)
    pb, _ = turned.axis_of(case.tool, cb)
    ra = tuple(pa[i]-O[i] for i in range(3))
    rb = tuple(pb[i]-O[i] for i in range(3))
    za, xb = dot(ra, Z), dot(rb, X)
    return Perpendicular(ca[2], dot(ra, X), dot(ra, E), za+case.obj.start, za+case.obj.end,
                         cb[2], dot(rb, E), dot(rb, Z), xb+case.tool.start, xb+case.tool.end, (O, X, E, Z))


def segment(r, d):
    """The part of a disc of radius `r` below a chord at `d` above its
    centre: its area and its first moment about the centre across the
    chord."""
    r, d = mp.mpf(r), mp.mpf(d)
    return r*r*mp.acos(-d/r)+d*mp.sqrt(r*r-d*d), -2*(r*r-d*d)**mp.mpf(1.5)/3


def band(R, y1, y2):
    """The area of a disc of radius `R` between the chords at `y1 < y2`
    (measured from its centre)."""
    R = mp.mpf(R)
    G = lambda u: u*mp.sqrt(R*R-u*u)+R*R*mp.asin(u/R)
    return G(mp.mpf(y2))-G(mp.mpf(y1))


def hole_rim_forms():
    """`hole_rim`: the box `[0, 10]^2 x [0, 4]` less a hole of radius 2
    about (5, 5), a pipe of radius 0.8 about `y = 5.5`, `z = 3.9` along `x`
    over [-3, 13]. The common is the pipe below `z = 4` between `x = 0` and
    `x = 10` (a segment times 10) less its common `H` with the hole's
    cylinder over `z` in [0, 4]; the box's boundary inside the pipe is its
    two sides' segments, its top face's strip less the hole's disc between
    the strip's chords, and the hole's wall inside the pipe (`H`'s part of
    the hole's boundary less that disc part, its top cap); the pipe's
    boundary inside the box is its wall below `z = 4` over the width less
    its part inside the hole's cylinder."""
    r, zc, yc, top = mp.mpf(0.8), mp.mpf(3.9), mp.mpf(5.5), 4
    d = top-zc
    H = Perpendicular(2, 5, 5, 0, top, r, yc, zc, -3, 13, procedural.WORLD)
    vh, mh, _ = H.volume()
    h_in, p_in = H.areas()
    area, moment = segment(r, d)
    vs = 10*area
    ms = (vs*5, vs*yc, 10*(zc*area+moment))
    vc = vs-vh
    mc = tuple(ms[i]-mh[i] for i in range(3))
    w = mp.sqrt(r*r-d*d)
    disc_part = band(2, yc-w-5, yc+w-5)
    a_in = 2*area+(10*2*w-disc_part)+(h_in-disc_part)
    b_in = 10*(2*mp.pi*r-2*r*mp.acos(d/r))-p_in
    A = box_with_hole(10, top, 5, 5, 2)
    B = cylinder(r, -3, 13, (5, yc, zc))
    return ops_from_common(A, B, (vc, mc, a_in, b_in))


def closed_forms():
    """{pair: {op: (volume, moments, area)}} in the ideal frames, with
    whether the pair lies in exact frames."""
    firsts = {}
    for c in cases():
        firsts.setdefault(c.pair_name, c)
    out = {}
    for name in PERPENDICULAR:
        c = firsts[name]
        P = perpendicular(c)
        vc, mc, _ = P.volume()
        a_in, b_in = P.areas()
        out[name] = (ops_from_common(whole_cylinder(c.obj), whole_cylinder(c.tool), (vc, mc, a_in, b_in)),
                     exact_pair(c))
    out['hole_rim'] = (hole_rim_forms(), True)
    return out


# ------------------------------------------------------------------ caps' circles

def cap_crossings(case):
    """Every cap's circle of one input against every cylinder of the other,
    in the stored models (floats). Returns (crossings, the least `|f|` of a
    circle without crossings), each crossing `(which input's cap, the cap's
    height, the circle's angle, the signed margin inside the other face
    along its axis (negative outside), |df/dt|)`, `f` the other cylinder's
    function over its radius squared."""
    A, B = ref.Prism(case.obj), ref.Prism(case.tool)
    crossings, least_clear = [], math.inf
    for tag, P, Q in (('obj', A, B), ('tool', B, A)):
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

                def local(X):
                    d = [X[i]-fq.o[i] for i in range(3)]
                    return [sum(fq.inv[k][i]*d[i] for i in range(3)) for k in range(3)]

                for el2 in Q.profile.elements:
                    if el2.kind != 'circle':
                        continue
                    qx, qy, qr = float(el2.c[0]), float(el2.c[1]), float(el2.r)

                    def f(t):
                        u, v, _ = local(point(t))
                        return ((u-qx)**2+(v-qy)**2-qr*qr)/(qr*qr)

                    ts = [2*math.pi*k/SAMPLES for k in range(SAMPLES+1)]
                    fs = [f(t) for t in ts]
                    roots = []
                    for k in range(SAMPLES):
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
                        margin = min(w-fq.lo, fq.hi-w)
                        slope = abs(f(t+1e-6)-f(t-1e-6))/2e-6
                        crossings.append((tag, h, t, margin, slope))
    return crossings, least_clear


def geometry_checks(case):
    """The pair's class (ideal and model intervals) and its caps' circles."""
    ideal, model = turned.intervals(case)
    return turned.classify(ideal), turned.classify(model), cap_crossings(case)


def geometry_job(case):
    return case.pair_name, geometry_checks(case)


def crossed(crossings):
    """The caps whose circles cross the other cylinder on its face:
    `obj lo`, `tool hi` and so on, with the count."""
    out = {}
    for tag, h, _, margin, _ in crossings:
        if margin > 0:
            key = f'{tag} {h:g}'
            out[key] = out.get(key, 0)+1
    return out


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
    for name, (forms, exact) in closed_forms().items():
        _, _, _, res, _, _, _ = by_pair[name]
        key = 'closed_forms_exact_frames' if exact else 'closed_forms_turned_frames'
        for op in ref.OPS:
            if res[op][0] == 0:
                continue
            note(key, procedural.deviation(res[op], forms[op]))
    for name, rows, _, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    margins = {'class_gap': mp.inf, 'cap_crossing_inside': math.inf, 'cap_crossing_outside': math.inf,
               'cap_crossing_slope': math.inf, 'cap_circle_clear': math.inf, 'cap_crossings_on_faces': 0,
               'cap_crossings': 0}
    for name, ((kind, gap), (mkind, mgap), (crossings, clear)) in geometry.items():
        want = CURVES[name]
        assert kind == want and mkind == want, (name, kind, mkind, want)
        assert gap > mp.mpf(1)/10 and mgap > mp.mpf(1)/10, (name, 'interval ends close', gap, mgap)
        margins['class_gap'] = min(margins['class_gap'], gap, mgap)
        inside = [c for c in crossings if c[3] > 0]
        assert inside, (name, "no cap's circle crosses the other cylinder on its face")
        for tag, h, t, margin, slope in crossings:
            assert abs(margin) > FACE_MARGIN, (name, tag, h, 'a crossing near the other face\'s end', margin)
            assert slope > SLOPE, (name, tag, h, 'a crossing nearly tangent', slope)
            key = 'cap_crossing_inside' if margin > 0 else 'cap_crossing_outside'
            margins[key] = min(margins[key], abs(margin))
            margins['cap_crossing_slope'] = min(margins['cap_crossing_slope'], slope)
        if clear < math.inf:
            assert clear > CLEAR, (name, 'a cap circle nearly tangent to the other cylinder', clear)
            margins['cap_circle_clear'] = min(margins['cap_circle_clear'], clear)
        margins['cap_crossings_on_faces'] += len(inside)
        margins['cap_crossings'] += len(crossings)
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
        assert all(f in EXACT for f in frames) or any(f in TURNED for f in frames), c.name
        _, _, _, na = stored_axes(c.obj.frame)
        _, _, _, nb = stored_axes(c.tool.frame)
        assert cross(na, nb) != (0.0, 0.0, 0.0), f'{c.name}: parallel axes'
        assert turned.cylinders(c.obj) and turned.cylinders(c.tool), f'{c.name}: not two cylinders'
        assert c.kind in ('solid', 'empty', 'degenerate') and (c.kind == 'degenerate') == (c.reason is not None)


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
    near = [f'{name}: near coincidences {n}' for name, _, _, _, _, n, _ in results if n]
    if near:
        raise SystemExit('\n'.join(near))
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames')
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('geometry:', ', '.join(f'{k} {v}' if isinstance(v, int) else f'{k} {float(v):.3g}'
                                 for k, v in margins.items()))
    for name, ((kind, gap), (mkind, mgap), (crossings, clear)) in geometry.items():
        on = ', '.join(f'{k} x{v}' for k, v in sorted(crossed(crossings).items()))
        print(f'  {name}: {kind} (ends {mp.nstr(gap, 3)} apart), cap circles crossing on the other face: {on}; '
              f'{len(crossings)} crossings in all, clear by {clear:.3g}')
    by_pair = {r[0]: r for r in results}
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
