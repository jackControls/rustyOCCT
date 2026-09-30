#!/usr/bin/env python3
"""Fixtures for S9e.1 of REVIEW_NOTES.md: a Boolean's result given to
another Boolean. The first Boolean is S9c.1's (two prisms of line, arc and
circle profiles in frames of different axes, at least one arc), its result
one solid; the second takes that result as its object, or as its tool
(`swapped`), with a third prism of lines, arcs and circles, every pair of
the three prisms' faces meeting in lines, circles or ellipses.

`boolean-chained-cases.txt` lists each case in the Boolean protocol with a
`then` row (`identity_reference.encode_chained_case`); only frames whose
stored axes the kernel's `Frame3::new` gives bit for bit as `stored_axes`
are used (`boolean-chained-frames.tsv` records them): the exact `XY` and
`SIDE` and the turned `TILT` (`generate_curved_boolean_fixtures.FRAMES`).

`boolean-chained-expected.tsv` gives per case, from
`chained_boolean_reference.py`:

* `expect KIND S9e.1`: the declared outcome (`solid`: one or more solids;
  `empty`: none; `degenerate`: a tangency between the first result and the
  third prism or an edge of one on the other's face, `Degenerate` in the
  decisions), then for a degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

The classes: a box less a tilted hole halved by a box whose wall holds the
hole's axis (the wall crossing the hole's ellipse edges); that result as
the tool of a slab crossing the hole (the cut two solids); a box fused with
a tilted pin, then cut by a box whose wall holds the pin's axis and whose
bottom crosses the box; the holed box under a box on its top face
(coplanar faces, edges crossing on them; the common empty); a box less a
groove across its top (its top face in two faces of the result), then
drilled through one of them; the holed box against a column of the hole's
frame whose cap cuts the hole in a circle; a quarter cylinder (a cylinder's
common with a box) against a coaxial cylinder (closed forms); the fused
pin's result as the tool of a slab it crosses; declared `degenerate`: a
third cylinder tangent to the groove's wall along a generatrix, a third box
whose edge lies on the groove's wall.

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size): closed forms (the groove and
the quarter cylinder, exact frames) within 1e-30; S9c.1's pair reference
for the first result `X = A op1 B` in the identities `V(X u C) + V(X n C)
= V(X) + V(C)`, `V(X - C) = V(X) - V(X n C)` (swapped: `V(C - X) = V(C) -
V(X n C)`), moments likewise, and `area(X u C) + area(X n C) = area(X) +
area(C)` where no face of the third lies on a face of the others, within
1e-30; every face's classes summing to its closed-form area; Monte-Carlo
estimates (200,000 uniform points per chain) of every volume and centre
within 5 standard errors. A scan for near coincidences (a class of a face or
a result's volume positive but below 1e-9 of the case's size, slicing
breakpoints closer than that) must find none outside the declared cases.
Chains run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Case, encode_chained_case
from curve_surface_reference import stored_axes
import curved_boolean_reference as cref
import chained_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, FloatPrism, at, disc, square

ROOT = Path(__file__).resolve().parents[1]
FIRST_OPERATION = 93
THIRD = 94
SECOND_OPERATION = 95


def prism(boundaries, frame, start, end, op):
    return Case('', 1e-7, op, frame, start, end, boundaries)


class Chained:
    def __init__(self, name, op1, obj, tool, then, third, swapped, kind='solid', reason=None):
        self.name, self.op1, self.then, self.swapped = name, op1, then, swapped
        self.kind, self.reason = kind, reason
        self.chain_name = name.rsplit('_', 1)[0]
        self.obj = prism(*obj, op=91)
        self.tool = prism(*tool, op=92)
        self.third = prism(*third, op=THIRD)
        for c in (self.obj, self.tool, self.third):
            c.name = name
        self.frames = (obj[1], tool[1], third[1])

    def encode(self):
        return encode_chained_case(self.obj, self.op1, self.tool, FIRST_OPERATION, self.then, self.third,
                                   SECOND_OPERATION, self.swapped)


def group(name, first, op1, third, ops, swapped=False):
    """Cases of one chain: `first` the first Boolean's object and tool,
    `ops` maps each second operation to its declared kind (or `(kind,
    reason)`)."""
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Chained(f'{name}_{op}', op1, first[0], first[1], op, third, swapped, kind, reason))
    return out


SOLID3 = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
TANGENT = 'a third cylinder tangent to the first result\'s cylinder along a generatrix'
EDGE_ON_FACE = 'an edge of the third prism on the first result\'s cylindrical face'

BOX = [square(0.0, 0.0, 10.0, 10.0)]
# A box less a tilted hole (the hole's cylinder leaves through the top and
# bottom faces in ellipses).
HOLED = ((BOX, at('XY', (0, 0, 0)), 0.0, 5.0), ([disc(0.0, 0.0, 1.5)], at('TILT', (5, 3, 0)), -2.0, 8.0))
# A box fused with a tilted pin sticking out of its top and bottom.
PINNED = ((BOX, at('XY', (0, 0, 0)), 0.0, 5.0), ([disc(0.0, 0.0, 0.8)], at('TILT', (5, 5, 0)), -4.0, 10.0))
# A box less a groove of radius 5/4 across its top along x (SIDE: u = y,
# v = z): the top face left in two faces.
GROOVE_R = 1.25
GROOVED = (([square(0.0, 0.0, 10.0, 6.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
           ([disc(3.0, 4.0, GROOVE_R)], at('SIDE', (0, 0, 0)), -1.0, 11.0))
# A cylinder of radius 3 in common with a box: a quarter cylinder.
QUARTER = (([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
           ([square(0.0, -1.0, 5.0, 6.0)], at('SIDE', (0, 0, 0)), 0.0, 5.0))


def cases():
    out = []
    out += group('hole_halved', HOLED, 'cut',
                 ([square(5.0, -1.0, 12.0, 11.0)], at('XY', (0, 0, 0)), -1.0, 6.0), SOLID3)
    out += group('hole_tool', HOLED, 'cut',
                 ([square(4.0, -2.0, 12.0, 12.0)], at('XY', (0, 0, 0)), 1.0, 4.0), SOLID3, swapped=True)
    out += group('pin_step', PINNED, 'fuse',
                 ([square(5.0, -5.0, 20.0, 15.0)], at('XY', (0, 0, 0)), 4.0, 12.0), SOLID3)
    out += group('hole_capped', HOLED, 'cut',
                 ([square(2.0, 2.0, 8.0, 8.0)], at('XY', (0, 0, 0)), 5.0, 7.5),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'empty'})
    out += group('groove_drilled', GROOVED, 'cut',
                 ([disc(5.0, 0.75, 0.5)], at('XY', (0, 0, 0)), -1.0, 5.0), SOLID3)
    out += group('hole_column', HOLED, 'cut',
                 ([square(-3.0, -3.0, 3.0, 3.0)], at('TILT', (5, 3, 0)), 3.0, 9.0), SOLID3)
    out += group('quarter_bored', QUARTER, 'common',
                 ([disc(0.0, 0.0, 1.5)], at('XY', (0, 0, 0)), 1.0, 7.0), SOLID3)
    out += group('pin_slab', PINNED, 'fuse',
                 ([square(-2.0, -2.0, 12.0, 14.0)], at('XY', (0, 0, 0)), 6.0, 7.0), SOLID3, swapped=True)
    out += group('groove_tangent', GROOVED, 'cut',
                 ([disc(3.0, 4.0-GROOVE_R-0.5, 0.5)], at('SIDE', (0, 0, 0)), -1.0, 12.0),
                 {op: ('degenerate', TANGENT) for op in ('fuse', 'cut', 'common')})
    # The groove's wall holds the line y = 2.25, z = 3 (0.75^2 + 1^2 =
    # 1.25^2): the third box's top edge along x lies on it.
    out += group('groove_edge', GROOVED, 'cut',
                 ([square(5.0, -1.0, 8.0, 2.25)], at('XY', (0, 0, 0)), -1.0, 3.0),
                 {op: ('degenerate', EDGE_ON_FACE) for op in ('fuse', 'cut', 'common')})
    return out


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert stored_axes(c.obj.frame)[1:] != stored_axes(c.tool.frame)[1:], \
            f'{c.name}: the first Boolean in frames with equal axes (S9a)'
        assert any(b.circle is not None or b.segments is not None for b in c.obj.boundaries+c.tool.boundaries), \
            f'{c.name}: no arc in the first Boolean'
        for frame in c.frames:
            assert frame_name(frame) in FRAMES, c.name


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


# ------------------------------------------------------------------ closed forms

def box_moments(x0, x1, y0, y1, z0, z1):
    v = (x1-x0)*(y1-y0)*(z1-z0)
    return v, (v*(x0+x1)/2, v*(y0+y1)/2, v*(z0+z1)/2)


def closed_forms():
    """{(chain, second operation): (volume, area, centre)} (mpmath)."""
    pi = mp.pi
    out = {}

    def add(key, parts, area):
        vol = sum(p[0] for p in parts)
        mom = [sum(p[1][i] for p in parts) for i in range(3)]
        out[key] = (vol, area, tuple(m/vol for m in mom))

    # The groove: box [0, 10] x [0, 6] x [0, 4] less the half disc of radius
    # r about (y, z) = (3, 4) below z = 4, along x in [0, 10].
    r = mp.mpf(GROOVE_R)
    half = pi*r*r/2
    groove = (-10*half, (-10*half*5, -10*half*3, -10*half*(4-4*r/(3*pi))))
    box = box_moments(0, 10, 0, 6, 0, 4)
    area_x = 2*(60+40+24)-10*2*r-pi*r*r+pi*r*10
    # The drill: radius 1/2 about (5, 0.75), z in [-1, 5]; in the result's
    # material for z in [0, 4].
    a = mp.mpf(1)/2
    inside = pi*a*a*4
    drill_in = (inside, (inside*5, inside*mp.mpf(0.75), inside*2))
    drill_all = pi*a*a*6
    stubs = (drill_all-inside, ((drill_all-inside)*5, (drill_all-inside)*mp.mpf(0.75), (drill_all-inside)*2))
    add(('groove_drilled', 'common'), [drill_in], 2*pi*a*a+2*pi*a*4)
    add(('groove_drilled', 'cut'), [box, groove, (-inside, tuple(-m for m in drill_in[1]))],
        area_x-2*pi*a*a+2*pi*a*4)
    add(('groove_drilled', 'fuse'), [box, groove, stubs], area_x+2*pi*a*2)
    # The quarter cylinder of radius 3 over [0, 5] and the coaxial cylinder
    # of radius 3/2 over [1, 7].
    R, h = mp.mpf(3), mp.mpf(5)
    b, lo, hi = mp.mpf(3)/2, mp.mpf(1), mp.mpf(7)

    def quarter(rad, z0, z1, sign=1):
        v = sign*pi*rad*rad*(z1-z0)/4
        g = 4*rad/(3*pi)
        return (v, (v*g, v*g, v*(z0+z1)/2))

    def whole(rad, z0, z1, sign=1):
        v = sign*pi*rad*rad*(z1-z0)
        return (v, (0*v, 0*v, v*(z0+z1)/2))

    area_q = 2*pi*R*R/4+2*pi*R*h/4+2*R*h
    area_c = 2*pi*b*b+2*pi*b*(hi-lo)
    common = quarter(b, lo, h)
    area_common = 2*pi*b*b/4+2*pi*b*(h-lo)/4+2*b*(h-lo)
    add(('quarter_bored', 'common'), [common], area_common)
    add(('quarter_bored', 'cut'), [quarter(R, 0, h), quarter(b, lo, h, -1)],
        area_q-2*b*(h-lo)+2*pi*b*(h-lo)/4)
    add(('quarter_bored', 'fuse'), [quarter(R, 0, h), whole(b, lo, hi), quarter(b, lo, h, -1)],
        area_q+area_c-area_common)
    return out


# ------------------------------------------------------------------ Monte Carlo

def monte_carlo(chain, n, seed):
    """Volume and centre estimates of the three second operations and their
    standard errors, from `n` uniform points in the prisms' box."""
    fps = [FloatPrism(p) for p in chain.prisms]
    lo, hi = [], []
    for p in chain.prisms:
        b0, b1 = p.bounds()
        lo.append([float(c) for c in b0])
        hi.append([float(c) for c in b1])
    box0 = [min(l[i] for l in lo) for i in range(3)]
    box1 = [max(h[i] for h in hi) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    fns = {op: chain.fn(op) for op in ref.OPS}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        m = [fp.contains(X) for fp in fps]
        for op in ref.OPS:
            if fns[op](*m):
                e = acc[op]
                e[0] += 1
                for i in range(3):
                    e[1][i] += X[i]
                    e[2][i] += X[i]*X[i]
    out = {}
    for op, (k, s1, s2) in acc.items():
        p = k/n
        vol = vbox*p
        sv = vbox*(max(p*(1-p), 1e-300)/n)**0.5
        if k > 1:
            c = [s/k for s in s1]
            sc = [(max(s2[i]/k-c[i]*c[i], 0)/k)**0.5 for i in range(3)]
        else:
            c, sc = None, None
        out[op] = (vol, sv, c, sc)
    return out


# ------------------------------------------------------------------ the chains' work

def evaluate(job):
    """One chain: the reference's rows for its second operations and every
    check's deviations (run in a worker process)."""
    name, obj, tool, third, op1, swapped, ops, mc_n = job
    chain = ref.Chain(obj, tool, third, op1, swapped)
    res = {op: chain.result(op) for op in ref.OPS}
    rows = {op: ref.rows(chain, op) for op in ops}
    size = chain.size
    vols = chain.volumes()
    checks = {}
    # S9c.1's pair reference for the first result.
    pair = cref.Pair(obj, tool)
    vx, mx = pair.volumes()[op1]
    ax = pair.area(op1)
    C = chain.prisms[2]
    vc, mc, ac = C.measures()
    vf, mf = vols['fuse']
    vn, mn = vols['common']
    vt, mt = vols['cut']
    dev = abs(vf+vn-vx-vc)/size**3
    for i in range(3):
        dev = max(dev, abs(mf[i]+mn[i]-mx[i]-mc[i])/size**4)
    if swapped:
        dev = max(dev, abs(vt-(vc-vn))/size**3, max(abs(mt[i]-(mc[i]-mn[i])) for i in range(3))/size**4)
    else:
        dev = max(dev, abs(vt-(vx-vn))/size**3, max(abs(mt[i]-(mx[i]-mn[i])) for i in range(3))/size**4)
    checks['pair_identities'] = dev
    # The first result's pair counts one solid (the second's argument).
    assert pair.result(op1)[0] == 1, f'{name}: the first result is not one solid'
    shared_with_third = False
    for i, f, cls, _ in chain.face_areas():
        others = [j for j in range(3) if j != i]
        for classes, a in cls.items():
            if a <= 0:
                continue
            for o, c in zip(others, classes):
                if c in ('same', 'opp') and 2 in (i, o):
                    shared_with_third = True
    if not shared_with_third:
        checks['area_identity'] = abs(chain.area('fuse')+chain.area('common')-ax-ac)/size**2
    # Faces: classes against closed-form areas.
    fdev = mp.mpf(0)
    for i, f, cls, _ in chain.face_areas():
        fdev = max(fdev, abs(sum(cls.values())-f.closed_area())/size**2)
    checks['face_classes'] = fdev
    # Monte Carlo.
    mc_ = monte_carlo(chain, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mc_[op]
        v = float(vols[op][0])
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/vols[op][0]) for x in vols[op][1]]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(chain)
    stats = {'volume_quadrature': chain.slicing.quad_error/size**4,
             'face_quadrature': max(sw.quad_error for _, _, _, sw in chain.face_areas())/size**2,
             'breaks': len(chain.slicing.breaks)}
    return name, rows, res, checks, near, stats


def near_coincidences(chain):
    size = chain.size
    out = []
    for i, f, cls, _ in chain.face_areas():
        for c, v in cls.items():
            if 0 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{"ABC"[i]} {cref.describe(f)} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = chain.volumes()[op][0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    b = chain.slicing.breaks
    for p, q in zip(b, b[1:]):
        if q-p < mp.mpf(10)**-9*size:
            out.append(f'breakpoints {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    by_chain = {r[0]: r for r in results}
    for (name, op), (V, S, Cn) in closed_forms().items():
        _, _, res, _, _, _ = by_chain[name]
        n, vol, area, centre = res[op]
        size = max(abs(V)**(mp.mpf(1)/3), 1)
        note('closed_forms', max(abs(vol-V)/abs(V), abs(area-S)/abs(S),
                                 max(abs(centre[i]-Cn[i]) for i in range(3))/size))
    for name, rows, _, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    return worst, covered


def generate(results):
    by_chain = {r[0]: r for r in results}
    blocks = []
    out = ['# case\trow (S9e.1, chained_boolean_reference.py: expect KIND S9e.1, reason TEXT for a degenerate '
           'case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, _, _, _ = by_chain[case.chain_name]
        row = rows[case.then]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9e.1')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (obj, tool, third) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in cases():
        for label, c in (('obj', case.obj), ('tool', case.tool), ('third', case.third)):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\t{label} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-chained-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-chained-expected.tsv': '\n'.join(out)+'\n',
            'boolean-chained-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n):
    chains = {}
    for c in cases():
        chains.setdefault(c.chain_name, [c.obj, c.tool, c.third, c.op1, c.swapped, []])[5].append(c.then)
    return [(name, obj, tool, third, op1, swapped, ops, mc_n)
            for name, (obj, tool, third, op1, swapped, ops) in chains.items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(6, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per chain')
    parser.add_argument('--only', help='one chain, printed (no files)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in run(todo, evaluate, 1):
            print(r)
        return
    results = run(todo, evaluate, args.workers)
    declared = {c.chain_name: c.kind == 'degenerate' for c in listed}
    for name, _, _, _, near, _ in results:
        if near and not declared[name]:
            raise SystemExit(f'{name}: near coincidences {near}')
    worst, covered = reference_checks(results)
    limits = {'closed_forms': 1e-30, 'pair_identities': 1e-30, 'area_identity': 1e-30, 'face_classes': 1e-30,
              'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-30, 'face_quadrature_estimate': 1e-30}
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
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.then == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    for name, _, _, _, near, _ in results:
        if near:
            print('declared near coincidences:', name, '; '.join(near[:4]))


if __name__ == '__main__':
    main()
