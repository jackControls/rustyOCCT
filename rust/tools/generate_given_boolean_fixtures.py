#!/usr/bin/env python3
"""Fixtures for S9e.2 of REVIEW_NOTES.md: a Boolean's result given to
another Boolean where the result's stored topology comes from another
assembly or holds part of its construction's region: a stack with an arc
wall (S9a.2), an S9b.1 result of two line prisms given with arcs, and one
solid of a first result of several. The second Boolean takes the given
solid as its object, or as its tool (`swapped`), with a third prism of
lines, arcs and circles, every pair of the three prisms' faces meeting in
lines, circles or ellipses (S9c.1's pairs).

`boolean-given-cases.txt` lists each case in the Boolean protocol with a
`then` row (`identity_reference.encode_chained_case`); a first result of
several solids picks one by a point inside it (`then OP ID [swapped] solid X
Y Z`). Only frames whose stored axes the kernel's `Frame3::new` gives bit
for bit as `stored_axes` are used (`boolean-given-frames.tsv` records
them): the exact `XY`, `SIDE` and `DOWN`, the turned `TILT` and `R125`
(`generate_curved_boolean_fixtures.FRAMES`).

`boolean-given-expected.tsv` gives per case, from
`chained_boolean_reference.py`:

* `expect KIND S9e.2 CLASS`: the declared outcome (`solid`: one or more
  solids; `empty`: none; `degenerate`: `Degenerate` in the decisions) and
  the class of the given solid (`stack`, `polyhedral`, `several`), then for
  a degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

The chains: DRAW's `bcut_simple/L3` to `L6` in the rollex's geometry (a
disc of radius 60 less a pocket of radius 40 across its rim over the top 6
of its height 20, a stack; then a cylinder of radius 30 whose bottom lies
on the pocket's floor, crossing the pocket's wall and the disc's rim, in the
turned frame `DOWN` as `L3` and `L4`, in the stack's frame as `L5` and
`L6`); a box with a boss (a stack) bored coaxially (closed forms), sliced by
a tilted slab, and as the tool of a tilted box; S9b.1 results: a box fused
with a turned box, bored through both; a box less a tilted box (a slanted
top), bored through the slanted face; several solids: a box severed by a
cylinder across it (S9e.1's class, two solids), the lower one drilled and
as the tool of a box; a box severed by a tilted slab (two S9b.1 solids),
the lower one bored through its slanted face; declared `degenerate`: a
third cylinder tangent to the rollex's rim along a generatrix, a third box
whose wall is tangent to the other solid of the severed box (the
arrangement holds both solids: the decisions' conservative refusal).

For a first result of several solids the reference adds a selector `D`, a
box holding the picked solid and no part of another
(`chained_boolean_reference.Chain(..., selector=D)`, the chain `op2(op1(A,
B) and D, C)`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size): closed forms (the bored boss,
exact frames) within 1e-30; for the first result `X` (S9c.1's pair
reference `curved_boolean_reference.Pair`, or for a picked solid the
three-prism chain `op1(A, B) and D`, one solid, while the pair counts
several) the identities `V(X u C) + V(X n C) = V(X) + V(C)`, `V(X - C) =
V(X) - V(X n C)` (swapped: `V(C - X) = V(C) - V(X n C)`), moments
likewise, and `area(X u C) + area(X n C) = area(X) + area(C)` where no face
of the third lies on a face of the others, within 1e-30; every face's
classes summing to its closed-form area; the selector's faces bounding no
part of the first result (zero area) and the pick point inside the picked
solid; Monte-Carlo estimates (200,000 uniform points per chain) of every
volume and centre within 5 standard errors. A scan for near coincidences (a
class of a face or a result's volume positive but below 1e-9 of the case's
size, slicing breakpoints closer than that) must find none outside the
declared cases. Chains run in worker processes (`--workers`).
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
CLASSES = ('stack', 'polyhedral', 'several')


def prism(boundaries, frame, start, end, op):
    return Case('', 1e-7, op, frame, start, end, boundaries)


class Given:
    def __init__(self, name, klass, op1, obj, tool, then, third, swapped, kind='solid', reason=None,
                 pick=None, selector=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.then, self.swapped = name, klass, op1, then, swapped
        self.kind, self.reason, self.pick = kind, reason, pick
        self.chain_name = name.rsplit('_', 1)[0]
        self.obj = prism(*obj, op=91)
        self.tool = prism(*tool, op=92)
        self.third = prism(*third, op=THIRD)
        # The reference's selector (not in the protocol: the pick point is).
        self.selector = None if selector is None else prism(*selector, op=96)
        for c in (self.obj, self.tool, self.third):
            c.name = name
        self.frames = (obj[1], tool[1], third[1])

    def encode(self):
        return encode_chained_case(self.obj, self.op1, self.tool, FIRST_OPERATION, self.then, self.third,
                                   SECOND_OPERATION, self.swapped, self.pick)


def group(name, klass, first, op1, third, ops, swapped=False, pick=None, selector=None):
    """Cases of one chain: `first` the first Boolean's object and tool,
    `ops` maps each second operation to its declared kind (or `(kind,
    reason)`)."""
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Given(f'{name}_{op}', klass, op1, first[0], first[1], op, third, swapped, kind, reason,
                         pick, selector))
    return out


SOLID3 = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
TANGENT = 'a third cylinder tangent to the stack\'s cylindrical rim along a generatrix'
TANGENT_OTHER = ('a third box\'s wall tangent to the other solid of the first result along a generatrix '
                 '(the arrangement holds every solid of it)')

# DRAW's rollex (bcut_simple/L3 to L6): a disc of radius 60 about (0, 60)
# less a pocket of radius 40 about (10, 20) over z in [14, 20] (the pocket
# crosses the rim): a stack.
ROLLEX = (([disc(0.0, 60.0, 60.0)], at('XY', (0, 0, 0)), 0.0, 20.0),
          ([disc(10.0, 20.0, 40.0)], at('XY', (0, 0, 0)), 14.0, 20.0))
# L3 and L4: the cylinder of radius 30 about (50, 40) over z in [14, 23] in
# the turned frame DOWN (local y = -y: the centre (50, -40) there).
ROLLEX_DOWN = ([disc(50.0, -40.0, 30.0)], at('DOWN', (0, 0, 23)), 0.0, 9.0)
# L5 and L6: the same cylinder in the stack's frame.
ROLLEX_XY = ([disc(50.0, 40.0, 30.0)], at('XY', (0, 0, 0)), 14.0, 23.0)
# A box with a boss of radius 2 over z in [2, 8] (above the box's top at
# 4): a stack.
BOSS = (([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
        ([disc(5.0, 5.0, 2.0)], at('XY', (0, 0, 0)), 2.0, 8.0))
# A slab above a tilted plane (TILT: w = 0.6 y + 0.8 (z - 10)) crossing
# the boss in an ellipse and the box's top.
TILTED_SLAB = ([square(-2.0, -15.0, 12.0, 15.0)], at('TILT', (0, 0, 10)), 0.0, 8.0)
# A box fused with a box turned about z (R125), S9b.1.
TURNED = (([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
          ([square(0.0, 0.0, 6.0, 6.0)], at('R125', (4, -2, 0)), 2.0, 8.0))
# A box less a tilted box above w = 0 (TILT at (0, 0, 7)): a slanted top,
# S9b.1.
SLANTED = (([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 6.0),
           ([square(-2.0, -15.0, 12.0, 15.0)], at('TILT', (0, 0, 7)), 0.0, 8.0))
# A box severed by a cylinder of radius 3 along x about (y, z) = (5, 2.5)
# (SIDE: u = y, v = z): two solids, y <= 5 - sqrt(2.75) and y >= 5 +
# sqrt(2.75) (S9e.1's class).
SEVERED = (([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 5.0),
           ([disc(5.0, 2.5, 3.0)], at('SIDE', (0, 0, 0)), -1.0, 11.0))
SEVERED_PICK = (5.0, 1.0, 2.5)
SEVERED_SELECTOR = ([square(-1.0, -1.0, 11.0, 5.0)], at('XY', (0, 0, 0)), -1.0, 6.0)
# A box severed by a tilted slab 3 <= w <= 5 (TILT at the origin: w = 0.6 y
# + 0.8 z): two S9b.1 solids.
SPLIT = (([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
         ([square(-2.0, -15.0, 12.0, 15.0)], at('TILT', (0, 0, 0)), 3.0, 5.0))
SPLIT_PICK = (5.0, 1.0, 1.0)
SPLIT_SELECTOR = ([square(-5.0, -20.0, 15.0, 20.0)], at('TILT', (0, 0, 0)), -10.0, 4.0)


def cases():
    out = []
    out += group('rollex_turned', 'stack', ROLLEX, 'cut', ROLLEX_DOWN, SOLID3)
    out += group('rollex_flat', 'stack', ROLLEX, 'cut', ROLLEX_XY, SOLID3)
    out += group('boss_bored', 'stack', BOSS, 'fuse',
                 ([disc(5.0, 5.0, 1.0)], at('XY', (0, 0, 0)), 1.0, 10.0), SOLID3)
    out += group('boss_sliced', 'stack', BOSS, 'fuse', TILTED_SLAB, SOLID3)
    out += group('boss_tool', 'stack', BOSS, 'fuse',
                 ([square(-1.0, -1.0, 11.0, 7.0)], at('TILT', (0, 0, 4)), -2.0, 2.0), SOLID3, swapped=True)
    out += group('turned_bored', 'polyhedral', TURNED, 'fuse',
                 ([disc(5.0, 5.0, 1.5)], at('XY', (0, 0, 0)), -1.0, 9.0), SOLID3)
    out += group('slanted_bored', 'polyhedral', SLANTED, 'cut',
                 ([disc(4.0, 6.0, 1.0)], at('XY', (0, 0, 0)), -1.0, 8.0), SOLID3)
    out += group('severed_drilled', 'several', SEVERED, 'cut',
                 ([disc(5.0, 1.0, 0.5)], at('XY', (0, 0, 0)), -1.0, 7.0), SOLID3,
                 pick=SEVERED_PICK, selector=SEVERED_SELECTOR)
    out += group('severed_tool', 'several', SEVERED, 'cut',
                 ([square(3.0, -1.0, 7.0, 4.0)], at('XY', (0, 0, 0)), 1.0, 7.0), SOLID3, swapped=True,
                 pick=SEVERED_PICK, selector=SEVERED_SELECTOR)
    out += group('split_bored', 'several', SPLIT, 'cut',
                 ([disc(5.0, 1.5, 1.0)], at('XY', (0, 0, 0)), -1.0, 5.5), SOLID3,
                 pick=SPLIT_PICK, selector=SPLIT_SELECTOR)
    # The rim's circle about (0, 60) of radius 60 holds (0, 0); the third
    # cylinder of radius 10 about (0, 10) touches it there from inside,
    # below the pocket.
    out += group('rollex_tangent', 'stack', ROLLEX, 'cut',
                 ([disc(0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 4.0, 12.0),
                 {op: ('degenerate', TANGENT) for op in ('fuse', 'cut', 'common')})
    # The severing cylinder reaches y = 8 at z = 2.5: the third box's wall
    # y = 8 touches the upper solid there, crossing the lower one.
    out += group('severed_tangent', 'several', SEVERED, 'cut',
                 ([square(3.0, 1.0, 7.0, 8.0)], at('XY', (0, 0, 0)), 1.0, 4.0),
                 {op: ('degenerate', TANGENT_OTHER) for op in ('fuse', 'cut', 'common')},
                 pick=SEVERED_PICK, selector=SEVERED_SELECTOR)
    return out


def lines_only(c):
    return all(b.circle is None and b.segments is None for b in c.boundaries)


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        same = stored_axes(c.obj.frame)[1:] == stored_axes(c.tool.frame)[1:]
        if c.klass == 'stack':
            assert same and (c.obj.start, c.obj.end) != (c.tool.start, c.tool.end), \
                f'{c.name}: a stack is a Boolean in one frame of different heights (S9a.2)'
        if c.klass == 'polyhedral':
            assert not same and lines_only(c.obj) and lines_only(c.tool), \
                f'{c.name}: an S9b.1 result is of line prisms in frames with different axes'
        assert (c.klass == 'several') == (c.pick is not None) == (c.selector is not None), c.name
        assert not lines_only(c.obj) or not lines_only(c.tool) or not lines_only(c.third), \
            f'{c.name}: no arc among the three prisms'
        for frame in c.frames:
            assert frame_name(frame) in FRAMES, c.name


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


# ------------------------------------------------------------------ closed forms

def closed_forms():
    """{(chain, second operation): (volume, area, centre)} (mpmath): the
    boss (a box [0, 10]^2 x [0, 4] with a disc of radius 2 about (5, 5) up
    to z = 8) bored by the coaxial cylinder of radius 1 over z in [1, 10]."""
    pi = mp.pi
    out = {}

    def add(key, parts, area):
        vol = sum(p[0] for p in parts)
        mom = [sum(p[1][i] for p in parts) for i in range(3)]
        out[key] = (vol, area, tuple(m/vol for m in mom))

    def cyl(r, z0, z1, sign=1):
        v = sign*pi*r*r*(z1-z0)
        return (v, (v*5, v*5, v*(z0+z1)/2))

    box = (mp.mpf(400), (mp.mpf(2000), mp.mpf(2000), mp.mpf(800)))
    area_x = 360+16*pi
    add(('boss_bored', 'common'), [cyl(1, 1, 8)], 2*pi+14*pi)
    add(('boss_bored', 'cut'), [box, cyl(2, 4, 8), cyl(1, 1, 8, -1)], area_x-pi+14*pi+pi)
    add(('boss_bored', 'fuse'), [box, cyl(2, 4, 8), cyl(1, 8, 10)], area_x-pi+4*pi+pi)
    return out


# ------------------------------------------------------------------ Monte Carlo

def monte_carlo(chain, n, seed):
    """Volume and centre estimates of the three second operations and their
    standard errors, from `n` uniform points in the prisms' box."""
    fps = [FloatPrism(p) for p in chain.prisms]
    lo, hi = [], []
    for p in chain.prisms[:3]:
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

def first_result(obj, tool, op1, selector):
    """The first result `X` given to the second Boolean: its volume,
    moments and area, and its solid count; for a picked solid the chain
    `op1(A, B) and D` (one solid, the pair counting several)."""
    pair = cref.Pair(obj, tool)
    if selector is None:
        n, _, _, _ = pair.result(op1)
        vx, mx = pair.volumes()[op1]
        return vx, mx, pair.area(op1), n, n
    part = ref.Chain(obj, tool, selector, op1, False)
    vx, mx = part.volumes()['common']
    return vx, mx, part.area('common'), part.solids('common'), pair.result(op1)[0]


def evaluate(job):
    """One chain: the reference's rows for its second operations and every
    check's deviations (run in a worker process)."""
    name, obj, tool, third, op1, swapped, selector, pick, ops, mc_n = job
    chain = ref.Chain(obj, tool, third, op1, swapped, selector=selector)
    res = {op: chain.result(op) for op in ref.OPS}
    rows = {op: ref.rows(chain, op) for op in ops}
    size = chain.size
    vols = chain.volumes()
    checks = {}
    vx, mx, ax, solids, whole = first_result(obj, tool, op1, selector)
    assert solids == 1, f'{name}: the given first result is not one solid'
    if selector is not None:
        assert whole > 1, f'{name}: a pick among one solid'
        checks['separation'] = ref.separation(chain)/size**2
        # The pick point well inside the picked solid.
        fps = [FloatPrism(p) for p in (chain.prisms[0], chain.prisms[1], chain.prisms[3])]
        first = ref.SET[op1]
        for k in range(7):
            X = list(pick)
            if k:
                X[(k-1)//2] += 0.01 if k % 2 else -0.01
            a, b, d = (fp.contains(X) for fp in fps)
            assert first(a, b) and d, f'{name}: the pick point is not inside the picked solid'
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
    shared_with_third = False
    for i, f, cls, _ in chain.face_areas():
        others = [j for j in range(len(chain.prisms)) if j != i]
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
                out.append(f'{"ABCD"[i]} {cref.describe(f)} {c} {mp.nstr(v, 3)}')
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
    out = ['# case\trow (S9e.2, chained_boolean_reference.py: expect KIND S9e.2 CLASS, reason TEXT for a '
           'degenerate case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, _, _, _ = by_chain[case.chain_name]
        row = rows[case.then]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9e.2 {case.klass}')
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
    return {'boolean-given-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-given-expected.tsv': '\n'.join(out)+'\n',
            'boolean-given-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n):
    chains = {}
    for c in cases():
        chains.setdefault(c.chain_name, [c.obj, c.tool, c.third, c.op1, c.swapped, c.selector, c.pick,
                                         []])[7].append(c.then)
    return [(name, *rest[:7], rest[7], mc_n) for name, rest in chains.items()]


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
              'separation': 0, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-30,
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
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.then == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{sum(1 for c in listed if c.klass == k)} {k}' for k in CLASSES))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    for name, _, _, _, near, _ in results:
        if near:
            print('declared near coincidences:', name, '; '.join(near[:4]))


if __name__ == '__main__':
    main()
