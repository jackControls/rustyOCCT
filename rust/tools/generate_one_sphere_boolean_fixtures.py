#!/usr/bin/env python3
"""Fixtures for S9e.4b.3c.1 of REVIEW_NOTES.md: two pieces of one sphere
given to a Boolean (faces of both inputs on one sphere, in general
position), and an imported sphere piece whose rim OCCT split into two arcs
(the DRAW survey's `so1` and `so4`).

The bodies are OCCT's own output, as S9e.4b.3a's
(`generate_imported_pieces_boolean_fixtures.py`): `boolean-one-sphere-
bodies.txt` holds one `write NAME imported/NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a sphere's rows, a `boolean common` row and a
box's or prism's rows, and for the split rims a `divide` row:
`ShapeUpgrade_ShapeDivideClosedEdges` divides every closed edge of the
written solid in two), which `compare_one_sphere_boolean.py
--write-bodies` runs to write `rust/fixtures/imported/NAME.brep`. Nothing
here reads the kernel; the files are read only to check that each is the
body it claims to be.

`boolean-one-sphere-cases.txt` lists each case in the Boolean protocol (an
imported input its one `brep imported/NAME.brep` row, a chain's further
`then` rows as S9e.1's); `boolean-one-sphere-expected.tsv` gives per case
`expect KIND S9e.4b.3c.1 CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a later
sub-step's; the class `pieces`: two imported pieces of one sphere,
`sphere`: an imported piece and a constructed ball of its sphere, `chain`:
their result given to another Boolean), a degenerate or unsupported case's
`reason TEXT`, then `result N volume area cx cy cz` or `empty`;
`boolean-one-sphere-frames.tsv` the stored axes of every solid the kernel
builds from its rows.

Every body is a piece of one sphere, the ball of radius 5 about `(5, 5, 4)`
of S9e.4b.3a's `octant`, each on a turned rational frame of S9e.4b.3a's
(none of whose axes is a world axis or normal to one: the reference's
sphere families are the world's meridians and parallels): `hemi`, the ball
above its equator's plane on `SKEW4`, its rim divided (as `so1`'s); `cap`,
above the parallel's plane at `3/2` on `SKEW4`, its rim divided (as
`so4`'s); `octant2`, the corner of `SKEW2`'s axes (`so2`'s wedge: two
meridian planes and the equator's, its pole vertex inside the cap's sphere
face); `octant4`, the corner of `SKEW4`'s axes; `wedge`, above the
parallel's plane at `1` between two meridian planes of `SKEW` (`so5`'s, in
another frame). In `so4` and `so2` the wedge's pole vertex is the cap's
stored pole; in turned rational frames it would lie within rounding of it
(two planes' line against a rounded axis), so `octant2`'s axis is another.
Cases (each the three operations, the cut either way where it differs):
`hemi_cap`, `cap_hemi` (`so1` and `so4`: the cap inside the hemisphere,
their rims parallel); `cap_octant`, `octant_cap` (`so4` and `so2`: the
octant's meridian arcs across the cap's rim); `ball_cap`, a ball of the
sphere on the world's axes and the cap (a whole sphere's split great
circles only); `wedge_cap`, the wedge and the cap; `chain_cap`, the cap
less the octant, then with a level slab across it (its cut two solids).
Declared `degenerate`: `ball_near`, a ball about the sphere's centre whose
radius is `5 + 2^-30`, within the resolution of the cap's sphere and not it
(its fuse and common; its cut a shell thinner than the reference's rays
resolve). `hemi_octant`, the hemisphere and `octant4` on one frame (the
octant's equator arc on the hemisphere's rim circle, both bases on one
plane, as `so1` and `so2`), declared `unsupported` until S9e.4b.3c.2's
kernel decided its exact incidences, is solid (its margins those of any
solid case: the least sine 1.0, the events' spacing 4.9e-3).
S9e.4b.3a's `one_sphere` (its octant and `octant_low`) and
S9e.4b.3b's (a zone's half and its own ball) are this step's too, kept in
their sets.

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each imported piece its
first Boolean, two pieces `(P common B) op (Q common D)` with `P` and `Q`
one sphere (the reference's faces of two inputs on one surface, which
leave the area identity), with S9e.4b.3a's checks (the two families, each
solid's closed form, the pair identities, Monte Carlo, solid counts by
rays at two resolutions, here with neighbouring rays' intervals joined
within two grid spacings (`count_solids`: the slivers two pieces of one
sphere leave taper to it), each piece one solid, the meetings' sines and
the events' spacing outside the declared cases) and each piece's closed form
(`corner_closed`, the hemisphere, `cap_closed`). With `--check`, where the
bodies' files exist, each is read: its faces one sphere and planes, every
stored vertex on the construction's surfaces within 1e-12 of the size, and
a divided body's rim two arcs (its stored vertices on the rim's circle).
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case, native_case
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
import generate_imported_pieces_boolean_fixtures as pieces
from generate_curved_boolean_fixtures import square

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'sphere', 'chain')
BODIES_FILE = 'boolean-one-sphere-bodies.txt'
STEP = 'S9e.4b.3c.1'
# Not `SKEW3`: its normal `(10, 11, 2)` lies in the split plane of one of
# the kernel's whole sphere's own rotations (rows `(1, 2, 2)` and `(2, 1,
# -2)`), a piece's pole vertex on its own split great circle within
# rounding.
FRAMES = pieces.FRAMES
prism, sphere, construction = base.prism, base.sphere, base.construction
ref = base.ref
box = pieces.box


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]

# The one sphere: S9e.4b.3a's octant's, radius 5 about (5, 5, 4).
RADIUS = 5.0
CENTRE = (5, 5, 4)
ON4 = at('SKEW4', CENTRE)
ON2 = at('SKEW2', CENTRE)
ON1 = at('SKEW', CENTRE)
CAP_HEIGHT = 1.5
WEDGE_HEIGHT = 1.0


class Body(pieces.Body):
    """A piece `sphere common box` on the one sphere, its rim divided where
    `divide` (OCCT's `divide` row: every closed edge in two)."""

    def __init__(self, name, klass, primitive, box, closed=None, divide=False, path=None):
        super().__init__(name, klass, primitive, box, closed=closed)
        self.divide = divide
        if path:
            self.path = path

    def native_rows(self):
        return super().native_rows()+(['divide'] if self.divide else [])


def cap_closed(r, h):
    """A ball of radius `r` above a plane at height `h` from its centre:
    its volume `pi t^2 (3 r - t) / 3` (`t = r - h`), its area the zone `2
    pi r t` and the disc `pi (r^2 - h^2)`."""
    R, H = mp.mpf(r), mp.mpf(h)
    t = R-H
    return mp.pi*t**2*(3*R-t)/3, 2*mp.pi*R*t+mp.pi*(R**2-H**2)


def stored_length(frame):
    """The length of a frame's stored normal (binary64 components, a unit
    within rounding): a parallel's plane at height `h` of the frame lies
    `h` times it from the centre."""
    n = stored_axes(frame)[3]
    q = lambda c: mp.mpf(F(c).numerator)/F(c).denominator
    return mp.sqrt(sum(q(c)**2 for c in n))


def make_bodies():
    out = [
        Body('sphere_hemi', 'sphere', sphere(RADIUS, ON4), prism([square(-6.0, -6.0, 6.0, 6.0)], ON4, 0.0, 6.0),
             closed=lambda: cap_closed(5, 0), divide=True),
        # Its box's walls and top apart from the hemisphere's (the
        # reference's rays meet no faces of two inputs on one plane).
        Body('sphere_cap', 'sphere', sphere(RADIUS, ON4),
             prism([square(-5.5, -5.5, 5.5, 5.5)], ON4, CAP_HEIGHT, 5.75),
             closed=lambda: cap_closed(5, CAP_HEIGHT*stored_length(ON4)), divide=True),
        Body('sphere_octant2', 'sphere', sphere(RADIUS, ON2), box(ON2, (8.0, 8.0, 8.0)),
             closed=lambda: pieces.corner_closed(5, ON2)),
        Body('sphere_octant4', 'sphere', sphere(RADIUS, ON4), box(ON4, (8.0, 8.0, 8.0)),
             closed=lambda: pieces.corner_closed(5, ON4)),
        Body('sphere_wedge', 'sphere', sphere(RADIUS, ON1), prism([square(0.0, 0.0, 8.0, 8.0)], ON1, WEDGE_HEIGHT,
                                                                   8.0)),
    ]
    bodies = {b.name: b for b in out}
    return bodies


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported piece its
    primitive and box) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.specs()
            exprs.append((b.op, k, k+1))
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(pieces.Imported):
    """One case, its imported inputs this step's bodies."""

    def __init__(self, name, klass, items, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.items = items
        self.solid_cases = []
        for k, s in enumerate(items):
            if s[0] == 'imported':
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=BODIES[s[1]].path)
            else:
                c = construction(s, OPERATIONS[k])
                c.name = name
            self.solid_cases.append(c)

    def expr(self):
        _, e = flatten(self.items)
        if not self.stages:
            return (self.op1, e[0], e[1])
        (op2, sw2) = self.stages[0]
        first = (self.op1, e[0], e[1])
        return (op2, e[2], first) if sw2 else (op2, first, e[2])


def group(name, klass, items, outcomes, first=None, reason=None, kind=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty)."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('degenerate' if reason else ('solid' if n else 'empty'))
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, reason))
    return out


NEAR = ('a ball about the sphere\'s centre whose radius is 5 + 2^-30: within the resolution of the cap\'s sphere and '
        'not it')
SLAB = prism([square(-4.0, -4.0, 14.0, 14.0)], at('XY', (0, 0, 6.375)), 0.0, 1.5)
BALL = sphere(RADIUS, at('XY', CENTRE))
NEAR_BALL = sphere(RADIUS+2.0**-30, at('XY', CENTRE))


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('hemi_cap', 'pieces', [imported('sphere_hemi'), imported('sphere_cap')], three)
    out += group('cap_hemi', 'pieces', [imported('sphere_cap'), imported('sphere_hemi')], {'cut': 0})
    out += group('cap_octant', 'pieces', [imported('sphere_cap'), imported('sphere_octant2')], three)
    out += group('octant_cap', 'pieces', [imported('sphere_octant2'), imported('sphere_cap')], {'cut': 1})
    out += group('ball_cap', 'sphere', [BALL, imported('sphere_cap')], three)
    out += group('wedge_cap', 'pieces', [imported('sphere_wedge'), imported('sphere_cap')], three)
    out += group('chain_cap', 'chain', [imported('sphere_cap'), imported('sphere_octant2'), SLAB],
                 {'fuse': 1, 'cut': 2, 'common': 1}, first='cut')
    out += group('ball_near', 'sphere', [NEAR_BALL, imported('sphere_cap')], {'fuse': 1, 'common': 1}, reason=NEAR)
    out += group('hemi_octant', 'pieces', [imported('sphere_hemi'), imported('sphere_octant4')], three)
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.items) == 2+len(c.stages), c.name
        assert c.last in OPS, c.name
        assert any(s[0] == 'imported' for s in c.items), c.name
        if c.klass == 'pieces':
            assert all(s[0] == 'imported' for s in c.items), c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        # The first Boolean's inputs on one sphere.
        specs, _ = flatten(c.items[:2])
        balls = {(s[1], tuple(s[2][:3])) for s in specs if s[0] == 'sphere'}
        assert len({(round(r, 6), o) for r, o in balls}) == 1, (c.name, balls)


# ------------------------------------------------------------------ the reference

def count_solids(chain, expr, n):
    """The chained reference's `count_solids` with each ray's intervals
    joined to a neighbouring ray's that overlap or lie within two grid
    spacings of it: two pieces of one sphere leave slivers tapering to its
    surface (a wedge less a cap), whose neighbouring rays' intervals are
    disjoint where the sliver is thinner than its slope across a spacing,
    which the reference's strict overlap counts as several solids."""
    used = sorted(ref.leaves(expr))
    fms = {j: ref.Float(chain.inputs[j]) for j in used}
    boxes = [chain.inputs[j].box() for j in used]
    lo = [min(b[0][i] for b in boxes) for i in range(3)]
    hi = [max(b[1][i] for b in boxes) for i in range(3)]
    pad = 1e-3*max(hi[i]-lo[i] for i in range(3))
    lo = [x-pad for x in lo]
    hi = [x+pad for x in hi]
    tmax = hi[2]-lo[2]
    xs = [lo[0]+(hi[0]-lo[0])*(i+0.3819660112501051)/n for i in range(n)]
    ys = [lo[1]+(hi[1]-lo[1])*(j+0.6180339887498949)/n for j in range(n)]
    reach = 2*max(hi[0]-lo[0], hi[1]-lo[1])/n
    rays = {}
    mem = [False]*len(chain.inputs)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            P, D = [x, y, lo[2]], [0.0, 0.0, 1.0]
            ivs = {k: fms[k].intervals(P, D, tmax) for k in used}
            pts = sorted(set([0.0, tmax]+[t for k in used for iv in ivs[k] for t in iv]))
            out = []
            for a, b in zip(pts, pts[1:]):
                m = (a+b)/2
                for k in used:
                    mem[k] = any(p < m < q for p, q in ivs[k])
                if ref.evaluate(expr, mem):
                    if out and out[-1][1] == a:
                        out[-1] = (out[-1][0], b)
                    else:
                        out.append((a, b))
            rays[(i, j)] = out
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for (i, j), ivs in rays.items():
        for k in range(len(ivs)):
            parent[(i, j, k)] = (i, j, k)
    for (i, j), ivs in rays.items():
        for di, dj in ((1, 0), (0, 1)):
            nb = rays.get((i+di, j+dj))
            if not nb:
                continue
            for k, (a, b) in enumerate(ivs):
                for kk, (a2, b2) in enumerate(nb):
                    if min(b, b2)-max(a, a2) > -reach:
                        ra, rb = find((i, j, k)), find((i+di, j+dj, kk))
                        if ra != rb:
                            parent[ra] = rb
    return len({find(a) for a in parent})


def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3a's, its bodies this step's)."""
    name, items, op1, swapped, chained, ops, degenerate, mc_n = job
    specs, sub = flatten(items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    size = chain.size
    first = (op1, sub[0], sub[1]) if chained else None

    def expr_of(op):
        if not chained:
            return (op, sub[0], sub[1])
        return (op, sub[2], first) if swapped else (op, first, sub[2])
    exprs = {op: (expr_of(op), solids) for op, solids in ops.items()}
    rows, res = {}, {}
    for op, (expr, solids) in exprs.items():
        rows[op] = ref.rows(chain, expr, solids)
        res[op] = chain.result(expr, solids)
    checks = {}
    fam = mp.mpf(0)
    for op, (expr, _) in exprs.items():
        a, b = chain.measures(expr, 0), chain.measures(expr, 1)
        fam = max(fam, abs(a[0]-b[0])/size**3, abs(a[2]-b[2])/size**2,
                  max(abs(a[1][i]-b[1][i]) for i in range(3))/size**4)
    checks['families'] = fam
    cdev = mp.mpf(0)
    for k, I in enumerate(chain.inputs):
        V, mom, A = I.closed()
        v, m, a = chain.measures(k)
        cdev = max(cdev, abs(v-V)/size**3, abs(a-A)/size**2, max(abs(m[i]-mom[i]) for i in range(3))/size**4)
    checks['solid_closed_forms'] = cdev
    for it, e in zip(items, sub):
        if it[0] != 'imported':
            continue
        assert count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
        closed = BODIES[it[1]].closed
        if closed is not None:
            V, A = closed()
            v, _, a = chain.measures(e)
            d = max(abs(v-V)/size**3, abs(a-A)/size**2)
            checks['piece_closed_forms'] = max(checks.get('piece_closed_forms', mp.mpf(0)), d)
    if chained:
        assert count_solids(chain, first, 48) == 1, f'{name}: the given result is not one solid'
        x_expr, c_expr = (sub[2], first) if swapped else (first, sub[2])
    else:
        x_expr, c_expr = sub[0], sub[1]
    vx, mx, ax = chain.measures(x_expr)
    vc, mc, ac = chain.measures(c_expr)
    m = {op: chain.measures(e) for op, (e, _) in exprs.items()}
    if set(m) == set(OPS):
        vf, mf, af = m['fuse']
        vn, mn, an = m['common']
        vt, mt, at_ = m['cut']
        dev = abs(vf+vn-vx-vc)/size**3
        for i in range(3):
            dev = max(dev, abs(mf[i]+mn[i]-mx[i]-mc[i])/size**4)
        dev = max(dev, abs(vt-(vx-vn))/size**3, max(abs(mt[i]-(mx[i]-mn[i])) for i in range(3))/size**4)
        checks['pair_identities'] = dev
        owner = pieces.owners(items, sub)
        shared = any(owner[i] != owner[sw.others.surfs[k][0]] for i, _, sw in chain.sweeps(0) for k in sw.skip)
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    for op, (expr, solids) in exprs.items():
        for n in (48, 71):
            got = count_solids(chain, expr, n)
            assert got == solids, f'{name} {op}: {got} solids by rays ({n}), {solids} declared'
    mc = ref.monte_carlo(chain, [e for e, _ in exprs.values()], mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op, (expr, _) in exprs.items():
        vol, c, sv, sc = mc[expr]
        v, mom, _ = chain.measures(expr)
        z = max(z, abs(vol-float(v))/max(sv, 1e-12))
        if c is not None and v > 0:
            cent = [float(mom[i]/v) for i in range(3)]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    sine, gap, _ = chain.margins()
    margins = {'sine': sine, 'gap': gap, 'spacing': base.spacing(chain)}
    stats = {'quadrature': chain.quad_error()/size**4, 'missed': chain.missed(), 'refined': chain.refined()}
    return name, rows, res, checks, margins, stats, degenerate


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, c.op1 if c.stages else None, bool(c.stages and c.stages[0][1]),
                                        bool(c.stages), {}, c.kind in ('degenerate', 'unsupported')])
        e[4][c.last] = c.solids
    return [(name, items, op1, swapped, chained, ops, degenerate, mc_n)
            for name, (items, op1, swapped, chained, ops, degenerate) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('ONE_SPHERE_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+repr(flatten(job[1])[0])+inspect.getsource(evaluate)
                          + Path(ref.__file__).read_text()).encode()).hexdigest()[:24]
    path = Path(cache)/f'{job[0]}-{key}.pickle'
    if path.exists():
        return pickle.loads(path.read_bytes())
    try:
        out = evaluate(job)
    except Exception as e:
        raise RuntimeError(f'{job[0]}: {e!r}') from e
    path.write_bytes(pickle.dumps(out))
    return out


LIMITS = pieces.LIMITS
reference_checks = pieces.reference_checks


# ------------------------------------------------------------------ the files

def bodies_text():
    """The writer's blocks: each new body's rows under `write NAME PATH`
    (S9e.4b.3a's reused bodies are its own)."""
    out = []
    for b in BODIES.values():
        if getattr(b, 'reused', False):
            continue
        out += [f'write {b.name} {b.path}', *b.native_rows(), 'end']
    return '\n'.join(out)+'\n'


def generate(results):
    by_group = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}; chained_curved_boolean_reference.py: expect KIND {STEP} CLASS, reason TEXT for a '
           'degenerate or unsupported case, then result N volume area cx cy cz or empty)']
    for case in all_cases():
        blocks.append(case.encode())
        row = by_group[case.group][1][case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True, 'unsupported': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP} {case.klass}')
        if case.kind in ('degenerate', 'unsupported'):
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (solid k) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits; the solids '
              'the kernel builds from rows, not the imported ones)']
    for case in all_cases():
        for k, (item, c) in enumerate(zip(case.items, case.solid_cases)):
            if item[0] == 'imported':
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-one-sphere-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-one-sphere-expected.tsv': '\n'.join(out)+'\n',
            'boolean-one-sphere-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_one_sphere_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its faces' kinds one sphere and planes, each
    stored vertex on the construction's surfaces within 1e-12 of the body's
    size; a divided body's rim two arcs: two stored vertices on its rim's
    circle (its sphere and its lowest plane), where the construction's ring
    has one. Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        curved = [k for k in kinds if k != 'plane']
        assert curved == ['sphere'], (b.name, kinds)
        surfs = []
        for spec in b.specs():
            surfs += base.construction_surfaces(spec)[1]
        size = max([1.0]+[abs(x) for v in vertices for x in v])
        dist = lambda s, X: abs(s.f(X))/max(base.tc.norm(s.grad(X)), mp.mpf(10)**-30)
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            worst = max(worst, min(dist(s, X) for s in surfs)/size)
        if b.divide:
            # The rim: the sphere's points on the box's lowest plane.
            low = b.box[3]
            frame = pieces.spec_frame(b.box)
            o, _, _, n = stored_axes(frame)
            q = lambda c: mp.mpf(F(c).numerator)/F(c).denominator
            nn = mp.sqrt(sum(q(c)**2 for c in n))
            height = lambda X: sum((X[i]-q(o[i]))*q(n[i]) for i in range(3))/nn-mp.mpf(low)
            ball = surfs[0]
            on_rim = [v for v in vertices if abs(height(tuple(mp.mpf(c) for c in v))) <= 1e-9*size
                      and dist(ball, tuple(mp.mpf(c) for c in v)) <= 1e-9*size]
            assert len(on_rim) == 2, (b.name, 'a divided rim of two arcs', len(on_rim))
    assert worst <= 1e-12, ('a stored vertex off its construction', mp.nstr(worst, 3))
    return worst, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=100000, help='Monte-Carlo points per group')
    parser.add_argument('--only', help='one group, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each group\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['ONE_SPHERE_CACHE'] = args.cache
    listed = all_cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in base.run(todo, cached, 1):
            print(r)
        return
    results = base.run(todo, cached, args.workers)
    worst, covered = reference_checks(results)
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
    deviation, read = check_bodies()
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.last == op)} {op}' for op in OPS),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{sum(1 for c in listed if c.klass == k)} {k}' for k in CLASSES))
    print('reference checks (largest deviation, relative to the case size; least_sine and event_spacing as '
          'reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print(f'bodies: {read} of {len(BODIES)} files read, every stored vertex within {mp.nstr(deviation, 3)} of the '
          'size of its construction\'s surfaces')


if __name__ == '__main__':
    main()
