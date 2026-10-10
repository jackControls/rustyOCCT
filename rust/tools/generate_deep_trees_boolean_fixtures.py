#!/usr/bin/env python3
"""Fixtures for S9e.4b.4c.2a of REVIEW_NOTES.md: imported bodies of one
cylinder face and plane faces that are a Boolean tree of their
primitive and convex hulls of their planes whose pockets hold pockets of
their own (a tooth or a post standing in a pocket, a U island whose notch is
a pocket again; bodies OCCT wrote to `.brep` files), given to Booleans.

The bodies are OCCT's own output: S9e.4b.3c.3b's `tooth`
(`generate_piece_trees_boolean_fixtures.py`, its file
`rust/fixtures/imported/form_tooth.brep` read again, not written again),
S9e.4a's `cyl` (`generate_imported_boolean_fixtures.py`) and this step's
own, `boolean-deep-trees-bodies.txt` holding one `write NAME
imported/deep_NAME.brep` block per body for `occt_boolean_oracle.cpp` (the
primitive's rows, a `boolean cut` row and a prism of a square with a hole:
the Boolean's one solid written), which `compare_deep_trees_boolean.py
--write-bodies` runs to write `rust/fixtures/imported/deep_NAME.brep` (format
version 1, no triangulations). Nothing here reads the kernel; the files are
read only to check that each is the body it claims to be (below).

`boolean-deep-trees-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/NAME.brep` row, a chain's further
`then` rows as S9e.1's; `boolean-deep-trees-expected.tsv` gives per case
`expect KIND S9e.4b.4c.2a CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a later step's; the
class `pieces`: an imported body and a construction, `both`: two imported
inputs, `chain`: an imported input's result given to another Boolean), a
degenerate or unsupported case's `reason TEXT`, then `result N volume area
cx cy cz` (totals over the N solids, world coordinates) or `empty`;
`boolean-deep-trees-frames.tsv` the stored axes of every solid the kernel
builds from its rows.

The bodies, each a primitive less a prism whose pockets nest (every section
a circle or a line: the prisms' walls meet the primitive's top cap alone, in
lines, or the tooth's cylinder along lines): `tooth` (S9e.4b.3c.3b's: a
cylinder of radius 3 about `(4, 4)` over `[0, 6]` less a U prism over `[2,
9/2]` entering its wall, its slot inside the cylinder: a tooth standing in
the pocket, two deep); `post`, a cylinder of radius 4 on `R815` over `[0,
6]` less a prism of the square `[-5/2, 5/2]^2` with the hole `[-1, 1]^2`
over `[3, 7]` (a square post standing in a square pocket from the cap, its
top on the cap's plane: two deep); `well`, a cylinder of radius 4 on
`R125` over `[0, 13/2]` less a prism of the square `[-5/2, 5/2]^2` with a
U-shaped hole (the square `[-3/2, 3/2]^2` less its notch `[-1/2, 3/2] x
[-1/2, 1/2]`, open toward `+x`) over `[5/2, 15/2]` (a U island standing in a
square pocket, its notch a pocket of the island's: three deep).

Cases (each the three operations): `tooth_pin`, the tooth and an upright rod
through the tooth's wall into the pocket; `post_rod`, the post and an
upright rod on its frame through the post's wall; `well_rod`, the well and
an upright rod on its frame through the island's arm into its notch;
`rod_post`, an upright rod on the post's frame less the post, through the
post's wall (the body the tool); `post_cyl`, the post and S9e.4a's imported
cylinder across the post, the pocket and the cylinder's wall (both
imported); `chain_post`, the post less an upright rod through the post,
then with an upright rod through the pocket's wall. Declared `degenerate`:
`post_touch`, the post and a ball resting on the post's top at its centre,
exactly (`Degenerate("a tangency between the inputs (S9c)")`). The
reference's cost grows with a case's faces and with oblique and curved
partners: a first set with the well on a frustum, a `TILT` rod less the
well, S9e.4a's ball across the post and the chain's `TILT` slab ran over 37
minutes of a worker for each of its three slowest groups (the four others
6 to 15), so every partner but the declared ball stands along the bodies'
axes and the well's primitive is a cylinder.

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`, which models convex
profiles): a prism of a square with a square hole the square's box less the
hole's box (past both caps), with a U-shaped hole the square's box less the
U's box less the notch's box (past the U's open side and both caps), the
tooth's U S9e.4b.3c.3b's box less its slot's box; each imported body its
Boolean, `(P cut X) op C` (swapped, the body the tool), two bodies, a chain
`((P cut X) op1 C) op2 D`, with S9e.4b.3c.3a's checks (the two families,
each solid's closed form, the pair identities on the last Boolean's
arguments, Monte Carlo, solid counts by rays at two resolutions, each body
one solid, the meetings' sines and the events' spacing outside the declared
group) and each body's closed form (times the stored axes' determinant).
With `--check`, where the bodies' files exist, each is read
(`stored_records`): its curved face one cylinder, every stored
vertex within 1e-12 of the size on the construction's surfaces.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, native_case
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_piece_trees_boolean_fixtures as trees
import generate_primitive_chains_boolean_fixtures as chains
from generate_curved_boolean_fixtures import square

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-deep-trees'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4c.2a'
# `R815`: the world's `z` as the normal, its `x` turned to `(8, 15, 0) / 17`
# (rounded): a cylinder's walls along its axis, their normals in the world's
# `xy` exactly.
FRAMES = dict(forms.FRAMES, R815=(0.0, 0.0, 1.0, 8.0, 15.0, 0.0))
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
q, point, moved = forms.q, forms.point, forms.moved
cylinder, block = forms.cylinder, forms.block
det = chains.det


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

def u_hole(x0, y0, x1, y1, notch):
    """A U: the square `[x0, x1] x [y0, y1]` less its notch `(s0, t0, t1)`,
    `[s0, x1] x [t0, t1]` open at `x = x1`, counter-clockwise."""
    s0, t0, t1 = notch
    return Boundary(points=[(x0, y0), (x1, y0), (x1, t0), (s0, t0), (s0, t1), (x1, t1), (x1, y1), (x0, y1)])


class Holed:
    """A prism of a square with a hole for OCCT (`native`: one prism of both
    boundaries) and the reference's same set (`parts`, `expr` over them): the
    square's box less the hole's, a square hole its box past both caps, a
    U-shaped hole (`notch`) its box less the notch's box past the U's open
    side and both caps."""

    def __init__(self, frame, h0, h1, rect, hole, notch=None):
        x0, y0, x1, y1 = rect
        a0, b0, a1, b1 = hole
        inner = square(a0, b0, a1, b1) if notch is None else u_hole(a0, b0, a1, b1, notch)
        self.native = prism([square(x0, y0, x1, y1), inner], frame, h0, h1)
        self.parts = [block(frame, x0, y0, x1, y1, h0, h1), block(frame, a0, b0, a1, b1, h0-1.0, h1+1.0)]
        self.expr = ('cut', 0, 1)
        if notch is not None:
            s0, t0, t1 = notch
            self.parts.append(block(frame, s0, t0, (a1+x1)/2, t1, h0-2.0, h1+2.0))
            self.expr = ('cut', 0, ('cut', 1, 2))


def shifted(expr, k):
    """An expression over a body's own solids, its leaves moved by `k`."""
    if isinstance(expr, int):
        return expr+k
    return (expr[0], shifted(expr[1], k), shifted(expr[2], k))


class Body:
    """An imported body `primitive cut holed`, OCCT's rows for it, its
    closed form `(volume, area)` (the area `None` where not checked) and how
    deep its pockets nest (`depth`)."""

    def __init__(self, name, primitive, holed, depth, closed):
        self.name, self.first, self.holed, self.depth, self.closed = name, primitive, holed, depth, closed
        self.kind, self.reason = None, None
        self.path = f'imported/deep_{name}.brep'
        self.ref_specs = [primitive]+holed.parts
        self.ref_expr = ('cut', 0, shifted(holed.expr, 1))

    def specs(self):
        return [self.first, self.holed.native]

    def primitive(self):
        return self.first

    def native_rows(self):
        rows = lambda spec: native_case(construction(spec, 0)).split('\n')[1:-1]
        return rows(self.first)+['boolean cut']+rows(self.holed.native)


# The post on `R815`, placed so that S9e.4a's cylinder (radius 3 about the
# world's `(5, 5)`, over `[-1, 5]`) has its axis at the chart's `(1.71,
# 0.77)`: across the post, the pocket's walls and the cylinder's wall,
# missing every other plane along it by at least 0.27, its top above the
# cap's plane by 1/4.
POST_FRAME = at('R815', (4.875, 3.125, -1.25))
# The well on `R125` (its `x` rounded), its cylinder about its frame's origin.
WELL_FRAME = at('R125', (1, -1, 0))
TOOTH_FRAME = trees.BITES_FRAME


def post_volume():
    """The cylinder `96 pi` less the ring `25 - 4` over the pocket's depth
    3, in the chart, times the stored axes' determinant."""
    return (96*mp.pi-63)*det(POST_FRAME)


def well_volume():
    """The cylinder `104 pi` less the square `25` less the U island `9 - 2`
    over the pocket's depth 4, in the chart, times the stored axes'
    determinant."""
    return (104*mp.pi-72)*det(WELL_FRAME)


def make_bodies():
    out = [
        Body('post', cylinder(POST_FRAME, 0.0, 0.0, 4.0, 0.0, 6.0),
             Holed(POST_FRAME, 3.0, 7.0, (-2.5, -2.5, 2.5, 2.5), (-1.0, -1.0, 1.0, 1.0)), 2,
             closed=lambda: (post_volume(), None)),
        Body('well', cylinder(WELL_FRAME, 0.0, 0.0, 4.0, 0.0, 6.5),
             Holed(WELL_FRAME, 2.5, 7.5, (-2.5, -2.5, 2.5, 2.5), (-1.5, -1.5, 1.5, 1.5), notch=(-0.5, -0.5, 0.5)), 3,
             closed=lambda: (well_volume(), None)),
    ]
    return {b.name: b for b in out}


# This step's own bodies, and S9e.4b.3c.3b's tooth read again.
OWN = make_bodies()
BODIES = dict(OWN, tooth=trees.BODIES['tooth'])


def imported(name):
    return ('imported', name)


def earlier(name):
    """S9e.4a's body `name` (its file read again)."""
    return ('base', name)


def path_of(item):
    return (BODIES[item[1]] if item[0] == 'imported' else base.BODIES[item[1]]).path


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its own
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(shifted(b.ref_expr, k))
        elif it[0] == 'base':
            exprs.append(len(specs))
            specs.append(base.BODIES[it[1]].spec)
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(pieces.Imported):
    """One case, its imported inputs this step's bodies, the tooth or S9e.4a's
    cylinder."""

    def __init__(self, name, klass, items, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.items = items
        self.solid_cases = []
        for k, s in enumerate(items):
            if s[0] in ('imported', 'base'):
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=path_of(s))
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


def group(name, klass, items, outcomes, first=None, kind=None, reason=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); the group's own refusal (`kind`, `reason`)
    declares its cases."""
    out = []
    for op, n in outcomes.items():
        if kind:
            k, why = kind, reason
        else:
            k, why = ('solid' if n else 'empty'), None
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, why))
    return out


TOOTH_PIN = cylinder(TOOTH_FRAME, 2.125, 3.75, 0.3125, -1.0, 7.0)
POST_ROD = cylinder(POST_FRAME, 1.25, 0.25, 0.5, -1.0, 7.0)
WELL_ROD = cylinder(WELL_FRAME, 0.5, 0.625, 0.375, -1.0, 7.5)
ROD_POST = cylinder(POST_FRAME, -1.125, 0.25, 0.375, -1.0, 7.0)
CHAIN_ROD = cylinder(POST_FRAME, 0.25, -0.25, 0.375, -1.0, 7.0)
CHAIN_PIN = cylinder(POST_FRAME, 2.25, -1.5, 0.625, 1.0, 8.0)
POST_TOUCH = sphere(1.0, at('XY', point(POST_FRAME, 0.0, 0.0, 7.0)))
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('well_rod', 'pieces', [imported('well'), WELL_ROD], three)
    out += group('post_cyl', 'both', [imported('post'), earlier('cyl')], three)
    out += group('chain_post', 'chain', [imported('post'), CHAIN_ROD, CHAIN_PIN], three, first='cut')
    out += group('rod_post', 'pieces', [ROD_POST, imported('post')], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('post_rod', 'pieces', [imported('post'), POST_ROD], three)
    out += group('tooth_pin', 'pieces', [imported('tooth'), TOOTH_PIN], three)
    out += group('post_touch', 'pieces', [imported('post'), POST_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
                 kind='degenerate', reason=TOUCH)
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
        if c.klass == 'both':
            assert sum(s[0] in ('imported', 'base') for s in c.items) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone'), c.name


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3c.3a's, each body's closed form where it has one)."""
    name, items, op1, swapped, chained, ops, declared, mc_n = job
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
    # Each imported body one solid, its closed form where it has one.
    for it, e in zip(items, sub):
        if it[0] != 'imported':
            continue
        assert forms.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
        closed = BODIES[it[1]].closed
        if closed is not None:
            V, A = closed()
            v, _, a = chain.measures(e)
            d = abs(v-V)/size**3
            if A is not None:
                d = max(d, abs(a-A)/size**2)
            checks['body_closed_forms'] = max(checks.get('body_closed_forms', mp.mpf(0)), d)
    # The pair identities on the last Boolean's arguments.
    if chained:
        assert forms.count_solids(chain, first, 48) == 1, f'{name}: the given result is not one solid'
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
        # Faces of two inputs of the case on one surface (not a body's own
        # solids') leave the area identity.
        owner = pieces.owners(items, sub)
        shared = any(owner[i] != owner[sw.others.surfs[k][0]] for i, _, sw in chain.sweeps(0) for k in sw.skip)
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    for op, (expr, _) in exprs.items():
        for n in (48, 71):
            got = forms.count_solids(chain, expr, n)
            assert got == ops[op], f'{name} {op}: {got} solids by rays ({n}), {ops[op]} declared'
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
    return name, rows, res, checks, margins, stats, declared


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, c.op1 if c.stages else None, bool(c.stages and c.stages[0][1]),
                                        bool(c.stages), {}, c.kind in ('degenerate', 'unsupported')])
        e[4][c.last] = c.solids
    return [(name, items, op1, swapped, chained, ops, declared, mc_n)
            for name, (items, op1, swapped, chained, ops, declared) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('DEEP_TREES_CACHE')
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


LIMITS = forms.LIMITS
# Groups whose events' spacing is not checked (their sines and gaps are).
SPACING_UNCHECKED = set()


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # A declared group's reference meets a tangency: its checks kept
        # apart, looser.
        for key, value in checks.items():
            note(('degenerate_'+key) if declared and key != 'monte_carlo_sigma' else key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not declared:
            note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
            if name not in SPACING_UNCHECKED:
                note('event_spacing', 1/max(margins['spacing'], mp.mpf(10)**-40))
    return worst, covered


# ------------------------------------------------------------------ the files

def bodies_text():
    """The writer's blocks: each of this step's bodies' rows under `write
    NAME PATH`."""
    out = []
    for b in OWN.values():
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
            if item[0] in ('imported', 'base'):
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q_).hex() for q_ in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_deep_trees_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file (this step's and the tooth's): its curved
    faces one cylinder, each stored vertex on
    the construction's surfaces within 1e-12 of the body's size. Returns the
    largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        curved = [k for k in kinds if k != 'plane']
        p = b.primitive()
        want = 'cylinder' if p[0] == 'prism' else p[0]
        assert curved == [want], (b.name, kinds)
        surfs = []
        for spec in b.ref_specs:
            surfs += base.construction_surfaces(spec)[1]
        size = max([1.0]+[abs(x) for v in vertices for x in v])
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            d = min(abs(s.f(X))/max(base.tc.norm(s.grad(X)), mp.mpf(10)**-30) for s in surfs)
            worst = max(worst, d/size)
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
        os.environ['DEEP_TREES_CACHE'] = args.cache
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
