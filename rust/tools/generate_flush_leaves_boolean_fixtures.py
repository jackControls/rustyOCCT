#!/usr/bin/env python3
"""Fixtures for S9e.4b.4b.2b.2b.2 of REVIEW_NOTES.md: imported bodies of
several cylinder faces and plane faces holding a part flush with a prism
leaf's walls (bodies OCCT wrote to `.brep` files), given to Booleans: a
prism within the leaf's heights fused with it and standing out of its walls
(its loops partly the leaf's walls), or cut from it across its walls (a
pocket open sideways, closed above by the part's other cap or open through
the leaf's cap).

The bodies are OCCT's own output, as S9e.4b.4b.2b.2b.1's
(`generate_rim_leaves_boolean_fixtures.py`): `boolean-flush-leaves-bodies.txt`
holds one `write NAME imported/flush_NAME.brep` block per body for
`occt_boolean_oracle.cpp` (the first prism's rows, a `boolean` row and the
second's: the Boolean's one solid written), which
`compare_flush_leaves_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/flush_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-flush-leaves-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/flush_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-flush-leaves-expected.tsv` gives
per case `expect KIND S9e.4b.4b.2b.2b.2 CLASS` (`solid`, `empty`,
`degenerate`: a `Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a
later step's; the class `pieces`: an imported body and a construction,
`both`: two imported inputs, `chain`: an imported input's result given to
another Boolean), a degenerate or unsupported case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-flush-leaves-frames.tsv` the stored axes
of every solid the kernel builds from its rows.

Each body is one Boolean OCCT was given of two prisms on one frame, every
section a circle or a line, the frames turned about the world's `z` (their
normals the world's `z` exactly), the leaf a stadium plate of half discs of
radius 2 about `(0, 0)` and `(6, 0)`: `arms`, the plate over `[0, 3]` on
`R724` fused with the box `[5, 15/2] x [-3, 3]` over `[1, 2]`, within its
heights: the box stands out of the plate's line walls and arc end in two
parts, each part's loop running along the plate's line wall, through its
joint and along its arc (the box's walls cross the lines at `x = 5` and the
arc at `x = 15/2`, `y = +-sqrt 7 / 2`); `slot`, the plate over `[0, 3]` on
`R815` less the square `[-3, 1] x [-3, 3]` with a hole of radius 1 about
`(-1/2, 1/4)` over `[1, 2]`, within its heights: a pocket open sideways
across the plate's arc end and both line walls, closed above and below,
the hole's post standing in it from its floor to its ceiling; `open`, the
plate over `[0, 5/2]` on `R125` less the square `[-4, 2] x [-3, 3]` with
the hole `[-5/2, 1] x [-1, 1]` over `[1, 4]`, from inside the plate past
its top: a pocket open sideways and through the plate's top, its post the
hole's square within the plate, flush with the plate's arc (its walls `y =
+-1` meet the arc at `x = -sqrt 3`).

Cases (each the three operations): `arms_rod`, the arms and an upright rod
across a part's wall `x = 5`, its top and the plate's line wall;
`slot_rod`, the slot and an upright rod across the pocket's wall `x = 1`;
`rod_slot`, an upright rod less the slot (the body the tool) across the
post and the plate's arc; `open_rod`, the open pocket and an upright rod
across its wall `x = 2`; `arms_slot`, the arms and the slot both imported
(their frames turned apart about `z`); `chain_arms`, the arms less an
upright rod through the plate, then with a box on its frame across a part
and the plate's top. Declared `degenerate`: `slot_touch`, the slot and a
ball of radius 1 resting on the plate's top cap at a point, exactly
(`Degenerate("a tangency between the inputs (S9c)")`). The bodies'
partners stand along their axes, on their frames (normals the world's `z`
exactly, so the import's leaves and the rods are parallel exactly).

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`, which models convex
profiles): each imported body its Boolean of two prisms, the square with a
round hole its box less its hole's cylinder (`Holed`), the square with a
square hole its box less its hole's box (`Ring`), `(X op Y) op2 C`
(swapped, the body the tool), two bodies `(X op Y) op2 (Z op W)`, a chain
`((X op Y) op2 C) op3 D`, with S9e.4b.3c.3a's checks (the two families,
each solid's closed form, the pair identities on the last Boolean's
arguments, Monte Carlo, solid counts by rays at two resolutions, each body
one solid, the meetings' sines and the events' spacing outside the declared
group) and each body's closed form (its prisms' sections in their chart
times their stored axes' determinant). With `--check`, where the bodies'
files exist, each is read (`stored_records`): its curved faces' kinds its
prisms' arcs', every stored vertex within 1e-12 of the size on the
constructions' surfaces.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, native_case
from generate_curved_boolean_fixtures import disc
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_piece_trees_boolean_fixtures as trees
import generate_primitive_chains_boolean_fixtures as chains
import generate_plane_parts_boolean_fixtures as parts
import generate_two_leaves_boolean_fixtures as two
import generate_rim_leaves_boolean_fixtures as rim

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-flush-leaves'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4b.2b.2b.2'
FRAMES = dict(two.FRAMES, R125=(0.0, 0.0, 1.0, 12.0, 5.0, 0.0))
prism, sphere, construction = base.prism, base.sphere, base.construction
ref = base.ref
q, point = forms.q, forms.point
cylinder, block = forms.cylinder, forms.block
det, Holed, Ring = chains.det, parts.Holed, rim.Ring
stadium_at, stadium_area = rim.stadium_at, rim.stadium_area


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Body(rim.Body):
    """An imported body `first op second` of two prisms (a `Holed` or `Ring`
    one the reference's set of its own), OCCT's rows for it, its closed form
    `(volume, area)` (the area `None` where not checked)."""

    def __init__(self, name, first, op, second, closed=None, kind=None, reason=None):
        super().__init__(name, first, op, second, closed, kind, reason)
        self.path = f'imported/flush_{name}.brep'


# The bodies' frames (each body's two prisms on one), the slot across the
# arms (`arms_slot`), their frames turned apart.
ARMS_FRAME = at('R724', (1, 0.5, 0.75))
SLOT_FRAME = at('R815', (2, 1, 0.5))
OPEN_FRAME = at('R125', (-1, 2, 0.25))


def plate():
    """The stadium's section: `24 + 4 pi`."""
    return stadium_area(6, 2)


def arms_volume():
    """The plate `3 (24 + 4 pi)` and the diamond over `[1, 2]` outside it:
    the diamond `21 / 2` less its part inside the plate. With `u = x - 6`,
    the diamond's half widths `7/4 (u + 3/2)` (`u <= 1/2`) and `7/2 (3/2 -
    u)`, the plate's `2` (`u <= 0`) and `sqrt(4 - u^2)`: inside, `16/7` for
    `x` in `[9/2, 9/2 + 8/7]` (the diamond's), `10/7` to `x = 6` (the
    plate's), the half disc's `u* sqrt(4 - u*^2) + 4 asin(u*/2)` to `u*` and
    the diamond's `7/2 (3/2 - u*)^2` beyond, `u* = (147 - 2 sqrt 407) / 106`
    the diamond's edge on the arc; in the chart, times the stored axes'
    determinant."""
    u = (147-2*mp.sqrt(407))/106
    inside = q(F(26, 7))+u*mp.sqrt(4-u**2)+4*mp.asin(u/2)+q(F(7, 2))*(q(F(3, 2))-u)**2
    return (3*plate()+q(F(21, 2))-inside)*det(ARMS_FRAME)


def slot_volume():
    """The plate `3 (24 + 4 pi)` less the pocket over its depth 1: the
    plate's part in the disc of radius 3 about `(-3/2, 0)` (its left half disc
    `2 pi`, and with `v = x + 3/2` the band `|y| <= 2` to `v = sqrt 5`, `4
    sqrt 5 - 6`, then the disc's `9 pi / 2 - 2 sqrt 5 - 9 asin(sqrt 5 / 3)`)
    less the post's disc `pi`."""
    pocket = 2*mp.pi+2*mp.sqrt(5)-6+9*mp.pi/2-9*mp.asin(mp.sqrt(5)/3)-mp.pi
    return (3*plate()-pocket)*det(SLOT_FRAME)


def open_volume():
    """The plate `5/2 (24 + 4 pi)` less the pocket over its depth 3/2: the
    plate's outline left of `x = 2` (`2 pi + 8`) less the post, the hole's
    square within the plate (`2` for `x` in `[0, 1]`, the half disc's strip
    `|y| <= 1`, `sqrt 3 + 2 pi / 3`)."""
    post = 2+mp.sqrt(3)+2*mp.pi/3
    return (q(F(5, 2))*plate()-q(F(3, 2))*(2*mp.pi+8-post))*det(OPEN_FRAME)


class Annulus:
    """A prism of a disc with a round hole for OCCT (`native`: one prism of
    both circles) and the reference's same set, the disc's cylinder less the
    hole's (past both caps)."""

    def __init__(self, frame, h0, h1, outer, hole):
        cx, cy, r = outer
        hx, hy, hr = hole
        self.native = prism([disc(cx, cy, r), disc(hx, hy, hr)], frame, h0, h1)
        self.parts = [cylinder(frame, cx, cy, r, h0, h1), cylinder(frame, hx, hy, hr, h0-1.0, h1+1.0)]
        self.op = 'cut'


def diamond():
    """The arms' diamond, counter-clockwise."""
    return Boundary(points=[(4.5, 0.0), (6.5, -3.5), (7.5, 0.0), (6.5, 3.5)])


def make_bodies():
    out = [
        Body('arms', prism([stadium_at(0, 0, 6, 2)], ARMS_FRAME, 0.0, 3.0), 'fuse',
             prism([diamond()], ARMS_FRAME, 1.0, 2.0),
             closed=lambda: (arms_volume(), None)),
        Body('slot', prism([stadium_at(0, 0, 6, 2)], SLOT_FRAME, 0.0, 3.0), 'cut',
             Annulus(SLOT_FRAME, 1.0, 2.0, (-1.5, 0.0, 3.0), (-0.5, 0.25, 1.0)),
             closed=lambda: (slot_volume(), None)),
        Body('open', prism([stadium_at(0, 0, 6, 2)], OPEN_FRAME, 0.0, 2.5), 'cut',
             Ring(OPEN_FRAME, 1.0, 4.0, (-4.0, -3.0, 2.0, 3.0), (-2.5, -1.0, 1.0, 1.0)),
             closed=lambda: (open_volume(), None)),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its own
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(trees.shifted(b.ref_expr, k))
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


ARMS_ROD = cylinder(ARMS_FRAME, 5.75, 2.25, 0.5, -1.0, 4.0)
SLOT_ROD = cylinder(SLOT_FRAME, 1.0, -1.0, 0.75, -1.0, 4.0)
ROD_SLOT = cylinder(SLOT_FRAME, -0.75, 1.0, 0.875, -1.0, 4.0)
OPEN_ROD = cylinder(OPEN_FRAME, 2.0, 1.0, 0.5, -1.0, 4.0)
CHAIN_ROD = cylinder(ARMS_FRAME, 3.0, 0.5, 0.75, -1.0, 4.0)
CHAIN_BOX = block(ARMS_FRAME, 5.5, 1.25, 6.75, 4.0, 1.5, 3.5)
SLOT_TOUCH = sphere(1.0, at('XY', point(SLOT_FRAME, 4.0, 0.0, 4.0)))
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('arms_slot', 'both', [imported('arms'), imported('slot')], three)
    out += group('chain_arms', 'chain', [imported('arms'), CHAIN_ROD, CHAIN_BOX], three, first='cut')
    out += group('arms_rod', 'pieces', [imported('arms'), ARMS_ROD], three)
    out += group('slot_rod', 'pieces', [imported('slot'), SLOT_ROD], three)
    out += group('rod_slot', 'pieces', [ROD_SLOT, imported('slot')], three)
    out += group('open_rod', 'pieces', [imported('open'), OPEN_ROD], three)
    out += group('slot_touch', 'pieces', [imported('slot'), SLOT_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
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
            assert sum(s[0] == 'imported' for s in c.items) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone'), c.name


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3c.3a's, each body's closed form)."""
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
    # Each imported body one solid, its closed form.
    for it, e in zip(items, sub):
        if it[0] != 'imported':
            continue
        assert forms.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
        V, A = BODIES[it[1]].closed()
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
    cache = os.environ.get('FLUSH_LEAVES_CACHE')
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
    """The writer's blocks: each body's rows under `write NAME PATH`."""
    out = []
    for b in BODIES.values():
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
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q_).hex() for q_ in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_flush_leaves_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its curved faces' kinds its prisms' arcs'
    (compared as sets, at least two curved faces), each stored vertex on the
    constructions' surfaces within 1e-12 of the body's size. Returns the
    largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        curved = sorted(set(k for k in kinds if k != 'plane'))
        assert curved == sorted(set(b.curved_kinds())), (b.name, kinds)
        assert len([k for k in kinds if k != 'plane']) >= 2, (b.name, kinds)
        surfs = []
        for spec in b.specs():
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
    parser.add_argument('--bodies', action='store_true', help='write the bodies\' blocks alone (no reference)')
    args = parser.parse_args()
    if args.cache:
        os.environ['FLUSH_LEAVES_CACHE'] = args.cache
    if args.bodies:
        (ROOT/'fixtures'/BODIES_FILE).write_text(bodies_text())
        return
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
