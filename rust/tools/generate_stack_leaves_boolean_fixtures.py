#!/usr/bin/env python3
"""Fixtures for S9e.4b.4b.2b.2b.3a of REVIEW_NOTES.md: imported bodies of
several cylinder faces and plane faces that are stacks of two prisms on one
frame, the upper one standing on the lower one's cap, both caps partly
exposed (two prism leaves touching along a face, their caps on one plane),
given to Booleans.

The bodies are OCCT's own output, as S9e.4b.4b.2b.2b.2's
(`generate_flush_leaves_boolean_fixtures.py`): `boolean-stack-leaves-bodies.txt`
holds one `write NAME imported/stack_NAME.brep` block per body for
`occt_boolean_oracle.cpp` (the first prism's rows, a `boolean` row and the
second's: the Boolean's one solid written), which
`compare_stack_leaves_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/stack_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-stack-leaves-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/stack_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-stack-leaves-expected.tsv` gives
per case `expect KIND S9e.4b.4b.2b.2b.3a CLASS` (`solid`, `empty`,
`degenerate`: a `Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a
later step's; the class `pieces`: an imported body and a construction,
`both`: two imported inputs, `chain`: an imported input's result given to
another Boolean), a degenerate or unsupported case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-stack-leaves-frames.tsv` the stored axes
of every solid the kernel builds from its rows.

Each body is one fuse OCCT was given of two prisms on one frame, the upper
one over heights starting at the lower one's top, every section a circle or
a line, every crossing of the two outlines on the shared plane rational, the
frames turned about the world's `z` (their normals the world's `z`
exactly): `step`, the box `[-3, 3] x [-2, 2]` over `[0, 3/2]` on `R724` and
the stadium of half discs of radius 5/4 about `(1, 1)` and `(5, 1)` over
`[3/2, 3]` standing on it across its corner (its outline crossing the box's
at `(3, -1/4)` and `(1/4, 2)`: the box's cap and the stadium's both partly
exposed on the plane `3/2`); `bar`, the stadium plate of half discs of
radius 25/8 about `(0, 0)` and `(6, 0)` over `[0, 2]` on `R125` and the box
`[5, 21/2] x [-15/8, 15/8]` over `[2, 3]` standing on it across its arc end
(the outlines crossing at `(17/2, +-15/8)`: the plate's cap exposed round
the box, the box's beyond the arc). A lower prism whose cap the upper one
hides wholly (a prism of one cap, its holes closed by the upper prism's cap)
is no evidence here: on such a frame OCCT writes its hole's cylinder capped
on the upper prism's cap plane exactly and the import takes it already; the
kernel's own such bodies on a tilted frame, whose caps its writer rounds
apart, are the kernel's tests' and the fuzz target's.

Cases (each the three operations): `step_rod`, the step and an upright rod
inside both outlines across the shared plane (its seam's middle on it),
through the box's bottom and the stadium's top; `step_wall`, the step and
an upright rod across the box's wall `x = 3` under the stadium; `bar_rod`,
the bar and an upright rod inside both outlines across the shared plane
(its seam's middle on it); `rod_bar`, an upright rod less the bar (the body
the tool) across its corner `(21/2, 15/8)`; `step_bar`, the step and the bar
both imported (their frames turned apart about `z`, the bar across the
step's box and stadium); `chain_step`, the step less an upright rod through
the box beside the stadium, then with a box on its frame across the box's
wall `x = 3`, the shared plane, the stadium's line `y = -1/4` and its right
arc. Declared `degenerate`: `bar_touch`, the bar and a ball of radius 1
resting on the box's top cap at a point, exactly (`Degenerate("a tangency
between the inputs (S9c)")`). The bodies' partners stand along their axes,
on their frames (normals the world's `z` exactly, so the import's leaves and
the rods are parallel exactly).

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`, which models convex
profiles): each imported body its fuse of two prisms, `(X op Y) op2 C`
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

from identity_reference import Case
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_piece_trees_boolean_fixtures as trees
import generate_primitive_chains_boolean_fixtures as chains
import generate_two_leaves_boolean_fixtures as two
import generate_rim_leaves_boolean_fixtures as rim

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-stack-leaves'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4b.2b.2b.3a'
FRAMES = dict(two.FRAMES, R125=(0.0, 0.0, 1.0, 12.0, 5.0, 0.0))
prism, sphere, construction = base.prism, base.sphere, base.construction
ref = base.ref
q, point = forms.q, forms.point
cylinder, block = forms.cylinder, forms.block
det = chains.det
stadium_at, stadium_area = rim.stadium_at, rim.stadium_area


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Body(rim.Body):
    """An imported body `first op second` of two prisms, OCCT's rows for it,
    its closed form `(volume, area)` (the area `None` where not checked)."""

    def __init__(self, name, first, op, second, closed=None, kind=None, reason=None):
        super().__init__(name, first, op, second, closed, kind, reason)
        self.path = f'imported/stack_{name}.brep'


# The bodies' frames (each body's two prisms on one), the bar across the
# step (`step_bar`), their frames turned apart.
STEP_FRAME = at('R724', (1, 0.5, 0.75))
BAR_FRAME = at('R125', (2, 4, 0.5))


def step_volume():
    """The box `24 * 3/2` and the stadium `3/2 (10 + 25 pi / 16)` (their
    common part a plane): `51 + 75 pi / 32`; in the chart, times the stored
    axes' determinant."""
    return (36+q(F(3, 2))*stadium_area(4, F(5, 4)))*det(STEP_FRAME)


def bar_volume():
    """The plate `2 (75/2 + 625 pi / 64)` and the box `11/2 * 15/4`: `75 +
    165/8 + 625 pi / 32`."""
    return (2*stadium_area(6, F(25, 8))+q(F(11, 2))*q(F(15, 4)))*det(BAR_FRAME)


def make_bodies():
    out = [
        Body('step', prism([forms.square(-3.0, -2.0, 3.0, 2.0)], STEP_FRAME, 0.0, 1.5), 'fuse',
             prism([stadium_at(1, 1, 4, F(5, 4))], STEP_FRAME, 1.5, 3.0),
             closed=lambda: (step_volume(), None)),
        Body('bar', prism([stadium_at(0, 0, 6, F(25, 8))], BAR_FRAME, 0.0, 2.0), 'fuse',
             prism([forms.square(5.0, -1.875, 10.5, 1.875)], BAR_FRAME, 2.0, 3.0),
             closed=lambda: (bar_volume(), None)),
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


STEP_ROD = cylinder(STEP_FRAME, 2.0, 1.0, 0.5, -1.0, 4.0)
STEP_WALL = cylinder(STEP_FRAME, 3.0, 0.75, 0.5, -1.0, 4.0)
BAR_ROD = cylinder(BAR_FRAME, 7.5, 1.0, 0.5, -1.0, 5.0)
ROD_BAR = cylinder(BAR_FRAME, 10.0, 1.5, 0.75, -1.0, 4.0)
CHAIN_ROD = cylinder(STEP_FRAME, -1.5, 0.0, 0.5, -1.0, 4.0)
CHAIN_BOX = block(STEP_FRAME, 2.5, -1.0, 6.0, 0.5, 1.0, 2.25)
BAR_TOUCH = sphere(1.0, at('XY', point(BAR_FRAME, 7.0, 0.0, 4.0)))
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('step_bar', 'both', [imported('step'), imported('bar')], three)
    out += group('chain_step', 'chain', [imported('step'), CHAIN_ROD, CHAIN_BOX], three, first='cut')
    out += group('step_rod', 'pieces', [imported('step'), STEP_ROD], three)
    out += group('step_wall', 'pieces', [imported('step'), STEP_WALL], three)
    out += group('bar_rod', 'pieces', [imported('bar'), BAR_ROD], three)
    out += group('rod_bar', 'pieces', [ROD_BAR, imported('bar')], three)
    out += group('bar_touch', 'pieces', [imported('bar'), BAR_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
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
    cache = os.environ.get('STACK_LEAVES_CACHE')
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
    """Every case's native rows (`compare_stack_leaves_boolean.py`)."""
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
        os.environ['STACK_LEAVES_CACHE'] = args.cache
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
