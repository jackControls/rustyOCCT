#!/usr/bin/env python3
"""Fixtures for S9e.4b.4b.2b.2b.3b of REVIEW_NOTES.md: imported bodies of
several cylinder faces and plane faces that are a prism leaf whose caps a
leaf crossing it splits (a stadium plate with a stadium slab through it on
crossing axes), given to Booleans.

The bodies are OCCT's own output, as S9e.4b.4c.2b.1's
(`generate_nested_leaves_boolean_fixtures.py`):
`boolean-crossed-leaves-bodies.txt` holds one `write NAME
imported/crossed_NAME.brep` block per body for `occt_boolean_oracle.cpp` (the
plate's rows, a `boolean fuse` row and the slab's: the Boolean's one solid
written), which `compare_crossed_leaves_boolean.py --write-bodies` runs to
write `rust/fixtures/imported/crossed_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-crossed-leaves-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/crossed_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-crossed-leaves-expected.tsv` gives
per case `expect KIND S9e.4b.4b.2b.2b.3b CLASS` (`solid`, `empty`,
`degenerate`: a `Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a
later step's; the class `pieces`: an imported body and a construction,
`both`: two imported inputs, `chain`: an imported input's result given to
another Boolean), a degenerate or unsupported case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-crossed-leaves-frames.tsv` the stored axes
of every solid the kernel builds from its rows.

Each body is one fuse OCCT was given of a stadium plate and a stadium slab
on crossing axes, every section a circle or a line: the plate on the
world's axes about its origin, the slab's normal the world's `x` and its `x`
the world's `y` (`ALONGX`: its chart the plate's `(y, z)`), so that the
plate's walls run along the slab's axis exactly and partners standing along
the plate's normal are parallel to its walls exactly. (On frames turned off
the world's axes the slab's in-plane axes are rounded: its stadium's joints
then lie a rounding apart from its arcs' extremes along the reference's
chords, each counted twice, `1e-18` of their range apart.) `joint`: the plate of half discs of radius 2 about
`(0, 0)` and `(4, 0)` over `[0, 3/2]` with the slab of half discs of radius 2
about the plate's `(y, z) = (-1, 3/4)` and `(3, 3/4)` over the plate's `x`
in `[-1/4, 3/4]` (the `boolean` fuzz target's measured body moved: the slab
across the plate's line-arc joints at `x = 0`, both plate
caps cut in two, the joints hidden in the slab, the slab whole); `waist`:
the plate of half discs of radius 2 about `(0, 0)` and `(8, 0)` with the
same slab over the plate's `x` in `[3, 4]` (across the plate's straight
middle: its line walls each one line through the slab).

Cases (each the three operations): `joint_rod`, the joint and an upright
rod on the plate's frame across the slab's cap `x = 3/4` and the plate's
right piece; `joint_arc`, the joint and an upright rod across the slab's
cap `x = -1/4` and the plate's arc piece; `waist_rod`, the waist and an
upright rod across the slab's cap `x = 4`; `rod_joint`, an upright rod
less the joint (the body the tool) across the slab and both plate pieces;
`joint_cyl`, the joint and S9e.4a's imported cylinder `cyl`
(`generate_imported_boolean_fixtures.py`: radius 3 about the world's `(5,
5)` over `[-1, 5]`) across the slab and the plate; `chain_waist`, the waist
less an upright rod across the slab, then with a box on the plate's frame
across the plate's right piece and the rod's hole. Declared `degenerate`:
`joint_touch`, the joint and a ball of radius 1 resting on the plate's top
at its point `(5/2, 0)`, exactly (`Degenerate("a tangency between the
inputs (S9c)")`).

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`, which models convex
profiles): each imported body the plate fused with the slab, `(X op Y) op2
C` (swapped, the body the tool), two bodies `(X op Y) op2 Z`, a chain `((X
op Y) op2 C) op3 D`, with S9e.4b.3c.3a's checks (the two families, each
solid's closed form, the pair identities on the last Boolean's arguments,
Monte Carlo, solid counts by rays at two resolutions, each body one solid,
the meetings' sines and the events' spacing outside the declared group) and
each body's closed form (the plate's and the slab's sections in their charts
times their stored axes' determinants, less their common part in the
plate's chart: the plate's section over the slab's range of its `x`, the
slab covering the plate's `(y, z)` there). With `--check`, where the
bodies' files exist, each is read (`stored_records`): its curved faces'
kinds its prisms' arcs', every stored vertex within 1e-12 of the size on the
constructions' surfaces.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case, native_case
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_piece_trees_boolean_fixtures as trees
import generate_primitive_chains_boolean_fixtures as chains
import generate_prism_leaves_boolean_fixtures as leaves
import generate_two_leaves_boolean_fixtures as two

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-crossed-leaves'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4b.2b.2b.3b'
FRAMES = two.FRAMES
prism, sphere, construction = base.prism, base.sphere, base.construction
ref = base.ref
q, point = forms.q, forms.point
cylinder, block = forms.cylinder, forms.block
det, stadium, stadium_area = chains.det, leaves.stadium, leaves.stadium_area


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Body(chains.Body):
    """An imported body `plate fuse slab`, OCCT's rows for it, its closed form
    `(volume, area)` (the area `None` where not checked)."""

    def __init__(self, name, plate, slab, closed):
        super().__init__(name, plate, 'fuse', slab, closed)
        self.path = f'imported/crossed_{name}.brep'


# The bodies' frames: each plate's, and its slab's at the plate's point `(x,
# -1, 3/4)` (its chart's origin the slab's half disc about the plate's `(y,
# z) = (-1, 3/4)`): the joint's slab over `[-1/4, 3/4]` about the plate's `x =
# 0`, the waist's over `[-1/2, 1/2]` about its `x = 7/2`.
JOINT_PLATE = at('XY', (5.25, 1.625, 0.5))
JOINT_SLAB = at('ALONGX', point(JOINT_PLATE, 0, -1, 0.75))
WAIST_PLATE = at('XY', (-2, 1, 0.25))
WAIST_SLAB = at('ALONGX', point(WAIST_PLATE, 3.5, -1, 0.75))
SLAB = stadium(4, 2.0)


def chord_area(x0, x1, r):
    """The disc of radius `r` about the origin between `x = x0` and `x =
    x1` (`-r <= x0 <= x1 <= r`): `[x sqrt(r^2 - x^2) + r^2 asin(x / r)]`."""
    R = q(r)

    def g(x):
        X = q(x)
        return X*mp.sqrt(R**2-X**2)+R**2*mp.asin(X/R)
    return g(x1)-g(x0)


def joint_volume():
    """The plate `3/2 (16 + 4 pi)` and the slab `16 + 4 pi`, less the
    plate's section over `x` in `[-1/4, 3/4]` (its half disc's from `-1/4` to
    `0`, the strip `3/4 x 4`) over `3/2`, each in its chart times its stored
    axes' determinant."""
    common = q(F(3, 2))*(chord_area(F(-1, 4), 0, 2)+3)
    return (q(F(3, 2))*stadium_area(4, 2)-common)*det(JOINT_PLATE)+stadium_area(4, 2)*det(JOINT_SLAB)


def waist_volume():
    """The plate `3/2 (32 + 4 pi)` and the slab `16 + 4 pi`, less the strip
    `1 x 4` over `3/2`."""
    return (q(F(3, 2))*stadium_area(8, 2)-6)*det(WAIST_PLATE)+stadium_area(4, 2)*det(WAIST_SLAB)


def make_bodies():
    out = [
        Body('joint', prism([stadium(4, 2.0)], JOINT_PLATE, 0.0, 1.5), prism([SLAB], JOINT_SLAB, -0.25, 0.75),
             closed=lambda: (joint_volume(), None)),
        Body('waist', prism([stadium(8, 2.0)], WAIST_PLATE, 0.0, 1.5), prism([SLAB], WAIST_SLAB, -0.5, 0.5),
             closed=lambda: (waist_volume(), None)),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


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
            exprs.append(trees.shifted(b.ref_expr, k))
        elif it[0] == 'base':
            exprs.append(len(specs))
            specs.append(base.BODIES[it[1]].spec)
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(pieces.Imported):
    """One case, its imported inputs this step's bodies or S9e.4a's
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


JOINT_ROD = cylinder(JOINT_PLATE, 1.0, 0.25, 0.5, -0.5, 2.5)
JOINT_ARC = cylinder(JOINT_PLATE, -0.625, 0.875, 0.5, -0.5, 2.5)
WAIST_ROD = cylinder(WAIST_PLATE, 4.625, -0.875, 0.75, -0.5, 2.5)
ROD_JOINT = cylinder(JOINT_PLATE, 0.25, 2.5, 1.0, 0.5, 1.0)
CHAIN_ROD = cylinder(WAIST_PLATE, 3.375, 0.875, 0.5, -0.5, 2.5)
CHAIN_BOX = block(WAIST_PLATE, 2.75, -0.25, 6.25, 1.25, 1.0, 2.0)
JOINT_TOUCH = sphere(1.0, at('XY', point(JOINT_PLATE, 2.5, 0.0, 2.5)))
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('joint_cyl', 'both', [imported('joint'), earlier('cyl')], three)
    out += group('chain_waist', 'chain', [imported('waist'), CHAIN_ROD, CHAIN_BOX], three, first='cut')
    out += group('joint_rod', 'pieces', [imported('joint'), JOINT_ROD], three)
    out += group('joint_arc', 'pieces', [imported('joint'), JOINT_ARC], three)
    out += group('waist_rod', 'pieces', [imported('waist'), WAIST_ROD], three)
    out += group('rod_joint', 'pieces', [ROD_JOINT, imported('joint')], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('joint_touch', 'pieces', [imported('joint'), JOINT_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
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
    cache = os.environ.get('CROSSED_LEAVES_CACHE')
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
    """Every case's native rows (`compare_crossed_leaves_boolean.py`)."""
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
        os.environ['CROSSED_LEAVES_CACHE'] = args.cache
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
