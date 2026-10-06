#!/usr/bin/env python3
"""Fixtures for S9e.4b.4a of REVIEW_NOTES.md: imported prisms whose arcs of
two circles meet at a joint (a profile's corner where two circles cross or
touch), given to Booleans and decided on S9e.4a's construction with each
such arc taken through its two ends.

The bodies are OCCT's own output, as S9e.4b.1's
(`generate_imported_arcs_boolean_fixtures.py`):
`boolean-imported-joints-bodies.txt` holds one `write NAME
imported/NAME.brep` block per body for `occt_boolean_oracle.cpp`
(`MakePrism` of a profile face whose world coordinates are its turned
frame's roundings), which `compare_imported_joints_boolean.py
--write-bodies` runs to write `rust/fixtures/imported/NAME.brep` (format
version 1, no triangulations). Nothing here reads the kernel; the files are
read only to check that each is the body it claims to be and that it is
this step's (below).

`boolean-imported-joints-cases.txt` lists each case in the Boolean
protocol, an imported input as its one `brep imported/NAME.brep` row, a
chain's further `then` rows as S9e.1's; `boolean-imported-joints-
expected.tsv` gives per case `expect KIND S9e.4b.4a CLASS` (`solid`,
`empty`, `degenerate`: a `Degenerate` of S9's rules; `unsupported`:
`OutOfDomain`, a later sub-step's; the class `prism`: an imported prism and
a construction, `both`: two imported inputs, `chain`: an imported input's
result given to another Boolean), a degenerate or unsupported case's
`reason TEXT`, then `result N volume area cx cy cz` (totals over the N
solids, world coordinates) or `empty`; `boolean-imported-joints-frames.tsv`
the stored axes of every solid the kernel builds from its rows.

The bodies, each a `MakePrism` of a profile on a frame turned about the
world's `z` (so its vertices' coordinates in its cap's frame are
roundings), every joint of two circles a rational common point of both in
the profile's frame (the construction exact, the stored joints its
roundings): `quad`, the common of four discs of radius 5 about `(+-1, 0)`
and `(0, +-1)`, its corners `(+-3, +-3)` (DRAW's `bcut_complex/E8` and
`bfuse_complex/D5` part, four discs' common), in the `R125` frame at `(5, 5,
0)`, height 4; `arch`, a base from `(-2, 0)` to `(2, 0)` and two arcs of
radius 5 about `(-3, 0)` and `(3, 0)` up to their common point `(0, 4)` (a
pointed arch: two crossing joints' arcs ending at lines too), in the
`TURN30` frame at `(5, 3, -1/2)`, height 3; `cam`, an arc about the origin
of radius 5 from `(-3, -4)` to `(5, 0)`, one about `(2, 0)` of radius 3
tangent inside it there on to `(2, 3)`, then lines by `(-3, 3)`
(`bfuse_complex/E1`'s internal fillet), in the `R125` frame at `(5, 5, 0)`, height 3; `blade`, an
arc about the origin of radius 25/2 from `(12, -7/2)` to `(12, 7/2)`, one
about `(6, 1)` of radius 13/2 crossing it there at 0.11 rad on to `(6,
15/2)`, then lines by `(0, 15/2)` and `(0, -7/2)` (`bugs/modalg_2/bug4993`'s
crossing at a small angle), in the `TURN30` frame at the origin, height 3;
and `split_lens`, S9e.4b.1's lens (arcs of radius 5 about `(-+3, 0)` meeting
at `(0, -+4)`) with each arc split in two at its middle point, in the
`TURN30` frame at `(5, 3, 1/2)`, height 4 (declared `unsupported`: a joint
of two circles each holding several arcs, S9e.4b.4's).

Cases (each the three operations): `quad_box`, the quad and a box across
its corner `(3, 3)`; `box_quad`, the box and the quad as the tool;
`quad_rod`, the quad and a vertical rod through its arc wall; `quad_ball`,
the quad and a ball about its corner's vertical edge; `arch_box`, the arch
and a box across its apex; `arch_slab`, the arch and a `TILT` slab;
`cam_box`, the cam and a box across its tangent joint; `cam_rod`, the cam
and a rod through its small arc's wall; `blade_box`, the blade and a box
across its crossing joint; `both`, the quad and the arch, both imported;
`chain_quad`, the quad less the rod, then with the box. Declared
`degenerate`: `quad_seat`, a rod in the quad's frame on the circle of its
arc about `(-1, 0)` below its bottom face (DRAW's `bfuse_complex/D5`: its
tool continuing the part's wall), within rounding of the circle the arc is
taken through. Declared `unsupported`: `split_lens_box`, the split lens and
a box.

Two circles touching from either side at a joint (an S curve, DRAW's
`bcut_complex/P4`) make a profile that is not convex, which the chained
reference does not take: the kernel's tests write such a body of its own
and the DRAW trial reads `P4`'s.

The reference is the construction OCCT was given, through S9e.3a's chained
reference with S9e.4a's checks (`generate_imported_boolean_fixtures.
evaluate_chain`: the two families, each solid's closed form, the pair
identities, Monte Carlo, solid counts by rays, the meetings' sines and the
events' spacing outside the declared cases). With `--check`, where the
bodies' files exist, each is read (`stored_records`): its faces' kinds the
construction's, every stored vertex within 1e-12 of the size on the
construction's surfaces, and each joint of two circles, its stored
vertices' local coordinates in the construction's frame (the exact affine
map of its binary64 axes) rounded once to binary64, off one of its circles
at least (as the kernel's S9e.4a construction takes them in the cap's
stored frame, a rounding of this one): the reason each body is this
step's.
"""
import argparse
from fractions import Fraction as F
import os
import struct

import mpmath as mp

from identity_reference import Boundary
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
import generate_imported_arcs_boolean_fixtures as arcs
from generate_curved_boolean_fixtures import square, disc

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('prism', 'both', 'chain')
BODIES_FILE = 'boolean-imported-joints-bodies.txt'
PREFIX = 'boolean-imported-joints'
STEP = 'S9e.4b.4a'
at, prism, sphere, construction = base.at, base.prism, base.sphere, base.construction
local_map = arcs.local_map


# ------------------------------------------------------------------ bodies

def quad():
    """The common of four discs of radius 5 about `(+-1, 0)` and `(0, +-1)`:
    its corners `(+-3, +-3)`, each the common point of two neighbours'
    circles, every side an arc of the disc about the opposite centre."""
    return Boundary(points=[(3.0, -3.0), (3.0, 3.0), (-3.0, 3.0), (-3.0, -3.0)],
                    segments=[(-1.0, 0.0, 5.0, True), (0.0, -1.0, 5.0, True), (1.0, 0.0, 5.0, True),
                              (0.0, 1.0, 5.0, True)])


def arch():
    """A pointed arch: a base from `(-2, 0)` to `(2, 0)`, an arc about `(-3,
    0)` of radius 5 up to `(0, 4)` and one about `(3, 0)` of radius 5 down
    to `(-2, 0)`, their common point the apex."""
    return Boundary(points=[(-2.0, 0.0), (2.0, 0.0), (0.0, 4.0)],
                    segments=[None, (-3.0, 0.0, 5.0, True), (3.0, 0.0, 5.0, True)])


def cam():
    """An arc about the origin of radius 5 from `(-3, -4)` to `(5, 0)`, one
    about `(2, 0)` of radius 3 tangent inside it there on to `(2, 3)`, then
    lines by `(-3, 3)` back."""
    return Boundary(points=[(-3.0, -4.0), (5.0, 0.0), (2.0, 3.0), (-3.0, 3.0)],
                    segments=[(0.0, 0.0, 5.0, True), (2.0, 0.0, 3.0, True), None, None])


def blade():
    """An arc about the origin of radius 25/2 from `(12, -7/2)` to `(12,
    7/2)`, one about `(6, 1)` of radius 13/2 crossing it there (radii along
    `(24, 7)` and `(12, 5)`: 0.11 rad apart) on to `(6, 15/2)`, then lines by
    `(0, 15/2)` and `(0, -7/2)` back."""
    return Boundary(points=[(12.0, -3.5), (12.0, 3.5), (6.0, 7.5), (0.0, 7.5), (0.0, -3.5)],
                    segments=[(0.0, 0.0, 12.5, True), (6.0, 1.0, 6.5, True), None, None, None])


def split_lens():
    """S9e.4b.1's lens (arcs of radius 5 about `(-3, 0)` and `(3, 0)` meeting
    at `(0, -4)` and `(0, 4)`), each arc split in two at its middle point."""
    return Boundary(points=[(0.0, -4.0), (2.0, 0.0), (0.0, 4.0), (-2.0, 0.0)],
                    segments=[(-3.0, 0.0, 5.0, True), (-3.0, 0.0, 5.0, True), (3.0, 0.0, 5.0, True),
                              (3.0, 0.0, 5.0, True)])


class Body(base.Body):
    pass


BODIES = {b.name: b for b in [
    Body('quad', prism([quad()], at('R125', (5, 5, 0)), 0.0, 4.0), 'prism'),
    Body('arch', prism([arch()], at('TURN30', (5, 3, -0.5)), 0.0, 3.0), 'prism'),
    Body('cam', prism([cam()], at('R125', (5, 5, 0)), 0.0, 3.0), 'prism'),
    Body('blade', prism([blade()], at('TURN30', (0, 0, 0)), 0.0, 3.0), 'prism'),
    Body('split_lens', prism([split_lens()], at('TURN30', (5, 3, 0.5)), 0.0, 4.0), 'prism'),
]}


def imported(name):
    return ('imported', name)


def resolve(spec):
    """A spec's construction: an imported body's, or the spec itself."""
    return BODIES[spec[1]].spec if spec[0] == 'imported' else spec


BOX = prism([square(6.0, 6.0, 12.0, 12.0)], at('XY', (0, 0, 1)), 0.0, 1.5)
ROD = prism([disc(8.5, 4.75, 1.0)], at('XY', (0, 0, -3)), 0.0, 11.0)


class Imported(arcs.Imported):
    """One case (S9e.4b.1's), its imported bodies this step's."""

    def __init__(self, name, klass, specs, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.specs = specs
        self.constructions = [construction(resolve(s), OPERATIONS[k]) for k, s in enumerate(specs)]
        self.solid_cases = []
        for k, s in enumerate(specs):
            c = construction(resolve(s), OPERATIONS[k])
            if s[0] == 'imported':
                c.brep = BODIES[s[1]].path
            c.name = name
            self.solid_cases.append(c)
        for c in self.constructions:
            c.name = name


def group(name, klass, specs, outcomes, first=None, reason=None, kind=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty)."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('degenerate' if reason else ('solid' if n else 'empty'))
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, specs, first or op, stages, k, n, reason))
    return out


SEAT = ('a rod in the quad\'s frame on the circle of its arc about (-1, 0) below its bottom face, within rounding of '
        'the circle the arc is taken through (two surfaces within the resolution of one)')
SPLIT = 'S9e.4b.4: the imported prism\'s joint of two circles each holding several arcs'


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('quad_box', 'prism', [imported('quad'), BOX], three)
    out += group('box_quad', 'prism', [BOX, imported('quad')], three)
    out += group('quad_rod', 'prism', [imported('quad'), ROD], three)
    out += group('quad_ball', 'prism', [imported('quad'), sphere(1.5, at('XY', (6.375, 8.625, 2.0)))], three)
    out += group('arch_box', 'prism', [imported('arch'), prism([square(2.0, 5.5, 4.0, 7.5)], at('XY', (0, 0, 0.5)),
                                                              0.0, 1.0)], three)
    out += group('arch_slab', 'prism', [imported('arch'), prism([square(-2.0, -15.0, 12.0, 15.0)],
                                                               at('TILT', (0, 0, 2.25)), 0.0, 1.0)],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('cam_box', 'prism', [imported('cam'), prism([square(8.5, 5.5, 12.0, 9.0)], at('XY', (0, 0, 1)),
                                                            0.0, 1.0)], three)
    out += group('cam_rod', 'prism', [imported('cam'), prism([disc(7.875, 8.375, 0.5)], at('XY', (0, 0, -2)),
                                                            0.0, 7.0)], three)
    out += group('blade_box', 'prism', [imported('blade'), prism([square(7.5, 8.0, 10.0, 10.5)],
                                                                at('XY', (0, 0, 1)), 0.0, 1.0)], three)
    out += group('both', 'both', [imported('quad'), imported('arch')], three)
    out += group('chain_quad', 'chain', [imported('quad'), ROD, BOX], three, first='cut')
    out += group('quad_seat', 'prism', [imported('quad'), prism([disc(-1.0, 0.0, 5.0)], at('R125', (5, 5, 0)),
                                                               -3.0, 0.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=SEAT)
    out += group('split_lens_box', 'prism', [imported('split_lens'), prism([square(4.0, -2.0, 6.0, 12.0)],
                                                                          at('XY', (0, 0, 1)), 0.0, 1.5)],
                 three, kind='unsupported', reason=SPLIT)
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.specs) == 2+len(c.stages), c.name
        assert c.last in OPS, c.name
        assert any(s[0] == 'imported' for s in c.specs), c.name
        if c.klass == 'both':
            assert sum(s[0] == 'imported' for s in c.specs) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        for s in c.specs:
            r = resolve(s)
            frame = r[{'sphere': 2, 'prism': 2}[r[0]]]
            assert tuple(frame[3:]) in base.FRAMES.values(), c.name


# ------------------------------------------------------------------ the reference

def jobs(mc_n):
    """One job per group for `base.evaluate_chain`, its specs resolved."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [[resolve(s) for s in c.specs], c.op1 if c.stages else None,
                                        bool(c.stages), {}, c.kind in ('degenerate', 'unsupported')])
        e[3][c.last] = c.solids
    return [(name, specs, op1, stages, ops, degenerate, mc_n, 'chain')
            for name, (specs, op1, stages, ops, degenerate) in groups.items()]


def cached(job):
    """`base.evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('IMPORTED_JOINTS_CACHE')
    if not cache:
        return base.evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+inspect.getsource(base.evaluate_chain)
                          + Path(base.ref.__file__).read_text()).encode()).hexdigest()[:24]
    path = Path(cache)/f'{job[0]}-{key}.pickle'
    if path.exists():
        return pickle.loads(path.read_bytes())
    try:
        out = base.evaluate(job)
    except Exception as e:
        raise RuntimeError(f'{job[0]}: {e!r}') from e
    path.write_bytes(pickle.dumps(out))
    return out


reference_checks = arcs.reference_checks


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
        for k, (spec, c) in enumerate(zip(case.specs, case.constructions)):
            if spec[0] == 'imported':
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_imported_joints_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


# ------------------------------------------------------------------ the files' check

def joints(spec):
    """A prism's profile points where arcs of two different circles meet,
    each with both circles `((cx, cy), r)` as binary64."""
    out = []
    for b in spec[1]:
        if b.circle is not None or b.segments is None:
            continue
        n = len(b.points)
        for j, p in enumerate(b.points):
            before, after = b.segments[(j-1) % n], b.segments[j]
            if before is not None and after is not None and tuple(before[:3]) != tuple(after[:3]):
                out.append((p, [tuple(before[:3]), tuple(after[:3])]))
    return out


def on_circle(q, circle):
    cx, cy, r = circle
    return (F(q[0])-F(cx))**2+(F(q[1])-F(cy))**2 == F(r)**2


def off_joints(body, vertices):
    """`(off, joints)`: of the stored vertices at the profile's joints of two
    circles, how many lie off one of their circles at least once their local
    coordinates in the construction's frame are rounded once to binary64
    (exactly, in rationals: S9e.4b.1's refusal), and how many there are."""
    local = local_map(body.spec[2])
    near = []
    for v in vertices:
        u = local(v)
        near.append((float(u[0]), float(u[1])))
    off = count = 0
    for p, circles in joints(body.spec):
        hits = [q for q in near if abs(q[0]-p[0]) <= 1e-9 and abs(q[1]-p[1]) <= 1e-9]
        assert len(hits) == 2, (body.name, p, hits)       # the bottom's and the top's
        for q in hits:
            count += 1
            off += not all(on_circle(q, c) for c in circles)
    return off, count


def check_bodies():
    """Each existing body file: its stored surfaces' kinds the construction's
    faces', each stored vertex on one of the construction's walls and one of
    its caps within 1e-12 of the body's size, and its joints of two circles
    off one of them at least. Returns the largest deviation, the files read and each body's
    `(off, joints)`."""
    worst, read, offs = mp.mpf(0), 0, {}
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        want, surfs = base.construction_surfaces(b.spec)
        assert sorted(kinds) == sorted(want), (b.name, kinds, want)
        I = construction(b.spec, 0)
        size = max([1.0]+[abs(v) for v in I.frame[:3]]+[abs(x) for v in vertices for x in v])
        walls = [s for s in surfs if isinstance(s.label, tuple) and s.label[0] == 'wall']
        caps = [s for s in surfs if s not in walls]
        for v in vertices:
            # On a wall and on a cap (a vertex on a cap plane exactly, as
            # most are, would hide its wall's distance in one minimum).
            X = tuple(mp.mpf(c) for c in v)
            d = max(min(abs(s.f(X))/max(base.tc.norm(s.grad(X)), mp.mpf(10)**-30) for s in group)
                    for group in (walls, caps))
            worst = max(worst, d/size)
        offs[b.name] = off_joints(b, vertices)
        assert offs[b.name][0] > 0, (b.name, 'every joint of two circles on both: not this step\'s body')
    assert worst <= 1e-12, ('a stored vertex off its construction', mp.nstr(worst, 3))
    return worst, read, offs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per group')
    parser.add_argument('--only', help='one group, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each group\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['IMPORTED_JOINTS_CACHE'] = args.cache
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
        assert value <= base.LIMITS[key], (key, mp.nstr(value, 3))
    files = generate(results)
    for name, contents in files.items():
        target = ROOT/'fixtures'/name
        if args.check:
            if target.read_text() != contents:
                raise SystemExit(f'{target} is stale')
        else:
            target.write_text(contents)
    deviation, read, offs = check_bodies()
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
          'size of its construction\'s surfaces; joints of two circles off one of them at least once rounded in the '
          'construction\'s frame: ' + ', '.join(f'{k} {o} of {e}' for k, (o, e) in offs.items()))


if __name__ == '__main__':
    main()
