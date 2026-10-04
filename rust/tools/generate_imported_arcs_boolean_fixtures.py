#!/usr/bin/env python3
"""Fixtures for S9e.4b.1 of REVIEW_NOTES.md: imported prisms whose arcs'
ends round off their circles in their caps' frames, given to Booleans and
decided on S9e.4a's construction with each arc's ends taken onto its circle.

The bodies are OCCT's own output, as S9e.4a's
(`generate_imported_boolean_fixtures.py`): `boolean-imported-arcs-bodies.txt`
holds one `write NAME imported/NAME.brep` block per body for
`occt_boolean_oracle.cpp` (`MakePrism` of a profile face whose world
coordinates are its turned frame's roundings), which
`compare_imported_arcs_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/NAME.brep` (format version 1, no triangulations).
Nothing here reads the kernel; the files are read only to check that each is
the body it claims to be and that it is this step's (below).

`boolean-imported-arcs-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-imported-arcs-expected.tsv` gives
per case `expect KIND S9e.4b.1 CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a later sub-step's;
the class `prism`: an imported prism and a construction, `both`: two
imported inputs, `chain`: an imported input's result given to another
Boolean), a degenerate or unsupported case's `reason TEXT`, then `result N
volume area cx cy cz` (totals over the N solids, world coordinates) or
`empty`; `boolean-imported-arcs-frames.tsv` the stored axes of every solid
the kernel builds from its rows.

The bodies (each a `MakePrism` of a profile on a frame turned from the
world's axes, so its vertices' coordinates in its cap's frame are
roundings): `slot`, a stadium of straight length 6 and radius 2 (its lines
tangent to its arcs) in the `TILT` frame at `(5, 3, 0)`, height 5;
`halves_turn`, a circle of radius 5/2 about `(1/2, 3/4)` as two arcs (a
joint of two arcs of one circle, a seam-split cylinder) in the `TURN30`
frame at `(5, 3, 0)`, height 5; `rounded`, the rectangle `[-4, 4] x [-5/2, 5/2]` with
corners rounded by fillets of radius 1 (eight tangent joints) in the `R125`
frame at `(5, 5, 0)`, height 3; `notch`, `(0, -4)` to `(8, -4)`, an arc
about `(5, 0)` of radius 5 to `(8, 4)` (crossing the lines at an angle), back
to `(0, 4)`, in the `TURN30` frame at `(9/4, 1/8, -1/2)`, height 4; `lens`, two
arcs of circles of radius 5 about `(-3, 0)` and `(3, 0)` meeting at `(0,
-4)` and `(0, 4)`, in the `TURN30` frame at `(5, 3, 1/2)`, height 4 (two
circles at a joint: S9e.4b.4's).

Cases (each the three operations): `slot_box`, the slot and a box across
one of its arc ends; `box_slot`, the box and the slot as the tool;
`slot_rod`, the slot and a vertical rod through it; `slot_ball`, the slot
and a ball across its arc wall; `halves_box`, the split circle and a box
across its seam; `rounded_slab`, the rounded rectangle and a `TILT` slab;
`notch_box`, the notch and a box across its arc's crossing joint; `both`,
the notch and the slot, both imported; `chain_slot`, the slot less the rod,
then with the box. Declared `degenerate`: `slot_flush`, a box in the `TILT`
frame on the slot's flat wall's plane (two faces within the resolution of
one plane); `slot_kiss`, a box in the `TILT` frame whose wall touches the
slot's arc wall along a generatrix. Declared `unsupported`: `lens_box`, the
lens (S9e.4b.4) and the box.

The reference is the construction OCCT was given, through S9e.3a's chained
reference with S9e.4a's checks (`generate_imported_boolean_fixtures.
evaluate_chain`: the two families, each solid's closed form, the pair
identities, Monte Carlo, solid counts by rays, the meetings' sines and the
events' spacing outside the declared cases). With `--check`, where the
bodies' files exist, each is read (`stored_records`): its faces' kinds the
construction's, every stored vertex within 1e-12 of the size on the
construction's surfaces, and its arcs' ends, their stored vertices' local
coordinates in the construction's frame (`stored_axes`, the exact affine
map of its binary64 axes) rounded once to binary64, off their circles
exactly (as the kernel's S9e.4a construction takes them in the cap's
stored frame, a rounding of this one): the reason each body is this
step's.
"""
import argparse
from fractions import Fraction as F
import os
import struct

import mpmath as mp

from identity_reference import Boundary, Case
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
from generate_curved_boolean_fixtures import square, disc, stadium

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('prism', 'both', 'chain')
BODIES_FILE = 'boolean-imported-arcs-bodies.txt'
STEP = 'S9e.4b.1'
at, prism, sphere, construction = base.at, base.prism, base.sphere, base.construction


# ------------------------------------------------------------------ bodies

def halves(cx, cy, r):
    """A circle as two arcs (two faces on one cylinder)."""
    return Boundary(points=[(cx+r, cy), (cx-r, cy)], segments=[(cx, cy, r, True), (cx, cy, r, True)])


def rounded(x0, y0, x1, y1, f):
    """The rectangle `[x0, x1] x [y0, y1]` with its corners rounded by
    fillets of radius `f` (each line tangent to its arcs)."""
    return Boundary(points=[(x0+f, y0), (x1-f, y0), (x1, y0+f), (x1, y1-f), (x1-f, y1), (x0+f, y1), (x0, y1-f),
                            (x0, y0+f)],
                    segments=[None, (x1-f, y0+f, f, True), None, (x1-f, y1-f, f, True), None, (x0+f, y1-f, f, True),
                              None, (x0+f, y0+f, f, True)])


def notch():
    """`(0, -4)` to `(8, -4)`, the arc about `(5, 0)` of radius 5 to `(8, 4)`
    (meeting the lines at an angle), back to `(0, 4)`."""
    return Boundary(points=[(0.0, -4.0), (8.0, -4.0), (8.0, 4.0), (0.0, 4.0)],
                    segments=[None, (5.0, 0.0, 5.0, True), None, None])


def lens():
    """Two arcs of radius 5 about `(-3, 0)` and `(3, 0)` meeting at `(0, -4)`
    and `(0, 4)`."""
    return Boundary(points=[(0.0, -4.0), (0.0, 4.0)], segments=[(-3.0, 0.0, 5.0, True), (3.0, 0.0, 5.0, True)])


class Body(base.Body):
    pass


BODIES = {b.name: b for b in [
    Body('slot', prism([stadium(3.0, 2.0)], at('TILT', (5, 3, 0)), 0.0, 5.0), 'prism'),
    Body('halves_turn', prism([halves(0.5, 0.75, 2.5)], at('TURN30', (5, 3, 0)), 0.0, 5.0), 'prism'),
    Body('rounded', prism([rounded(-4.0, -2.5, 4.0, 2.5, 1.0)], at('R125', (5, 5, 0)), 0.0, 3.0), 'prism'),
    Body('notch', prism([notch()], at('TURN30', (2.25, 0.125, -0.5)), 0.0, 4.0), 'prism'),
    Body('lens', prism([lens()], at('TURN30', (5, 3, 0.5)), 0.0, 4.0), 'prism'),
]}


def imported(name):
    return ('imported', name)


def resolve(spec):
    """A spec's construction: an imported body's, or the spec itself."""
    return BODIES[spec[1]].spec if spec[0] == 'imported' else spec


BOX = prism([square(7.0, 0.0, 12.0, 10.0)], at('XY', (0, 0, 1)), 0.0, 1.5)
ROD = prism([disc(3.5, 4.5, 0.8)], at('XY', (0, 0, -3)), 0.0, 11.0)


class Imported:
    """One case: its solids (`imported(NAME)` or a construction), the first
    Boolean, the fixed stages and the last stage's operation (as S9e.4a's)."""

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

    @property
    def last(self):
        return self.stages[-1][0] if self.stages else self.op1

    def encode(self):
        from identity_reference import encode_boolean_case, encode_chained_case
        c = self.solid_cases
        if not self.stages:
            return encode_boolean_case(c[0], self.op1, c[1], BOOLEANS[0])
        (op2, sw2) = self.stages[0]
        return encode_chained_case(c[0], self.op1, c[1], BOOLEANS[0], op2, c[2], BOOLEANS[1], sw2, None)

    def native(self):
        from identity_reference import native_boolean_case, native_chained_case
        c = self.solid_cases
        if not self.stages:
            return native_boolean_case(c[0], self.op1, c[1])
        (op2, sw2) = self.stages[0]
        return native_chained_case(c[0], self.op1, c[1], op2, c[2], sw2, None)


def group(name, klass, specs, outcomes, first=None, reason=None, kind=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty)."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('degenerate' if reason else ('solid' if n else 'empty'))
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, specs, first or op, stages, k, n, reason))
    return out


FLUSH = ('a box in the slot\'s frame on the plane of its flat wall (coplanar in the construction, within the '
         'resolution of one plane in the imported body)')
KISS = 'a box\'s wall touching the imported slot\'s arc wall along a generatrix'
LENS = 'S9e.4b.4: the imported prism\'s arcs of two circles meeting at a joint'


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('slot_box', 'prism', [imported('slot'), BOX], three)
    out += group('box_slot', 'prism', [BOX, imported('slot')], three)
    out += group('slot_rod', 'prism', [imported('slot'), ROD], three)
    out += group('slot_ball', 'prism', [imported('slot'), sphere(1.5, at('XY', (10.25, 4.625, 2.125)))], three)
    out += group('halves_box', 'prism', [imported('halves_turn'), prism([square(6.5, 0.0, 9.0, 10.0)],
                                                                       at('XY', (0, 0, 1)), 0.0, 2.0)], three)
    out += group('rounded_slab', 'prism', [imported('rounded'), prism([square(-2.0, -15.0, 12.0, 15.0)],
                                                                      at('TILT', (0, 0, 2.5)), 0.0, 1.0)],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('notch_box', 'prism', [imported('notch'), prism([square(4.0, 5.0, 8.0, 9.0)],
                                                                 at('XY', (0, 0, 1)), 0.0, 2.0)], three)
    out += group('both', 'both', [imported('notch'), imported('slot')], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('chain_slot', 'chain', [imported('slot'), ROD, BOX], three, first='cut')
    out += group('slot_flush', 'prism', [imported('slot'), prism([square(-1.0, 2.0, 1.0, 4.0)],
                                                                 at('TILT', (5, 3, 0)), 1.0, 3.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=FLUSH)
    out += group('slot_kiss', 'prism', [imported('slot'), prism([square(5.0, -1.0, 8.0, 1.0)],
                                                                at('TILT', (5, 3, 0)), 1.0, 3.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=KISS)
    out += group('lens_box', 'prism', [imported('lens'), prism([square(4.0, -2.0, 6.0, 12.0)], at('XY', (0, 0, 1)),
                                                               0.0, 1.5)],
                 three, kind='unsupported', reason=LENS)
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
    cache = os.environ.get('IMPORTED_ARCS_CACHE')
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


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, degenerate in results:
        for key, value in checks.items():
            note(key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not degenerate:
            note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
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
        for k, (spec, c) in enumerate(zip(case.specs, case.constructions)):
            if spec[0] == 'imported':
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-imported-arcs-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-imported-arcs-expected.tsv': '\n'.join(out)+'\n',
            'boolean-imported-arcs-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_imported_arcs_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


# ------------------------------------------------------------------ the files' check

def local_map(frame):
    """The construction's exact affine map on its stored axes: a world point's
    local coordinates `(u, v, w)`, `p = o + u x + v y + w n` (Cramer's rule in
    rationals)."""
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))

    def det(a, b, c):
        return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))
    d = det(x, y, n)

    def local(p):
        q = tuple(F(p[i])-o[i] for i in range(3))
        return (det(q, y, n)/d, det(x, q, n)/d, det(x, y, q)/d)
    return local


def arc_ends(spec):
    """A prism's profile points ending or starting an arc, each with its arcs'
    circles `((cx, cy), r)` as binary64."""
    out = []
    for b in spec[1]:
        if b.circle is not None or b.segments is None:
            continue
        n = len(b.points)
        for j, p in enumerate(b.points):
            circles = [(s[0], s[1], s[2]) for s in (b.segments[(j-1) % n], b.segments[j]) if s is not None]
            if circles:
                out.append((p, circles))
    return out


def off_circles(body, vertices):
    """`(off, ends)`: of the stored vertices at the profile's arcs' ends, how
    many lie off one of their circles once their local coordinates in the
    construction's frame are rounded once to binary64 (exactly, in
    rationals), and how many there are."""
    local = local_map(body.spec[2])
    off = ends = 0
    near = []
    for v in vertices:
        u = local(v)
        near.append((float(u[0]), float(u[1])))
    for p, circles in arc_ends(body.spec):
        hits = [q for q in near if abs(q[0]-p[0]) <= 1e-9 and abs(q[1]-p[1]) <= 1e-9]
        assert len(hits) == 2, (body.name, p, hits)       # the bottom's and the top's
        for q in hits:
            ends += 1
            on = all((F(q[0])-F(cx))**2+(F(q[1])-F(cy))**2 == F(r)**2 for cx, cy, r in circles)
            off += not on
    return off, ends


def check_bodies():
    """Each existing body file: its stored surfaces' kinds the construction's
    faces', each stored vertex on the construction's surfaces within 1e-12
    of the body's size, and its arcs' ends off their circles. Returns the
    largest deviation, the files read and each body's `(off, ends)`."""
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
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            d = min(abs(s.f(X))/max(base.tc.norm(s.grad(X)), mp.mpf(10)**-30) for s in surfs)
            worst = max(worst, d/size)
        offs[b.name] = off_circles(b, vertices)
        assert offs[b.name][0] > 0, (b.name, 'every arc end on its circle: not this step\'s body')
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
        os.environ['IMPORTED_ARCS_CACHE'] = args.cache
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
          'size of its construction\'s surfaces; arcs\' ends off their circles in the construction\'s frame: '
          + ', '.join(f'{k} {o} of {e}' for k, (o, e) in offs.items()))


if __name__ == '__main__':
    main()
