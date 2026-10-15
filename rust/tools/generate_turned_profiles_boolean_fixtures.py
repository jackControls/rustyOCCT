#!/usr/bin/env python3
"""Fixtures for S9e.4b.4d.2a of REVIEW_NOTES.md: imported turned bodies of
cylinders, cones and planes about one axis meeting at corners only (no
smooth join, no cone within the resolution of a plane or a cylinder), given
to Booleans (`bugs/modalg_6/bug28626`'s kind: a point, a waist, steps).

The bodies are OCCT's own output, as S9e.4b.4d.1's
(`write_turned_bodies.py`): `boolean-turned-profiles-bodies.txt` holds one
`write NAME imported/profile_NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a `plane` row, the profile's `wire` of lines from
the axis out and back to it, and a `revolve` row turning it a whole turn
about the axis: BRepPrimAPI_MakeRevol), which
`compare_turned_profiles_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/profile_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-turned-profiles-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/profile_NAME.brep` row, a
chain's further `then` rows as S9e.1's; `boolean-turned-profiles-
expected.tsv` gives per case `expect KIND S9e.4b.4d.2a CLASS` (`solid`,
`empty`, `degenerate`: a `Degenerate` of S9's rules; `unsupported`:
`OutOfDomain`, a later step's; the class `pieces`: an imported body and a
construction, `both`: two imported inputs, `chain`: an imported input's
result given to another Boolean), a degenerate or unsupported case's
`reason TEXT`, then `result N volume area cx cy cz` (totals over the N
solids, world coordinates) or `empty`; `boolean-turned-profiles-frames.tsv`
the stored axes of every solid the kernel builds from its rows.

Each body's profile is a path of dyadic points `(r, w)` (`r` the radius,
`w` the height along the axis) from a point on the axis out and back to it,
its heights never falling, every segment a line: rising, a cylinder or a
cone; level, a plane (a disc at an end, an annulus at a step); its ends an
apex or a disc. The bodies, each on the world's axes: `point`
(`bug28626_case1_shapes`' kind, about `z`): an apex, a cone of half-angle
pi/4 to radius 3 at 3, a cylinder to 6, a cone widening by 2/5 to radius 4
at 17/2, a cylinder to 16 and its disc; `point_x`, the same profile about
the world's `x` through `(-3/2, 33/8, 3/8)`; `waist` (about `z`): a disc of
radius 8, a cylinder to 10, a cone narrowing by 3/4 to radius 5 at 14 (a
convex corner below it, a concave one above), a cylinder to 20 and its disc;
`steps` (about `z`): a disc of radius 6, a cylinder to 4, a chamfer
narrowing by 1 to radius 5 at 5, a step in to radius 3 (an annulus facing
up), a cylinder to 9, a cone widening by 1/2 to radius 4 at 11 (a
concave corner), a cylinder to 13 and its disc; `crown` (about `z` through
`(1, 2, 0)`): a disc of radius 3, a cylinder to 4, a step out to radius 6
(an annulus facing down), a cylinder to 7, a cone narrowing by 1 to radius
4 at 9, a cylinder to 11 and its disc. No two of a body's cones share an apex or have
rulings of one direction (two cones of one half-angle opening apart), which
S9's rules refuse.

Cases (each the three operations): `point_box`, the point and a `TILT` box
across its first corner; `waist_rod`, the waist and a rod along the world's
`x` through its cone; `steps_ball`, the steps and a ball across the
annulus; `box_crown`, a `TILT` box less the crown (the body the tool) across
its step out; `point_tilt`, the point along `x` and S9e.4a's imported box
on the `TILT` frame both imported (a second turned body on a crossing axis,
`bug28626`'s configuration, or S9e.4a's imported cylinder crossing it, ran
over 30 worker-minutes in the reference); `chain_waist`, the waist less a box
biting its convex corner, then with a `TILT` slab across its neck. Declared `degenerate`: `waist_touch`, the waist and a ball resting on
its top disc at a point, exactly (`Degenerate("a tangency between the inputs
(S9c)")`).

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each body the max-min
form of its profile's lines (`Lattice`, this module's own: for each rising
segment the common of the primitives whose lines lie on or outside its own
over its heights, each such term then made as small as the form allows,
the fuse of those commons, each primitive a cylinder or a
frustum over a range of dyadic heights at which its radius is dyadic:
capped at a step's or an end's plane on its wider side, past a step's plane
into the wider side, to its apex, or past a corner over the heights it
bounds), checked exactly against the profile at every height, then
`(X) op C` (swapped, the body the tool), two bodies `(X) op (Y)`, a chain
`((X) op C) op2 D`, with S9e.4b.3c.3a's checks (the two families, each
solid's closed form, the pair identities on the last Boolean's arguments,
Monte Carlo, solid counts by rays at two resolutions, each body one solid,
the meetings' sines and the events' spacing outside the declared group) and
each body's closed form (Pappus: its frusta's volumes `pi h (a^2 + a b +
b^2) / 3` and its lateral, disc and annulus areas, in exact arithmetic but
for the slant heights' square roots). With `--check`, where the bodies'
files exist, each is read (`stored_records`): its faces' kinds its
profile's, every stored vertex within 1e-12 of the size on the profile's
surfaces.
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

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-turned-profiles'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2a'
FRAMES = forms.FRAMES
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
q = forms.q
cylinder, block = forms.cylinder, forms.block
# A whole turn, binary64 (`repr(2 * math.pi)`).
TURN = '6.283185307179586'


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def num(x):
    """A coordinate's shortest round-trip text (exact: every coordinate here
    is dyadic)."""
    x = F(x)
    assert float(x) == x, x
    return str(x.numerator) if x.denominator == 1 else repr(float(x))


# ------------------------------------------------------------------ the lattice

class Lattice:
    """A profile's max-min form: `points`, the dyadic `(r, w)` from the axis
    out and back to it (the first and last on the axis), as Fractions.
    `parts` are its rising segments' lines `(r0, k)` (`r = r0 + k w`), their
    heights and ranges; `terms` the max-min form's commons (indices into
    `parts`)."""

    def __init__(self, points):
        self.points = [(F(r), F(w)) for r, w in points]
        pts = self.points
        assert pts[0][0] == 0 and pts[-1][0] == 0, 'a profile from the axis to the axis'
        for (r0, w0), (r1, w1) in zip(pts, pts[1:]):
            assert w1 >= w0, 'heights never falling'
        self.w_min, self.w_max = pts[0][1], pts[-1][1]
        segs = list(zip(pts, pts[1:]))
        self.segs = segs
        rising = [i for i, (a, b) in enumerate(segs) if b[1] > a[1]]
        self.parts = []
        for i in rising:
            (ra, wa), (rb, wb) = segs[i]
            k = (rb-ra)/(wb-wa)
            r0 = ra-k*wa
            reach = [self.reach(i, 0), self.reach(i, 1)]
            self.parts.append({'seg': i, 'r0': r0, 'k': k, 'interval': (wa, wb), 'reach': reach})
        # Each rising segment's term over potential ranges (a corner's
        # anywhere).
        big = (self.w_max-self.w_min)+max(r for r, _ in pts)
        for p in self.parts:
            p['range'] = self.range_of(p, self.w_min-big, self.w_max+big)
        self.every = []
        for i, pi in enumerate(self.parts):
            a, b = pi['interval']
            term = [j for j, pj in enumerate(self.parts)
                    if j == i or (self.covers(pj, a, b) and self.at(pj, a) >= self.at(pi, a)
                                  and self.at(pj, b) >= self.at(pi, b))]
            self.every.append(term)
        # Each term then made as small as the form allows (a member dropped
        # where the form, its ranges taken again, is still the profile; in the
        # parts' order): fewer and shorter primitives for the reference to
        # sweep.
        self.settle()
        for i in range(len(self.every)):
            for j in sorted(self.every[i], reverse=True):
                if j == i:
                    continue
                kept = self.every[i]
                self.every[i] = [x for x in kept if x != j]
                self.settle()
                if not self.holds():
                    self.every[i] = kept
                    self.settle()
        self.check()

    def settle(self):
        """The terms (`every`'s, each once, a term holding another's dropped)
        and each part's range: past a corner the heights of the faces whose
        terms hold it and a quarter of its own beyond, on its radius's
        dyadic grid, moved on along it while within a sixteenth of a node's
        height, a cap's or an earlier part's end (two primitives' caps on one
        plane or nearly make events of the reference's sweeps a rounding
        apart)."""
        terms = []
        for t in self.every:
            if t not in terms:
                terms.append(t)
        self.terms = [t for t in terms if not any(u != t and set(u) <= set(t) for u in terms)]
        taken = {w for _, w in self.points}
        for p in self.parts:
            taken |= {h[1] for h in p['reach'] if h[0] == 'cap'}
        for j, p in enumerate(self.parts):
            lo, hi = p['interval']
            for i, t in enumerate(self.every):
                if j in t:
                    lo, hi = min(lo, self.parts[i]['interval'][0]), max(hi, self.parts[i]['interval'][1])
            own = p['interval'][1]-p['interval'][0]
            ends = [self.grid(p, lo-own/4, -1), self.grid(p, hi+own/4, 1)]
            step = F(1, 8)*(p['k'].denominator if p['k'] != 0 else 1)
            for end, sign in ((0, -1), (1, 1)):
                if p['reach'][end][0] in ('cap', 'apex'):
                    continue
                while any(abs(ends[end]-x) < F(1, 16) for x in taken):
                    ends[end] += sign*step
            p['range'] = self.range_of(p, ends[0], ends[1])
            taken |= set(p['range'])

    def holds(self):
        try:
            self.check()
        except AssertionError:
            return False
        return True

    @staticmethod
    def at(p, w):
        return p['r0']+p['k']*w

    def reach(self, i, end):
        """How a rising segment's end reaches: ('apex',), ('cap', h),
        ('near', h) or ('corner',)."""
        segs = self.segs
        (ra, wa), (rb, wb) = segs[i]
        node = segs[i][end]
        if node[0] == 0:
            return ('apex',)
        j = i-1 if end == 0 else i+1
        (sa, ha), (sb, hb) = segs[j]
        if ha != hb:
            return ('corner',)
        # A plane: this segment's end its wider side's (its cap), or the
        # narrower side's (past it into the wider side by a quarter of its
        # height).
        other = sb if end == 1 else sa
        if node[0] > other:
            return ('cap', ha)
        wider = segs[j-1] if end == 0 else segs[j+1]
        span = wider[1][1]-wider[0][1]
        return ('near', ha-span/4 if end == 0 else ha+span/4)

    def grid(self, p, w, sign):
        """The height at or beyond `w` (away from the segment by `sign`) on
        its line's dyadic grid (its radius there dyadic, eighths)."""
        k = p['k']
        step = F(1, 8)*(k.denominator if k != 0 else 1)
        base_w = p['interval'][0]
        n = (w-base_w)/step
        n = -((-n.numerator)//n.denominator) if sign > 0 else n.numerator//n.denominator
        return base_w+n*step

    def range_of(self, p, lo, hi):
        """A part's range: its ends by their reach, a cone clamped at its
        apex."""
        ends = []
        for end, want in ((0, lo), (1, hi)):
            how = p['reach'][end]
            if how[0] in ('cap', 'near'):
                ends.append(how[1])
            else:
                ends.append(want)
        if p['k'] != 0:
            apex = -p['r0']/p['k']
            if p['k'] > 0 and (ends[0] <= apex or p['reach'][0][0] == 'apex'):
                ends[0] = apex
            if p['k'] < 0 and (ends[1] >= apex or p['reach'][1][0] == 'apex'):
                ends[1] = apex
        return tuple(ends)

    def covers(self, p, a, b):
        lo, hi = p['range']
        return lo <= a and hi >= b

    def value(self, p, w):
        lo, hi = p['range']
        if w < lo or w > hi:
            return None
        r = self.at(p, w)
        return r if r > 0 else None

    def form(self, w):
        best = None
        for t in self.terms:
            vals = [self.value(self.parts[j], w) for j in t]
            if any(v is None for v in vals):
                continue
            m = min(vals)
            best = m if best is None else max(best, m)
        return best

    def profile(self, w):
        for (ra, wa), (rb, wb) in self.segs:
            if wa < w < wb:
                return ra+(rb-ra)*(w-wa)/(wb-wa)
        return None

    def check(self):
        """The form is the profile at every height (two points between every
        two breakpoints, and beyond them)."""
        brk = {w for _, w in self.points}
        for p in self.parts:
            brk |= set(p['range'])
            if p['k'] != 0:
                brk.add(-p['r0']/p['k'])
        for a, x in enumerate(self.parts):
            for y in self.parts[a+1:]:
                if x['k'] != y['k']:
                    brk.add((y['r0']-x['r0'])/(x['k']-y['k']))
        brk = sorted(brk)
        pts = [brk[0]-1, brk[-1]+1]
        for a, b in zip(brk, brk[1:]):
            pts += [a+(b-a)/3, a+2*(b-a)/3]
        for w in pts:
            assert self.form(w) == self.profile(w), ('the max-min form off the profile', w)

    def specs(self, frame_at):
        """Each part's reference solid: a cylinder (the prism of its circle)
        or a frustum on the axis frame from its lower height (`frame_at(w)`
        the axis frame at height `w`; `frame_at(None)` at the origin)."""
        out = []
        for p in self.parts:
            lo, hi = p['range']
            if p['k'] == 0:
                out.append(prism([forms.disc(0.0, 0.0, float(p['r0']))], frame_at(None), float(lo), float(hi)))
            else:
                r_lo, r_hi = self.at(p, lo), self.at(p, hi)
                for x in (r_lo, r_hi, hi-lo, lo):
                    assert float(x) == x, ('a frustum off the dyadic grid', x)
                out.append(cone(float(r_lo), float(r_hi), float(hi-lo), frame_at(lo)))
        return out

    def expr(self, k=0):
        """The max-min form over the specs from index `k`: the fuse of the
        terms' commons."""
        def chain(op, xs):
            e = xs[0]
            for x in xs[1:]:
                e = (op, e, x)
            return e
        return chain('fuse', [chain('common', [k+j for j in t]) for t in self.terms])

    def closed(self):
        """Pappus: `(volume, area)` exactly but for the slant heights."""
        vol, area = mp.mpf(0), mp.mpf(0)
        for (ra, wa), (rb, wb) in self.segs:
            if wb > wa:
                vol += q(F(1, 3)*(wb-wa)*(ra*ra+ra*rb+rb*rb))
                area += q(ra+rb)*mp.sqrt(q((rb-ra)**2+(wb-wa)**2))
            else:
                area += q(abs(rb*rb-ra*ra))
        return mp.pi*vol, mp.pi*area


# ------------------------------------------------------------------ bodies

class Body:
    """An imported turned body: its profile about an axis through `origin`
    along the world's `z` or `x`, OCCT's rows for it, its lattice, its closed
    form."""

    def __init__(self, name, points, axis, origin):
        self.name, self.axis, self.origin = name, axis, tuple(F(c) for c in origin)
        self.path = f'imported/profile_{name}.brep'
        self.lattice = Lattice(points)
        self.ref_specs = self.lattice.specs(self.frame_at)
        self.ref_expr = self.lattice.expr()

    def place(self, r, w):
        """The world point at radius `r` (along the frame's `x`) and height
        `w`."""
        o = self.origin
        if self.axis == 'z':
            return (o[0]+F(r), o[1], o[2]+F(w))
        return (o[0]+F(w), o[1], o[2]+F(r))

    def frame_at(self, w):
        p = self.place(0, 0 if w is None else w)
        name = 'XY' if self.axis == 'z' else 'SIDE'
        return at(name, p)

    def native_rows(self):
        o = self.origin
        pts = self.lattice.points
        wire = [f'wire S {len(pts)}']
        for r, w in pts:
            wire.append(' '.join(num(c) for c in self.place(r, w)))
            wire.append('L')
        direction = '0 0 1' if self.axis == 'z' else '1 0 0'
        return [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} {direction} {TURN}']

    def closed(self):
        return self.lattice.closed()

    def kinds(self):
        """Its faces' surface kinds: a rising segment's cylinder or cone, a
        level one's plane."""
        out = []
        for (ra, wa), (rb, wb) in self.lattice.segs:
            out.append('plane' if wa == wb else ('cylinder' if ra == rb else 'cone'))
        return sorted(out)


def make_bodies():
    point = [(0, 0), (3, 3), (3, 6), (4, F(17, 2)), (4, 16), (0, 16)]
    out = [
        Body('point', point, 'z', (0, 0, 0)),
        Body('point_x', point, 'x', (F(-3, 2), F(33, 8), F(3, 8))),
        Body('waist', [(0, 0), (8, 0), (8, 10), (5, 14), (5, 20), (0, 20)], 'z', (0, 0, 0)),
        Body('steps', [(0, 0), (6, 0), (6, 4), (5, 5), (3, 5), (3, 9), (4, 11), (4, 13), (0, 13)], 'z',
             (0, 0, 0)),
        Body('crown', [(0, 0), (3, 0), (3, 4), (6, 4), (6, 7), (4, 9), (4, 11), (0, 11)], 'z', (1, 2, 0)),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its
    lattice's, S9e.4a's imported body its construction) and each input's
    expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(b.lattice.expr(k))
        elif it[0] == 'base':
            exprs.append(len(specs))
            specs.append(base.BODIES[it[1]].spec)
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
            if s[0] in ('imported', 'base'):
                path = (BODIES if s[0] == 'imported' else base.BODIES)[s[1]].path
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=path)
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


# The point's box across its first corner (radius 3 at height 3), the corner's
# point (3, 0, 3) at the box's chart point (1, 1/5, 53/20) inside it.
POINT_BOX = block(at('TILT', (2, -1.75, 1)), 0.0, 0.0, 3.0, 3.0, 0.0, 3.0)
# The waist's rod along the world's x about (y, z) = (3/4, 12), through its cone.
WAIST_ROD = cylinder(at('SIDE', (-12, 0.75, 0)), 0.0, 12.0, 1.5, 0.0, 24.0)
STEPS_BALL = sphere(1.75, at('XY', (3.5, 1.0, 5.25)))
# A box across the crown's step out (height 4, radii 3 to 6).
BOX_CROWN = block(at('TILT', (3.5, 1.25, 2.5)), 0.0, 0.0, 4.0, 3.0, 0.0, 2.5)
# The chain's box biting the waist's convex corner (radius 8 at height 10)
# from the side, then a `TILT` slab across its neck and cone.
CHAIN_BOX = block(at('XY', (0, 0, 0)), 5.75, -2.5, 10.0, 2.5, 7.5, 12.5)
CHAIN_SLAB = block(at('TILT', (0, 0, 9.5)), -9.0, -9.0, 9.0, 9.0, 0.0, 1.25)
WAIST_TOUCH = sphere(2.0, at('XY', (1.5, -1.0, 22.0)))
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('point_tilt', 'both', [imported('point_x'), ('base', 'box_tilt')], three)
    out += group('chain_waist', 'chain', [imported('waist'), CHAIN_BOX, CHAIN_SLAB], three, first='cut')
    out += group('point_box', 'pieces', [imported('point'), POINT_BOX], three)
    out += group('waist_rod', 'pieces', [imported('waist'), WAIST_ROD], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('steps_ball', 'pieces', [imported('steps'), STEPS_BALL], three)
    out += group('box_crown', 'pieces', [BOX_CROWN, imported('crown')], three)
    out += group('waist_touch', 'pieces', [imported('waist'), WAIST_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
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
        d = max(abs(v-V)/size**3, abs(a-A)/size**2)
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
    cache = os.environ.get('TURNED_PROFILES_CACHE')
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
# Groups whose events' spacing is not checked (their sines and gaps are):
# a body's corner circle (two of its reference solids' surfaces meeting)
# through the partner box's walls, where the reference finds each point as
# both surfaces' event, a rounding apart.
SPACING_UNCHECKED = {'point_box', 'box_crown'}


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
    """Every case's native rows (`compare_turned_profiles_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def profile_distance(body, X):
    """A world point's distance from the body's profile's surfaces (each
    segment's cone, cylinder or plane, unbounded), in its axis frame."""
    o = [q(c) for c in body.origin]
    d = [X[i]-o[i] for i in range(3)]
    if body.axis == 'z':
        r, w = mp.sqrt(d[0]**2+d[1]**2), d[2]
    else:
        r, w = mp.sqrt(d[1]**2+d[2]**2), d[0]
    best = None
    for (ra, wa), (rb, wb) in body.lattice.segs:
        if wa == wb:
            dist = abs(w-q(wa))
        else:
            # The line through the two points, `(rb - ra) (w - wa) - (wb -
            # wa) (r - ra) = 0`.
            a, b = q(rb-ra), -q(wb-wa)
            dist = abs(a*(w-q(wa))+b*(r-q(ra)))/mp.sqrt(a*a+b*b)
        best = dist if best is None else min(best, dist)
    return best


def check_bodies():
    """Each existing body file: its faces' kinds its profile's, each stored
    vertex on the profile's surfaces within 1e-12 of the body's size.
    Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        assert sorted(kinds) == b.kinds(), (b.name, kinds)
        size = max([1.0]+[abs(x) for v in vertices for x in v])
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            worst = max(worst, profile_distance(b, X)/size)
    assert worst <= 1e-12, ('a stored vertex off its profile', mp.nstr(worst, 3))
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
        os.environ['TURNED_PROFILES_CACHE'] = args.cache
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
          'size of its profile\'s surfaces')


if __name__ == '__main__':
    main()
