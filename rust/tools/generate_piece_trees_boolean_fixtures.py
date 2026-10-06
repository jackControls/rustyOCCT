#!/usr/bin/env python3
"""Fixtures for S9e.4b.3c.3b of REVIEW_NOTES.md: imported bodies of one
sphere or cylinder face and plane faces that are a Boolean tree of their
primitive and several convex hulls of their planes (bodies OCCT wrote to
`.brep` files) given to Booleans.

The bodies are OCCT's own output, as S9e.4b.3c.3a's
(`generate_piece_forms_boolean_fixtures.py`): `boolean-piece-trees-
bodies.txt` holds one `write NAME imported/form_NAME.brep` block per body
for `occt_boolean_oracle.cpp` (one solid's rows, a `boolean` row and
another solid's rows: the Boolean's one solid written), which
`compare_piece_trees_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/form_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-piece-trees-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/form_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-piece-trees-expected.tsv` gives
per case `expect KIND S9e.4b.3c.3b CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a later step's;
the class `pieces`: an imported body and a construction, `both`: two
imported inputs, `chain`: an imported input's result given to another
Boolean), a degenerate or unsupported case's `reason TEXT`, then `result N
volume area cx cy cz` (totals over the N solids, world coordinates) or
`empty`; `boolean-piece-trees-frames.tsv` the stored axes of every solid the
kernel builds from its rows.

Each body is one Boolean OCCT was given, of a primitive (a sphere's
hemisphere or whole ball, or a cylinder: a prism of a disc or of a profile
of lines and an arc) and a prism of a profile of lines (a box, or a U: a
body not convex in its planes), whose sections by each other's planes are
circles and lines (the `.brep` reader takes no ellipse, and a sphere's
section off its stored frame's parallels and meridians has a pcurve the
converter does not certify, S9e.4b.3a's evidence), written to
`imported/form_NAME.brep`; each needs more than one of S9e.4b.3c.3a's hulls
(`why`): `u_scoop`, a U prism on `SKEW4` less a ball of radius 5/2 centred
on the top face of one arm (the ball's frame normal along the U's `y`, its
`x` the U's normal: the top face a meridian plane, the arm's walls
parallels' planes, the slot's floor, the ends and the bottom out of reach),
the groove crossing the arm into the slot (the fuzz target's ball's groove
in a U prism: the hull of the U's outer planes less its slot, less the
ball); `u_boss`, a U prism on the world's axes fused with a cylindrical
boss of radius 3/2 on its base's top face from that face to a cap (G4's
boss on a body not convex); `bites`, a cylinder of radius 3 on the world's
axes less a U prism along its axis between its caps whose two arms bite
its wall from outside (a primitive bitten twice, by a body not convex;
the seam between the bites); `cap_boss`, a box on `SKEW` fused with an
upper hemisphere of radius 3 whose disc lies 1 below the box (the ball
meeting the box's bottom alone, a parallel's plane: the primitive common
its disc's half-space fused with the box, the fuzz target's ball's half
with a box); `cap_pocket`, a box on `SKEW2` less an upper hemisphere of
radius 3 whose disc lies 2 inside it (a dish with a flat floor: the box
less the primitive common its disc's half-space); `dee_boss`, a box on the
world's axes fused with a prism along `x` of a disc of radius 5/2 less a
segment, its flat along the axis (a plane holding the cylinder's axis
direction, through no seam: its frame's `x` beyond the flat), the box's
top face through its axis (a cylinder's end not normal to its axis);
declared `degenerate`: `u_notch`, a U prism less a cylinder of radius 3
along `z` about the middle of its base's end face, tangent to the two
faces it meets (`bcut_complex/I6`'s notch on a body not convex; its frame's
`x` reversed: its seam off the body); declared `unsupported`: `tooth`, a
cylinder less a U prism whose base and slot lie inside it, a tooth left
between the arms (a pocket within a pocket, S9e.4b.4's).

Cases (each the three operations): `u_scoop_rod`, the U's groove and an
upright rod through it; `u_boss_slab`, the U's boss and a `TILT` slab
across the boss and the slot; `bites_box`, the bitten cylinder and a
`TILT` box across one bite; `cap_boss_rod`, the hemisphere's boss and a rod
along `x` through its disc; `rod_cap_pocket`, a `TILT` rod less the dish
(the body the tool); `dee_boss_ball`, the flattened boss and a ball across
its flat; `pieces`, the hemisphere's boss and the flattened boss, both
imported; `chain_dee`, the flattened boss less a rod, then with a `TILT`
slab. Declared `degenerate`: `u_notch_rod`, the U's notch (its wall
tangent to its own faces) and a rod; declared `unsupported`: `tooth_rod`,
the tooth and a rod. (The partners are rods where they can be: the
reference's cost grows with the faces of a case.)

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`), which models convex
profiles and whole spheres: a U prism as the box of its outer planes less
its slot's box (`Slotted`), a hemisphere as its ball common a box above its
equator's plane (`Hemisphere`), the same sets; each imported body its
first Boolean, `(X op Y) op2 C` (swapped, the body the tool), two bodies
`(X op Y) op2 (Z op W)`, a chain `((X op Y) op2 C) op3 D`, with S9e.4b.3c.3a's
checks (the two families, each solid's closed form, the pair identities on
the last Boolean's arguments, Monte Carlo, solid counts by rays at two
resolutions, each body one solid, the meetings' sines and the events'
spacing outside the declared groups) and each body's closed form where it
has one (`closed`). With `--check`, where the bodies' files exist, each is
read (`stored_records`): its faces' kinds one sphere or cylinder and
planes, every stored vertex within 1e-12 of the size on the construction's
surfaces.
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
from generate_curved_boolean_fixtures import square, disc

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-piece-trees'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.3c.3b'
FRAMES = forms.FRAMES
HP = base.HP
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
at, q, point, moved, reversed_x = forms.at, forms.q, forms.point, forms.moved, forms.reversed_x
cylinder, block = forms.cylinder, forms.block
segment_area, cap_volume, box_volume, plane_distance = (forms.segment_area, forms.cap_volume, forms.box_volume,
                                                        forms.plane_distance)


# ------------------------------------------------------------------ bodies

class Body(forms.Body):
    """An imported body `first op second`, OCCT's rows for it, which of its
    two solids is the curved primitive (`curved`: 0 or 1), why no single
    form of S9e.4b.3c.3a's is it (`why`), its closed form `(volume, area)`
    where it has one, and its declared refusal (`kind` `degenerate` or
    `unsupported`, its reason `reason`). The reference takes a U prism (a
    profile not convex, which it does not model) as the box of its hull
    less its slot's box (`Slotted`), and an upper hemisphere (it models
    whole spheres) as its ball common a box above its equator's plane
    (`Hemisphere`): the same sets, its solids `ref_specs` and its expression
    over them `ref_expr`."""

    def __init__(self, name, first, op, second, curved, why, closed=None, kind=None, reason=None):
        super().__init__(name, first, op, second, None, closed, reason)
        self.curved, self.why, self.kind = curved, why, kind
        specs, exprs = [], []
        for s in (first, second):
            if isinstance(s, (Slotted, Hemisphere)):
                k = len(specs)
                specs += s.parts
                exprs.append((s.op, k, k+1))
            else:
                exprs.append(len(specs))
                specs.append(s)
        self.ref_specs, self.ref_expr = specs, (op, exprs[0], exprs[1])

    def native_rows(self):
        rows = lambda spec: native_case(construction(spec, 0)).split('\n')[1:-1]
        first, second = (s.native if isinstance(s, (Slotted, Hemisphere)) else s
                         for s in (self.first, self.second))
        return rows(first)+[f'boolean {self.op}']+rows(second)

    def primitive(self):
        p = (self.first, self.second)[self.curved]
        return p.native if isinstance(p, Hemisphere) else p


class Slotted:
    """A U prism for OCCT (`prism`: a profile of lines) and the reference's
    same set, a box less a slot's box reaching past the U's open end and
    past both its caps."""

    def __init__(self, frame, h0, h1, rect, slot):
        x0, y0, x1, y1 = rect
        s0, s1, t0, t1 = slot
        assert s1 == x1 and x0 < s0 < x1 and y0 < t0 < t1 < y1
        self.rect, self.slot_rect = rect, slot
        self.prism = prism([u_profile(x0, y0, x1, y1, slot)], frame, h0, h1)
        self.box = block(frame, x0, y0, x1, y1, h0, h1)
        self.slot = block(frame, s0, t0, x1+1.0, t1, h0-1.0, h1+1.0)
        self.native, self.parts, self.op = self.prism, [self.box, self.slot], 'cut'


class Hemisphere:
    """An upper hemisphere for OCCT (`native`: a sphere's row from its
    equator to its pole) and the reference's same set, its whole ball common
    a cylinder on its frame from its equator's plane, of twice its radius,
    four radii high (its wall and top apart from the ball, its top above the
    body's box; three faces where a box has six, the reference's cost)."""

    def __init__(self, r, frame):
        self.native = sphere(r, frame, 0.0, HP)
        self.ball = sphere(r, frame)
        self.half = cylinder(frame, 0.0, 0.0, 2.0*r, 0.0, 4.0*r)
        self.parts, self.op = [self.ball, self.half], 'common'


def u_profile(x0, y0, x1, y1, slot):
    """A U: the rectangle `[x0, x1] x [y0, y1]` less the slot `(s0, s1, t0,
    t1)`, `[s0, x1] x [t0, t1]` open at `x = x1`, counter-clockwise."""
    s0, _, t0, t1 = slot
    return Boundary(points=[(x0, y0), (x1, y0), (x1, t0), (s0, t0), (s0, t1), (x1, t1), (x1, y1), (x0, y1)])


def dee(cx, cy, r, flat):
    """The disc of radius `r` about `(cx, cy)` on the side `x <= flat` of a
    chord at distance 3 r / 5 (its ends rational: a 3-4-5 triangle)."""
    h = F(4, 5)*F(r)
    assert F(flat)-F(cx) == F(3, 5)*F(r)
    return Boundary(points=[(flat, float(F(cy)-h)), (flat, float(F(cy)+h))],
                    segments=[None, (cx, cy, r, True)])


# The U's groove: `forms`'s scoop slab, a U, its ball on one arm's top face.
U_SCOOP_FRAME = at('SKEW4', (1, 2, 0.5))
U_SCOOP_CENTRE = point(U_SCOOP_FRAME, 8, 1.5, 6)
U_SCOOP_BALL = U_SCOOP_CENTRE+(4.0, -8.0, -1.0, 4.0, 1.0, 8.0)
U_BOSS_FRAME = at('XY', (0, 0, 0))
BITES_FRAME = at('XY', (0, 0, 0))
CAP_BOSS_FRAME = at('SKEW', (1, 2, 0.5))
CAP_BOSS_CENTRE = point(CAP_BOSS_FRAME, 4, 4, -1)
CAP_POCKET_FRAME = at('SKEW2', (0.5, 1, 1.5))
CAP_POCKET_CENTRE = point(CAP_POCKET_FRAME, 4, 4, 2)
# The flattened boss along the world's `x`: a `SIDE` frame's local `(u, v,
# w)` the world's `(y, z, x)`; its circle's `x` the frame's, beyond the flat.
DEE_FRAME = at('SIDE', (0, 0, 0))
U_NOTCH_AXIS = at('XY', (0, 0, 0))[:6]+(-1.0, 0.0, 0.0)


def u_area(x0, y0, x1, y1, slot):
    s0, s1, t0, t1 = slot
    return (x1-x0)*(y1-y0)-(s1-s0)*(t1-t0)


def make_bodies():
    out = [
        # The groove across one arm into the slot: the half ball between
        # the arm's walls `y = 0` and `y = 3` (a slab of the ball about its
        # centre, 3/2 either side).
        Body('u_scoop', Slotted(U_SCOOP_FRAME, 0.0, 6.0, (0.0, 0.0, 12.0, 8.0), (4.0, 12.0, 3.0, 5.0)), 'cut',
             sphere(2.5, U_SCOOP_BALL), 1,
             'the hull of the U\'s outer planes less its slot (a pocket), less the ball'),
        Body('u_boss', Slotted(U_BOSS_FRAME, 0.0, 4.0, (0.0, 0.0, 11.0, 8.0), (5.0, 11.0, 3.0, 5.0)), 'fuse',
             cylinder(U_BOSS_FRAME, 2.5, 4.0, 1.5, 4.0, 7.0), 1,
             'the capped cylinder fused with the hull of the U\'s outer planes less its slot (a pocket)',
             closed=lambda: (u_area(0, 0, 11, 8, (5, 11, 3, 5))*4+mp.pi*q(1.5)**2*3, None)),
        Body('bites', cylinder(BITES_FRAME, 4.0, 4.0, 3.0, 0.0, 6.0), 'cut',
             Slotted(BITES_FRAME, 2.0, 4.5, (-3.0, 1.5, 3.5, 6.5), (-1.0, 3.5, 2.5, 5.5)), 0,
             'the capped cylinder less two pockets, one for each arm'),
        Body('cap_boss', block(CAP_BOSS_FRAME, 0.0, 0.0, 8.0, 8.0, 0.0, 6.0), 'fuse',
             Hemisphere(3.0, moved(CAP_BOSS_FRAME, CAP_BOSS_CENTRE)), 1,
             'the ball common its disc\'s half-space, fused with the box',
             closed=lambda: (box_volume(CAP_BOSS_FRAME, 384)+2*mp.pi*27/3
                             - cap_volume(3, -plane_distance(CAP_BOSS_FRAME, 0, CAP_BOSS_CENTRE)), None)),
        Body('cap_pocket', block(CAP_POCKET_FRAME, 0.0, 0.0, 8.0, 8.0, 0.0, 4.0), 'cut',
             Hemisphere(3.0, moved(CAP_POCKET_FRAME, CAP_POCKET_CENTRE)), 1,
             'the box less the ball common its disc\'s half-space',
             closed=lambda: (box_volume(CAP_POCKET_FRAME, 256)-2*mp.pi*27/3
                             + cap_volume(3, -plane_distance(CAP_POCKET_FRAME, 4, CAP_POCKET_CENTRE)), None)),
        # The flat at `y = 21/8` (local `u`), the box's top face `z = 3`
        # through the axis: the half disc above it less half the segment.
        Body('dee_boss', block(BITES_FRAME, 3.25, -2.875, 13.25, 4.625, -1.0, 3.0), 'fuse',
             prism([dee(1.125, 3.0, 2.5, 2.625)], DEE_FRAME, 5.625, 11.625), 1,
             'the capped cylinder common its flat\'s half-space, fused with the box',
             closed=lambda: (10*q(7.5)*4+6*(mp.pi*q(2.5)**2/2-segment_area(2.5, 1.5)/2), None)),
        # I6's notch on a U: the wall tangent to `y = 0` and `y = 6`.
        Body('u_notch', Slotted(at('XY', (0, 0, 0)), 0.0, 3.0, (0.0, 0.0, 8.0, 6.0), (5.0, 8.0, 2.5, 3.5)), 'cut',
             cylinder(U_NOTCH_AXIS, 0.0, -3.0, 3.0, -1.0, 4.0), 1,
             'the hull of the U\'s outer planes less its slot (a pocket), less the cylinder',
             kind='degenerate', reason='an imported plane piece whose curved face is tangent to its plane faces'),
        Body('tooth', cylinder(BITES_FRAME, 4.0, 4.0, 3.0, 0.0, 6.0), 'cut',
             Slotted(BITES_FRAME, 2.0, 4.5, (-3.0, 2.5, 3.0, 5.5), (2.0, 3.0, 3.5, 4.5)), 0,
             'the capped cylinder less a U-shaped pocket holding a tooth (a pocket within a pocket)',
             kind='unsupported',
             reason='an imported plane piece other than a Boolean tree of its primitive and its planes\' hulls '
                    '(S9e.4b.4)'),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def shifted(expr, k):
    """An expression over a body's own solids, its leaves moved by `k`."""
    if isinstance(expr, int):
        return expr+k
    return (expr[0], shifted(expr[1], k), shifted(expr[2], k))


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its two
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(shifted(b.ref_expr, k))
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


def group(name, klass, items, outcomes, first=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); a declared body's refusal declares its cases."""
    out = []
    declared = [BODIES[s[1]] for s in items if s[0] == 'imported' and BODIES[s[1]].kind]
    for op, n in outcomes.items():
        if declared:
            kind, reason = declared[0].kind, declared[0].reason
        else:
            kind, reason = ('solid' if n else 'empty'), None
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, kind, n, reason))
    return out


U_SCOOP_ROD = cylinder(at('XY', point(U_SCOOP_FRAME, 7.25, 2.25, 0)), 0.0, 0.0, 0.625, -4.0, 12.0)
U_BOSS_SLAB = block(at('TILT', (0, 0, 0)), -1.0, -4.0, 12.0, 8.0, 4.0, 6.0)
BITES_BOX = block(at('TILT', (0, 0, 0)), 1.3125, 1.4375, 6.5625, 4.5625, 5.0625, 6.9375)
CAP_BOSS_ROD = cylinder(at('SIDE', point(CAP_BOSS_FRAME, 4.25, 3.75, -0.5)), 0.0, 0.0, 0.75, -6.0, 6.0)
CAP_POCKET_ROD = cylinder(at('TILT', point(CAP_POCKET_FRAME, 4.5, 3.75, 2.5)), 0.0, 0.0, 1.0, -6.0, 6.0)
ROD_CAP_POCKET = {'fuse': 1, 'cut': 2, 'common': 1}
DEE_BALL = sphere(1.75, at('XY', (7.75, 2.75, 4.25)))
CHAIN_ROD = cylinder(BITES_FRAME, 8.125, 1.875, 0.625, -2.0, 8.0)
CHAIN_SLAB = block(at('TILT', (0, 0, 0)), 2.0, -6.0, 15.0, 8.0, 4.75, 6.5)
U_NOTCH_ROD = cylinder(at('XY', (0, 0, 0)), 3.25, 1.25, 0.5, -1.0, 4.0)
TOOTH_ROD = cylinder(BITES_FRAME, 5.0, 4.0, 0.75, -1.0, 7.0)


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('u_scoop_rod', 'pieces', [imported('u_scoop'), U_SCOOP_ROD], three)
    out += group('pieces', 'both', [imported('cap_boss'), imported('dee_boss')], three)
    out += group('bites_box', 'pieces', [imported('bites'), BITES_BOX], three)
    out += group('rod_cap_pocket', 'pieces', [CAP_POCKET_ROD, imported('cap_pocket')], ROD_CAP_POCKET)
    out += group('u_boss_slab', 'pieces', [imported('u_boss'), U_BOSS_SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('cap_boss_rod', 'pieces', [imported('cap_boss'), CAP_BOSS_ROD], three)
    out += group('chain_dee', 'chain', [imported('dee_boss'), CHAIN_ROD, CHAIN_SLAB], {'fuse': 1, 'cut': 1, 'common': 2},
                 first='cut')
    out += group('dee_boss_ball', 'pieces', [imported('dee_boss'), DEE_BALL], three)
    out += group('tooth_rod', 'pieces', [imported('tooth'), TOOTH_ROD], three)
    out += group('u_notch_rod', 'pieces', [imported('u_notch'), U_NOTCH_ROD], three)
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
        for s in c.items:
            if s[0] == 'imported' and BODIES[s[1]].kind:
                assert c.kind == BODIES[s[1]].kind and c.reason == BODIES[s[1]].reason, c.name
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
        # two solids') leave the area identity.
        owner = pieces.owners(items, sub)
        shared = any(owner[i] != owner[sw.others.surfs[k][0]] for i, _, sw in chain.sweeps(0) for k in sw.skip)
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    for op, (expr, solids) in exprs.items():
        for n in (48, 71):
            got = forms.count_solids(chain, expr, n)
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
    cache = os.environ.get('PIECE_TREES_CACHE')
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
# the dish's reference half box has its base, the hemisphere's equator's
# plane through the ball's rounded centre, within rounding of a whole family
# curve of the dish's box's walls (lines at that height), so the reference
# finds one event there twice, rounding apart.
SPACING_UNCHECKED = {'rod_cap_pocket'}


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # A declared group's reference meets its body's own tangency or its
        # pocket within a pocket: its checks kept apart, looser.
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
    """Every case's native rows (`compare_piece_trees_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def reflex(boundary):
    """Whether a polygon has a reflex corner (a body not convex in its
    planes)."""
    p = boundary.points
    n = len(p)
    turns = [(p[(i+1) % n][0]-p[i][0])*(p[(i+2) % n][1]-p[(i+1) % n][1])
             - (p[(i+1) % n][1]-p[i][1])*(p[(i+2) % n][0]-p[(i+1) % n][0]) for i in range(n)]
    return any(t < 0 for t in turns) and any(t > 0 for t in turns)


def check_bodies():
    """Each existing body file: its faces' kinds one sphere or cylinder and
    planes, each stored vertex on the construction's surfaces within 1e-12
    of the body's size. Returns the largest deviation and the files read."""
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
        os.environ['PIECE_TREES_CACHE'] = args.cache
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
