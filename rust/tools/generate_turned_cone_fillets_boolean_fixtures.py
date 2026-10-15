#!/usr/bin/env python3
"""Fixtures for S9e.4b.4d.2c.1b of REVIEW_NOTES.md: imported turned bodies
whose fillets and rounds (tori) join a cone smoothly (the kind of
`OCC496b`'s third torus, between a cylinder and a cone), each torus snapped
tangent to its two lines, given to Booleans whose partners cross no smooth
edge.

The bodies are OCCT's own output, as S9e.4b.4d.2c.1a's
(`generate_turned_fillets_boolean_fixtures.py`):
`boolean-turned-cone-fillets-bodies.txt` holds one `write NAME
imported/cone_fillet_NAME.brep` block per body for `occt_boolean_oracle.cpp`
(a `plane` row, the profile's `wire` of lines and arcs, a `revolve` row
turning it a whole turn about the world's `z`, and for the refitted body a
`tolerance` row setting every tolerance to 1e-3), which
`compare_turned_cone_fillets_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/cone_fillet_NAME.brep`. Nothing here reads the
kernel; the files are read only to check that each is the body it claims to
be (below). The cases' files are S9e.4b.4d.2c.1a's, under the prefix
`boolean-turned-cone-fillets`, the class `pieces`.

Each body has two profiles `(r, w)` (`r` the radius, `w` the height along
the axis), from the axis out and back to it, each segment a line (a plane
where level, a cylinder where upright, else a cone) or an arc (a torus) to
the next point: the one OCCT is given (`occt`, binary64 points: the model's
rounded once where it is no binary64 number) and the model's (`model`,
exact): each arc a fillet (concave, its material outside the tube) or a
round (convex) between a cone and a plane or a cylinder, tangent to both,
the direction from its centre to its point on the cone `((1 - s^2), 2 s) /
(1 + s^2)` for a dyadic half-angle tangent `s` (the kernel's snap keeps it:
REVIEW_NOTES' "S9e.4b.4d.2c.1b refined"), the cone its tangent there. The
bodies, each about the world's `z` through the origin: `round` (`OCC496b`'s
third torus's kind): a disc of radius 6, a cylinder to 4, a round of radius 2
about `(4, 4)` turning by `s = 1/4` (the direction `(15/17, 8/17)`) to the
cone narrowing by 8/15 to the top's disc at 8; `fillet`: a disc of radius
8, a cylinder to 2, a step in to radius 5 (an annulus facing up), a fillet of
radius 1 about `(5, 3)` turning to `s = -4` (the direction `(-15/17,
-8/17)`) to a boss's cone narrowing by 8/15 to its disc at 7; `top`: a disc
of radius 3, a cone widening by 4/3 to a round of radius 1 about `(8, 5)`
from `s = -1/2` (the direction `(3/5, -4/5)`) to the top's disc at 6;
`refit` (`OCC496b`'s kind at a resolution of 1e-3): `round`'s profile whose
round OCCT is given with its centre raised by `D` = 2^-20 above the
cylinder's top at 4, its radius the distance to that point (rounded once),
and its cone through `(17/4, 8)` tangent to that circle (the tangent point
rounded once; no Pythagorean direction): the model's round is tangent to the
cylinder at `4 + D`, its direction's half-angle tangent OCCT's rounded to
the multiples of 2^-40 (the kernel's snap), its cone that direction's
tangent through the round's point on it.

The kernel's construction (`Model.construction_inside`, checked here against
the profile at points of the meridian, not swept): the sharp profile, each
arc replaced by the corner its two lines meet at, is S9e.4b.4d.2a's lattice
(`profiles.Lattice`); each fillet's piece fused to it and each round's taken
from it, a piece its box less the whole torus: the parallelogram in the
meridian from the tube's centre along the two directions to its tangent
points (beside a cone a coaxial cone through the centre normal to it, beside
a plane a cylinder, beside a cylinder a plane) past the corner by half the
tube (halved up to four times while another segment of the sharp profile
would meet it), its far plane or cylinder and its far cone's point on the
near one on the multiples of 2^-24 of the tube's power of two, outward.

Cases (each the three operations): `round_rod`, the round body and a rod
along the world's `x` across its round alone (between its tangent circles);
`fillet_ball`, the fillet body and a ball across its fillet alone;
`top_box`, the top body and a skew box across its bottom's sharp corner;
`refit_box`, the refitted body and a `TILT` box across its top's sharp
corner; and declared `degenerate`, `round_touch`, the round body and a ball
touching its round at a point (`Degenerate("a torus tangent to the other
input's surface (S9d.4b.2)")`).

The reference is S9e.4b.4d.2c.1a's: S9e.3a's chained reference with each
body one input turned from its model profile (`TurnedInput`, its faces the
profile's own: each plane's annulus or disc, each cylinder's and each
cone's band, each arc's band of its whole torus's tube, each covered by two
families of the reference's circles and lines; its membership the profile
with arcs in the meridian), with S9e.4b.3c.3a's checks and each body's
closed form by Pappus over its model profile (its lines' frusta and annuli,
each arc's by Green's theorem in the meridian, exactly but for pi, the
slant heights' square roots and the arcs' angles). With `--check`, where the
bodies' files exist, each is read (`stored_records`): its faces' kinds its
OCCT profile's, every stored vertex within 1e-12 of the size on the OCCT
profile's lines and circles.
"""
import argparse
from fractions import Fraction as F
import math
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case
from curve_surface_reference import stored_axes
import generate_turned_fillets_boolean_fixtures as fillets

base = fillets.base
pieces = fillets.pieces
forms = fillets.forms
profiles = fillets.profiles
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS = base.OPERATIONS
CLASSES = ('pieces',)
PREFIX = 'boolean-turned-cone-fillets'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2c.1b'
sphere, construction, torus = base.sphere, base.construction, base.torus
ref = base.ref
q = forms.q
cylinder, block = forms.cylinder, forms.block
Lattice, at, TURN = profiles.Lattice, profiles.at, profiles.TURN
L, arc, points, meet, segment_meets = fillets.L, fillets.arc, fillets.points, fillets.meet, fillets.segment_meets
# The refitted round's offset.
D = F(1, 2**20)


# ------------------------------------------------------------------ profiles

def num(x):
    """A coordinate's text: a Fraction rounded once to binary64 (the
    shortest round-trip text), or a binary64 number's text as given."""
    if isinstance(x, str):
        return x
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else repr(float(x))


def kind_of(a, b):
    """A meridian segment's surface."""
    if a[1] == b[1]:
        return 'plane'
    return 'cylinder' if a[0] == b[0] else 'cone'


def foot(c, a, b):
    """The foot of `c` on the line through `a` and `b` (exact)."""
    d = (b[0]-a[0], b[1]-a[1])
    s = ((c[0]-a[0])*d[0]+(c[1]-a[1])*d[1])/(d[0]*d[0]+d[1]*d[1])
    return (a[0]+s*d[0], a[1]+s*d[1])


def solve(ea, eb, v):
    """`(x, y)` with `v = x ea + y eb` (exact)."""
    det = ea[0]*eb[1]-ea[1]*eb[0]
    assert det != 0, 'parallel directions'
    return ((v[0]*eb[1]-v[1]*eb[0])/det, (ea[0]*v[1]-ea[1]*v[0])/det)


def direction(s):
    """The rational unit direction of half-angle tangent `s`."""
    s = F(s)
    return ((1-s*s)/(1+s*s), 2*s/(1+s*s))


def half_tangent(e):
    """The half-angle tangent of a unit direction `(c, s)`."""
    c, s = e
    return s/(1+c) if c >= 0 else (1-c)/s


def power_above(x):
    """The least power of two at or above `x > 0` (exact)."""
    x = F(x)
    p = F(1)
    while p < x:
        p *= 2
    while p/2 >= x:
        p /= 2
    return p


def outward(x, g, sign):
    """`x` on the multiples of `g`, moved away by `sign` (exact)."""
    n = x/g
    n = -((-n.numerator)//n.denominator) if sign > 0 else n.numerator//n.denominator
    return n*g


def in_convex(poly, p):
    """Whether `p` lies strictly inside the convex polygon `poly` (either
    turn), in any ordered field."""
    sides = []
    n = len(poly)
    for k in range(n):
        (x0, y0), (x1, y1) = poly[k], poly[(k+1) % n]
        sides.append((x1-x0)*(p[1]-y0)-(y1-y0)*(p[0]-x0))
    return all(s > 0 for s in sides) or all(s < 0 for s in sides)


class Model:
    """A model profile with arcs: its sharp profile's lattice, each arc's
    corner (`centre`, `t`, `round`, `corner`, its tangent points `ends`, the
    unit directions `dirs` from its centre to them, the neighbours' `kinds`)
    and its piece's box (`box`: the parallelogram `C + x e_c + y e_f`, `0 <=
    x <= X`, `0 <= y <= Y`, `e_c` the direction toward the cone)."""

    def __init__(self, profile):
        self.profile = profile
        pts = points(profile)
        segs = [s for _, _, s in profile]
        self.pts, self.segs = pts, segs
        assert pts[0][0] == 0 and pts[-1][0] == 0, 'a profile from the axis to the axis'
        sharp, corners = [], []
        i = 0
        while i < len(pts):
            if segs[i][0] == 'A':
                _, cr, cw, t, ccw = segs[i]
                a0, a1, b0, b1 = pts[i-1], pts[i], pts[i+1], pts[i+2]
                kinds = (kind_of(a0, a1), kind_of(b0, b1))
                # A cone and a plane or a cylinder (S9e.4b.4d.2c.1a's plane
                # and cylinder, two cones no turned body's).
                assert sorted(kinds) in (['cone', 'cylinder'], ['cone', 'plane']), ('a torus not beside a cone', kinds)
                corner = meet(a0, a1, b0, b1)
                dirs = []
                for p, s, end in ((a0, a1, a1), (b0, b1, b0)):
                    assert foot((cr, cw), p, s) == end, ('an arc not tangent to its line', end)
                    assert (end[0]-cr)**2+(end[1]-cw)**2 == t*t, 'an end off its arc'
                    dirs.append(((end[0]-cr)/t, (end[1]-cw)/t))
                cone_side = kinds.index('cone')
                corners.append({'centre': (cr, cw), 't': t, 'round': ccw, 'corner': corner, 'ends': (a1, b0),
                                'dirs': dirs, 'kinds': kinds, 'cone': cone_side})
                sharp.append(corner)
                i += 2
            else:
                sharp.append(pts[i])
                i += 1
        self.sharp, self.corners = sharp, corners
        self.lattice = Lattice(sharp)
        for c in corners:
            c['box'] = self.box(c)

    @staticmethod
    def oblique(c, p):
        """A meridian point's coordinates `(x, y)` on the corner's box's
        axes (along `e_c` and `e_f` from the centre)."""
        e_c, e_f = c['dirs'][c['cone']], c['dirs'][1-c['cone']]
        return solve(e_c, e_f, (p[0]-c['centre'][0], p[1]-c['centre'][1]))

    def box(self, c):
        """The piece's box `(X, Y)`: past the corner by half the tube (halved
        up to four times), its far plane or cylinder and its far cone's point
        on the near plane or cylinder on the grid, outward; clear of the sharp
        profile's other segments, the lines at the corner leaving it."""
        C, t = c['centre'], c['t']
        e_c, e_f = c['dirs'][c['cone']], c['dirs'][1-c['cone']]
        xP, yP = self.oblique(c, c['corner'])
        assert xP > 0 and yP > 0, 'a corner off its box'
        g = power_above(t)/2**24
        # The far side along `e_f` crosses the axis-aligned coordinate `k`
        # of `e_c` (a plane's height beside a cylinder, a cylinder's radius
        # beside a plane); the far cone's point on the near side the other.
        k = 1 if e_f[1] == 0 else 0
        d = t/2
        for _ in range(5):
            far = C[k]+(xP+d)*e_c[k]
            X = (outward(far, g, 1 if e_c[k] > 0 else -1)-C[k])/e_c[k]
            j = 1-k
            near = C[j]+(yP+d)*e_f[j]
            Y = (outward(near, g, 1 if e_f[j] > 0 else -1)-C[j])/e_f[j]
            inside = lambda x: 0 <= x[0] <= X and 0 <= x[1] <= Y
            clear = True
            for a, b in zip(self.sharp, self.sharp[1:]):
                oa, ob = self.oblique(c, a), self.oblique(c, b)
                if c['corner'] in (a, b):
                    other = ob if a == c['corner'] else oa
                    clear = clear and not inside(other)
                else:
                    clear = clear and not segment_meets(oa, ob, (F(0), F(0)), (X, Y))
            if clear:
                return X, Y
            d /= 2
        raise AssertionError('no box clear of the profile')

    def in_box(self, c, r, w):
        x, y = self.oblique(c, (r, w))
        X, Y = c['box']
        return 0 <= x <= X and 0 <= y <= Y

    def construction_inside(self, r, w):
        """The kernel's construction's membership at a meridian point (exact
        rationals): the lattice's, with each fillet's box less its tube and
        less each round's."""
        f = self.lattice.form(w)
        inside = f is not None and r < f
        for c in self.corners:
            (cr, cw), t = c['centre'], c['t']
            beyond = (r-cr)**2+(w-cw)**2 > t*t
            if self.in_box(c, r, w) and beyond:
                inside = not c['round']
        return inside

    def kite(self, c):
        """The corner's kite: its centre, tangent points and corner, in
        order about it."""
        a, b = c['ends']
        return [c['centre'], a, c['corner'], b]

    def inside(self, r, w):
        """Whether a meridian point lies inside the profile with arcs: inside
        the sharp polygon, less each round's kite beyond its arc, with each
        fillet's kite beyond its arc."""
        pts = [(q(a), q(b)) for a, b in self.sharp]
        wind = False
        for (x0, y0), (x1, y1) in zip(pts, pts[1:]+pts[:1]):
            if (y0 > w) != (y1 > w) and r < x0+(w-y0)*(x1-x0)/(y1-y0):
                wind = not wind
        for c in self.corners:
            kite = [(q(a), q(b)) for a, b in self.kite(c)]
            (cr, cw), t = c['centre'], c['t']
            if in_convex(kite, (r, w)) and (r-q(cr))**2+(w-q(cw))**2 > q(t)**2:
                wind = not c['round']
        return wind

    moment = fillets.Model.moment
    closed = fillets.Model.closed


# ------------------------------------------------------------------ bodies

class Body(fillets.Body):
    """An imported turned body about the world's `z` through `origin`: OCCT's
    profile and the model's, the relation between them checked."""

    def __init__(self, name, occt, model=None, tolerance=None, origin=(0, 0, 0)):
        self.name, self.axis, self.origin = name, 'z', tuple(F(c) for c in origin)
        self.path = f'imported/cone_fillet_{name}.brep'
        self.occt, self.tolerance = occt, tolerance
        self.model = Model(model or occt)
        if model:
            self.relate()
        self.lattice = self.model.lattice
        self.ref_specs = [('turned', name)]
        self.ref_expr = 0

    def native_rows(self):
        o = self.origin
        wire = [f'wire S {len(self.occt)}']
        for r, w, seg in self.occt:
            wire.append(' '.join(num(c) for c in self.place(r, w)))
            if seg[0] == 'L':
                wire.append('L')
            else:
                _, cr, cw, t, ccw = seg
                c = self.place(cr, cw)
                wire.append(f'A {num(c[0])} {num(c[1])} {num(c[2])} {num(t)} {1 if ccw else 0}')
        rows = [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} 0 0 1 {TURN}']
        if self.tolerance:
            rows.append(f'tolerance {self.tolerance}')
        return rows



def rounded(profile):
    """A model profile's points rounded once to binary64 (OCCT's)."""
    return [(F(float(F(r))), F(float(F(w))), seg) for r, w, seg in profile]


def refit():
    """`round`'s profile with its round refitted (OCCT's rows) and the
    model's (tangent): OCCT's round about `(4, 4 + D)` through the
    cylinder's top `(6, 4)`, its radius rounded once, its cone through `(17/4,
    8)` tangent to that circle at a point rounded once; the model's round
    tangent to the cylinder at `4 + D`, its direction's half-angle tangent
    OCCT's on the multiples of 2^-40."""
    with mp.workdps(60):
        radius = float(mp.sqrt(q(4+D*D)))
        C = (q(4), q(4+D))
        Q = (q(F(17, 4)), q(8))
        rho = mp.mpf(radius)
        d = mp.sqrt((Q[0]-C[0])**2+(Q[1]-C[1])**2)
        u = ((Q[0]-C[0])/d, (Q[1]-C[1])/d)
        cb, sb = rho/d, mp.sqrt(1-(rho/d)**2)
        # The tangent point on the side of larger radius.
        cands = [(C[0]+rho*(u[0]*cb+sg*u[1]*sb), C[1]+rho*(u[1]*cb-sg*u[0]*sb)) for sg in (1, -1)]
        T1 = max(cands, key=lambda p: p[0])
        T1 = (F(float(T1[0])), F(float(T1[1])))
        # The stored cone's direction from the centre: normal to OCCT's line
        # through `T1` and `Q`, toward it from the centre.
        dr, dw = q(F(17, 4)-T1[0]), q(8-T1[1])
        n = (dw, -dr)
        length = mp.sqrt(n[0]**2+n[1]**2)
        n = (n[0]/length, n[1]/length)
        if n[0]*(q(T1[0])-C[0])+n[1]*(q(T1[1])-C[1]) < 0:
            n = (-n[0], -n[1])
        s = half_tangent(n)
        scaled = s*2**40
        nearest = mp.nint(scaled)
        assert abs(scaled-nearest) < mp.mpf(1)/4, 'a direction near a half step of the grid'
        s = F(int(nearest), 2**40)
    occt = [(0, 0, L), (6, 0, L), (6, 4, arc(4, 4+D, repr(radius), True)), (T1[0], T1[1], L),
            (F(17, 4), 8, L), (0, 8, L)]
    e = direction(s)
    centre = (F(4), 4+D)
    t1 = (centre[0]+2*e[0], centre[1]+2*e[1])
    top = t1[0]-e[1]/e[0]*(8-t1[1])
    model = [(0, 0, L), (6, 0, L), (6, 4+D, arc(4, 4+D, 2, True)), (t1[0], t1[1], L), (top, 8, L), (0, 8, L)]
    return occt, model


def tangent_profile(start, centre, t, s, ccw, rest):
    """A profile's arc about `centre` of radius `t` from `start` to its point
    in the direction of half-angle tangent `s`, then `rest`."""
    e = direction(s)
    end = (centre[0]+t*e[0], centre[1]+t*e[1])
    return [(start[0], start[1], arc(centre[0], centre[1], t, ccw)), (end[0], end[1], L)]+rest


def make_bodies():
    # The round: its cone through its tangent point `(98/17, 84/17)` normal
    # to `(15/17, 8/17)`, at 8 of radius 1054/255.
    round_ = [(0, 0, L), (6, 0, L)]+tangent_profile((6, 4), (4, 4), 2, F(1, 4), True,
                                                     [(F(1054, 255), 8, L), (0, 8, L)])
    # The fillet: its boss's cone through `(70/17, 43/17)` normal to
    # `(-15/17, -8/17)`, at 7 of radius 26/15.
    fillet = [(0, 0, L), (8, 0, L), (8, 2, L)]+tangent_profile((5, 2), (5, 3), 1, -4, False,
                                                               [(F(26, 15), 7, L), (0, 7, L)])
    # The top: its cone from radius 3 at 0 through `(43/5, 21/5)` normal to
    # `(3/5, -4/5)`.
    e = direction(F(-1, 2))
    t1 = (8+e[0], 5+e[1])
    top = [(0, 0, L), (3, 0, L), (t1[0], t1[1], arc(8, 5, 1, True)), (8, 6, L), (0, 6, L)]
    occt, model = refit()
    out = [
        Body('round', rounded(round_), round_),
        Body('fillet', rounded(fillet), fillet),
        Body('top', rounded(top), top),
        Body('refit', occt, model, tolerance='0.001'),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


flatten = fillets.flatten


class Imported(fillets.Imported):
    """One case, its imported input this step's body."""

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


def group(name, items, outcomes, kind=None, reason=None):
    """Cases of one group: the operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); the group's own refusal declares its cases."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('solid' if n else 'empty')
        out.append(Imported(f'{name}_{op}', 'pieces', items, op, [], k, n, reason if kind else None))
    return out


# The round body's rod along the world's x about (y, z) = (1/2, 9/2): across
# its round's torus alone (between its tangent circles at heights 4 and
# 84/17).
ROUND_ROD = cylinder(at('SIDE', (-10, 0, 0)), 0.5, 4.5, 0.25, 0.0, 20.0)
# A ball across the fillet body's torus alone (clear of its step's plane and
# its boss's cone).
FILLET_BALL = sphere(0.125, at('XY', (4.5, 0, F(69, 32))))
# A box across the top body's bottom corner (radius 3 at height 0), on a
# skew frame (its faces on `TILT`'s would hold the cone's rulings).
TOP_BOX = block(at('SKEW', (3, 0, 0)), -1.5, -1.0, 1.5, 1.0, -1.0, 1.0)
# A box across the refitted body's top corner (radius about 17/4 at 8).
REFIT_BOX = block(at('TILT', (0, 0, 8)), 3.0, -1.0, 5.5, 1.0, -1.0, 1.0)
# A ball touching the round body's round at (148/25, 107/25): its centre
# (4, 4) plus 2 (24/25, 7/25); the ball's centre 25/8 along that direction.
ROUND_TOUCH = sphere(1.125, at('XY', (7, 0, 4.875)))
TOUCH = "a torus tangent to the other input's surface (S9d.4b.2)"


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out = []
    out += group('round_rod', [imported('round'), ROUND_ROD], three)
    out += group('top_box', [imported('top'), TOP_BOX], three)
    out += group('refit_box', [imported('refit'), REFIT_BOX], three)
    out += group('fillet_ball', [imported('fillet'), FILLET_BALL], three)
    out += group('round_touch', [imported('round'), ROUND_TOUCH], {'fuse': 2, 'cut': 1, 'common': 0},
                 kind='degenerate', reason=TOUCH)
    return out


def all_cases():
    return cases()


validate = fillets.validate


# ------------------------------------------------------------------ the reference's turned input

tc = ref.tc
Z = mp.mpf(0)
mq = fillets.mq


class TurnedInput(fillets.TurnedInput):
    """S9e.4b.4d.2c.1a's turned input with cone faces: each cone's band, its
    surface `x^2 + y^2 - (r0 + k w)^2` of exact data, covered by its
    generatrices and its rows."""

    def __init__(self, body):
        self.body = body
        self.model = m = body.model
        self.fr = tc.Frame(body.frame_at(None))
        self.surfs, self.faces = [], []
        pts, segs = m.pts, m.segs
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if (ra, wa) == (rb, wb):
                continue
            if segs[i][0] == 'A':
                self.arc_face(i, segs[i], (ra, wa), (rb, wb))
            elif wa == wb:
                self.plane_face(i, wa, min(ra, rb), max(ra, rb), 1 if rb < ra else -1)
            elif ra == rb:
                self.cylinder_face(i, ra, wa, wb, 1 if wb > wa else -1)
            else:
                self.cone_face(i, (ra, wa), (rb, wb), 1 if wb > wa else -1)

    def cone_face(self, i, a, b, out):
        """The cone through meridian points `a` and `b` over their heights
        (either way), its outward normal `out` times its radial side's."""
        fr = self.fr
        k = (b[0]-a[0])/(b[1]-a[1])
        r0 = a[0]-k*a[1]
        K, R0 = mq(k), mq(r0)
        lo, hi = mq(min(a[1], b[1])), mq(max(a[1], b[1]))

        def f(X):
            p = fr.chart(X)
            rr = R0+K*p[2]
            return p[0]*p[0]+p[1]*p[1]-rr*rr

        def grad(X):
            p = fr.chart(X)
            rr = R0+K*p[2]
            return fr.covector((2*p[0], 2*p[1], -2*K*rr))
        self.surfs.append(tc.Surf(('cone', i), 2, f, grad))
        zero = (Z, Z, Z)

        def generatrix(t):
            c, s = mp.cos(t), mp.sin(t)
            return tc.line_curve(fr.world(R0*c, R0*s, Z), fr.vec(K*c, K*s, 1), fr.vec(-R0*s, R0*c, Z),
                                 fr.vec(-K*s, K*c, Z), lo, hi)

        def row(w):
            rr = R0+K*w
            return tc.circle_curve(fr.world(Z, Z, w), fr.vec(rr, Z, Z), fr.vec(Z, rr, Z), fr.nm, fr.vec(K, Z, Z),
                                   fr.vec(Z, K, Z))
        self.faces.append(tc.Face(('cone', i), [tc.Family('generatrices', Z, 2*mp.pi, True, generatrix, out),
                                                 tc.Family('rows', lo, hi, False, row, -out)]))


class FloatTurned(fillets.FloatTurned):
    """`TurnedInput` in binary64 for the rays and Monte Carlo: its surfaces'
    polynomials along a ray (a cone's too) and its membership (the profile
    with arcs in binary64, each corner's kite)."""

    def __init__(self, S):
        super().__init__(S)
        m = S.model
        self.kites = [([(float(a), float(b)) for a, b in m.kite(c)], (float(c['centre'][0]), float(c['centre'][1])),
                       float(c['t']), c['round']) for c in m.corners]

    def polys(self, P, D_):
        p0, p1 = self.chart(P), self.chart_vec(D_)
        lin_ = lambda i: [p0[i], p1[i]]
        sq = lambda i: tc.fpoly_mul(lin_(i), lin_(i))
        out = []
        for i in range(len(self.pts)-1):
            (ra, wa), (rb, wb) = self.pts[i], self.pts[i+1]
            seg = self.segs[i]
            if seg[0] == 'A':
                _, cr, cw, t, _ = seg
                R, W, T = float(cr), float(cw), float(t)
                zz = [p0[2]-W, p1[2]]
                s = tc.fpoly_add(tc.fpoly_add(sq(0), sq(1)), tc.fpoly_mul(zz, zz))
                s[0] += R*R-T*T
                h = tc.fpoly_add(sq(0), sq(1))
                out.append(tc.fpoly_add(tc.fpoly_mul(s, s), [-4*R*R*c for c in h]))
            elif wa == wb:
                out.append([p0[2]-wa, p1[2]])
            else:
                k = (rb-ra)/(wb-wa)
                r0 = ra-k*wa
                h = tc.fpoly_add(sq(0), sq(1))
                g = [r0+k*p0[2], k*p1[2]]
                out.append(tc.fpoly_add(h, [-c for c in tc.fpoly_mul(g, g)]))
        return out

    def contains(self, X):
        p = self.chart(X)
        r, w = math.sqrt(p[0]*p[0]+p[1]*p[1]), p[2]
        wind = False
        n = len(self.sharp)
        for k in range(n):
            (x0, y0), (x1, y1) = self.sharp[k], self.sharp[(k+1) % n]
            if (y0 > w) != (y1 > w) and r < x0+(w-y0)*(x1-x0)/(y1-y0):
                wind = not wind
        for kite, (cr, cw), t, rnd in self.kites:
            if in_convex(kite, (r, w)) and (r-cr)*(r-cr)+(w-cw)*(w-cw) > t*t:
                wind = not rnd
        return wind


# S9e.4b.4d.2c.1a's (its turned bodies' inputs; each module's own bodies
# its own inputs, both loaded in one process).
_make_input, _Float = ref.make_input, ref.Float


def make_input(case):
    """The chained reference's inputs, this step's turned body's
    `TurnedInput`."""
    body = getattr(case, 'turned', None)
    return TurnedInput(body) if isinstance(body, Body) else _make_input(case)


def make_float(S):
    return FloatTurned(S) if isinstance(S, TurnedInput) else _Float(S)


ref.make_input = make_input
ref.Float = make_float
TurnedCase = fillets.TurnedCase


def reference_case(spec):
    return TurnedCase(BODIES[spec[1]]) if spec[0] == 'turned' else construction(spec, 0)


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3c.3a's, each body's closed form)."""
    name, items, ops, declared, mc_n = job
    specs, sub = flatten(items)
    chain = ref.Chain([reference_case(s) for s in specs])
    size = chain.size
    exprs = {op: ((op, sub[0], sub[1]), solids) for op, solids in ops.items()}
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
    # The imported body one solid, its closed form.
    it, e = items[0], sub[0]
    assert forms.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
    V, A = BODIES[it[1]].closed()
    v, _, a = chain.measures(e)
    checks['body_closed_forms'] = max(abs(v-V)/size**3, abs(a-A)/size**2)
    # The pair identities.
    vx, mx, ax = chain.measures(sub[0])
    vc, mc, ac = chain.measures(sub[1])
    m = {op: chain.measures(e) for op, (e, _) in exprs.items()}
    vf, mf, af = m['fuse']
    vn, mn, an = m['common']
    vt, mt, _ = m['cut']
    dev = abs(vf+vn-vx-vc)/size**3
    for i in range(3):
        dev = max(dev, abs(mf[i]+mn[i]-mx[i]-mc[i])/size**4)
    dev = max(dev, abs(vt-(vx-vn))/size**3, max(abs(mt[i]-(mx[i]-mn[i])) for i in range(3))/size**4)
    checks['pair_identities'] = dev
    checks['area_identity'] = abs(af+an-ax-ac)/size**2
    # Solid counts by rays (not the declared group's: its ball touches the
    # round at a point; nor a result too thin for them, `THIN`).
    for op, (expr, _) in exprs.items():
        for n in () if declared or (name, op) in THIN else (48, 71):
            got = forms.count_solids(chain, expr, n)
            assert got == ops[op], f'{name} {op}: {got} solids by rays ({n}), {ops[op]} declared'
    mcr = ref.monte_carlo(chain, [e for e, _ in exprs.values()], mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op, (expr, _) in exprs.items():
        if (name, op) in THIN:
            continue
        vol, c, sv, sc = mcr[expr]
        v, mom, _ = chain.measures(expr)
        z = max(z, abs(vol-float(v))/max(sv, 1e-12))
        if c is not None and v > 0:
            cent = [float(mom[i]/v) for i in range(3)]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    sine, _, _ = chain.margins()
    margins = {'sine': sine, 'spacing': base.spacing(chain)}
    stats = {'quadrature': chain.quad_error()/size**4, 'missed': chain.missed(), 'refined': chain.refined()}
    return name, rows, res, checks, margins, stats, declared


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, {}, c.kind == 'degenerate'])
        e[1][c.last] = c.solids
    return [(name, items, ops, declared, mc_n) for name, (items, ops, declared) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('TURNED_CONE_FILLETS_CACHE')
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


LIMITS = fillets.LIMITS
# Results thinner than the rays' spacing and smaller than Monte Carlo's
# resolution over the case's box (their counts and Monte Carlo unchecked):
# the ball in the fillet's common.
THIN = {('fillet_ball', 'common')}
reference_checks = fillets.reference_checks


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
           'degenerate case, then result N volume area cx cy cz or empty)']
    for case in all_cases():
        blocks.append(case.encode())
        row = by_group[case.group][1][case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP} {case.klass}')
        if case.kind == 'degenerate':
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
    """Every case's native rows (`compare_turned_cone_fillets_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([reference_case(s) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def profile_distance(body, X):
    """A world point's distance from the body's OCCT profile's lines and
    circles (unbounded), in its axis frame."""
    o = [q(c) for c in body.origin]
    d = [X[i]-o[i] for i in range(3)]
    r, w = mp.sqrt(d[0]**2+d[1]**2), d[2]
    best = None
    pts = [(F(a), F(b)) for a, b, _ in body.occt]
    for ((ra, wa), (rb, wb)), (_, _, seg) in zip(zip(pts, pts[1:]), body.occt):
        if seg[0] == 'A':
            _, cr, cw, t, _ = seg
            radius = mp.mpf(t) if isinstance(t, str) else q(t)
            dist = abs(mp.sqrt((r-q(cr))**2+(w-q(cw))**2)-radius)
        else:
            a, b = q(rb-ra), -q(wb-wa)
            dist = abs(a*(w-q(wa))+b*(r-q(ra)))/mp.sqrt(a*a+b*b)
        best = dist if best is None else min(best, dist)
    return best


def check_bodies():
    """Each existing body file: its faces' kinds its OCCT profile's, each
    stored vertex on that profile's lines and circles within 1e-12 of the
    body's size. Returns the largest deviation and the files read."""
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
        os.environ['TURNED_CONE_FILLETS_CACHE'] = args.cache
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))
    print('reference checks (largest deviation, relative to the case size; least_sine and event_spacing as '
          'reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print(f'bodies: {read} of {len(BODIES)} files read, every stored vertex within {mp.nstr(deviation, 3)} of the '
          'size of its profile\'s lines and circles')


if __name__ == '__main__':
    main()
