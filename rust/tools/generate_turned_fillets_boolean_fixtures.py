#!/usr/bin/env python3
"""Fixtures for S9e.4b.4d.2c.1a of REVIEW_NOTES.md: imported turned bodies
whose fillets and rounds (tori) join a plane and a cylinder smoothly (the
kind of `OCC485a`'s fillet under its step and its round at its top), each
torus snapped tangent to its two lines, given to Booleans whose partners
cross no smooth edge.

The bodies are OCCT's own output, as S9e.4b.4d.2a's
(`generate_turned_profiles_boolean_fixtures.py`):
`boolean-turned-fillets-bodies.txt` holds one `write NAME
imported/fillet_NAME.brep` block per body for `occt_boolean_oracle.cpp` (a
`plane` row, the profile's `wire` of lines and arcs (`A cx 0 cz r ccw` to
the next point, the frame's normal `-y` so that a counter-clockwise arc in
`(x, z)` is the flag 1), a `revolve` row turning it a whole turn about the
world's `z` (BRepPrimAPI_MakeRevol), and for the refitted body a `tolerance`
row setting every tolerance to 1e-3 as a translated part stores its own),
which `compare_turned_fillets_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/fillet_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-turned-fillets-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/fillet_NAME.brep` row;
`boolean-turned-fillets-expected.tsv` gives per case `expect KIND
S9e.4b.4d.2c.1a CLASS` (`solid`, `empty`, `degenerate`: a `Degenerate` of
S9's rules; the class `pieces`: an imported body and a construction), a
degenerate case's `reason TEXT`, then `result N volume area cx cy cz`
(totals over the N solids, world coordinates) or `empty`;
`boolean-turned-fillets-frames.tsv` the stored axes of every solid the
kernel builds from its rows.

Each body has two profiles of dyadic points `(r, w)` (`r` the radius, `w`
the height along the axis), from the axis out and back to it, each segment
a line (a plane where level, else a cylinder) or an arc (a torus) to the
next point: the one OCCT is given (`occt`) and the model's (`model`): each
arc a fillet (concave, its material outside the tube) or a round (convex)
between a plane and a cylinder, tangent to both. The bodies, each about the
world's `z` through the origin: `step` (`turned_fillet`'s and `OCC485a`'s
kind): a disc of radius 8, a cylinder to 3, a step in to radius 9/2, a
fillet of radius 1/2 about `(9/2, 7/2)` to the shaft of radius 4, the shaft
to 7, a round of radius 5/4 about `(11/4, 7)` to the top's disc at 33/4;
`collar`: a round of radius 1 about `(2, 1)` from the bottom's disc to a
shaft of radius 3, the shaft to 9/2, a fillet of radius 1/2 about `(7/2,
9/2)` to a step out at 5 (an annulus facing down), a collar of radius 6 to
7, a round of radius 1 about `(5, 7)` to the top's disc at 8; `refit`
(`OCC485a`'s round's kind): `step`'s profile whose round OCCT is given with
its centre raised by `D` = 2^-20 (9.5e-7) above the shaft's top at 7, its
radius the distance to that point (`sqrt((5/4)^2 + D^2)`, rounded once:
5/4 + 3.6e-13) and the top at that radius above its centre (rounded once),
at a resolution of 1e-3: its tube crosses the shaft by 3.6e-13 and its
stored edge with the shaft lies `D` below the model's tangency, within the
resolution; the model's round is tangent to the shaft and the top's plane at
`7 + D` and `33/4 + D`.

The kernel's construction (as REVIEW_NOTES' "S9e.4b.4d.2c.1 refined"
decides; `Model.construction_specs`, checked here against the profile at
points of the meridian, not swept): the sharp profile, each arc replaced by
the corner its two lines meet at, is S9e.4b.4d.2a's lattice
(`profiles.Lattice`); each fillet's piece fused to it and each round's taken
from it, a piece the box (an annular prism about the axis) from the tube's
centre past the corner by half the tube less the whole torus.

Cases (each the three operations): `step_rod`, the step and a rod along the
world's `x` across its round alone (between its tangent circles); `step_ball`,
the step and a ball across its fillet alone; `collar_box`, the collar and a
`TILT` box across the sharp corner of its step out and its collar (away
from its fillet's and rounds' circles); `refit_box`, the refitted body and a
`TILT` box across its bottom's corner; and declared `degenerate`,
`step_touch`, the step and a ball touching its round at a point
(`Degenerate("a torus tangent to the other input's surface (S9d.4b.2)")`).

The reference is S9e.3a's chained reference
(`chained_curved_boolean_reference.py`) with each body one input turned
from its model profile (`TurnedInput`, none of the kernel's construction:
its faces the profile's own, each plane's annulus or disc, each cylinder's
band and each arc's band of its whole torus's tube, each covered by two
families of the reference's circles and lines; its membership the profile
with arcs in the meridian; the construction swept instead, its lattice's,
boxes' and tori's faces, ran 50 minutes on five workers
without finishing a case), then `(X) op C`, with
S9e.4b.3c.3a's checks (the two families, each solid's closed form, the pair
identities, Monte Carlo, the body one solid, and outside the declared group
solid counts by rays at two resolutions, the meetings' sines and the events'
spacing; a result too thin for the rays and Monte Carlo, `THIN`, unchecked
by them) and each body's closed form by Pappus over its model profile (its
lines' frusta and annuli, and each arc's by Green's theorem in the meridian,
exactly but for pi and the slant heights' square roots, its moment along
the axis likewise). With `--check`, where the bodies' files exist, each is read
(`stored_records`): its faces' kinds its OCCT profile's, every stored vertex
within 1e-12 of the size on the OCCT profile's lines and circles.
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
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_turned_profiles_boolean_fixtures as profiles

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces',)
PREFIX = 'boolean-turned-fillets'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2c.1a'
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
torus = base.torus
ref = base.ref
q = forms.q
cylinder, block = forms.cylinder, forms.block
Lattice, num, at, TURN = profiles.Lattice, profiles.num, profiles.at, profiles.TURN
# The refitted round's offset.
D = F(1, 2**20)


# ------------------------------------------------------------------ profiles

L = ('L',)


def arc(cr, cw, t, ccw):
    """An arc to the next point about `(cr, cw)` of radius `t` (a Fraction,
    or for OCCT a binary64 number's text), counter-clockwise in `(r, w)`
    where `ccw`."""
    return ('A', F(cr), F(cw), t if isinstance(t, str) else F(t), ccw)


def points(profile):
    return [(F(r), F(w)) for r, w, _ in profile]


def meet(a0, a1, b0, b1):
    """The common point of the line through `a0`, `a1` and that through
    `b0`, `b1` (exact)."""
    (x1, y1), (x2, y2) = a0, a1
    (x3, y3), (x4, y4) = b0, b1
    den = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
    assert den != 0, 'parallel lines at an arc'
    a, b = x1*y2-y1*x2, x3*y4-y3*x4
    return ((a*(x3-x4)-(x1-x2)*b)/den, (a*(y3-y4)-(y1-y2)*b)/den)


def segment_meets(a, b, lo, hi):
    """Whether the segment from `a` to `b` meets the closed box `[lo, hi]`
    (exactly)."""
    t0, t1 = F(0), F(1)
    for k in range(2):
        d = b[k]-a[k]
        if d == 0:
            if a[k] < lo[k] or a[k] > hi[k]:
                return False
            continue
        e0, e1 = (lo[k]-a[k])/d, (hi[k]-a[k])/d
        if e0 > e1:
            e0, e1 = e1, e0
        t0, t1 = max(t0, e0), min(t1, e1)
        if t0 > t1:
            return False
    return True


class Model:
    """A model profile with arcs: its sharp profile's lattice, each arc's
    corner (`centre`, `t`, `round`, `corner`) and piece's box `[(r, w) low,
    (r, w) high]`."""

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
                # Each neighbour a plane (level) or a cylinder (upright).
                for p, s in ((a0, a1), (b0, b1)):
                    assert p[0] == s[0] or p[1] == s[1], 'a torus beside a cone'
                corner = meet(a0, a1, b0, b1)
                # Tangent to both at its ends (the centre's feet).
                for p, s, end in ((a0, a1, a1), (b0, b1, b0)):
                    foot = (p[0], cw) if p[0] == s[0] else (cr, p[1])
                    assert foot == end, ('an arc not tangent to its line', foot, end)
                    assert (end[0]-cr)**2+(end[1]-cw)**2 == t*t, 'an end off its arc'
                corners.append({'centre': (cr, cw), 't': t, 'round': ccw, 'corner': corner})
                sharp.append(corner)
                i += 2
            else:
                sharp.append(pts[i])
                i += 1
        self.sharp, self.corners = sharp, corners
        self.lattice = Lattice(sharp)
        for c in corners:
            c['box'] = self.box(c)

    def box(self, c):
        """The piece's box: from the tube's centre past the corner by half
        the tube, clear of the sharp profile's other segments and of the
        other boxes."""
        centre, corner, d = c['centre'], c['corner'], c['t']/2
        toward = [1 if corner[k] > centre[k] else -1 for k in range(2)]
        far = tuple(corner[k]+toward[k]*d for k in range(2))
        lo = tuple(min(centre[k], far[k]) for k in range(2))
        hi = tuple(max(centre[k], far[k]) for k in range(2))
        inside = lambda x: lo[0] <= x[0] <= hi[0] and lo[1] <= x[1] <= hi[1]
        for a, b in zip(self.sharp, self.sharp[1:]):
            if corner in (a, b):
                other = b if a == corner else a
                assert not inside(other), 'a line at the corner ending inside its box'
            else:
                assert not segment_meets(a, b, lo, hi), 'a box across another segment'
        assert lo[0] > 0, 'a box reaching the axis'
        return lo, hi

    def construction_specs(self, frame_at):
        """The kernel's construction's solids: the lattice's, then each corner's box
        (the kernel's annular prism as the cylinder of its outer radius less
        that of its inner one, the latter past its ends by an eighth of the
        tube: the reference's prisms are of one convex profile) and its whole
        torus."""
        out = list(self.lattice.specs(frame_at))
        for c in self.corners:
            lo, hi = c['box']
            past = c['t']/8
            out.append(cylinder(frame_at(None), 0.0, 0.0, float(hi[0]), float(lo[1]), float(hi[1])))
            out.append(cylinder(frame_at(None), 0.0, 0.0, float(lo[0]), float(lo[1]-past), float(hi[1]+past)))
            out.append(torus(float(c['centre'][0]), float(c['t']), frame_at(c['centre'][1])))
            for x in (*lo, *hi, *c['centre'], c['t'], lo[1]-past, hi[1]+past):
                assert float(x) == x, ('a piece off the binary64 numbers', x)
        return out

    def construction_expr(self, k=0):
        """The lattice fused with each fillet's piece (its box less its
        torus) and less each round's."""
        e = self.lattice.expr(k)
        j = k+len(self.lattice.parts)
        for c in self.corners:
            piece = ('cut', ('cut', j, j+1), j+2)
            e = ('cut', e, piece) if c['round'] else ('fuse', e, piece)
            j += 3
        return e

    def inside(self, r, w):
        """Whether a meridian point lies inside the profile with arcs: inside
        the sharp polygon, less each round's corner beyond its arc, with each
        fillet's corner beyond its arc."""
        pts = [(q(a), q(b)) for a, b in self.sharp]
        wind = False
        for (x0, y0), (x1, y1) in zip(pts, pts[1:]+pts[:1]):
            if (y0 > w) != (y1 > w) and r < x0+(w-y0)*(x1-x0)/(y1-y0):
                wind = not wind
        for c in self.corners:
            (cr, cw), (pr, pw), t = c['centre'], c['corner'], c['t']
            if q(min(cr, pr)) < r < q(max(cr, pr)) and q(min(cw, pw)) < w < q(max(cw, pw)):
                if (r-q(cr))**2+(w-q(cw))**2 > q(t)**2:
                    wind = not c['round']
        return wind

    def moment(self):
        """The profile's `int r w dA = int r^2 w / 2 dw` along its boundary
        (counter-clockwise): a line's exactly (Simpson's rule, exact for its
        cubic), an arc's by quadrature."""
        out = mp.mpf(0)
        pts, segs = self.pts, self.segs
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if segs[i][0] == 'A':
                _, cr, cw, t, ccw = segs[i]
                R, W, T = q(cr), q(cw), q(t)
                th0 = mp.atan2(q(wa-cw), q(ra-cr))
                th1 = mp.atan2(q(wb-cw), q(rb-cr))
                sweep = th1-th0
                if ccw and sweep < 0:
                    sweep += 2*mp.pi
                if not ccw and sweep > 0:
                    sweep -= 2*mp.pi
                out += mp.quad(lambda th: (R+T*mp.cos(th))**2*(W+T*mp.sin(th))/2*T*mp.cos(th),
                               [th0, th0+sweep])
            elif wb != wa:
                f = lambda s: (ra+(rb-ra)*s)**2*(wa+(wb-wa)*s)/2
                out += q((wb-wa)*(f(F(0))+4*f(F(1, 2))+f(F(1)))/6)
        return out

    def closed(self):
        """Pappus: `(volume, area)` of the profile turned about the axis: the
        volume `2 pi` times the profile's `int r dA = int r^2 / 2 dw` along
        its boundary (counter-clockwise; the axis adds nothing), the area `2
        pi int r ds` along its path, each arc's in closed form."""
        vol, area = mp.mpf(0), mp.mpf(0)
        pts, segs = self.pts, self.segs
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if segs[i][0] == 'A':
                _, cr, cw, t, ccw = segs[i]
                th0 = mp.atan2(q(wa-cw), q(ra-cr))
                th1 = mp.atan2(q(wb-cw), q(rb-cr))
                sweep = th1-th0
                # The short way about the centre, as the arc turns.
                if ccw and sweep < 0:
                    sweep += 2*mp.pi
                if not ccw and sweep > 0:
                    sweep -= 2*mp.pi
                th1 = th0+sweep
                R, W, T = q(cr), q(cw), q(t)
                # int r^2 / 2 dw, r = R + T cos, w = W + T sin, dw = T cos.
                def prim(th):
                    c, s = mp.cos(th), mp.sin(th)
                    return T*(R*R*s/2+R*T*(th/2+s*c/2)+T*T*(s-s**3/3)/2)
                vol += prim(th1)-prim(th0)
                area += abs(T*(R*sweep+T*(mp.sin(th1)-mp.sin(th0))))
            elif wb != wa:
                vol += q(F(1, 6)*(wb-wa)*(ra*ra+ra*rb+rb*rb))
                area += q(ra+rb)/2*mp.sqrt(q((rb-ra)**2+(wb-wa)**2))
            else:
                area += q(abs(rb*rb-ra*ra))/2
        return 2*mp.pi*vol, 2*mp.pi*area


# ------------------------------------------------------------------ bodies

class Body:
    """An imported turned body about the world's `z` through `origin`: OCCT's
    profile and the model's, the relation between them checked, the model's
    construction."""

    def __init__(self, name, occt, model=None, tolerance=None, origin=(0, 0, 0)):
        self.name, self.axis, self.origin = name, 'z', tuple(F(c) for c in origin)
        self.path = f'imported/fillet_{name}.brep'
        self.occt, self.tolerance = occt, tolerance
        self.model = Model(model or occt)
        if model:
            self.relate()
        self.lattice = self.model.lattice
        # The reference's input: the body turned from its model profile
        # (`TurnedInput`), one input.
        self.ref_specs = [('turned', name)]
        self.ref_expr = 0

    def relate(self):
        """The model's profile is OCCT's within `4 D` point by point, its
        segments of the same kinds."""
        assert len(self.occt) == len(self.model.profile), self.name
        for (ra, wa, sa), (rb, wb, sb) in zip(self.occt, self.model.profile):
            assert abs(F(ra)-rb) <= 4*D and abs(F(wa)-wb) <= 4*D, (self.name, ra, wa)
            assert sa[0] == sb[0], self.name

    def place(self, r, w):
        o = self.origin
        return (o[0]+F(r), o[1], o[2]+F(w))

    def frame_at(self, w):
        return at('XY', self.place(0, 0 if w is None else w))

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
                radius = t if isinstance(t, str) else num(t)
                wire.append(f'A {num(c[0])} {num(c[1])} {num(c[2])} {radius} {1 if ccw else 0}')
        rows = [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} 0 0 1 {TURN}']
        if self.tolerance:
            rows.append(f'tolerance {self.tolerance}')
        return rows

    def closed(self):
        return self.model.closed()

    def kinds(self):
        """Its stored faces' surface kinds: OCCT's profile's segments'."""
        out = []
        pts = [(F(r), F(w)) for r, w, _ in self.occt]
        for (a, b), (_, _, seg) in zip(zip(pts, pts[1:]), self.occt):
            if seg[0] == 'A':
                out.append('torus')
            else:
                out.append('plane' if a[1] == b[1] else ('cylinder' if a[0] == b[0] else 'cone'))
        return sorted(out)


def refit():
    """`step`'s profile with its round refitted (OCCT's rows) and the model's
    (tangent): OCCT's round about `(11/4, 7 + D)` through the shaft's top `(4,
    7)`, its radius and the top's height rounded once."""
    radius = float(mp.sqrt(q(F(5, 4)**2+D*D)))
    top = float(7+D+F(radius))
    occt = [(0, 0, L), (8, 0, L), (8, 3, L), (F(9, 2), 3, arc(F(9, 2), F(7, 2), F(1, 2), False)), (4, F(7, 2), L),
            (4, 7, arc(F(11, 4), 7+D, repr(radius), True)), (F(11, 4), top, L), (0, top, L)]
    model = [(0, 0, L), (8, 0, L), (8, 3, L), (F(9, 2), 3, arc(F(9, 2), F(7, 2), F(1, 2), False)), (4, F(7, 2), L),
             (4, 7+D, arc(F(11, 4), 7+D, F(5, 4), True)), (F(11, 4), F(33, 4)+D, L), (0, F(33, 4)+D, L)]
    return occt, model


def make_bodies():
    step = [(0, 0, L), (8, 0, L), (8, 3, L), (F(9, 2), 3, arc(F(9, 2), F(7, 2), F(1, 2), False)),
            (4, F(7, 2), L), (4, 7, arc(F(11, 4), 7, F(5, 4), True)), (F(11, 4), F(33, 4), L),
            (0, F(33, 4), L)]
    collar = [(0, 0, L), (2, 0, arc(2, 1, 1, True)), (3, 1, L), (3, F(9, 2), arc(F(7, 2), F(9, 2), F(1, 2), False)),
              (F(7, 2), 5, L), (6, 5, L), (6, 7, arc(5, 7, 1, True)), (5, 8, L), (0, 8, L)]
    occt, model = refit()
    out = [
        Body('step', step),
        Body('collar', collar),
        Body('refit', occt, model, tolerance='0.001'),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its
    construction's) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            exprs.append(len(specs))
            specs.append(('turned', it[1]))
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(pieces.Imported):
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

    def expr(self):
        _, e = flatten(self.items)
        return (self.op1, e[0], e[1])


def group(name, items, outcomes, kind=None, reason=None):
    """Cases of one group: the operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); the group's own refusal declares its cases."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('solid' if n else 'empty')
        out.append(Imported(f'{name}_{op}', 'pieces', items, op, [], k, n, reason if kind else None))
    return out


# The step's rod along the world's x about (y, z) = (1/2, 61/8): across its
# round's torus alone (between its tangent circles at heights 7 and 33/4).
STEP_ROD = cylinder(at('SIDE', (-10, 0, 0)), 0.5, 7.625, 0.25, 0.0, 20.0)
# A ball across the step's fillet's torus alone (1/4 off its plane and its
# shaft, 2^(-1/2)/2 from the tube's centre circle).
STEP_BALL = sphere(0.1875, at('XY', (4.25, 0, 3.25)))
# A box across the collar's sharp corner (radius 6 at height 5), its step
# out's annulus from radius 5 and its collar to 6.4, clear of its fillet's
# circle (radius 7/2) and its round's (height 7).
COLLAR_BOX = block(at('TILT', (0, 0, 5)), 5.0, -1.0, 7.5, 1.0, -1.0, 1.0)
# A box across the refitted body's bottom corner (radius 8 at height 0).
REFIT_BOX = block(at('TILT', (0, 0, 0)), 6.5, -1.0, 9.5, 1.0, -1.0, 1.0)
# A ball touching the step's round at (7/2, 0, 8): the round's centre (11/4,
# 7) plus 5/4 (3/5, 4/5); the ball's centre 15/8 along that direction.
STEP_TOUCH = sphere(0.625, at('XY', (3.875, 0, 8.5)))
TOUCH = "a torus tangent to the other input's surface (S9d.4b.2)"


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out = []
    out += group('step_rod', [imported('step'), STEP_ROD], three)
    out += group('collar_box', [imported('collar'), COLLAR_BOX], three)
    out += group('refit_box', [imported('refit'), REFIT_BOX], three)
    out += group('step_ball', [imported('step'), STEP_BALL], three)
    out += group('step_touch', [imported('step'), STEP_TOUCH], {'fuse': 2, 'cut': 1, 'common': 0},
                 kind='degenerate', reason=TOUCH)
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.items) == 2 and not c.stages, c.name
        assert c.last in OPS, c.name
        assert c.items[0][0] == 'imported', c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone', 'turned'), c.name


# ------------------------------------------------------------------ the reference's turned input

tc = ref.tc
Z = mp.mpf(0)


def mq(x):
    return q(x)


class TurnedInput:
    """A body turned about its axis as one input of the chained reference
    (S9e.3a's `Chain`): its faces the model profile's own (each plane's
    annulus or disc, each cylinder's band, each arc's band of its whole
    torus's tube), its surfaces theirs, its membership the profile with
    arcs in the meridian (`Model.inside`), its closed form by Pappus. Every
    surface is the reference's own of exact data (`tc`'s families of circles
    and lines), none of the kernel's construction."""
    kind = 'turned'

    def __init__(self, body):
        self.body = body
        self.model = m = body.model
        self.fr = fr = tc.Frame(body.frame_at(None))
        self.surfs, self.faces = [], []
        pts, segs = m.pts, m.segs
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if (ra, wa) == (rb, wb):
                continue
            # The outward normal in the meridian `(r, w)`: the profile runs
            # counter-clockwise, its material on the left.
            if segs[i][0] == 'A':
                self.arc_face(i, segs[i], (ra, wa), (rb, wb))
            elif wa == wb:
                self.plane_face(i, wa, min(ra, rb), max(ra, rb), 1 if rb < ra else -1)
            else:
                assert ra == rb, 'a cone in a fillets profile'
                self.cylinder_face(i, ra, wa, wb, 1 if wb > wa else -1)

    def plane_face(self, i, h, r0, r1, out):
        """The plane `w = h` over radii `r0` to `r1`, its outward normal `out`
        times the axis."""
        fr = self.fr
        H = mq(h)
        self.surfs.append(tc.plane_surf(('plane', i), fr, (F(0), F(0), F(1)), F(h)))
        zero = (Z, Z, Z)

        def ring(a):
            return tc.circle_curve(fr.world(Z, Z, H), fr.vec(a, Z, Z), fr.vec(Z, a, Z), zero, fr.vec(1, Z, Z),
                                   fr.vec(Z, 1, Z))

        def spoke(t):
            c, s = mp.cos(t), mp.sin(t)
            return tc.line_curve(fr.world(Z, Z, H), fr.vec(c, s, Z), zero, fr.vec(-s, c, Z), mq(r0), mq(r1))
        self.faces.append(tc.Face(('plane', i), [tc.Family('rings', mq(r0), mq(r1), False, ring, out),
                                                  tc.Family('spokes', Z, 2*mp.pi, True, spoke, -out)]))

    def cylinder_face(self, i, r, w0, w1, out):
        """The cylinder of radius `r` over heights `w0` to `w1` (either way),
        its outward normal `out` times the radial direction."""
        fr = self.fr
        R = mq(r)
        lo, hi = mq(min(w0, w1)), mq(max(w0, w1))

        def f(X):
            p = fr.chart(X)
            return p[0]*p[0]+p[1]*p[1]-R*R

        def grad(X):
            p = fr.chart(X)
            return fr.covector((2*p[0], 2*p[1], Z))
        self.surfs.append(tc.Surf(('cylinder', i), 2, f, grad))
        zero = (Z, Z, Z)

        def generatrix(t):
            c, s = mp.cos(t), mp.sin(t)
            return tc.line_curve(fr.world(R*c, R*s, Z), fr.nm, fr.vec(-R*s, R*c, Z), zero, lo, hi)

        def row(w):
            return tc.circle_curve(fr.world(Z, Z, w), fr.vec(R, Z, Z), fr.vec(Z, R, Z), fr.nm, zero, zero)
        self.faces.append(tc.Face(('cylinder', i), [tc.Family('generatrices', Z, 2*mp.pi, True, generatrix, out),
                                                     tc.Family('rows', lo, hi, False, row, -out)]))

    def arc_face(self, i, seg, a, b):
        """The band of the whole torus about the arc's centre over the arc's
        angles in the meridian (`v` from the radial direction toward the
        axis's), its outward normal the tube's for a round, inward for a
        fillet."""
        fr = self.fr
        _, cr, cw, t, ccw = seg
        R, W, T = mq(cr), mq(cw), mq(t)
        th0 = mp.atan2(mq(a[1]-cw), mq(a[0]-cr))
        th1 = mp.atan2(mq(b[1]-cw), mq(b[0]-cr))
        lo, hi = (th0, th1) if ccw else (th1, th0)
        if hi < lo:
            hi += 2*mp.pi
        assert hi-lo < mp.pi, 'an arc of half a turn or more'
        k = R*R-T*T

        def f(X):
            p = fr.chart(X)
            z = p[2]-W
            s = p[0]*p[0]+p[1]*p[1]+z*z
            return (s+k)**2-4*R*R*(p[0]*p[0]+p[1]*p[1])

        def grad(X):
            p = fr.chart(X)
            z = p[2]-W
            s = p[0]*p[0]+p[1]*p[1]+z*z
            g = 4*(s+k)
            return fr.covector((g*p[0]-8*R*R*p[0], g*p[1]-8*R*R*p[1], g*z))
        self.surfs.append(tc.Surf(('torus', i), 4, f, grad))
        out = 1 if ccw else -1
        zero = (Z, Z, Z)

        def tube(u):
            c, s = mp.cos(u), mp.sin(u)
            e, de = fr.vec(c, s, Z), fr.vec(-s, c, Z)
            return tc.circle_curve(tc.lin((1, fr.world(Z, Z, W)), (R, e)), tc.scale(e, T), tc.scale(fr.nm, T),
                                   tc.scale(de, R), tc.scale(de, T), zero, rng=(lo, hi))

        def parallel(v):
            c, s = mp.cos(v), mp.sin(v)
            rho = R+T*c
            return tc.circle_curve(fr.world(Z, Z, W+T*s), fr.vec(rho, Z, Z), fr.vec(Z, rho, Z),
                                   tc.scale(fr.nm, T*c), fr.vec(-T*s, Z, Z), fr.vec(Z, -T*s, Z))
        self.faces.append(tc.Face(('torus', i), [tc.Family('tube circles', Z, 2*mp.pi, True, tube, out),
                                                  tc.Family('parallels', lo, hi, False, parallel, -out)]))

    def contains(self, X):
        p = self.fr.chart(X)
        return self.model.inside(mp.sqrt(p[0]*p[0]+p[1]*p[1]), p[2])

    def closed(self):
        """Volume, world moments and area by Pappus: `2 pi` times the
        profile's `int r dA` and `int r w dA` (Green's theorem along its
        boundary), its area `2 pi int r ds`."""
        V, A = self.model.closed()
        rw = self.model.moment()
        wbar = rw/(V/(2*mp.pi))
        mom = tuple(c*V for c in self.fr.world(Z, Z, wbar))
        return V, mom, A

    def box(self):
        m = self.model
        r = float(max(x for x, _ in m.pts))
        ws = [float(w) for _, w in m.pts]
        return tc.frame_box(self.fr, [(u, v, w) for u in (-r, r) for v in (-r, r) for w in (min(ws), max(ws))])


class FloatTurned(tc.FloatModel):
    """`TurnedInput` in binary64 for the rays and Monte Carlo: its surfaces'
    polynomials along a ray and its membership (the profile with arcs in
    binary64, `math.sqrt` alone)."""

    def __init__(self, S):
        self.S = S
        fr = S.fr
        self.o = [float(v) for v in fr.om]
        self.inv = [[float(v) for v in row] for row in fr.invm]
        self.kind = 'turned'
        m = S.model
        self.pts = [(float(r), float(w)) for r, w in m.pts]
        self.segs = m.segs
        self.sharp = [(float(r), float(w)) for r, w in m.sharp]
        self.corners = [((float(c['centre'][0]), float(c['centre'][1])), (float(c['corner'][0]),
                         float(c['corner'][1])), float(c['t']), c['round']) for c in m.corners]

    def polys(self, P, D):
        p0, p1 = self.chart(P), self.chart_vec(D)
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
                h = tc.fpoly_add(sq(0), sq(1))
                h[0] -= ra*ra
                out.append(h)
        return out

    def contains(self, X):
        p = self.chart(X)
        return inside_float(self.sharp, self.corners, math.sqrt(p[0]*p[0]+p[1]*p[1]), p[2])


def inside_float(sharp, corners, r, w):
    wind = False
    n = len(sharp)
    for k in range(n):
        (x0, y0), (x1, y1) = sharp[k], sharp[(k+1) % n]
        if (y0 > w) != (y1 > w) and r < x0+(w-y0)*(x1-x0)/(y1-y0):
            wind = not wind
    for (cr, cw), (pr, pw), t, rnd in corners:
        if min(cr, pr) < r < max(cr, pr) and min(cw, pw) < w < max(cw, pw):
            if (r-cr)*(r-cr)+(w-cw)*(w-cw) > t*t:
                wind = not rnd
    return wind


_make_input, _Float = ref.make_input, ref.Float


def make_input(case):
    """The chained reference's inputs, a turned body's `TurnedInput`."""
    body = getattr(case, 'turned', None)
    return TurnedInput(body) if body is not None else _make_input(case)


def make_float(S):
    return FloatTurned(S) if S.kind == 'turned' else _Float(S)


ref.make_input = make_input
ref.Float = make_float


class TurnedCase:
    """A turned body as the chained reference's case."""

    def __init__(self, body):
        self.turned = body


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
    cache = os.environ.get('TURNED_FILLETS_CACHE')
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


LIMITS = profiles.LIMITS
# Results thinner than the rays' spacing and smaller than Monte Carlo's
# resolution over the case's box (their counts and Monte Carlo unchecked):
# the ball in the fillet's common, a lens about 0.05 thick of volume 1.1e-3.
THIN = {('step_ball', 'common')}


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # The declared group's ball touches the round: its checks kept
        # apart, looser.
        for key, value in checks.items():
            note(('degenerate_'+key) if declared and key != 'monte_carlo_sigma' else key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not declared:
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
    """Every case's native rows (`compare_turned_fillets_boolean.py`)."""
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
        os.environ['TURNED_FILLETS_CACHE'] = args.cache
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
