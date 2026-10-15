#!/usr/bin/env python3
"""Fixtures for S9e.4b.4d.2c.1c of REVIEW_NOTES.md: imported turned bodies
whose rounds (tori) are tangent to one line of the meridian and meet the
other at a corner (`OCC496b`'s first two tori: a round tangent to the wider
of two coaxial cylinders, its tube crossing the narrower before its top),
given to Booleans whose partners cross the corner rings.

The bodies are OCCT's own output, as S9e.4b.4d.2c.1a's to 2c.2's:
`boolean-turned-corners-bodies.txt` holds one `write NAME
imported/corner_NAME.brep` block per body for `occt_boolean_oracle.cpp` (a
`plane` row, the profile's `wire` of lines and arcs, a `revolve` row turning
it a whole turn about the body's axis, and for `shaft` a `tolerance` row
setting every tolerance to 3e-3, `OCC496b`'s), which
`compare_turned_corners_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/corner_NAME.brep`. Nothing here reads the kernel;
the files are read only to check that each is the body it claims to be
(below). The cases' files are S9e.4b.4d.2c.1a's, under the prefix
`boolean-turned-corners`, the class `pieces`.

Each body has two profiles `(r, w)` (`r` the radius, `w` the height along
the axis), from the axis out and back to it, each segment a line (a plane
where level, a cylinder where upright, else a cone) or an arc (a torus) to
the next point: the one OCCT is given (`occt`, binary64 points: the model's
rounded once where it is no binary64 number) and the model's (`model`,
exact; a corner's height `a + b sqrt(d)`, a `Surd`). A corner's round is
tangent to the wider cylinder at the height of its centre and its tube
crosses the narrower one at `sin v = sqrt(1 - cos^2 v)`, `cos v = (r_Y -
r_c) / t` (REVIEW_NOTES' "S9e.4b.4d.2c.1c refined"). The bodies:

`shaft` (`OCC496b`'s profile about its own axis, the world's `x` through `(0,
-2^-21, OZ)`, in its own numbers on the multiples of `2^-30`, at its
resolution 3e-3): a disc at `W0`, the cylinder of radius `15 + 2^-21`, a
round of radius 2 about `(13 + 2^-21, C1)` crossing the cylinder of radius
`R14` a sixth of a turn round (within 5e-7), that cylinder, a round of
radius 2 about `(R14 - 2, C2)` (OCCT's arc from `9 2^-19` below its centre,
`OCC496b`'s stored edge, the tube 7.4e-11 off that point) crossing the
cylinder of radius `13 + 2^-21`, that cylinder, a round of radius 9/2 about
`(17/2 + 2^-21, C3)` (`OCC496b`'s centre moved 9.5e-7 onto the tangent: an
arc's ends must lie on its circle) to `OCC496b`'s stored angle `V3`, a cone
from there to the top circle `(RT, WP)` (`OCC496b`'s: its line 1.8e-7 rad
off the arc's tangent), the top's disc; the model's third round turned to
the direction whose half-angle tangent is OCCT's cone's normal's on the
multiples of 2^-40 (2c.1b's snap), its cone that direction's tangent;
`step`: a disc of radius 4, the cylinder to 3, a round of radius 1 about
`(3, 3)` crossing the shaft of radius 7/2 at `3 + sqrt(3)/2`, the shaft to
6, the top's disc; `collar` (tangent at its upper end, a Pythagorean
corner): a disc of radius 7/2, the shaft to the corner at 4, a round of
radius 5/4 about `(11/4, 5)` from there (the direction `(3/5, -4/5)`) to the
collar of radius 4, tangent at 5, the collar to 8, the top's disc.

The kernel's construction (`Model.construction_radius`, checked here
against the profile at heights of the meridian, not swept): the sharp
profile, each corner's round replaced by a step between its two cylinders
at the height halfway from the corner to the nearer of the tube's top (its
bottom) and the narrower cylinder's far end on the multiples of `2^-24` of
the tube's power of two, and each tangent round by its lines' corner, is
S9e.4b.4d.2a's lattice (`profiles.Lattice`); each corner's piece (the box
from the centre past the step at the wider cylinder by half the tube, less
the whole torus) taken from the wider cylinder's part alone.

Cases (each the three operations): `shaft_box`, a box across the shaft's
first corner ring (its sides `y = +-1`); `shaft_rod`, a rod along the
world's `y` across that ring near its top; `shaft_skew`, a cube on the skew
frame `SKEW` about the second corner ring's point (evaluating turned);
`step_box`, `step_rod`, `step_skew`: a box, a rod along `x` and a skew cube
across `step`'s corner ring; `collar_box` and `collar_ball`: a box and a
ball (about a point 1/8 outside the ring) across `collar`'s; and declared `degenerate`, `step_touch`, a ball
touching `step`'s round between its tangent circle and its corner
(`Degenerate("a torus tangent to the other input's surface (S9d.4b.2)")`).

The reference is S9e.4b.4d.2c.1a's: S9e.3a's chained reference with each
body one input turned from its model profile (`TurnedInput`, its faces the
profile's own: each plane's annulus or disc, each cylinder's and cone's
band, each arc's band of its whole torus's tube, a corner's height exact to
the working precision; its membership the profile's radius at each
height), with S9e.4b.3c.3a's checks and each body's closed form by Pappus
over its model profile (its lines' frusta and annuli, each arc's band by
its angles). With `--check`, where the bodies' files exist, each is read
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
import generate_turned_fillets_boolean_fixtures as fillets
import generate_turned_cone_fillets_boolean_fixtures as cones

base = fillets.base
forms = fillets.forms
profiles = fillets.profiles
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS = base.OPERATIONS
CLASSES = ('pieces',)
PREFIX = 'boolean-turned-corners'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2c.1c'
sphere, construction = base.sphere, base.construction
ref = base.ref
q = forms.q
cylinder, block = forms.cylinder, forms.block
Lattice, at, TURN = profiles.Lattice, profiles.at, profiles.TURN
L, arc = fillets.L, fillets.arc
direction, half_tangent, power_above = cones.direction, cones.half_tangent, cones.power_above


# ------------------------------------------------------------------ numbers

class Surd:
    """`a + b sqrt(d)` exactly (Fractions, `d` positive): a corner's
    height."""

    def __init__(self, a, b, d):
        self.a, self.b, self.d = F(a), F(b), F(d)
        assert self.d > 0

    def mp(self):
        return q(self.a)+q(self.b)*mp.sqrt(q(self.d))

    def __float__(self):
        with mp.workdps(60):
            return float(self.mp())

    def __repr__(self):
        return f'Surd({self.a}, {self.b}, {self.d})'

    def __eq__(self, other):
        return isinstance(other, Surd) and (self.a, self.b, self.d) == (other.a, other.b, other.d)

    def __hash__(self):
        return hash((self.a, self.b, self.d))


def val(x):
    """A number as an mpf at the working precision: a Fraction, a Surd, a
    binary64 number's text or an mpf."""
    if isinstance(x, mp.mpf):
        return x
    if isinstance(x, Surd):
        return x.mp()
    if isinstance(x, str):
        return mp.mpf(x)
    return q(x)


def rounded(x):
    """A model number rounded once to binary64, as a Fraction."""
    return F(float(x)) if isinstance(x, Surd) else F(float(F(x)))


def num(x):
    """A coordinate's text: a number rounded once to binary64 (the shortest
    round-trip text)."""
    if isinstance(x, str):
        return x
    x = rounded(x)
    return str(x.numerator) if x.denominator == 1 else repr(float(x))


def on_grid(x, g):
    """`x` (a Fraction) on the multiples of `g`, the nearest (a half step
    refused)."""
    n = F(x)/g
    k = math.floor(n)
    # Within a millionth of a step of a half step refused: the kernel's
    # binary64 node may round it the other way.
    assert abs(n-k-F(1, 2)) > F(1, 10**6), 'a half step of the grid'
    return (k+1 if n-k > F(1, 2) else k)*g


def less(a, b):
    """`a < b` for Fractions and Surds (at 60 digits: the numbers compared
    differ by far more)."""
    with mp.workdps(60):
        return val(a) < val(b)


# ------------------------------------------------------------------ the model

def offset_square(w, cw):
    """`(w - cw)^2` exactly for a number or a Surd `cw + b sqrt(d)`."""
    if isinstance(w, Surd):
        assert w.a == cw, 'a corner not level with its centre plus a surd'
        return w.b*w.b*w.d
    return (F(w)-cw)**2


class Model:
    """A model profile with arcs: each arc a round or a fillet tangent to
    both lines (2c.1a's and 2c.1b's, `kind` 'tangent') or a round tangent to
    the wider of two cylinders and meeting the narrower at a corner
    (`kind` 'corner'); its sharp profile's lattice (each corner's round a
    step, each tangent arc its lines' corner) and each corner's piece."""

    def __init__(self, profile):
        self.profile = profile
        self.pts = [(F(r), w if isinstance(w, Surd) else F(w)) for r, w, _ in profile]
        self.segs = [s for _, _, s in profile]
        pts, segs = self.pts, self.segs
        assert pts[0][0] == 0 and pts[-1][0] == 0, 'a profile from the axis to the axis'
        for (_, wa), (_, wb) in zip(pts, pts[1:]):
            assert not less(wb, wa), 'heights never falling'
        self.corners = []
        for i, seg in enumerate(segs):
            if seg[0] != 'A':
                continue
            _, cr, cw, t, ccw = seg
            a, b = pts[i], pts[i+1]
            for e in (a, b):
                assert (e[0]-cr)**2+offset_square(e[1], cw) == t*t, ('an end off its arc', e)
            lines = (pts[i-1], a), (b, pts[i+2])
            kinds = ['cylinder' if ln[0][0] == ln[1][0] else 'plane' if ln[0][1] == ln[1][1] else 'cone'
                     for ln in lines]
            tangent = []
            for (p, s), e, k in zip(lines, (a, b), kinds):
                if k == 'cylinder':
                    tangent.append(e[1] == cw and abs(e[0]-cr) == t)
                elif k == 'plane':
                    tangent.append(e[0] == cr and abs(F(e[1])-cw) == t)
                else:
                    tangent.append(cones.foot((cr, cw), p, s) == e)
            c = {'centre': (cr, cw), 't': t, 'round': ccw, 'at': i, 'ends': (a, b), 'kinds': kinds}
            if all(tangent):
                c['kind'] = 'tangent'
            else:
                # S9e.4b.4d.2c.1c: a round between two cylinders, tangent to
                # the wider at one end, the narrower crossing the tube.
                assert ccw and kinds == ['cylinder', 'cylinder'] and sum(tangent) == 1, ('a corner', i)
                x = tangent.index(True)
                ry = (a, b)[1-x][0]
                rx = (a, b)[x][0]
                assert rx-t == cr and 0 < ry-cr < t, ('a corner of the wrong form', i)
                c.update(kind='corner', tangent=x, rx=rx, ry=ry, rise=t*t-(ry-cr)**2, up=(x == 0))
            self.corners.append(c)
        self.sharp = self.sharp_profile()
        self.lattice = Lattice(self.sharp)
        for c in self.corners:
            if c['kind'] == 'corner':
                c['box'] = self.corner_box(c)

    def step(self, c):
        """A corner's step: halfway from the corner to the nearer of the
        tube's top (bottom) and the narrower cylinder's far end, on the
        multiples of `2^-24` of the tube's power of two (the corner's height
        rounded once first, as the kernel's node), checked exactly past the
        corner and short of both."""
        i, (cr, cw), t, up = c['at'], c['centre'], c['t'], c['up']
        if up:
            corner, far, top = self.pts[i+1][1], self.pts[i+2][1], cw+t
            limit = min(top, F(far))
        else:
            corner, far, top = self.pts[i][1], self.pts[i-1][1], cw-t
            limit = max(top, F(far))
        g = power_above(t)/2**24
        mid = F(float((rounded(corner)+limit)/2))
        h = on_grid(mid, g)
        d = h-cw
        assert (d > 0) == up and d*d > c['rise'], 'a step not past its corner'
        assert (h < limit) == up and h != limit, 'a step not short of the tube and the line'
        return h

    def sharp_profile(self):
        """Each tangent arc its lines' corner, each corner's round a step
        between its two cylinders."""
        pts, segs = self.pts, self.segs
        out = []
        i = 0
        while i < len(pts):
            c = next((c for c in self.corners if c['at'] == i), None)
            if c is None:
                out.append(pts[i])
                i += 1
                continue
            if c['kind'] == 'corner':
                h = self.step(c)
                c['step'] = h
                a, b = c['ends']
                out += [(a[0], h), (b[0], h)]
            else:
                a0, a1, b0, b1 = pts[i-1], pts[i], pts[i+1], pts[i+2]
                # A cylinder's far end at a corner's surd height: another
                # point of its line.
                if isinstance(a0[1], Surd):
                    a0 = (a0[0], a1[1]-1)
                if isinstance(b1[1], Surd):
                    b1 = (b1[0], b0[1]+1)
                c['corner'] = fillets.meet(a0, a1, b0, b1)
                out.append(c['corner'])
            i += 2
        assert not any(isinstance(w, Surd) for _, w in out)
        return out

    def corner_box(self, c):
        """A corner's piece's box `[(r, w) low, (r, w) high]`: from the
        centre past `(r_X, H)` by half the tube."""
        (cr, cw), t, h, rx = c['centre'], c['t'], c['step'], c['rx']
        d = t/2
        far = (rx+d, h+d if c['up'] else h-d)
        return (cr, min(cw, far[1])), (far[0], max(cw, far[1]))

    def radius(self, w):
        """The profile's radius at a height (an mpf; at a node between two
        curved faces both's; None outside it)."""
        pts, segs = self.pts, self.segs
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            lo, hi = val(wa), val(wb)
            if wa == wb or not lo <= w <= hi:
                continue
            seg = segs[i]
            if seg[0] == 'A':
                _, cr, cw, t, _ = seg
                side = 1 if val(ra)+val(rb) > 2*val(cr) else -1
                return val(cr)+side*mp.sqrt(val(t)**2-(w-val(cw))**2)
            return val(ra)+(val(rb)-val(ra))*(w-lo)/(hi-lo)
        return None

    def inside(self, r, w):
        R = self.radius(w)
        return R is not None and r < R

    def arcs(self):
        """Each arc's `(R, W, T, th0, th1)` (its angles in the meridian from
        its start to its end, the short way it turns)."""
        out = {}
        pts, segs = self.pts, self.segs
        for i, seg in enumerate(segs):
            if seg[0] != 'A':
                continue
            _, cr, cw, t, ccw = seg
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            R, W, T = val(cr), val(cw), val(t)
            th0 = mp.atan2(val(wa)-W, val(ra)-R)
            th1 = mp.atan2(val(wb)-W, val(rb)-R)
            sweep = th1-th0
            if ccw and sweep < 0:
                sweep += 2*mp.pi
            if not ccw and sweep > 0:
                sweep -= 2*mp.pi
            out[i] = (R, W, T, th0, th0+sweep)
        return out

    def closed(self):
        """Pappus: `(volume, area)` of the profile turned about the axis,
        each arc's band by its angles."""
        vol, area = mp.mpf(0), mp.mpf(0)
        pts, segs = self.pts, self.segs
        arcs = self.arcs()
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if i in arcs:
                R, W, T, th0, th1 = arcs[i]

                def prim(th):
                    c, s = mp.cos(th), mp.sin(th)
                    return T*(R*R*s/2+R*T*(th/2+s*c/2)+T*T*(s-s**3/3)/2)
                vol += prim(th1)-prim(th0)
                area += abs(T*(R*(th1-th0)+T*(mp.sin(th1)-mp.sin(th0))))
            elif wb != wa:
                Ra, Rb, h = val(ra), val(rb), val(wb)-val(wa)
                vol += h*(Ra*Ra+Ra*Rb+Rb*Rb)/6
                area += (Ra+Rb)/2*mp.sqrt((Rb-Ra)**2+h*h)
            else:
                area += abs(val(rb)**2-val(ra)**2)/2
        return 2*mp.pi*vol, 2*mp.pi*area

    def moment(self):
        """The profile's `int r w dA = int r^2 w / 2 dw` along its boundary
        (counter-clockwise): a line's exactly (Simpson's rule, exact for its
        cubic), an arc's by quadrature."""
        out = mp.mpf(0)
        pts = self.pts
        arcs = self.arcs()
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if i in arcs:
                R, W, T, th0, th1 = arcs[i]
                out += mp.quad(lambda th: (R+T*mp.cos(th))**2*(W+T*mp.sin(th))/2*T*mp.cos(th), [th0, th1])
            elif wb != wa:
                Ra, Rb, Wa, Wb = val(ra), val(rb), val(wa), val(wb)
                f = lambda s: (Ra+(Rb-Ra)*s)**2*(Wa+(Wb-Wa)*s)/2
                out += (Wb-Wa)*(f(0)+4*f(mp.mpf(1)/2)+f(1))/6
        return out

    def construction_radius(self, w):
        """The kernel's construction's radius at a height (an mpf): the sharp
        profile's lattice's form, each corner's wider part less its piece
        (the box from the centre past the step, less the tube): that part's
        radius where the box holds the height the tube's outer radius there
        (the centre's radius off the tube)."""
        L = self.lattice
        cut = {}
        for c in self.corners:
            if c['kind'] != 'corner':
                continue
            j = next(j for j, p in enumerate(L.parts) if p['k'] == 0 and p['r0'] == c['rx'])
            cut[j] = c
        best = None
        for t in L.terms:
            vals = []
            for j in t:
                p = L.parts[j]
                lo, hi = (val(x) for x in p['range'])
                if not lo <= w <= hi:
                    vals = None
                    break
                v = val(p['r0'])+val(p['k'])*w
                if j in cut:
                    c = cut[j]
                    (blo, bhi) = c['box']
                    if val(blo[1]) < w < val(bhi[1]):
                        (cr, cw), T = c['centre'], val(c['t'])
                        dz = w-val(cw)
                        outer = val(cr)+(mp.sqrt(T*T-dz*dz) if abs(dz) < T else 0)
                        v = min(v, outer)
                vals.append(v)
            if vals:
                m = min(vals)
                best = m if best is None else max(best, m)
        return best


# ------------------------------------------------------------------ bodies

class Body:
    """An imported turned body about the world's `z` (or `x`) through
    `origin`: OCCT's profile and the model's, point by point within `2e-5`
    (the shaft's second round's edge `9 2^-19` below its centre) and of the
    same segments' kinds."""

    def __init__(self, name, occt, model, tolerance=None, axis='z', origin=(0, 0, 0)):
        self.name, self.axis, self.origin = name, axis, tuple(F(c) for c in origin)
        self.path = f'imported/corner_{name}.brep'
        self.occt, self.tolerance = occt, tolerance
        self.model = Model(model)
        assert len(occt) == len(model), name
        for (ra, wa, sa), (rb, wb, sb) in zip(occt, model):
            assert abs(val(ra)-val(rb)) <= 2e-5 and abs(val(wa)-val(wb)) <= 2e-5, (name, ra, wa)
            assert sa[0] == sb[0], name
        self.lattice = self.model.lattice
        self.ref_specs = [('turned', name)]
        self.ref_expr = 0

    def place(self, r, w):
        o = self.origin
        if self.axis == 'z':
            return (o[0]+rounded(r), o[1], o[2]+rounded(w))
        return (o[0]+rounded(w), o[1], o[2]+rounded(r))

    def frame_at(self, w):
        o = self.origin
        h = F(0) if w is None else F(w)
        if self.axis == 'z':
            return at('XY', (o[0], o[1], o[2]+h))
        return at('SIDE', (o[0]+h, o[1], o[2]))

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
                # About the world's x the profile's plane's `(x, z)` is `(w,
                # r)`: an arc counter-clockwise in `(r, w)` turns clockwise
                # there.
                flag = ccw if self.axis == 'z' else not ccw
                wire.append(f'A {num(c[0])} {num(c[1])} {num(c[2])} {num(t)} {1 if flag else 0}')
        n = '0 0 1' if self.axis == 'z' else '1 0 0'
        rows = [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} {n} {TURN}']
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
            out.append('torus' if seg[0] == 'A' else cones.kind_of(a, b))
        return sorted(out)

    def frame_point(self, r, w):
        """The world point at `(r, w)` in the profile's plane (mpf)."""
        o = [q(c) for c in self.origin]
        if self.axis == 'z':
            return (o[0]+r, o[1], o[2]+w)
        return (o[0]+w, o[1], o[2]+r)

    def chart(self, X):
        """A world point's `(r, w)` about the axis (mpf)."""
        o = [q(c) for c in self.origin]
        d = [X[i]-o[i] for i in range(3)]
        if self.axis == 'z':
            return mp.sqrt(d[0]**2+d[1]**2), d[2]
        return mp.sqrt(d[1]**2+d[2]**2), d[0]


def corner_height(cw, cr, t, ry, up):
    """The tube's crossing of the cylinder of radius `ry`: `cw +- sqrt(t^2 -
    (ry - cr)^2)` (a Fraction where it is a square)."""
    rise = t*t-(ry-cr)**2
    s = F(1 if up else -1)
    num_, den = rise.numerator, rise.denominator
    rn, rd = math.isqrt(num_), math.isqrt(den)
    if rn*rn == num_ and rd*rd == den:
        return cw+s*F(rn, rd)
    return Surd(cw, s, rise)


def grid30(x):
    """A number of `OCC496b`'s on the multiples of `2^-30`."""
    return on_grid(F(x), F(1, 2**30))


# `OCC496b`'s numbers (its file's), on the multiples of 2^-30.
OY = grid30(-4.7683715820310001e-07)
OZ = grid30(38.500011444089999)
W0 = grid30(-38.250450134280001)
R15 = grid30(15.000000476837158)
C1 = grid30(-33.750469207759998)
R14 = grid30(13.999999523167158)
C2 = grid30(-31.35044670105)
E2 = C2-F(9, 2**19)
R13 = grid30(13.000000476837158)
C3 = grid30(-20.23610687256)
V3 = '0.4432409130372571'
WP = grid30(-3.0000028610229998)
RT = grid30(5.2986021041871583)


def shaft():
    """`OCC496b`'s profile (OCCT's rows) and the model's: OCCT's second round
    from its stored edge `9 2^-19` below its centre, its third round to
    `OCC496b`'s stored angle and a cone from there to the top circle; the
    model's second round tangent at its centre's height, its third turned to
    the direction whose half-angle tangent is OCCT's cone's normal's on the
    multiples of 2^-40 (2c.1b's snap), its cone that direction's tangent."""
    corner1 = corner_height(C1, R15-2, F(2), R14, True)
    corner2 = corner_height(C2, R14-2, F(2), R13, True)
    c3 = (R13-F(9, 2), C3)
    with mp.workdps(60):
        v = mp.mpf(V3)
        p3 = (F(float(q(c3[0])+F(9, 2)*mp.cos(v))), F(float(q(c3[1])+F(9, 2)*mp.sin(v))))
        # The stored cone's normal toward it from the centre: normal to
        # OCCT's line from `p3` to the top circle.
        dr, dw = q(RT-p3[0]), q(WP-p3[1])
        n = (dw, -dr)
        length = mp.sqrt(n[0]**2+n[1]**2)
        n = (n[0]/length, n[1]/length)
        if n[0]*(q(p3[0])-q(c3[0]))+n[1]*(q(p3[1])-q(c3[1])) < 0:
            n = (-n[0], -n[1])
        s = half_tangent(n)
        scaled = s*2**40
        nearest = mp.nint(scaled)
        assert abs(scaled-nearest) < mp.mpf(1)/4, 'a direction near a half step of the grid'
        s = F(int(nearest), 2**40)
    e = direction(s)
    t3 = (c3[0]+F(9, 2)*e[0], c3[1]+F(9, 2)*e[1])
    top = t3[0]-e[1]/e[0]*(WP-t3[1])
    occt = [(0, W0, L), (R15, W0, L), (R15, C1, arc(R15-2, C1, 2, True)), (R14, rounded(corner1), L),
            (R14, E2, arc(R14-2, C2, 2, True)), (R13, rounded(corner2), L), (R13, C3, arc(*c3, F(9, 2), True)),
            (p3[0], p3[1], L), (RT, WP, L), (0, WP, L)]
    model = [(0, W0, L), (R15, W0, L), (R15, C1, arc(R15-2, C1, 2, True)), (R14, corner1, L),
             (R14, C2, arc(R14-2, C2, 2, True)), (R13, corner2, L), (R13, C3, arc(*c3, F(9, 2), True)),
             (t3[0], t3[1], L), (top, WP, L), (0, WP, L)]
    return occt, model


def as_occt(model):
    """A model profile's points rounded once to binary64 (OCCT's)."""
    return [(rounded(r), rounded(w), seg) for r, w, seg in model]


def make_bodies():
    corner = corner_height(F(3), F(3), F(1), F(7, 2), True)
    step = [(0, 0, L), (4, 0, L), (4, 3, arc(3, 3, 1, True)), (F(7, 2), corner, L), (F(7, 2), 6, L), (0, 6, L)]
    collar = [(0, 0, L), (F(7, 2), 0, L), (F(7, 2), 4, arc(F(11, 4), 5, F(5, 4), True)), (4, 5, L), (4, 8, L),
              (0, 8, L)]
    occt, model = shaft()
    out = [
        Body('shaft', occt, model, tolerance='0.003', axis='x', origin=(0, OY, OZ)),
        Body('step', as_occt(step), step),
        Body('collar', as_occt(collar), collar),
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


# The shaft's first corner ring: radius R14 at about -32.018 about the
# world's x through (0, OY, OZ). The box across it (its sides y = +-1, its
# bottom z = 51 below the ring's top at OZ + R14), the second round's
# tangent circle at C2 inside it too.
SHAFT_BOX = block(at('XY', (0, 0, 0)), -33.0, -1.0, -31.0, 1.0, 51.0, 55.0)
# A rod along the world's y about (x, z) = (-32, 419/8) of radius 1/4: the
# ring's top inside it, the ring leaving it on either side.
SHAFT_ROD = cylinder(at('ALONGY', (0, 0, 0)), 52.375, -32.0, 0.25, -3.0, 3.0)
# A cube of side 1 on `SKEW` about (-237/8, 0, 103/2), the second corner
# ring's point (radius 13 + 2^-21 at about -29.618) within it.
SHAFT_SKEW = block(at('SKEW', (-29.625, 0, 51.5)), -0.5, -0.5, 0.5, 0.5, -0.5, 0.5)
# `step`'s corner ring: radius 7/2 at 3 + sqrt(3)/2. A box across it (its
# sides y = +-1/2), a rod along x about (0, 31/8) of radius 1/4, a cube of
# side 1 on `SKEW` about (7/2, 0, 4).
STEP_BOX = block(at('XY', (0, 0, 0)), 3.0, -0.5, 4.5, 0.5, 3.25, 4.5)
STEP_ROD = cylinder(at('SIDE', (-10, 0, 0)), 0.0, 3.875, 0.25, 0.0, 20.0)
STEP_SKEW = block(at('SKEW', (3.5, 0, 4)), -0.5, -0.5, 0.5, 0.5, -0.5, 0.5)
# `collar`'s corner ring: radius 7/2 at 4. A box across it, and a ball of
# radius 5/16 about (29/8, 0, 4), 1/8 off the ring's point (7/2, 0, 4) (its
# poles off the narrower cylinder's line, the ball clear of the collar's
# cylinder of radius 4).
COLLAR_BOX = block(at('XY', (0, 0, 0)), 3.0, -0.5, 4.5, 0.5, 3.0, 4.5)
COLLAR_BALL = sphere(0.3125, at('XY', (3.625, 0, 4)))
# A ball touching `step`'s round at (19/5, 0, 18/5): its centre (3, 3) plus 1
# (4/5, 3/5), between the tangent circle (angle 0) and the corner (a sixth of
# a turn); the ball's centre 5/4 along that direction.
STEP_TOUCH = sphere(0.25, at('XY', (4, 0, 3.75)))
TOUCH = "a torus tangent to the other input's surface (S9d.4b.2)"


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out = []
    out += group('shaft_box', [imported('shaft'), SHAFT_BOX], three)
    out += group('shaft_rod', [imported('shaft'), SHAFT_ROD], three)
    out += group('shaft_skew', [imported('shaft'), SHAFT_SKEW], three)
    out += group('step_rod', [imported('step'), STEP_ROD], three)
    out += group('step_box', [imported('step'), STEP_BOX], three)
    out += group('step_skew', [imported('step'), STEP_SKEW], three)
    out += group('collar_box', [imported('collar'), COLLAR_BOX], three)
    out += group('collar_ball', [imported('collar'), COLLAR_BALL], three)
    out += group('step_touch', [imported('step'), STEP_TOUCH], {'fuse': 2, 'cut': 1, 'common': 0},
                 kind='degenerate', reason=TOUCH)
    return out


def all_cases():
    return cases()


validate = fillets.validate


# ------------------------------------------------------------------ the reference's turned input

tc = ref.tc
Z = mp.mpf(0)


class TurnedInput(cones.TurnedInput):
    """S9e.4b.4d.2c.1b's turned input over a profile whose corners' heights
    may be surds: each face's numbers at the working precision, its
    membership the profile's radius at each height."""

    def __init__(self, body):
        self.body = body
        self.model = m = body.model
        self.fr = tc.Frame(body.frame_at(None))
        self.surfs, self.faces = [], []
        pts, segs = m.pts, m.segs
        for i in range(len(pts)-1):
            (ra, wa), (rb, wb) = pts[i], pts[i+1]
            if segs[i][0] == 'A':
                self.arc_face_c(i, segs[i], (ra, wa), (rb, wb))
            elif wa == wb:
                self.plane_face(i, wa, min(ra, rb), max(ra, rb), 1 if rb < ra else -1)
            elif ra == rb:
                self.cylinder_face_c(i, ra, wa, wb, 1 if less(wa, wb) else -1)
            else:
                self.cone_face(i, (ra, wa), (rb, wb), 1 if wb > wa else -1)

    def cylinder_face_c(self, i, r, w0, w1, out):
        """The cylinder of radius `r` over heights `w0` to `w1` (either way;
        a surd's at the working precision)."""
        fr = self.fr
        R = val(r)
        lo, hi = sorted((val(w0), val(w1)))

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

    def arc_face_c(self, i, seg, a, b):
        """The band of the whole torus about the arc's centre over the arc's
        angles in the meridian (`th0` to `th1` the way it turns), its
        outward normal the tube's for a round, inward for a fillet."""
        fr = self.fr
        _, cr, cw, t, ccw = seg
        R, W, T, th0, th1 = self.model.arcs()[i]
        lo, hi = (th0, th1) if ccw else (th1, th0)
        assert hi-lo < mp.pi and hi > lo, 'an arc of half a turn or more'
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

    def box(self):
        m = self.model
        r = float(max(x for x, _ in m.pts))
        ws = [float(w) for _, w in m.pts]
        return tc.frame_box(self.fr, [(u, v, w) for u in (-r, r) for v in (-r, r) for w in (min(ws), max(ws))])


class FloatTurned(tc.FloatModel):
    """`TurnedInput` in binary64 for the rays and Monte Carlo: its surfaces'
    polynomials along a ray and its membership (the profile's radius at a
    height, `math.sqrt` alone)."""

    def __init__(self, S):
        self.S = S
        fr = S.fr
        self.o = [float(v) for v in fr.om]
        self.inv = [[float(v) for v in row] for row in fr.invm]
        self.kind = 'turned'
        m = S.model
        self.pts = [(float(r), float(w)) for r, w in m.pts]
        self.segs = m.segs

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
        for i in range(len(self.pts)-1):
            (ra, wa), (rb, wb) = self.pts[i], self.pts[i+1]
            if wa == wb or not wa <= w <= wb:
                continue
            seg = self.segs[i]
            if seg[0] == 'A':
                _, cr, cw, t, _ = seg
                cr, cw, t = float(cr), float(cw), float(t)
                side = 1.0 if ra+rb > 2*cr else -1.0
                return r < cr+side*math.sqrt(max(0.0, t*t-(w-cw)*(w-cw)))
            return r < ra+(rb-ra)*(w-wa)/(wb-wa)
        return False


# 2c.1b's hooks (its bodies' and 2c.1a's inputs; each module's own bodies
# its own inputs, all loaded in one process).
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
    # round at a point).
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
    cache = os.environ.get('TURNED_CORNERS_CACHE')
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
# Results smaller than Monte Carlo's resolution over the case's box (their
# counts and Monte Carlo unchecked): the shaft's rod's and skew cube's
# commons (about 0.8 and 0.65 in a box of some 30,000, a few points of the
# 100,000).
THIN = {('shaft_rod', 'common'), ('shaft_skew', 'common')}
# Groups whose events' spacing is not checked (their sines and gaps are): a
# partner across a corner ring (the round's torus and the narrower
# cylinder meeting there) or a tangent circle, where the reference finds
# each point of the ring as both surfaces' event, a rounding apart.
SPACING_UNCHECKED = {'shaft_box', 'shaft_rod', 'step_rod'}


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
    """Every case's native rows (`compare_turned_corners_boolean.py`)."""
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
    r, w = body.chart(X)
    best = None
    pts = [(F(a), F(b)) for a, b, _ in body.occt]
    for ((ra, wa), (rb, wb)), (_, _, seg) in zip(zip(pts, pts[1:]), body.occt):
        if seg[0] == 'A':
            _, cr, cw, t, _ = seg
            dist = abs(mp.sqrt((r-q(cr))**2+(w-q(cw))**2)-val(t))
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
        os.environ['TURNED_CORNERS_CACHE'] = args.cache
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
