#!/usr/bin/env python3
"""Fixtures for S9d.4b.1 of REVIEW_NOTES.md: Booleans of a torus v-segment or
wedge (`Solid::torus_with` other than a whole torus) against a polyhedral
prism, before any of its kernel code.

`boolean-torus-segment-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a torus's block is its `frame` and
`torus MAJOR MINOR LOW HIGH ANGLE` rows); `boolean-torus-segment-expected.tsv`
gives per case `expect KIND S9d.4b.1` (`solid`, `empty` or `degenerate`,
then `reason TEXT`) and `result N volume area cx cy cz` or `empty` from
`torus_segment_boolean_reference.py`, as S9d.4a's; `boolean-torus-segment-
frames.tsv` the stored axes' bits. Frames are the curved generator's.

Parts, all of radii 5/2 and 3/2 on `XY` at the origin unless said (S3's
construction: a segment the region between the tube's arc and the axis,
revolved, its ends planar discs normal to the axis; a wedge the whole tube
over a partial turn, its ends the tube's discs): `OH` the outer half
(latitudes `-pi/2..pi/2`: a barrel of radius up to 4, its end discs of
radius 5/2 at `w = +-3/2` tangent to the wall along their rings), `IH` the
inner half (`pi/2..3pi/2`: a spool, its waist of radius 1, inside out as
OCCT builds it), `UB` an upper band (`0.5..2.25`: from the axis to the
tube's outer side between its end heights 0.719 and 1.167, the tube's cap
above; its end discs of radii 3.816 below and 1.557 above), `HW` a half
turn and `QW` a quarter turn of the whole tube (their ends in the half-
planes of `x` and of the rounded `-x`, `y`). Pairs: a bar through the axis
across `OH` (`oh_bar`), a box across its top end disc and wall
(`oh_side`); a slab through the spool's waist (`ih_slab`: two solids), a
column through its waist, its corners out of the spool there (`ih_column`,
and `column_ih` the column first: six solids); a box across `UB`'s lower
end disc and wall (`ub_side`), a box standing on its upper end disc
(`ub_lid`: coplanar, opposite), a box through its axis (`ub_hole`); a bar
along both of `HW`'s end discs (`hw_bar`: its common two solids), a slab
across it (`slab_hw`, the slab first); a box across `QW`'s start disc
(`qw_box`), a box in the hole cutting a cap from the tube's inside over
both ends (`qw_hole`), a box beyond its turn (`qw_far`: the common empty),
a box against its start disc (`qw_back`: coplanar, opposite); turned
frames: `OH` in `TILT` above a plane (`oh_tilt`), `HW` in `TILT2` across a
slab (`hw_tilt`); `degenerate` with reasons: a face on `OH`'s end plane
tangent to its wall along the ring (`oh_lid`), a face tangent to `UB`'s wall
along its top circle (`ub_top`), a face tangent to `QW` at its end circle's
outermost point (`qw_touch`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms, in the part's ideal frame by one quadrature (mpmath's
  tanh-sinh between the kinks) of exact 2D forms along the axis, each
  slice's section as the segment's shape gives it by hand (the outer half a
  disc to `R + q`, the inner half to `R - q`, the band a disc between its
  end heights and the tube's annulus above, a wedge the annulus) and a
  wedge's box clipped to its sector (the clipped faces its end discs): an
  aligned box by S9d.1's rectangle-in-disc antiderivatives (the wall by its
  latitudes, `r rho` times the circle's angle in the rectangle; the end
  discs by rectangles in their discs; the box's faces by chords), a
  half-space by circular segments, inputs meeting on a face or at a point
  by their sums; within 1e-30 in exact frames, 1e-15 in turned ones and
  where a wedge's rounded end (off its ideal plane by 1.2e-16 or 6.1e-17)
  cuts the box;
* the operations by inclusion and exclusion in the meridian sweep, every
  operation's volume and moments two ways (normal slices clipping the
  common, meridian half-planes classifying boundaries), both inputs'
  slicings against their closed forms (Pappus), the wall two ways, the end
  discs two ways (exactly at their planes, and by the meridians' radii or
  the normal slices' chords), the area identity `area(fuse) + area(common)
  + 2 opp = area(A) + area(B)`, every face's classes summing to its area,
  within 1e-30; the fuse's solid count by the sweep as by the rule;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, slicing breakpoints of either
  direction closer than that) must find none but in the declared degenerate
  pairs; every other pair's margins (a vertex's distance from the part's
  faces and rims, a face's plane's from the wall's tangent planes whose
  points of contact lie on the face and the wall, the sine at which an edge
  crosses the wall or an end plane, a rim's crossing of a face or its gap)
  are reported.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import math
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from brep_reference import sin_rn
from identity_reference import Boundary, Case, encode_boolean_case
from curve_surface_reference import stored_axes
import torus_segment_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES
from generate_sphere_boolean_fixtures import FloatPrism, arc_in_rect, frame_name, ideal_axes, rect_disc
from generate_cone_boolean_fixtures import prism_measures, segment

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-torus-segment'
STEP = 'S9d.4b.1'
BOOLEAN_OPERATION = 93
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
TWO_PI = ref.TWO_PI
PI = 3.141592653589793
HALF_PI = 1.5707963267948966
THREE_HALVES_PI = 4.71238898038469


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def square(x0, y0, x1, y1):
    return Boundary(points=[(float(x0), float(y0)), (float(x1), float(y0)), (float(x1), float(y1)),
                            (float(x0), float(y1))])


def make(spec, op):
    if spec[0] == 'part':
        _, torus, frame, _ = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], torus=tuple(float(v) for v in torus))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def part(shape, frame, R=2.5, r=1.5):
    """A segment or wedge: `shape` ('outer_half' | 'inner_half' | ('band',
    low, high) | ('wedge', angle))."""
    if shape == 'outer_half':
        torus = (R, r, -HALF_PI, HALF_PI, TWO_PI)
    elif shape == 'inner_half':
        torus = (R, r, HALF_PI, THREE_HALVES_PI, TWO_PI)
    elif shape[0] == 'band':
        torus = (R, r, shape[1], shape[2], TWO_PI)
    else:
        torus = (R, r, 0.0, TWO_PI, shape[1])
    return ('part', torus, frame, shape)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def spec_frame(spec):
    return spec[2]


class Boolean:
    def __init__(self, name, operation, obj, tool, form, kind='solid', reason=None):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.form = form
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.frames = (spec_frame(obj), spec_frame(tool))
        self.shape = obj[3] if obj[0] == 'part' else tool[3]

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, form):
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, form, kind, reason))
    return out


LID_TANGENT = "a face of the prism on an end disc's plane, tangent to the wall along the disc's ring"
TOP_TANGENT = 'a face of the prism tangent to the wall along a circle'
POINT_TANGENT = "a face of the prism tangent to the wall at a point of an end disc's circle"

UB = ('band', 0.5, 2.25)
Z_HI_UB = 1.5*sin_rn(2.25)


def cases():
    """S9d.4b.1's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    XY = at('XY', O)
    OH, IH, UB_ = part('outer_half', XY), part('inner_half', XY), part(UB, XY)
    HW, QW = part(('wedge', PI), XY), part(('wedge', HALF_PI), XY)
    BIG = square(-5, -5, 5, 5)
    all3 = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
    box = ('box', None)
    out = []
    out += group('oh_bar', OH, prism([square(-0.5, -5, 0.5, 5)], XY, -1, 1), all3, box)
    out += group('oh_side', OH, prism([square(1, -6, 6, 6)], XY, 0.5, 3), {'common': 'solid'}, box)
    out += group('oh_lid', OH, prism([BIG], XY, 1.5, 3), {'fuse': ('degenerate', LID_TANGENT)},
                 ('touch', 'oh_lid'))
    out += group('ih_slab', IH, prism([BIG], XY, -0.5, 0.5), {'cut': 'solid'}, box)
    out += group('ih_column', IH, prism([square(-0.75, -0.75, 0.75, 0.75)], XY, -3, 3), {'common': 'solid'}, box)
    out += group('column_ih', prism([square(-0.75, -0.75, 0.75, 0.75)], XY, -3, 3), IH, {'cut': 'solid'}, box)
    out += group('ub_side', UB_, prism([square(2, -6, 6, 6)], XY, -1, 1), all3, box)
    out += group('ub_lid', UB_, prism([square(-2, -2, 2, 2)], XY, Z_HI_UB, 2), {'fuse': 'solid', 'common': 'solid'},
                 None)
    out += group('ub_hole', UB_, prism([square(-1.25, -1.25, 1.25, 1.25)], XY, 0, 2), {'cut': 'solid'}, box)
    out += group('ub_top', UB_, prism([BIG], XY, 1.5, 3), {'fuse': ('degenerate', TOP_TANGENT)}, ('touch', None))
    out += group('hw_bar', HW, prism([square(-5, -0.5, 5, 0.5)], XY, -2, 2), all3, ('box', 'rounded'))
    out += group('slab_hw', prism([BIG], XY, -0.5, 0.5), HW, {'cut': 'solid'}, ('box', 'rounded'))
    out += group('qw_box', QW, prism([square(1.75, -1, 5, 1.25)], XY, -2, 0.5), {'cut': 'solid', 'common': 'solid'},
                 box)
    out += group('qw_hole', QW, prism([square(-1.25, -1.25, 1.25, 1.25)], XY, -2, 2), {'common': 'solid'},
                 ('box', 'rounded'))
    out += group('qw_far', QW, prism([square(-5, -5, -1.5, 5)], XY, -2, 2), {'fuse': 'solid'}, box)
    out += group('qw_back', QW, prism([square(0.5, -3, 5, 0)], XY, -2, 2), {'fuse': 'solid', 'common': 'empty'},
                 ('touch', 'qw_back'))
    out += group('qw_touch', QW, prism([square(4, -3, 6, 3)], XY, -2, 2), {'fuse': ('degenerate', POINT_TANGENT)},
                 ('touch', None))
    out += group('oh_tilt', part('outer_half', at('TILT', O)), prism([square(-6, -6, 6, 6)], XY, 0.75, 5),
                 {'common': 'solid'}, ('planes', [('cap', 0)]))
    out += group('hw_tilt', part(('wedge', PI), at('TILT2', (0.25, 0.5, 0))),
                 prism([square(-6, -6, 6, 6)], XY, -0.25, 0.25), {'cut': 'solid', 'common': 'solid'}, None)
    return out


def exact_pair(case):
    return all(frame_name(f) in EXACT for f in case.frames)


# ------------------------------------------------------------------ closed forms

class Shape:
    """A part's sections by hand, in its ideal frame: `structure` [(s0, s1,
    kind)] (`out` the disc to `R + q`, `in` to `R - q`, `annulus`), the wall's
    latitudes, its end discs (height, radius) or a wedge's turn."""

    def __init__(self, shape, torus):
        R, r = mp.mpf(torus[0]), mp.mpf(torus[1])
        self.R, self.r = R, r
        self.wedge = shape[0] == 'wedge'
        if shape == 'outer_half':
            self.structure = [(-r, r, 'out')]
            self.lats = (-mp.pi/2, mp.pi/2)
            self.ends = [(-r, R), (r, R)]
        elif shape == 'inner_half':
            self.structure = [(-r, r, 'in')]
            self.lats = (mp.pi/2, 3*mp.pi/2)
            self.ends = [(r, R), (-r, R)]
        elif shape[0] == 'band':
            # Between its end heights the band reaches the axis and the
            # tube's outer side; above the upper one it is the tube's cap.
            zl, zh = mp.mpf(torus[1]*sin_rn(shape[1])), mp.mpf(torus[1]*sin_rn(shape[2]))
            assert 0 < zl < zh < r
            ql, qh = mp.sqrt(r*r-zl*zl), mp.sqrt(r*r-zh*zh)
            self.structure = [(zl, zh, 'out'), (zh, r, 'annulus')]
            self.lats = (mp.asin(zl/r), mp.pi-mp.asin(zh/r))
            self.ends = [(zl, R+ql), (zh, R-qh)]
        else:
            self.structure = [(-r, r, 'annulus')]
            self.lats = (mp.mpf(0), 2*mp.pi)
            self.ends = []
            c, s = mp.mpf(ref.cos_rn(torus[4])), mp.mpf(ref.sin_rn(torus[4]))
            U = mp.atan2(s, c)
            self.U = U if U > 0 else U+2*mp.pi

    def q(self, s):
        return mp.sqrt(max(self.r*self.r-s*s, mp.mpf(0)))

    def discs(self, s):
        for s0, s1, kind in self.structure:
            if s0 < s < s1:
                q = self.q(s)
                return {'out': [(self.R+q, 1)], 'in': [(self.R-q, 1)],
                        'annulus': [(self.R+q, 1), (self.R-q, -1)]}[kind]
        return []

    def kinks(self):
        return sorted(set([self.r, -self.r]+[x for s0, s1, _ in self.structure for x in (s0, s1)]))

    def measures(self):
        """Volume, moments (u, v, w), wall area and end areas (Pappus)."""
        R, r = self.R, self.r
        if self.wedge:
            U = self.U
            I2 = mp.pi*r*r*(R*R+r*r/4)
            return (U*mp.pi*r*r*R, (I2*mp.sin(U), I2*(1-mp.cos(U)), mp.mpf(0)), U*2*mp.pi*R*r,
                    2*mp.pi*r*r)
        ks = self.kinks()
        V = mp.quad(lambda s: sum(sg*mp.pi*rho*rho for rho, sg in self.discs(s)), ks)
        Mw = mp.quad(lambda s: s*sum(sg*mp.pi*rho*rho for rho, sg in self.discs(s)), ks)
        a, b = self.lats
        wall = 2*mp.pi*r*(R*(b-a)+r*(mp.sin(b)-mp.sin(a)))
        return V, (mp.mpf(0), mp.mpf(0), Mw), wall, sum(mp.pi*rho*rho for _, rho in self.ends)


def box_form(sh, box, clipped=()):
    """The part (ideal frame) and an aligned box `[u0, u1] x [v0, v1] x [w0,
    w1]` (a wedge's clipped to its sector): the common's volume and moments,
    the wall inside, the end discs inside (a wedge's: the clipped faces
    `clipped`, ('u', 0) or ('v', 0), inside the tube), the box's other faces
    inside."""
    R, r = sh.R, sh.r
    (u0, u1), (v0, v1), (w0, w1) = [tuple(mp.mpf(x) for x in rng) for rng in box]
    zero = {'volume': mp.mpf(0), 'moments': (mp.mpf(0),)*3, 'wall': mp.mpf(0), 'ends': mp.mpf(0),
            'faces': mp.mpf(0)}
    if u1 <= u0 or v1 <= v0 or w1 <= w0:
        return zero
    wa, wb = max(w0, -r), min(w1, r)
    if wb <= wa:
        return zero
    kinks = set([wa, wb]+[k for k in sh.kinks() if wa < k < wb])
    vals = [abs(u0), abs(u1), abs(v0), abs(v1)]+[mp.sqrt(x*x+y*y) for x in (u0, u1) for y in (v0, v1)]
    for val in vals:
        dq = abs(val-R)
        if dq <= r:
            for s in (-mp.sqrt(r*r-dq*dq), mp.sqrt(r*r-dq*dq)):
                if wa < s < wb:
                    kinks.add(s)
    pts = sorted(kinks)

    def rd(s):
        out = (mp.mpf(0),)*3
        for rho, sg in sh.discs(s):
            d = rect_disc(u0, u1, v0, v1, rho)
            out = tuple(x+sg*y for x, y in zip(out, d))
        return out
    V = mp.quad(lambda s: rd(s)[0], pts)
    Mu = mp.quad(lambda s: rd(s)[1], pts)
    Mv_ = mp.quad(lambda s: rd(s)[2], pts)
    Mw = mp.quad(lambda s: s*rd(s)[0], pts)
    out = {'volume': V, 'moments': (Mu, Mv_, Mw)}
    a, b = sh.lats
    # The wall by its latitudes phi (s = r sin phi), the element r rho dphi.
    phis = sorted(set([a, b]+[p for s in pts for p in (mp.asin(max(min(s/r, 1), -1)),
                                                       mp.pi-mp.asin(max(min(s/r, 1), -1)))
                              for p in (p, p+2*mp.pi, p-2*mp.pi) if a < p < b]))

    def wall(phi):
        s = r*mp.sin(phi)
        if not w0 < s < w1:
            return mp.mpf(0)
        rho = R+r*mp.cos(phi)
        return r*rho*arc_in_rect(u0, u1, v0, v1, rho)
    out['wall'] = mp.quad(wall, phis)

    def chord(c, lo, hi, s):
        total = mp.mpf(0)
        for rho, sg in sh.discs(s):
            Q2 = rho*rho-c*c
            if Q2 <= 0:
                continue
            Q = mp.sqrt(Q2)
            total += sg*max(mp.mpf(0), min(hi, Q)-max(lo, -Q))
        return total
    faces = {('u', 0): (u0, v0, v1), ('u', 1): (u1, v0, v1), ('v', 0): (v0, u0, u1), ('v', 1): (v1, u0, u1)}
    face_in = ends_in = mp.mpf(0)
    for key, (c, lo, hi) in faces.items():
        val = mp.quad(lambda s: chord(c, lo, hi, s), pts)
        if key in clipped:
            ends_in += val
        else:
            face_in += val
    for w in (w0, w1):
        if -r < w < r:
            face_in += rd(w)[0]
    for z, rho in sh.ends:
        if w0 < z < w1:
            ends_in += rect_disc(u0, u1, v0, v1, rho)[0]
    out['ends'], out['faces'] = ends_in, face_in
    return out


def halfspace_form(sh, a, b):
    """A segment (ideal frame) and the half-space `a . p >= b` (not normal to
    the axis): the common's volume and moments, the wall inside, the end
    discs inside, the plane's face inside."""
    R, r = sh.R, sh.r
    au, av, aw = (mp.mpf(c) for c in a)
    b = mp.mpf(b)
    nuv = mp.sqrt(au*au+av*av)
    mu, mv = au/nuv, av/nuv
    d = lambda s: (b-aw*s)/nuv
    kinks = list(sh.kinks())
    N = nuv*nuv
    D2 = [aw*aw, -2*b*aw, b*b]
    E = [D2[0]+N, D2[1], D2[2]+N*(R*R-r*r)]
    poly = [E[0]*E[0], 2*E[0]*E[1], E[1]*E[1]+2*E[0]*E[2], 2*E[1]*E[2], E[2]*E[2]]
    poly = [poly[0], poly[1], poly[2]-4*R*R*N*D2[0], poly[3]-4*R*R*N*D2[1], poly[4]-4*R*R*N*D2[2]]
    for s in mp.polyroots(poly, maxsteps=400, extraprec=400):
        s = mp.mpc(s)
        if abs(s.imag) < mp.mpf(10)**-15 and -r < s.real < r:
            kinks.append(s.real)
    kinks = sorted(kinks)

    def part(k):
        def g(s):
            out = mp.mpf(0)
            for rho, sg in sh.discs(s):
                out += sg*segment(rho, d(s))[k]
            return out
        return g
    V = mp.quad(part(0), kinks)
    Mn = mp.quad(part(1), kinks)
    Mw = mp.quad(lambda s: s*part(0)(s), kinks)
    lo, hi = sh.lats

    def wall(phi):
        s = r*mp.sin(phi)
        rho = R+r*mp.cos(phi)
        return r*rho*segment(rho, d(s))[2]
    phis = sorted(set([lo, hi]+[p for s in kinks for p0 in (mp.asin(max(min(s/r, 1), -1)),)
                                for p in (p0, mp.pi-p0, p0+2*mp.pi, p0-2*mp.pi) if lo < p < hi]))
    return {'volume': V, 'moments': (mu*Mn, mv*Mn, Mw), 'wall': mp.quad(wall, phis),
            'ends': sum((segment(rho, d(z))[0] for z, rho in sh.ends), mp.mpf(0)),
            'faces': mp.sqrt(au*au+av*av+aw*aw)/nuv*mp.quad(part(3), kinks)}


def touch_area(case, sh):
    """The area on which the inputs meet with opposite sides (by hand)."""
    if case.form[1] == 'oh_lid':
        return mp.pi*sh.R*sh.R
    if case.form[1] == 'qw_back':
        return mp.pi*sh.r*sh.r
    return mp.mpf(0)


def pair_forms(case):
    """{op: (volume, area, centre)} in closed form for a pair (see the
    module's docstring), in the part's ideal frame; None for a pair
    without one."""
    if case.form is None:
        return None
    K_role = 'A' if case.obj.torus is not None else 'B'
    kc, pc = (case.obj, case.tool) if K_role == 'A' else (case.tool, case.obj)
    sh = Shape(case.shape, kc.torus)
    X, Y, Nn = ideal_axes(frame_name(kc.frame))
    ko = tuple(mp.mpf(v) for v in kc.frame[:3])
    local = lambda P: tuple(ref.dot(ref.sub(P, ko), ax) for ax in (X, Y, Nn))
    world_m = lambda V, m: tuple(ko[j]*V+sum(m[i]*(X, Y, Nn)[i][j] for i in range(3)) for j in range(3))
    VK, mK, AW, AE = sh.measures()
    AK = AW+AE
    MK = world_m(VK, mK)
    VP, MP, AP = prism_measures(pc)
    px, py, pn = ideal_axes(frame_name(pc.frame))
    po = tuple(mp.mpf(v) for v in pc.frame[:3])
    pworld = lambda u, v, w: tuple(po[i]+u*px[i]+v*py[i]+w*pn[i] for i in range(3))
    kind = case.form[0]
    opp = mp.mpf(0)
    if kind == 'touch':
        Vc, Mc, K_in, P_in = mp.mpf(0), (mp.mpf(0),)*3, mp.mpf(0), mp.mpf(0)
        opp = touch_area(case, sh)
    elif kind == 'box':
        pts = pc.boundaries[0].points
        corners = [local(pworld(mp.mpf(u), mp.mpf(v), mp.mpf(w))) for u, v in pts for w in (pc.start, pc.end)]
        box = []
        for i in range(3):
            vals = sorted(c[i] for c in corners)
            assert vals[3]-vals[0] < mp.mpf(10)**-14 and vals[7]-vals[4] < mp.mpf(10)**-14, 'an aligned box'
            box.append([vals[0], vals[7]])
        clipped = []
        if sh.wedge:
            # The sector: v >= 0 (a half turn), u, v >= 0 (a quarter).
            assert abs(sh.U-mp.pi) < 1e-15 or abs(sh.U-mp.pi/2) < 1e-15
            for axis, key in ((1, 'v'),)+(((0, 'u'),) if sh.U < 2 else ()):
                if box[axis][0] < 0:
                    box[axis][0] = mp.mpf(0)
                    clipped.append((key, 0))
        f = box_form(sh, box, clipped)
        Vc, Mc, K_in, P_in = f['volume'], f['moments'], f['wall']+f['ends'], f['faces']
    else:
        # The part less the parts beyond each listed face.
        _, faces = case.form
        Vc, Mc, K_in, P_in = VK, mK, AK, mp.mpf(0)
        pts = [(mp.mpf(u), mp.mpf(v)) for u, v in pc.boundaries[0].points]
        for tag, idx in faces:
            if tag == 'wall':
                p, q = pts[idx], pts[(idx+1) % len(pts)]
                e = (q[0]-p[0], q[1]-p[1])
                a = tuple(e[1]*px[i]-e[0]*py[i] for i in range(3))
                P0 = pworld(p[0], p[1], mp.mpf(0))
            else:
                sg = -1 if idx == 0 else 1
                a = tuple(sg*c for c in pn)
                P0 = pworld(mp.mpf(0), mp.mpf(0), mp.mpf(pc.start if idx == 0 else pc.end))
            ac = tuple(ref.dot(a, ax) for ax in (X, Y, Nn))
            f = halfspace_form(sh, ac, ref.dot(a, P0)-ref.dot(a, ko))
            Vc -= f['volume']
            Mc = tuple(m-x for m, x in zip(Mc, f['moments']))
            K_in -= f['wall']+f['ends']
            P_in += f['faces']
    McW = world_m(Vc, Mc)
    parts = {'K': (VK, MK, AK, K_in), 'P': (VP, MP, AP, P_in)}
    first, second = ('K', 'P') if K_role == 'A' else ('P', 'K')
    VA, MA, AA, A_in = parts[first]
    VB, MB, AB, B_in = parts[second]
    ops = {'common': (Vc, McW, A_in+B_in),
           'fuse': (VA+VB-Vc, tuple(p+q-c for p, q, c in zip(MA, MB, McW)), AA-A_in+AB-B_in-2*opp),
           'cut': (VA-Vc, tuple(p-c for p, c in zip(MA, McW)), AA-A_in+B_in)}
    out = {}
    for op, (V, m, A) in ops.items():
        if V <= mp.mpf(10)**-30:
            out[op] = None
            continue
        out[op] = (V, A, tuple(v/V for v in m))
    return out


# ------------------------------------------------------------------ Monte Carlo

class FloatPart:
    def __init__(self, K):
        T = K.T
        self.K = K
        self.o = [float(v) for v in T.o]
        self.inv = [[float(c) for c in row] for row in T.inv]
        self.R, self.r = float(T.R), float(T.r)
        if K.kind == 'segment':
            self.vlo, self.vhi = float(K.vlo), float(K.vhi)
            self.zb = [float(z) for z in K.zband]
        else:
            self.cu, self.su = (float(v) for v in K.eu)
            self.big = float(K.Um) > math.pi
        x, y, n = T.x, T.y, T.n
        e = self.R+self.r
        pts = []
        for su in (-1, 1):
            for sv in (-1, 1):
                for sw in (-1, 1):
                    pts.append([self.o[i]+su*e*float(x[i])+sv*e*float(y[i])+sw*self.r*float(n[i]) for i in range(3)])
        self.box = ([min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)])

    def in_arc(self, v):
        k = math.ceil((self.vlo-v)/(2*math.pi))
        return v+2*math.pi*k <= self.vhi

    def contains(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        u, v, w = (sum(row[i]*d[i] for i in range(3)) for row in self.inv)
        t = math.hypot(u, v)
        R, r = self.R, self.r
        if abs(w) >= r:
            return False
        q = math.sqrt(r*r-w*w)
        if self.K.kind == 'wedge':
            if not (R-q < t < R+q):
                return False
            a, b = v > 0, u*self.su-v*self.cu > 0
            return (a or b) if self.big else (a and b)
        a = math.asin(w/r)
        n = 0
        for side, ang in ((1, a), (-1, math.pi-a)):
            if self.in_arc(ang) and R+side*q > t:
                n += 1
        return n % 2 == 1


def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fk, fp = FloatPart(pair.K), FloatPrism(pair.P)
    box0 = [min(fk.box[0][i], fp.box[0][i]) for i in range(3)]
    box1 = [max(fk.box[1][i], fp.box[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        s, p = fk.contains(X), fp.contains(X)
        a, b = (s, p) if pair.part_first else (p, s)
        for op in ref.OPS:
            if sets[op](a, b):
                e = acc[op]
                e[0] += 1
                for i in range(3):
                    e[1][i] += X[i]
                    e[2][i] += X[i]*X[i]
    out = {}
    for op, (k, s1, s2) in acc.items():
        p = k/n
        vol = vbox*p
        sv = vbox*math.sqrt(max(p*(1-p), 1e-300)/n)
        if k > 1:
            c = [s/k for s in s1]
            sc = [math.sqrt(max(s2[i]/k-c[i]*c[i], 0)/k) for i in range(3)]
        else:
            c, sc = None, None
        out[op] = (vol, sv, c, sc)
    return out


# ------------------------------------------------------------------ margins

def rims(K):
    """The part's rims: (centre, radius, unit normal, unit x, unit y) in the
    chart (mpf; orthonormal coordinates of an ideal frame)."""
    out = []
    if K.kind == 'segment':
        for e in K.ends:
            out.append(((mp.mpf(0), mp.mpf(0), e['zm']), e['rho'], (0, 0, 1), (1, 0, 0), (0, 1, 0)))
    else:
        for c, s in K.end_units():
            out.append(((K.Rm*c, K.Rm*s, mp.mpf(0)), K.rm, (-s, c, 0), (c, s, 0), (0, 0, 1)))
    return out


def on_polygon(X, loop, N):
    n = len(loop)
    sides = [ref.dot(ref.cross(ref.sub(loop[(k+1) % n], loop[k]), ref.sub(X, loop[k])), N) for k in range(n)]
    return all(v >= 0 for v in sides) or all(v <= 0 for v in sides)


def margins(pair):
    """The pair's distances from tangency and incidence, relative to its
    size (see the module's docstring)."""
    K, P, size = pair.K, pair.P, pair.size
    R, r = K.Rm, K.rm
    pc = {k: K.chart(X) for k, X in P.P.items()}
    pm = {k: ref.Mv(p) for k, p in pc.items()}

    def on_wall(X):
        phi = mp.atan2(X[2], mp.sqrt(X[0]**2+X[1]**2)-R)
        if K.kind == 'segment':
            return ref.in_range(phi, K.vlo, K.vhi)
        return K.in_turn(mp.atan2(X[1], X[0]))
    vertex = mp.inf
    for p in pm.values():
        rho = mp.sqrt(p[0]**2+p[1]**2)
        if on_wall(p):
            vertex = min(vertex, abs(mp.sqrt((rho-R)**2+p[2]**2)-r))
        for c, rad, nrm, ex, ey in rims(K):
            d = ref.sub(p, c)
            h = ref.dot(d, nrm)
            inplane = mp.sqrt(ref.dot(d, ex)**2+ref.dot(d, ey)**2)
            if K.kind == 'segment':
                if inplane <= rad:
                    vertex = min(vertex, abs(h))
            else:
                if ref.dot(p, (ex[0], ex[1], 0)) >= 0 and inplane <= rad:
                    vertex = min(vertex, abs(h))
            vertex = min(vertex, mp.sqrt((inplane-rad)**2+h*h))
    face = rim = mp.inf
    edge = mp.inf
    for piece in P.pieces:
        for f in piece.faces:
            if f['tag'][0] == 'internal':
                continue
            loop = [pm[k] for k in f['loop']]
            a = ref.cross(ref.sub(loop[1], loop[0]), ref.sub(loop[2], loop[0]))
            na = mp.sqrt(ref.dot(a, a))
            N = tuple(c/na for c in a)
            h = ref.dot(N, loop[0])
            sb = mp.sqrt(N[0]**2+N[1]**2)
            if sb < mp.mpf(10)**-30:
                # Normal to the axis: from the tube's top and bottom circles
                # (on the wall) and the end planes (unless on one).
                for c0 in (r, -r):
                    if on_wall((R, mp.mpf(0), c0)):
                        face = min(face, abs(h*N[2]-c0))
            else:
                g = mp.atan2(N[1], N[0])
                cb = N[2]
                for th in (g, g+mp.pi):
                    ct = mp.cos(th-g)
                    for sg in (1, -1):
                        cphi, sphi = sg*sb*ct, sg*cb
                        X = ((R+r*cphi)*mp.cos(th), (R+r*cphi)*mp.sin(th), r*sphi)
                        if not on_wall(X):
                            continue
                        dist = ref.dot(N, X)-h
                        Xp = tuple(X[i]-dist*N[i] for i in range(3))
                        if on_polygon(Xp, loop, N):
                            face = min(face, abs(dist))
            # The rims against the face: their crossings' sines, or gaps.
            for c, rad, nrm, ex, ey in rims(K):
                A1, A2 = rad*ref.dot(N, ex), rad*ref.dot(N, ey)
                A = mp.sqrt(A1*A1+A2*A2)
                d0 = ref.dot(N, c)-h
                if A < mp.mpf(10)**-30:
                    if abs(d0) > mp.mpf(10)**-30:
                        rim = min(rim, abs(d0))
                    continue
                g0 = mp.atan2(A2, A1)
                if abs(d0) >= A:
                    X = tuple(c[i]+rad*(mp.cos(g0)*ex[i]+mp.sin(g0)*ey[i])*(-1 if d0 > 0 else 1) for i in range(3))
                    Xp = tuple(X[i]-(ref.dot(N, X)-h)*N[i] for i in range(3))
                    if on_polygon(Xp, loop, N):
                        rim = min(rim, abs(d0)-A)
                    continue
                w = mp.acos(-d0/A)
                for t in (g0+w, g0-w):
                    X = tuple(c[i]+rad*(mp.cos(t)*ex[i]+mp.sin(t)*ey[i]) for i in range(3))
                    if K.kind == 'wedge' and ref.dot(X, (ex[0], ex[1], 0)) < 0:
                        continue
                    if on_polygon(X, loop, N):
                        rim = min(rim, mp.sqrt(A*A-d0*d0)/rad)
    edges = set(e for piece in P.pieces for e in piece.edges)
    for i, j in edges:
        vi, vj = pc[i], pc[j]
        e = ref.sub(vj, vi)
        em = ref.Mv(e)
        ne = mp.sqrt(ref.dot(em, em))
        for t in ref.real_roots(K.T.line_quartic(vi, e)):
            if not 0 <= t <= 1:
                continue
            X = tuple(pm[i][k]+t*em[k] for k in range(3))
            if not on_wall(X):
                continue
            s2 = X[0]**2+X[1]**2+X[2]**2+R*R-r*r
            grad = (4*s2*X[0]-8*R*R*X[0], 4*s2*X[1]-8*R*R*X[1], 4*s2*X[2])
            ng = mp.sqrt(ref.dot(grad, grad))
            edge = min(edge, abs(ref.dot(grad, em))/(ng*ne))
        for c, rad, nrm, ex, ey in rims(K):
            den = ref.dot(em, nrm)
            if den == 0:
                continue
            t = ref.dot(ref.sub(c, pm[i]), nrm)/den
            if not 0 < t < 1:
                continue
            X = tuple(pm[i][k]+t*em[k] for k in range(3))
            d = ref.sub(X, c)
            if K.kind == 'wedge' and ref.dot(X, (ex[0], ex[1], 0)) < 0:
                continue
            if ref.dot(d, d) <= rad*rad:
                edge = min(edge, abs(den)/ne)
    return {'vertex': vertex/size, 'face': face/size, 'edge_crossing': edge, 'rim': rim}


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, first, ops, mc_n = job
    obj, tool = first.obj, first.tool
    pair = ref.Pair(obj, tool)
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    size = pair.size
    r = pair.sliced()
    m = pair.swept()
    VK, MK, AK = pair.K.closed()
    VP, MP, AP = pair.P.closed()
    checks = {}
    dev = mp.mpf(0)
    for res_ in (r, m):
        for (v, mm), (V, Mm) in ((res_['K'], (VK, MK)), (res_['P'], (VP, MP))):
            dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(mm, Mm))/size**4)
    checks['inputs_sliced'] = dev
    dev = mp.mpf(0)
    vc, mc = m['common']
    for key, (V, Mm) in (('fuse', (VK+VP-vc, [a+b-c for a, b, c in zip(MK, MP, mc)])),
                         ('K-P', (VK-vc, [a-c for a, c in zip(MK, mc)])),
                         ('P-K', (VP-vc, [a-c for a, c in zip(MP, mc)]))):
        v, mm = m[key]
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(mm, Mm))/size**4)
    checks['inclusion_exclusion'] = dev
    faces = pair.faces()
    checks['face_classes'] = max(abs(sum(cls.values())-area) for _, _, cls, area in faces)/size**2
    opp = pair.classes('K')['opp']+pair.classes('P')['opp']
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')+opp-AK-AP)/size**2
    checks['input_areas'] = max(abs(sum(sum(cls.values()) for w, _, cls, _ in faces if w == k)-A)
                                for k, A in (('K', AK), ('P', AP)))/size**2
    dev = mp.mpf(0)
    for key in ('K', 'P', 'common', 'fuse', 'K-P', 'P-K'):
        dev = max(dev, abs(m[key][0]-r[key][0])/size**3,
                  max(abs(x-y) for x, y in zip(m[key][1], r[key][1]))/size**4)
    checks['second_direction'] = dev
    checks['wall_two_ways'] = max(abs(r['wall'][k]-m['wall'][k]) for k in ('in', 'out'))/size**2
    dev = mp.mpf(0)
    ends = [cls for w, tag, cls, _ in faces if tag[0] == 'end']
    if pair.K.kind == 'segment':
        f = pair.K.face_factor((0, 0, 1))
        for cls, (i_, s_, o_, all_) in zip(ends, m['ends']):
            dev = max(dev, abs(cls['in']-i_*f), abs(cls['same']-s_*f), abs(cls['opp']-o_*f))
    else:
        for cls, ((c, s), _), (i_, all_) in zip(ends, pair.K.end_planes, r['ends']):
            f = pair.K.face_factor((-s, c, 0))
            dev = max(dev, abs(cls['in']+cls['same']+cls['opp']-i_*f))
    checks['ends_two_ways'] = dev/size**2
    mcr = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mcr[op]
        V, mm = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in mm]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(pair)
    fuse_swept = pair.swept_count('fuse')
    assert fuse_swept == pair.solids('fuse'), (name, 'the fuse swept', fuse_swept, 'by the rule', pair.solids('fuse'))
    stats = {'volume_quadrature': max(pair.normal.quad_error, pair.meridian.quad_error)/size**4,
             'breaks': (len(pair.normal.breaks), len(pair.meridian.breaks)), 'margins': margins(pair),
             'solids': {k: pair.swept_count(k) for k in ('common', 'K-P', 'P-K', 'fuse')}}
    return name, rows, res, checks, near, stats


def near_coincidences(pair):
    """A face class, or a result's volume, positive but thinner than 1e-9 of
    the case's size; slicing breakpoints of either direction closer than
    that (and farther than 1e-30)."""
    size = pair.size
    out = []
    for w, tag, cls, _ in pair.faces():
        for c, v in cls.items():
            if mp.mpf(10)**-30*size**2 < abs(v) < mp.mpf(10)**-9*size**2:
                out.append(f'{w} {tag} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volume(op)[0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    # A wedge's rounded end (6.1e-17 off pi/2, 1.2e-16 off pi) beside the
    # angle where the half-plane is parallel to a face normal to the ideal
    # end's plane: no meeting (the planes cross 5e15 away), left out.
    Um = pair.K.Um if pair.K.kind == 'wedge' else None
    rounded = lambda x: Um is not None and abs(x-(((Um+mp.pi) % (2*mp.pi))-mp.pi)) < mp.mpf(10)**-35
    for label, b, unit in (('levels', pair.normal.raw_breaks, size), ('angles', pair.meridian.raw_breaks, 1)):
        b = sorted(b)
        for p, q in zip(b, b[1:]):
            if label == 'angles' and (rounded(p) or rounded(q)) and q-p < mp.mpf(10)**-15:
                continue
            if mp.mpf(10)**-30*unit < q-p < mp.mpf(10)**-9*unit:
                out.append(f'{label} {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, pair_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15, 'closed_forms_rounded_ends': 1e-15,
          'closed_forms_degenerate': 1e-25, 'inputs_sliced': 1e-30, 'inclusion_exclusion': 1e-30,
          'face_classes': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30, 'second_direction': 1e-30,
          'wall_two_ways': 1e-30, 'ends_two_ways': 1e-30, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-30}


def form_kind(first):
    if first.form is not None and first.form[1] == 'rounded':
        return 'rounded_ends'
    return 'exact_frames' if exact_pair(first) else 'turned_frames'


def reference_checks(results, forms):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1
    by_pair = {r[0]: r for r in results}
    firsts = {}
    for c in cases():
        firsts.setdefault(c.pair_name, c)
    for name, form in forms.items():
        if form is None:
            continue
        _, _, res, _, _, _ = by_pair[name]
        kind = form_kind(firsts[name])
        if any(c.kind == 'degenerate' for c in cases() if c.pair_name == name):
            kind = 'degenerate'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}' if kind != 'degenerate' else 'closed_forms_degenerate', dev)
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, torus_segment_boolean_reference.py: expect KIND {STEP}, reason TEXT for a '
           'degenerate case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, _, _, _ = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (obj, tool) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in cases():
        for label, c in (('obj', case.obj), ('tool', case.tool)):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\t{label} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n, only=None):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items() if only is None or name in only]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert (c.obj.torus is None) != (c.tool.torus is None), f'{c.name}: a torus part and a prism'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', nargs='*', help='evaluate these pairs only (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    if args.only:
        results = run(jobs(args.samples, set(args.only)), evaluate, args.workers)
        firsts = {c.pair_name: c for c in listed}
        for name, rows, res, checks, near, stats in results:
            print(name, rows, {k: mp.nstr(v, 3) for k, v in checks.items()}, near,
                  {k: (mp.nstr(v, 3) if not isinstance(v, (dict, tuple)) else v) for k, v in stats.items()})
            form = pair_forms(firsts[name])
            if form:
                for op, want in form.items():
                    n, vol, area, centre = res[op]
                    if want is None:
                        print('  form', op, 'empty; reference', n)
                    else:
                        print('  form', op, mp.nstr(vol-want[0], 3), mp.nstr(area-want[1], 3),
                              [mp.nstr(centre[i]-want[2][i], 3) for i in range(3)])
        return
    results = run(jobs(args.samples), evaluate, args.workers)
    degenerate = {c.pair_name for c in listed if c.kind == 'degenerate'}
    near = [f'{name}: near coincidences {n}' for name, _, _, _, n, _ in results if n and name not in degenerate]
    if near:
        raise SystemExit('\n'.join(near))
    firsts = {}
    for c in listed:
        firsts.setdefault(c.pair_name, c)
    forms = dict(run(list(firsts.values()), closed_job, args.workers))
    worst, covered = reference_checks(results, forms)
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
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames')
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('closed forms:', sum(1 for f in forms.values() if f is not None), 'of', len(forms), 'pairs')
    by_pair = {r[0]: r for r in results}
    least = {}
    for name, _, _, _, n, stats in results:
        if n:
            print('declared near coincidences:', name, '; '.join(n[:4]))
        if name not in degenerate:
            for k, v in stats['margins'].items():
                if v < least.get(k, (mp.inf, None))[0]:
                    least[k] = (v, name)
    print('least margins (non-degenerate pairs):', ', '.join(f'{k} {mp.nstr(v, 3)} ({n})'
                                                             for k, (v, n) in sorted(least.items())))
    for name, _, _, _, _, stats in results:
        print('solids', name, stats['solids'])
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
