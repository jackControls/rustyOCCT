#!/usr/bin/env python3
"""Fixtures for S9f.1 of REVIEW_NOTES.md: Booleans of spline prisms against
polyhedral prisms in any relative position (a profile of lines, arcs and
S8b's nonrational splines in one frame, a profile of lines in another,
turned, leaning or tilted, same-axis pairs whose offset rounds included).

`boolean-spline-any-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); in every case exactly one
profile holds a spline and the other lines only, and the two frames differ
(their axes, or, for `bulge_offset`, the origins' offset, which is not
binary64 in the object's frame: S9a's domain excludes it). The profiles are
S9a.2's (`generate_boolean_fixtures.spline_profiles`: the bulge, dome,
blob, wave, capsule and lens hole) and `kink` (R4 of the decisions: a
quadratic whose interior knot at its apex has multiplicity two, C1 exactly
in its frame, used in the turned `TILT`, where its lifted poles stay C1:
the kernel's extrusion refuses a spline whose lifted poles lose C1 at such
a knot, `InvalidTopology("edge_not_c1")`, before any Boolean). Frames are S9c.1's
(`generate_curved_boolean_fixtures.FRAMES`, stored bit for bit by the
kernel's `Frame3::new`; `boolean-spline-any-frames.tsv` records them).

`boolean-spline-any-r4-*` hold the fixtures added with S9f.1's kernel,
after the native capture (the compared set is the 38 above, unchanged):
`knot` (R4 of the decisions: a quadratic from (10, 2) to (0, 6) whose
interior knot of multiplicity two at (5, 4), C1 exactly there along (-2,
1), loses C1 lifted into `TILT`, the knot's pole 2^-53 off its neighbours'
midpoint: the deferred rounded-knot fixture, under a box whose planes crease
its wall across the knot and cut it in generatrices of planes exactly
parallel to its axis), and the blob in `TILT` against a `SIDE` box (`n_A .
m` exactly zero for its caps: generatrices, where `TILTX`'s is -8.9e-17,
`blob_rounding`).

`boolean-spline-any-expected.tsv` gives per case, from
`curved_boolean_reference.py` with S9f.1's spline walls:

* `expect KIND S9f.1`: the declared outcome (`solid`: one or more solids;
  `empty`: none; `degenerate`: a plane tangent to a spline wall along a
  generatrix, tangent to it at a knot, or within rounding of its axis,
  `Degenerate` in the decisions), then for a degenerate case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

Before writing, `reference_checks` compares the reference with
independent results (limits relative to the case's size):

* every operation two ways: the slicing's volume and first moments against
  the divergence theorem over the face sweeps' kept pieces
  (`Pair.divergence_volumes`); `fuse = A + B - common` and `cut = A -
  common` (volumes and moments, `A` and `B` in closed form, each
  operation's slices apart); `area(fuse) + area(common) = area(A) +
  area(B)` (pairs without faces shared with opposite orientations); every
  face's classes summing to its closed-form area; both sides' shared areas
  equal;
* closed forms: each spline profile's prism cut by a plane (a half-space
  box: one face across the prism, the others clear of it), oblique
  (creases) and parallel to its axis (generatrices), against S8b's split
  reference (`split_reference.rows`, an independent slicing of the
  profile's exact Bezier pieces), its sides the common and the cut; and
  `bulge_offset`, whose frames share their axes, against S9a.2's
  `SplinePair` (atoms by Green's theorem, meetings by root finding) with
  the rounding offset taken exactly (turned frames: S9a's frame
  coordinates take the stored axes as orthonormal, hence 1e-15);
* S9a.2's own fixtures: every one of its spline cases where one profile
  holds a spline and the other lines only (same frames, prisms and
  stacks), their 25 printed digits within 1e-24 in `XY` and 1e-15 in
  `TILT`, with S9a.2's solid counts;
* solid counts: each case's declared count, the reference's by its
  slicing's union-find;
* margins that flag near coincidences: the scan of S9c.1 (a face's class
  or a result's volume positive but below 1e-9 of the size, slicing
  breakpoints closer than that) and a geometric one: every plane of the
  polyhedral prism against the spline walls (a plane parallel to the axis
  kept at least 1e-3 from tangency to every span, knots included; an
  oblique one at a sine of at least 1e-3 from the axis, its creases'
  extreme heights at least 1e-3 from the caps), every edge of it crossing
  a spline wall at a sine of at least 1e-3, and every vertex of either
  prism at least 1e-3 from the other's faces (a vertex on a face it shares
  the plane of excepted). The declared degenerate pairs must fail their
  margin (below 1e-12) and are exempt from the rest.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import os
from pathlib import Path
import struct

import mpmath as mp

from identity_reference import Boundary, Case, Spline, encode_boolean_case
from curve_surface_reference import stored_axes
import curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, at, square
from generate_boolean_fixtures import spline, spline_profiles

ROOT = Path(__file__).resolve().parents[1]
BOOLEAN_OPERATION = 93
STEP = 'S9f.1'

TANGENT_PLANE = 'a plane tangent to a spline wall along a generatrix'
KNOT_TANGENT = 'a plane tangent to a spline wall along the generatrix of a knot'
AXIS_ROUNDING = "a plane within rounding of a spline wall's axis"


def profiles():
    """S9a.2's spline profiles, `kink`: a quadratic from (8, 0) over its
    apex (4, 4) back to the origin, closed by the base line, its interior
    knot of multiplicity two (the degree) at the apex, C1 there exactly (the
    knot's pole the midpoint of its neighbours, equal spans); and `knot`
    (R4's of the evidence): a quadratic from (10, 2) to (0, 6) through the
    knot of multiplicity two at (5, 4), along (-2, 1), poles (7, 3), (5, 4),
    (3, 5), C1 there exactly, closed by lines through (0, 0) and (10, 0)."""
    p = spline_profiles()
    p['kink'] = [Boundary(points=[(0.0, 0.0), (8.0, 0.0)], segments=[
        None, spline(2, [(8.0, 0.0), (6.0, 4.0), (4.0, 4.0), (2.0, 4.0), (0.0, 0.0)], (0.0, 1.0, 2.0), (3, 2, 3))])]
    p['knot'] = [Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 2.0), (0.0, 6.0)], segments=[
        None, None,
        spline(2, [(10.0, 2.0), (7.0, 3.0), (5.0, 4.0), (3.0, 5.0), (0.0, 6.0)], (0.0, 1.0, 2.0), (3, 2, 3)),
        None])]
    return p


def prism(boundaries, frame, start, end, op):
    return Case('', 1e-7, op, frame, start, end, boundaries)


class Boolean:
    def __init__(self, name, operation, obj, tool, kind='solid', solids=1, reason=None, opposite=False,
                 near=False):
        self.name, self.operation, self.kind, self.solids, self.reason = name, operation, kind, solids, reason
        self.opposite, self.near = opposite, near
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = prism(*obj, op=91)
        self.obj.name = name
        self.tool = prism(*tool, op=92)
        self.tool.name = name
        self.frames = (obj[1], tool[1])

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, **kw):
    """Cases of one pair: `ops` maps each operation to its solid count, or to
    `('empty',)`, or to `('degenerate', reason)`."""
    out = []
    for op, want in ops.items():
        if isinstance(want, int):
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'solid', want, **kw))
        elif want[0] == 'empty':
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'empty', 0, **kw))
        else:
            out.append(Boolean(f'{name}_{op}', op, obj, tool, 'degenerate', None, want[1], **kw))
    return out


ONE3 = {'fuse': 1, 'cut': 1, 'common': 1}


def cases():
    """Every class of S9f.1's decisions: oblique planes creasing spline walls
    (the bulge under a tilted box, a turned dome under a box, the wave in a
    tilted slab), planes parallel to the axis cutting them in generatrices
    (the dome between two side planes, the wave's slab walls, the turned
    dome's), a leaning box whose edges pierce the blob's wall, a tilted pin
    across a lens hole's walls, a turned box on the capsule's coplanar base
    and one standing on its top (a face shared with the opposite
    orientation), a same-axis pair whose offset rounds, the blob prism as a
    leaning tool, R4's knot of multiplicity two in a turned frame, and the
    degenerate classes: a plane tangent to the dome's wall along its apex's
    generatrix, one tangent to the kink's wall at its knot, and a plane
    within rounding of a turned blob's axis."""
    p = profiles()
    out = []
    out += group('bulge_tilt', (p['bulge'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(0.0, -8.0, 8.0, 8.0)], at('TILT', (6, 3, 2)), 0.0, 6.0), ONE3)
    out += group('dome_side', (p['dome'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(-1.0, 1.0, 3.0, 6.0)], at('SIDE', (1, 0, 0)), 0.0, 2.0), ONE3)
    out += group('blob_lean', (p['blob'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(-1.0, -1.0, 1.0, 1.0)], at('LEAN', (6.5, 3, 1)), 0.0, 6.0), ONE3)
    out += group('wave_tiltx', (p['wave'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(0.0, -2.0, 4.0, 1.5)], at('TILTX', (3, 6, 2.5)), -3.0, 3.0), ONE3)
    out += group('capsule_turned', (p['capsule'], at('XY', (0, 0, 0)), 0.0, 2.0),
                 ([square(-2.0, -0.75, 1.0, 0.75)], at('R125', (-0.5, 0, 0)), 0.0, 1.0), ONE3)
    out += group('capsule_stand', (p['capsule'], at('XY', (0, 0, 0)), 0.0, 2.0),
                 ([square(-2.0, -0.75, 1.0, 0.75)], at('R125', (-0.5, 0, 0)), 2.0, 3.0),
                 {'fuse': 1, 'cut': 1, 'common': ('empty',)}, opposite=True)
    out += group('lens_pin', (p['lens_hole'], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([square(-0.5, -2.5, 0.5, 2.5)], at('TILT2', (5, 5, 2)), -4.0, 4.0), {'fuse': 1, 'cut': 1, 'common': 2})
    out += group('bulge_offset', (p['bulge'], at('TILT', (1, -2, 0.5)), 0.0, 5.0),
                 ([square(7.0, 1.0, 13.0, 4.0)], at('TILT', (1.1, -2.3, 0.7)), -1.0, 3.0), ONE3)
    out += group('dome_turned', (p['dome'], at('TILT', (0, 0, 0)), 0.0, 5.0),
                 ([square(1.0, -5.0, 3.0, 8.0)], at('XY', (0, 0, 1)), 0.0, 8.0), ONE3)
    out += group('kink_lean', (p['kink'], at('TILT', (0, 0, 0)), 0.0, 4.0),
                 ([square(-1.0, -1.0, 1.0, 1.0)], at('LEAN', (4, 4.4, 0.8)), -3.0, 3.0),
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('blob_tool', ([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 4.0),
                 (p['blob'], at('LEAN', (6, 1, 1)), 0.0, 6.0), ONE3)
    out += group('dome_tangent', (p['dome'], at('XY', (0, 0, 0)), 0.0, 5.0),
                 ([square(2.0, -5.0, 5.0, 1.0)], at('TURN', (0, 0, -1)), 0.0, 7.0),
                 {op: ('degenerate', TANGENT_PLANE) for op in ('fuse', 'common')})
    out += group('kink_knot', (p['kink'], at('XY', (0, 0, 0)), 0.0, 4.0),
                 ([square(4.0, -9.0, 6.0, 1.0)], at('TURN', (0, 0, -1)), 0.0, 6.0),
                 {'cut': ('degenerate', KNOT_TANGENT)})
    out += group('blob_rounding', (p['blob'], at('TILT', (0, 0, 0)), 0.0, 5.0),
                 ([square(-6.0, -6.0, 12.0, 6.0)], at('TILTX', (0, 3.2, 2.4)), 0.0, 6.0),
                 {op: ('degenerate', AXIS_ROUNDING) for op in ('cut', 'common')}, near=True)
    return out


def r4_cases():
    """The fixtures added with S9f.1's kernel (`boolean-spline-any-r4-*`):
    R4's rounded knot in `TILT` (`knot`, which the extrusion refused before
    the kernel removed the knot exactly before lifting) under an `XY` box
    whose bottom plane creases its wall across the knot (its top plane's
    crease leaving through the top cap) and whose `x` planes, exactly
    parallel to `TILT`'s axis, cut it in generatrices (the box off the
    plane `z = 0`, which holds the prism's base edge: an edge in the plane
    of a face of the other, off that face, is refused, S9c.1's rule);
    and the blob in `TILT` against a `SIDE` box, its caps exactly parallel
    to the blob's axis (generatrices; `TILTX`'s, -8.9e-17 off, is
    `blob_rounding`'s degenerate rounding)."""
    p = profiles()
    out = []
    out += group('knot_tilt', (p['knot'], at('TILT', (0, 0, 0)), 0.0, 5.0),
                 ([square(0.0, 0.0, 4.0, 6.0)], at('XY', (3, 1, 0.25)), 0.0, 1.5), ONE3)
    out += group('blob_side_turned', (p['blob'], at('TILT', (0, 0, 0)), 0.0, 5.0),
                 ([square(1.0, -3.0, 6.0, 2.0)], at('SIDE', (5, 0, 0)), 0.0, 2.0), ONE3)
    return out


# ------------------------------------------------------------------ closed forms: half-spaces

def half_spaces():
    """(name, spline prism, half-space box, exact frame): each spline
    profile's prism and a box whose low cap's plane crosses it, every other
    face of the box clear of it: oblique planes (creases: `TILT`, `LEAN` and
    `TILTX` boxes) and planes parallel to the axis (generatrices: `SIDE`
    boxes), the spline prism in `XY` and in the turned `TILT`."""
    p = profiles()
    big = [square(-40.0, -40.0, 40.0, 40.0)]
    out = []
    for name, prof, fa, oa, fb, ob, h in (
            ('bulge_oblique', 'bulge', 'XY', (0, 0, 0), 'TILT', (0, 3, 2), 5.0),
            ('bulge_parallel', 'bulge', 'XY', (0, 0, 0), 'SIDE', (10.5, 0, 0), 5.0),
            ('dome_oblique', 'dome', 'XY', (0, 0, 0), 'LEAN', (2, 0, 2.5), 5.0),
            ('blob_parallel_turned', 'blob', 'TILT', (0, 0, 0), 'SIDE', (5, 0, 0), 5.0),
            ('wave_oblique', 'wave', 'XY', (0, 0, 0), 'TILTX', (4, 6, 2), 5.0),
            ('capsule_parallel', 'capsule', 'XY', (0, 0, 0), 'SIDE', (-1, 0, 0), 2.0),
            ('lens_hole_oblique_turned', 'lens_hole', 'TILT', (0, 0, 0), 'LEAN', (5, 3, 4), 4.0),
            ('kink_parallel_turned', 'kink', 'TILT', (0, 0, 0), 'SIDE', (3, 0, 0), 4.0)):
        A = prism(p[prof], at(fa, oa), 0.0, h, 91)
        B = prism(big, at(fb, ob), 0.0, 40.0, 92)
        out.append((name, A, B, fa == 'XY'))
    return out


def half_space_job(job):
    """The pair's common and cut against S8b's split reference: the sides
    of the box's low cap's plane (above: along the box's axis)."""
    import split_reference
    name, A, B, exact = job
    pair = ref.Pair(A, B)
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(B.frame))
    sides = {r[0]: r for r in split_reference.rows(A, o+ref.cross(x, y))}
    assert set(sides) == {'below', 'above'}, (name, set(sides))
    dev = mp.mpf(0)
    size = pair.size
    for op, side in (('common', 'above'), ('cut', 'below')):
        _, vol, area, centre = pair.result(op)
        _, V, S, C = sides[side]
        dev = max(dev, abs(vol-V)/size**3, abs(area-S)/size**2, max(abs(centre[i]-C[i]) for i in range(3))/size)
    return name, exact, dev


# ------------------------------------------------------------------ S9a.2's SplinePair

def exact_offset(obj, tool):
    """The origins' offset `(a, b, c)` in the object's stored axes, exactly
    (frames sharing their axes)."""
    axes = tuple(tuple(F(c) for c in v) for v in stored_axes(obj.frame))
    taxes = tuple(tuple(F(c) for c in v) for v in stored_axes(tool.frame))
    assert axes[1:] == taxes[1:], 'frames with different axes'
    o, x, y, n = axes
    D = tuple(taxes[0][i]-o[i] for i in range(3))
    den = ref.det3(x, y, n)
    return ref.det3(D, y, n)/den, ref.det3(x, D, n)/den, ref.det3(x, y, D)/den


def s9a2_result(obj, operation, tool):
    """S9a.2's reference (`boolean_reference.SplinePair`: atoms by Green's
    theorem over classified pieces, meetings by root finding) with the
    origins' offset taken exactly, binary64 or not (S9a's own domain
    requires binary64): (solids, volume, area, centre)."""
    import boolean_reference as s9a

    class ExactOffsetPair(s9a.SplinePair):
        def __init__(self, obj, tool):
            self.axes = tuple(tuple(F(c) for c in v) for v in stored_axes(obj.frame))
            a, b, c = exact_offset(obj, tool)
            self.offset = (a, b, c)
            self.A = s9a.profile_elements(obj.boundaries, obj.tolerance)
            self.B = s9a.profile_elements(tool.boundaries, tool.tolerance, a, b)
            self.heights = {'A': tuple(sorted((F(obj.start), F(obj.end)))),
                            'B': tuple(sorted((F(tool.start)+c, F(tool.end)+c)))}
            coords = [abs(s9a.M(v)) for e in self.A+self.B for v in (*e.p, *e.q)]
            coords += [abs(s9a.M(e.c[i]))+s9a.M(e.r) for e in self.A+self.B if e.kind == 'A' for i in range(2)]
            coords += [abs(s9a.M(v)) for e in self.A+self.B if e.kind == 'S' for ctrl in e.ctrls for p in ctrl
                       for v in p]
            self.scale = max([s9a.M(1)]+coords)
            self.eps = s9a.M(10)**-25*self.scale
            self._meetings = {}
            self._bands()
            self.slice_atoms = self.atoms
            self.classes, self.atoms = self._green()

    solids, V, S, centre, _ = ExactOffsetPair(obj, tool).result(operation)
    return solids, V, S, centre


def s9a2_cases():
    """S9a.2's spline fixtures whose one profile holds a spline and the
    other lines only (no tangency: S9a.2's `dome_tangent_*` and
    `dome_touch_*` are left out)."""
    import generate_boolean_fixtures as s9a
    out = []
    for c in s9a.spline_cases():
        sa, sb = (any(b.segments is not None and any(isinstance(s, Spline) for s in b.segments)
                      for b in x.boundaries) for x in (c.obj, c.tool))
        arcs = any(b.circle is not None or (b.segments and any(isinstance(s, tuple) for s in b.segments))
                   for b in (c.tool if sa else c.obj).boundaries)
        if sa != sb and not arcs and 'tangent' not in c.name and 'touch' not in c.name:
            out.append(c)
    return out


def s9a2_job(case):
    rows, _ = ref.rows(case.obj, case.operation, case.tool)
    return case.name, rows[0]


# ------------------------------------------------------------------ margins

INF = mp.mpf('inf')


def surface_distance(prism, X):
    """The distance, in the prism's frame coordinates, from the exact point
    `X` to the prism's surface; a height term exactly zero (a point in the
    plane of a cap, as coplanar caps put it) is left out."""
    u, v, w = ref.apply(prism.inv, ref.sub(X, prism.o))
    um, vm, wm = ref.M(u), ref.M(v), ref.M(w)
    lo, hi = prism.lo, prism.hi
    heights = [ref.M(abs(w-h)) for h in (lo, hi) if w != h]
    dh = mp.mpf(0) if lo <= w <= hi else min(heights)
    bd = prism.profile.boundary_distance(um, vm)
    inside = prism.profile.inside(um, vm) and bd > 0
    if inside and lo <= w <= hi:
        return min([bd]+heights)
    if inside:
        return dh
    if lo <= w <= hi:
        return bd
    return mp.sqrt(bd*bd+dh*dh)


def margins(pair):
    """Distances from the degeneracies of the decisions (see the module's
    docstring): {'tangent', 'axis', 'crease', 'edge', 'vertex'}."""
    A, B = pair.A, pair.B
    S, Q = (A, B) if A.profile.splines else (B, A)
    assert not Q.profile.splines, 'S9f.1: splines in one profile'
    out = {k: INF for k in ('tangent', 'axis', 'crease', 'edge', 'vertex')}
    spans = [el for el in S.profile.elements if el.kind == 'spline']
    nS = S.n
    nlen = mp.sqrt(ref.M(ref.dot(nS, nS)))
    lo, hi = ref.M(S.lo), ref.M(S.hi)
    for f in Q.faces:
        N, P0 = f.normal, f.point
        a0, ax, ay, an = ref.dot(N, ref.sub(S.o, P0)), ref.dot(N, S.x), ref.dot(N, S.y), ref.dot(N, nS)
        Nlen = mp.sqrt(ref.M(ref.dot(N, N)))
        for el in spans:
            g = ref.padd(ref.padd(ref.pscale(el.X, ax), ref.pscale(el.Y, ay)), [a0])
            dg = ref.pder(g)
            if ref.pzero(dg):
                continue
            crit = [t for t in ref.real_roots(dg) if -mp.mpf(10)**-30 <= t <= 1+mp.mpf(10)**-30]
            for tau in crit:
                value = ref.peval([ref.M(c) for c in g], tau)
                if an == 0:
                    out['tangent'] = min(out['tangent'], abs(value)/mp.sqrt(ref.M(ax*ax+ay*ay)))
                else:
                    wstar = -value/ref.M(an)
                    out['crease'] = min(out['crease'], min(abs(wstar-lo), abs(wstar-hi))*nlen)
        if an != 0:
            out['axis'] = min(out['axis'], abs(ref.M(an))/(Nlen*nlen))
            if ax == 0 and ay == 0:
                w0 = -a0/an
                if w0 not in (S.lo, S.hi):
                    out['crease'] = min(out['crease'], min(abs(ref.M(w0)-lo), abs(ref.M(w0)-hi))*nlen)
    # Edges of the polyhedral prism crossing spline walls.
    edges, vertices = [], []
    for vtx in Q.profile.vertices:
        P, R = Q.world_exact(vtx[0], vtx[1], Q.lo), Q.world_exact(vtx[0], vtx[1], Q.hi)
        edges.append((P, R))
        vertices += [P, R]
    for el in Q.profile.elements:
        for h in (Q.lo, Q.hi):
            edges.append((Q.world_exact(el.p[0], el.p[1], h), Q.world_exact(el.q[0], el.q[1], h)))
    for P, R in edges:
        lp, lr = ref.apply(S.inv, ref.sub(P, S.o)), ref.apply(S.inv, ref.sub(R, S.o))
        D = ref.Mv(ref.sub(R, P))
        for el in spans:
            if (lr[0]-lp[0], lr[1]-lp[1]) == (0, 0):
                continue
            for tau, t in el.point_events(lp[:2], (lr[0]-lp[0], lr[1]-lp[1])):
                w = ref.M(lp[2])+t*ref.M(lr[2]-lp[2])
                if not (0 <= tau <= 1 and 0 <= t <= 1 and lo <= w <= hi):
                    continue
                dx, dy = el.tangent(tau)
                T = tuple(dx*S.xm[i]+dy*S.ym[i] for i in range(3))
                Nw = ref.cross(T, S.nm)
                sine = abs(ref.dot(D, Nw))/mp.sqrt(ref.dot(D, D)*ref.dot(Nw, Nw))
                out['edge'] = min(out['edge'], sine)
    # Vertices against the other's surface (a spline prism's knots included).
    for X in vertices:
        out['vertex'] = min(out['vertex'], surface_distance(S, X))
    for vtx in S.profile.vertices:
        for h in (S.lo, S.hi):
            out['vertex'] = min(out['vertex'], surface_distance(Q, S.world_exact(vtx[0], vtx[1], h)))
    return out


MARGIN = mp.mpf(10)**-3
DEGENERATE_MARGIN = {TANGENT_PLANE: 'tangent', KNOT_TANGENT: 'tangent', AXIS_ROUNDING: 'axis'}


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, obj, tool, ops, opposite = job
    pair = ref.Pair(obj, tool)
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    A, B = pair.A, pair.B
    size = pair.size
    va, ma, aa = A.measures()
    vb, mb, ab = B.measures()
    vols = pair.volumes()
    checks = {}
    # Inclusion and exclusion (operations computed apart).
    dev = mp.mpf(0)
    vf, mf = vols['fuse']
    vc, mc = vols['common']
    vt, mt = vols['cut']
    for x, y in ((vf, va+vb-vc), (vt, va-vc)):
        dev = max(dev, abs(x-y)/size**3)
    for i in range(3):
        dev = max(dev, abs(mf[i]-(ma[i]+mb[i]-mc[i]))/size**4, abs(mt[i]-(ma[i]-mc[i]))/size**4)
    checks['inclusion_exclusion'] = dev
    # Every operation a second way: the divergence theorem over the pieces.
    div = pair.divergence_volumes()
    dev = mp.mpf(0)
    for op in ref.OPS:
        dev = max(dev, abs(div[op][0]-vols[op][0])/size**3,
                  max(abs(x-y) for x, y in zip(div[op][1], vols[op][1]))/size**4)
    checks['divergence'] = dev
    # Faces: classes against closed-form areas; the area identity.
    fdev = mp.mpf(0)
    for tag, f, cls, _ in pair.face_areas():
        fdev = max(fdev, abs(sum(cls.values())-f.closed_area())/size**2)
    checks['face_classes'] = fdev
    total = lambda t, c: sum((cls[c] for tag, f, cls, _ in pair.face_areas() if tag == t), mp.mpf(0))
    checks['shared_both_sides'] = max(abs(total('A', 'same')-total('B', 'same')),
                                      abs(total('A', 'opp')-total('B', 'opp')))/size**2
    if not opposite:
        checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')-aa-ab)/size**2
    # A second slicing direction (parallel axes).
    if pair.slicing.parallel:
        axes = [(F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))]
        used = min(axes, key=lambda e: abs(ref.dot(A.n, e)))
        other = [e for e in axes if e != used and ref.cross(A.n, e) != (0, 0, 0)][0]
        second = ref.Slicing(A, B, other).measure()
        dev = mp.mpf(0)
        for op in ref.OPS:
            dev = max(dev, abs(second[op][0]-vols[op][0])/size**3)
            dev = max(dev, max(abs(x-y) for x, y in zip(second[op][1], vols[op][1]))/size**4)
        checks['second_axis'] = dev
    near = near_coincidences(pair)
    stats = {'volume_quadrature': pair.slicing.quad_error/size**4,
             'face_quadrature': max(sw.quad_error for _, _, _, sw in pair.face_areas())/size**2,
             'breaks': len(pair.slicing.breaks)}
    return name, rows, {op: res[op][0] for op in ref.OPS}, res, checks, near, margins(pair), stats


def near_coincidences(pair):
    """S9c.1's scan: a face class, or a result's volume, positive but thinner
    than 1e-9 of the case's size; slicing breakpoints closer than that."""
    size = pair.size
    out = []
    for tag, f, cls, _ in pair.face_areas():
        for c, v in cls.items():
            if 0 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{tag} {ref.describe(f)} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volumes()[op][0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    b = pair.slicing.breaks
    for p, q in zip(b, b[1:]):
        if q-p < mp.mpf(10)**-9*size:
            out.append(f'breakpoints {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def expected_file(path):
    out = {}
    for line in path.read_text().splitlines()[1:]:
        name, row = line.split('\t')
        if row.split()[0] in ('result', 'empty'):
            out[name] = row
    return out


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


def reference_checks(results, workers):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, res, checks, near, margin, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
        note('face_quadrature_estimate', stats['face_quadrature'])
    # Closed forms: half-spaces against S8b's split reference.
    for name, exact, dev in run(half_spaces(), half_space_job, workers):
        note(f'split_reference_{"exact" if exact else "turned"}_frames', dev)
    # The same-axis pair whose offset rounds against S9a.2's SplinePair.
    by_pair = {r[0]: r for r in results}
    for c in cases():
        if c.pair_name != 'bulge_offset':
            continue
        a, b, cc = exact_offset(c.obj, c.tool)
        assert any(F(float(v)) != v for v in (a, b, cc)), 'bulge_offset: the offset is binary64'
        n, V, S, C = s9a2_result(c.obj, c.operation, c.tool)
        got = by_pair[c.pair_name][3][c.operation]
        assert n == got[0], (c.name, n, got[0])
        size = max(abs(V)**(mp.mpf(1)/3), 1)
        note('s9a2_rounded_offset', max(abs(got[1]-V)/abs(V), abs(got[2]-S)/abs(S),
                                        max(abs(got[3][i]-C[i]) for i in range(3))/size))
    # S9a.2's own spline fixtures (one spline profile, the other lines).
    want = expected_file(ROOT/'fixtures'/'boolean-spline-expected.tsv')
    listed = s9a2_cases()
    frames = {c.name: frame_name(c.obj.frame) for c in listed}
    for case_name, row in run(listed, s9a2_job, workers):
        got, exp = row.split(), want[case_name].split()
        assert got[0] == exp[0] and got[1:2] == exp[1:2], (case_name, got, exp)
        dev = mp.mpf(0)
        for x, y in zip(got[2:], exp[2:]):
            dev = max(dev, abs(mp.mpf(x)-mp.mpf(y))/max(1, abs(mp.mpf(y))))
        note(f's9a2_fixtures_{"exact" if frames[case_name] == "XY" else "turned"}_frames', dev)
    return worst, covered


def generate(results, listed, prefix):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, curved_boolean_reference.py with spline walls: expect KIND {STEP}, reason TEXT '
           'for a degenerate case, then result N volume area cx cy cz or empty)']
    for case in listed:
        blocks.append(case.encode())
        _, rows, _, res, _, near, _, _ = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        if case.kind == 'solid':
            assert n == case.solids, f'{case.name}: declared {case.solids} solids, the reference gives {row[0]}'
        elif case.kind == 'empty':
            assert n == 0, f'{case.name}: declared empty, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (obj, tool) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in listed:
        for label, c in (('obj', case.obj), ('tool', case.tool)):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\t{label} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{prefix}-cases.txt': '\n'.join(blocks)+'\n',
            f'{prefix}-expected.tsv': '\n'.join(out)+'\n',
            f'{prefix}-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(listed):
    pairs = {}
    for c in listed:
        pairs.setdefault(c.pair_name, [c.obj, c.tool, [], c.opposite])[2].append(c.operation)
    return [(name, obj, tool, ops, opposite) for name, (obj, tool, ops, opposite) in pairs.items()]


def validate(listed):
    """The fixture list's own rules (see the module's docstring)."""
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        has = [any(b.segments is not None and any(isinstance(s, Spline) for s in b.segments) for b in x.boundaries)
               for x in (c.obj, c.tool)]
        assert sum(has) == 1, f'{c.name}: a spline in exactly one profile'
        other = c.tool if has[0] else c.obj
        assert all(b.circle is None and not (b.segments and any(s is not None for s in b.segments))
                   for b in other.boundaries), f'{c.name}: the other profile of lines only'
        if stored_axes(c.obj.frame)[1:] == stored_axes(c.tool.frame)[1:]:
            a, b, cc = exact_offset(c.obj, c.tool)
            assert any(F(float(v)) != v for v in (a, b, cc)), f'{c.name}: an S9a pair (a binary64 offset)'
        frame_name(c.obj.frame), frame_name(c.tool.frame)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    args = parser.parse_args()
    listed = cases()+r4_cases()
    validate(listed)
    results = run(jobs(listed), evaluate, args.workers)
    declared = {c.pair_name: c for c in listed}
    worst_margin, degenerate_margin = {}, {}
    for name, _, _, _, _, near, margin, _ in results:
        first = declared[name]
        kinds = {c.kind for c in listed if c.pair_name == name}
        if 'degenerate' in kinds:
            key = DEGENERATE_MARGIN[first.reason]
            assert margin[key] < mp.mpf(10)**-12, f'{name}: declared degenerate, {key} margin {margin[key]}'
            degenerate_margin[key] = max(degenerate_margin.get(key, mp.mpf(0)), margin[key])
            continue
        if near and not first.near:
            raise SystemExit(f'{name}: near coincidences {near}')
        for key, value in margin.items():
            if value < MARGIN:
                raise SystemExit(f'{name}: {key} margin {mp.nstr(value, 3)}')
            worst_margin[key] = min(worst_margin.get(key, INF), value)
    worst, covered = reference_checks(results, args.workers)
    limits = {'inclusion_exclusion': 1e-30, 'divergence': 1e-30, 'face_classes': 1e-30,
              'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'second_axis': 1e-30,
              'volume_quadrature_estimate': 1e-30, 'face_quadrature_estimate': 1e-30,
              'split_reference_exact_frames': 1e-30, 'split_reference_turned_frames': 1e-15,
              's9a2_rounded_offset': 1e-15, 's9a2_fixtures_exact_frames': 1e-24,
              's9a2_fixtures_turned_frames': 1e-15}
    for key, value in worst.items():
        assert value <= limits[key], (key, mp.nstr(value, 3))
    files = {**generate(results, cases(), 'boolean-spline-any'),
             **generate(results, r4_cases(), 'boolean-spline-any-r4')}
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('smallest margins (non-degenerate pairs):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst_margin.items())))
    print('declared degenerate margins (largest):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(degenerate_margin.items())))
    for name, _, _, _, _, near, _, _ in results:
        if near:
            print('declared near coincidences:', name, '; '.join(near[:4]))


if __name__ == '__main__':
    main()
