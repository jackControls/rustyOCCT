#!/usr/bin/env python3
"""Fixtures for S8a of REVIEW_NOTES.md: prisms of line and arc profiles
split by a plane.

`split-cases.txt` lists each case in the protocol of `identity-cases.txt`
(a prism: `case`, `op`, `frame`, `offsets`, `boundary` rows) with a `split ox
oy oz nx ny nz` row (the plane through a point with a normal) before `end`;
`split-expected.tsv` gives each case's sides from `split_reference.py`
(`side below|above|whole volume area cx cy cz`: the totals of the pieces on
each side), and `split-frames.tsv` each
prism's stored frame axes (`stored_axes`, the kernel's `Frame3::new`; its tests
check them bit for bit). Every class has a case: planes crossing the caps and
the walls at angles, normal to the axis (M3's height split), parallel to it
across line and arc walls, through two vertical edges, through one vertex
only, tangent to an arc wall along a ruling, lying in a cap, missing, and
crossing a hole, a non-convex profile in two pieces of section, and a
circle's whole wall in an ellipse. No Rust result supplies an expectation.

S8c: `split-primitive-cases.txt` gives one cone, frustum, sphere or zone per
line (`cone NAME tol ox oy oz nx ny nz xx xy xz r1 r2 h split ...` or
`sphere NAME tol ox oy oz nx ny nz xx xy xz R a1 a2 split ...`, as the
kernel's builders and `BRepPrimAPI_MakeCone`/`MakeSphere` take them) and
`split-primitive-expected.tsv` its sides (`split_reference.revolved_rows`):
planes normal to the axis (through a frustum, an apex cone, a zone, a whole
sphere, at an apex, in a cap, missing, touching a sphere), a whole sphere
by oblique planes, planes containing the axis, and the conic and circle
sections that wait for procedural edges (S8d). S8d: `split-torus-cases.txt`
the same with `torus NAME tol ... major minor angle` lines (whole tori, the
angle a full turn) and `split-torus-expected.tsv` (`split_reference.
torus_rows`): planes normal to the axis (through the equator, above it,
touching, missing, in a tilted frame), containing it, parallel to it
through the tube, and oblique.

S8b: `split-spline-cases.txt` holds prisms whose profiles have spline
segments (the `B` segment of an `S` boundary row, `identity_reference.
encode_case`), quadratic and cubic, with interior knots, beside lines and an
arc, in holes, two of them given clockwise; `split-spline-expected.tsv`
their sides from the same slicing (`spline_cases`). Before writing,
`reference_checks` compares the spline slicing with closed forms: straight
splines give the polygon's rows, and every profile of lines and splines its
Green's-theorem area and moments in exact Fractions, both within 1e-30.

S8e: `split-sheet-cases.txt` holds face and wire bodies (S6: a `make face`
or `make wire` row in place of `offsets`) with their `split` rows, and
`split-sheet-expected.tsv` their sides (`split_reference.planar_rows`) in
the rows of `split-expected.tsv` with the measures one dimension down: `side
S area perimeter cx cy cz` for a sheet (its area, the total length of its
boundary, chords included, and its centroid), `side S length 0 cx cy cz`
for a wire (its length, 0 for its ends' measure, and its centroid).
`planar_reference_checks` compares them with closed forms (a square's
halves, a U's and a holed square's perimeters, circles' arcs and
segments, a quadratic's length), exact half-plane clipping of polygons in
Fractions, straight splines, Green's theorem and the sides' sums, within
1e-30.
"""
import argparse
import dataclasses
from fractions import Fraction as F
import math
from pathlib import Path
import struct

from identity_reference import Boundary, Case, Spline, encode_case, reversed_segment, stored_points
from curve_surface_reference import stored_axes
import split_reference as ref

ROOT = Path(__file__).resolve().parents[1]
XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
TILT = (1.0, -2.0, 0.5, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0)


def arc(cx, cy, r, ccw=True):
    return (cx, cy, r, ccw)


def path(points, segments):
    return Boundary(points=points, segments=segments)


def prism(name, boundaries, frame=XY, start=0.0, end=5.0):
    return Case(name, 1e-7, 81, frame, start, end, boundaries)


def cases():
    square = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    hole = Boundary(points=[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)])
    ushape = Boundary(points=[(0.0, 0.0), (9.0, 0.0), (9.0, 6.0), (6.0, 6.0), (6.0, 2.0), (3.0, 2.0),
                              (3.0, 6.0), (0.0, 6.0)])
    stadium = path([(0.0, -1.0), (3.0, -1.0), (3.0, 1.0), (0.0, 1.0)],
                   [None, arc(3.0, 0.0, 1.0), None, arc(0.0, 0.0, 1.0)])
    rounded = path([(1.0, 0.0), (7.0, 0.0), (8.0, 1.0), (8.0, 5.0), (7.0, 6.0), (1.0, 6.0), (0.0, 5.0), (0.0, 1.0)],
                   [None, arc(7.0, 1.0, 1.0), None, arc(7.0, 5.0, 1.0), None, arc(1.0, 5.0, 1.0), None,
                    arc(1.0, 1.0, 1.0)])
    disc = Boundary(circle=(0.0, 0.0, 2.0))
    slot = path([(-2.0, -1.0), (2.0, -1.0), (2.0, 1.0), (-2.0, 1.0)],
                [None, arc(2.0, 0.0, 1.0), None, arc(-2.0, 0.0, 1.0)])
    frame10 = Boundary(points=[(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)])
    return [
        (prism('sq_height', [square]), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        (prism('sq_tilted', [square]), (5.0, 5.0, 2.5, 1.0, 0.0, 2.0)),
        (prism('sq_skew', [square]), (4.0, 6.0, 2.0, 1.0, 2.0, 3.0)),
        (prism('sq_vertical', [square]), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('sq_diagonal_edges', [square]), (5.0, 5.0, 0.0, 1.0, 1.0, 0.0)),
        (prism('sq_through_corner', [square]), (10.0, 10.0, 5.0, 1.0, 1.0, 1.0)),
        (prism('sq_in_cap', [square]), (3.0, 4.0, 5.0, 0.0, 0.0, 1.0)),
        (prism('sq_corner_cut', [square]), (0.0, 0.0, 7.0, 1.0, 1.0, 4.0)),
        (prism('sq_miss', [square]), (0.0, 0.0, 7.0, 0.25, 0.25, 4.0)),
        (prism('sq_hole_crossed', [square, hole]), (5.0, 5.0, 2.5, 1.0, 0.5, 1.5)),
        (prism('sq_hole_vertical', [square, hole]), (5.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('u_two_sections', [ushape]), (0.0, 4.0, 0.0, 0.0, 1.0, 0.0)),
        (prism('u_tilted', [ushape], start=0.0, end=3.0), (4.5, 3.0, 1.5, 0.25, 1.0, 0.5)),
        (prism('stadium_height', [stadium], end=2.0), (0.0, 0.0, 0.75, 0.0, 0.0, 1.0)),
        (prism('stadium_tilted', [stadium], end=2.0), (1.5, 0.0, 1.0, 1.0, 0.0, 2.0)),
        (prism('stadium_arc_only', [stadium], end=2.0), (3.5, 0.0, 1.0, 1.0, 0.0, 0.25)),
        (prism('stadium_vertical_arc', [stadium], end=2.0), (3.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('stadium_tangent', [stadium], end=2.0), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('rounded_skew', [rounded], end=3.0), (4.0, 3.0, 1.5, 1.0, 1.0, 2.0)),
        (prism('disc_wall_ellipse', [disc], end=4.0), (0.0, 0.0, 2.0, 1.0, 0.0, 4.0)),
        (prism('disc_through_caps', [disc], end=4.0), (0.0, 0.0, 2.0, 1.0, 0.0, 1.0)),
        (prism('disc_vertical', [disc], end=4.0), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('slot_hole_tilted', [frame10, slot], frame=TILT, start=0.0, end=3.0),
         (1.0, -2.0, 2.0, 0.25, 0.5, 1.0)),
        (prism('tilted_frame_in_cap', [square], frame=TILT, start=0.0, end=4.0),
         (1.0, -2.0, 0.5, 0.0, 3.0, 4.0)),
        (prism('tilted_frame_height', [square], frame=TILT, start=0.0, end=4.0),
         (1.0, -0.8, 2.1, 0.0, 3.0, 4.0)),
        (prism('reversed_offsets', [stadium], start=2.0, end=-1.0), (1.0, 0.0, 0.5, 1.0, 1.0, 2.0)),
    ]


HALF_PI = math.pi/2


def primitive_cases():
    """(kind, name, frame, params, plane)."""
    cone, sphere = 'cone', 'sphere'
    return [
        (cone, 'frustum_height', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 1.25, 0.0, 0.0, 1.0)),
        (cone, 'apex_height', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.5, 0.0, 0.0, 1.0)),
        (cone, 'apex_up_height', XY, (0.0, 1.5, 2.0), (0.0, 0.0, 0.5, 0.0, 0.0, -1.0)),
        (cone, 'widening_height', XY, (1.0, 2.5, 0.5), (0.0, 0.0, 0.375, 0.0, 0.0, 1.0)),
        (cone, 'apex_touch', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 3.0, 0.0, 0.0, 1.0)),
        (cone, 'frustum_in_cap', XY, (2.0, 1.0, 3.0), (0.5, 0.0, 0.0, 0.0, 0.0, 1.0)),
        (cone, 'frustum_miss', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 4.0, 0.0, 0.0, 1.0)),
        (cone, 'tilted_height', TILT, (2.0, 0.5, 3.0), (1.0, -0.8, 2.1, 0.0, 3.0, 4.0)),
        (cone, 'frustum_meridian', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (cone, 'apex_meridian', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.0, 1.0, 1.0, 0.0)),
        (cone, 'apex_oblique', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.0, 1.0, 0.0, 2.0)),
        (cone, 'frustum_parallel', XY, (2.0, 1.0, 3.0), (0.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (sphere, 'sphere_equator', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
        (sphere, 'sphere_height', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 0.0, 0.75, 0.0, 0.0, 1.0)),
        (sphere, 'sphere_oblique', XY, (2.0, -HALF_PI, HALF_PI), (0.3, 0.2, 0.1, 1.0, 2.0, 2.0)),
        (sphere, 'sphere_tilted', TILT, (1.5, -HALF_PI, HALF_PI), (1.25, -2.0, 0.5, 1.0, -1.0, 0.5)),
        (sphere, 'sphere_touch', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        (sphere, 'sphere_miss', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 3.0, 0.0, 0.0, 1.0, 0.0)),
        (sphere, 'zone_height', XY, (2.0, -0.5, 1.0), (0.0, 0.0, 0.25, 0.0, 0.0, 1.0)),
        (sphere, 'cap_height', XY, (2.0, 0.25, HALF_PI), (0.0, 0.0, 1.5, 0.0, 0.0, 1.0)),
        (sphere, 'zone_meridian', XY, (2.0, -0.5, 1.0), (0.0, 0.0, 0.0, 0.0, 1.0, 0.0)),
        (sphere, 'cap_meridian', XY, (2.0, 0.25, HALF_PI), (0.0, 0.0, 0.0, 1.0, 1.0, 0.0)),
        (sphere, 'zone_oblique', XY, (2.0, -0.5, 1.0), (0.0, 0.0, 0.5, 1.0, 0.0, 1.0)),
    ]


def conic_cases():
    """(kind, name, frame, params, plane): cones, frusta, zones and caps by
    planes neither normal to their axes nor containing them (S8d.2), beside
    S8c's apex_oblique, frustum_parallel and zone_oblique: a closed ellipse,
    a parabola, hyperbolas, a tongue at the top, the rulings through a
    frustum's virtual apex, an ellipse touching both rims, a zone's side cap
    and circle arcs across both rims."""
    cone, sphere = 'cone', 'sphere'
    return [
        (cone, 'apex_ellipse', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.0, -0.2, 0.0, 1.0)),
        (cone, 'apex_parabola', XY, (2.0, 0.0, 2.0), (0.5, 0.0, 0.0, 1.0, 0.0, 1.0)),
        (cone, 'apex_hyperbola', XY, (2.0, 0.0, 3.0), (0.8, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (cone, 'frustum_top_tongue', XY, (2.0, 1.0, 3.0), (0.6, 0.0, 3.0, 1.0, 0.0, 1.0)),
        (cone, 'frustum_rulings', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 6.0, 1.0, 0.0, 0.1)),
        (cone, 'frustum_touch_both', XY, (2.0, 1.0, 3.0), (2.0, 0.0, 0.0, 1.0, 0.0, 1.0)),
        (cone, 'widening_oblique', XY, (1.0, 2.5, 0.5), (0.0, 0.0, 0.25, 0.3, 0.2, 1.0)),
        (cone, 'tilted_oblique', TILT, (2.0, 0.5, 3.0), (1.0, -0.8, 2.1, 1.0, 3.0, 4.0)),
        (sphere, 'zone_side_cap', XY, (2.0, -0.5, 1.0), (1.8, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (sphere, 'cap_oblique', XY, (2.0, 0.25, HALF_PI), (0.0, 0.0, 1.2, 1.0, 0.0, 1.0)),
        (sphere, 'zone_both_rims', XY, (2.0, -0.5, 1.0), (0.3, 0.0, 0.0, 1.0, 0.0, 0.2)),
        (sphere, 'tilted_zone_oblique', TILT, (1.5, -0.7, 0.9), (1.25, -2.0, 0.5, 1.0, -1.0, 0.5)),
    ]


def torus_cases():
    """(name, frame, (major, minor), plane): whole tori (S8d)."""
    return [
        ('torus_equator', XY, (3.0, 1.0), (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
        ('torus_height', XY, (3.0, 1.0), (0.0, 0.0, 0.5, 0.0, 0.0, 1.0)),
        ('torus_high', XY, (3.0, 1.0), (0.0, 0.0, 0.875, 0.0, 0.0, -1.0)),
        ('torus_touch', XY, (3.0, 1.0), (0.0, 0.0, 1.0, 0.0, 0.0, 1.0)),
        ('torus_miss', XY, (3.0, 1.0), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        ('torus_meridian', XY, (3.0, 1.0), (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_meridian_skew', XY, (3.0, 1.0), (0.0, 0.0, 0.5, 1.0, 1.0, 0.0)),
        ('torus_tilted_height', TILT, (2.5, 0.75), (1.0, -1.82, 0.74, 0.0, 3.0, 4.0)),
        ('torus_parallel_outer', XY, (3.0, 1.0), (3.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_parallel_inner', XY, (3.0, 1.0), (2.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_oblique', XY, (3.0, 1.0), (0.0, 0.0, 0.25, 0.25, 0.0, 1.0)),
        ('torus_oblique_tube', XY, (3.0, 1.0), (3.0, 0.0, 0.0, 1.0, 0.0, 1.0)),
    ]


def spiric_cases():
    """(name, frame, (major, minor), plane): whole tori cut in spiric
    sections (S8d.3), beside S8d's four: two loops winding about the axis
    (bands), one contractible loop (a cap, the rest the torus less a disc)
    and two loops winding about the tube (C-shaped pieces), in both frames."""
    return [
        ('torus_gentle', XY, (3.0, 1.0), (0.0, 0.0, 0.0, 0.1, 0.0, 1.0)),
        ('torus_band_offset', XY, (3.0, 1.0), (0.0, 0.0, 0.3, 0.15, 0.05, 1.0)),
        ('torus_top_band', XY, (3.0, 1.0), (0.0, 0.0, 0.8, 0.05, 0.0, 1.0)),
        ('torus_top_cap', XY, (3.0, 1.0), (0.0, 0.0, 0.9, 0.3, 0.0, 1.0)),
        ('torus_cap_outer', XY, (3.0, 1.0), (3.6, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_peanut', XY, (3.0, 1.0), (2.2, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_hole_slice', XY, (3.0, 1.0), (2.05, 0.0, 0.0, 1.0, 0.0, 0.3)),
        ('torus_two_ovals', XY, (3.0, 1.0), (1.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_steep', XY, (3.0, 1.0), (0.0, 0.0, 0.0, 1.0, 0.0, 1.0)),
        ('torus_skew_ovals', XY, (3.0, 1.0), (0.5, 0.2, 0.0, 1.0, 0.3, 0.4)),
        ('torus_small_tube', XY, (4.0, 0.5), (0.0, 0.0, 0.1, 0.2, 0.0, 1.0)),
        ('torus_tilted_spiric', TILT, (2.5, 0.75), (1.0, -1.82, 0.74, 1.0, 3.0, 4.0)),
        ('torus_tilted_band', TILT, (2.5, 0.75), (1.0, -1.82, 0.74, 0.3, 3.0, 4.0)),
    ]


def spline(degree, poles, knots=(0.0, 1.0), mults=None):
    """A clamped nonrational spline segment (S8b); a single Bezier span by
    default."""
    if mults is None:
        mults = (degree+1,)+(1,)*(len(knots)-2)+(degree+1,)
    return Spline(degree, tuple(poles), tuple(knots), tuple(mults))


def clockwise(points, segments):
    """The same closed path given the other way round (its first point
    kept), each spline reversed: the stored order must restore it."""
    n = len(points)
    return path([points[0]]+points[:0:-1], [reversed_segment(segments[n-1-j]) for j in range(n)])


def framed(frame, point, normal):
    """A plane given in a frame's stored coordinates, in world coordinates
    (binary64)."""
    o, x, y, n = stored_axes(frame)
    at = tuple(o[i]+point[0]*x[i]+point[1]*y[i]+point[2]*n[i] for i in range(3))
    return at+tuple(normal[0]*x[i]+normal[1]*y[i]+normal[2]*n[i] for i in range(3))


def spline_profiles():
    """S8b's profiles: a rectangle with a quadratic bulge, a closed blob of
    four cubics (given clockwise), a rectangle under a quadratic wave with
    two interior knots, a stadium whose left end is a cubic with a double
    interior knot, and a square with a lens-shaped hole of two cubics (given
    clockwise)."""
    rect = [(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)]
    bulge = path(rect, [None, spline(2, [(10.0, 0.0), (12.0, 3.0), (10.0, 6.0)]), None, None])
    blob = clockwise([(4.0, 0.0), (8.0, 4.0), (4.0, 8.0), (0.0, 4.0)], [
        spline(3, [(4.0, 0.0), (6.5, 0.0), (8.0, 1.5), (8.0, 4.0)]),
        spline(3, [(8.0, 4.0), (8.0, 6.5), (6.5, 8.0), (4.0, 8.0)]),
        spline(3, [(4.0, 8.0), (1.5, 8.0), (0.0, 6.5), (0.0, 4.0)]),
        spline(3, [(0.0, 4.0), (0.0, 1.5), (1.5, 0.0), (4.0, 0.0)])])
    wave = path(rect, [None, None, spline(2, [(10.0, 6.0), (8.0, 8.0), (5.0, 5.0), (2.0, 8.0), (0.0, 6.0)],
                                          (0.0, 1.0, 2.0, 3.0)), None])
    capsule = path([(0.0, -1.0), (3.0, -1.0), (3.0, 1.0), (0.0, 1.0)],
                   [None, arc(3.0, 0.0, 1.0), None,
                    spline(3, [(0.0, 1.0), (-1.0, 1.0), (-1.5, 0.5), (-1.5, -0.5), (-1.0, -1.0), (0.0, -1.0)],
                           (0.0, 1.0, 2.0), (4, 2, 4))])
    square = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    lens = clockwise([(3.0, 5.0), (7.0, 5.0)], [
        spline(3, [(3.0, 5.0), (4.0, 2.0), (6.0, 2.0), (7.0, 5.0)]),
        spline(3, [(7.0, 5.0), (6.0, 8.0), (4.0, 8.0), (3.0, 5.0)])])
    return {'bulge': [bulge], 'blob': [blob], 'wave': [wave], 'capsule': [capsule], 'lens': [square, lens]}


def spline_cases():
    """(case, plane): prisms of S8b's profiles (`spline_profiles`) by planes
    normal to the axis, parallel to it (crossing a spline once, twice, four
    times, through an interior knot's point, through two joins of splines,
    tangent to a spline's extreme, beyond a spline but inside its control
    polygon), oblique (across spline walls and caps, across every wall and
    no cap, across a spline beside an arc, touching one join), and in a
    cap, in both frames."""
    p = spline_profiles()
    return [
        (prism('bulge_height', p['bulge']), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        (prism('bulge_parallel_once', p['bulge']), (0.0, 2.0, 0.0, 0.0, 1.0, 0.0)),
        (prism('bulge_parallel_twice', p['bulge']), (10.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('bulge_tangent', p['bulge']), (11.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('bulge_miss_hull', p['bulge']), (12.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('bulge_oblique', p['bulge']), (10.0, 3.0, 2.5, 1.0, 0.25, 1.5)),
        (prism('blob_walls_only', p['blob']), (4.0, 4.0, 2.5, 0.125, 0.0625, 1.0)),
        (prism('blob_joins', p['blob']), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('blob_join_touch', p['blob']), (8.0, 4.0, 5.0, 1.0, 0.0, 1.0)),
        (prism('blob_tilted_in_cap', p['blob'], frame=TILT, end=4.0), (1.0, -2.0, 0.5, 0.0, 3.0, 4.0)),
        (prism('wave_knot', p['wave']), (6.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('wave_parallel_four', p['wave']), (0.0, 6.8, 0.0, 0.0, 1.0, 0.0)),
        (prism('wave_tilted', p['wave'], frame=TILT, end=4.0), framed(TILT, (5.0, 6.0, 2.0), (0.5, 1.0, 1.0))),
        (prism('capsule_parallel', p['capsule'], end=2.0), (-1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('capsule_steep', p['capsule'], end=2.0), (1.5, 0.0, 1.0, 1.0, 0.25, 2.75)),
        (prism('lens_parallel', p['lens']), (0.0, 4.0, 0.0, 0.0, 1.0, 0.0)),
        (prism('lens_tilted', p['lens'], frame=TILT, end=3.0), framed(TILT, (5.0, 4.0, 1.5), (0.25, 1.0, 0.5))),
    ]


def reference_checks():
    """The spline reference against closed forms: a square whose sides are
    straight splines (a quadratic and a cubic with an interior knot, poles
    collinear and evenly spaced) gives the polygon's rows for S8a's square
    planes; and every spline profile of lines and splines has, by slicing,
    the area and first moments Green's theorem gives in exact Fractions.
    Returns the largest relative deviations; raises beyond 1e-30."""
    square = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    straight = path([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)], [
        spline(2, [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]),
        spline(3, [(10.0, 0.0), (10.0, 2.5), (10.0, 5.0), (10.0, 7.5), (10.0, 10.0)], (0.0, 1.0, 2.0)),
        None, None])
    worst = [ref.M(0), ref.M(0)]
    planes = [plane for case, plane in cases() if case.name in
              ('sq_height', 'sq_tilted', 'sq_skew', 'sq_vertical', 'sq_diagonal_edges', 'sq_corner_cut')]
    for plane in planes:
        want = ref.rows(prism('polygon', [square]), plane)
        got = ref.rows(prism('straight', [straight]), plane)
        assert [r[0] for r in want] == [r[0] for r in got]
        for w, g in zip(want, got):
            for a, b in zip([w[1], w[2], *w[3]], [g[1], g[2], *g[3]]):
                worst[0] = max(worst[0], abs(a-b)/max(1, abs(a)))
    seen = set()
    for case, _ in spline_cases():
        if any(b.circle is not None or (b.segments and any(isinstance(s, tuple) for s in b.segments))
               for b in case.boundaries) or case.name.split('_')[0] in seen:
            continue
        seen.add(case.name.split('_')[0])
        area, mx, my = ref.green_moments(case)
        h = ref.M(case.end-case.start)
        prism_ = ref.Prism(case, (0.0, 0.0, 0.0, 0.0, 0.0, 1.0))
        V, mu, mv, _ = prism_.volume_moments(None)
        cap = prism_.region_measure([], lambda u, v: True)
        for a, b in ((V/h, area), (mu/h, mx), (mv/h, my), (cap, area)):
            worst[1] = max(worst[1], abs(a-ref.M(b))/max(1, abs(ref.M(b))))
    if max(worst) > ref.M(10)**-30:
        raise SystemExit(f'spline reference checks failed: {worst}')
    return worst


def planar(name, boundaries, make, frame=XY):
    """A face (`Body::face_from_profile`) or wire (`Body::wire_from_boundary`)
    body of S6 on the frame's plane (S8e)."""
    return Case(name, 1e-7, 83, frame, 0.0, 0.0, boundaries, make=make)


def planar_cases():
    """(case, plane): S8e's sheets and wires, lines, arcs, circles, holes and
    S8b's splines, in both frames: planes crossing (at angles to the body's
    plane, through vertices, through a spline's knot point and joins, four
    times across a spline), along an edge (a notch's bottom between two
    prongs; a wire's edge between runs on both sides), tangent to an arc, a
    hole and a spline, missing, parallel off the body and containing it."""
    square = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    ushape = Boundary(points=[(0.0, 0.0), (9.0, 0.0), (9.0, 6.0), (6.0, 6.0), (6.0, 2.0), (3.0, 2.0),
                              (3.0, 6.0), (0.0, 6.0)])
    ell = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 2.0), (6.0, 2.0), (6.0, 6.0), (0.0, 6.0)])
    stadium = path([(0.0, -1.0), (3.0, -1.0), (3.0, 1.0), (0.0, 1.0)],
                   [None, arc(3.0, 0.0, 1.0), None, arc(0.0, 0.0, 1.0)])
    disc = Boundary(circle=(0.0, 0.0, 2.0))
    ring_hole = Boundary(circle=(5.0, 5.0, 2.0))
    slot = path([(-2.0, -1.0), (2.0, -1.0), (2.0, 1.0), (-2.0, 1.0)],
                [None, arc(2.0, 0.0, 1.0), None, arc(-2.0, 0.0, 1.0)])
    frame10 = Boundary(points=[(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)])
    p = spline_profiles()
    face, wire = 'face', 'wire'
    return [
        (planar('sheet_square_skew', [square], face), (4.0, 6.0, 2.0, 1.0, 2.0, 3.0)),
        (planar('sheet_square_diagonal', [square], face), (5.0, 5.0, 0.0, 1.0, 1.0, 0.0)),
        (planar('sheet_square_corner', [square], face), (10.0, 10.0, 0.0, 1.0, 1.0, 1.0)),
        (planar('sheet_u_along_edge', [ushape], face), (0.0, 2.0, 0.0, 0.0, 1.0, 0.0)),
        (planar('sheet_stadium_tangent', [stadium], face), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (planar('sheet_disc_chord', [disc], face), (1.0, 0.0, 0.0, 1.0, 0.0, 1.0)),
        (planar('sheet_disc_miss', [disc], face), (3.0, 0.0, 0.0, 1.0, 0.0, 0.5)),
        (planar('sheet_hole_tangent', [square, ring_hole], face), (0.0, 7.0, 0.0, 0.0, 1.0, 0.0)),
        (planar('sheet_slot_hole_tilted', [frame10, slot], face, TILT),
         framed(TILT, (1.0, 0.5, 0.0), (1.0, 0.5, 0.75))),
        (planar('sheet_parallel_off', [square], face), (0.0, 0.0, 1.0, 0.0, 0.0, 1.0)),
        (planar('sheet_in_plane', [square], face), (3.0, 4.0, 0.0, 0.0, 0.0, -1.0)),
        (planar('sheet_bulge_twice', p['bulge'], face), (10.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (planar('sheet_wave_knot', p['wave'], face), (6.5, 0.0, 0.0, 1.0, 0.0, 1.0)),
        (planar('sheet_lens_tilted', p['lens'], face, TILT), framed(TILT, (5.0, 4.0, 0.0), (0.25, 1.0, 0.5))),
        (planar('wire_square_skew_tilted', [square], wire, TILT), framed(TILT, (4.0, 6.0, 0.0), (1.0, 2.0, 3.0))),
        (planar('wire_square_diagonal', [square], wire), (5.0, 5.0, 0.0, 1.0, 1.0, 0.0)),
        (planar('wire_l_along_edge', [ell], wire), (0.0, 2.0, 0.0, 0.0, 1.0, 0.0)),
        (planar('wire_circle_chord_tilted', [disc], wire, TILT), framed(TILT, (0.5, 0.0, 0.0), (1.0, 1.0, 1.0))),
        (planar('wire_stadium_arc', [stadium], wire), (3.5, 0.0, 0.0, 1.0, 0.0, 0.25)),
        (planar('wire_parallel_off_tilted', [square], wire, TILT), framed(TILT, (0.0, 0.0, 1.0), (0.0, 0.0, 1.0))),
        (planar('wire_in_plane_tilted', [square], wire, TILT), (1.0, -2.0, 0.5, 0.0, 3.0, 4.0)),
        (planar('wire_blob_joins', p['blob'], wire), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (planar('wire_wave_four', p['wave'], wire), (0.0, 6.8, 0.0, 0.0, 1.0, 0.0)),
        (planar('wire_capsule_tilted', p['capsule'], wire, TILT), framed(TILT, (1.5, 0.0, 0.0), (1.0, 0.25, 2.75))),
        (planar('wire_bulge_tangent', p['bulge'], wire), (11.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
    ]


def _clip_moments(points, a, b, d, below):
    """Area and first moments (exact Fractions) of a polygon clipped to the
    half-plane `a u + b v + d <= 0` (or `>= 0`), by Sutherland-Hodgman and
    the shoelace sums (exact for one half-plane, whatever the polygon's
    shape)."""
    pts = [(F(x), F(y)) for x, y in points]
    g = lambda p: a*p[0]+b*p[1]+d
    keep = (lambda f: f <= 0) if below else (lambda f: f >= 0)
    out = []
    for i, p in enumerate(pts):
        q = pts[(i+1) % len(pts)]
        fp, fq = g(p), g(q)
        if keep(fp):
            out.append(p)
        if fp*fq < 0:
            t = fp/(fp-fq)
            out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    area = mx = my = F(0)
    for i, p in enumerate(out):
        q = out[(i+1) % len(out)]
        cross = p[0]*q[1]-q[0]*p[1]
        area += cross/2
        mx += (p[0]+q[0])*cross/6
        my += (p[1]+q[1])*cross/6
    return area, mx, my


def _quadratic_length(p0, p1, p2):
    """The closed-form length of a quadratic Bezier curve."""
    A = [ref.M(p1[i])-p0[i] for i in range(2)]
    C = [ref.M(p2[i])-2*ref.M(p1[i])+p0[i] for i in range(2)]
    # |B'(t)| = 2 sqrt(al t^2 + be t + ga).
    al = C[0]**2+C[1]**2
    be = 2*(A[0]*C[0]+A[1]*C[1])
    ga = A[0]**2+A[1]**2
    def antiderivative(t):
        f = ref.mp.sqrt(al*t*t+be*t+ga)
        return (2*al*t+be)*f/(4*al)+(4*al*ga-be*be)/(8*al**ref.M(1.5))*ref.mp.log(2*al*t+be+2*ref.mp.sqrt(al)*f)
    return 2*(antiderivative(ref.M(1))-antiderivative(ref.M(0)))


def planar_reference_checks():
    """S8e's reference against closed forms and exact sums: a square's
    halves (sheet: area, perimeter, centre; wire: length, centre) by a
    vertical plane and a diagonal one, a U's perimeters along its notch, a
    square with a tangent hole; every polygon sheet's sides against exact
    half-plane clipping in Fractions; every circle's sides against `r
    theta` arcs and circular segments; a square of straight splines against
    the polygon; spline sheets' whole area and moments against Green's
    theorem in Fractions and their sides' sums; every wire's sides' lengths
    against its whole length; and the bulge's quadratic Bezier against its
    closed-form length. Returns the largest relative deviation; raises
    beyond 1e-30."""
    M, mp = ref.M, ref.mp
    worst = [M(0)]
    def near(got, want):
        worst[0] = max(worst[0], abs(M(got)-M(want))/max(1, abs(M(want))))
    def same_rows(got, want):
        assert [r[0] for r in got] == [r[0] for r in want], (got, want)
        for g, w in zip(got, want):
            for x, y in zip([g[1], g[2], *g[3]], [w[1], w[2], *w[3]]):
                near(x, y)
    square = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    vertical, diagonal = (4.0, 0.0, 0.0, 1.0, 0.0, 0.0), (5.0, 5.0, 0.0, 1.0, 1.0, 0.0)
    third, root2 = M(10)/3, mp.sqrt(2)
    closed = [
        ('face', vertical, [('below', 40, 28, (2, 5, 0)), ('above', 60, 32, (7, 5, 0))]),
        ('face', diagonal, [('below', 50, 20+10*root2, (third, third, 0)),
                            ('above', 50, 20+10*root2, (2*third, 2*third, 0))]),
        ('wire', vertical, [('below', 18, 0, (M(8)/9, 5, 0)), ('above', 22, 0, (M(184)/22, 5, 0))]),
        ('wire', diagonal, [('below', 20, 0, (M(5)/2, M(5)/2, 0)), ('above', 20, 0, (M(15)/2, M(15)/2, 0))]),
    ]
    for make, plane, want in closed:
        same_rows(ref.planar_rows(planar('check', [square], make), plane), want)
    by_name = {case.name: (case, plane) for case, plane in planar_cases()}
    for name, want in (('sheet_u_along_edge', (22, 28)), ('sheet_hole_tangent', (34+4*mp.pi, 26))):
        rows = ref.planar_rows(*by_name[name])
        near(rows[0][2], want[0])
        near(rows[1][2], want[1])
    straight = path([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)], [
        spline(2, [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]),
        spline(3, [(10.0, 0.0), (10.0, 2.5), (10.0, 5.0), (10.0, 7.5), (10.0, 10.0)], (0.0, 1.0, 2.0)),
        None, None])
    for make in ('face', 'wire'):
        for plane in (vertical, diagonal, (4.0, 6.0, 2.0, 1.0, 2.0, 3.0), (10.0, 10.0, 0.0, 1.0, 1.0, 1.0)):
            same_rows(ref.planar_rows(planar('straight', [straight], make), plane),
                      ref.planar_rows(planar('polygon', [square], make), plane))
    for case, plane in planar_cases():
        p = ref.Planar(case, plane)
        rows = p.rows()
        segments = [s for b in case.boundaries for s in (b.segments or [])]
        lines_and_splines = all(b.circle is None for b in case.boundaries) and \
            not any(isinstance(s, tuple) for s in segments)
        if case.make == 'face':
            V, mu, mv, _ = p.prism.volume_moments(None)
            if len(rows) == 2:
                # The sides add up to the whole.
                parts = [p.prism.volume_moments(below) for below in (True, False)]
                for k, total in enumerate((V, mu, mv)):
                    near(parts[0][k]+parts[1][k], total)
            polygons = all(b.circle is None and b.segments is None for b in case.boundaries)
            if polygons and len(rows) == 2:
                for below, row in zip((True, False), rows):
                    area, mx, my = F(0), F(0), F(0)
                    for k, b in enumerate(case.boundaries):
                        s = 1 if k == 0 else -1
                        A_, X_, Y_ = _clip_moments(stored_points(b, case.tolerance), p.a, p.b, p.d, below)
                        area, mx, my = area+s*A_, mx+s*X_, my+s*Y_
                    got = p.prism.volume_moments(below)
                    for x, y in zip(got[:3], (area, mx, my)):
                        near(x, y)
            if lines_and_splines and any(isinstance(s, Spline) for s in segments):
                area, mx, my = ref.green_moments(dataclasses.replace(case, make=None))
                for x, y in zip((V, mu, mv), (area, mx, my)):
                    near(x, y)
        else:
            total = sum(q[1] for q in p.wire_pieces())
            near(sum(r[1] for r in rows), total)
        if all(b.circle is not None for b in case.boundaries):
            cx, cy, r = (M(x) for x in case.boundaries[0].circle)
            a, b, d = M(p.a), M(p.b), M(p.d)
            ab = mp.sqrt(a*a+b*b)
            h = -(a*cx+b*cy+d)/ab
            if abs(h) < r:
                al = mp.acos(h/r)
                n = (a/ab, b/ab)
                at = lambda s: p.world(cx+n[0]*s, cy+n[1]*s)
                if case.make == 'wire':
                    want = [('below', 2*r*(mp.pi-al), 0, at(-r*mp.sin(al)/(mp.pi-al))),
                            ('above', 2*r*al, 0, at(r*mp.sin(al)/al))]
                else:
                    seg = r*r*(al-mp.sin(al)*mp.cos(al))
                    arm = M(2)/3*r**3*mp.sin(al)**3
                    rest = mp.pi*r*r-seg
                    want = [('below', rest, 2*r*(mp.pi-al)+2*r*mp.sin(al), at(-arm/rest)),
                            ('above', seg, 2*r*al+2*r*mp.sin(al), at(arm/seg))]
                same_rows(rows, want)
    bulge = spline_profiles()['bulge']
    whole = ref.planar_rows(planar('bulge', bulge, 'wire'), (0.0, 0.0, 1.0, 0.0, 0.0, 1.0))
    near(whole[0][1], 26+_quadratic_length((10, 0), (12, 3), (10, 6)))
    if worst[0] > M(10)**-30:
        raise SystemExit(f'sheet and wire reference checks failed: {worst[0]}')
    return worst[0]


def primitive_line(kind, name, frame, params, plane):
    words = [kind, name, '1e-07'] + [repr(float(v)) for v in frame+params] + ['split'] + \
        [repr(float(v)) for v in plane]
    return ' '.join(words)


def encode(case, plane):
    text = encode_case(case)
    body, end = text.rsplit('\nend', 1)
    return body+'\nsplit '+' '.join(repr(float(v)) for v in plane)+'\nend'


def generate():
    blocks, out = [], ['# case\trow (S8a, split_reference.py)']
    frames = ['# case\taxis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case, plane in cases():
        blocks.append(encode(case, plane))
        for row in ref.rows(case, plane):
            out.append(f'{case.name}\t{ref.text(row)}')
        _, x, y, n = stored_axes(case.frame)
        for key, v in (('n', n), ('x', x), ('y', y)):
            frames.append(f'{case.name}\t{key}\t'+' '.join(struct.pack('>d', c).hex() for c in v))
    lines, prim = [], ['# case\trow (S8c, split_reference.revolved_rows)']
    for kind, name, frame, params, plane in primitive_cases():
        lines.append(primitive_line(kind, name, frame, params, plane))
        for row in ref.revolved_rows(kind, frame, params, plane):
            prim.append(f'{name}\t{ref.text(row)}')
    clines, conics = [], ['# case\trow (S8d.2, split_reference.revolved_rows)']
    for kind, name, frame, params, plane in conic_cases():
        clines.append(primitive_line(kind, name, frame, params, plane))
        for row in ref.revolved_rows(kind, frame, params, plane):
            conics.append(f'{name}\t{ref.text(row)}')
    tlines, tori = [], ['# case\trow (S8d, split_reference.torus_rows)']
    for name, frame, params, plane in torus_cases():
        tlines.append(primitive_line('torus', name, frame, params+(2*math.pi,), plane))
        for row in ref.torus_rows(frame, params, plane):
            tori.append(f'{name}\t{ref.text(row)}')
    slines, spirics = [], ['# case\trow (S8d.3, split_reference.torus_rows)']
    for name, frame, params, plane in spiric_cases():
        slines.append(primitive_line('torus', name, frame, params+(2*math.pi,), plane))
        for row in ref.torus_rows(frame, params, plane):
            spirics.append(f'{name}\t{ref.text(row)}')
    bblocks, bsplines = [], ['# case\trow (S8b, split_reference.py)']
    for case, plane in spline_cases():
        bblocks.append(encode(case, plane))
        for row in ref.rows(case, plane):
            bsplines.append(f'{case.name}\t{ref.text(row)}')
    pblocks, planars = [], ['# case\trow (S8e, split_reference.planar_rows: side S area perimeter cx cy cz '
                            'of a sheet, side S length 0 cx cy cz of a wire)']
    for case, plane in planar_cases():
        pblocks.append(encode(case, plane))
        for row in ref.planar_rows(case, plane):
            planars.append(f'{case.name}\t{ref.text(row)}')
    return {'split-cases.txt': '\n'.join(blocks)+'\n', 'split-expected.tsv': '\n'.join(out)+'\n',
            'split-sheet-cases.txt': '\n'.join(pblocks)+'\n',
            'split-sheet-expected.tsv': '\n'.join(planars)+'\n',
            'split-frames.tsv': '\n'.join(frames)+'\n',
            'split-primitive-cases.txt': '\n'.join(lines)+'\n',
            'split-primitive-expected.tsv': '\n'.join(prim)+'\n',
            'split-conic-cases.txt': '\n'.join(clines)+'\n',
            'split-conic-expected.tsv': '\n'.join(conics)+'\n',
            'split-torus-cases.txt': '\n'.join(tlines)+'\n',
            'split-torus-expected.tsv': '\n'.join(tori)+'\n',
            'split-spiric-cases.txt': '\n'.join(slines)+'\n',
            'split-spiric-expected.tsv': '\n'.join(spirics)+'\n',
            'split-spline-cases.txt': '\n'.join(bblocks)+'\n',
            'split-spline-expected.tsv': '\n'.join(bsplines)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    straight, green = reference_checks()
    planar_worst = planar_reference_checks()
    files = generate()
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                raise SystemExit(f'{path} is stale')
        else:
            path.write_text(contents)
    print(len(cases()), 'cases,', len(primitive_cases()), 'primitive cases,', len(conic_cases()),
          'conic cases,', len(torus_cases()), 'tori,', len(spiric_cases()), 'spiric tori,',
          len(spline_cases()), 'spline prisms')
    print('spline reference: straight splines within', ref.mp.nstr(straight, 3),
          'of the polygon, slicing within', ref.mp.nstr(green, 3), 'of Green (relative)')
    print(len(planar_cases()), 'sheets and wires; their reference within', ref.mp.nstr(planar_worst, 3),
          'of closed forms, exact clipping, Green and their sums (relative)')


if __name__ == '__main__':
    main()
