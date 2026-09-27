# Spline cells (S4a): native observations before any kernel code certifies spline geometry

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code certified spline geometry: the kernel had
the spline variants and the R4 continuity check, and every other check of a
spline reported its uncertified issue. `inputs.txt` holds the explicit OCCT
constructions of `generate_brep_fixtures.spline_models()`, `oracle.cpp` is
`occt_brep_check_oracle.cpp` as captured (B-spline curves, 2D curves and
surfaces added), and `native.txt` its three rows per case: BRepCheck
statuses, distinct subshape counts, and tolerances with OCCT's own measured
gaps and deviations.

## The models

The independent prism builder gained a spline side: a clamped B-spline from
one profile point to the next, its edges at both caps, its pcurves on the
caps (the same poles in the cap's frame), and its wall the ruled spline
surface `C(u) + v z` over `v` in `[0, h]`, whose pcurves are lines. Every
relation is exact, so a kernel that composes the rational pieces must find
no deviation at all.

* `spline_bulge` (a quadratic side), `spline_cubic_bulge` (a cubic with a
  simple knot), `spline_c0_bulge` (a quadratic with a knot of multiplicity
  2 that is a corner: every relation still exact) and
  `spline_rounded_corner` (the rational quadratic of a quarter circle).
* Mutations of the bulge: a cap's spline pcurve shifted by `1e-3`, a vertex
  moved by `1e-3`, and the wall reversed.
* A stadium whose vertical edge's pcurve on the cylinder is a degree-1
  spline (`spline_stadium_pcurve`), the same shifted by `1e-3` in `u`, and
  the edge itself a degree-1 spline (`spline_stadium_edge`): spline geometry
  on a periodic analytic surface, the Taylor path of S4b.

## Observations

* The four sides and the two stadium variants are BRepCheck-valid, with the
  counts of their seamed encodings (8 12 6 6 1 1, or 10 15 7 7 1 1 for the
  rounded corner). `spline_c0_bulge` is valid natively: BRepCheck has no
  continuity status, so the kernel's `edge_not_c1`, `pcurve_not_c1` and
  `face_not_c1` there have no native counterpart.
* The shifted pcurves report `BRepCheck_InvalidCurveOnSurface` (8) and
  `BRepCheck_InvalidSameParameterFlag` (11) on the edge and an unclosed wire
  (27) on the face, measuring `1.00001e-3` and `1.0000099583e-3`. The moved
  vertex reports `BRepCheck_InvalidPointOnCurve` (1), measuring
  `9.9999999999989e-4`. The reversed wall reports a bad orientation of a
  subshape (32) on its face.
* Every other measured deviation is OCCT's evaluation rounding, at most
  `9.2e-16` (the rational corner's wall), and every vertex gap at most
  `1.8e-16` (the stadium's circles). Tolerances are the requested `1e-7`.

A kernel enclosure of a valid spline model must not be below these
measurements (up to `2^-46` of the case's size, as for M5) nor above the
tolerance. `compare_brep.py --family spline` requires every later run to
reproduce these rows exactly, the measurements within that allowance.
