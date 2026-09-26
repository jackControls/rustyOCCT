# Spheres (S3): native observations before any kernel sphere code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
(`build_pinned_occt.py`). `capture.json` records the Rust revision and that
no kernel sphere code existed; the kernel tree had no changes. `inputs.txt`
is `sphere-cases.txt` (22 spheres and zones from
`generate_primitive_fixtures.py`), `oracle.cpp` is `occt_primitive_oracle.cpp`
as captured, `native.txt` its output.

## OCCT source review

* `BRepPrimAPI_MakeSphere(gp_Ax2, R, angle1, angle2)` builds a
  `BRepPrim_Sphere`, a `BRepPrim_Revolution` of the meridian circle about
  `-y` through `x` (`SetMeridian`, offset `2π`) between the latitudes
  `angle1` and `angle2`. The lateral face is a `Geom_SphericalSurface` on
  the frame; its seam lies along `x` (angle 0).
* An end at latitude `±π/2` (the binary64 value DRAW's `±90` degrees give)
  is a pole: one vertex and a degenerated edge. Any other end is a circle
  closed at its seam vertex, bounding a planar disc.
* So every sphere or zone has 2 vertices, 3 edges (the seam and two ends),
  `1 + caps` wires and faces, one shell and one solid: a full sphere
  2/3/1/1/1/1, a hemisphere 2/3/2/2/1/1, a zone 2/3/3/3/1/1.

## Observations

All 22 solids are BRepCheck-valid. Counts, faces (type, area, centre),
edges (degenerated, closed, length, point at the middle parameter) and
vertices equal `primitive_reference.py`'s, which derives them from the
specification alone. Volume, area, centre and inertia agree to `6.7e-15`
of the case's scale (bound `1e-9`). As on the cone, a degenerated edge
reports a length of about `7.7e-16` from `BRepGProp`, not 0.
