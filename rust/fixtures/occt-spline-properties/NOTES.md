# Spline mass properties (S4d): native observations before any kernel code integrates them

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code computed mass properties with a spline
surface or a spline pcurve on a cylinder or cone: those were `None`, and
only spline pcurves on planes, spheres and tori were integrated. The inputs
are the spline models of the S4a capture (`inputs.txt`, identical to
`occt-spline-preimplementation/inputs.txt`); `oracle.cpp` is
`occt_brep_check_oracle.cpp` as captured, with `BREP_ORACLE_PROPERTIES` set,
and `native.txt` its fourth row per case: `BRepGProp::VolumeProperties` and
`SurfaceProperties` at a requested relative precision of `1e-12`, as volume,
its relative error estimate, area, its relative error estimate, centre of
mass, and the matrix of inertia about the centre (row-major).

## Observations

* The volumes are the profiles' areas times the height: `20/3` for the
  quadratic bulge (the rectangle `3 × 2` and the parabolic segment's
  `2/3`), `5.75 + π/16` for the rounded corner, `6 + π/2` for the stadiums
  (the spline pcurve and the degree-1 spline edge change nothing).
* Every error estimate is at most `4.5e-16`, but the estimates are not
  bounds: the quadratic bulge's area is `40/3 + 8 + √2 + asinh 1 =
  23.628920482725971` (the caps, the three flat walls and the arc length of
  the spline side), and OCCT reports `23.6289205238958`, `1.7e-9` relative
  above it; the rounded corner's volume is `5.75 + π/16 =
  5.9463495408493621` and OCCT's is `5.7e-12` relative below it. The
  cubic bulge, the stadiums and every volume of a nonrational wall agree
  to rounding.
* The reversed wall has the same volume but an area smaller by twice the
  wall's (`2 × 2.2956`): `SurfaceProperties` counts the reversed face
  negatively. The kernel reports that case invalid (`loop_winding`), so it
  is never compared.
* The shifted pcurve and moved vertex change the properties by less than
  `1e-3` relative; they are invalid on both sides and not compared either.

The bridge (`compare_brep.py --family spline`) requires these rows to
reproduce within `1e-9` relative plus OCCT's error estimates, and, for every
model valid on both sides, each kernel enclosure (volume, area, centroid,
inertia) to contain OCCT's value within OCCT's error estimate plus `1e-8`
relative, the allowance the bulge's area needs.
