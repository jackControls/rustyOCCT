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
* Every error estimate is at most `4.5e-16`: the integrands are polynomial
  on the ruled spline walls, and OCCT's Gauss integration of them is exact
  up to rounding.
* The reversed wall has the same volume but an area smaller by twice the
  wall's (`2 × 2.2956`): `SurfaceProperties` counts the reversed face
  negatively. The kernel reports that case invalid (`loop_winding`), so it
  is never compared.
* The shifted pcurve and moved vertex change the properties by less than
  `1e-3` relative; they are invalid on both sides and not compared either.

The bridge (`compare_brep.py --family spline`) requires these rows to
reproduce within `1e-9` relative plus OCCT's error estimates, and, for every
model valid on both sides, each kernel enclosure (volume, area, centroid,
inertia) to contain OCCT's value within the same allowance.
