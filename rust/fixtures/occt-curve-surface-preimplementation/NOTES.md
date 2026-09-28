# Lines and circles against analytic surfaces (S7c.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these pairs existed
(`rust/kernel/src/intersection/curve_surface.rs` absent). `inputs.txt` holds
the 36 cases of `curve-surface-cases.txt` with the kernel's stored frame
axes (a line by a point and a direction); `oracle.cpp` is
`occt_curve_surface_oracle.cpp` as captured, and `native.txt` its rows:
`GeomAPI_IntCS`, each point with its curve parameter, each segment by its
parameters.

## Observations

* Against the independent reference (`curve_surface_reference.py`: exact
  polynomials along lines and exact Groebner bases for circles by sympy,
  80-digit roots for cones), 29 cases agree: every reference point matched
  within `1e-6` of the case's size and no other native point (a tangency may
  be reported as a cluster of points at it, as on `c_torus_tangent`).
* Reviewed (`occt-curve-surface-divergences.json`): a curve lying on the
  surface is never reported as a segment: nothing for a line in a plane or
  along a cylinder, a circle in a plane, on a sphere or on a cylinder, and
  hundreds of separate points along a torus's meridian or parallel.

`compare_curve_surface.py` requires every later run to reproduce these rows
(on another platform, its reviewed record) and compares the kernel's results
with both.
