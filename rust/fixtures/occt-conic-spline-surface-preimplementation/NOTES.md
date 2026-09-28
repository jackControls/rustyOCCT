# Ellipses, hyperbolas and splines against analytic surfaces (S7c.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these curves existed
(`rust/kernel/src/intersection/conic_surface.rs` and `spline_revolved.rs`
absent). `inputs.txt` holds the 35 cases of `curve-surface-cases.txt` with
prefixes `e_`, `h_` and `s_`: ellipses and hyperbolas with the kernel's
stored frame axes (`Geom_Ellipse`, `Geom_Hyperbola`), and clamped rational
B-splines (`Geom_BSplineCurve`) against tori and cones; `oracle.cpp` is
`occt_curve_surface_oracle.cpp` as captured, and `native.txt` its rows:
`GeomAPI_IntCS`, each point with its curve parameter, each segment by its
parameters.

## Observations

* Against the independent reference (`curve_surface_reference.py`: the
  conics' rational parameterisations `tan(t / 2)` and `tanh(t / 2)`
  substituted exactly, the splines' exact span polynomials, sympy's
  `real_roots` with multiplicities; 80 digits for cones), 30 cases agree:
  every reference point matched within `1e-6` of the case's size and no
  other native point (a tangency may be reported as a cluster of points at
  it, as on `e_torus_tangent`).
* Reviewed (`occt-curve-surface-divergences.json`): a curve lying on the
  surface is never reported as a segment: nothing for an ellipse or a
  hyperbola in a plane or a spline on a cone's reference circle, and dozens
  of separate points along a spline on a torus's parallel.
* The rational-parameterisation reference reproduces the rows an earlier
  Groebner-basis draft of it gave on every conic case the two share.

`compare_curve_surface.py` requires every later run to reproduce these rows
(on another platform, its reviewed record) and compares the kernel's results
with both.
