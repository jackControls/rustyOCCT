# Splines against circles, ellipses and hyperbolas (S7d.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these pairs existed
(`rust/kernel/src/intersection/spline_curve.rs` absent). `inputs.txt` holds
the 13 cases of `curve-curve-cases.txt` whose names start with `s`: a clamped
rational B-spline (`Geom_BSplineCurve`, its whole domain) first, a circle, an
ellipse or a hyperbola second, with the kernel's stored frame axes;
`oracle.cpp` is `occt_curve_curve_oracle.cpp` as captured, and `native.txt`
its `IntTools_EdgeEdge` common parts.

## Observations

* Against the independent reference (`curve_curve_reference.py`,
  `spline_rows`: the spline's exact span polynomials substituted into the
  conic's implicit equations, sympy's gcd and `real_roots`, tangency by
  parallel tangent directions), 11 cases agree, including the exact rational
  quarter circle reported as an edge part.
* Reviewed (`occt-curve-curve-divergences.json`): the same quarter circle as
  the first span of a two-span spline gives two separate vertices on it and
  no edge part (`sc_partial`); a straight spline touching an ellipse's
  vertex gives two vertices 5.5e-5 either side of the tangency
  (`se_tangent`).

`compare_curve_curve.py` requires every later run to reproduce these rows
(on another platform, its reviewed record) and compares the kernel's results
with both.
