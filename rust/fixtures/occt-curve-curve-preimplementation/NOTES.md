# Pairs of lines, circles, ellipses and hyperbolas (S7d.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these pairs existed
(`rust/kernel/src/intersection/curve_curve.rs` absent). `inputs.txt` holds
the 41 cases of `curve-curve-cases.txt` with the kernel's stored frame axes
(a line by a point and a direction); `oracle.cpp` is
`occt_curve_curve_oracle.cpp` as captured: `IntTools_EdgeEdge` on edges of
the two curves (lines over `[-10, 10]`, circles and ellipses whole,
hyperbolas over `[-3, 3]`, every reference point inside those ranges), and
`native.txt` its common parts: vertices with both parameters, edge parts by
their range on the first edge.

## Observations

* Against the independent reference (`curve_curve_reference.py`: both
  curves' implicit equations solved together by sympy's Groebner basis,
  coincidence by reduction modulo the first curve's basis, tangency by
  parallel tangent directions), all 41 cases agree: every reference point
  matched by a native vertex within `1e-6` of the case's size, no other
  native part, and every coincidence an edge part.

`compare_curve_curve.py` requires every later run to reproduce these rows
(on another platform, its reviewed record) and compares the kernel's results
with both.
