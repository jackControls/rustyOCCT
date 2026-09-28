# Intersections with a torus (S7b.3a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel torus intersection code existed
(`rust/kernel/src/intersection/toroidal.rs` absent). `inputs.txt` holds the
58 `tp_`, `ts_` and `tx_` cases of `procedural-intersection-cases.txt` (a
torus and a plane, a torus and a sphere, coaxial pairs) with the kernel's
stored frame axes; `oracle.cpp` is `occt_procedural_intersection_oracle.cpp`
as captured (planes and tori added), and `native.txt` its rows:
`GeomInt_IntSS` with tolerance `1e-7`, each line's closedness and 17
samples, and each isolated point.

## Observations

* Against the independent reference (`procedural_intersection_reference.py`
  with the meridian parameterisation, exact classes, 80-digit loop ends and
  points), 55 cases agree: every native sample within `1e-6` of the case's
  size from the exact curve or circles (by closest-point distance, loop ends
  found by bisection on the sign of `D`), every component carrying native
  samples, the empty, same and not-yet-parameterised cases matched (for the
  last, the samples lie on both surfaces). OCCT splits a closed curve into
  several lines, and returns the special and coaxial cases as circles.
* Reviewed (`occt-procedural-intersection-divergences.json`): a plane
  touching the outer equator, a sphere in the hole touching the inner
  equator and a sphere touching the outer equator from outside give no
  result natively (the isolated tangency is missed).
* The first capture of these cases used a plane normal `(0, 1, sqrt 15)`
  for the near-Villarceau case, whose stored normal the reference's
  emulation of `Frame3::new` and the kernel round one unit in the last place
  apart; the case now uses the Pythagorean normal `(0, 3, 4)` with a torus of
  ratio near `3/5`, the fixture checks stored normals bit for bit, and this
  capture was retaken with the kernel's torus code set aside.

`compare_procedural_intersections.py` requires every later run to reproduce
these rows and compares the kernel's results with both.
