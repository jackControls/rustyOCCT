# Procedural intersections of cylinders and spheres (S7b.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these pairs existed
(`rust/kernel/src/intersection/procedural.rs` absent). `inputs.txt` holds the
20 cases of `procedural-intersection-cases.txt` (two cylinders with crossing
axes, a cylinder and a sphere off its axis) with the kernel's stored frame
axes; `oracle.cpp` is `occt_procedural_intersection_oracle.cpp` as captured,
and `native.txt` its rows: `GeomInt_IntSS` (the engine under
`GeomAPI_IntSS`, which also reports isolated points) with tolerance `1e-7`,
each line's closedness and 17 samples, and each isolated point.

## Observations

* Against the independent reference (`procedural_intersection_reference.py`,
  exact classes and 80-digit roots and points), 17 cases agree: every native
  sample within `1e-6` of the case's size from the exact curve, every
  component (a loop, each ring, a figure-eight) carrying native samples, and
  the empty and tangent-point cases matched. OCCT splits a closed curve into
  several lines (two to four here); the comparison is by component, not by
  line.
* Reviewed (`occt-procedural-intersection-divergences.json`): a sphere
  tangent to a cylinder from outside or from inside gives no result natively
  (the isolated tangency is missed), and near a Viviani configuration one of
  OCCT's 51 approximated samples is `3.3e-5` from the exact curve, at the
  loop's nearly pinched neck.

`compare_procedural_intersections.py` requires every later run to reproduce
these rows and, once the kernel intersects these pairs, compares its
certified curves with both.
