# A torus with a cylinder or a cone (S7b.3b.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these pairs existed
(`rust/kernel/src/intersection/torus_curves.rs` absent). `inputs.txt` holds
the 24 cases of `torus-curve-cases.txt` with the kernel's stored frame axes;
`oracle.cpp` is `occt_procedural_intersection_oracle.cpp` as captured, and
`native.txt` its rows: `GeomInt_IntSS` with tolerance `1e-7`, each line's
closedness and 17 samples, and each isolated point.

## Observations

* Against the independent reference (`torus_curve_reference.py`: every
  meridian's roots from 80-digit polynomial roots, the critical meridians
  from a resultant, components by continuity, tangencies by exact algebra),
  18 cases agree: every native sample within `1e-6` of both surfaces, every
  component of the reference curve carrying native samples, the empty cases
  empty.
* Reviewed (`occt-torus-curve-divergences.json`): the three cases of isolated
  tangency points give no result natively; one line near a nearly singular
  saddle is `4.9e-6` off, and on two cone cases single approximated lines
  are up to `2.6e-4` off the cone (walking-line approximation), their end
  points on both surfaces.

`compare_torus_curves.py` requires every later run to reproduce these rows
and compares the kernel's results with both.
