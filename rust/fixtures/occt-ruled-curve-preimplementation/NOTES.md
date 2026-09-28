# Two cones, and a cone's apex on another surface (S7b.4): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for these pairs existed
(`rust/kernel/src/intersection/ruled_curves.rs` absent). `inputs.txt` holds
the 11 cases of `ruled-curve-cases.txt` with the kernel's stored frame axes;
`oracle.cpp` is `occt_procedural_intersection_oracle.cpp` as captured, and
`native.txt` its rows: `GeomInt_IntSS` with tolerance `1e-7`, each line's
closedness and 17 samples, and each isolated point.

## Observations

* Against the independent reference (`ruled_curve_reference.py`: the first
  cone's rulings with `v = tan(t / 2)`, so the curve and its points at
  infinity are the zero set of a smooth function on a torus; the scan and
  components of `torus_curve_reference.py`; the apex and parallel twin cones
  by their factors), 9 cases agree: every native sample within 20 of the
  apex on both surfaces within `1e-6`, every component covered (cut 20 from
  the apex), the empty case empty. On `kk_unbounded` the native lines run to
  about `1.3e5`, their approximation errors growing along them (hundreds
  there); the comparison holds them to the surfaces near the apex.
* Reviewed (`occt-ruled-curve-divergences.json`): the apex lying on the
  sphere as an isolated point is missed, and the tiny loop beside the apex
  when the sphere is `2^-20` larger.

`compare_ruled_curves.py` requires every later run to reproduce these rows
(on another platform, its reviewed record) and compares the kernel's results
with both.
