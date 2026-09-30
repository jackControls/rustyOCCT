# Booleans of a sphere against cylinders where a frame is turned (S9d.2c): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-30 by
`compare_spheres_turned_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9d.2c's fixtures). `capture.json` records the Rust
revision and that the kernel's S9d.2c code did not exist
(`rust_spheres_turned_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/spheres_turned.rs` absent). The
kernel's probe reports every case `unsupported`: a turned hemisphere's
circle against a cylinder is refused with `OutOfDomain("a sphere's circle
of unequal axes against a cylinder (S9d.2b)")` (8 cases: `dome_tilt_pipe`,
`dome_lean_bite`, `dome_tiltx_rings`, `rim_tangent`), a sphere meeting a
cylinder in a turned frame in a loop with `OutOfDomain("a sphere meeting a
turned cylinder in a loop (S9d.2b)")` (10: `bite_lean`, `bite_tilt`,
`graze_tiltx`, `bite_r125`, and `dome_bite_turned`, whose loop is refused
before its cap's circles), the messages read from `protocol::run` by a
throwaway example. No protocol or oracle change was needed. `inputs.txt`
holds the 18 cases of `boolean-spheres-turned-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
sphere or hemisphere is `BRepPrimAPI_MakeSphere(gp_Ax2, R, low, high)`, each
prism `BRepPrimAPI_MakePrism` of its circle in its own frame. Frames:
hemispheres in `TILT`, `LEAN` and `TILTX` against prisms in `XY`; prisms
in `LEAN`, `TILT`, `TILTX` and `R125` against whole spheres; a hemisphere
in `TILT` against a prism in `LEAN`.

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 18 (`spheres_boolean_reference.py`,
  each result one solid).
* 1 case matches within the 2e-8 allowance: `dome_tilt_pipe_common` (the
  pipe's circles and the tilted disc's ellipse, no B-spline edge) within
  5.9e-9.
* 17 cases are reviewed (`occt-boolean-spheres-turned-divergences.json`),
  all `measure` (and `centre` but for `dome_tilt_pipe_cut` and
  `dome_lean_bite_cut`): BRepGProp's default integration on faces bounded
  by the Boolean's approximated intersection curves misses by up to 1.2e-5
  in volume (`rim_tangent_common`), 1.4e-5 in area and 1.8e-6 of the case's
  size in centre (`graze_tiltx_cut`). The same native results measured by
  BRepGProp with an adaptive accuracy of 1e-10 (a diagnostic build of the
  oracle, not committed; unchanged at 1e-12, so converged) are within
  1.5e-8 in volume, 1.1e-8 in area and 8.6e-9 in centre but for
  `dome_lean_bite_common` (4.0e-8 in volume and 2.1e-8 in area of a
  result of volume 0.26: 1.0e-8 absolute, against its area times its edges'
  tolerance, 2.9e-7): their edges' tolerances at most 7.5e-7 (`bite_tilt`),
  up to 4 of their edges B-splines. The geometry agrees with the reference
  within the native tolerances; the default measure does not reach the
  allowance.
* The degenerate fixture is a valid result natively: `rim_tangent_common`
  one solid (the hemisphere's rim touching the cylinder's wall, a crossing
  of two points 4.4e-22 of the radius squared from tangency in the
  reference's exact model; reviewed for its measure). The decisions refuse
  it (`Degenerate`).
* Five solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `dome_tilt_pipe_cut` 6/12/8 to
  4/7/5, `dome_tilt_pipe_common` 5/12/8 to 4/10/7, `bite_tilt_common` 3/3/2
  to 2/2/2, `graze_tiltx_common` 3/3/2 to 2/2/2 and `dome_bite_turned_cut`
  4/6/4 to 3/4/3 (faces of one surface split at seams merged). The
  kernel's results will be compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_spheres_turned_boolean_fixtures.py --check`): closed forms (a
  whole sphere against a turned cylinder by S9d.2's lens of two discs along
  its axis; the hemisphere against the coaxial pipe by its sections along
  the pipe's axis, circular segments of one quadrature, its faces by
  central symmetry and the disc inside the pipe a circle and an ellipse)
  within 8.1e-17 (ideal axes against stored ones), a cap and its
  complement against the whole sphere sliced along the prism's axis within
  1.7e-41, inclusion and exclusion 1.7e-41, the area identity 1.4e-40,
  every face's classes 9.2e-41, a second slicing direction (whole spheres)
  1.4e-40, Monte Carlo 2.2 standard errors; its scan for near coincidences
  found none but in the declared degenerate pair, and every other cap
  circle is at least 0.16 of the radius squared from tangency.
