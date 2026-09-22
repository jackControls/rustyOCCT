# Native spline/quadric observations

The 92-case corpus captures `GeomAPI_IntCS` against sphere and infinite-cylinder
surfaces before the Rust implementation. It includes rational Bézier and
B-spline curves through degree 25, exact contacts, contained arcs/generators,
unclamped knots, corners, explicit trims and multiple periodic turns. The
source reference is `3d097a0328e71b826377d4814ab05ec3c3d23871`; local OCCT 7.9.3
and the distribution runtime used by Linux CI are separate builds.

Local OCCT 7.9.3 has **66 matches and 26 reviewed differences**. Every Rust
answer also matches an independent exact basis-polynomial and continued-fraction
calculation. The differences are:

- Native OCCT omits overlap intervals for constant points on the surface,
  exactly rational quarter-circle arcs, and cylinder generators. These have an
  identically zero implicit polynomial on the entire queried interval.
- For partial overlap, native OCCT omits the interval and can return its endpoint
  as an isolated point. Rust returns the maximal closed overlap and excludes its
  endpoints from isolated events.
- A periodic polyline queried on `[3.5,8.5]` has six exact parameter events,
  including both endpoints. Native OCCT returns four, ending at 6.5. The source
  `ComputeAppendPoint` normalizes periodic parameters into one period; Rust
  retains all events in the explicitly requested interval.

`occt-spline-quadric-divergences.json` pins each runtime version, complete input
hash, every native output bit and the independent certificate hash. The shared
comparison rejects an altered Rust result even when a native review exists.
All native points and intervals remain in the report; no numerical budget was
widened and reviewed cases are not counted as parity matches. `--strict-native`
fails on these reviewed differences as well.

Run `target/math-oracle-venv/bin/python rust/tools/compare_spline_quadric.py`
against an installed SDK. The ordinary 118 exact fixtures additionally cover
dense degree-25 equations, order-50 contacts and extreme exponents. Neither
corpus establishes general surface intersection or full OCCT equivalence.
