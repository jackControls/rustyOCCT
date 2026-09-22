# Native tensor patch editing observations

The 691 native inputs were captured from OCCT 7.9.3 before the Rust surface
editing implementation. They cover 4,023 output patches and 324 isocurves.
Every complete Rust result matches independent exact Cox tensor polynomials
and binomial affine substitutions. The native bridge currently has **637
matching cases and 54 reviewed differences** on macOS OCCT 7.9.3. All output
counts, types, degrees and original U/V ranges agree. Native runtime versions
are separate from the pinned source reference; this is not a runtime build of
`3d097a0328e71b826377d4814ab05ec3c3d23871`.

All 54 differing cases involve a degree-25 direction and an operation that
calls `Geom_BezierSurface::Segment`: 12 rectangular trims, 12 U splits, 12 V
splits and 18 composed edits. The budget remains `1e-10 + 2e-12*abs(exact)`
per domain, Euclidean pole or normalized-weight component. The bridge permits
one common weight scale per item, not independent pole rescaling. No tolerance
was widened, and reviewed native results never exempt Rust from the exact
homogeneous coefficient contract.

For `d25_25_p00_w0_op2`, the last patch's upper corner must retain the original
clamped surface corner `(24,25,5)`. The native U split returns
`(24.000000211690278,25.000000239022278,5.000000221933909)`.
The V split (`...op3`) returns
`(23.999999962659242,24.999999950903764,4.999999995247982)`.
An endpoint of a clamped rational Bézier patch equals that Euclidean control
point, so these examples demonstrate changes in geometry, not just an
alternative projective representation.

The source path is `Geom_BezierSurface::Segment` → `BSplSLib::BuildCache` →
`PLib::{UTrimming,VTrimming,CoefficientsPoles}`. It edits floating power
coefficients and reconstructs Bernstein controls. The observed high-degree
errors are consistent with cancellation along that conversion path; the
individual floating instruction responsible has not been isolated. Rust uses
exact homogeneous de Casteljau restriction and retains rational controls.

`fixtures/occt-surface-editing-divergences.json` pins the native runtime, input
SHA-256, exact-output SHA-256, every differing field and a SHA-256 covering all
native binary64 domain/pole/weight bits, types and degrees. Consecutive field
indices are encoded as ranges without dropping fields. Raw native observations
and complete exact Rust outputs remain in each CI artifact. Any changed pin,
new discrepancy or exact Rust error fails the comparison. `--strict-native`
also fails these reviewed numerical differences. Native exceptions and invalid
control values are preserved as failures, rather than silently filtered out.

The extra source GTest input families are adaptations, not new unchanged DRAW
passes. Arbitrary face trimming, sewing, intersection curves, generic topology
history, Booleans and STEP remain outside this capability.
