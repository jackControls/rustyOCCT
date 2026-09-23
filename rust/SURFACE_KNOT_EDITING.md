# Exact B-spline surface knot editing

The surface representation retains two exact rational knot vectors and a
U-major grid of positive homogeneous controls `(wx, wy, wz, w)`. Each axis
supports degrees 1 through 25 and independent periodicity. At most 4096
controls are permitted in the complete tensor grid. Conversion from existing
binary64 surfaces or exact Bezier patches preserves every represented value.

Refinement requests total multiplicities. Unordered requests, duplicate knots
and the two periodic seam aliases use the same contract as exact curve knot
editing. A combined U/V batch validates both axes and the complete resulting
grid before control arithmetic. Invalid requests are errors even if another
request is a no-op. All operations are immutable and atomic.

Removal succeeds only if every transverse homogeneous control curve has an
exact positive-weight representation in the candidate basis. A valid but
inexact removal returns no surface. No tolerance, sampled deviation test or
partial row update is used. Deleting the last occurrence of a periodic origin
advances that axis's origin and preserves its period. The untouched axis and
the complete homogeneous tensor function are preserved, including raw support
outside a nonperiodic unclamped fundamental domain.

Exact evaluation retains position and all partials through total order two,
with explicit quadrant selection and exact continuity decisions. Isocurves
retain an exact B-spline representation and can feed the existing certified
intersection APIs. Rational rectangle extraction returns complete exact Bezier
patches with the existing Cartesian patch/control limits. Binary64 enclosures
are optional and can fail without discarding exact results.

These count limits bound storage and iteration counts, not arbitrary-precision
operand sizes or wall-clock time. General surface degree elevation, arbitrary
face trim loops, surface/surface intersection and B-rep topology remain separate
capabilities.

Source reference: OCCT `3d097a0328e71b826377d4814ab05ec3c3d23871`:
`Geom_BSplineSurface::{InsertUKnots,InsertVKnots,RemoveUKnot,RemoveVKnot}`,
`BSplSLib::{SetPoles,GetPoles,InsertKnots,RemoveKnot}`, and the called
`BSplCLib` insertion/inverse-insertion routines. OCCT packs transverse controls
as a higher-dimensional curve; Rust shares a sparse exact Boehm insertion map
across homogeneous rows/columns and uses the existing exact curve inverse for
removal. Native removal uses an explicit tolerance and
can differ from exact homogeneous removal. Differences must be reported and
independently justified, not silently treated as native matches.

The original `Geom_BSplineSurface_Test.cxx` SetUp and U/V knot insertion/removal
fixtures supply comparison inputs. These are adapted native API comparisons,
not unchanged upstream GTest passes. Native observations are captured before
the Rust implementation. Independent Cox coefficient equations verify
all 1,756 complete control grids and removal decisions, including eight extreme
range cases beyond the native corpus. Separate Rust full-support coefficient
equations, tensor identities and closed quotient formulas check operation
sequences and every structured fuzz input. All 791 existing surface jet fixtures
also compare the new retained representation with the prior exact de Boor path.

The evaluator constructs differentiated de Boor pole weights once per parameter
axis, then applies exact integer dot products to every homogeneous row. This
preserves the same derivatives and quadrant contract while avoiding repeated
arbitrary-precision interpolation through each row. The independent checker
shares Cox matrices and fraction-free elimination across transverse fields,
while retaining every coefficient equation and every residual check.


Refinement builds its sparse map on unit controls once per axis, using the same
five-period unrolling/cropping convention as exact curves. Integer dot products
then apply that map to the entire grid. This avoids repeated control rounding
and repeated fraction reductions without approximating any tensor coefficient.
The sustained surface-knot campaign uses an explicit 60-second instrumented
verification budget per input and the common 2 GiB process limit; those are
harness resource bounds, not production latency guarantees.
