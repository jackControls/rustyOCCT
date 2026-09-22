# Complete linear intersection contract

`intersection::linear_intersection` accepts all 25 ordered pairings of the
existing `LinearPrimitive3` points, lines, closed segments, planes and closed
triangles. Finite binary64 values define exact rational geometry. No tolerance,
snapping or normalization to rounded unit vectors changes that geometry.
Collapsed segments are singletons; collapsed lines and collinear plane/triangle
definitions are invalid. Nonfinite coordinates are errors.

The result is the **entire closed-set intersection**: empty, a point, a segment,
a convex polygon (three to six vertices), an infinite line or an infinite plane.
All construction coordinates remain rational, even beyond finite binary64.
Converting a construction to finite bounds is a separate fallible operation.
Rounded representatives do not certify incidence.

Results have a canonical representation independent of operand order and
defining-point permutations: lexicographically ordered segment endpoints;
polygon vertices in cyclic order beginning with the lexicographic minimum,
with the lexicographically smaller of the two cycle directions; a line with its
first nonzero direction coordinate equal to one and the corresponding origin
coordinate zero; a plane equation with its first nonzero normal coordinate one.
No duplicate or redundant collinear polygon vertices remain.

This is a geometric construction API, not a Boolean operation on B-rep solids.
It does not create edge/face identifiers, pcurves, topology or history maps.
The existing specialized intersection APIs retain their contracts.

Source baseline: `3d097a0328e71b826377d4814ab05ec3c3d23871`.
Read `IntAna_QuadQuadGeo::Perform(gp_Pln,gp_Pln)` (cross-normal direction,
parallel/coincident cases), `IntTools_EdgeEdge` (line incidence and interval
overlap), `IntTools_EdgeFace` / `IntCurveSurface_HInter` (surface intersections
within edge parameters), `BOPAlgo_BOP::BuildRC` (COMMON dimension filtering)
and `BOPAlgo_Section::PerformInternal1` / `BuildSection` (contact edges/vertices).
Rust uses independently implemented exact affine equations and halfspaces;
it does not port OCCT's floating tolerance or general B-rep algorithms.

The native harness captures COMMON and SECTION separately before the Rust
implementation is written. COMMON deliberately filters results below the
minimum input dimension; those omitted contacts are not Rust failures or
matching COMMON results. Both native outputs remain available for inspection.
Independent exact references, complete-set checks, invariants and sustained
mutation fuzzing are required in addition to native agreement.
