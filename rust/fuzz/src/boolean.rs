//! Booleans of two prisms in one frame (S9a of REVIEW_NOTES.md): the split
//! target's line and arc profiles (rectangles, regular polygons, stadiums,
//! U shapes, squares with round or square holes) on dyadic sizes, the tool
//! offset by a dyadic vector in the axis-aligned frame (its origin moved) or
//! sharing the tilted frame's origin, its heights equal to the object's,
//! spanning them, overlapping them, disjoint from them, inside them or on
//! them, their profiles lines and arcs or the split target's splines
//! (S9a.2). Fuse, cut and common never panic and fail only as documented (a
//! cavity among several solids, S9a.2's, or arcs in frames with different
//! axes, S9c's; S9b turns, leans or tilts the tool's frame, S9c.2a stands it
//! on its side (perpendicular cylinders), S9d.1 makes it a sphere, a cap or
//! a zone (S9d.2c, with no decode change, meets a turned cap's circles
//! against the object's arcs and a sphere in the shared tilted frame in
//! loops, S9d.2b's refusals), S9d.3a a cone or frustum, S9d.4a a whole
//! torus, S9d.4b.1 a torus v-segment or wedge, S9d.4b.2a makes the object a
//! sphere or a cone against a whole torus, S9d.4b.2b a whole torus, S9d.3c
//! a sphere against a cone or sphere tool and a cone against a sphere, cap
//! or zone tool or a cone tool (`CONE_PAIRS`), S9d.4c a sphere,
//! cap, zone or cone against a torus v-segment or wedge tool and a cap or
//! zone against a whole torus; a
//! result thinner than the resolution or touching itself; an undecided
//! comparison); each result
//! validates as it is built and its history passes the independent check
//! (debug builds); when all three succeed their volumes agree,
//! `V(A ∪ B) = V(A) + V(B) - V(A ∩ B)` and `V(A - B) = V(A) - V(A ∩ B)`
//! (two tori: one operation, its volume within its bounds);
//! every result moves rigidly with its ids, and each of its vertices
//! classifies on its boundary; each operation's first result is an input
//! again (S9b.2) against a turned box, with the same identities; one with
//! curved faces (S9e.1: a result of prisms with arcs in frames with
//! different axes) too, on its construction's exact model (`GIVEN_CURVED`);
//! S9e.2: a stack with arc walls and one solid of several are given too,
//! and the chained byte's upper bits may turn the box into a turned
//! cylinder (`GIVEN_ROUND`: a stack or an S9b.1 result given with arcs);
//! S9e.3a: first results of spheres, cones and tori are given too, and the
//! chained byte's next bit may make the partner a sphere (`GIVEN_BALL`: a
//! given result against a sphere); S9e.3b: their meetings of two curved
//! faces met by the partner, which the chained byte's next bit centres on
//! the first such edge (`GIVEN_MET`). S9f.1: a spline prism
//! against a line prism in a turned, leaning, tilted or side frame is
//! decided by the curved engine's spline walls (`SPLINE_WALLS`; before it,
//! refused), and the byte after the chained one's place, at or above 128,
//! makes the object's spline profile R4's knot (`knot_profile`: a quadratic
//! whose interior knot of multiplicity two is C1 exactly, removed before
//! lifting). S9f.2a: a spline prism against a prism with arcs or splines
//! in the frame turned about the axis meets it on an exactly parallel
//! axis (`SPLINE_PARALLEL`), and in the tilted frame the byte after R4's,
//! at or above 128, offsets a spline variant's tool by an amount that
//! rounds (the curved engine instead of S9a.2's one frame). S9f.2b: a
//! spline prism against a prism with arcs on crossing axes (the tool
//! leaning, tilted or on its side) meets it along the walls' meetings with
//! the cylinders (`SPLINE_CROSSING`), their loops round a turning point
//! inside the faces and a cap circle on a wall in a plane along its axis
//! too (S9f.2b.2); spline walls against spline walls on crossing axes, a
//! torus, a cone or a given result stay refused (S9f's). S9f.3a: a spline
//! prism object against the sphere, cap or zone tool (the spline byte from
//! 192, its low bit set, the flags' bits 5 and 6 keeping the object's
//! prism) meets it along the walls' meetings with the sphere, their loops
//! and the rims' and split's circles on the walls (`SPLINE_SPHERE`).
//! S9f.3b: a spline prism object against the cone or frustum tool (the
//! spline byte in 160..192, its low bit set, the flags' bits 5 and 6
//! keeping the object's prism) meets it along the walls' meetings with the
//! cone, both nappes or its loops, and the rims on the walls
//! (`SPLINE_CONE`).
//! S9e.4a: the object written by the kernel's `.brep` writer and read back
//! is imported (`Solid::imported_with`) and given the chosen operation with
//! the tool again, its volume the object's own result's (`IMPORTED`).
//! S9e.4b.2: the chained stage's first result of plane faces likewise,
//! imported (a polyhedron on its stored vertices where it is no prism) and
//! cut by the turned box again, its volume the chained cut's. S9e.4b.3a:
//! a first result of one sphere, cylinder or cone face and plane faces too
//! (a plane piece, its primitive common its planes' half-spaces).
//! S9e.4b.4b.1: one of several sphere, cylinder and cone faces too (a
//! Boolean chain of their primitives; S9e.4b.4b.2a: one led by a prism leaf,
//! its walls' arcs tangent to its lines or not; S9e.4b.4b.2b.1: one with
//! other plane faces, a primitive's flat, a prism of one cap or a pocket).
//! S9e.4b.4a: a lens prism in the object's frame given the chosen operation
//! with the tool, and again once written, read back and imported, its arcs
//! of two circles taken through their ends where their joints round off
//! both (`JOINTS`). S9e.4b.4c.1: the chained cut's first solid of plane
//! faces written, read back and imported (a polyhedron on its stored
//! vertices in the curved engine) and the kernel's own cut each less a ball
//! about the turned box's axis, the same volume (`MESHES`). S9e.4b.4c.2a:
//! the kernel's own cylinder on the object's frame less a prism of a square
//! with a square or a U-shaped hole from inside it past its top (a post or
//! a U island standing in a pocket, pockets two or three deep) and its
//! import, each cut by the turned box, the same volume (`NESTED`).
//! S9e.4b.4b.2b.2a: by the input's 17th byte (past those it decodes), the
//! kernel's own box on the object's frame fused with a stadium through both
//! its caps (two prism leaves, the stadium's walls split by the box) or a
//! stadium plate less a square turned about its axis from inside it past its
//! top (a lone leaf with a pocket), and its import, each cut by the chained
//! stage's partner, the same volume (`LEAVES`). S9e.4b.4b.2b.2b.1: by the
//! input's 18th byte, the kernel's own stadium plate less a square with a
//! round hole across its rim from inside it past its top (a leaf cut at its
//! rim with a post) or a stadium plate fused with a square with a square hole
//! over its arc end (a well), and its import, each cut by the chained stage's
//! partner, the same volume (`RIMS`). S9e.4b.4b.2b.2b.2: by the input's 19th
//! byte, the kernel's own stadium plate fused with a diamond within its
//! heights standing out of its walls (two parts flush with them) or less a
//! disc with a round hole within its heights across its end (a pocket open
//! sideways round a post), and its import, each cut by the chained stage's
//! partner, the same volume (`FLUSH`). S9e.4b.4b.2b.2b.3a: by the input's
//! 20th byte, the kernel's own box fused with a stadium standing on its top
//! across its corner (two leaves touching along their caps) or a square with
//! a round hole under a stadium with a round hole hiding its top, and its
//! import, each cut by the chained stage's partner, the same volume
//! (`STACKS`). S9e.4b.4c.2b.1: by the input's 21st byte, the kernel's own
//! stadium plate less a square with a square hole (a post standing in a
//! pocket) or with a U-shaped hole (a U island, its notch a pocket of its own)
//! from inside it past its top, and its import, each cut by the chained
//! stage's partner, the same volume (`NESTED_LEAVES`). Of `NESTED`,
//! `LEAVES`, `RIMS`, `FLUSH`, `STACKS` and `NESTED_LEAVES` one runs an input
//! at most, the first whose key is set, and none where `JOINTS` takes the
//! lens (`own` in `check_boolean`).
use crate::analytic_intersections::Bytes;
use crate::split::{profile, spline_profile};
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

/// Whether two whole tori are decoded (S9d.4b.2b). On since the torus
/// pairs' and turned parts' track (REVIEW_NOTES.md, "The boolean target's
/// last two off switches"): of 1,203 replayed variants (every corpus input
/// moved onto two tori, debug assertions) none fails, 3.6 G instructions
/// at the median, 18.9 G at the ninth decile and 54.6 G at the slowest
/// (3.8, 25 and 82 G before: a torus meeting's root sampling, its jets
/// along both faces' pcurves, the moved result's second validation), and
/// under the sanitizer the slowest 22 s a run on the Mac (30 s before, side
/// by side), the 40 heaviest 8 to 22 s of CPU.
const TORUS_PAIRS: bool = true;

/// Whether a cone object against a cone tool is decoded (S9d.3c). On since
/// the off switches' track (REVIEW_NOTES.md, "The boolean target's off
/// switches"): two coaxial cones' ring took thousands of pieces a turn in
/// its certified integrals (159 G instructions at the slowest of the
/// variants with debug assertions), and a cone section of a given result
/// was rounded on the other cone's frame (`vertex_off_curve`, 29 inputs);
/// of 1,203 replayed variants (every corpus input moved onto a cone object
/// against the cone tool, debug assertions) none fails, 3.3 G at the
/// median, 12.5 G at the ninth decile and 48 G at the slowest, and under the
/// sanitizer the slowest 28 s a run on the Mac, the next 21 to 23 s, where
/// the corpus's slowest input takes 19 s.
const CONE_PAIRS: bool = true;

/// Whether a torus v-segment or wedge tool meets a sphere or cone object,
/// and a cap or zone object a torus or a part, in the tilted or a turned
/// frame too (S9d.4c): their rims' and circles' crossings lie in fields
/// over the rounded frames' inverses. On since the torus pairs' and turned
/// parts' track (REVIEW_NOTES.md, "The boolean target's last two off
/// switches"): of 752 replayed variants (every corpus input moved onto a
/// part or a cap against a torus in those frames, debug assertions) none
/// fails, 10.3 G instructions at the median, 23.3 G at the ninth decile and
/// 43.8 G at the slowest (16, 36 and 80 G before: a circle's resultants
/// with a torus found again for every face, seam and operation, the
/// fields' signs by Sturm-Tarski), and under the sanitizer the slowest 19
/// to 22 s a run on the Mac (37 to 53 s before, side by side), the 40
/// heaviest 13 to 22 s of CPU.
const TURNED_PARTS: bool = true;

/// Whether an operation's first result with curved faces is given to the
/// chained cut and common with the turned box (S9e.1: its model the first
/// arrangement run again, cached). On: the stage keeps its limit of 12
/// faces.
const GIVEN_CURVED: bool = true;

/// Whether the chained stage's partner may be a turned cylinder instead of
/// the turned box (S9e.2: a stack or an S9b.1 result of line prisms given
/// with arcs, decided on its construction's curved arrangement matched to
/// its stored topology), by the chained byte's upper bits.
const GIVEN_ROUND: bool = true;

/// Whether the chained stage's partner may be a sphere of radius 1.25 about
/// the turned box's centre instead (S9e.3a: a given result against a
/// sphere), by the chained byte's next bit.
const GIVEN_BALL: bool = true;

/// Whether the chained stage's partner may reach the first result's
/// meetings of two curved faces (`Meet`, `Rise`, `Toric`) or torus sections
/// (S9e.3b: a given result's meeting met by the partner), and by the
/// chained byte's next bit stand about the middle of the first such edge.
/// On since the arithmetic of fields of degree eight was sped up
/// (REVIEW_NOTES.md's track of that name): the corpus's slowest chained
/// operations reaching them took 60 to 71 s an input under
/// AddressSanitizer on the Mac before, about a third of that now, below the
/// corpus's slowest input (`e36969f1`, a torus meeting). Off, a partner
/// whose bounds meet such an edge's sampled box is not given that result
/// (S9e.3a's results whose meetings the partner does not reach still are).
const GIVEN_MET: bool = true;

/// Whether the object written by the kernel's `.brep` writer and read back
/// is imported and given the chosen operation again (S9e.4a: an imported
/// solid decided on the construction its stored surfaces give), and the
/// chained stage's first result of plane faces cut by the turned box again
/// (S9e.4b.2: an imported polyhedron on its stored vertices), or of one
/// sphere, cylinder or cone face and planes (S9e.4b.3a: a plane piece). On.
const IMPORTED: bool = true;

/// Whether a lens prism (two arcs of circles about `(+-3k, 0)` of radius
/// `5k` meeting at `(0, +-4k)`) in the object's frame is given the chosen
/// operation with the tool, then written, read back, imported and given it
/// again, by the chained byte's next bit (S9e.4b.4a: in the tilted frame its
/// joints round off both circles once read back, each arc taken through its
/// ends). On.
const JOINTS: bool = true;

/// Whether the chained cut's first solid of plane faces (no prism: the
/// turned box's hole in it), written, read back and imported (S9e.4b.4c.1:
/// a polyhedron on its stored vertices in the curved engine), and the
/// kernel's own cut are each cut by a ball about the box's axis, by the
/// chained byte's next bit. On.
const MESHES: bool = true;

/// Whether the kernel's own cylinder on the object's frame less a prism of
/// a square with a square hole (a post in a pocket) or a U-shaped hole (a U
/// island, its notch a pocket of its own), from inside it past its top,
/// written, read back and imported (S9e.4b.4c.2a: a tree whose pockets
/// nest), and the kernel's body are each cut by the chained stage's
/// partner, by the chained byte's next bit (the U by `MESHES`'s). On.
const NESTED: bool = true;

/// Whether the kernel's own box on the object's frame fused with a stadium
/// through both its caps (two prism leaves) or, by the byte's next bit, a
/// stadium plate less a turned square from inside it past its top (a lone
/// leaf with a pocket), written, read back and imported (S9e.4b.4b.2b.2a),
/// and the kernel's body are each cut by the chained stage's partner, by the
/// input's 17th byte (at or above 128: earlier inputs decode as before). On.
const LEAVES: bool = true;

/// Whether the kernel's own stadium plate on the object's frame less a square
/// with a round hole across its rim from inside it past its top (a leaf cut
/// at its rim, the hole's post a tooth to its top) or, by the byte's bit 6,
/// fused with a square with a square hole over its arc end (its hole a well
/// on the plate's top), written, read back and imported (S9e.4b.4b.2b.2b.1),
/// and the kernel's body are each cut by the chained stage's partner, by the
/// input's 18th byte (at or above 128: earlier inputs decode as before). On.
const RIMS: bool = true;

/// Whether the kernel's own stadium plate on the object's frame fused with a
/// diamond within its heights standing out of its walls (two parts flush with
/// them) or, by the byte's bit 6, less a disc with a round hole within its
/// heights across its end (a pocket open sideways round a post), written,
/// read back and imported (S9e.4b.4b.2b.2b.2), and the kernel's body are each
/// cut by the chained stage's partner, by the input's 19th byte (at or above
/// 128: earlier inputs decode as before). On.
const FLUSH: bool = true;

/// Whether the kernel's own box on the object's frame fused with a stadium
/// standing on its top across its corner (two leaves touching along their
/// caps) or, by the byte's bit 6, a square with a round hole under a stadium
/// with a round hole hiding its top (a prism of one cap under a leaf),
/// written, read back and imported (S9e.4b.4b.2b.2b.3a), and the kernel's
/// body are each cut by the chained stage's partner, by the input's 20th
/// byte (at or above 128: earlier inputs decode as before). On.
const STACKS: bool = true;

/// Whether the kernel's own stadium plate on the object's frame less a
/// square with a square hole from inside it past its top (a post standing in
/// a pocket) or, by the byte's bit 6, with a U-shaped hole (a U island, its
/// notch a pocket of its own), written, read back and imported
/// (S9e.4b.4c.2b.1), and the kernel's body are each cut by the chained
/// stage's partner, by the input's 21st byte (at or above 128: earlier inputs
/// decode as before), in the axis-aligned frames alone (on the tilted one the
/// import refuses such a body as before the step). On.
const NESTED_LEAVES: bool = true;

/// Whether a spline prism meets a line prism in frames with different
/// axes (S9f.1's spline walls in the curved engine: creases, generatrices
/// and vertices in the arcs' fields `Q(alpha)`). On: 2,852 spline variants
/// of the corpus replayed with debug assertions in 0.29 s at the median and
/// 2.95 s at the slowest (a lens hole as a tilted tool), that one 23 s
/// under AddressSanitizer where the corpus's slowest input takes 29 s.
const SPLINE_WALLS: bool = true;

/// Whether a spline prism meets a prism with arcs, circles or splines on an
/// exactly parallel axis (S9f.2a's walls: generatrices over the profiles'
/// exact crossings, in the spline arcs' fields): the tool turned about the
/// axis, or in the tilted frame offset by an amount that rounds (the byte
/// after R4's, at or above 128).
const SPLINE_PARALLEL: bool = true;

/// Whether a spline prism meets a prism with arcs or circles (the split
/// target's stadium and round hole) on crossing axes: the tool leaning,
/// tilted or on its side (S9f.2b's meetings of spline walls with
/// cylinders, in the curved engine; S9f.2b.2's loops round a turning point
/// inside the faces, graphs over the height, and cap circles on a wall in a
/// plane along its axis among them: of 1,431 replayed variants the 150
/// reaching them took 0.85 s at the median and 6.8 s at the slowest with
/// debug assertions on a host at load 30 to 48, half the corpus's slowest
/// input there; 32 s under AddressSanitizer at load 6 to 15, the corpus's
/// slowest 36 s).
const SPLINE_CROSSING: bool = true;

/// Whether a spline prism meets a sphere, a cap or a zone (S9f.3a's
/// meetings of spline walls with the sphere, S9f.2b's quadratic over its
/// three rows: branches over the run, loops' graphs over the height, rims
/// and the split's great circle on the walls along their creases or in
/// tower fields). On since the certified integrals beside the loops'
/// turning points were sped up (REVIEW_NOTES.md's track of that name): of
/// 477 replayed variants (debug assertions) none failed, 4.9 G
/// instructions at the median, 11.1 G at the ninth decile and 36.8 G at the
/// slowest (5.5, 14.2 and 472 G before), the corpus's own `d7599dbe` (a
/// lens hole under a leaning sphere whose loops turn back near the lens's
/// tips) the slowest, 472 s under AddressSanitizer before and 24 to 38 s
/// now on a host at load 8 to 22, where the corpus's slowest input,
/// unrelated to the switch, took 87 to 138 s (`e36969f1`, 45 s in the last
/// campaign at load 3 to 7).
const SPLINE_SPHERE: bool = true;

/// Whether a spline prism meets a cone or a frustum (S9f.3b's meetings of
/// spline walls with the cone, S9f.2b's quadratic with its radius row:
/// both nappes for a ruling within its half angle, loops and branches
/// beyond it, its rims on the walls along their creases or in the rims'
/// half-angle charts). On since the certified integrals beside the loops'
/// turning points were sped up (as `SPLINE_SPHERE`): of 402 replayed
/// variants (debug assertions) none failed, 3.3 G instructions at the
/// median, 8.6 G at the ninth decile and 27.1 G at the slowest (3.5, 11.7
/// and 87.7 G before); under AddressSanitizer the slowest (a cone on its
/// side across a lens hole) 23 to 33 s on a host at load 8 to 22, 97 s
/// before.
const SPLINE_CONE: bool = true;

/// R4's knot on the target's sizes: a rectangle `2s` by `t` under a
/// quadratic from `(2s, t)` to `(0, 2t)` whose interior knot of
/// multiplicity two (the degree) at `(s, 2t)` is C1 exactly (its pole the
/// midpoint of `(7s/4, 3t/2)` and `(s/4, 5t/2)`, equal spans), its lifted
/// poles off C1 by rounding in turned frames unless removed first.
fn knot_profile(s: f64, t: f64) -> Option<Profile> {
    let tol = Tolerance::default();
    let poles = [
        (2.0 * s, t),
        (1.75 * s, 1.5 * t),
        (s, 2.0 * t),
        (0.25 * s, 2.5 * t),
        (0.0, 2.0 * t),
    ]
    .iter()
    .map(|(x, y)| Point2::new(*x, *y))
    .collect();
    let curve = BSplineCurve2::new(2, poles, None, vec![0.0, 1.0, 2.0], vec![3, 2, 3]).ok()?;
    let outer = Boundary::path(
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0 * s, 0.0),
            Point2::new(2.0 * s, t),
            Point2::new(0.0, 2.0 * t),
        ],
        vec![
            Segment::Line,
            Segment::Line,
            Segment::Spline(SplineSpan::whole(curve)),
            Segment::Line,
        ],
        tol,
    )
    .ok()?;
    Profile::new(outer, vec![], tol).ok()
}

pub fn check_boolean(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let (ka, kb, flags) = (b.next(), b.next(), b.next());
    let (s1, t1) = (
        1.0 + f64::from(b.next() % 16) / 4.0,
        0.5 + f64::from(b.next() % 16) / 8.0,
    );
    let (s2, t2) = (
        1.0 + f64::from(b.next() % 16) / 4.0,
        0.5 + f64::from(b.next() % 16) / 8.0,
    );
    let (dx, dy) = (
        f64::from(b.next() % 33) / 8.0 - 2.0,
        f64::from(b.next() % 33) / 8.0 - 2.0,
    );
    let h = 0.5 + f64::from(b.next() % 8) / 4.0;
    let (Some(mut pa), Some(mut pb)) = (profile(ka, s1, t1), profile(kb, s2, t2)) else {
        return;
    };
    let tolerance = Tolerance::default();
    let tilted = flags & 1 == 1;
    let (fa, fb) = if tilted {
        let Ok(f) = Frame3::new(
            Point3::new(1.0, -2.0, 0.5),
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            tolerance,
        ) else {
            return;
        };
        (f, f)
    } else {
        let Ok(f) = Frame3::new(
            Point3::new(dx, dy, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tolerance,
        ) else {
            return;
        };
        (Frame3::xy(), f)
    };
    // S9a.2 added heights inside the object's (pockets, cavities) and on
    // its top (touching stacks), chosen by a byte after the others.
    let heights = b.next();
    let (lo, hi) = match (heights % 3, (flags >> 1) % 4) {
        (1, _) => (h / 4.0, h * 0.75),
        (2, _) => (h, h + 1.0),
        (_, 0) => (0.0, h),
        (_, 1) => (-1.0, h + 1.0),
        (_, 2) => (h / 2.0, h * 1.5),
        _ => (h + 1.0, h + 2.0),
    };
    // S9b: the tool's frame turned about the axis, leaning or tilted
    // (frames with different axes), chosen by a byte after the others.
    // S9c.2a: or stood on its side (its axis along x, an exact frame),
    // by the byte's top bit.
    let pick = b.next();
    let offset_fb = fb;
    let fb = match (tilted, pick % 4) {
        (false, 0) if pick >= 128 => {
            let Ok(f) = Frame3::new(
                Point3::new(dx - h / 2.0, dy, h / 2.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                tolerance,
            ) else {
                return;
            };
            f
        }
        (false, turn @ 1..=3) => {
            let (normal, x) = match turn {
                1 => (Vec3::new(0.0, 0.0, 1.0), Vec3::new(3.0, 4.0, 0.0)),
                2 => (Vec3::new(3.0, 0.0, 4.0), Vec3::new(0.0, 1.0, 0.0)),
                _ => (Vec3::new(0.0, 3.0, 4.0), Vec3::new(1.0, 0.0, 0.0)),
            };
            let Ok(f) = Frame3::new(Point3::new(dx, dy, h / 2.0), normal, x, tolerance) else {
                return;
            };
            f
        }
        _ => fb,
    };
    // S9a.2: the split target's spline profiles for the object, the tool
    // or both, chosen by a byte after the others.
    let spline_byte = b.next();
    let splines = spline_byte % 4;
    if splines & 1 == 1 {
        // S9f.1: R4's knot by the byte after the chained one's place (read
        // there, so earlier inputs decode as before).
        let knot = data.get(14).is_some_and(|k| *k >= 128);
        let Some(p) = (if knot {
            knot_profile(s1, t1)
        } else {
            spline_profile(ka, s1, t1)
        }) else {
            return;
        };
        pa = p;
    }
    // S9f.1: a spline prism against a line prism in frames with different
    // axes (off: refused as before its kernel).
    if !SPLINE_WALLS && splines != 0 && !tilted && (pick % 4 != 0 || pick >= 128) {
        return;
    }
    if splines & 2 == 2 {
        let Some(p) = spline_profile(kb, s2, t2) else {
            return;
        };
        pb = p;
    }
    // S9f.2a: a spline prism against a prism with arcs (the split target's
    // stadium and round hole) or splines on an exactly parallel axis: the
    // tool turned about it, or in the tilted frame offset by an amount that
    // rounds (S9a.2's one frame otherwise).
    let curved = |spline: bool, kind: u8| spline || kind % 6 == 2 || kind % 6 == 4;
    let curved_pair = splines != 0
        && spline_byte < 144
        && curved(splines & 1 == 1, ka)
        && curved(splines & 2 == 2, kb);
    let offset_tilt = tilted && splines != 0 && data.get(15).is_some_and(|k| *k >= 128);
    if !SPLINE_PARALLEL && curved_pair && (offset_tilt || (!tilted && pick % 4 == 1)) {
        return;
    }
    // S9f.2b: a spline prism against a prism with arcs on crossing axes
    // (the tool leaning, tilted or on its side; off: refused as before
    // its kernel). Spline walls against spline walls there stay refused
    // (by design, S9f).
    let crossing = !tilted && (pick % 4 >= 2 || (pick % 4 == 0 && pick >= 128));
    if !SPLINE_CROSSING && curved_pair && splines != 3 && crossing {
        return;
    }
    // S9f.3a: a spline prism object against the sphere, cap or zone tool
    // (off: refused as before its kernel).
    if !SPLINE_SPHERE && splines & 1 == 1 && spline_byte >= 192 && matches!((flags >> 5) % 4, 0 | 3)
    {
        return;
    }
    // S9f.3b: a spline prism object against the cone or frustum tool (off:
    // refused as before its kernel).
    let prism_object = match (flags >> 5) % 4 {
        0 | 3 => true,
        2 => !CONE_PAIRS,
        _ => false,
    };
    if !SPLINE_CONE && splines & 1 == 1 && (160..192).contains(&spline_byte) && prism_object {
        return;
    }
    let fb = if offset_tilt {
        // The tilted frame's own normal and `x` (normalized again they
        // could turn by an ulp: no longer exactly parallel).
        let Ok(f) = Frame3::new(
            Point3::new(1.0 + dx / 3.0, -2.0 + dy / 3.0, 0.5 + h / 5.0),
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            tolerance,
        ) else {
            return;
        };
        f
    } else {
        fb
    };
    // S9d.4b.2a: against a whole torus tool (the spline byte in 144..148),
    // the flags' bits 5 and 6 make the object a sphere (1) or a
    // cone (2) in its frame instead of its prism (0 and 3 keep it);
    // S9d.4b.2b: with both clear, the flags' and the heights byte's top
    // bits a whole torus, its radii by the object's scale and its kind
    // byte (no corpus input of 1,435 decodes to one). Two tori's Booleans
    // are the target's slowest under the sanitizer (their exact meetings
    // and their faces' certified integrals along them; in a turned frame,
    // of degree eight, longer): the tool then keeps its offset frame, and
    // one operation (the chained one's byte) is checked by its volume's
    // bounds, its history and its rigid motion.
    let whole_torus = (144..148).contains(&spline_byte);
    // S9d.4c: a torus v-segment or wedge tool (148..160) takes the same
    // curved objects as a whole torus, and a sphere object against either
    // is a cap or zone by the heights byte's top bit, in the axis-aligned
    // frames (the tool offset or on its side; in the tilted or a turned
    // frame only with `TURNED_PARTS`).
    let torus_tool = (144..160).contains(&spline_byte);
    let exact = !tilted && pick % 4 == 0;
    let part_curved = (148..160).contains(&spline_byte) && (exact || TURNED_PARTS);
    let caps = torus_tool && heights >= 128 && (exact || TURNED_PARTS);
    let tori =
        TORUS_PAIRS && whole_torus && (flags >> 5) % 4 == 0 && flags >= 128 && heights >= 128;
    let fb = if tori && !tilted { offset_fb } else { fb };
    let a = match (whole_torus, (flags >> 5) % 4) {
        (true, 0) if tori => {
            let big = 0.75 * s1;
            let small = big * [0.25, 0.375, 0.5, 0.625][usize::from(ka % 4)];
            let turn = std::f64::consts::TAU;
            let Ok((a, _)) =
                Solid::torus_with(OperationId(1), fa, big, small, 0.0, turn, turn, tolerance)
            else {
                return;
            };
            a
        }
        // S9d.3c: against a cone or a sphere tool (the spline byte from
        // 160) the same bits make the object a sphere or a cone: a turned
        // cone's loops with a sphere, a turned cap against a cone, two
        // cones' loops (`CONE_PAIRS`).
        (_, 1) if whole_torus || part_curved || spline_byte >= 160 => {
            // S9d.4c: against a torus or a part, the heights byte's top bit
            // makes it a cap or zone (its next two which), its circles of
            // surd radii.
            let half = std::f64::consts::FRAC_PI_2;
            let (low, high) = if caps {
                [(-half, 0.0), (-0.5, 0.75), (0.25, half), (-half, 0.5)]
                    [usize::from((heights >> 5) % 4)]
            } else {
                (-half, half)
            };
            let Ok((a, _)) =
                Solid::sphere_with(OperationId(1), fa, 0.75 * s1, low, high, tolerance)
            else {
                return;
            };
            a
        }
        (_, 2)
            if whole_torus
                || part_curved
                || spline_byte >= 192
                || (CONE_PAIRS && spline_byte >= 160) =>
        {
            let r = 0.75 * s1;
            let Ok((a, _)) = Solid::cone_with(OperationId(1), fa, r, r / 2.0, h, tolerance) else {
                return;
            };
            a
        }
        _ => {
            let Ok((a, _)) = Solid::extrude_with(OperationId(1), pa, fa, 0.0, h) else {
                return;
            };
            a
        }
    };
    // S9d.1: a sphere, a cap or a zone for the tool in its frame, by the
    // spline byte's two top bits.
    let tool = if spline_byte >= 192 {
        let half = std::f64::consts::FRAC_PI_2;
        let (low, high) = [(-half, half), (-half, 0.0), (-0.5, 0.75), (0.25, half)]
            [usize::from((spline_byte >> 2) % 4)];
        let Ok((tool, _)) = Solid::sphere_with(OperationId(2), fb, 0.75 * s2, low, high, tolerance)
        else {
            return;
        };
        tool
    } else if (144..160).contains(&spline_byte) {
        // S9d.4a: a whole torus about the tool's frame, its radii by the
        // byte's low bits; S9d.4b.1: a v-segment or a wedge by its next two
        // (0 whole, as before), which of them by the flags' fourth and fifth.
        let big = 0.75 * s2;
        let small = big * [0.25, 0.375, 0.5, 0.625][usize::from(spline_byte % 4)];
        let turn = std::f64::consts::TAU;
        let (half, pi) = (std::f64::consts::FRAC_PI_2, std::f64::consts::PI);
        let pick = usize::from((flags >> 3) % 4);
        let (low, high, angle) = match (spline_byte >> 2) % 4 {
            0 => (0.0, turn, turn),
            // The outer and inner halves (their end planes tangent to the
            // torus along their rings, the inner inside out).
            1 => {
                let (low, high) = [(-half, half), (half, 3.0 * half)][pick % 2];
                (low, high, turn)
            }
            2 => (0.0, turn, [half, pi, 3.0 * half, 2.0][pick]),
            // Bands from the axis to the tube's outer side with its cap,
            // and from its inner side below.
            _ => {
                let (low, high) = [(0.5, 2.25), (-2.5, -0.25)][pick % 2];
                (low, high, turn)
            }
        };
        let Ok((tool, _)) =
            Solid::torus_with(OperationId(2), fb, big, small, low, high, angle, tolerance)
        else {
            return;
        };
        tool
    } else if spline_byte >= 160 {
        // S9d.3a: a cone or frustum over the tool's heights, its radii by
        // the byte's next bits.
        let r = 0.75 * s2;
        let (bottom, top) =
            [(r, 0.0), (0.0, r), (r, r / 2.0), (r / 2.0, r)][usize::from((spline_byte >> 2) % 4)];
        let Ok(base) = Frame3::new(
            fb.point(rusty_occt::Point2::new(0.0, 0.0), lo),
            fb.normal(),
            fb.x(),
            tolerance,
        ) else {
            return;
        };
        let Ok((tool, _)) = Solid::cone_with(OperationId(2), base, bottom, top, hi - lo, tolerance)
        else {
            return;
        };
        tool
    } else {
        let Ok((tool, _)) = Solid::extrude_with(OperationId(2), pb, fb, lo, hi) else {
            return;
        };
        tool
    };
    let run = |r: Result<(Vec<Solid>, rusty_occt::history::History), Error>| -> Option<Vec<Solid>> {
        match r {
            Ok((out, _)) => Some(out),
            Err(Error::Degenerate(_) | Error::ComputationLimit(_)) => None,
            // A cavity among several solids (S9a.2's), or arcs in frames
            // with different axes (S9c's).
            // Splines along one curve of different forms, or a spline span
            // along a line (S9a.2's).
            Err(Error::OutOfDomain(m))
                if m.contains("cavity")
                    || m.contains("S9c")
                    || m.contains("S9d")
                    || m.contains("S9e")
                    || m.contains("S9f")
                    || m.contains("different forms")
                    || m.contains("along the plane") =>
            {
                None
            }
            Err(e) => panic!("unexpected error {e}"),
        }
    };
    let volume = |out: &[Solid]| -> f64 { out.iter().map(|s| s.mass_properties().volume).sum() };
    let (va, vb) = (a.mass_properties().volume, tool.mass_properties().volume);
    let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * x.abs().max(y.abs()).max(1.0);
    let motion = rusty_occt::RigidTransform::rotation(
        Point3::new(0.5, -1.0, 2.0),
        Vec3::new(1.0, 2.0, 2.0),
        0.5,
    )
    .expect("a rotation");
    let moves = |out: &[Solid]| {
        for piece in out {
            let (moved, _) = piece
                .transform_with(OperationId(6), motion)
                .expect("a result moves rigidly");
            let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
            assert_eq!(ids(piece), ids(&moved), "a moved result keeps its ids");
            for v in piece.topology().vertices() {
                assert_eq!(
                    piece.classify(v.position).expect("a vertex classifies"),
                    rusty_occt::Location::Boundary,
                    "a result's vertex on its boundary"
                );
            }
        }
    };
    if tori {
        // One operation: fuse between the larger input and their sum,
        // cut between the object less the tool and the object, common at
        // most the smaller (within 1e-9).
        let slack = 1e-9 * (va + vb);
        let op = b.next() % 3;
        let r = match op {
            0 => a.fuse(OperationId(3), &tool),
            1 => a.cut(OperationId(4), &tool),
            _ => a.common(OperationId(5), &tool),
        };
        if let Some(out) = run(r) {
            let v = volume(&out);
            let (lo, hi) = match op {
                0 => (va.max(vb), va + vb),
                1 => (va - vb, va),
                _ => (0.0, va.min(vb)),
            };
            assert!(v >= lo - slack && v <= hi + slack, "{v} outside {lo}..{hi}");
            moves(&out);
        }
        return;
    }
    let fused = run(a.fuse(OperationId(3), &tool));
    let cut = run(a.cut(OperationId(4), &tool));
    let common = run(a.common(OperationId(5), &tool));
    if let (Some(f), Some(m)) = (&fused, &common) {
        assert!(
            near(volume(f), va + vb - volume(m)),
            "fuse {} for {va} + {vb} - {}",
            volume(f),
            volume(m)
        );
    }
    if let (Some(c), Some(m)) = (&cut, &common) {
        assert!(
            near(volume(c), va - volume(m)),
            "cut {} for {va} - {}",
            volume(c),
            volume(m)
        );
    }
    // One operation's first result, chosen by a byte after the others, of
    // at most 12 faces, cut by the box and in common with it: exact
    // fragments of larger stored models took up to 165 s an input under
    // ASan (fuzz/regressions/README.md); the kernel's tests take them.
    // S9e.2: the byte's upper bits may make the partner a cylinder of radius
    // 1.25 on the box's frame.
    let chained = b.next();
    let round = GIVEN_ROUND && (chained / 3) % 2 == 1;
    let ball = GIVEN_BALL && (chained / 6) % 2 == 1;
    let met = GIVEN_MET && (chained / 12) % 2 == 1;
    // One stage of the kernel's own bodies an input (REVIEW_NOTES.md, "One own
    // body an input"): of `NESTED` (the chained byte's bit 96), `LEAVES`,
    // `RIMS`, `FLUSH`, `STACKS` and `NESTED_LEAVES` (the input's 17th to
    // 21st bytes at or above 128, the last in the axis-aligned frames
    // alone) the first whose key is set, and none where the chained byte
    // gives the lens to `JOINTS`. Each is a body of several Booleans written,
    // read back, imported (its chains tried in turn) and cut twice, and an
    // input past 20 bytes ran all whose keys it held, after the lens's
    // Booleans with the tool: CI's scheduled replay at `bbfdee16` timed out
    // on a corpus input running four (`1608293d`, 86 s under
    // AddressSanitizer against the 60 s limit).
    let joints = IMPORTED && JOINTS && (chained / 24) % 2 == 1;
    let keys = [
        (NESTED, ((chained / 96) % 2 == 1).then_some(chained)),
        (LEAVES, data.get(16).copied().filter(|k| *k >= 128)),
        (RIMS, data.get(17).copied().filter(|k| *k >= 128)),
        (FLUSH, data.get(18).copied().filter(|k| *k >= 128)),
        (STACKS, data.get(19).copied().filter(|k| *k >= 128)),
        (
            NESTED_LEAVES,
            data.get(20).copied().filter(|k| *k >= 128 && !tilted),
        ),
    ];
    let own = keys
        .iter()
        .position(|(on, key)| *on && key.is_some())
        .filter(|_| IMPORTED && !joints);
    let own_key = |i: usize| keys[i].1.filter(|_| own == Some(i));
    let chosen = [&fused, &cut, &common][usize::from(chained % 3)];
    let small = |s: &&Solid| {
        s.topology().faces().len() <= 12
            && (GIVEN_CURVED
                || s.topology()
                    .faces()
                    .iter()
                    .all(|f| matches!(f.surface, rusty_occt::topology::Surface::Plane(_))))
    };
    let first = chosen.as_ref().and_then(|out| out.first()).filter(small);
    // S9b.2: an operation's first result is an input again, against a
    // turned box about the object's origin (a fresh operation's ids). S9e.3b:
    // by the chained byte's next bit, about the middle of the first result's
    // first meeting of two curved faces or torus section instead (the box
    // and the cylinder half above it), so the partner crosses that edge.
    let meeting = first.filter(|_| met).and_then(|s| {
        use rusty_occt::topology::Curve3 as C;
        s.topology().edges().iter().find_map(|e| {
            matches!(
                e.curve,
                C::Meet(_) | C::Rise(_) | C::Toric(_) | C::Section(_)
            )
            .then(|| e.curve.point(0.5))
        })
    });
    let origin = match meeting {
        Some(p) if ball => p,
        Some(p) => p + fa.normal() * (-0.5 * h),
        None => fa.point(rusty_occt::Point2::new(0.5, 0.25), h / 3.0),
    };
    let turned = Frame3::new(origin, fa.normal(), fa.x() * 3.0 + fa.y() * 4.0, tolerance).ok();
    let turned = turned.and_then(|f| {
        if ball {
            let half = std::f64::consts::FRAC_PI_2;
            return Solid::sphere_with(OperationId(7), f, 1.25, -half, half, tolerance).ok();
        }
        let outline = if round {
            rusty_occt::Boundary::circle(rusty_occt::Point2::new(0.0, 0.0), 1.25, tolerance)
        } else {
            rusty_occt::Boundary::polygon(
                vec![
                    rusty_occt::Point2::new(-1.0, -1.0),
                    rusty_occt::Point2::new(1.0, -1.0),
                    rusty_occt::Point2::new(1.0, 1.0),
                    rusty_occt::Point2::new(-1.0, 1.0),
                ],
                tolerance,
            )
        }
        .ok()?;
        let profile = rusty_occt::Profile::new(outline, vec![], tolerance).ok()?;
        Solid::extrude_with(OperationId(7), profile, f, 0.0, h).ok()
    });
    // Off (`GIVEN_MET`), a partner reaching a meeting of two curved faces
    // or a torus section of the first result is not given it.
    let reaches = |first: &Solid, partner: &Solid| {
        use rusty_occt::topology::Curve3 as C;
        let b = partner.bounds();
        first.topology().edges().iter().any(|e| {
            if !matches!(
                e.curve,
                C::Meet(_) | C::Rise(_) | C::Toric(_) | C::Section(_)
            ) {
                return false;
            }
            let ps: Vec<[f64; 3]> = (0..=16)
                .map(|k| e.curve.point(f64::from(k) / 16.0).to_array())
                .collect();
            (0..3).all(|i| {
                let lo = ps.iter().map(|p| p[i]).fold(f64::INFINITY, f64::min);
                let hi = ps.iter().map(|p| p[i]).fold(f64::NEG_INFINITY, f64::max);
                let pad = 0.25 * (hi - lo) + 1e-6;
                let (bl, bh) = (b.min.to_array()[i], b.max.to_array()[i]);
                hi + pad >= bl && lo - pad <= bh
            })
        })
    };
    if let Some((box_, _)) = turned {
        if let Some(first) = first.filter(|f| GIVEN_MET || !reaches(f, &box_)) {
            let (c, m) = (
                run(first.cut(OperationId(9), &box_)),
                run(first.common(OperationId(10), &box_)),
            );
            if let (Some(c), Some(m)) = (&c, &m) {
                let v1 = first.mass_properties().volume;
                assert!(near(volume(c), v1 - volume(m)), "chained cut");
            }
            // S9e.4b.2: the first result of plane faces written, read back
            // and imported, cut by the box again: the chained cut's volume.
            // S9e.4b.3a: one of a sphere, cylinder or cone face and plane
            // faces likewise (an imported plane piece).
            use rusty_occt::topology::Surface as S;
            let curved: Vec<&S> = first
                .topology()
                .faces()
                .iter()
                .map(|f| &f.surface)
                .filter(|s| !matches!(s, S::Plane(_)))
                .collect();
            // S9e.4b.4b.1: or of several such faces (a Boolean chain of
            // their primitives; the writer writes circles alone, so their
            // coaxial meetings and plane sections).
            let planar = curved
                .iter()
                .all(|s| matches!(s, S::Sphere { .. } | S::Cylinder { .. } | S::Cone { .. }));
            if let (true, true, Some(c)) = (IMPORTED, planar, &c) {
                if let Some(again) = reimported(first) {
                    if let Some(out) = run(again.cut(OperationId(12), &box_)) {
                        assert!(
                            near(volume(&out), volume(c)),
                            "imported result cut {} for {}",
                            volume(&out),
                            volume(c)
                        );
                    }
                }
            }
            // S9e.4b.4c.1: by the chained byte's next bit, the chained cut's
            // first solid of plane faces (its turned box's hole in it: no
            // prism) written, read back and imported (a polyhedron on its
            // stored vertices, a cavity's too) and the kernel's own cut each
            // less a ball about the box's axis crossing the hole's walls:
            // where both evaluate, the same volume (`MESHES`).
            let cut = c.as_ref().and_then(|c| c.first()).filter(|s| {
                s.topology().faces().len() <= 16
                    && s.topology()
                        .faces()
                        .iter()
                        .all(|f| matches!(f.surface, S::Plane(_)))
            });
            if let (true, Some(cut)) = (IMPORTED && MESHES && (chained / 48) % 2 == 1, cut) {
                let centre = origin + fa.normal() * (0.5 * h);
                let ball = Frame3::new(centre, fa.normal(), fa.x() * 3.0 + fa.y() * 4.0, tolerance)
                    .ok()
                    .and_then(|f| {
                        let half = std::f64::consts::FRAC_PI_2;
                        Solid::sphere_with(OperationId(14), f, 1.25, -half, half, tolerance).ok()
                    });
                if let (Some((ball, _)), Some(again)) = (ball, reimported(cut)) {
                    if let (Some(x), Some(y)) = (
                        run(cut.cut(OperationId(15), &ball)),
                        run(again.cut(OperationId(16), &ball)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported polyhedron less a ball {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
        // S9e.4b.4c.2a: by the chained byte's next bit, the kernel's own
        // cylinder on the object's frame less a square with a square hole
        // (a post in a pocket) or, by `MESHES`'s bit, a U-shaped hole (a U
        // island, its notch a pocket of its own) from inside it past its
        // top, written, read back and imported, each cut by the chained
        // stage's partner: where both evaluate, the same volume (`NESTED`).
        if own_key(0).is_some() {
            let u = (chained / 48) % 2 == 1;
            if let Some(body) = nested_body(fa, s1, h, u, tolerance) {
                if let Some(again) = reimported(&body) {
                    if let (Some(x), Some(y)) = (
                        run(body.cut(OperationId(23), &box_)),
                        run(again.cut(OperationId(24), &box_)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported nested pockets {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
        // S9e.4b.4b.2b.2a: by the input's 17th byte (read there, so earlier
        // inputs decode as before), the kernel's own two prism leaves (a box
        // with a stadium through it) or, by its next bit, a lone leaf with a
        // pocket (a stadium plate less a turned square from its top),
        // written, read back and imported, each cut by the chained stage's
        // partner: where both evaluate, the same volume (`LEAVES`).
        let leaves = own_key(1);
        if let (true, true, Some(k)) = (IMPORTED, LEAVES, leaves) {
            if let Some(body) = leaves_body(fa, s1, h, k >= 192, tolerance) {
                if let Some(again) = reimported(&body) {
                    if let (Some(x), Some(y)) = (
                        run(body.cut(OperationId(26), &box_)),
                        run(again.cut(OperationId(27), &box_)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported prism leaves {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
        // S9e.4b.4b.2b.2b.1: by the input's 18th byte (read there, so earlier
        // inputs decode as before), the kernel's own leaf cut at its rim by a
        // pocket with a post or, by its bit 6, a plate with a well standing on
        // it (`rims_body`), written, read back and imported, each cut by the
        // chained stage's partner: where both evaluate, the same volume
        // (`RIMS`).
        let rims = own_key(2);
        if let (true, true, Some(k)) = (IMPORTED, RIMS, rims) {
            if let Some(body) = rims_body(fa, s1, h, k >= 192, tolerance) {
                if let Some(again) = reimported(&body) {
                    if let (Some(x), Some(y)) = (
                        run(body.cut(OperationId(34), &box_)),
                        run(again.cut(OperationId(35), &box_)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported rim leaves {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
        // S9e.4b.4b.2b.2b.2: by the input's 19th byte (read there, so earlier
        // inputs decode as before), the kernel's own plate with parts flush
        // with its walls: a diamond within its heights standing out of them or,
        // by its bit 6, a pocket open sideways round a post (`flush_body`),
        // written, read back and imported, each cut by the chained stage's
        // partner: where both evaluate, the same volume (`FLUSH`).
        let flush = own_key(3);
        if let (true, true, Some(k)) = (IMPORTED, FLUSH, flush) {
            if let Some(body) = flush_body(fa, s1, h, k >= 192, tolerance) {
                if let Some(again) = reimported(&body) {
                    if let (Some(x), Some(y)) = (
                        run(body.cut(OperationId(39), &box_)),
                        run(again.cut(OperationId(40), &box_)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported flush parts {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
        // S9e.4b.4b.2b.2b.3a: by the input's 20th byte (read there, so earlier
        // inputs decode as before), the kernel's own stack of two prisms whose
        // caps touch: a box with a stadium standing on its top across its
        // corner or, by its bit 6, a holed square under a holed stadium hiding
        // its top (`stacks_body`), written, read back and imported, each cut
        // by the chained stage's partner: where both evaluate, the same volume
        // (`STACKS`).
        let stacks = own_key(4);
        if let (true, true, Some(k)) = (IMPORTED, STACKS, stacks) {
            if let Some(body) = stacks_body(fa, s1, h, k >= 192, tolerance) {
                if let Some(again) = reimported(&body) {
                    if let (Some(x), Some(y)) = (
                        run(body.cut(OperationId(43), &box_)),
                        run(again.cut(OperationId(44), &box_)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported stacks {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
        // S9e.4b.4c.2b.1: by the input's 21st byte (read there, so earlier
        // inputs decode as before), the kernel's own stadium plate less a
        // square with a square hole (a post in a pocket) or, by its bit 6, a
        // U-shaped hole (a U island, its notch a pocket of its own) from inside
        // it past its top (`nested_leaves_body`), written, read back and
        // imported, each cut by the chained stage's partner: where both
        // evaluate, the same volume (`NESTED_LEAVES`).
        // On the tilted frame the import refuses such a body (its pocket walls
        // within rounding of the plate's cylinders' direction, S9's rule in
        // the nested chains' Booleans), so it is built on the world's axes
        // alone.
        let nested = own_key(5);
        if let (true, true, Some(k)) = (IMPORTED, NESTED_LEAVES, nested) {
            if let Some(body) = nested_leaves_body(fa, s1, h, k >= 192, tolerance) {
                if let Some(again) = reimported(&body) {
                    if let (Some(x), Some(y)) = (
                        run(body.cut(OperationId(47), &box_)),
                        run(again.cut(OperationId(48), &box_)),
                    ) {
                        assert!(
                            near(volume(&x), volume(&y)),
                            "imported nested pockets of a leaf {} for {}",
                            volume(&y),
                            volume(&x)
                        );
                    }
                }
            }
        }
    }
    // S9e.4a: the object written by the kernel's writer, read back and
    // imported (`IMPORTED`), given the chosen operation with the same tool:
    // where both evaluate, the same volume within 1e-9 (its construction
    // read off the written surfaces is the object's within rounding).
    if IMPORTED {
        if let (Some(imported), Some(direct)) = (reimported(&a), chosen) {
            let r = match chained % 3 {
                0 => imported.fuse(OperationId(11), &tool),
                1 => imported.cut(OperationId(11), &tool),
                _ => imported.common(OperationId(11), &tool),
            };
            if let Some(out) = run(r) {
                assert!(
                    near(volume(&out), volume(direct)),
                    "imported {} for {}",
                    volume(&out),
                    volume(direct)
                );
            }
        }
    }
    // S9e.4b.4a: a lens in the object's frame given the chosen operation
    // with the tool, and once written, read back and imported (`JOINTS`):
    // where both evaluate, the same volume within 1e-9.
    if joints {
        if let Some(lens) = lens_prism(fa, s1, h, tolerance) {
            if let Some(again) = reimported(&lens) {
                let op = |s: &Solid| match chained % 3 {
                    0 => s.fuse(OperationId(13), &tool),
                    1 => s.cut(OperationId(13), &tool),
                    _ => s.common(OperationId(13), &tool),
                };
                if let (Some(x), Some(y)) = (run(op(&lens)), run(op(&again))) {
                    assert!(
                        near(volume(&x), volume(&y)),
                        "imported lens {} for {}",
                        volume(&y),
                        volume(&x)
                    );
                }
            }
        }
    }
    for out in [&fused, &cut, &common].into_iter().flatten() {
        moves(out);
    }
}

/// A lens of two arcs of circles about `(+-3k, 0)` of radius `5k` meeting at
/// `(0, +-4k)`, `k = s / 8` (dyadic, every number exact), over `[0, h]` in
/// `frame` (S9e.4b.4a).
fn lens_prism(frame: Frame3, s: f64, h: f64, tolerance: Tolerance) -> Option<Solid> {
    let k = s / 8.0;
    let arc = |cx: f64| Segment::Arc {
        center: Point2::new(cx, 0.0),
        radius: 5.0 * k,
        ccw: true,
    };
    let outline = Boundary::path(
        vec![Point2::new(0.0, -4.0 * k), Point2::new(0.0, 4.0 * k)],
        vec![arc(-3.0 * k), arc(3.0 * k)],
        tolerance,
    )
    .ok()?;
    let profile = Profile::new(outline, vec![], tolerance).ok()?;
    Solid::extrude_with(OperationId(12), profile, frame, 0.0, h)
        .ok()
        .map(|(s, _)| s)
}

/// A cylinder of radius `s` over `[0, h + 1]` on `frame` less a prism on it
/// over `[(h + 1) / 2, h + 2]` of the square of half side `5k` (`k = s / 8`,
/// dyadic) with a square hole of half side `2k` (a post in a pocket) or,
/// where `u`, a U-shaped hole, the square of half side `3k` less its notch
/// `[-k, 3k] x [-k, k]` (a U island, its notch a pocket of its own)
/// (S9e.4b.4c.2a).
fn nested_body(frame: Frame3, s: f64, h: f64, u: bool, tolerance: Tolerance) -> Option<Solid> {
    let k = s / 8.0;
    let square = |a: f64| {
        Boundary::polygon(
            vec![
                Point2::new(-a, -a),
                Point2::new(a, -a),
                Point2::new(a, a),
                Point2::new(-a, a),
            ],
            tolerance,
        )
    };
    let hole = if u {
        let p = |x: f64, y: f64| Point2::new(x * k, y * k);
        Boundary::polygon(
            vec![
                p(-3.0, -3.0),
                p(3.0, -3.0),
                p(3.0, -1.0),
                p(-1.0, -1.0),
                p(-1.0, 1.0),
                p(3.0, 1.0),
                p(3.0, 3.0),
                p(-3.0, 3.0),
            ],
            tolerance,
        )
    } else {
        square(2.0 * k)
    }
    .ok()?;
    let profile = Profile::new(square(5.0 * k).ok()?, vec![hole], tolerance).ok()?;
    let top = h + 1.0;
    let (cylinder, _) =
        Solid::cylinder_with(OperationId(20), frame, s, 0.0, top, tolerance).ok()?;
    let (ring, _) =
        Solid::extrude_with(OperationId(21), profile, frame, top / 2.0, top + 1.0).ok()?;
    let (out, _) = cylinder.cut(OperationId(22), &ring).ok()?;
    let [body] = <[Solid; 1]>::try_from(out).ok()?;
    Some(body)
}

/// The kernel's own two prism leaves on `frame` (`k = s / 8`, dyadic): the
/// box of the square of half side `5k` over `[0, h]` fused with a stadium
/// (half discs of radius `2k` about `(-k, k)` and `(k, k)`) over `[-h / 2,
/// 3h / 2]`, through both its caps; or, where `pocket`, a stadium (half discs
/// of radius `3k` about `(-2k, 0)` and `(2k, 0)`) over `[0, h]` less a
/// square of half side `3k / 2` on the frame turned about its axis (its `x`
/// along `3 x + 4 y`) over `[h / 2, 2h]`, from inside it past its top
/// (S9e.4b.4b.2b.2a).
fn leaves_body(frame: Frame3, s: f64, h: f64, pocket: bool, tolerance: Tolerance) -> Option<Solid> {
    let k = s / 8.0;
    let stadium = |c: f64, y: f64, r: f64| {
        let arc = |x: f64| Segment::Arc {
            center: Point2::new(x, y),
            radius: r,
            ccw: true,
        };
        Boundary::path(
            vec![
                Point2::new(-c, y - r),
                Point2::new(c, y - r),
                Point2::new(c, y + r),
                Point2::new(-c, y + r),
            ],
            vec![Segment::Line, arc(c), Segment::Line, arc(-c)],
            tolerance,
        )
    };
    let square = |a: f64| {
        Boundary::polygon(
            vec![
                Point2::new(-a, -a),
                Point2::new(a, -a),
                Point2::new(a, a),
                Point2::new(-a, a),
            ],
            tolerance,
        )
    };
    let prism = |b: Boundary, f: Frame3, h0: f64, h1: f64, op: u64| {
        let profile = Profile::new(b, vec![], tolerance).ok()?;
        Solid::extrude_with(OperationId(op), profile, f, h0, h1)
            .ok()
            .map(|(s, _)| s)
    };
    let (out, _) = if pocket {
        let plate = prism(stadium(2.0 * k, 0.0, 3.0 * k).ok()?, frame, 0.0, h, 28)?;
        let turned = Frame3::new(
            frame.origin(),
            frame.normal(),
            frame.x() * 3.0 + frame.y() * 4.0,
            tolerance,
        )
        .ok()?;
        let hole = prism(square(1.5 * k).ok()?, turned, h / 2.0, 2.0 * h, 29)?;
        plate.cut(OperationId(30), &hole).ok()?
    } else {
        let block = prism(square(5.0 * k).ok()?, frame, 0.0, h, 28)?;
        let rod = prism(stadium(k, k, 2.0 * k).ok()?, frame, -h / 2.0, 1.5 * h, 29)?;
        block.fuse(OperationId(30), &rod).ok()?
    };
    let [body] = <[Solid; 1]>::try_from(out).ok()?;
    Some(body)
}

/// The kernel's own bodies of S9e.4b.4b.2b.2b.1 on `frame` (`k = s / 8`,
/// dyadic): a stadium plate (half discs of radius `3k` about `(0, 0)` and
/// `(4k, 0)`) over `[0, h]` less the square `[-6k, 2k] x [-6k, 6k]` with a
/// hole of radius `3k / 2` about `(-k, k / 2)` over `[h / 2, 2h]` (the pocket
/// across the plate's rim, the hole's post standing to its top); or, where
/// `well`, on the frame turned about its axis (its `x` along `8 x + 15 y`) a
/// stadium plate (half discs of radius `5k / 2` about `(0, 0)` and
/// `(4k, 0)`) over `[0, h]` fused with the square `[5k / 4, 9k] x [-3k / 2,
/// 3k / 2]` with the hole `[7k / 2, 11k / 2] x [-k / 2, k / 2]` over `[h / 2,
/// 3h / 2]` (the ring's near wall and hole standing on the plate's top, its
/// long walls across the plate's arc end at `x = 6k`: a rational crossing,
/// every vertical edge rational; no edge on an arc's circle or its tangent,
/// no wall on an arc's tangent plane).
fn rims_body(frame: Frame3, s: f64, h: f64, well: bool, tolerance: Tolerance) -> Option<Solid> {
    let k = s / 8.0;
    let stadium = |cx: f64, cy: f64, length: f64, r: f64| {
        let arc = |x: f64| Segment::Arc {
            center: Point2::new(x, cy),
            radius: r,
            ccw: true,
        };
        Boundary::path(
            vec![
                Point2::new(cx, cy - r),
                Point2::new(cx + length, cy - r),
                Point2::new(cx + length, cy + r),
                Point2::new(cx, cy + r),
            ],
            vec![Segment::Line, arc(cx + length), Segment::Line, arc(cx)],
            tolerance,
        )
    };
    let square = |x0: f64, y0: f64, x1: f64, y1: f64| {
        Boundary::polygon(
            vec![
                Point2::new(x0, y0),
                Point2::new(x1, y0),
                Point2::new(x1, y1),
                Point2::new(x0, y1),
            ],
            tolerance,
        )
    };
    // The well on the frame turned about its axis (on the world's axes its
    // walls' exact incidences are edges along the other leaf's surfaces in
    // the import's model, and the kernel's own body's Booleans with a ball
    // meet its arc walls' crossings in fields its sphere's sides do not mix).
    let on = if well {
        Frame3::new(
            frame.origin(),
            frame.normal(),
            frame.x() * 8.0 + frame.y() * 15.0,
            tolerance,
        )
        .ok()?
    } else {
        frame
    };
    let prism = |outer: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64| {
        let profile = Profile::new(outer, holes, tolerance).ok()?;
        Solid::extrude_with(OperationId(op), profile, on, h0, h1)
            .ok()
            .map(|(s, _)| s)
    };
    let (out, _) = if well {
        let plate = prism(
            stadium(0.0, 0.0, 4.0 * k, 2.5 * k).ok()?,
            vec![],
            0.0,
            h,
            31,
        )?;
        let hole = square(3.5 * k, -0.5 * k, 5.5 * k, 0.5 * k).ok()?;
        let outer = square(1.25 * k, -1.5 * k, 9.0 * k, 1.5 * k).ok()?;
        let ring = prism(outer, vec![hole], h / 2.0, 1.5 * h, 32)?;
        plate.fuse(OperationId(33), &ring).ok()?
    } else {
        let plate = prism(
            stadium(0.0, 0.0, 4.0 * k, 3.0 * k).ok()?,
            vec![],
            0.0,
            h,
            31,
        )?;
        let hole = Boundary::circle(Point2::new(-k, 0.5 * k), 1.5 * k, tolerance).ok()?;
        let outer = square(-6.0 * k, -6.0 * k, 2.0 * k, 6.0 * k).ok()?;
        let tool = prism(outer, vec![hole], h / 2.0, 2.0 * h, 32)?;
        plate.cut(OperationId(33), &tool).ok()?
    };
    let [body] = <[Solid; 1]>::try_from(out).ok()?;
    Some(body)
}

/// The kernel's own bodies of S9e.4b.4b.2b.2b.2 on `frame` (`k = s / 8`,
/// dyadic): a stadium plate (half discs of radius `2k` about `(0, 0)` and
/// `(6k, 0)`) over `[0, h]` fused with the diamond of corners `(9k / 2, 0)`,
/// `(13k / 2, -7k / 2)`, `(15k / 2, 0)` and `(13k / 2, 7k / 2)` over `[h / 4,
/// 3h / 4]` (two parts standing out of the plate's line walls and arc end);
/// or, where `void`, less the disc of radius `3k` about `(-3k / 2, 0)` with a
/// hole of radius `k` about `(-k / 2, k / 4)` over `[h / 4, 3h / 4]` (a pocket
/// open sideways across the plate's arc end round a post).
fn flush_body(frame: Frame3, s: f64, h: f64, void: bool, tolerance: Tolerance) -> Option<Solid> {
    let k = s / 8.0;
    let arc = |x: f64| Segment::Arc {
        center: Point2::new(x, 0.0),
        radius: 2.0 * k,
        ccw: true,
    };
    let stadium = Boundary::path(
        vec![
            Point2::new(0.0, -2.0 * k),
            Point2::new(6.0 * k, -2.0 * k),
            Point2::new(6.0 * k, 2.0 * k),
            Point2::new(0.0, 2.0 * k),
        ],
        vec![Segment::Line, arc(6.0 * k), Segment::Line, arc(0.0)],
        tolerance,
    )
    .ok()?;
    let prism = |outer: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64| {
        let profile = Profile::new(outer, holes, tolerance).ok()?;
        Solid::extrude_with(OperationId(op), profile, frame, h0, h1)
            .ok()
            .map(|(s, _)| s)
    };
    let plate = prism(stadium, vec![], 0.0, h, 36)?;
    let (out, _) = if void {
        let disc = Boundary::circle(Point2::new(-1.5 * k, 0.0), 3.0 * k, tolerance).ok()?;
        let hole = Boundary::circle(Point2::new(-0.5 * k, 0.25 * k), k, tolerance).ok()?;
        let tool = prism(disc, vec![hole], h / 4.0, 0.75 * h, 37)?;
        plate.cut(OperationId(38), &tool).ok()?
    } else {
        let points = [(4.5, 0.0), (6.5, -3.5), (7.5, 0.0), (6.5, 3.5)]
            .map(|(x, y)| Point2::new(x * k, y * k))
            .to_vec();
        let part = prism(
            Boundary::polygon(points, tolerance).ok()?,
            vec![],
            h / 4.0,
            0.75 * h,
            37,
        )?;
        plate.fuse(OperationId(38), &part).ok()?
    };
    let [body] = <[Solid; 1]>::try_from(out).ok()?;
    Some(body)
}

/// The kernel's own stacks of S9e.4b.4b.2b.2b.3a on the frame turned about
/// `frame`'s axis (its `x` along `8 x + 15 y`: on the world's axes two walls
/// of the leaves cross on the shared plane at a point where their faces'
/// boxes only touch, which S9's plane rule takes for two faces within the
/// resolution of one plane) (`k = s / 8`, dyadic, every crossing of the two
/// outlines rational, no wall on another's circle's tangent): the box `[-3k, 3k] x [-2k, 2k]` over `[0, h]` fused
/// with the stadium of half discs of radius `5k / 4` about `(k, k)` and `(5k,
/// k)` over `[h, 2h]` (standing on its top across its corner, both caps
/// partly exposed); or, where `cover`, the square `[-2k, 2k]^2` with a hole
/// of radius `k` about the origin over `[0, h]` fused with the stadium of half
/// discs of radius `3k` about `(0, 0)` and `(3k, 0)` with a hole of radius `k
/// / 2` about `(k / 4, k / 4)` over `[h, 2h]` (hiding the square's top, its
/// own hole opening into the square's).
fn stacks_body(frame: Frame3, s: f64, h: f64, cover: bool, tolerance: Tolerance) -> Option<Solid> {
    let k = s / 8.0;
    let stadium = |cx: f64, cy: f64, length: f64, r: f64| {
        let arc = |x: f64| Segment::Arc {
            center: Point2::new(x, cy),
            radius: r,
            ccw: true,
        };
        Boundary::path(
            vec![
                Point2::new(cx, cy - r),
                Point2::new(cx + length, cy - r),
                Point2::new(cx + length, cy + r),
                Point2::new(cx, cy + r),
            ],
            vec![Segment::Line, arc(cx + length), Segment::Line, arc(cx)],
            tolerance,
        )
    };
    let square = |x0: f64, y0: f64, x1: f64, y1: f64| {
        Boundary::polygon(
            vec![
                Point2::new(x0, y0),
                Point2::new(x1, y0),
                Point2::new(x1, y1),
                Point2::new(x0, y1),
            ],
            tolerance,
        )
    };
    let on = Frame3::new(
        frame.origin(),
        frame.normal(),
        frame.x() * 8.0 + frame.y() * 15.0,
        tolerance,
    )
    .ok()?;
    let prism = |outer: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64| {
        let profile = Profile::new(outer, holes, tolerance).ok()?;
        Solid::extrude_with(OperationId(op), profile, on, h0, h1)
            .ok()
            .map(|(s, _)| s)
    };
    let (lower, upper) = if cover {
        let hole = Boundary::circle(Point2::new(0.0, 0.0), k, tolerance).ok()?;
        let lower = prism(
            square(-2.0 * k, -2.0 * k, 2.0 * k, 2.0 * k).ok()?,
            vec![hole],
            0.0,
            h,
            41,
        )?;
        let hole = Boundary::circle(Point2::new(0.25 * k, 0.25 * k), 0.5 * k, tolerance).ok()?;
        let upper = prism(
            stadium(0.0, 0.0, 3.0 * k, 3.0 * k).ok()?,
            vec![hole],
            h,
            2.0 * h,
            42,
        )?;
        (lower, upper)
    } else {
        let lower = prism(
            square(-3.0 * k, -2.0 * k, 3.0 * k, 2.0 * k).ok()?,
            vec![],
            0.0,
            h,
            41,
        )?;
        let upper = prism(
            stadium(k, k, 4.0 * k, 1.25 * k).ok()?,
            vec![],
            h,
            2.0 * h,
            42,
        )?;
        (lower, upper)
    };
    let (out, _) = lower.fuse(OperationId(43), &upper).ok()?;
    let [body] = <[Solid; 1]>::try_from(out).ok()?;
    Some(body)
}

/// The kernel's own stadium plate with a pocket holding pockets of its own on
/// `frame` turned about its axis (its `x` along `8 x + 15 y`; `k = s / 8`,
/// dyadic): the plate of half discs of radius `3k` about `(0, 0)` and `(6k,
/// 0)` over `[0, h]` less, over `[h / 2, h + 1]`, the square `[k, 5k] x [-2k,
/// 2k]` with the hole `[2k, 4k] x [-k, k]` (a post standing in the pocket) or,
/// where `u`, the square `[k / 2, 11k / 2] x [-5k / 2, 5k / 2]` with a U-shaped
/// hole (the square `[3k / 2, 9k / 2] x [-3k / 2, 3k / 2]` less its notch
/// `[5k / 2, 9k / 2] x [-k / 2, k / 2]`: a U island, its notch a pocket of its
/// own; S9e.4b.4c.2b.1).
fn nested_leaves_body(
    frame: Frame3,
    s: f64,
    h: f64,
    u: bool,
    tolerance: Tolerance,
) -> Option<Solid> {
    let k = s / 8.0;
    let polygon = |points: &[(f64, f64)]| {
        Boundary::polygon(
            points
                .iter()
                .map(|&(x, y)| Point2::new(x * k, y * k))
                .collect(),
            tolerance,
        )
    };
    let square =
        |x0: f64, y0: f64, x1: f64, y1: f64| polygon(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]);
    let arc = |x: f64| Segment::Arc {
        center: Point2::new(x, 0.0),
        radius: 3.0 * k,
        ccw: true,
    };
    let plate = Boundary::path(
        vec![
            Point2::new(0.0, -3.0 * k),
            Point2::new(6.0 * k, -3.0 * k),
            Point2::new(6.0 * k, 3.0 * k),
            Point2::new(0.0, 3.0 * k),
        ],
        vec![Segment::Line, arc(6.0 * k), Segment::Line, arc(0.0)],
        tolerance,
    )
    .ok()?;
    let on = Frame3::new(
        frame.origin(),
        frame.normal(),
        frame.x() * 8.0 + frame.y() * 15.0,
        tolerance,
    )
    .ok()?;
    let (outer, hole) = if u {
        (
            square(0.5, -2.5, 5.5, 2.5).ok()?,
            polygon(&[
                (1.5, -1.5),
                (4.5, -1.5),
                (4.5, -0.5),
                (2.5, -0.5),
                (2.5, 0.5),
                (4.5, 0.5),
                (4.5, 1.5),
                (1.5, 1.5),
            ])
            .ok()?,
        )
    } else {
        (
            square(1.0, -2.0, 5.0, 2.0).ok()?,
            square(2.0, -1.0, 4.0, 1.0).ok()?,
        )
    };
    let (plate, _) = Solid::extrude_with(
        OperationId(52),
        Profile::new(plate, vec![], tolerance).ok()?,
        on,
        0.0,
        h,
    )
    .ok()?;
    let (ring, _) = Solid::extrude_with(
        OperationId(53),
        Profile::new(outer, vec![hole], tolerance).ok()?,
        on,
        h / 2.0,
        h + 1.0,
    )
    .ok()?;
    let (out, _) = plate.cut(OperationId(54), &ring).ok()?;
    let [body] = <[Solid; 1]>::try_from(out).ok()?;
    Some(body)
}

/// A solid written by the kernel's `.brep` writer, read back and imported
/// (S9e.4a); none where the writer, the reader or the recognition refuses
/// it (a spline prism, S9f). Arcs off their circles once rounded are taken
/// onto them by the Boolean's model (S9e.4b.1); arcs of different circles
/// meeting at a joint rounded off either are taken through their ends
/// there (S9e.4b.4a).
fn reimported(s: &Solid) -> Option<Solid> {
    use rusty_occt::occt_brep::{import, read, write};
    let text = write(s.topology(), s.resolution().linear()).ok()?;
    let doc = read(&text).ok()?;
    let solid = import(&doc).solids.into_iter().next()?;
    let topology = solid.result.ok()?;
    Solid::imported_with(OperationId(11), topology, solid.tolerance)
        .ok()
        .map(|(s, _)| s)
}
