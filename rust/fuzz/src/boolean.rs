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
//! or zone tool (against a cone tool off, `CONE_PAIRS`), S9d.4c a sphere,
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
//! the first such edge (`GIVEN_MET`, off for time). S9f.1: a spline prism
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
//! the cylinders (`SPLINE_CROSSING`), their loops round a cylinder
//! refused (S9f.2b.2's); spline walls against spline walls on crossing
//! axes, a curved solid or a given result stay refused (S9f's). S9e.4a: the
//! object written by the kernel's `.brep` writer and read back is imported
//! (`Solid::imported_with`) and given the chosen operation with the tool
//! again, its volume the object's own result's (`IMPORTED`).
use crate::analytic_intersections::Bytes;
use crate::split::{profile, spline_profile};
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

/// Whether two whole tori are decoded (S9d.4b.2b). Off: since the certified
/// integrals along their meetings were sped up (REVIEW_NOTES.md's track of
/// that name) intersecting variants take about 3 s under the sanitizer at
/// the median, 22 s at the ninth decile and 51 s at the slowest of 44 on
/// the Mac (46, 111 and several minutes before), too close to the target's
/// 60 seconds at the Linux runners' 2.6 times; the kernel's tests and the
/// replayed variants cover them until the degree-eight arrangement's
/// arithmetic is faster too.
const TORUS_PAIRS: bool = false;

/// Whether a cone object against a cone tool is decoded (S9d.3c). Off
/// until two cones' certified integrals are fast enough for the target's 60
/// seconds under the sanitizer (up to 80 s with debug assertions alone, 14
/// of 711 variants above 20 s, the time in the result's mass and
/// validation): the kernel's tests and the replays of 711 cone-pair
/// variants cover them meanwhile. A sphere object against a cone tool and a
/// cone object against a sphere, cap or zone tool stay on (26 s at most).
const CONE_PAIRS: bool = false;

/// Whether a torus v-segment or wedge tool meets a sphere or cone object,
/// and a cap or zone object a torus or a part, in the tilted or a turned
/// frame too (S9d.4c). Off: their rims' and circles' crossings lie in
/// fields over the rounded frames' inverses, and of 3,000 replayed variants
/// (debug assertions, no sanitizer) the 1,201 caps and 592 parts in turned
/// frames took 2.1 s and 0.6 s at the median, 7.5 s and 4.6 s at the ninth
/// decile and 26 s and 50 s at the slowest (about twelve times that under
/// the sanitizer), against 2.3 s at the slowest of 457 in the axis-aligned
/// frames, which stay on; parts against prisms with arcs are on in every
/// frame (4.5 s at the slowest of 750).
const TURNED_PARTS: bool = false;

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
/// Off: the corpus's slowest chained operations reaching them take 60 to
/// 71 s an input under AddressSanitizer on the Mac (8 to 10 s with debug
/// assertions: algebraic vertices in fields of degree eight through the
/// second arrangement), past the target's 60 s and the Linux runners' 2.6
/// times; off, a partner whose bounds meet such an edge's sampled box is
/// not given that result (S9e.3a's results whose meetings the partner does
/// not reach still are). The kernel's tests and the corpus replayed with it
/// on (debug assertions, every input's partner centred on its meeting
/// where it has one) cover them.
const GIVEN_MET: bool = false;

/// Whether the object written by the kernel's `.brep` writer and read back
/// is imported and given the chosen operation again (S9e.4a: an imported
/// solid decided on the construction its stored surfaces give). On.
const IMPORTED: bool = true;

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
/// cylinders, in the curved engine; their loops round the cylinder,
/// S9f.2b.2's, stay refused).
const SPLINE_CROSSING: bool = true;

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
    // (S9f.3).
    let crossing = !tilted && (pick % 4 >= 2 || (pick % 4 == 0 && pick >= 128));
    if !SPLINE_CROSSING && curved_pair && splines != 3 && crossing {
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
    for out in [&fused, &cut, &common].into_iter().flatten() {
        moves(out);
    }
}

/// A solid written by the kernel's `.brep` writer, read back and imported
/// (S9e.4a); none where the writer, the reader or the recognition refuses
/// it (a spline prism, S9f; arcs off their circles once rounded are
/// refused by the Boolean, S9e.4b).
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
