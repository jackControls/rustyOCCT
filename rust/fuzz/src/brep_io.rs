//! OCCT `.brep` interop (T2 of TOPOLOGY_MODEL.md). A base document (a prism
//! from the `identity` target's structure as the kernel writes it, or a small
//! upstream `data/occ` file) takes structure-aware mutations of its tokens and
//! lines. Reading never panics: malformed text is a typed `BrepError`. Every
//! imported solid validates, and when the writer expresses it, it round-trips
//! to the same synthesized counts, cell counts and bit-identical vertices.
//! An unmutated prism always writes and round-trips. The same bytes also make
//! a cone (S3 of REVIEW_NOTES.md), written with OCCT's seam and degenerated
//! apex edge: it always round-trips, and its mutated text reads, imports and
//! validates or fails cleanly. So do the valid spline fixtures (S4e): spline
//! edges, pcurves and surfaces with B-spline records. And (S6) the valid
//! sheet, shell, wire and acorn fixtures and the face and wire bodies of a
//! profile: written as free shapes, they round-trip, and every free shape
//! of mutated text validates and round-trips or fails cleanly.
use crate::identity::{build, cone_spec, spec, sphere_spec, torus_spec};
use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::occt_brep::{import, read, write};
use rusty_occt::topology::Topology;
use rusty_occt::Tolerance;

#[path = "../../kernel/tests/support/brep_protocol.rs"]
#[allow(dead_code)]
mod brep_protocol;

/// Valid spline cases of the B-rep fixtures (S4).
const SPLINES: [&str; 6] = [
    "case spline_bulge",
    "case spline_cubic_bulge",
    "case spline_rounded_corner",
    "case spline_stadium_pcurve",
    "case spline_stadium_edge",
    "case spline_face_c1",
];

/// A valid spline fixture written with B-spline records always round-trips;
/// its mutated text reads, imports and validates or fails cleanly.
fn check_spline(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(pick) = u.choose_index(SPLINES.len()) else {
        return;
    };
    let block = include_str!("../../fixtures/brep-cases.txt")
        .split("\nend")
        .find(|b| b.trim().lines().next() == Some(SPLINES[pick]))
        .expect("a fixture block");
    let (_, tol, parts) = brep_protocol::parse(block.trim());
    let tolerance = Tolerance::new(tol, 1e-12).unwrap();
    let t = Topology::from_parts(parts, tolerance).expect("a valid spline fixture");
    assert!(round_trip(&t, tol), "an unmutated spline fixture writes");
    let base = write(&t, tol).unwrap();
    if let Ok(text) = mutate(&mut u, &base) {
        check_text(&text);
    }
}

const UPSTREAM: [&str; 3] = [
    include_str!("../../../data/occ/wedge_ok.brep"),
    include_str!("../../../data/occ/mal_ecrou.brep"),
    include_str!("../../../data/occ/solid.brep"),
];

/// Tokens a mutation may insert: record kinds, flags, references,
/// continuity, special numbers.
const WORDS: [&str; 22] = [
    "Ve", "Ed", "Wi", "Fa", "Sh", "So", "Co", "CS", "*", "+1", "-1", "+2", "0", "1", "2", "3", "7",
    "CN", "6CN", "nan", "-inf", "1e308",
];

fn vertices(t: &Topology) -> Vec<[u64; 3]> {
    let mut v: Vec<[u64; 3]> = t
        .vertices()
        .iter()
        .map(|x| x.position.to_array().map(f64::to_bits))
        .collect();
    v.sort_unstable();
    v
}

fn mutate(u: &mut Unstructured, text: &str) -> Result<String> {
    let mut lines: Vec<Vec<String>> = text
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect();
    for _ in 0..u.int_in_range(1..=6)? {
        if lines.is_empty() {
            break;
        }
        let at = u.choose_index(lines.len())?;
        match u.int_in_range(0..=7)? {
            // A number moves, flips sign, becomes an integer or a word.
            0 | 1 => {
                if lines[at].is_empty() {
                    continue;
                }
                let k = u.choose_index(lines[at].len())?;
                let word = &mut lines[at][k];
                *word = match (word.parse::<f64>(), u.int_in_range(0..=3)?) {
                    (Ok(x), 0) => format!("{:?}", x * (1.0 + 1e-9)),
                    (Ok(x), 1) => format!("{:?}", -x),
                    // `as` saturates a huge or infinite number at the end of
                    // i64's range; the nudge saturates there too.
                    (Ok(x), 2) => format!(
                        "{}",
                        (x as i64).saturating_add(i64::from(u.int_in_range(-2..=2)?))
                    ),
                    _ => u.choose(&WORDS)?.to_string(),
                };
            }
            2 => {
                if !lines[at].is_empty() {
                    let k = u.choose_index(lines[at].len())?;
                    lines[at].remove(k);
                }
            }
            3 => {
                let copy = lines[at].clone();
                lines.insert(at, copy);
            }
            4 => {
                lines.remove(at);
            }
            5 => {
                let other = u.choose_index(lines.len())?;
                lines.swap(at, other);
            }
            // Orientation of a shape reference: flipped, or made internal or
            // external (OCCT's `bug21246.brep` has internal edges and faces).
            6 => {
                let mark = *u.choose(&["", "i", "e"])?;
                for word in &mut lines[at] {
                    let Some(rest) = word.strip_prefix(['+', '-']) else {
                        continue;
                    };
                    if !rest.chars().all(|c| c.is_ascii_digit()) {
                        continue;
                    }
                    *word = match (mark, word.starts_with('+')) {
                        ("", true) => format!("-{rest}"),
                        ("", false) => format!("+{rest}"),
                        (m, _) => format!("{m}{rest}"),
                    };
                }
            }
            _ => lines.truncate(at),
        }
    }
    Ok(lines
        .iter()
        .map(|l| l.join(" "))
        .collect::<Vec<_>>()
        .join("\n"))
}

/// Written again, read and imported, a topology comes back with the same
/// counts and vertices.
fn round_trip(t: &Topology, tolerance: f64) -> bool {
    let Ok(text) = write(t, tolerance) else {
        return false;
    };
    let doc = read(&text).expect("the writer's text reads");
    let back = import(&doc);
    assert!(back.unsupported.is_empty(), "{:?}", back.unsupported);
    let solid = t.class().name() == "solid";
    assert_eq!(
        (back.solids.len(), back.free.len()),
        if solid { (1, 0) } else { (0, 1) }
    );
    let again = if solid {
        back.solids[0].result.as_ref()
    } else {
        back.free[0].result.as_ref()
    }
    .expect("a written body imports");
    assert_eq!(again.class(), t.class());
    assert_eq!(again.occt_counts(), t.occt_counts());
    assert_eq!(again.faces().len(), t.faces().len());
    assert_eq!(again.edges().len(), t.edges().len());
    // A ring edge is closed at a vertex in the file and may come back a ring
    // (S6: a disc's or a circle wire's).
    if t.edges().iter().all(|e| !e.is_ring()) && again.edges().iter().all(|e| !e.is_ring()) {
        assert_eq!(vertices(again), vertices(t));
    }
    true
}

/// Mutated text never panics; what imports validates and round-trips.
fn check_text(text: &str) {
    let Ok(doc) = read(text) else { return };
    let im = import(&doc);
    let bodies = im
        .solids
        .iter()
        .map(|s| (&s.result, s.tolerance))
        .chain(im.free.iter().map(|f| (&f.result, f.tolerance)));
    for (result, resolution) in bodies {
        if let Ok(t) = result {
            assert!(t.check(resolution).is_empty());
            let tolerance = resolution.linear().max(Tolerance::default().linear());
            round_trip(t, tolerance);
        }
    }
}

fn check_cone(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = cone_spec(&mut u) else {
        return;
    };
    let Some((solid, _)) = s.build(1.0) else {
        return;
    };
    let tolerance = solid.resolution().linear();
    assert!(
        round_trip(solid.topology(), tolerance),
        "an unmutated cone writes"
    );
    let base = write(solid.topology(), tolerance).unwrap();
    if let Ok(text) = mutate(&mut u, &base) {
        check_text(&text);
    }
}

fn check_sphere(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = sphere_spec(&mut u) else {
        return;
    };
    let Some((solid, _)) = s.build(1.0) else {
        return;
    };
    let tolerance = solid.resolution().linear();
    assert!(
        round_trip(solid.topology(), tolerance),
        "an unmutated sphere writes"
    );
    let base = write(solid.topology(), tolerance).unwrap();
    if let Ok(text) = mutate(&mut u, &base) {
        check_text(&text);
    }
}

fn check_torus(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = torus_spec(&mut u) else {
        return;
    };
    let Some((solid, _)) = s.build(1.0) else {
        return;
    };
    let tolerance = solid.resolution().linear();
    assert!(
        round_trip(solid.topology(), tolerance),
        "an unmutated torus writes"
    );
    let base = write(solid.topology(), tolerance).unwrap();
    if let Ok(text) = mutate(&mut u, &base) {
        check_text(&text);
    }
}

/// S6: a valid sheet, shell, wire or acorn fixture, or the face or wire
/// body of a profile, written as a free shape, round-trips; its mutated
/// text reads, imports and validates or fails cleanly.
fn check_free(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let (t, tol) = if u.ratio(1, 2).unwrap_or(true) {
        // The S6 models (generate_brep_fixtures.sheet_models); the three
        // mutations among them are invalid and skipped.
        let blocks: Vec<&str> = include_str!("../../fixtures/brep-cases.txt")
            .split("\nend")
            .map(str::trim)
            .filter(|b| {
                ["case sheet_", "case shell_", "case wire_", "case acorn_"]
                    .iter()
                    .any(|p| b.starts_with(p))
            })
            .collect();
        let Ok(pick) = u.choose_index(blocks.len()) else {
            return;
        };
        let (_, tol, parts) = brep_protocol::parse(blocks[pick]);
        let Ok(t) = Topology::from_parts(parts, Tolerance::new(tol, 1e-12).unwrap()) else {
            return;
        };
        (t, tol)
    } else {
        let Ok(Some(s)) = spec(&mut u) else { return };
        let Some(solid) = build(&s, None, 1.0, false) else {
            return;
        };
        let profile = solid.profile().unwrap().clone();
        let (op, frame, tolerance) = (solid.operation(), solid.frame(), profile.tolerance());
        let body = if u.ratio(1, 2).unwrap_or(true) {
            rusty_occt::Body::face_from_profile_with(op, profile, frame)
        } else {
            rusty_occt::Body::wire_from_boundary_with(op, profile.outer().clone(), frame, tolerance)
        };
        let Ok((body, _)) = body else { return };
        (body.topology().clone(), tolerance.linear())
    };
    assert!(round_trip(&t, tol), "an unmutated free shape writes");
    let base = write(&t, tol).unwrap();
    if let Ok(text) = mutate(&mut u, &base) {
        check_text(&text);
    }
}

pub fn check_brep_io(data: &[u8]) {
    if data.first() == Some(&0x53) {
        // 'S': the spline family.
        check_spline(&data[1..]);
        return;
    }
    if data.first() == Some(&0x46) {
        // 'F': free shapes (S6).
        check_free(&data[1..]);
        return;
    }
    check_cone(data);
    check_sphere(data);
    check_torus(data);
    let mut u = Unstructured::new(data);
    let base = if u.ratio(1, 3).unwrap_or(false) {
        UPSTREAM[u.choose_index(UPSTREAM.len()).unwrap_or(0)].to_string()
    } else {
        let Ok(Some(s)) = spec(&mut u) else { return };
        let Some(solid) = build(&s, None, 1.0, false) else {
            return;
        };
        let tolerance = solid.resolution().linear();
        assert!(
            round_trip(solid.topology(), tolerance),
            "an unmutated prism writes"
        );
        write(solid.topology(), tolerance).unwrap()
    };
    let text = if u.ratio(1, 8).unwrap_or(true) {
        base
    } else {
        match mutate(&mut u, &base) {
            Ok(t) => t,
            Err(_) => return,
        }
    };
    check_text(&text);
}
