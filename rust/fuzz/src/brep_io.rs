//! OCCT `.brep` interop (T2 of TOPOLOGY_MODEL.md). A base document (a prism
//! from the `identity` target's structure as the kernel writes it, or a small
//! upstream `data/occ` file) takes structure-aware mutations of its tokens and
//! lines. Reading never panics: malformed text is a typed `BrepError`. Every
//! imported solid validates, and when the writer expresses it, it round-trips
//! to the same synthesized counts, cell counts and bit-identical vertices.
//! An unmutated prism always writes and round-trips.
use crate::identity::{build, spec};
use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::occt_brep::{import, read, write};
use rusty_occt::topology::Topology;
use rusty_occt::Tolerance;

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
                    (Ok(x), 2) => format!("{}", x as i64 + i64::from(u.int_in_range(-2..=2)?)),
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
            // Orientation of a shape reference.
            6 => {
                for word in &mut lines[at] {
                    if let Some(rest) = word.strip_prefix('+') {
                        *word = format!("-{rest}");
                    } else if let Some(rest) = word.strip_prefix('-') {
                        if rest.chars().all(|c| c.is_ascii_digit()) {
                            *word = format!("+{rest}");
                        }
                    }
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
    assert_eq!(back.solids.len(), 1);
    let again = back.solids[0]
        .result
        .as_ref()
        .expect("a written solid imports");
    assert_eq!(again.occt_counts(), t.occt_counts());
    assert_eq!(again.faces().len(), t.faces().len());
    assert_eq!(again.edges().len(), t.edges().len());
    assert_eq!(vertices(again), vertices(t));
    true
}

pub fn check_brep_io(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let base = if u.ratio(1, 3).unwrap_or(false) {
        UPSTREAM[u.choose_index(UPSTREAM.len()).unwrap_or(0)].to_string()
    } else {
        let Ok(Some(s)) = spec(&mut u) else { return };
        let Some(solid) = build(&s, None, 1.0, false) else {
            return;
        };
        let tolerance = solid.profile().tolerance().linear();
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
    let Ok(doc) = read(&text) else { return };
    for solid in import(&doc).solids {
        if let Ok(t) = &solid.result {
            assert!(t.check(solid.tolerance).is_empty());
            let tolerance = solid.tolerance.linear().max(Tolerance::default().linear());
            round_trip(t, tolerance);
        }
    }
}
