//! The generic B-rep validator against independently generated expectations.
#[path = "support/brep_protocol.rs"]
mod brep_protocol;
use brep_protocol::parse;
use rusty_occt::Tolerance;
use std::collections::BTreeMap;

#[test]
fn complete_issue_sets_match_the_independent_oracle() {
    let expected: BTreeMap<&str, Vec<&str>> = include_str!("../../fixtures/brep-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let (name, issues) = l.split_once('\t').unwrap();
            (name, issues.split(';').filter(|s| !s.is_empty()).collect())
        })
        .collect();
    let cases = include_str!("../../fixtures/brep-cases.txt");
    let mut checked = 0;
    let mut valid = 0;
    let mut failures = Vec::new();
    for block in cases.split("\nend").filter(|b| !b.trim().is_empty()) {
        let (name, tolerance, parts) = parse(block.trim());
        let tolerance = Tolerance::new(tolerance, 1e-12).unwrap();
        let mut actual: Vec<String> = parts
            .check(tolerance)
            .iter()
            .map(|i| i.to_string())
            .collect();
        actual.sort();
        let mut want: Vec<String> = expected[name.as_str()]
            .iter()
            .map(|s| s.to_string())
            .collect();
        want.sort();
        if actual != want {
            failures.push(format!(
                "{name}\n  expected {want:?}\n  actual   {actual:?}"
            ));
        }
        if want.is_empty() {
            valid += 1;
        }
        checked += 1;
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!((checked, valid), (expected.len(), 78));
}

/// Measured enclosures (M5) of every valid case: never below the reference's
/// certain lower value of the gap (sound), never above the bound the
/// reference declared (tight enough for the checker to accept), and the
/// measured parts validate.
///
/// The reference builds each frame's axes as `Frame3::new` does but exactly,
/// while the kernel stores them rounded, so the two gaps are of surfaces a
/// few ulps apart; and the kernel measures in binary64 intervals, whose
/// outward rounding grows with the coordinates. Both comparisons allow
/// `2^-46` of the case's largest coordinate, pcurve coordinate or radius:
/// vacuous for rounding-level gaps, strict for the cases moved half the
/// resolution.
#[test]
fn measured_enclosures_lie_between_the_reference_gap_and_its_declared_bound() {
    let mut lows: BTreeMap<&str, Vec<(&str, usize, f64)>> = BTreeMap::new();
    for line in include_str!("../../fixtures/brep-enclosure-lows.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let mut w = line.split('\t');
        let (name, key, low) = (w.next().unwrap(), w.next().unwrap(), w.next().unwrap());
        let (kind, index) = key.split_once(' ').unwrap();
        lows.entry(name)
            .or_default()
            .push((kind, index.parse().unwrap(), low.parse().unwrap()));
    }
    let (mut compared, mut strict) = (0, 0);
    for block in include_str!("../../fixtures/brep-cases.txt")
        .split("\nend")
        .filter(|b| !b.trim().is_empty())
    {
        let (name, tolerance, declared) = parse(block.trim());
        let Some(keys) = lows.get(name.as_str()) else {
            continue;
        };
        let tolerance = Tolerance::new(tolerance, 1e-12).unwrap();
        let mut bare = declared.clone();
        bare.vertices.iter_mut().for_each(|v| v.enclosure = None);
        bare.fins.iter_mut().for_each(|f| f.enclosure = None);
        bare.faces.iter_mut().for_each(|f| f.enclosure = None);
        let measured = bare.with_measured_enclosures();
        let mut size: f64 = 0.0;
        for v in &declared.vertices {
            size = size.max(
                v.position
                    .to_array()
                    .iter()
                    .fold(0.0, |m, x| m.max(x.abs())),
            );
        }
        for f in &declared.fins {
            for t in [0.0, 1.0] {
                let p = f.pcurve.point(t);
                size = size.max(p.x.abs()).max(p.y.abs());
            }
        }
        for e in &declared.edges {
            if let rusty_occt::topology::Curve3::CircularArc { radius, frame, .. } = e.curve {
                size = size.max(radius).max(
                    frame
                        .origin()
                        .to_array()
                        .iter()
                        .fold(0.0, |m, x| m.max(x.abs())),
                );
            }
        }
        let slack = size * 2f64.powi(-46);
        assert!(measured.check(tolerance).is_empty(), "{name}");
        for &(kind, i, low) in keys {
            let pick = |p: &rusty_occt::topology::TopologyParts| match kind {
                "v" => p.vertices[i].enclosure,
                "u" => p.fins[i].enclosure,
                _ => p.faces[i].enclosure,
            };
            let m = pick(&measured).expect("measured").bound;
            let d = pick(&declared).expect("declared").bound;
            assert!(
                low - slack <= m && m <= d + slack,
                "{name} {kind} {i}: {low} <= {m} <= {d}"
            );
            if low > 1e3 * slack {
                strict += 1;
            }
            compared += 1;
        }
    }
    assert_eq!(compared, lows.values().map(Vec::len).sum::<usize>());
    assert_eq!(lows.len(), 78);
    // The gaps moved half the resolution compare without the allowance
    // mattering: a vertex, a cap fin, and a side fin with its face.
    assert!(strict >= 4, "{strict}");
}

/// Mass properties integrate spline pcurves on planes (S4d): the box with a
/// cap's line pcurve given as an equivalent quadratic spline keeps the box's
/// volume, area and centroid, enclosed.
#[test]
fn spline_pcurves_on_planes_integrate() {
    let case = |name: &str| {
        let block = include_str!("../../fixtures/brep-cases.txt")
            .split("\nend")
            .find(|b| b.trim().lines().next() == Some(&format!("case {name}")))
            .unwrap();
        let (_, tolerance, parts) = parse(block.trim());
        rusty_occt::topology::Topology::from_parts(parts, Tolerance::new(tolerance, 1e-12).unwrap())
            .expect("a valid case")
    };
    let m = case("spline_pcurve_c1")
        .mass_enclosure()
        .expect("integrated");
    let within = |x: [f64; 2], want: f64| x[0] <= want && want <= x[1] && x[1] - x[0] <= 1e-9;
    assert!(within(m.volume, 6.0), "{:?}", m.volume);
    assert!(within(m.surface_area, 22.0), "{:?}", m.surface_area);
    for (c, want) in m.centroid.iter().zip([1.5, 1.0, 0.5]) {
        assert!(within(*c, want), "{c:?}");
    }
}

/// Every valid spline case's certified mass enclosure (S4d) contains the
/// independent reference's quadrature (`brep-spline-mass.tsv`): volume,
/// area, centroid and the inertia about it.
#[test]
fn spline_mass_encloses_the_reference() {
    let cases = include_str!("../../fixtures/brep-cases.txt");
    let mut checked = 0;
    for line in include_str!("../../fixtures/brep-spline-mass.tsv")
        .lines()
        .skip(1)
    {
        let (name, values) = line.split_once('\t').unwrap();
        let want: Vec<f64> = values.split(' ').map(|x| x.parse().unwrap()).collect();
        let block = cases
            .split("\nend")
            .find(|b| b.trim().lines().next() == Some(&format!("case {name}")))
            .unwrap();
        let (_, tolerance, parts) = parse(block.trim());
        let t = rusty_occt::topology::Topology::from_parts(
            parts,
            Tolerance::new(tolerance, 1e-12).unwrap(),
        )
        .expect("a valid case");
        let m = t
            .mass_enclosure()
            .unwrap_or_else(|| panic!("{name}: integrated"));
        let mut got = vec![m.volume, m.surface_area];
        got.extend(m.centroid);
        got.extend(m.inertia.iter().flatten().copied());
        assert_eq!(got.len(), want.len(), "{name}");
        for (k, ([lo, hi], x)) in got.iter().zip(&want).enumerate() {
            let slack = 1e-12 * x.abs().max(1.0);
            assert!(
                lo - slack <= *x && *x <= hi + slack,
                "{name} value {k}: {x} outside [{lo}, {hi}]"
            );
        }
        checked += 1;
    }
    assert_eq!(checked, 16);
}

fn topology_of(name: &str) -> rusty_occt::topology::Topology {
    let block = include_str!("../../fixtures/brep-cases.txt")
        .split("\nend")
        .find(|b| b.trim().lines().next() == Some(&format!("case {name}")))
        .unwrap();
    let (_, tolerance, parts) = parse(block.trim());
    rusty_occt::topology::Topology::from_parts(parts, Tolerance::new(tolerance, 1e-12).unwrap())
        .expect("a valid case")
}

/// D9 (S6): every valid case's computed class is the independent
/// reference's (`brep-classes.tsv`).
#[test]
fn classes_match_the_independent_reference() {
    let mut checked = 0;
    for line in include_str!("../../fixtures/brep-classes.tsv")
        .lines()
        .skip(1)
    {
        let (name, class) = line.split_once('\t').unwrap();
        assert_eq!(topology_of(name).class().name(), class, "{name}");
        checked += 1;
    }
    assert_eq!(checked, 78);
}

/// S6: each valid sheet's, wire's and acorn's certified area or length and
/// centre contains OCCT's `BRepGProp` value from the pre-implementation
/// capture, up to `1e-9` of the row's largest magnitude (OCCT's own
/// integration error; the circle's centre is `7e-15` off its origin).
#[test]
fn sheet_and_wire_measures_contain_the_native_properties() {
    let mut checked = 0;
    for line in include_str!("../../fixtures/occt-sheet-preimplementation/native.txt").lines() {
        let words: Vec<&str> = line.split(' ').collect();
        if words[1] != "G" {
            continue;
        }
        let name = words[0];
        // The three mutations are invalid.
        if [
            "sheet_square_vertex_moved",
            "sheet_open_box_pcurve_shift",
            "wire_disconnected",
        ]
        .contains(&name)
        {
            continue;
        }
        let t = topology_of(name);
        let want: Vec<f64> = words[2..].iter().map(|x| x.parse().unwrap()).collect();
        let m = t
            .measure_enclosure()
            .unwrap_or_else(|| panic!("{name}: measured"));
        let slack = 1e-9 * want.iter().fold(1.0_f64, |a, x| a.max(x.abs()));
        for (k, ([lo, hi], x)) in std::iter::once(m.measure)
            .chain(m.centre)
            .zip(&want)
            .enumerate()
        {
            assert!(
                lo - slack <= *x && *x <= hi + slack,
                "{name} value {k}: {x} outside [{lo}, {hi}]"
            );
        }
        checked += 1;
    }
    assert_eq!(checked, 15);
}
