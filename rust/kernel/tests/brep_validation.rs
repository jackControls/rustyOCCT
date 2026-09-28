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
    assert_eq!((checked, valid), (expected.len(), 81));
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
    assert_eq!(lows.len(), 81);
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
/// area, centroid and the inertia about it. The reference agrees with
/// itself on halved pieces within `1e-20` of each property's scale
/// (`generate_brep_fixtures.spline_mass`), so the only slack is that and
/// the rounding of its twenty printed digits.
///
/// F8: every enclosure's width is within `1e-12` of its property's scale,
/// the analytic family's (the first-order enclosures of S4d were 2e-3 to 5e-3
/// on nonrational walls' areas and 5e-2 to 0.6 on the rational corner):
/// the volume and area relative to themselves, the centroid to the body's
/// length `L = V^(1/3)` (plus the rounding of its world coordinates), the
/// inertia to `V L^2`.
#[test]
fn spline_mass_encloses_the_reference() {
    let mut checked = 0;
    let mut worst = [0.0_f64; 4];
    for line in include_str!("../../fixtures/brep-spline-mass.tsv")
        .lines()
        .skip(1)
    {
        let (name, values) = line.split_once('\t').unwrap();
        let want: Vec<f64> = values.split(' ').map(|x| x.parse().unwrap()).collect();
        let m = topology_of(name)
            .mass_enclosure()
            .unwrap_or_else(|| panic!("{name}: integrated"));
        let mut got = vec![m.volume, m.surface_area];
        got.extend(m.centroid);
        got.extend(m.inertia.iter().flatten().copied());
        assert_eq!(got.len(), want.len(), "{name}");
        // Volume, area, centroid, inertia: the reference's scale of each.
        let group = |k: usize| match k {
            0 => 0..1,
            1 => 1..2,
            2..=4 => 2..5,
            _ => 5..14,
        };
        let scale = |k: usize| want[group(k)].iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        for (k, ([lo, hi], x)) in got.iter().zip(&want).enumerate() {
            let slack = 1e-20 * scale(k) + 4.0 * f64::EPSILON * x.abs() + 1e-25;
            assert!(
                lo - slack <= *x && *x <= hi + slack,
                "{name} value {k}: {x} outside [{lo}, {hi}]"
            );
        }
        let length = want[0].cbrt();
        let relative = |k: usize| {
            let [lo, hi] = got[k];
            let (width, x) = (hi - lo, want[k].abs());
            match k {
                0 | 1 => width / x,
                2..=4 => (width - 8.0 * f64::EPSILON * x).max(0.0) / length,
                _ => width / (want[0] * length * length),
            }
        };
        for (k, w) in (0..14).map(|k| (k, relative(k))) {
            let kind = [0, 1, 2, 2, 2, 3][k.min(5)];
            worst[kind] = worst[kind].max(w);
            assert!(w <= 1e-12, "{name} value {k}: relative width {w:e}");
        }
        checked += 1;
    }
    assert_eq!(checked, 19);
    eprintln!(
        "relative widths: volume {:.1e}, area {:.1e}, centroid {:.1e}, inertia {:.1e}",
        worst[0], worst[1], worst[2], worst[3]
    );
}

/// F8: a parallel of a cone, sphere or torus written as a degree-1 spline
/// pcurve bounds the same face as the line it replaces (both are linear in
/// the fraction), so the mass enclosures along the spline, by the certified
/// quadrature of `F(u, v) du`, overlap the closed forms' along the line and
/// are as narrow (within `1e-12` of each property's scale). The spline's
/// use is certified by Taylor bounds, which need a coarser enclosure than
/// the line's harmonic one: the case runs at tolerance `1e-4`.
#[test]
fn spline_parallels_integrate_as_their_lines() {
    for name in [
        "cone_frustum",
        "sphere_zone",
        "torus_segment",
        "torus_outer_half",
    ] {
        let block = include_str!("../../fixtures/brep-cases.txt")
            .split("\nend")
            .find(|b| b.trim().lines().next() == Some(&format!("case {name}")))
            .unwrap()
            .trim();
        let mut replaced = 0;
        let text: Vec<String> = block
            .lines()
            .map(|line| {
                let w: Vec<&str> = line.split_whitespace().collect();
                if line.starts_with("tolerance") {
                    return "tolerance 0.0001".to_string();
                }
                // A use `u e S line u0 v0 u1 v1 enc x` along a parallel.
                if w.len() == 10 && w[0] == "u" && w[3] == "line" && w[5] == w[7] {
                    replaced += 1;
                    let poles = w[4..8].join(" ");
                    return format!(
                        "u {} {} bspline 1 0 2 0.0 1.0 2 2 2 {poles} 1.0 1.0 enc 1e-5",
                        w[1], w[2]
                    );
                }
                line.to_string()
            })
            .collect();
        assert!(replaced > 0, "{name}");
        let (_, tolerance, parts) = parse(&text.join("\n"));
        let spline = rusty_occt::topology::Topology::from_parts(
            parts,
            Tolerance::new(tolerance, 1e-12).unwrap(),
        )
        .unwrap_or_else(|issues| panic!("{name}: {issues:?}"))
        .mass_enclosure()
        .expect("integrated");
        let line = topology_of(name).mass_enclosure().expect("integrated");
        let pair = |m: &rusty_occt::topology::MassEnclosure| {
            let mut out = vec![m.volume, m.surface_area];
            out.extend(m.centroid);
            out.extend(m.inertia.iter().flatten().copied());
            out
        };
        let volume = line.volume[1];
        let length = volume.cbrt();
        for (k, ([a, b], [c, d])) in pair(&spline).into_iter().zip(pair(&line)).enumerate() {
            assert!(
                a <= d && c <= b,
                "{name} value {k}: [{a}, {b}] and [{c}, {d}]"
            );
            let scale = match k {
                0 => volume,
                1 => line.surface_area[1],
                2..=4 => length,
                _ => volume * length * length,
            };
            assert!(b - a <= 1e-12 * scale, "{name} value {k}: width {}", b - a);
        }
    }
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
    assert_eq!(checked, 81);
}

/// S6: each valid sheet's, wire's and acorn's certified area or length and
/// centre contains OCCT's `BRepGProp` value from the pre-implementation
/// capture, up to `1e-9` of the row's largest magnitude (OCCT's own
/// integration error; the circle's centre is `7e-15` off its origin). F8:
/// the bulge's spline wall is now enclosed to `1e-15` around its closed
/// form `√2 + asinh 1`, and OCCT's area is `1.8e-8` relative above it (its
/// centre `4.6e-9` off), so that row allows `2e-8`, as S4d allows `1e-8`
/// for the bulge's total area.
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
        let relative = if name == "sheet_spline_wall" {
            2e-8
        } else {
            1e-9
        };
        let slack = relative * want.iter().fold(1.0_f64, |a, x| a.max(x.abs()));
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
