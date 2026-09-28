//! S7b.4: two cones, and a cone whose rational apex lies on a sphere or a
//! cylinder, as traced curves or lines, against the independent reference
//! (`fixtures/ruled-curve-*.txt|tsv` from
//! `tools/generate_ruled_curve_fixtures.py`).
use rusty_occt::intersection::{
    surface_surface, AnalyticItem, SurfaceIntersection, TracedComponent, TracedCurve,
};
use rusty_occt::topology::Surface;
use rusty_occt::{Frame3, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;
use std::f64::consts::PI;

fn surface(words: &[&str]) -> Surface {
    let v: Vec<f64> = words[1..].iter().map(|w| w.parse().unwrap()).collect();
    let frame = Frame3::new(
        Point3::new(v[0], v[1], v[2]),
        Vec3::new(v[3], v[4], v[5]),
        Vec3::new(v[6], v[7], v[8]),
        Tolerance::default(),
    )
    .unwrap();
    match words[0] {
        "sphere" => Surface::Sphere {
            frame,
            radius: v[9],
        },
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        _ => Surface::Cone {
            frame,
            radius: v[9],
            half_angle: v[10],
        },
    }
}

fn cases() -> Vec<(String, Surface, Surface)> {
    include_str!("../../fixtures/ruled-curve-cases.txt")
        .split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|b| {
            let lines: Vec<Vec<&str>> = b
                .trim()
                .lines()
                .map(|l| l.split_whitespace().collect())
                .collect();
            (
                lines[0][1].to_string(),
                surface(&lines[1][1..]),
                surface(&lines[2][1..]),
            )
        })
        .collect()
}

/// Each case's rows as words, in the reference's order.
fn expected() -> BTreeMap<String, Vec<Vec<String>>> {
    let mut out: BTreeMap<String, Vec<Vec<String>>> = BTreeMap::new();
    for line in include_str!("../../fixtures/ruled-curve-expected.tsv")
        .lines()
        .skip(1)
    {
        let (name, row) = line.split_once('\t').unwrap();
        out.entry(name.into())
            .or_default()
            .push(row.split(' ').map(String::from).collect());
    }
    out
}

/// An angle's enclosure contains `x` up to whole turns.
fn contains_angle(e: [f64; 2], x: f64) -> bool {
    [-2.0 * PI, 0.0, 2.0 * PI]
        .iter()
        .any(|s| contains(e, x + s))
}
fn contains(e: [f64; 2], x: f64) -> bool {
    let slack = 1e-25 * x.abs();
    e[0] - slack <= x && x <= e[1] + slack
}
fn numbers(row: &[String]) -> Vec<f64> {
    row[1..].iter().filter_map(|w| w.parse().ok()).collect()
}

fn check(name: &str, c: &TracedCurve, rows: &[Vec<String>]) -> Result<(), String> {
    let want = |kind: &'static str| rows.iter().filter(move |r| r[0] == kind);
    let folds: Vec<_> = want("fold").collect();
    if folds.len() != c.folds().len() {
        return Err(format!("{} folds for {}", c.folds().len(), folds.len()));
    }
    for r in &folds {
        let x = numbers(r);
        let hits = c
            .folds()
            .iter()
            .filter(|f| {
                contains_angle(f.phi, x[0])
                    && contains_angle(f.t, x[1])
                    && (0..3).all(|k| contains(f.point[k], x[2 + k]))
            })
            .count();
        if hits != 1 {
            return Err(format!("{name}: fold {x:?} in {hits}"));
        }
    }
    // The apex, when it lies on the other surface.
    let apexes: Vec<_> = want("apex").collect();
    if apexes.len() != c.nodes().len() {
        return Err(format!("{} apexes for {}", c.nodes().len(), apexes.len()));
    }
    for (r, n) in apexes.iter().zip(c.nodes()) {
        let x = numbers(r);
        let crossing = r.last().unwrap() == "crossing";
        if n.crossing != crossing || !(0..3).all(|k| contains(n.point[k], x[k])) {
            return Err(format!("{name}: apex {x:?}"));
        }
    }
    // Components: folds, winding numbers and crossings of infinity.
    let mut got: Vec<String> = Vec::new();
    for comp in c.components() {
        match comp {
            TracedComponent::Smooth {
                folds,
                winding,
                infinite,
                ..
            } => got.push(format!(
                "component {folds} {} {} {infinite}",
                winding[0], winding[1]
            )),
            TracedComponent::Crossing {
                tracks, infinite, ..
            } => got.push(format!(
                "component 0 1 {} {infinite}",
                winding_of(c, tracks)
            )),
            TracedComponent::Isolated { .. } => {}
        }
    }
    got.sort();
    let mut comps: Vec<String> = want("component").map(|r| r.join(" ")).collect();
    comps.sort();
    if got != comps {
        return Err(format!("{got:?} for {comps:?}"));
    }
    // Rings: each reference point (at u = 0, or u = 1 on a factor's chart)
    // on some track.
    for r in want("ring") {
        let x = numbers(r);
        let hit = (0..c.tracks().len()).any(|t| {
            c.point_at(t, x[0])
                .is_ok_and(|p| (0..3).all(|k| contains(p[k], x[2 + k])))
        });
        if !hit {
            return Err(format!("{name}: ring {x:?}"));
        }
    }
    Ok(())
}

/// The winding in the chart's second parameter of the component through the
/// apex (its tracks cover one turn in u).
fn winding_of(c: &TracedCurve, tracks: &[usize]) -> i64 {
    let mut total = 0.0;
    for &k in tracks {
        let tr = &c.tracks()[k];
        let (a, b) = (c.t_at(k, tr.phi[0]).unwrap(), c.t_at(k, tr.phi[1]).unwrap());
        total += 0.5 * (b[0] + b[1]) - 0.5 * (a[0] + a[1]);
    }
    (total / c.period()).round() as i64
}

/// Lines and points (two cones with one apex), each reference row inside
/// one of the kernel's items of its kind.
fn check_items(name: &str, items: &[AnalyticItem], rows: &[Vec<String>]) -> Result<(), String> {
    if items.len() != rows.len() {
        return Err(format!("{name}: {} items for {}", items.len(), rows.len()));
    }
    let mut left: Vec<&AnalyticItem> = items.iter().collect();
    for r in rows {
        let x = numbers(r);
        let Some(k) = left.iter().position(|item| {
            item.kind() == r[0]
                && item.values().len() == x.len()
                && item.values().iter().zip(&x).all(|(e, v)| contains(*e, *v))
        }) else {
            return Err(format!("{name}: {} {x:?}", r[0]));
        };
        left.remove(k);
    }
    Ok(())
}

/// The kernel stores the frames the reference intersected, bit for bit.
#[test]
fn stored_normals_are_the_reference_inputs() {
    let all = cases();
    let rows: Vec<&str> = include_str!("../../fixtures/ruled-curve-frames.tsv")
        .lines()
        .skip(1)
        .collect();
    assert_eq!(rows.len(), 2 * all.len());
    for ((name, a, b), pair) in all.iter().zip(rows.chunks(2)) {
        for (s, row) in [a, b].into_iter().zip(pair) {
            let frame = match s {
                Surface::Cylinder { frame: f, .. }
                | Surface::Cone { frame: f, .. }
                | Surface::Sphere { frame: f, .. } => f,
                _ => unreachable!(),
            };
            let want: Vec<u64> = row
                .split('\t')
                .nth(2)
                .unwrap()
                .split(' ')
                .map(|h| u64::from_str_radix(h, 16).unwrap())
                .collect();
            let got: Vec<u64> = frame.normal().to_array().map(f64::to_bits).to_vec();
            assert_eq!(got, want, "{name}");
        }
    }
}

#[test]
fn every_case_matches_the_exact_reference() {
    let want = expected();
    let all = cases();
    let mut failures = Vec::new();
    for (name, a, b) in &all {
        let got = surface_surface(a, b).unwrap_or_else(|e| panic!("{name}: {e}"));
        let rows = &want[name];
        let verdict = match &got {
            SurfaceIntersection::Empty if rows[0][0] == "empty" => Ok(()),
            SurfaceIntersection::Traced(c) => check(name, c, rows),
            SurfaceIntersection::Items(items) => check_items(name, items, rows),
            other => Err(format!("{other:?}")),
        };
        if let Err(why) = verdict {
            failures.push(format!("{name}: {why}"));
        }
        assert_eq!(surface_surface(b, a).unwrap(), got, "{name}: order");
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(all.len(), 12);
}

/// Points along every track lie on both surfaces, enclosed tightly.
#[test]
fn sampled_points_lie_on_both_surfaces() {
    let mut count = 0;
    for (name, a, b) in cases() {
        let Ok(SurfaceIntersection::Traced(c)) = surface_surface(&a, &b) else {
            continue;
        };
        for (k, track) in c.tracks().iter().enumerate() {
            let [lo, hi] = track.phi;
            for j in 0..=4 {
                let phi = lo + (hi - lo) * f64::from(j) / 4.0;
                // Points at infinity (a cone's unbounded branches) and far
                // ones are skipped.
                let Ok(e) = c.point_at(k, phi) else {
                    continue;
                };
                if e.iter().any(|[l, h]| l.abs().max(h.abs()) > 1e3) {
                    continue;
                }
                for [l, h] in e {
                    assert!(
                        h - l <= 1e-12 * l.abs().max(h.abs()).max(1.0),
                        "{name}: wide {e:?}"
                    );
                }
                let mid = |k: usize| 0.5 * e[k][0] + 0.5 * e[k][1];
                let p = Vec3::new(mid(0), mid(1), mid(2));
                let scale = p.length().max(1.0);
                for s in [&a, &b] {
                    let gap = distance(s, p);
                    assert!(gap <= 4e-12 * scale, "{name}: {gap} off {s:?} at {phi}");
                }
                count += 1;
            }
        }
    }
    assert!(count > 40, "{count} points");
}

fn distance(s: &Surface, p: Vec3) -> f64 {
    match s {
        Surface::Cylinder { frame, radius } => {
            let rel = p - (frame.origin() - Point3::ORIGIN);
            let h = rel.dot(frame.normal());
            ((rel - frame.normal() * h).length() - radius).abs()
        }
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            let rel = p - (frame.origin() - Point3::ORIGIN);
            let h = rel.dot(frame.normal());
            let rho = (rel - frame.normal() * h).length();
            let (s, c) = half_angle.sin_cos();
            let shift = radius * c + h * s;
            (rho * c - shift).abs().min((rho * c + shift).abs())
        }
        Surface::Sphere { frame, radius } => {
            ((p - (frame.origin() - Point3::ORIGIN)).length() - radius).abs()
        }
        Surface::Torus {
            frame,
            major,
            minor,
        } => {
            let rel = p - (frame.origin() - Point3::ORIGIN);
            let h = rel.dot(frame.normal());
            let rho = (rel - frame.normal() * h).length();
            ((rho - major).hypot(h) - minor).abs()
        }
        _ => unreachable!(),
    }
}
