//! S7b.3b.1: a torus and a cylinder or a cone off its axis, as traced curves,
//! against the independent reference (`fixtures/torus-curve-*.txt|tsv` from
//! `tools/generate_torus_curve_fixtures.py`).
use rusty_occt::intersection::{
    surface_surface, SurfaceIntersection, TracedComponent, TracedCurve,
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
        "torus" => Surface::Torus {
            frame,
            major: v[9],
            minor: v[10],
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
    include_str!("../../fixtures/torus-curve-cases.txt")
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
    for line in include_str!("../../fixtures/torus-curve-expected.tsv")
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
    // Folds: each reference fold inside exactly one of the kernel's.
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
    let nodes: Vec<_> = want("node").collect();
    if nodes.len() != c.nodes().len() {
        return Err(format!("{} nodes for {}", c.nodes().len(), nodes.len()));
    }
    for r in &nodes {
        let x = numbers(r);
        let crossing = r.last().unwrap() == "crossing";
        let hits = c
            .nodes()
            .iter()
            .filter(|n| {
                n.crossing == crossing
                    && contains_angle(n.phi, x[0])
                    && contains_angle(n.t, x[1])
                    && (0..3).all(|k| contains(n.point[k], x[2 + k]))
            })
            .count();
        if hits != 1 {
            return Err(format!("{name}: node {x:?} in {hits}"));
        }
    }
    // Components as sorted words.
    let mut got: Vec<String> = Vec::new();
    for comp in c.components() {
        match comp {
            TracedComponent::Smooth { folds, winding, .. } => {
                got.push(format!("component {folds} {} {}", winding[0], winding[1]))
            }
            TracedComponent::Crossing { folds, nodes, .. } => {
                got.push(format!("cluster {folds} {}", nodes.len()))
            }
            TracedComponent::Isolated { .. } => {}
        }
    }
    got.sort();
    let mut comps: Vec<String> = rows
        .iter()
        .filter(|r| r[0] == "component" || r[0] == "cluster")
        .map(|r| r.join(" "))
        .collect();
    comps.sort();
    if got != comps {
        return Err(format!("{got:?} for {comps:?}"));
    }
    // Rings: each reference point at phi = 0 on some fold-free component.
    for r in want("ring") {
        let x = numbers(r);
        let hit = c.components().iter().any(|comp| match comp {
            TracedComponent::Smooth {
                tracks, folds: 0, ..
            } => tracks.iter().any(|t| {
                c.point_at(*t, 0.0)
                    .is_ok_and(|p| (0..3).all(|k| contains(p[k], x[2 + k])))
            }),
            _ => false,
        });
        if !hit {
            return Err(format!("{name}: ring {x:?}"));
        }
    }
    Ok(())
}

/// The kernel stores the frames the reference intersected, bit for bit.
#[test]
fn stored_normals_are_the_reference_inputs() {
    let all = cases();
    let rows: Vec<&str> = include_str!("../../fixtures/torus-curve-frames.tsv")
        .lines()
        .skip(1)
        .collect();
    assert_eq!(rows.len(), 2 * all.len());
    for ((name, a, b), pair) in all.iter().zip(rows.chunks(2)) {
        for (s, row) in [a, b].into_iter().zip(pair) {
            let frame = match s {
                Surface::Cylinder { frame: f, .. }
                | Surface::Cone { frame: f, .. }
                | Surface::Torus { frame: f, .. } => f,
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
            other => Err(format!("{other:?}")),
        };
        if let Err(why) = verdict {
            failures.push(format!("{name}: {why}"));
        }
        assert_eq!(surface_surface(b, a).unwrap(), got, "{name}: order");
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(all.len(), 24);
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
                let e = c.point_at(k, phi).unwrap();
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
    assert!(count > 250, "{count} points");
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
