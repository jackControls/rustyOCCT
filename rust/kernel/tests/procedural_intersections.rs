//! S7b.1 and S7b.2: procedural intersection curves of two cylinders with
//! crossing axes, a cylinder and a sphere off its axis, a cylinder and a cone
//! (axes not coaxial) and a sphere and a cone (the centre off the axis),
//! against the independent reference (`fixtures/procedural-intersection-*.txt|tsv` from
//! `tools/generate_procedural_intersection_fixtures.py`).
use rusty_occt::intersection::{
    surface_surface, AnalyticItem, Branch, Component, ProceduralCurve, SurfaceIntersection,
};
use rusty_occt::topology::Surface;
use rusty_occt::{Frame3, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;

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
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        "cone" => Surface::Cone {
            frame,
            radius: v[9],
            half_angle: v[10],
        },
        _ => Surface::Sphere {
            frame,
            radius: v[9],
        },
    }
}

fn cases() -> Vec<(String, Surface, Surface)> {
    include_str!("../../fixtures/procedural-intersection-cases.txt")
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

/// Each case's rows: kind and numbers, in the reference's order.
fn expected() -> BTreeMap<String, Vec<(String, Vec<f64>)>> {
    let mut out: BTreeMap<String, Vec<(String, Vec<f64>)>> = BTreeMap::new();
    for line in include_str!("../../fixtures/procedural-intersection-expected.tsv")
        .lines()
        .skip(1)
    {
        let (name, row) = line.split_once('\t').unwrap();
        let mut w = row.split(' ');
        let kind = w.next().unwrap().to_string();
        let numbers: Vec<f64> = w.map(|x| x.parse().unwrap()).collect();
        out.entry(name.into()).or_default().push((kind, numbers));
    }
    out
}

fn contains(e: [f64; 2], x: f64) -> bool {
    let slack = 1e-25 * x.abs();
    e[0] - slack <= x && x <= e[1] + slack
}
fn contains3(e: [[f64; 2]; 3], x: &[f64]) -> bool {
    (0..3).all(|k| contains(e[k], x[k]))
}
fn close3(e: [[f64; 2]; 3], x: &[f64], scale: f64) -> bool {
    (0..3).all(|k| (0.5 * e[k][0] + 0.5 * e[k][1] - x[k]).abs() <= 1e-12 * scale)
}

/// The enclosure of pi as binary64 bounds.
const PI: [f64; 2] = [std::f64::consts::PI, 3.1415926535897936];

fn check(got: &SurfaceIntersection, rows: &[(String, Vec<f64>)]) -> Result<(), String> {
    let kinds: Vec<&str> = rows.iter().map(|(k, _)| k.as_str()).collect();
    match got {
        SurfaceIntersection::Empty if kinds == ["empty"] => Ok(()),
        SurfaceIntersection::Items(items) => match (items.as_slice(), rows) {
            ([AnalyticItem::Point(e)], [(k, p)]) if k == "point" && contains3(*e, p) => Ok(()),
            (other, _) => Err(format!("{other:?} for {kinds:?}")),
        },
        SurfaceIntersection::Procedural(c) => check_curve(c, rows),
        other => Err(format!("{other:?} for {kinds:?}")),
    }
}

fn check_curve(c: &ProceduralCurve, rows: &[(String, Vec<f64>)]) -> Result<(), String> {
    let at = |u: [f64; 2], b| c.point_at(u, b).map_err(|e| e.to_string());
    let row_kinds: Vec<&str> = rows.iter().map(|(k, _)| k.as_str()).collect();
    match c.components() {
        // Two rings are one row, by their points at 0 and pi.
        [Component::Ring {
            branch: Branch::Plus,
        }, Component::Ring {
            branch: Branch::Minus,
        }] => {
            let [(k, r)] = rows else {
                return Err(format!("rings for {row_kinds:?}"));
            };
            if k != "rings" {
                return Err(format!("rings for {k}"));
            }
            for (n, (u, b)) in [
                ([0.0, 0.0], Branch::Plus),
                ([0.0, 0.0], Branch::Minus),
                (PI, Branch::Plus),
                (PI, Branch::Minus),
            ]
            .into_iter()
            .enumerate()
            {
                if !contains3(at(u, b)?, &r[3 * n..3 * n + 3]) {
                    return Err(format!("ring point {n}"));
                }
            }
            Ok(())
        }
        [Component::FigureEight { node }] => {
            let [(k, r)] = rows else {
                return Err(format!("a figure-eight for {row_kinds:?}"));
            };
            if k != "figure_eight" {
                return Err(format!("a figure-eight for {k}"));
            }
            if !contains3(at(*node, Branch::Plus)?, &r[0..3]) {
                return Err("node".into());
            }
            // Opposite the node.
            if !contains3(at([0.0, 0.0], Branch::Plus)?, &r[3..6])
                || !contains3(at([0.0, 0.0], Branch::Minus)?, &r[6..9])
            {
                return Err("opposite the node".into());
            }
            Ok(())
        }
        // Each loop its own row, in order.
        loops => {
            if loops.len() != rows.len() || row_kinds.iter().any(|k| *k != "loop") {
                return Err(format!("{loops:?} for {row_kinds:?}"));
            }
            for (comp, (_, r)) in loops.iter().zip(rows) {
                let Component::Loop { u } = comp else {
                    return Err(format!("{comp:?} among loops"));
                };
                let (u0, u1) = (r[0], r[1]);
                if !contains(u[0], u0) || !contains(u[1], u1) {
                    return Err(format!("loop range {u:?} for {u0}, {u1}"));
                }
                // Tight: a loop's ends are narrowed to a few units in the
                // last place.
                if u.iter()
                    .any(|[lo, hi]| hi - lo > 1e-12 * lo.abs().max(hi.abs()).max(1.0))
                {
                    return Err(format!("wide loop range {u:?}"));
                }
                // The ends: the curve's point over each root's enclosure.
                if !contains3(at(u[0], Branch::Plus)?, &r[2..5])
                    || !contains3(at(u[1], Branch::Plus)?, &r[5..8])
                {
                    return Err("loop ends".into());
                }
                // The middle, at the reference's parameter rounded to binary64.
                let mid = 0.5 * u0 + 0.5 * u1;
                let scale = r.iter().fold(1.0f64, |a, x| a.max(x.abs()));
                if !close3(at([mid, mid], Branch::Plus)?, &r[8..11], scale)
                    || !close3(at([mid, mid], Branch::Minus)?, &r[11..14], scale)
                {
                    return Err("loop middle".into());
                }
            }
            Ok(())
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
        if let Err(why) = check(&got, &want[name]) {
            failures.push(format!("{name}: {why}"));
        }
        assert_eq!(surface_surface(b, a).unwrap(), got, "{name}: order");
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(all.len(), 35);
}

/// Points along every component lie on both surfaces.
#[test]
fn sampled_points_lie_on_both_surfaces() {
    for (name, a, b) in cases() {
        let SurfaceIntersection::Procedural(c) = surface_surface(&a, &b).unwrap() else {
            continue;
        };
        let mut params = Vec::new();
        for comp in c.components() {
            match comp {
                Component::Loop { u } => {
                    let (lo, hi) = (u[0][1], u[1][0]);
                    for k in 1..16 {
                        let t = lo + (hi - lo) * f64::from(k) / 16.0;
                        params.push((t, Branch::Plus));
                        params.push((t, Branch::Minus));
                    }
                }
                Component::Ring { branch } => {
                    for k in 0..16 {
                        params.push((-3.0 + 0.375 * f64::from(k), *branch));
                    }
                }
                Component::FigureEight { .. } => {
                    for k in 0..16 {
                        let t = -3.0 + 0.375 * f64::from(k);
                        params.push((t, Branch::Plus));
                        params.push((t, Branch::Minus));
                    }
                }
            }
        }
        for (t, branch) in params {
            let e = c.point_at([t, t], branch).unwrap();
            for k in 0..3 {
                assert!(e[k][1] - e[k][0] <= 1e-12, "{name}: wide {e:?}");
            }
            let p = Vec3::new(e[0][0], e[1][0], e[2][0]);
            for s in [c.ruled(), c.other()] {
                let gap = distance(s, p);
                assert!(gap <= 1e-12, "{name}: {gap} off {s:?} at {t}");
            }
        }
    }
}

fn distance(s: &Surface, p: Vec3) -> f64 {
    match s {
        Surface::Cylinder { frame, radius } => {
            let rel = p - (frame.origin() - Point3::ORIGIN);
            let h = rel.dot(frame.normal());
            ((rel - frame.normal() * h).length() - radius).abs()
        }
        Surface::Sphere { frame, radius } => {
            ((p - (frame.origin() - Point3::ORIGIN)).length() - radius).abs()
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
        _ => unreachable!(),
    }
}
