//! S7a: intersections of analytic surfaces against the independent exact
//! reference (`fixtures/analytic-intersection-*.txt|tsv` from
//! `tools/generate_analytic_intersection_fixtures.py`).
use rusty_occt::intersection::{surface_surface, AnalyticItem, SurfaceIntersection};
use rusty_occt::topology::Surface;
use rusty_occt::{Frame3, Point3, RigidTransform, Tolerance, Vec3};
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
        "plane" => Surface::Plane(frame),
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        "sphere" => Surface::Sphere {
            frame,
            radius: v[9],
        },
        "cone" => Surface::Cone {
            frame,
            radius: v[9],
            half_angle: v[10],
        },
        other => panic!("surface kind {other}"),
    }
}

fn cases() -> Vec<(String, Surface, Surface)> {
    include_str!("../../fixtures/analytic-intersection-cases.txt")
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

fn expected() -> BTreeMap<String, Vec<String>> {
    let mut out: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in include_str!("../../fixtures/analytic-intersection-expected.tsv")
        .lines()
        .skip(1)
    {
        let (name, item) = line.split_once('\t').unwrap();
        out.entry(name.into()).or_default().push(item.into());
    }
    out
}

/// Whether the kernel's result has the reference's items, each reference
/// number inside the kernel's enclosure (up to 1e-25 relative, the
/// reference's printing).
fn matches(got: &SurfaceIntersection, want: &[String]) -> Result<(), String> {
    let words: Vec<&str> = want.iter().map(|w| w.as_str()).collect();
    match (got, words.as_slice()) {
        (SurfaceIntersection::Empty, ["empty"])
        | (SurfaceIntersection::Same, ["same"])
        | (SurfaceIntersection::NotConic, ["not_conic"])
        // A procedural curve (S7b) is not a conic.
        | (SurfaceIntersection::Procedural(_), ["not_conic"]) => Ok(()),
        (SurfaceIntersection::Items(items), _) if items.len() == want.len() => {
            for (item, row) in items.iter().zip(want) {
                let mut w = row.split(' ');
                let kind = w.next().unwrap();
                if item.kind() != kind {
                    return Err(format!("{} for {kind}", item.kind()));
                }
                let numbers: Vec<f64> = w.map(|x| x.parse().unwrap()).collect();
                let values = item.values();
                if values.len() != numbers.len() {
                    return Err(format!("{} numbers for {}", values.len(), numbers.len()));
                }
                for (k, ([lo, hi], x)) in values.iter().zip(&numbers).enumerate() {
                    let slack = 1e-25 * x.abs().max(1e-300);
                    if !(lo - slack <= *x && *x <= hi + slack) {
                        return Err(format!("{kind} value {k}: {x} outside [{lo}, {hi}]"));
                    }
                }
            }
            Ok(())
        }
        _ => Err(format!("{got:?} for {want:?}")),
    }
}

/// The kernel stores the frames the reference intersected, bit for bit.
#[test]
fn stored_normals_are_the_reference_inputs() {
    let all = cases();
    let rows: Vec<&str> = include_str!("../../fixtures/analytic-intersection-frames.tsv")
        .lines()
        .skip(1)
        .collect();
    assert_eq!(rows.len(), 2 * all.len());
    for ((name, a, b), pair) in all.iter().zip(rows.chunks(2)) {
        for (s, row) in [a, b].into_iter().zip(pair) {
            let frame = match s {
                Surface::Plane(f)
                | Surface::Cylinder { frame: f, .. }
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
    let mut failures = Vec::new();
    let all = cases();
    for (name, a, b) in &all {
        let got = surface_surface(a, b).unwrap_or_else(|e| panic!("{name}: {e}"));
        if let Err(why) = matches(&got, &want[name]) {
            failures.push(format!("{name}: {why}"));
        }
        // Symmetric in its arguments.
        assert_eq!(surface_surface(b, a).unwrap(), got, "{name}: order");
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(all.len(), 67);
}

/// Enclosures are tight: every bound is within a few units of its
/// midpoint's last place, except the huge near-degenerate parameters.
#[test]
fn enclosures_are_a_few_ulps_wide() {
    for (name, a, b) in cases() {
        if let SurfaceIntersection::Items(items) = surface_surface(&a, &b).unwrap() {
            for item in items {
                for [lo, hi] in item.values() {
                    let ulp = f64::EPSILON * lo.abs().max(hi.abs()).max(f64::MIN_POSITIVE);
                    assert!(hi - lo <= 4.0 * ulp, "{name}: [{lo}, {hi}]");
                }
            }
        }
    }
}

/// Every point of a returned curve (at a few parameters) lies on both
/// surfaces within a relative 1e-12.
#[test]
fn returned_curves_lie_on_both_surfaces() {
    let mid = |[lo, hi]: [f64; 2]| 0.5 * lo + 0.5 * hi;
    let v = |x: &[[f64; 2]; 3]| Vec3::new(mid(x[0]), mid(x[1]), mid(x[2]));
    for (name, a, b) in cases() {
        let SurfaceIntersection::Items(items) = surface_surface(&a, &b).unwrap() else {
            continue;
        };
        for item in items {
            let points: Vec<Vec3> = match &item {
                AnalyticItem::Point(p) => vec![v(p)],
                AnalyticItem::Line { point, direction } => (-2..=2)
                    .map(|k| v(point) + v(direction) * f64::from(k))
                    .collect(),
                AnalyticItem::Circle {
                    centre,
                    normal,
                    radius,
                } => {
                    let n = v(normal);
                    let x = any_perpendicular(n);
                    let y = n.cross(x);
                    (0..8)
                        .map(|k| {
                            let t = f64::from(k) * 0.785;
                            v(centre) + (x * t.cos() + y * t.sin()) * mid(*radius)
                        })
                        .collect()
                }
                AnalyticItem::Ellipse {
                    centre,
                    normal,
                    major,
                    semi_major,
                    semi_minor,
                } => {
                    let (u, w) = (v(major), v(normal).cross(v(major)));
                    (0..8)
                        .map(|k| {
                            let t = f64::from(k) * 0.785;
                            v(centre)
                                + u * (mid(*semi_major) * t.cos())
                                + w * (mid(*semi_minor) * t.sin())
                        })
                        .collect()
                }
                AnalyticItem::Hyperbola {
                    centre,
                    normal,
                    transverse,
                    semi_transverse,
                    semi_conjugate,
                } => {
                    let (u, w) = (v(transverse), v(normal).cross(v(transverse)));
                    (-2..=2)
                        .flat_map(|k| {
                            let t = f64::from(k) * 0.5;
                            [1.0, -1.0].map(|s| {
                                v(centre)
                                    + u * (s * mid(*semi_transverse) * t.cosh())
                                    + w * (mid(*semi_conjugate) * t.sinh())
                            })
                        })
                        .collect()
                }
            };
            let scale = points.iter().map(|p| p.length()).fold(1.0, f64::max);
            for p in points {
                for s in [&a, &b] {
                    let gap = distance(s, p);
                    assert!(gap <= 1e-12 * scale, "{name}: {gap} off {s:?}");
                }
            }
        }
    }
}

/// Rigid motions move every item with the surfaces (translation by an
/// exact dyadic vector keeps every exact predicate).
#[test]
fn translations_move_the_items() {
    let shift = Vec3::new(0.5, -0.25, 2.0);
    let t = RigidTransform::translation(shift).unwrap();
    let mv = |s: &Surface| -> Surface {
        let tol = Tolerance::default();
        match s {
            Surface::Plane(f) => Surface::Plane(f.transformed(t, tol).unwrap()),
            Surface::Cylinder { frame, radius } => Surface::Cylinder {
                frame: frame.transformed(t, tol).unwrap(),
                radius: *radius,
            },
            Surface::Sphere { frame, radius } => Surface::Sphere {
                frame: frame.transformed(t, tol).unwrap(),
                radius: *radius,
            },
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => Surface::Cone {
                frame: frame.transformed(t, tol).unwrap(),
                radius: *radius,
                half_angle: *half_angle,
            },
            _ => unreachable!(),
        }
    };
    for (name, a, b) in cases() {
        let before = surface_surface(&a, &b).unwrap();
        let after = surface_surface(&mv(&a), &mv(&b)).unwrap();
        match (&before, &after) {
            (SurfaceIntersection::Items(x), SurfaceIntersection::Items(y)) => {
                assert_eq!(x.len(), y.len(), "{name}");
                for (p, q) in x.iter().zip(y) {
                    assert_eq!(p.kind(), q.kind(), "{name}");
                }
            }
            // A procedural curve carries its (moved) surfaces: its
            // components are the same.
            (SurfaceIntersection::Procedural(x), SurfaceIntersection::Procedural(y)) => {
                assert_eq!(x.components().len(), y.components().len(), "{name}");
            }
            _ => assert_eq!(before, after, "{name}"),
        }
    }
}

fn any_perpendicular(n: Vec3) -> Vec3 {
    let e = if n.x.abs() < 0.5 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let x = n.cross(e);
    x * (1.0 / x.length())
}

/// Binary64 distance from a point to a surface (a test's check only).
fn distance(s: &Surface, p: Vec3) -> f64 {
    let (f, rho) = match s {
        Surface::Plane(f) => return (p - (f.origin() - Point3::ORIGIN)).dot(f.normal()).abs(),
        Surface::Cylinder { frame, radius } => (frame, Some((*radius, 0.0))),
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => (frame, Some((*radius, *half_angle))),
        Surface::Sphere { frame, radius } => {
            return ((p - (frame.origin() - Point3::ORIGIN)).length() - radius).abs()
        }
        _ => unreachable!(),
    };
    let (r, angle) = rho.unwrap();
    let rel = p - (f.origin() - Point3::ORIGIN);
    let h = rel.dot(f.normal());
    let radial = (rel - f.normal() * h).length();
    if angle == 0.0 {
        return (radial - r).abs();
    }
    let (s, c) = angle.sin_cos();
    (radial * c - (r * c + h * s))
        .abs()
        .min((radial * c + (r * c + h * s)).abs())
}
