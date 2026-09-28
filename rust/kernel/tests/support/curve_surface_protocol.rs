//! The case protocol of `curve-surface-cases.txt` (S7c) and the kernel's
//! rows for it, shared by `curve_surface_probe` and `tests/curve_surface.rs`:
//! `empty`, `contained`, `limit`, `overlap s0 s1`, or `point s_lo s_hi x_lo
//! x_hi y_lo y_hi z_lo z_hi crossing|tangent`, sorted by parameter.
use rusty_occt::intersection::{
    conic_surface, curve_surface, spline_cone, spline_torus, Conic, CurvePoint,
    CurveSurfaceIntersection,
};
use rusty_occt::topology::{Curve3, Surface};
use rusty_occt::{BSplineCurve3, Error, Frame3, Point3, Tolerance, Vec3};

pub enum Curve {
    Edge(Curve3),
    Conic(Conic),
    Spline(BSplineCurve3),
}

pub struct Case {
    pub name: String,
    pub curve: Curve,
    pub surface: Surface,
}

fn frame(v: &[f64]) -> Frame3 {
    Frame3::new(
        Point3::new(v[0], v[1], v[2]),
        Vec3::new(v[3], v[4], v[5]),
        Vec3::new(v[6], v[7], v[8]),
        Tolerance::default(),
    )
    .unwrap()
}

fn numbers(words: &[&str]) -> Vec<f64> {
    words.iter().map(|w| w.parse().unwrap()).collect()
}

fn curve(words: &[&str]) -> Curve {
    let v = numbers(&words[1..]);
    match words[0] {
        "line" => Curve::Edge(Curve3::LineSegment {
            start: Point3::new(v[0], v[1], v[2]),
            end: Point3::new(v[3], v[4], v[5]),
        }),
        "circle" => Curve::Edge(Curve3::Circle {
            frame: frame(&v),
            radius: v[9],
        }),
        "ellipse" => Curve::Conic(Conic::Ellipse {
            frame: frame(&v),
            major: v[9],
            minor: v[10],
        }),
        "hyperbola" => Curve::Conic(Conic::Hyperbola {
            frame: frame(&v),
            major: v[9],
            minor: v[10],
        }),
        _ => {
            let (degree, n) = (v[0] as usize, v[1] as usize);
            let poles: Vec<[f64; 4]> = (0..n)
                .map(|i| std::array::from_fn(|c| v[2 + 4 * i + c]))
                .collect();
            let at = 2 + 4 * n;
            let k = v[at] as usize;
            let knots = (0..k).map(|i| v[at + 1 + 2 * i]).collect();
            let mults = (0..k).map(|i| v[at + 2 + 2 * i] as usize).collect();
            Curve::Spline(
                BSplineCurve3::new(
                    degree,
                    poles
                        .iter()
                        .map(|p| Point3::new(p[0], p[1], p[2]))
                        .collect(),
                    Some(poles.iter().map(|p| p[3]).collect()),
                    knots,
                    mults,
                )
                .unwrap(),
            )
        }
    }
}

pub fn surface(words: &[&str]) -> Surface {
    let v = numbers(&words[1..]);
    let frame = frame(&v);
    match words[0] {
        "plane" => Surface::Plane(frame),
        "sphere" => Surface::Sphere {
            frame,
            radius: v[9],
        },
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        "torus" => Surface::Torus {
            frame,
            major: v[9],
            minor: v[10],
        },
        _ => Surface::Cone {
            frame,
            radius: v[9],
            half_angle: v[10],
        },
    }
}

pub fn cases(text: &str) -> Vec<Case> {
    text.split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|block| {
            let lines: Vec<Vec<&str>> = block
                .trim()
                .lines()
                .map(|l| l.split_whitespace().collect())
                .collect();
            Case {
                name: lines[0][1].to_string(),
                curve: curve(&lines[1][1..]),
                surface: surface(&lines[2][1..]),
            }
        })
        .collect()
}

fn point_row(p: &CurvePoint) -> String {
    let mut words = vec![format!("{:?} {:?}", p.parameter[0], p.parameter[1])];
    words.extend(p.point.iter().map(|[lo, hi]| format!("{lo:?} {hi:?}")));
    let contact = if p.tangent { "tangent" } else { "crossing" };
    format!("point {} {contact}", words.join(" "))
}

/// The kernel's rows for a case; `Err` for an unexpected error.
pub fn rows(case: &Case) -> Result<Vec<String>, Error> {
    let result = match &case.curve {
        Curve::Edge(c) => curve_surface(c, &case.surface),
        Curve::Conic(c) => conic_surface(c, &case.surface),
        Curve::Spline(c) => return spline_rows(c, &case.surface),
    };
    Ok(match result {
        Ok(CurveSurfaceIntersection::Empty) => vec!["empty".into()],
        Ok(CurveSurfaceIntersection::Contained) => vec!["contained".into()],
        Ok(CurveSurfaceIntersection::Points(points)) => points.iter().map(point_row).collect(),
        Err(Error::ComputationLimit(_)) => vec!["limit".into()],
        Err(e) => return Err(e),
    })
}

fn spline_rows(c: &BSplineCurve3, s: &Surface) -> Result<Vec<String>, Error> {
    let mut rows: Vec<(f64, String)> = Vec::new();
    match s {
        Surface::Torus { .. } => {
            let found = match spline_torus(c, s) {
                Err(Error::ComputationLimit(_)) => return Ok(vec!["limit".into()]),
                other => other?,
            };
            for o in found.overlaps() {
                let (a, b) = o.parameters();
                rows.push((a, format!("overlap {a:?} {b:?}")));
            }
            for p in found.points() {
                let t = p.parameter();
                let point = CurvePoint {
                    parameter: [t.lower(), t.upper()],
                    point: p.coordinate_bounds().map(|x| [x.lower(), x.upper()]),
                    tangent: p.multiplicities().iter().flatten().any(|m| *m > 1),
                };
                rows.push((t.lower(), point_row(&point)));
            }
        }
        _ => {
            let found = match spline_cone(c, s) {
                Err(Error::ComputationLimit(_)) => return Ok(vec!["limit".into()]),
                other => other?,
            };
            for [a, b] in found.overlaps {
                rows.push((a, format!("overlap {a:?} {b:?}")));
            }
            for p in &found.points {
                rows.push((p.parameter[0], point_row(p)));
            }
        }
    }
    rows.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(if rows.is_empty() {
        vec!["empty".into()]
    } else {
        rows.into_iter().map(|r| r.1).collect()
    })
}
