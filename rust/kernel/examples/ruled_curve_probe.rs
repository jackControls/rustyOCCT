//! Test-only input for compare_ruled_curves.py: reads the case protocol of
//! `ruled-curve-cases.txt` on stdin and prints per case the reference's rows
//! (`ruled_curve_reference.py`), every number as `lo hi`: `NAME empty`; `line`
//! and `point` items; `fold` rows (the angles and the point); `apex` rows
//! (with `crossing` or `isolated`); `component folds w_u w_t infinite` rows;
//! and `ring` rows (the point at `u = 0` of each component without folds, at
//! `u = 1` on a factor's chart).
use rusty_occt::intersection::{surface_surface, SurfaceIntersection, TracedComponent};
use rusty_occt::topology::Surface;
use rusty_occt::{Frame3, Point3, Tolerance, Vec3};
use std::f64::consts::PI;
use std::io::Read;

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

/// An angle's enclosure moved by whole periods into `(-p/2, p/2]` by its
/// middle.
fn canonical([lo, hi]: [f64; 2], period: f64) -> [f64; 2] {
    let m = 0.5 * lo + 0.5 * hi;
    let k = -((m - 0.5 * period) / period).ceil();
    let k = if m + k * period <= -0.5 * period + 1e-12 {
        k + 1.0
    } else {
        k
    };
    [lo + k * period, hi + k * period]
}

fn text(values: &[[f64; 2]]) -> String {
    values
        .iter()
        .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for block in input.split("\nend").filter(|b| !b.trim().is_empty()) {
        let lines: Vec<Vec<&str>> = block
            .trim()
            .lines()
            .map(|l| l.split_whitespace().collect())
            .collect();
        let name = lines[0][1];
        let (a, b) = (surface(&lines[1][1..]), surface(&lines[2][1..]));
        let c = match surface_surface(&a, &b) {
            Ok(SurfaceIntersection::Empty) => {
                println!("{name} empty");
                continue;
            }
            Ok(SurfaceIntersection::Items(items)) => {
                for item in items {
                    println!("{name} {} {}", item.kind(), text(&item.values()));
                }
                continue;
            }
            Ok(SurfaceIntersection::Traced(c)) => c,
            other => {
                println!("{name} error {other:?}");
                continue;
            }
        };
        let period = c.period();
        let mut folds: Vec<_> = c
            .folds()
            .iter()
            .map(|f| {
                (
                    [canonical(f.phi, 2.0 * PI), canonical(f.t, period)],
                    f.point,
                )
            })
            .collect();
        // By u, then by t where the u enclosures overlap.
        folds.sort_by(|x, y| {
            let (p, q) = (x.0, y.0);
            if p[0][1] < q[0][0] {
                std::cmp::Ordering::Less
            } else if q[0][1] < p[0][0] {
                std::cmp::Ordering::Greater
            } else {
                (p[1][0] + p[1][1]).total_cmp(&(q[1][0] + q[1][1]))
            }
        });
        for (angles, p) in folds {
            let mut v = angles.to_vec();
            v.extend(p);
            println!("{name} fold {}", text(&v));
        }
        for n in c.nodes() {
            let kind = if n.crossing { "crossing" } else { "isolated" };
            println!("{name} apex {} {kind}", text(&n.point));
        }
        let at = if period < 4.0 { 1.0 } else { 0.0 };
        let mut comps = Vec::new();
        let mut rings = Vec::new();
        for comp in c.components() {
            let (tracks, folds, winding, infinite) = match comp {
                TracedComponent::Smooth {
                    tracks,
                    folds,
                    winding,
                    infinite,
                } => (tracks, *folds, *winding, *infinite),
                // Through the apex: a factor's ring, one turn in u.
                TracedComponent::Crossing {
                    tracks, infinite, ..
                } => (tracks, 0, [1, winding_of(&c, tracks)], *infinite),
                TracedComponent::Isolated { .. } => continue,
            };
            comps.push(vec![folds as i64, winding[0], winding[1], infinite as i64]);
            if folds == 0 {
                let track = tracks
                    .iter()
                    .copied()
                    .find(|t| c.t_at(*t, at).is_ok())
                    .unwrap();
                let t = canonical(c.t_at(track, at).unwrap(), period);
                rings.push((t, c.point_at(track, at).unwrap()));
            }
        }
        comps.sort();
        for v in comps {
            let words: Vec<String> = v.iter().map(|x| x.to_string()).collect();
            println!("{name} component {}", words.join(" "));
        }
        rings.sort_by(|x, y| (x.0[0] + x.0[1]).total_cmp(&(y.0[0] + y.0[1])));
        for (t, p) in rings {
            let mut v = vec![[at, at], t];
            v.extend(p);
            println!("{name} ring {}", text(&v));
        }
    }
}

/// The winding in the chart's second parameter of a component through the
/// apex (its tracks cover one turn in u): the change of t over them.
fn winding_of(c: &rusty_occt::intersection::TracedCurve, tracks: &[usize]) -> i64 {
    let mut total = 0.0;
    for &k in tracks {
        let tr = &c.tracks()[k];
        let (a, b) = (c.t_at(k, tr.phi[0]).unwrap(), c.t_at(k, tr.phi[1]).unwrap());
        total += 0.5 * (b[0] + b[1]) - 0.5 * (a[0] + a[1]);
    }
    (total / c.period()).round() as i64
}
