//! Test-only input for compare_torus_curves.py: reads the case protocol of
//! `torus-curve-cases.txt` on stdin and prints per case the reference's rows
//! (`torus_curve_reference.py`), every number as `lo hi`: `NAME empty`, or
//! `fold` rows (the angles and the point), `node` rows (with `crossing` or
//! `isolated`), `component folds w_phi w_t` and `cluster folds nodes` rows,
//! and `ring` rows (the point at `phi = 0` of each component without folds).
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

/// An angle's enclosure moved by whole turns into `(-pi, pi]` by its middle.
fn canonical([lo, hi]: [f64; 2]) -> [f64; 2] {
    let m = 0.5 * lo + 0.5 * hi;
    let k = if m <= -PI + 1e-12 {
        1.0
    } else if m > PI + 1e-12 {
        -((m - PI) / (2.0 * PI)).ceil()
    } else {
        0.0
    };
    [lo + k * 2.0 * PI, hi + k * 2.0 * PI]
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
            Ok(SurfaceIntersection::Traced(c)) => c,
            other => {
                println!("{name} error {other:?}");
                continue;
            }
        };
        // By phi, then by t where the phi enclosures overlap (mirror images).
        let order = |x: &[[f64; 2]; 2], y: &[[f64; 2]; 2]| {
            if x[0][1] < y[0][0] {
                std::cmp::Ordering::Less
            } else if y[0][1] < x[0][0] {
                std::cmp::Ordering::Greater
            } else {
                (x[1][0] + x[1][1]).total_cmp(&(y[1][0] + y[1][1]))
            }
        };
        let mut folds: Vec<_> = c
            .folds()
            .iter()
            .map(|f| ([canonical(f.phi), canonical(f.t)], f.point))
            .collect();
        folds.sort_by(|x, y| order(&x.0, &y.0));
        for (angles, p) in folds {
            let mut v = angles.to_vec();
            v.extend(p);
            println!("{name} fold {}", text(&v));
        }
        let mut nodes: Vec<_> = c
            .nodes()
            .iter()
            .map(|n| ([canonical(n.phi), canonical(n.t)], n.point, n.crossing))
            .collect();
        nodes.sort_by(|x, y| order(&x.0, &y.0));
        for (angles, p, crossing) in nodes {
            let mut v = angles.to_vec();
            v.extend(p);
            let kind = if crossing { "crossing" } else { "isolated" };
            println!("{name} node {} {kind}", text(&v));
        }
        let mut comps = Vec::new();
        let mut rings = Vec::new();
        for comp in c.components() {
            match comp {
                TracedComponent::Smooth {
                    tracks,
                    folds,
                    winding,
                    ..
                } => {
                    comps.push(("component", vec![*folds as i64, winding[0], winding[1]]));
                    if *folds == 0 {
                        let track = tracks[0];
                        let t = canonical(c.t_at(track, 0.0).unwrap());
                        let p = c.point_at(track, 0.0).unwrap();
                        rings.push((t, p));
                    }
                }
                TracedComponent::Crossing { folds, nodes, .. } => {
                    comps.push(("cluster", vec![*folds as i64, nodes.len() as i64]));
                }
                TracedComponent::Isolated { .. } => {}
            }
        }
        comps.sort();
        for (kind, v) in comps {
            let words: Vec<String> = v.iter().map(|x| x.to_string()).collect();
            println!("{name} {kind} {}", words.join(" "));
        }
        rings.sort_by(|x, y| (x.0[0] + x.0[1]).partial_cmp(&(y.0[0] + y.0[1])).unwrap());
        for (t, p) in rings {
            let mut v = vec![[0.0, 0.0], t];
            v.extend(p);
            println!("{name} ring {}", text(&v));
        }
    }
}
