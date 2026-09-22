//! Test-only text worker for native OCCT comparisons of roots and intersections.
use rusty_occt::intersection::{line_circle, line_cylinder, line_sphere, CurvedIntersection};
use rusty_occt::polynomial::{quadratic_roots, QuadraticRoots};
use rusty_occt::{Circle3, Cylinder3, Point3, Sphere3, Vec3};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for row in io::stdin().lock().lines() {
        let row = row?;
        let words: Vec<_> = row.split_whitespace().collect();
        if words.len() < 2 {
            return Err("missing fixture kind/label".into());
        }
        let values = words[2..]
            .iter()
            .map(|x| x.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()?;
        print!("{}", words[1]);
        if words[0] == "q" {
            if values.len() != 3 {
                return Err("quadratic requires three coefficients".into());
            }
            let roots = match quadratic_roots(values[0], values[1], values[2])? {
                QuadraticRoots::All => {
                    println!(" A");
                    continue;
                }
                QuadraticRoots::None => vec![],
                QuadraticRoots::One(root) => vec![root],
                QuadraticRoots::Two { lower, upper } => vec![lower, upper],
            };
            let count: u8 = roots.iter().map(|r| r.multiplicity()).sum();
            print!(" {count}");
            for root in roots {
                let value = root.bounds()?.representative();
                for _ in 0..root.multiplicity() {
                    print!(" {value:.17e}");
                }
            }
        } else {
            if values.len() != 13 {
                return Err("curved primitive requires 13 values".into());
            }
            let center = Point3::new(values[0], values[1], values[2]);
            let axis = Vec3::new(values[3], values[4], values[5]);
            let radius = values[6];
            let p = Point3::new(values[7], values[8], values[9]);
            let q = Point3::new(values[10], values[11], values[12]);
            let result = match words[0] {
                "s" => line_sphere(p, q, &Sphere3::new(center, radius)?)?,
                "y" => line_cylinder(p, q, &Cylinder3::new(center, axis, radius)?)?,
                "c" => line_circle(p, q, &Circle3::new(center, axis, radius)?)?,
                _ => return Err("unknown curved primitive".into()),
            };
            let hits = match result {
                CurvedIntersection::Contained => {
                    println!(" A");
                    continue;
                }
                CurvedIntersection::Disjoint => vec![],
                CurvedIntersection::One(hit) => vec![hit],
                CurvedIntersection::Two { first, second } => vec![first, second],
            };
            print!(" {}", hits.len());
            for hit in hits {
                let p = hit.point.position();
                print!(
                    " {:.17e} {:.17e} {:.17e} {:.17e}",
                    hit.point.parameter().representative(),
                    p.x,
                    p.y,
                    p.z
                );
            }
        }
        println!();
    }
    Ok(())
}
