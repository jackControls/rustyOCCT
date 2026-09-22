//! Line/plane differential-test worker. Input: label + five xyz points.
use rusty_occt::intersection::{line_plane, LinePlaneIntersection};
use rusty_occt::{Plane3, Point3};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for row in io::stdin().lock().lines() {
        let row = row?;
        let words: Vec<_> = row.split_whitespace().collect();
        if words.len() != 16 {
            return Err("expected label and five xyz points".into());
        }
        let coordinates = words[1..]
            .iter()
            .map(|x| x.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()?;
        let p: Vec<_> = coordinates
            .chunks_exact(3)
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect();
        let plane = Plane3::through_points(p[0], p[1], p[2])?;
        match line_plane(p[3], p[4], &plane)? {
            LinePlaneIntersection::Disjoint => println!("{} N", words[0]),
            LinePlaneIntersection::Contained => println!("{} C", words[0]),
            LinePlaneIntersection::Point(point) => {
                let p = point.position();
                println!(
                    "{} P {:.17e} {:.17e} {:.17e} {:.17e}",
                    words[0],
                    point.parameter().representative(),
                    p.x,
                    p.y,
                    p.z
                );
            }
        }
    }
    Ok(())
}
