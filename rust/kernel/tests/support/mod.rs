//! Test-only text protocol shared with rust/tools/occt_oracle.cpp.
use rusty_occt::identity::OperationId;
use rusty_occt::{Boundary, Frame3, Location, Point2, Point3, Profile, Solid, Tolerance, Vec3};
use std::fmt;
use std::str::{FromStr, SplitWhitespace};

struct Input<'a>(SplitWhitespace<'a>);
impl Input<'_> {
    fn token<T: FromStr>(&mut self) -> Result<T, Box<dyn std::error::Error>> {
        self.0
            .next()
            .ok_or("missing fixture token")?
            .parse()
            .map_err(|_| "invalid fixture token".into())
    }
    fn point(&mut self) -> Result<Point3, Box<dyn std::error::Error>> {
        Ok(Point3::new(self.token()?, self.token()?, self.token()?))
    }
    fn vector(&mut self) -> Result<Vec3, Box<dyn std::error::Error>> {
        Ok(Vec3::new(self.token()?, self.token()?, self.token()?))
    }
    fn boundary(&mut self, tolerance: Tolerance) -> Result<Boundary, Box<dyn std::error::Error>> {
        match self.token::<String>()?.as_str() {
            "P" => {
                let count = self.token::<usize>()?;
                let mut points = Vec::new();
                for _ in 0..count {
                    points.push(Point2::new(self.token()?, self.token()?));
                }
                Ok(Boundary::polygon(points, tolerance)?)
            }
            "C" => Ok(Boundary::circle(
                Point2::new(self.token()?, self.token()?),
                self.token()?,
                tolerance,
            )?),
            _ => Err("unknown boundary fixture kind".into()),
        }
    }
}

#[derive(Debug)]
pub struct Observation {
    pub name: String,
    pub values: Vec<f64>,
}

impl fmt::Display for Observation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name)?;
        for value in &self.values {
            write!(f, " {value:.17e}")?;
        }
        Ok(())
    }
}

/// Every fixture solid by name, as built by `evaluate`.
#[allow(dead_code)]
pub fn solids(text: &str) -> Result<Vec<(String, Solid)>, Box<dyn std::error::Error>> {
    let mut input = Input(text.split_whitespace());
    let count = input.token::<usize>()?;
    let tolerance = Tolerance::default();
    let mut out = Vec::new();
    for _ in 0..count {
        let name: String = input.token()?;
        let boundary_count: usize = input.token()?;
        let frame = Frame3::new(input.point()?, input.vector()?, input.vector()?, tolerance)?;
        let (start, end) = (input.token()?, input.token()?);
        let outer = input.boundary(tolerance)?;
        let mut holes = Vec::new();
        for _ in 1..boundary_count {
            holes.push(input.boundary(tolerance)?);
        }
        let profile = Profile::new(outer, holes, tolerance)?;
        let solid = Solid::extrude_with(OperationId::UNSPECIFIED, profile, frame, start, end)?.0;
        let queries: usize = input.token()?;
        for _ in 0..queries {
            input.point()?;
        }
        out.push((name, solid));
    }
    Ok(out)
}

pub fn evaluate(text: &str) -> Result<Vec<Observation>, Box<dyn std::error::Error>> {
    let mut input = Input(text.split_whitespace());
    let count = input.token::<usize>()?;
    let tolerance = Tolerance::default();
    let mut observations = Vec::new();
    for _ in 0..count {
        let name: String = input.token()?;
        let boundary_count: usize = input.token()?;
        if boundary_count == 0 {
            return Err("fixture needs an outer boundary".into());
        }
        let frame = Frame3::new(input.point()?, input.vector()?, input.vector()?, tolerance)?;
        let (start, end) = (input.token()?, input.token()?);
        let outer = input.boundary(tolerance)?;
        let mut holes = Vec::new();
        for _ in 1..boundary_count {
            holes.push(input.boundary(tolerance)?);
        }
        let solid = Solid::extrude_with(
            OperationId::UNSPECIFIED,
            Profile::new(outer, holes, tolerance)?,
            frame,
            start,
            end,
        )
        .map(|(s, _)| s)?;
        let mass = solid.mass_properties();
        let bounds = solid.bounds();
        let counts = solid.topology().occt_counts();
        let mut values = vec![mass.volume, mass.surface_area];
        values.extend(
            mass.centroid
                .to_array()
                .into_iter()
                .chain(bounds.min.to_array())
                .chain(bounds.max.to_array())
                .chain(mass.inertia.into_iter().flatten()),
        );
        values.extend([
            // The counts OCCT reports for the same body, seams synthesized.
            counts.vertices as f64,
            counts.edges as f64,
            counts.faces as f64,
        ]);
        let queries: usize = input.token()?;
        for _ in 0..queries {
            let local = input.point()?;
            let location = solid.classify(frame.point(Point2::new(local.x, local.y), local.z))?;
            values.push(match location {
                Location::Outside => 0.0,
                Location::Boundary => 1.0,
                Location::Inside => 2.0,
            });
        }
        observations.push(Observation { name, values });
    }
    if input.0.next().is_some() {
        return Err("unexpected trailing fixture tokens".into());
    }
    Ok(observations)
}
