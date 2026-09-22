//! Test-only DRAW command adapter. Production code never depends on Tcl or OCCT.
//! One persistent shape table per test process; no expected values live here.
use rusty_occt::{Point3, RigidTransform, Solid, Tolerance, Vec3};
use std::collections::BTreeMap;
use std::io::{self, BufRead, Write};

enum Failure {
    Error(String),
    Unsupported(String),
}
impl From<rusty_occt::Error> for Failure {
    fn from(value: rusty_occt::Error) -> Self {
        Self::Error(value.to_string())
    }
}
type Result<T> = std::result::Result<T, Failure>;

fn unsupported(args: &[String]) -> Failure {
    Failure::Unsupported(format!("unsupported DRAW signature: {}", args.join(" ")))
}
fn numbers(args: &[String]) -> Result<Vec<f64>> {
    args.iter()
        .map(|s| {
            s.parse::<f64>()
                .ok()
                .filter(|n| n.is_finite())
                .ok_or_else(|| Failure::Error(format!("invalid finite number: {s}")))
        })
        .collect()
}
fn get<'a>(shapes: &'a BTreeMap<String, Solid>, name: &str) -> Result<&'a Solid> {
    shapes
        .get(name)
        .ok_or_else(|| Failure::Error(format!("unknown shape: {name}")))
}

fn dispatch(shapes: &mut BTreeMap<String, Solid>, args: &[String]) -> Result<String> {
    let t = Tolerance::default();
    let command = args.first().map(String::as_str).unwrap_or("");
    match command {
        "box" if args.len() == 5 || args.len() == 8 => {
            let n = numbers(&args[2..])?;
            let (origin, size) = if n.len() == 3 {
                (Point3::ORIGIN, Vec3::new(n[0], n[1], n[2]))
            } else {
                (Point3::new(n[0], n[1], n[2]), Vec3::new(n[3], n[4], n[5]))
            };
            let solid = Solid::box_at(origin, size, t)?;
            shapes.insert(args[1].clone(), solid);
            Ok(String::new())
        }
        "copy" if args.len() == 3 => {
            let solid = get(shapes, &args[1])?.clone();
            shapes.insert(args[2].clone(), solid);
            Ok(String::new())
        }
        "ttranslate" | "trotate" if args.len() == (if command == "ttranslate" { 5 } else { 9 }) => {
            let n = numbers(&args[2..])?;
            let transform = if command == "ttranslate" {
                RigidTransform::translation(Vec3::new(n[0], n[1], n[2]))?
            } else {
                RigidTransform::rotation(
                    Point3::new(n[0], n[1], n[2]),
                    Vec3::new(n[3], n[4], n[5]),
                    n[6].to_radians(),
                )?
            };
            let solid = get(shapes, &args[1])?.transformed(transform)?;
            shapes.insert(args[1].clone(), solid);
            Ok(String::new())
        }
        "isdraw" if args.len() == 2 => Ok(if shapes.contains_key(&args[1]) {
            "1"
        } else {
            "0"
        }
        .into()),
        "whatis" if args.len() == 2 => {
            get(shapes, &args[1])?;
            Ok(format!(
                "{} is a shape SOLID FORWARD Free Modified",
                args[1]
            ))
        }
        "checkshape" if args.len() == 2 => {
            get(shapes, &args[1])?.topology().validate(t)?;
            Ok("This shape seems to be valid".into())
        }
        "nbshapes" if args.len() == 2 => {
            let topology = get(shapes, &args[1])?.topology();
            let counts = [
                ("VERTEX", topology.vertices().len()),
                ("EDGE", topology.edges().len()),
                ("WIRE", topology.faces().iter().map(|f| f.loops.len()).sum()),
                ("FACE", topology.faces().len()),
                ("SHELL", 1),
                ("SOLID", 1),
                ("COMPSOLID", 0),
                ("COMPOUND", 0),
            ];
            let mut result = format!("Number of shapes in {}\n", args[1]);
            for (kind, count) in counts {
                result.push_str(&format!(" {kind:<10}: {count}\n"));
            }
            result.push_str(&format!(
                " SHAPE     : {}\n",
                counts.iter().map(|(_, n)| n).sum::<usize>()
            ));
            Ok(result)
        }
        "vprops" if args.len() == 2 || args.len() == 3 => {
            if args.len() == 3 && numbers(&args[2..])?[0] <= 0.0 {
                return Err(unsupported(args));
            }
            // Epsilon controls OCCT quadrature. These box properties are analytic.
            let m = get(shapes, &args[1])?.mass_properties();
            Ok(format!("Mass : {:.17e}\n\nCenter of gravity :\nX = {:.17e}\nY = {:.17e}\nZ = {:.17e}\nMatrix of Inertia :\n{:.17e} {:.17e} {:.17e}\n{:.17e} {:.17e} {:.17e}\n{:.17e} {:.17e} {:.17e}\n",
                m.volume, m.centroid.x, m.centroid.y, m.centroid.z,
                m.inertia[0][0], m.inertia[0][1], m.inertia[0][2],
                m.inertia[1][0], m.inertia[1][1], m.inertia[1][2],
                m.inertia[2][0], m.inertia[2][1], m.inertia[2][2]))
        }
        "isbbinterf" if args.len() == 3 => {
            let a = get(shapes, &args[1])?.bounds();
            let b = get(shapes, &args[2])?.bounds();
            Ok(if a.intersects(b, t)? {
                "The shapes are interfered by AABB.\n"
            } else {
                "The shapes are NOT interfered by AABB.\n"
            }
            .into())
        }
        _ => Err(unsupported(args)),
    }
}

// Hex-encoded UTF-8 tokens keep Tcl names, spaces and newlines out of framing.
fn unhex(s: &str) -> Result<String> {
    let bytes = s.as_bytes();
    if bytes.len() % 2 != 0 || !bytes.is_ascii() {
        return Err(Failure::Error("malformed wire token".into()));
    }
    let decoded = (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|_| Failure::Error("malformed wire hex".into()))?;
    String::from_utf8(decoded).map_err(|_| Failure::Error("malformed wire UTF-8".into()))
}
fn hex(s: &str) -> String {
    s.as_bytes().iter().map(|b| format!("{b:02x}")).collect()
}
fn main() -> io::Result<()> {
    let mut shapes = BTreeMap::new();
    let mut output = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let result = line?
            .split(' ')
            .map(unhex)
            .collect::<Result<Vec<_>>>()
            .and_then(|args| dispatch(&mut shapes, &args));
        let (status, message) = match result {
            Ok(value) => ("OK", value),
            Err(Failure::Error(value)) => ("ERROR", value),
            Err(Failure::Unsupported(value)) => ("UNSUPPORTED", value),
        };
        writeln!(output, "{status} {}", hex(&message))?;
        output.flush()?;
    }
    Ok(())
}
