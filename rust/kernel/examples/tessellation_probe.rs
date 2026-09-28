//! Test-only input for compare_tessellation.py (T-a of REVIEW_NOTES.md).
//!
//! Reads `tessellation-cases.txt` on stdin: identity-protocol case blocks
//! with `mesh SETTING DEFLECTION ANGLE` rows. Builds each body and prints
//! each setting's mesh as `mesh CASE SETTING deflection angle` (the mesh's
//! certified maxima), `v x y z` per node, `t a b c face deflection angle`
//! per triangle, `e edge deflection angle nodes...` per edge polyline (its
//! largest segment bounds), `end`; or `error CASE SETTING message`.
#[path = "../tests/support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use rusty_occt::tessellation::{Mesh, Parameters};
use std::fmt::Write;
use std::io::Read;

fn print(name: &str, setting: &str, mesh: &Mesh) {
    let mut out = String::new();
    writeln!(
        out,
        "mesh {name} {setting} {:?} {:?}",
        mesh.deflection, mesh.angle
    )
    .unwrap();
    for p in &mesh.nodes {
        writeln!(out, "v {:?} {:?} {:?}", p.x, p.y, p.z).unwrap();
    }
    for f in &mesh.faces {
        for k in f.triangles.clone() {
            let [a, b, c] = mesh.triangles[k];
            let bound = mesh.triangle_bounds[k];
            writeln!(
                out,
                "t {a} {b} {c} {} {:?} {:?}",
                f.face.index(),
                bound.deflection,
                bound.angle
            )
            .unwrap();
        }
    }
    for e in &mesh.edges {
        let d = e.segments.iter().map(|s| s.deflection).fold(0.0, f64::max);
        let a = e.segments.iter().map(|s| s.angle).fold(0.0, f64::max);
        let nodes: Vec<String> = e.nodes.iter().map(|n| n.to_string()).collect();
        writeln!(out, "e {} {d:?} {a:?} {}", e.edge.index(), nodes.join(" ")).unwrap();
    }
    out += "end\n";
    print!("{out}");
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for block in input.split("\nend").filter(|b| !b.trim().is_empty()) {
        let rows: Vec<&str> = block.trim().lines().collect();
        let settings: Vec<Vec<&str>> = rows
            .iter()
            .filter(|r| r.starts_with("mesh "))
            .map(|r| r.split_whitespace().collect())
            .collect();
        let body: String = rows
            .iter()
            .filter(|r| !r.starts_with("mesh "))
            .map(|r| format!("{r}\n"))
            .collect();
        let spec = identity_protocol::parse(&body);
        for s in settings {
            let (deflection, angle): (f64, f64) = (s[2].parse().unwrap(), s[3].parse().unwrap());
            let result = Parameters::new(deflection, angle).and_then(|p| match spec.make {
                Some(_) => identity_protocol::build_body(&spec).0.tessellate(p),
                None => identity_protocol::build(&spec).tessellate(p),
            });
            match result {
                Ok(mesh) => print(&spec.name, s[1], &mesh),
                Err(e) => println!("error {} {} {e}", spec.name, s[1]),
            }
        }
    }
}
