//! Test-only protocol for compare_revolve_history.py. Reads
//! `primitive-cases.txt` lines (`cone NAME ox oy oz nx ny nz xx xy xz r1 r2
//! h`), builds each with `Solid::cone_with` and prints the synthesized OCCT
//! counts (`C`), every construction relation with its meridian parent, role
//! and a geometric signature of its target (`G`) in the format of
//! rust/tools/occt_revolve_oracle.cpp, and the body (`B`).
#[path = "../tests/support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
#[path = "../tests/support/signatures.rs"]
mod signatures;
use identity_protocol::{parent_text, role_name};
use rusty_occt::history::Relation;
use rusty_occt::identity::OperationId;
use rusty_occt::{Frame3, Point3, Solid, Tolerance, Vec3};
use signatures::{body, counts, signature};
use std::io::Read;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for line in input.lines().filter(|l| !l.trim().is_empty()) {
        let w: Vec<&str> = line.split_whitespace().collect();
        let n: Vec<f64> = w[2..].iter().map(|x| x.parse().unwrap()).collect();
        let tol = Tolerance::default();
        let frame = Frame3::new(
            Point3::new(n[0], n[1], n[2]),
            Vec3::new(n[3], n[4], n[5]),
            Vec3::new(n[6], n[7], n[8]),
            tol,
        )
        .unwrap();
        let (solid, history) =
            Solid::cone_with(OperationId(1), frame, n[9], n[10], n[11], tol).unwrap();
        println!("case {}", w[1]);
        println!("{}", counts(&solid));
        for r in &history.relations {
            let Relation::Generated { from, to, role } = r else {
                panic!("a construction only generates");
            };
            let slot = solid.topology().slot_of(*to).expect("a live target");
            let parents: Vec<String> = from.iter().map(parent_text).collect();
            println!(
                "G {} {} {}",
                parents.join(","),
                role_name(*role),
                signature(&solid, slot)
            );
        }
        println!("{}", body(&solid));
        println!("end");
    }
}
