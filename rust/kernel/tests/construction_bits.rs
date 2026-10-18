//! R13 (the review of 18 October 2026): every binary64 the kernel
//! constructs is the same on every host. Each family of constructions that
//! called the platform's libm before is pinned here by a digest of every
//! number it builds (`{:?}` prints each binary64's shortest round trip, so
//! a digest moves with any bit): frames and lengths (`Frame3::new`,
//! `Vec3::length`), rotations, the identity cases' solids (prisms with
//! arcs, cones, spheres, tori, their transformed copies and mass
//! properties), every imported fixture body (its recognized model: turned
//! bodies' slopes, pieces' bounds), the STEP fixtures' imports, the `.brep`
//! texts written back and meshes. `fixtures/construction-bits.tsv` was
//! written on macOS (`RUSTY_WRITE_CONSTRUCTION_BITS=1`) and is reproduced
//! on a Linux host (glibc 2.39) and by every CI job.
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use rusty_occt::identity::OperationId;
use rusty_occt::tessellation::{tessellate, Parameters};
use rusty_occt::{occt_brep, step, Frame3, Point3, RigidTransform, Solid, Tolerance, Vec3};
use std::path::PathBuf;

/// FNV-1a of a text: stable across Rust releases, unlike `DefaultHasher`.
fn fnv(text: &str) -> String {
    let h = text.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    });
    format!("{h:016x}")
}

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../fixtures")
}

/// Deterministic directions: a linear congruential generator's bits.
fn directions(n: usize) -> Vec<Vec3> {
    let mut s: u64 = 0x0c0f_fee5;
    let mut next = || {
        s = s
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (s >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    let mut out = vec![
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(0.0, 1.0, 15f64.sqrt()),
        Vec3::new(4.0, 1.0, 8.0),
        Vec3::new(-7.0, -4.0, 4.0),
    ];
    while out.len() < n {
        out.push(Vec3::new(next(), next(), next()));
    }
    out
}

fn rows() -> Vec<String> {
    let mut out = Vec::new();
    let tol = Tolerance::new(1e-7, 1e-12).unwrap();
    let dirs = directions(48);
    for (i, w) in dirs.iter().enumerate() {
        let hint = dirs[(i + 1) % dirs.len()];
        let frame = Frame3::new(Point3::new(1.5, -2.25, 0.75), *w, hint, tol);
        let text = format!("{:?} {:?} {:?}", w.length(), w.normalized(), frame);
        out.push(format!("frame\t{i}\t{}", fnv(&text)));
    }
    for (i, w) in dirs.iter().take(24).enumerate() {
        let angle = (i as f64 - 11.5) * 0.37;
        let r = RigidTransform::rotation(Point3::new(0.5, -1.0, 2.0), *w, angle);
        out.push(format!("rotation\t{i}\t{}", fnv(&format!("{r:?}"))));
    }
    for spec in identity_protocol::cases(include_str!("../../fixtures/identity-cases.txt")) {
        let text = if spec.make.is_some() {
            format!("{:?}", identity_protocol::build_body(&spec))
        } else {
            let mut solid = identity_protocol::build(&spec);
            let mut text = format!("{solid:?}");
            for t in &spec.transforms {
                solid = solid.transform_with(OperationId(1), *t).unwrap().0;
                text += &format!("{solid:?}");
            }
            text
        };
        out.push(format!("identity\t{}\t{}", spec.name, fnv(&text)));
    }
    let mut breps: Vec<_> = std::fs::read_dir(fixtures().join("imported"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "brep"))
        .collect();
    breps.sort();
    for path in breps {
        let name = path.file_stem().unwrap().to_string_lossy().to_string();
        let doc = occt_brep::read(&std::fs::read_to_string(&path).unwrap()).unwrap();
        let imported = occt_brep::import(&doc);
        let mut text = String::new();
        for s in imported.solids {
            match s.result {
                Ok(t) => {
                    text += &occt_brep::write(&t, s.tolerance.linear())
                        .unwrap_or_else(|e| format!("{e:?}"));
                    let solid = Solid::imported_with(OperationId(91), t, s.tolerance);
                    text += &format!("{solid:?}");
                }
                Err(e) => text += &format!("{e:?}"),
            }
        }
        out.push(format!("imported\t{name}\t{}", fnv(&text)));
    }
    let mut steps: Vec<_> = std::fs::read_dir(fixtures().join("step"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "stp"))
        .collect();
    steps.sort();
    for path in steps {
        let name = path.file_stem().unwrap().to_string_lossy().to_string();
        let bytes = std::fs::read(&path).unwrap();
        let text = match step::read(&bytes) {
            Ok(file) => format!("{:?}", step::import(&file)),
            Err(e) => format!("{e:?}"),
        };
        out.push(format!("step\t{name}\t{}", fnv(&text)));
    }
    let parameters = Parameters::new(0.05, 0.3).unwrap();
    for spec in identity_protocol::cases(include_str!("../../fixtures/identity-cases.txt"))
        .into_iter()
        .filter(|s| {
            s.make.is_none() && (s.cone.is_some() || s.sphere.is_some() || s.torus.is_some())
        })
    {
        let solid = identity_protocol::build(&spec);
        let text = format!("{:?}", tessellate(solid.topology(), parameters));
        out.push(format!("mesh\t{}\t{}", spec.name, fnv(&text)));
    }
    out
}

#[test]
fn every_construction_is_bitwise_the_fixtures() {
    let got = rows();
    if std::env::var_os("RUSTY_WRITE_CONSTRUCTION_BITS").is_some() {
        let path = fixtures().join("construction-bits.tsv");
        let header = "# Digests of every binary64 the kernel constructs, by family: \
                      tests/construction_bits.rs\n";
        std::fs::write(path, header.to_string() + &got.join("\n") + "\n").unwrap();
    }
    let want: Vec<&str> = include_str!("../../fixtures/construction-bits.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    let differ: Vec<String> = got
        .iter()
        .zip(&want)
        .filter(|(g, w)| g != w)
        .map(|(g, w)| format!("{g} for {w}"))
        .collect();
    assert_eq!(got.len(), want.len(), "the fixture's rows");
    assert!(differ.is_empty(), "{differ:#?}");
}
