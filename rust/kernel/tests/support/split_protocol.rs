//! The case protocol of `split-cases.txt` (S8) and the kernel's rows for it,
//! shared by `split_probe` and `tests/split.rs`: `limit`, `unsupported`, or
//! per piece `piece below|above vol_lo vol_hi area_lo area_hi cx_lo cx_hi
//! cy_lo cy_hi cz_lo cz_hi faces edges vertices` (OCCT's counts of the piece,
//! `Topology::occt_counts`).
#[path = "identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
pub use identity_protocol::{build, parse, CaseSpec};
use rusty_occt::history::History;
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Frame3, Point3, Side, Solid, Tolerance, Vec3};

pub struct Case {
    pub spec: CaseSpec,
    /// The plane: a point and its normal.
    pub plane: [f64; 6],
}

pub fn cases(text: &str) -> Vec<Case> {
    text.split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|block| {
            let (mut rest, mut plane) = (Vec::new(), [0.0; 6]);
            for line in block.trim().lines() {
                match line.strip_prefix("split ") {
                    Some(p) => {
                        let v: Vec<f64> =
                            p.split_whitespace().map(|w| w.parse().unwrap()).collect();
                        plane.copy_from_slice(&v);
                    }
                    None => rest.push(line),
                }
            }
            Case {
                spec: parse(&rest.join("\n")),
                plane,
            }
        })
        .collect()
}

pub fn plane_frame(p: [f64; 6]) -> Frame3 {
    let n = Vec3::new(p[3], p[4], p[5]);
    let hint = if n.x.abs() < 0.5 * n.length() {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    Frame3::new(Point3::new(p[0], p[1], p[2]), n, hint, Tolerance::default()).unwrap()
}

/// A case's solid, its pieces with their sides, and the history.
pub type Split = (Solid, Vec<(Side, Solid)>, History);

/// The pieces and history of a case.
pub fn split(case: &Case) -> Result<Split, Error> {
    let solid = build(&case.spec);
    let (pieces, history) = solid.split_by_plane(OperationId(900), plane_frame(case.plane))?;
    Ok((solid, pieces, history))
}

/// The kernel's rows for a case; `Err` for an unexpected error.
pub fn rows(case: &Case) -> Result<Vec<String>, Error> {
    let pieces = match split(case) {
        Ok((_, pieces, _)) => pieces,
        Err(Error::ComputationLimit(_)) => return Ok(vec!["limit".into()]),
        Err(Error::OutOfDomain(_)) => return Ok(vec!["unsupported".into()]),
        Err(e) => return Err(e),
    };
    let mut out = Vec::new();
    for (side, piece) in pieces {
        let t = piece.topology();
        let m = t
            .mass_enclosure()
            .ok_or(Error::ComputationLimit("a piece's mass"))?;
        let c = t.occt_counts();
        let side = match side {
            Side::Below => "below",
            Side::Above => "above",
        };
        let mut words = vec![format!("piece {side}")];
        words.push(format!("{:?} {:?}", m.volume[0], m.volume[1]));
        words.push(format!("{:?} {:?}", m.surface_area[0], m.surface_area[1]));
        for [lo, hi] in m.centroid {
            words.push(format!("{lo:?} {hi:?}"));
        }
        words.push(format!("{} {} {}", c.faces, c.edges, c.vertices));
        out.push(words.join(" "));
    }
    Ok(out)
}
