//! The fix of S9e.4b.4d.2c.1a's consistency (REVIEW_NOTES.md, "S9e.4b.4d.2c.1a's
//! consistency fixed"): an imported turned body of `OCC485a`'s form (a
//! bottom and a step within the resolution of planes, a fillet and a round
//! off their snapped tangencies, an axis tilted by 2^-27 rad, resolution
//! 3e-3; `fixtures/imported/consistency_ribbed.brep`, written by OCCT from
//! `tools/generate_turned_consistency_fixtures.py`) is stored as its model,
//! so that its Booleans with a rib through its step, fillet, shaft and round
//! and with boxes across its cylinders, its fillet and its round evaluate,
//! agree with one another within their enclosures (the fuse and the common
//! in both orders, each cut and the common summing to its input, the fuse
//! to both inputs less the common), keep `y = 0` where both inputs are
//! symmetric in `y`, and give the model's measures
//! (`fixtures/boolean-turned-consistency-expected.tsv`). The other turned
//! fixtures are checked the same way against a box across each.
#[path = "support/turned_consistency.rs"]
#[allow(dead_code)]
mod consistency;
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use consistency::{consistent, imported};
use rusty_occt::identity::OperationId;

use rusty_occt::{Point3, Solid, Tolerance, Vec3};

/// The body is its model (the reference's volume by Pappus within 1e-9),
/// and every partner's six Booleans with it evaluate, agree with one another
/// and keep `y = 0` where the partner is symmetric in `y`, the common the
/// model's.
#[test]
fn the_ribbed_body_s_booleans_are_consistent() {
    let failures = consistency::check_partners(include_str!(
        "../../fixtures/boolean-turned-consistency-expected.tsv"
    ));
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The other imported turned fixtures, each against a box across it (from
/// inside past its side, across its middle heights): their six Booleans
/// evaluate and agree with one another.
#[test]
fn other_turned_bodies_booleans_are_consistent() {
    let tol = Tolerance::default();
    let mut failures = Vec::new();
    for name in [
        "turned_fillet",
        "turned_pin",
        "turned_step",
        "fillet_step",
        "fillet_collar",
        "fillet_refit",
        "flat_pin",
        "flat_collar",
        "profile_crown",
        "profile_point",
        "profile_steps",
        "profile_waist",
        "tilted_waist",
        "tilted_pin",
        "cone_fillet_fillet",
        "cone_fillet_refit",
        "cone_fillet_round",
        "cone_fillet_top",
    ] {
        let body = match imported(&format!("imported/{name}.brep")) {
            Ok(b) => b,
            Err(e) => {
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let b = body.bounds();
        let (lo, hi) = (b.min.to_array(), b.max.to_array());
        let at = |k: usize, f: f64| lo[k] + f * (hi[k] - lo[k]);
        let start = Point3::new(at(0, 0.5625), at(1, 0.3125), at(2, 0.296875));
        let size = Vec3::new(
            hi[0] - lo[0],
            (hi[1] - lo[1]) * 0.40625,
            (hi[2] - lo[2]) * 0.46875,
        );
        let partner = Solid::box_at_with(OperationId(92), start, size, tol)
            .unwrap()
            .0;
        if let Err(e) = consistent(&body, &partner, false) {
            failures.push(format!("{name}: {e}"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
