//! S9e.4b.4d.1 (REVIEW_NOTES.md, "S9e.4b.4d refined"): turned bodies OCCT
//! wrote (`tools/write_turned_bodies.py`, `fixtures/imported/turned_*.brep`)
//! refused by name (S9e.4b.4d.2c.1a: but a fillet's and a round's body
//! snapped within its resolution), a cylinder OCCT stores as two faces on
//! frames apart no longer the body's tangency, and every other imported
//! fixture's outcome as before the step (`eefb2294`).

#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;

use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Solid};

const SMOOTH: &str =
    "an imported body whose torus faces join their neighbours smoothly (S9e.4b.4d)";
const NEAR_FLAT: &str =
    "an imported cone within the resolution of a plane or a cylinder (S9e.4b.4d)";
const GENERAL: &str = "an imported solid other than a prism, a sphere, a cone, a torus or a \
                       Boolean of its primitives (S9e.4b)";
const TANGENT: &str =
    "an imported body of several primitives whose faces are tangent along an edge";
const NO_CHAIN: &str =
    "an imported body of several primitives other than a Boolean chain of them (S9e.4b.4c)";

fn body(name: &str) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(91), topology, resolution).map(|(s, _)| s)
}

fn refusal(name: &str) -> Error {
    match body(name) {
        Ok(_) => panic!("{name}: imported"),
        Err(e) => e,
    }
}

/// A turned part's fillet and round (tori tangent to the step's plane, the
/// shaft and the top's plane) were refused by name; S9e.4b.4d.2c.1a: such a
/// body imports, its tori snapped tangent to their lines (its volume its
/// profile's by Pappus: its sharp profile's `pi (400 + 150)` with its
/// fillet's corner region and less its round's), while one whose round's
/// stored edge with the shaft lies 1e-4 off the snapped tangency (its tube
/// crossing the shaft by 2.5e-9, the faces' normals at that edge 5e-5 rad
/// apart: a thousand times its resolution of 1e-7) stays refused by name; a
/// torus meeting cones at corners keeps the general refusal.
#[test]
fn turned_fillets_and_rounds_are_refused_by_name() {
    let fillet = body("turned_fillet").unwrap_or_else(|e| panic!("turned_fillet: {e}"));
    let pi = std::f64::consts::PI;
    // A corner region's `int r dA`: the square of side `t` from the centre
    // `cr` toward the corner (`sr` its side) less the quarter disc.
    let corner = |cr: f64, t: f64, sr: f64| {
        t * t * (cr + sr * t / 2.0) - (cr * pi * t * t / 4.0 + sr * t * t * t / 3.0)
    };
    let want = pi * 550.0 + 2.0 * pi * (corner(5.5, 0.5, -1.0) - corner(3.0, 2.0, 1.0));
    let v = fillet.mass_properties().volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} for {want}");
    assert_eq!(refusal("turned_round_off"), Error::OutOfDomain(SMOOTH));
    assert_eq!(refusal("turned_corner"), Error::OutOfDomain(GENERAL));
}

/// A pin whose step is a cone rising 1e-5 over its face (5e-6 rad from
/// flat) in a body of resolution 1e-3 was refused as within the resolution
/// of a plane (as `OCC497c`'s steps, refused before as two cylinders'
/// section near a node); S9e.4b.4d.2b: a turned body, it imports with the
/// step stored as the plane through its wider rim (its volume the flat
/// profile's), and a rod whose top cone rises 2^-17 to its apex (one curved
/// face beside it, no lattice) stays refused so. The same step rising 1e-3,
/// ten thousand resolutions of 1e-7, imports as before, its volume the
/// profile's.
#[test]
fn a_cone_within_the_resolution_of_a_plane_is_degenerate() {
    assert_eq!(refusal("flat_rod"), Error::Degenerate(NEAR_FLAT));
    let pin = body("turned_pin").unwrap_or_else(|e| panic!("turned_pin: {e}"));
    let flat = std::f64::consts::PI * (180.0 + 16.0 * 7.0);
    let v = pin.mass_properties().volume;
    assert!((v - flat).abs() <= 1e-9 * flat, "{v} for {flat}");
    let s = body("turned_step").unwrap_or_else(|e| panic!("turned_step: {e}"));
    // Radius 6 to height 5, a frustum from radius 6 to 4 over 1e-3, radius 4
    // to height 12.
    let volume = std::f64::consts::PI * (180.0 + 0.001 * 76.0 / 3.0 + 16.0 * 6.999);
    let v = s.mass_properties().volume;
    assert!((v - volume).abs() <= 1e-9 * volume, "{v} for {volume}");
}

/// A cylinder OCCT stores as two faces on frames 3 apart along its axis
/// (`bug25966_Revolution`'s configuration: a turned profile with a vertex on
/// its wall's line) is no longer the body's tangency: the two faces lie on
/// one surface across a whole ring. The body is a turned profile no chain of
/// its primitives builds (S9e.4b.4c.2b.2).
#[test]
fn faces_on_one_cylinder_across_a_ring_are_no_tangency() {
    let e = refusal("turned_split");
    assert_ne!(e, Error::Degenerate(TANGENT));
    assert_eq!(e, Error::OutOfDomain(NO_CHAIN));
}

/// Every other body OCCT wrote for the import's steps: its outcome as
/// before S9e.4b.4d.1 (`None` imports; else the refusal's kind and text).
/// `leaf_post`'s fillet face and post on one cylinder across an arc stay the
/// body's tangency.
#[test]
fn every_other_imported_fixture_is_as_before() {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(8);
    let chunks: Vec<_> = TABLE.chunks(TABLE.len().div_ceil(threads)).collect();
    std::thread::scope(|scope| {
        for chunk in chunks {
            scope.spawn(move || {
                for &(name, want) in chunk {
                    let got = match body(name) {
                        Ok(_) => None,
                        Err(Error::Degenerate(m)) => Some(("Degenerate", m)),
                        Err(Error::OutOfDomain(m)) => Some(("OutOfDomain", m)),
                        Err(e) => panic!("{name}: {e}"),
                    };
                    assert_eq!(got, want, "{name}");
                }
            });
        }
    });
}

type Outcome = Option<(&'static str, &'static str)>;

const TABLE: [(&str, Outcome); 113] = [
    ("arch", None),
    ("ball", None),
    ("bitten", None),
    ("blade", None),
    ("box", None),
    ("box_tilt", None),
    ("cam", None),
    ("chain_bead", None),
    ("chain_capsule", Some(("Degenerate", "an imported body of several primitives whose faces are tangent along an edge"))),
    ("chain_cup", None),
    ("chain_dome", None),
    ("chain_knob", None),
    ("chain_pin", None),
    ("chain_rounded", None),
    ("chain_shaft", None),
    ("cyl", None),
    ("cyl_side", None),
    ("dee", None),
    ("dee_turn", None),
    ("deep_post", None),
    ("deep_well", None),
    ("dome", None),
    ("draft", None),
    ("ell", None),
    ("flush_arms", None),
    ("flush_open", None),
    ("flush_slot", None),
    ("form_ball_boss", None),
    ("form_bite", None),
    ("form_bites", None),
    ("form_boss", None),
    ("form_cap_boss", None),
    ("form_cap_pocket", None),
    ("form_dee_boss", None),
    ("form_dimple", None),
    ("form_notch", Some(("Degenerate", "an imported plane piece whose curved face is tangent to its plane faces"))),
    ("form_quarter", Some(("Degenerate", "a plane through a cone's apex"))),
    ("form_scoop", None),
    ("form_sink", None),
    ("form_slot", None),
    ("form_tooth", None),
    ("form_u_boss", None),
    ("form_u_notch", Some(("Degenerate", "an imported plane piece whose curved face is tangent to its plane faces"))),
    ("form_u_scoop", None),
    ("frustum", None),
    ("half", None),
    ("halves", None),
    ("halves_turn", None),
    ("hollow", None),
    ("incidence_back", None),
    ("incidence_half", None),
    ("incidence_hemi", None),
    ("incidence_high23", None),
    ("incidence_lifted45", None),
    ("incidence_octant", None),
    ("incidence_turned45", None),
    ("incidence_wedge23", None),
    ("incidence_wedge45", None),
    ("incidence_x23", None),
    ("incidence_x45", None),
    ("leaf_bore", None),
    ("leaf_boss", None),
    ("leaf_dimple", None),
    ("leaf_dome", None),
    ("leaf_notch", Some(("OutOfDomain", "an imported body of several primitives with plane faces other than their ends, one or two prisms', a flat or a pocket, its pockets nested at most three deep (S9e.4b.4b.2b.2b.3b)"))),
    ("leaf_pocket", None),
    ("leaf_post", Some(("Degenerate", "an imported body of several primitives whose faces are tangent along an edge"))),
    ("lens", None),
    ("lune", None),
    ("nested_island", None),
    ("nested_post", None),
    ("notch", None),
    ("notched", None),
    ("octa", None),
    ("octant", None),
    ("octant_low", None),
    ("octant_tilt", None),
    ("part_cake", None),
    ("part_flat", None),
    ("part_slot", None),
    ("part_stack", None),
    ("pedestal", None),
    ("plate", None),
    ("poly_cavity", None),
    ("pyramid", None),
    ("quad", None),
    ("ridge", None),
    ("rim_post", None),
    ("rim_well", None),
    ("ring", None),
    ("rounded", None),
    ("slot", None),
    ("sphere_cap", None),
    ("sphere_hemi", None),
    ("sphere_octant2", None),
    ("sphere_octant4", None),
    ("sphere_wedge", None),
    ("spike", None),
    ("split_lens", None),
    ("stack_bar", None),
    ("stack_step", None),
    ("steps_high", None),
    ("steps_low", None),
    ("tetra", None),
    ("truncated", None),
    ("two_cross", None),
    ("two_pocketed", None),
    ("two_through", None),
    ("upper", None),
    ("vane_down", None),
    ("vane_up", None),
    ("wedge", None),
    ("zone_wedge", None),
];
