//! Headless analytic solid geometry with explicit numerical contracts.
//!
//! This first implementation supports normal extrusions of simple polygons and
//! circles, including disjoint polygonal/circular holes. It retains analytic
//! curves, surfaces, oriented topology and face parameter curves. It does not
//! yet implement general Boolean operations, blends, STEP or tessellation.
//!
//! All distances are in millimetres, angles in radians. Geometric constructors
//! validate their inputs and return errors; there is no C++ runtime dependency.
//!
//! ```
//! use rusty_occt::{Boundary, Frame3, Point2, Profile, Solid, Tolerance};
//! let tolerance = Tolerance::default();
//! let outer = Boundary::rectangle(40.0, 20.0, tolerance)?;
//! let hole = Boundary::circle(Point2::new(10.0, 10.0), 3.0, tolerance)?;
//! let profile = Profile::new(outer, vec![hole], tolerance)?;
//! let plate = Solid::extrude(profile, Frame3::xy(), 0.0, 5.0)?;
//! plate.topology().validate(tolerance)?;
//! assert_eq!(plate.topology().faces().len(), 7);
//! # Ok::<(), rusty_occt::Error>(())
//! ```

mod algebraic;
pub mod curve;
mod error;
mod exact;
pub mod intersection;
mod interval;
pub mod math;
pub mod polynomial;
pub mod predicates;
mod profile;
pub mod proximity;
mod solid;
pub mod spline;
pub mod surface;
pub mod topology;

pub use curve::{BSplineCurve3, BezierCurve3};
pub use error::{Error, Result};
pub use intersection::{Circle3, Cylinder3, Plane3, Sphere3, Triangle3};
pub use interval::ScalarInterval;
pub use math::{Bounds3, Frame3, Point2, Point3, RigidTransform, Tolerance, Vec3};
pub use profile::{Boundary, Location, Profile};
pub use solid::{MassProperties, Solid};
pub use spline::KnotVector;
pub use surface::{BSplineSurface3, BezierSurface3};
