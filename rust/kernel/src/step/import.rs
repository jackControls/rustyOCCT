//! STEP bodies into the cell model (STEP-a and STEP-b of the STEP import
//! track in `REVIEW_NOTES.md`). Each solid or surface-model shell becomes
//! OCCT's shape structure, an `occt_brep::Document` built as `StepToTopoDS`
//! builds its `TopoDS` shapes, and the `.brep` converter turns that into
//! cells: seams merge into windings, degenerated edges into poles, and the
//! result passes `Topology::from_parts` or is rejected with every validation
//! issue.
//!
//! What STEP does not carry is supplied here: lengths in millimetres, the
//! file's uncertainty as every entity's tolerance, and the pcurves. On a
//! plane the converter derives them from the edges (`CurveOnPlane`); on a
//! cylinder, cone, sphere or torus each use of a ruling, a parallel or a
//! meridian is a straight segment in `(u, v)`, placed to continue the last
//! in the universal cover, and a loop's passage through a pole gets the
//! degenerated edge OCCT's reader adds, running the way that keeps the face
//! on the loop's left. An ellipse that is a cylinder's plane section gets the
//! exact sinusoid over the same cover (STEP-b). On a B-spline surface the
//! file's own pcurve is taken (`spline.rs`), over the range where its image
//! meets the edge's vertices. The validator certifies every such pcurve
//! against its edge.
use super::part21::{Exchange, Instance, Parameter};
use super::spline;
use super::StepError;
use crate::occt_brep::read::{self as brep, Data, Document, EdgeRep, Kind, Orient, Shape, Sub};
use crate::occt_brep::{self, Rejected};
use crate::topology::Topology;
use crate::Tolerance;
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// Which STEP item a body comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Item {
    /// A `MANIFOLD_SOLID_BREP` or `BREP_WITH_VOIDS`.
    Solid,
    /// One shell of a `SHELL_BASED_SURFACE_MODEL`: a sheet or a closed shell
    /// bounding a void (S6).
    Shell,
}

/// One body of a file.
#[derive(Debug, Clone, PartialEq)]
pub struct StepBody {
    /// The entity number of the solid or surface model.
    pub entity: u64,
    pub item: Item,
    /// The body resolution: the context's uncertainty in millimetres,
    /// floored at `1e-7` and capped at `1`.
    pub tolerance: Tolerance,
    /// The cells, or why there are none: unsupported constructs by name, or
    /// the validation issues of the converted cells.
    pub result: Result<Topology, Rejected>,
}

/// Everything a file holds, as far as the kernel represents it.
#[derive(Debug, Clone, PartialEq)]
pub struct StepImport {
    /// Solids, then surface-model shells in order, by entity number.
    pub bodies: Vec<StepBody>,
    /// Every construct the kernel cannot represent, by name, with the
    /// number of bodies (or files, for `AssemblyPlacement`) that carry it.
    pub unsupported: BTreeMap<&'static str, usize>,
}

/// The schemas of AP203 (both editions), AP214 and AP242.
const SCHEMAS: [&str; 4] = [
    "CONFIG_CONTROL_DESIGN",
    "AP203_CONFIGURATION_CONTROLLED_3D_DESIGN_OF_MECHANICAL_PARTS_AND_ASSEMBLIES_MIM_LF",
    "AUTOMOTIVE_DESIGN",
    "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF",
];

/// Entity names reported when a body reaches them, most specific first
/// (complex instances list several).
const NAMES: [&str; 45] = [
    "RATIONAL_B_SPLINE_SURFACE",
    "RATIONAL_B_SPLINE_CURVE",
    "B_SPLINE_SURFACE_WITH_KNOTS",
    "B_SPLINE_CURVE_WITH_KNOTS",
    "BEZIER_SURFACE",
    "BEZIER_CURVE",
    "UNIFORM_SURFACE",
    "UNIFORM_CURVE",
    "QUASI_UNIFORM_SURFACE",
    "QUASI_UNIFORM_CURVE",
    "B_SPLINE_SURFACE",
    "B_SPLINE_CURVE",
    "ELLIPSE",
    "HYPERBOLA",
    "PARABOLA",
    "POLYLINE",
    "OFFSET_CURVE_3D",
    "COMPOSITE_CURVE",
    "INTERSECTION_CURVE",
    "PCURVE",
    "DEGENERATE_PCURVE",
    "SURFACE_OF_LINEAR_EXTRUSION",
    "SURFACE_OF_REVOLUTION",
    "OFFSET_SURFACE",
    "RECTANGULAR_TRIMMED_SURFACE",
    "CURVE_BOUNDED_SURFACE",
    "DEGENERATE_TOROIDAL_SURFACE",
    "POLY_LOOP",
    "VERTEX_LOOP",
    "SUBEDGE",
    "SUBFACE",
    "ORIENTED_EDGE",
    "ORIENTED_FACE",
    "FACETED_BREP",
    "AXIS2_PLACEMENT_2D",
    "CARTESIAN_POINT",
    "DIRECTION",
    "VECTOR",
    "LINE",
    "CIRCLE",
    "PLANE",
    "EDGE_CURVE",
    "VERTEX_POINT",
    "ADVANCED_FACE",
    "CLOSED_SHELL",
];

/// A construct the kernel cannot (yet) represent, or malformed data, by
/// name.
pub(super) type Named<T> = Result<T, &'static str>;

/// The name reported for an instance of an unexpected type.
fn name_of(i: &Instance) -> &'static str {
    NAMES
        .iter()
        .find(|n| i.record(n).is_some())
        .copied()
        .unwrap_or("UnknownEntity")
}

pub(super) fn real(p: &Parameter) -> Named<f64> {
    match p {
        Parameter::Real(x) => Ok(*x),
        Parameter::Integer(i) => Ok(*i as f64),
        Parameter::Typed(_, v) => real(v),
        _ => Err("MalformedEntity"),
    }
}

pub(super) fn reference(p: &Parameter) -> Named<u64> {
    match p {
        Parameter::Reference(n) => Ok(*n),
        _ => Err("MalformedEntity"),
    }
}

fn boolean(p: &Parameter) -> Named<bool> {
    match p {
        Parameter::Enumeration(e) if e == "T" => Ok(true),
        Parameter::Enumeration(e) if e == "F" => Ok(false),
        _ => Err("MalformedEntity"),
    }
}

pub(super) fn list(p: &Parameter) -> Named<&[Parameter]> {
    match p {
        Parameter::List(items) => Ok(items),
        _ => Err("MalformedEntity"),
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn along(p: [f64; 3], d: [f64; 3], t: f64) -> [f64; 3] {
    [p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: [f64; 3]) -> Named<[f64; 3]> {
    let n = norm(a);
    if n > 0.0 && n.is_finite() {
        Ok(a.map(|c| c / n))
    } else {
        Err("DegenerateDirection")
    }
}

/// A right-handed orthonormal placement.
#[derive(Debug, Clone, Copy)]
struct Axes {
    o: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    z: [f64; 3],
}

impl Axes {
    fn local(&self, p: [f64; 3]) -> [f64; 3] {
        let q = sub(p, self.o);
        [dot(q, self.x), dot(q, self.y), dot(q, self.z)]
    }
    fn angle(&self, p: [f64; 3]) -> f64 {
        let [a, b, _] = self.local(p);
        b.atan2(a)
    }
    fn circle(&self, r: f64, t: f64) -> [f64; 3] {
        let (s, c) = t.sin_cos();
        std::array::from_fn(|i| self.o[i] + r * (c * self.x[i] + s * self.y[i]))
    }
}

/// The units of a representation context.
#[derive(Debug, Clone, Copy)]
struct Units {
    /// Millimetres per length unit.
    length: f64,
    /// Radians per plane angle unit.
    angle: f64,
    /// Every entity's tolerance, in millimetres.
    tolerance: f64,
}

const DEFAULT_UNITS: Units = Units {
    length: 1.0,
    angle: 1.0,
    tolerance: 1e-7,
};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Measure {
    Length,
    Angle,
}

/// The factor of a unit to millimetres or radians.
fn unit_factor(x: &Exchange, n: u64, measure: Measure, depth: usize) -> Named<f64> {
    let i = x.instances.get(&n).ok_or("UnitsUndetermined")?;
    if depth > 8 {
        return Err("UnitsUndetermined");
    }
    if let Some(si) = i.record("SI_UNIT") {
        let [prefix, name] = si.parameters.as_slice() else {
            return Err("UnitsUndetermined");
        };
        let prefix = match prefix {
            Parameter::Enumeration(p) => Some(p.as_str()),
            _ => None,
        };
        let base = match name {
            Parameter::Enumeration(b) => b.as_str(),
            _ => "",
        };
        return match (measure, base, prefix) {
            (Measure::Length, "METRE", None) => Ok(1000.0),
            (Measure::Length, "METRE", Some("MILLI")) => Ok(1.0),
            (Measure::Length, "METRE", Some("CENTI")) => Ok(10.0),
            (Measure::Length, "METRE", Some("DECI")) => Ok(100.0),
            (Measure::Length, "METRE", Some("KILO")) => Ok(1e6),
            (Measure::Length, "METRE", Some("MICRO")) => Ok(1e-3),
            (Measure::Length, "METRE", Some("NANO")) => Ok(1e-6),
            (Measure::Angle, "RADIAN", None) => Ok(1.0),
            _ => Err("UnitsUndetermined"),
        };
    }
    if let Some(cb) = i.record("CONVERSION_BASED_UNIT") {
        let [_, m] = cb.parameters.as_slice() else {
            return Err("UnitsUndetermined");
        };
        let m = x.instances.get(&reference(m)?).ok_or("UnitsUndetermined")?;
        let (value, base) = measure_with_unit(m)?;
        let factor = real(value)? * unit_factor(x, base, measure, depth + 1)?;
        return if factor.is_finite() && factor > 0.0 {
            Ok(factor)
        } else {
            Err("UnitsUndetermined")
        };
    }
    Err("UnitsUndetermined")
}

/// A measure's value and unit: the attributes of `MEASURE_WITH_UNIT`, in its
/// own record when the instance is complex (`(LENGTH_MEASURE_WITH_UNIT()
/// MEASURE_WITH_UNIT(v, #u) UNCERTAINTY_MEASURE_WITH_UNIT('n', 'd'))`), else
/// the first two of the subtype's.
fn measure_with_unit(i: &Instance) -> Named<(&Parameter, u64)> {
    let record = i
        .record("MEASURE_WITH_UNIT")
        .or_else(|| {
            i.records
                .iter()
                .find(|r| r.name.ends_with("MEASURE_WITH_UNIT") && r.parameters.len() >= 2)
        })
        .ok_or("UnitsUndetermined")?;
    match record.parameters.as_slice() {
        [value, unit, ..] => Ok((value, reference(unit)?)),
        _ => Err("UnitsUndetermined"),
    }
}

/// The units of the context of the representation listing `item`.
fn units(x: &Exchange, contexts: &BTreeMap<u64, u64>, item: u64) -> Named<Units> {
    let Some(context) = contexts.get(&item).and_then(|c| x.instances.get(c)) else {
        return Ok(DEFAULT_UNITS);
    };
    let mut out = DEFAULT_UNITS;
    let kind = |n: u64, name: &str| {
        x.instances
            .get(&n)
            .is_some_and(|i| i.record(name).is_some())
    };
    if let Some(assigned) = context.record("GLOBAL_UNIT_ASSIGNED_CONTEXT") {
        for u in list(assigned.parameters.first().ok_or("UnitsUndetermined")?)? {
            let u = reference(u)?;
            if kind(u, "LENGTH_UNIT") {
                out.length = unit_factor(x, u, Measure::Length, 0)?;
            } else if kind(u, "PLANE_ANGLE_UNIT") {
                out.angle = unit_factor(x, u, Measure::Angle, 0)?;
            }
        }
    }
    if let Some(assigned) = context.record("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT") {
        for u in list(assigned.parameters.first().ok_or("UnitsUndetermined")?)? {
            let Some(measure) = x
                .instances
                .get(&reference(u)?)
                .filter(|i| i.record("UNCERTAINTY_MEASURE_WITH_UNIT").is_some())
            else {
                continue;
            };
            let (value, unit) = measure_with_unit(measure)?;
            if kind(unit, "LENGTH_UNIT") {
                let t = real(value)? * unit_factor(x, unit, Measure::Length, 0)?;
                if !t.is_finite() || t <= 0.0 {
                    return Err("UnitsUndetermined");
                }
                // OCCT's read.precision.mode "File", read.maxprecision.val 1.
                out.tolerance = t.clamp(1e-7, 1.0);
                break;
            }
        }
    }
    Ok(out)
}

/// A curve an edge lies on, in millimetres.
#[derive(Debug, Clone, Copy)]
enum Curve {
    Line {
        p: [f64; 3],
        d: [f64; 3],
    },
    Circle {
        axes: Axes,
        r: f64,
    },
    /// `o + a1 cos t x + a2 sin t y` (STEP-b), the semi-axes in either
    /// order.
    Ellipse {
        axes: Axes,
        a1: f64,
        a2: f64,
    },
    /// A B-spline curve, the entity of its data (STEP-b).
    BSpline(u64),
}

/// A surface a face lies on, in millimetres and radians.
#[derive(Debug, Clone, Copy)]
enum Carrier {
    Plane(Axes),
    Cylinder(Axes, f64),
    Cone(Axes, f64, f64),
    Sphere(Axes, f64),
    Torus(Axes, f64, f64),
}

impl Carrier {
    fn axes(&self) -> Axes {
        match *self {
            Self::Plane(a)
            | Self::Cylinder(a, _)
            | Self::Cone(a, ..)
            | Self::Sphere(a, _)
            | Self::Torus(a, ..) => a,
        }
    }
    /// The surface's inverse map (OCCT's and the kernel's parameterisation),
    /// `u` in `(-pi, pi]`.
    fn uv(&self, p: [f64; 3]) -> [f64; 2] {
        let [a, b, h] = self.axes().local(p);
        match *self {
            Self::Plane(_) => [a, b],
            Self::Cylinder(..) => [b.atan2(a), h],
            Self::Cone(_, r, angle) => {
                let v = h / angle.cos();
                if r + v * angle.sin() < 0.0 {
                    [(-b).atan2(-a), v]
                } else {
                    [b.atan2(a), v]
                }
            }
            Self::Sphere(..) => [b.atan2(a), h.atan2(a.hypot(b))],
            Self::Torus(_, major, _) => [b.atan2(a), h.atan2(a.hypot(b) - major)],
        }
    }
    /// The `v` of a pole at `p`: a sphere's poles and a cone's apex, where
    /// every `u` meets.
    fn pole(&self, p: [f64; 3], tol: f64) -> Option<f64> {
        let [a, b, h] = self.axes().local(p);
        match *self {
            Self::Sphere(..) if a.hypot(b) <= tol => Some(FRAC_PI_2.copysign(h)),
            Self::Cone(_, r, angle) => {
                let v = h / angle.cos();
                ((r + v * angle.sin()).abs() <= tol).then(|| -r / angle.sin())
            }
            _ => None,
        }
    }
    fn periodic_v(&self) -> bool {
        matches!(self, Self::Torus(..))
    }
    /// The surface's record in OCCT's structure.
    fn record(&self) -> brep::Surface {
        let Axes { o, x, y, z } = self.axes();
        match *self {
            Self::Plane(_) => brep::Surface::Plane { p: o, n: z, x, y },
            Self::Cylinder(_, r) => brep::Surface::Cylinder {
                p: o,
                n: z,
                x,
                y,
                r,
            },
            Self::Cone(_, r, a) => brep::Surface::Cone {
                p: o,
                n: z,
                x,
                y,
                r,
                a,
            },
            Self::Sphere(_, r) => brep::Surface::Sphere {
                p: o,
                n: z,
                x,
                y,
                r,
            },
            Self::Torus(_, major, minor) => brep::Surface::Torus {
                p: o,
                n: z,
                x,
                y,
                major,
                minor,
            },
        }
    }
}

/// A use's pcurve in the walking direction: a segment in `(u, v)`, with
/// whether each end is at a pole; or, with `sinus`, the sinusoid `v = a0 +
/// a1 cos u + a2 sin u` from `a` to `b` (an ellipse on a cylinder, STEP-b),
/// which a shift by whole turns in `u` leaves unchanged.
#[derive(Debug, Clone, Copy)]
struct Segment {
    a: [f64; 2],
    b: [f64; 2],
    pole_a: bool,
    pole_b: bool,
    sinus: Option<[f64; 3]>,
}

/// An edge record already made, by its `EDGE_CURVE`.
#[derive(Debug, Clone, Copy)]
struct EdgeInfo {
    record: usize,
    /// The edge runs along its basis curve from `EDGE_CURVE`'s start.
    along_basis: bool,
    start: u64,
    end: u64,
    curve: Curve,
    range: [f64; 2],
    /// The `EDGE_CURVE`'s geometry, which may carry the file's pcurves.
    geometry: u64,
    /// The record's first and last points (its vertices), and whether they
    /// are one vertex.
    from: [f64; 3],
    to: [f64; 3],
    closed: bool,
}

/// One `ORIENTED_EDGE` of a loop.
struct Use {
    edge: EdgeInfo,
    /// Traversed along the edge's curve.
    along: bool,
    /// The vertex it ends at, walking the loop as written.
    head: u64,
}

pub(super) struct Build<'a> {
    x: &'a Exchange,
    units: Units,
    doc: Document,
    vertices: BTreeMap<u64, usize>,
    edges: BTreeMap<u64, EdgeInfo>,
}

fn orient(forward: bool) -> Orient {
    if forward {
        Orient::Forward
    } else {
        Orient::Reversed
    }
}

impl<'a> Build<'a> {
    pub(super) fn instance(&self, n: u64) -> Named<&'a Instance> {
        // The reader checked every reference.
        let x: &'a Exchange = self.x;
        x.instances.get(&n).ok_or("MalformedEntity")
    }

    /// The parameters of `n`'s record `name`, exactly `count` of them;
    /// `n`'s own name when it is another entity.
    pub(super) fn params(&self, n: u64, name: &str, count: usize) -> Named<&'a [Parameter]> {
        let i = self.instance(n)?;
        let r = i.record(name).ok_or_else(|| name_of(i))?;
        if r.parameters.len() == count {
            Ok(&r.parameters)
        } else {
            Err("MalformedEntity")
        }
    }

    fn push(&mut self, kind: Kind, data: Data, subs: Vec<Sub>) -> usize {
        self.doc.shapes.push(Shape { kind, data, subs });
        self.doc.shapes.len() - 1
    }

    fn finite(&self, v: [f64; 3]) -> Named<[f64; 3]> {
        if v.iter().all(|c| c.is_finite()) {
            Ok(v)
        } else {
            Err("NonFiniteGeometry")
        }
    }

    pub(super) fn point(&self, n: u64) -> Named<[f64; 3]> {
        let p = self.params(n, "CARTESIAN_POINT", 2)?;
        let [a, b, c] = list(&p[1])? else {
            return Err("MalformedEntity");
        };
        let s = self.units.length;
        self.finite([real(a)? * s, real(b)? * s, real(c)? * s])
    }

    /// A point in a surface's parameters (a pcurve's), unscaled.
    pub(super) fn point2(&self, n: u64) -> Named<[f64; 2]> {
        let p = self.params(n, "CARTESIAN_POINT", 2)?;
        let [a, b] = list(&p[1])? else {
            return Err("MalformedEntity");
        };
        let q = [real(a)?, real(b)?];
        if q.iter().all(|c| c.is_finite()) {
            Ok(q)
        } else {
            Err("NonFiniteGeometry")
        }
    }

    fn direction(&self, n: u64) -> Named<[f64; 3]> {
        let p = self.params(n, "DIRECTION", 2)?;
        let [a, b, c] = list(&p[1])? else {
            return Err("MalformedEntity");
        };
        unit([real(a)?, real(b)?, real(c)?])
    }

    fn length(&self, p: &Parameter) -> Named<f64> {
        let r = real(p)? * self.units.length;
        if r.is_finite() && r > 0.0 {
            Ok(r)
        } else {
            Err("DegenerateRadius")
        }
    }

    /// `AXIS2_PLACEMENT_3D` by ISO 10303-42's `build_axes`: `z` the axis or
    /// `(0, 0, 1)`, `x` the reference direction (or `(1, 0, 0)`, `(0, 1, 0)`
    /// when the axis is `+-(1, 0, 0)`) made normal to `z`; a reference
    /// direction along the axis is `DegeneratePlacement`.
    fn placement(&self, n: u64) -> Named<Axes> {
        let p = self.params(n, "AXIS2_PLACEMENT_3D", 4)?;
        let o = self.point(reference(&p[1])?)?;
        let z = match &p[2] {
            Parameter::Unset => [0.0, 0.0, 1.0],
            q => self.direction(reference(q)?)?,
        };
        let v = match &p[3] {
            Parameter::Unset if z == [1.0, 0.0, 0.0] || z == [-1.0, 0.0, 0.0] => [0.0, 1.0, 0.0],
            Parameter::Unset => [1.0, 0.0, 0.0],
            q => self.direction(reference(q)?)?,
        };
        let k = dot(v, z);
        let x =
            unit(std::array::from_fn(|i| v[i] - k * z[i])).map_err(|_| "DegeneratePlacement")?;
        Ok(Axes {
            o,
            x,
            y: cross(z, x),
            z,
        })
    }

    /// An edge's curve and whether it runs against its basis (a trimmed
    /// curve of opposite sense).
    fn curve(&self, n: u64, depth: usize) -> Named<(Curve, bool)> {
        if depth > 8 {
            return Err("MalformedEntity");
        }
        let i = self.instance(n)?;
        if i.record("LINE").is_some() {
            let p = self.params(n, "LINE", 3)?;
            let origin = self.point(reference(&p[1])?)?;
            let v = self.params(reference(&p[2])?, "VECTOR", 3)?;
            if real(&v[2])? <= 0.0 {
                return Err("DegenerateDirection");
            }
            let d = self.direction(reference(&v[1])?)?;
            return Ok((Curve::Line { p: origin, d }, false));
        }
        if i.record("CIRCLE").is_some() {
            let p = self.params(n, "CIRCLE", 3)?;
            let axes = self.placement(reference(&p[1])?)?;
            return Ok((
                Curve::Circle {
                    axes,
                    r: self.length(&p[2])?,
                },
                false,
            ));
        }
        if i.record("ELLIPSE").is_some() {
            let p = self.params(n, "ELLIPSE", 4)?;
            let axes = self.placement(reference(&p[1])?)?;
            return Ok((
                Curve::Ellipse {
                    axes,
                    a1: self.length(&p[2])?,
                    a2: self.length(&p[3])?,
                },
                false,
            ));
        }
        if spline::is_curve(i) {
            return Ok((Curve::BSpline(n), false));
        }
        if let Some(name) = spline::knotless(i) {
            return Err(name);
        }
        for name in ["SURFACE_CURVE", "SEAM_CURVE"] {
            if i.record(name).is_some() {
                let p = self.params(n, name, 4)?;
                return self.curve(reference(&p[1])?, depth + 1);
            }
        }
        if i.record("TRIMMED_CURVE").is_some() {
            let p = self.params(n, "TRIMMED_CURVE", 6)?;
            let (c, reversed) = self.curve(reference(&p[1])?, depth + 1)?;
            return Ok((c, reversed == boolean(&p[4])?));
        }
        Err(name_of(i))
    }

    fn surface(&self, n: u64) -> Named<Carrier> {
        let i = self.instance(n)?;
        let one = |name: &str, count: usize| -> Named<Option<&[Parameter]>> {
            Ok(match i.record(name) {
                Some(_) => Some(self.params(n, name, count)?),
                None => None,
            })
        };
        if let Some(p) = one("PLANE", 2)? {
            return Ok(Carrier::Plane(self.placement(reference(&p[1])?)?));
        }
        if let Some(p) = one("CYLINDRICAL_SURFACE", 3)? {
            return Ok(Carrier::Cylinder(
                self.placement(reference(&p[1])?)?,
                self.length(&p[2])?,
            ));
        }
        if let Some(p) = one("CONICAL_SURFACE", 4)? {
            let r = real(&p[2])? * self.units.length;
            let a = real(&p[3])? * self.units.angle;
            if !(r >= 0.0 && r.is_finite() && a > 0.0 && a < FRAC_PI_2) {
                return Err("InvalidConicalSurface");
            }
            return Ok(Carrier::Cone(self.placement(reference(&p[1])?)?, r, a));
        }
        if let Some(p) = one("SPHERICAL_SURFACE", 3)? {
            return Ok(Carrier::Sphere(
                self.placement(reference(&p[1])?)?,
                self.length(&p[2])?,
            ));
        }
        if let Some(p) = one("TOROIDAL_SURFACE", 4)? {
            if real(&p[2])? < 0.0 {
                return Err("NegativeMajorRadius");
            }
            let (major, minor) = (self.length(&p[2])?, self.length(&p[3])?);
            if major <= minor {
                return Err("NonRingToroidalSurface");
            }
            return Ok(Carrier::Torus(
                self.placement(reference(&p[1])?)?,
                major,
                minor,
            ));
        }
        Err(name_of(i))
    }

    fn vertex(&mut self, n: u64) -> Named<usize> {
        if let Some(v) = self.vertices.get(&n) {
            return Ok(*v);
        }
        let p = self.params(n, "VERTEX_POINT", 2)?;
        let point = self.point(reference(&p[1])?)?;
        let record = self.push(
            Kind::Vertex,
            Data::Vertex {
                tolerance: self.units.tolerance,
                point,
            },
            Vec::new(),
        );
        self.vertices.insert(n, record);
        Ok(record)
    }

    /// The edge record of an `EDGE_CURVE`: along its basis curve, the
    /// curve's range between the vertices (a closed circle one turn).
    fn edge(&mut self, n: u64) -> Named<EdgeInfo> {
        if let Some(e) = self.edges.get(&n) {
            return Ok(*e);
        }
        let p = self.params(n, "EDGE_CURVE", 5)?;
        let (v1, v2) = (reference(&p[1])?, reference(&p[2])?);
        let geometry = reference(&p[3])?;
        let (curve, reversed) = self.curve(geometry, 0)?;
        let along_basis = boolean(&p[4])? != reversed;
        let (start, end) = if along_basis { (v1, v2) } else { (v2, v1) };
        let (a, b) = (self.vertex(start)?, self.vertex(end)?);
        let at = |record: usize| match self.doc.shapes[record].data {
            Data::Vertex { point, .. } => point,
            _ => unreachable!("a vertex record"),
        };
        let (pa, pb) = (at(a), at(b));
        let (record_curve, range) = match curve {
            Curve::Line { p, d } => {
                let (f, l) = (dot(sub(pa, p), d), dot(sub(pb, p), d));
                if f == l {
                    return Err("ZeroLengthEdge");
                }
                (brep::Curve3::Line { p, d }, [f, l])
            }
            Curve::Circle { axes, r } => {
                let f = axes.angle(pa);
                let sweep = if start == end {
                    TAU
                } else {
                    (axes.angle(pb) - f).rem_euclid(TAU)
                };
                if sweep <= 0.0 {
                    return Err("ZeroLengthEdge");
                }
                (
                    brep::Curve3::Circle {
                        p: axes.o,
                        n: axes.z,
                        x: axes.x,
                        y: axes.y,
                        r,
                    },
                    [f, f + sweep],
                )
            }
            Curve::Ellipse { axes, a1, a2 } => {
                let angle = |q: [f64; 3]| {
                    let [x, y, _] = axes.local(q);
                    (y / a2).atan2(x / a1)
                };
                let f = angle(pa);
                let sweep = if start == end {
                    TAU
                } else {
                    (angle(pb) - f).rem_euclid(TAU)
                };
                if sweep <= 0.0 {
                    return Err("ZeroLengthEdge");
                }
                (
                    brep::Curve3::Ellipse {
                        p: axes.o,
                        n: axes.z,
                        x: axes.x,
                        y: axes.y,
                        major: a1,
                        minor: a2,
                    },
                    [f, f + sweep],
                )
            }
            Curve::BSpline(entity) => {
                let (record, kernel) = spline::curve3(self, entity)?;
                let range =
                    spline::curve_range(&kernel, pa, pb, start == end, self.units.tolerance)?;
                (brep::Curve3::BSpline(record), range)
            }
        };
        self.doc.curves.push(record_curve);
        let rep = EdgeRep::Curve {
            curve: self.doc.curves.len(),
            location: 0,
            range,
        };
        let record = self.push(
            Kind::Edge,
            Data::Edge {
                tolerance: self.units.tolerance,
                degenerated: false,
                reps: vec![rep],
            },
            vec![
                Sub {
                    orient: Orient::Forward,
                    shape: a,
                    location: 0,
                },
                Sub {
                    orient: Orient::Reversed,
                    shape: b,
                    location: 0,
                },
            ],
        );
        let info = EdgeInfo {
            record,
            along_basis,
            start: v1,
            end: v2,
            curve,
            range,
            geometry,
            from: pa,
            to: pb,
            closed: start == end,
        };
        self.edges.insert(n, info);
        Ok(info)
    }

    /// A use's pcurve on a curved surface, in the walking direction: the
    /// image of a ruling, parallel or meridian under the inverse map.
    fn segment(&self, carrier: &Carrier, u: &Use) -> Named<Segment> {
        let tol = self.units.tolerance;
        let [f, l] = u.edge.range;
        let angle_tol = 1e-9;
        let axis = carrier.axes();
        // `a` against `b` modulo a turn.
        let same = |a: f64, b: f64, radius: f64| {
            let d = (a - b + PI).rem_euclid(TAU) - PI;
            d.abs() * radius.max(tol) <= 10.0 * tol
        };
        let mut s = match (u.edge.curve, carrier) {
            (Curve::Line { p, d }, Carrier::Cylinder(..) | Carrier::Cone(..)) => {
                let (p0, p1, pm) = (along(p, d, f), along(p, d, l), along(p, d, 0.5 * (f + l)));
                let um = carrier.uv(pm)[0];
                let radius = |q: [f64; 3]| {
                    let [a, b, _] = axis.local(q);
                    a.hypot(b)
                };
                for q in [p0, p1] {
                    if carrier.pole(q, tol).is_none() && !same(carrier.uv(q)[0], um, radius(q)) {
                        return Err("PCurveNotDerived");
                    }
                }
                let v = |q: [f64; 3]| carrier.pole(q, tol).unwrap_or(carrier.uv(q)[1]);
                Segment {
                    a: [um, v(p0)],
                    b: [um, v(p1)],
                    pole_a: carrier.pole(p0, tol).is_some(),
                    pole_b: carrier.pole(p1, tol).is_some(),
                    sinus: None,
                }
            }
            (Curve::Circle { axes, r }, _) => {
                let sweep = l - f;
                let (p0, p1) = (axes.circle(r, f), axes.circle(r, l));
                let centre = axis.local(axes.o);
                let parallel = norm(cross(axes.z, axis.z)) <= angle_tol
                    && centre[0].hypot(centre[1]) <= 10.0 * tol;
                let meridian = dot(axes.z, axis.z).abs() <= angle_tol
                    && match *carrier {
                        Carrier::Sphere(..) => norm(centre) <= 10.0 * tol,
                        Carrier::Torus(_, major, _) => {
                            centre[2].abs() <= 10.0 * tol
                                && (centre[0].hypot(centre[1]) - major).abs() <= 10.0 * tol
                        }
                        _ => false,
                    };
                if parallel {
                    let [u0, v0] = carrier.uv(p0);
                    let du = sweep * dot(axes.z, axis.z).signum();
                    Segment {
                        a: [u0, v0],
                        b: [u0 + du, v0],
                        pole_a: false,
                        pole_b: false,
                        sinus: None,
                    }
                } else if meridian {
                    // The meridian plane's longitude: the centre's on a
                    // torus, the arc's middle on a sphere (an arc across a
                    // pole changes longitude and is not derived).
                    let um = match *carrier {
                        Carrier::Torus(..) => centre[1].atan2(centre[0]),
                        _ => carrier.uv(axes.circle(r, 0.5 * (f + l)))[0],
                    };
                    let e = [
                        um.cos() * axis.x[0] + um.sin() * axis.y[0],
                        um.cos() * axis.x[1] + um.sin() * axis.y[1],
                        um.cos() * axis.x[2] + um.sin() * axis.y[2],
                    ];
                    let slope = dot(axes.z, cross(e, axis.z)).signum();
                    let v0 = carrier.uv(p0)[1];
                    let v1 = v0 + slope * sweep;
                    let radius = match *carrier {
                        Carrier::Torus(_, major, _) => major,
                        _ => r,
                    };
                    for q in [p0, p1] {
                        if carrier.pole(q, tol).is_none() && !same(carrier.uv(q)[0], um, radius) {
                            return Err("PCurveNotDerived");
                        }
                    }
                    if matches!(carrier, Carrier::Sphere(..)) && v1.abs() > FRAC_PI_2 + 1e-9 {
                        return Err("PCurveNotDerived");
                    }
                    Segment {
                        a: [um, v0],
                        b: [um, v1],
                        pole_a: carrier.pole(p0, tol).is_some(),
                        pole_b: carrier.pole(p1, tol).is_some(),
                        sinus: None,
                    }
                } else {
                    return Err("PCurveNotDerived");
                }
            }
            (Curve::Ellipse { axes, a1, a2 }, &Carrier::Cylinder(_, r)) => {
                // The section of the cylinder by the ellipse's plane: its
                // centre on the axis, its axes projecting onto the axis's
                // normal plane as two perpendicular radii. Then the angle
                // about the axis is `phi +- t` and the height the sinusoid
                // of the plane `n . (p - c) = 0` over it.
                let c = axis.local(axes.o);
                let flat = |d: [f64; 3], s: f64| [s * dot(d, axis.x), s * dot(d, axis.y)];
                let (px, py) = (flat(axes.x, a1), flat(axes.y, a2));
                let n = [
                    dot(axes.z, axis.x),
                    dot(axes.z, axis.y),
                    dot(axes.z, axis.z),
                ];
                let near = |x: f64, y: f64| (x - y).abs() <= 10.0 * tol;
                if n[2].abs() <= angle_tol
                    || !near(c[0].hypot(c[1]), 0.0)
                    || !near(px[0].hypot(px[1]), r)
                    || !near(py[0].hypot(py[1]), r)
                    || !near((px[0] * py[0] + px[1] * py[1]) / r, 0.0)
                {
                    return Err("PCurveNotDerived");
                }
                let a = [
                    (n[0] * c[0] + n[1] * c[1] + n[2] * c[2]) / n[2],
                    -r * n[0] / n[2],
                    -r * n[1] / n[2],
                ];
                let v = |u: f64| a[0] + a[1] * u.cos() + a[2] * u.sin();
                let (sine, cosine) = f.sin_cos();
                let p0 = std::array::from_fn(|k| {
                    axes.o[k] + a1 * cosine * axes.x[k] + a2 * sine * axes.y[k]
                });
                let u0 = carrier.uv(p0)[0];
                let u1 = u0 + n[2].signum() * (l - f);
                Segment {
                    a: [u0, v(u0)],
                    b: [u1, v(u1)],
                    pole_a: false,
                    pole_b: false,
                    sinus: Some(a),
                }
            }
            _ => return Err("PCurveNotDerived"),
        };
        if !u.along {
            s = Segment {
                a: s.b,
                b: s.a,
                pole_a: s.pole_b,
                pole_b: s.pole_a,
                sinus: s.sinus,
            };
        }
        Ok(s)
    }
}

/// Place a loop's segments in the universal cover: each continues the last
/// (the nearest branch in `u`, and in `v` on a torus); at a pole the next
/// use's `u` is the nearest branch beyond the last in the direction that
/// keeps the face on the loop's left (`-u` at a pole above the face, `+u`
/// below it, reversed when the face lies on the right), and the gap is a
/// degenerated edge: `(after use, from, to)`.
fn walk(
    carrier: &Carrier,
    segments: &mut [Segment],
    left: bool,
) -> Vec<(usize, [f64; 2], [f64; 2])> {
    let mut degenerated = Vec::new();
    let nearest = |x: f64, target: f64| x + TAU * ((target - x) / TAU).round();
    for i in 1..segments.len() {
        let prev = segments[i - 1];
        let s = &mut segments[i];
        let shift_u = if prev.pole_b && s.pole_a {
            let above = prev.a[1] < s.a[1];
            let step = if above == left { -1.0 } else { 1.0 };
            let gap = (step * (s.a[0] - prev.b[0])).rem_euclid(TAU);
            let gap = if gap <= 1e-12 { TAU } else { gap };
            let to = prev.b[0] + step * gap;
            degenerated.push((i - 1, prev.b, [to, s.a[1]]));
            to - s.a[0]
        } else {
            nearest(s.a[0], prev.b[0]) - s.a[0]
        };
        let shift_v = if carrier.periodic_v() {
            nearest(s.a[1], prev.b[1]) - s.a[1]
        } else {
            0.0
        };
        s.a = [s.a[0] + shift_u, s.a[1] + shift_v];
        s.b = [s.b[0] + shift_u, s.b[1] + shift_v];
    }
    if let (Some(first), Some(last)) = (segments.first(), segments.last()) {
        if last.pole_b && first.pole_a && (last.b[0] - first.a[0]).abs() > 1e-12 {
            degenerated.push((segments.len() - 1, last.b, first.a));
        }
    }
    degenerated
}

/// A pcurve record (1-based) and its range.
type Placed = (usize, [f64; 2]);

fn line2(doc: &mut Document, a: [f64; 2], b: [f64; 2]) -> usize {
    doc.curves2d.push(brep::Curve2::Line {
        p: a,
        d: [b[0] - a[0], b[1] - a[1]],
    });
    doc.curves2d.len()
}

/// A use's derived pcurve over `[0, 1]`: a segment, or a sinusoid from `a`
/// to `b` in `u`.
fn segment2(doc: &mut Document, a: [f64; 2], b: [f64; 2], sinus: Option<[f64; 3]>) -> usize {
    match sinus {
        None => line2(doc, a, b),
        Some(c) => {
            doc.curves2d.push(brep::Curve2::Sinusoid {
                u0: a[0],
                du: b[0] - a[0],
                a: c,
            });
            doc.curves2d.len()
        }
    }
}

impl<'a> Build<'a> {
    /// The file's pcurve of an edge geometry on the surface `surface`: the
    /// first `PCURVE` of its `SURFACE_CURVE` or `SEAM_CURVE` (through
    /// trimmed curves) whose basis is that surface (STEP-b).
    fn pcurve_on(&self, geometry: u64, surface: u64) -> Named<spline::Pcurve> {
        let mut n = geometry;
        for _ in 0..8 {
            let i = self.instance(n)?;
            if i.record("TRIMMED_CURVE").is_some() {
                n = reference(&self.params(n, "TRIMMED_CURVE", 6)?[1])?;
                continue;
            }
            let Some(name) = ["SURFACE_CURVE", "SEAM_CURVE"]
                .into_iter()
                .find(|k| i.record(k).is_some())
            else {
                return Err("PCurveNotDerived");
            };
            for a in list(&self.params(n, name, 4)?[2])? {
                let a = reference(a)?;
                if self.instance(a)?.record("PCURVE").is_some()
                    && reference(&self.params(a, "PCURVE", 3)?[1])? == surface
                {
                    return Ok(spline::pcurve(self, a)?.1);
                }
            }
            return Err("PCurveNotDerived");
        }
        Err("MalformedEntity")
    }

    /// A face: its surface record, its wires with the derived pcurves (or,
    /// on a spline surface, the file's) and degenerated edges, its
    /// orientation in the shell (`same_sense`).
    fn face(&mut self, n: u64) -> Named<(usize, bool)> {
        let i = self.instance(n)?;
        let name = ["ADVANCED_FACE", "FACE_SURFACE"]
            .into_iter()
            .find(|k| i.record(k).is_some())
            .ok_or_else(|| name_of(i))?;
        let p = self.params(n, name, 4)?;
        let same_sense = boolean(&p[3])?;
        let surface_entity = reference(&p[2])?;
        let si = self.instance(surface_entity)?;
        // An elementary surface, or a spline's kernel surface (STEP-b).
        let (carrier, spline_surface) = if spline::is_surface(si) {
            let (record, kernel) = spline::surface(self, surface_entity)?;
            self.doc
                .surfaces
                .push(brep::Surface::BSpline(Box::new(record)));
            (None, Some(kernel))
        } else if let Some(name) = spline::knotless(si) {
            return Err(name);
        } else {
            let carrier = self.surface(surface_entity)?;
            self.doc.surfaces.push(carrier.record());
            (Some(carrier), None)
        };
        let surface = self.doc.surfaces.len();
        let tol = self.units.tolerance;
        let bounds: Vec<u64> = list(&p[1])?.iter().map(reference).collect::<Named<_>>()?;
        let mut wires = Vec::new();
        // Per edge record, the pcurves of its uses in this face by stored
        // orientation (OCCT's first and second pcurve on a closed surface).
        let mut pcurves: BTreeMap<usize, [Option<Placed>; 2]> = BTreeMap::new();
        let mut uses_in_face: BTreeMap<usize, usize> = BTreeMap::new();
        let mut winds = false;
        for b in bounds {
            let bi = self.instance(b)?;
            let kind = ["FACE_OUTER_BOUND", "FACE_BOUND"]
                .into_iter()
                .find(|k| bi.record(k).is_some())
                .ok_or_else(|| name_of(bi))?;
            let bp = self.params(b, kind, 3)?;
            let (lp, bound_sense) = (reference(&bp[1])?, boolean(&bp[2])?);
            let li = self.instance(lp)?;
            if li.record("VERTEX_LOOP").is_some() {
                // The whole sphere or torus (ISO 10303-42), as OCCT reads it.
                if matches!(carrier, Some(Carrier::Sphere(..) | Carrier::Torus(..)))
                    && list(&p[1])?.len() == 1
                {
                    continue;
                }
                return Err("VERTEX_LOOP");
            }
            let lparams = self.params(lp, "EDGE_LOOP", 2)?;
            let mut uses = Vec::new();
            for oe in list(&lparams[1])? {
                let oe = reference(oe)?;
                let q = self.params(oe, "ORIENTED_EDGE", 5)?;
                let (element, sense) = (reference(&q[3])?, boolean(&q[4])?);
                let edge = self.edge(element)?;
                uses.push(Use {
                    edge,
                    along: sense == edge.along_basis,
                    head: if sense { edge.end } else { edge.start },
                });
            }
            if uses.is_empty() {
                return Err("EmptyLoop");
            }
            let wire_forward = bound_sense == same_sense;
            // Pcurves on curved surfaces, placed in the universal cover.
            let mut segments = Vec::new();
            let mut degenerated = Vec::new();
            if let Some(carrier) = carrier.filter(|c| !matches!(c, Carrier::Plane(_))) {
                for u in &uses {
                    segments.push(self.segment(&carrier, u)?);
                }
                degenerated = walk(&carrier, &mut segments, bound_sense == same_sense);
                let (first, last) = (segments[0].a, segments[segments.len() - 1].b);
                winds |= !segments[0].pole_a
                    && ((last[0] - first[0]).abs() > 1.0
                        || carrier.periodic_v() && (last[1] - first[1]).abs() > 1.0);
            }
            // On a spline surface, the file's pcurves, each over the range
            // where its image meets the edge's vertices (STEP-b).
            let mut placed = Vec::new();
            if let Some(kernel) = &spline_surface {
                for u in &uses {
                    let pcurve = self.pcurve_on(u.edge.geometry, surface_entity)?;
                    let (from, to) = (u.edge.from, u.edge.to);
                    placed.push(Some(spline::pcurve_range(
                        kernel,
                        pcurve,
                        from,
                        to,
                        u.edge.closed,
                        tol,
                    )?));
                }
            }
            let mut subs = Vec::new();
            for (k, u) in uses.iter().enumerate() {
                subs.push(Sub {
                    orient: orient(u.along),
                    shape: u.edge.record,
                    location: 0,
                });
                *uses_in_face.entry(u.edge.record).or_default() += 1;
                let made = if let Some(s) = segments.get(k) {
                    let (a, b) = if u.along { (s.a, s.b) } else { (s.b, s.a) };
                    Some((segment2(&mut self.doc, a, b, s.sinus), [0.0, 1.0]))
                } else if let Some((record, range)) = placed.get_mut(k).and_then(Option::take) {
                    self.doc.curves2d.push(record);
                    Some((self.doc.curves2d.len(), range))
                } else {
                    None
                };
                if let Some(c2) = made {
                    // BRep_Tool::CurveOnSurface: the second pcurve serves the
                    // use stored reversed in the unoriented face.
                    let stored = u.along == wire_forward;
                    let slot = &mut pcurves.entry(u.edge.record).or_default()[usize::from(!stored)];
                    if slot.is_some() {
                        return Err("EdgeUsedTwiceAlike");
                    }
                    *slot = Some(c2);
                }
                for &(after, from, to) in degenerated.iter().filter(|d| d.0 == k) {
                    let vertex = self.vertex(u.head)?;
                    let c2 = line2(&mut self.doc, from, to);
                    let rep = EdgeRep::OnSurface {
                        pcurves: vec![c2],
                        surface,
                        location: 0,
                        range: [0.0, 1.0],
                    };
                    debug_assert_eq!(after, k);
                    let record = self.push(
                        Kind::Edge,
                        Data::Edge {
                            tolerance: tol,
                            degenerated: true,
                            reps: vec![rep],
                        },
                        vec![
                            Sub {
                                orient: Orient::Forward,
                                shape: vertex,
                                location: 0,
                            },
                            Sub {
                                orient: Orient::Reversed,
                                shape: vertex,
                                location: 0,
                            },
                        ],
                    );
                    // Walking the loop as written, whatever the wire's sense.
                    subs.push(Sub {
                        orient: Orient::Forward,
                        shape: record,
                        location: 0,
                    });
                }
            }
            let wire = self.push(Kind::Wire, Data::None, subs);
            wires.push(Sub {
                orient: orient(wire_forward),
                shape: wire,
                location: 0,
            });
        }
        if winds && uses_in_face.values().all(|n| *n < 2) {
            return Err("PeriodicFaceWithoutSeam");
        }
        for (record, slots) in pcurves {
            // Derived pcurves share `[0, 1]`; a spline face's edge has one.
            let (pcurves, range) = match slots {
                [Some(a), Some(b)] => (vec![a.0, b.0], a.1),
                [Some(a), None] | [None, Some(a)] => (vec![a.0], a.1),
                [None, None] => continue,
            };
            if let Data::Edge { reps, .. } = &mut self.doc.shapes[record].data {
                reps.push(EdgeRep::OnSurface {
                    pcurves,
                    surface,
                    location: 0,
                    range,
                });
            }
        }
        let face = self.push(
            Kind::Face,
            Data::Face {
                tolerance: tol,
                surface,
                location: 0,
            },
            wires,
        );
        Ok((face, same_sense))
    }

    /// A shell of faces; `ORIENTED_CLOSED_SHELL`'s orientation composes.
    fn shell(&mut self, n: u64) -> Named<(usize, bool)> {
        let i = self.instance(n)?;
        if i.record("ORIENTED_CLOSED_SHELL").is_some() {
            let p = self.params(n, "ORIENTED_CLOSED_SHELL", 4)?;
            let element = reference(&p[2])?;
            if self
                .instance(element)?
                .record("ORIENTED_CLOSED_SHELL")
                .is_some()
            {
                return Err("MalformedEntity");
            }
            let (shell, forward) = self.shell(element)?;
            return Ok((shell, forward == boolean(&p[3])?));
        }
        let name = ["CLOSED_SHELL", "OPEN_SHELL"]
            .into_iter()
            .find(|k| i.record(k).is_some())
            .ok_or_else(|| name_of(i))?;
        let p = self.params(n, name, 2)?;
        let mut subs = Vec::new();
        for f in list(&p[1])? {
            let f = reference(f)?;
            let fi = self.instance(f)?;
            let (face, forward) = if fi.record("ORIENTED_FACE").is_some() {
                let q = self.params(f, "ORIENTED_FACE", 4)?;
                let (face, same) = self.face(reference(&q[2])?)?;
                (face, same == boolean(&q[3])?)
            } else {
                self.face(f)?
            };
            subs.push(Sub {
                orient: orient(forward),
                shape: face,
                location: 0,
            });
        }
        Ok((self.push(Kind::Shell, Data::None, subs), true))
    }
}

/// OCCT's shape structure of one body: a solid of its shells, or a shell.
fn document(x: &Exchange, units: Units, item: u64, shell: Option<u64>) -> Named<Document> {
    let mut b = Build {
        x,
        units,
        doc: Document {
            version: 2,
            locations: Vec::new(),
            curves2d: Vec::new(),
            curves: Vec::new(),
            surfaces: Vec::new(),
            triangulations: 0,
            shapes: Vec::new(),
            root: Sub {
                orient: Orient::Forward,
                shape: 0,
                location: 0,
            },
        },
        vertices: BTreeMap::new(),
        edges: BTreeMap::new(),
    };
    let root = if let Some(s) = shell {
        let (record, forward) = b.shell(s)?;
        b.doc.root.orient = orient(forward);
        record
    } else {
        let i = b.instance(item)?;
        let shells: Vec<(u64, bool)> = if i.record("MANIFOLD_SOLID_BREP").is_some() {
            vec![(
                reference(&b.params(item, "MANIFOLD_SOLID_BREP", 2)?[1])?,
                false,
            )]
        } else {
            let p = b.params(item, "BREP_WITH_VOIDS", 3)?;
            let mut out = vec![(reference(&p[1])?, false)];
            for v in list(&p[2])? {
                out.push((reference(v)?, true));
            }
            out
        };
        let mut subs = Vec::new();
        for (s, void) in shells {
            if void && b.instance(s)?.record("ORIENTED_CLOSED_SHELL").is_none() {
                return Err("MalformedEntity");
            }
            let (record, forward) = b.shell(s)?;
            subs.push(Sub {
                orient: orient(forward),
                shape: record,
                location: 0,
            });
        }
        b.push(Kind::Solid, Data::None, subs)
    };
    b.doc.root.shape = root;
    Ok(b.doc)
}

/// Every solid and surface-model shell of an exchange structure in the cell
/// model, and every construct the kernel cannot represent, by name.
pub fn import(x: &Exchange) -> Result<StepImport, StepError> {
    let schemas = x.schemas();
    let known = |s: &&str| {
        let name = s
            .split('{')
            .next()
            .unwrap_or("")
            .trim()
            .to_ascii_uppercase();
        // AP214's committee drafts (`AUTOMOTIVE_DESIGN_CC2` and the like)
        // have its geometry and topology.
        SCHEMAS.contains(&name.as_str()) || name.starts_with("AUTOMOTIVE_DESIGN_CC")
    };
    if !schemas.iter().any(known) {
        return Err(StepError::Schema(schemas.join(", ")));
    }
    let mut unsupported: BTreeMap<&'static str, usize> = BTreeMap::new();
    // Each item's representation context: the first representation listing
    // it, by entity number.
    let mut contexts = BTreeMap::new();
    for i in x.instances.values() {
        for r in &i.records {
            if !r.name.ends_with("REPRESENTATION") {
                continue;
            }
            if let [_, Parameter::List(items), Parameter::Reference(context)] =
                r.parameters.as_slice()
            {
                for item in items {
                    if let Parameter::Reference(item) = item {
                        contexts.entry(*item).or_insert(*context);
                    }
                }
            }
        }
    }
    let placed = x.instances.values().any(|i| {
        i.records.iter().any(|r| {
            matches!(
                r.name.as_str(),
                "MAPPED_ITEM"
                    | "ITEM_DEFINED_TRANSFORMATION"
                    | "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
                    | "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"
            )
        })
    });
    if placed {
        unsupported.insert("AssemblyPlacement", 1);
    }
    let mut items = Vec::new();
    for (&n, i) in &x.instances {
        if i.record("MANIFOLD_SOLID_BREP").is_some() || i.record("BREP_WITH_VOIDS").is_some() {
            items.push((n, Item::Solid, None));
        } else if let Some(model) = i.record("SHELL_BASED_SURFACE_MODEL") {
            match model.parameters.as_slice() {
                [_, Parameter::List(shells)] => {
                    for s in shells {
                        items.push((n, Item::Shell, Some(reference(s).unwrap_or(0))));
                    }
                }
                _ => items.push((n, Item::Shell, Some(0))),
            }
        }
    }
    let mut bodies = Vec::new();
    for (entity, item, shell) in items {
        let units = units(x, &contexts, entity);
        let tolerance =
            Tolerance::new(units.map_or(1e-7, |u| u.tolerance), 1e-12).unwrap_or_default();
        let converted = units
            .and_then(|u| document(x, u, entity, shell))
            .map(|doc| {
                let mut imported = occt_brep::import(&doc);
                for (name, count) in &imported.unsupported {
                    *unsupported.entry(name).or_default() += count;
                }
                match item {
                    Item::Solid => imported.solids.pop().map(|s| (s.tolerance, s.result)),
                    Item::Shell => imported.free.pop().map(|s| (s.tolerance, s.result)),
                }
            });
        let (tolerance, result) = match converted {
            Ok(Some(done)) => done,
            Ok(None) => (tolerance, Err(Rejected::Unsupported(vec!["EmptyBody"]))),
            Err(name) => {
                *unsupported.entry(name).or_default() += 1;
                (tolerance, Err(Rejected::Unsupported(vec![name])))
            }
        };
        bodies.push(StepBody {
            entity,
            item,
            tolerance,
            result,
        });
    }
    Ok(StepImport {
        bodies,
        unsupported,
    })
}
