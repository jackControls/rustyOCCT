use crate::decide;
use crate::history::{self, History, Relation};
use crate::identity::{AlgorithmLevel, EntityId, OperationId, OperationKind};
use crate::math::finite;
use crate::profile::BoundaryKind;
use crate::topology::Topology;

mod attrs;
mod boolean;
mod enclose;
pub(crate) mod split;
mod stack;
use crate::{
    Boundary, Bounds3, Error, Frame3, Location, Point2, Point3, Profile, Result, RigidTransform,
    Tolerance,
};
pub use attrs::Context;
pub use split::Side;

/// Geometric properties at unit density, evaluated from analytic geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    /// Cubic millimetres.
    pub volume: f64,
    /// Square millimetres; includes hole walls and both caps.
    pub surface_area: f64,
    pub centroid: Point3,
    /// Central inertia tensor in world axes, at unit density (mm^5).
    pub inertia: [[f64; 3]; 3],
}

impl MassProperties {
    /// The same properties after a rigid motion: the centroid moved, the
    /// inertia turned (`Q I Q^T`).
    pub(crate) fn moved(&self, transform: RigidTransform) -> Self {
        let q: [[f64; 3]; 3] = [crate::Vec3::X, crate::Vec3::Y, crate::Vec3::Z]
            .map(|e| transform.vector(e).to_array());
        // q[j][i] is the i-th component of the j-th rotated axis: Q = q^T.
        let inertia = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                let mut s = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        s += q[i][a] * self.inertia[i][j] * q[j][b];
                    }
                }
                s
            })
        });
        Self {
            volume: self.volume,
            surface_area: self.surface_area,
            centroid: transform.point(self.centroid),
            inertia,
        }
    }
}

/// How a solid was made: a normal extrusion of a profile, a right circular
/// cone or frustum, a sphere or zone, or a torus, v-segment or wedge (S3 of
/// REVIEW_NOTES.md).
#[derive(Debug, Clone, PartialEq)]
enum Construction {
    Prism(Box<Profile>),
    Cone {
        bottom: f64,
        top: f64,
        tolerance: Tolerance,
    },
    /// The latitudes of its ends; `start` and `end` are their heights.
    Sphere {
        radius: f64,
        low: f64,
        high: f64,
        tolerance: Tolerance,
    },
    /// The tube's latitudes and the turn; `start` and `end` are the heights
    /// of a v-segment's ends.
    Torus {
        major: f64,
        minor: f64,
        low: f64,
        high: f64,
        angle: f64,
        tolerance: Tolerance,
    },
    /// A piece of a prism split by a plane oblique to its axis (S8a.2).
    Clipped(Box<split::Clipped>),
    /// A half of a cone or sphere zone split by a plane containing its
    /// axis (S8c.2).
    Half(Box<split::Half>),
    /// A solid of a Boolean of two prisms that is a stack of slabs of
    /// different regions (S9a.2).
    Stack(Box<boolean::stack::Stack>),
}

/// An immutable, validated normal extrusion of one planar material region,
/// or a cone or frustum.
/// The retained construction supports exact queries without a triangle mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Solid {
    construction: Construction,
    frame: Frame3,
    start: f64,
    end: f64,
    topology: Topology,
    mass: MassProperties,
    bounds: Bounds3,
    operation: OperationId,
}

impl Solid {
    /// Extrude along the profile frame's normal between two signed offsets.
    /// Ids derive from [`OperationId::UNSPECIFIED`]; the history is discarded.
    #[deprecated(note = "use Solid::extrude_with, which also returns the history")]
    pub fn extrude(profile: Profile, frame: Frame3, start: f64, end: f64) -> Result<Self> {
        Self::extrude_with(OperationId::UNSPECIFIED, profile, frame, start, end).map(|(s, _)| s)
    }

    /// Extrude along the profile frame's normal between two signed offsets.
    /// Either direction is allowed; zero or sub-tolerance thickness is rejected.
    ///
    /// Every entity id derives from `operation` and the profile's labels (or
    /// stored indices). The history is a construction: every vertex, edge and
    /// face is `Generated` from its profile element(s) with its role.
    pub fn extrude_with(
        operation: OperationId,
        profile: Profile,
        frame: Frame3,
        start: f64,
        end: f64,
    ) -> Result<(Self, History)> {
        Self::extrude_at(
            AlgorithmLevel::CURRENT,
            operation,
            profile,
            frame,
            start,
            end,
        )
    }

    /// [`Solid::extrude_with`] at a recorded algorithm level (H8), so a
    /// stored history replays identically; a level this build does not
    /// provide is an error.
    pub fn extrude_at(
        level: AlgorithmLevel,
        operation: OperationId,
        profile: Profile,
        frame: Frame3,
        start: f64,
        end: f64,
    ) -> Result<(Self, History)> {
        Self::extrude_in(
            &Context::new(operation).at(level),
            profile,
            frame,
            start,
            end,
        )
    }

    /// [`Solid::extrude_with`] in an operation context. A construction has
    /// no input entities, so it carries no attributes (generation never
    /// creates any).
    pub fn extrude_in(
        context: &Context,
        profile: Profile,
        frame: Frame3,
        start: f64,
        end: f64,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let solid = Self::build(operation, profile, frame, start, end)?;
        let t = &solid.topology;
        let relations = t
            .ids()
            .map(|(id, _)| {
                let d = t.derivation(id).expect("every id has a derivation");
                Relation::Generated {
                    from: d.parents.clone(),
                    to: id,
                    role: d.role,
                }
            })
            .collect();
        let history = History::new(
            operation,
            OperationKind::Extrude,
            Vec::new(),
            vec![t.body_id()],
            relations,
            Vec::new(),
        );
        let history = history.at_level(level);
        solid.debug_check(&[], &history);
        Ok((solid, history))
    }

    /// Every operation checks its history independently in debug builds.
    fn debug_check(&self, inputs: &[history::EntitySet], history: &History) {
        if cfg!(debug_assertions) {
            let outputs = [self.topology.entity_set(self.resolution())];
            let issues = history::check(inputs, &outputs, history);
            assert!(
                issues.is_empty(),
                "operation history is invalid: {issues:?}"
            );
        }
    }

    fn build(
        operation: OperationId,
        profile: Profile,
        frame: Frame3,
        start: f64,
        end: f64,
    ) -> Result<Self> {
        let tolerance = profile.tolerance();
        tolerance.resolve(&[start, end])?;
        let low = start.min(end);
        let high = start.max(end);
        let height = finite(high - low, "extrusion height")?;
        if decide::sum_le(&[high, -low], &[tolerance.linear()]) {
            return Err(Error::Degenerate("extrusion height"));
        }
        let bounds = bounds(&profile, frame, low, high);
        bounds.min.checked(tolerance)?;
        bounds.max.checked(tolerance)?;
        let mass = properties(&profile, frame, low, height)?;
        let topology = Topology::prism(&profile, frame, low, high, start < end, operation)?;
        Ok(Self {
            construction: Construction::Prism(Box::new(profile)),
            frame,
            start,
            end,
            topology,
            mass,
            bounds,
            operation,
        })
    }

    fn build_cone(
        operation: OperationId,
        frame: Frame3,
        bottom: f64,
        top: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<Self> {
        tolerance.resolve(&[bottom, top, height])?;
        let topology = Topology::cone(frame, bottom, top, height, tolerance, operation)?;
        let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let (x, y) = (frame.x().to_array(), frame.y().to_array());
        for (radius, z) in [(bottom, 0.0), (top, height)] {
            let centre = frame.point(Point2::default(), z).to_array();
            for i in 0..3 {
                let extent = radius * x[i].hypot(y[i]);
                min[i] = min[i].min(centre[i] - extent);
                max[i] = max[i].max(centre[i] + extent);
            }
        }
        let bounds = Bounds3 {
            min: Point3::new(min[0], min[1], min[2]),
            max: Point3::new(max[0], max[1], max[2]),
        };
        bounds.min.checked(tolerance)?;
        bounds.max.checked(tolerance)?;
        let mass = topology
            .mass_enclosure()
            .ok_or(Error::Unrepresentable("cone mass properties"))?
            .midpoints();
        Ok(Self {
            construction: Construction::Cone {
                bottom,
                top,
                tolerance,
            },
            frame,
            start: 0.0,
            end: height,
            topology,
            mass,
            bounds,
            operation,
        })
    }

    /// The same construction in another frame.
    fn rebuilt(&self, operation: OperationId, frame: Frame3) -> Result<Self> {
        match &self.construction {
            Construction::Prism(profile) => {
                Self::build(operation, (**profile).clone(), frame, self.start, self.end)
            }
            Construction::Cone {
                bottom,
                top,
                tolerance,
            } => Self::build_cone(operation, frame, *bottom, *top, self.end, *tolerance),
            Construction::Sphere {
                radius,
                low,
                high,
                tolerance,
            } => Self::build_sphere(operation, frame, *radius, *low, *high, *tolerance),
            Construction::Torus {
                major,
                minor,
                low,
                high,
                angle,
                tolerance,
            } => Self::build_torus(
                operation, frame, *major, *minor, *low, *high, *angle, *tolerance,
            ),
            Construction::Clipped(clipped) => {
                clipped.rebuilt(operation, frame, self.start, self.end)
            }
            Construction::Half(half) => half.rebuilt(operation, frame),
            Construction::Stack(stack) => stack.rebuilt_with(operation, frame, None),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build_torus(
        operation: OperationId,
        frame: Frame3,
        major: f64,
        minor: f64,
        low: f64,
        high: f64,
        angle: f64,
        tolerance: Tolerance,
    ) -> Result<Self> {
        tolerance.resolve(&[major, minor])?;
        let topology =
            Topology::torus(frame, major, minor, low, high, angle, tolerance, operation)?;
        let (start, end) = if high - low == std::f64::consts::TAU {
            (0.0, 0.0)
        } else {
            (
                crate::math::scaled_sin(minor, low),
                crate::math::scaled_sin(minor, high),
            )
        };
        // Within the whole torus's box: the discs of radius major + minor at
        // heights -minor and minor.
        let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let (x, y) = (frame.x().to_array(), frame.y().to_array());
        for z in [-minor, minor] {
            let centre = frame.point(Point2::default(), z).to_array();
            for i in 0..3 {
                let extent = (major + minor) * x[i].hypot(y[i]);
                min[i] = min[i].min(centre[i] - extent);
                max[i] = max[i].max(centre[i] + extent);
            }
        }
        let bounds = Bounds3 {
            min: Point3::new(min[0], min[1], min[2]),
            max: Point3::new(max[0], max[1], max[2]),
        };
        bounds.min.checked(tolerance)?;
        bounds.max.checked(tolerance)?;
        let mass = topology
            .mass_enclosure()
            .ok_or(Error::Unrepresentable("torus mass properties"))?
            .midpoints();
        Ok(Self {
            construction: Construction::Torus {
                major,
                minor,
                low,
                high,
                angle,
                tolerance,
            },
            frame,
            start,
            end,
            topology,
            mass,
            bounds,
            operation,
        })
    }

    /// A torus, v-segment or wedge on the frame (S3 of REVIEW_NOTES.md):
    /// the tube of radius `minor` about the circle of radius `major` round
    /// the frame's normal, its latitudes `low..high` (radians, from the
    /// outer equator towards the normal) revolved by `angle` from the frame's
    /// x, as `BRepPrimAPI_MakeTorus(gp_Ax2, major, minor, low, high, angle)`.
    /// A whole torus is one face without loops; a v-segment is a full turn,
    /// a wedge the whole tube from `v = 0`; both at once are out of domain.
    /// Mass properties are the general certified ones.
    #[allow(clippy::too_many_arguments)]
    pub fn torus_with(
        operation: OperationId,
        frame: Frame3,
        major: f64,
        minor: f64,
        low: f64,
        high: f64,
        angle: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::torus_in(
            &Context::new(operation),
            frame,
            major,
            minor,
            low,
            high,
            angle,
            tolerance,
        )
    }

    /// [`Solid::torus_with`] at a recorded algorithm level (H8).
    #[allow(clippy::too_many_arguments)]
    pub fn torus_at(
        level: AlgorithmLevel,
        operation: OperationId,
        frame: Frame3,
        major: f64,
        minor: f64,
        low: f64,
        high: f64,
        angle: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::torus_in(
            &Context::new(operation).at(level),
            frame,
            major,
            minor,
            low,
            high,
            angle,
            tolerance,
        )
    }

    /// [`Solid::torus_with`] in an operation context. Every entity is
    /// generated from its meridian element; nothing carries an attribute.
    #[allow(clippy::too_many_arguments)]
    pub fn torus_in(
        context: &Context,
        frame: Frame3,
        major: f64,
        minor: f64,
        low: f64,
        high: f64,
        angle: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let solid = Self::build_torus(operation, frame, major, minor, low, high, angle, tolerance)?;
        let history = solid.revolve_history(level);
        solid.debug_check(&[], &history);
        Ok((solid, history))
    }

    fn build_sphere(
        operation: OperationId,
        frame: Frame3,
        radius: f64,
        low: f64,
        high: f64,
        tolerance: Tolerance,
    ) -> Result<Self> {
        tolerance.resolve(&[radius])?;
        let topology = Topology::sphere(frame, radius, low, high, tolerance, operation)?;
        let (start, end) = (
            crate::math::scaled_sin(radius, low),
            crate::math::scaled_sin(radius, high),
        );
        let (start, end) = (
            if low == -std::f64::consts::FRAC_PI_2 {
                -radius
            } else {
                start
            },
            if high == std::f64::consts::FRAC_PI_2 {
                radius
            } else {
                end
            },
        );
        // The zone lies within the hull of discs of its largest radius at
        // both ends (the radius is largest at the equator when it is in).
        let widest = if low <= 0.0 && 0.0 <= high {
            radius
        } else {
            crate::math::scaled_cos(radius, low).max(crate::math::scaled_cos(radius, high))
        };
        let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let (x, y) = (frame.x().to_array(), frame.y().to_array());
        for z in [start, end] {
            let centre = frame.point(Point2::default(), z).to_array();
            for i in 0..3 {
                let extent = widest * x[i].hypot(y[i]);
                min[i] = min[i].min(centre[i] - extent);
                max[i] = max[i].max(centre[i] + extent);
            }
        }
        let bounds = Bounds3 {
            min: Point3::new(min[0], min[1], min[2]),
            max: Point3::new(max[0], max[1], max[2]),
        };
        bounds.min.checked(tolerance)?;
        bounds.max.checked(tolerance)?;
        let mass = topology
            .mass_enclosure()
            .ok_or(Error::Unrepresentable("sphere mass properties"))?
            .midpoints();
        Ok(Self {
            construction: Construction::Sphere {
                radius,
                low,
                high,
                tolerance,
            },
            frame,
            start,
            end,
            topology,
            mass,
            bounds,
            operation,
        })
    }

    /// A sphere or spherical zone on the frame (S3 of REVIEW_NOTES.md):
    /// centre at the frame origin, latitudes `low < high` in
    /// `[-pi/2, pi/2]` (radians) about its normal, an end at `+-pi/2` a pole,
    /// as `BRepPrimAPI_MakeSphere(gp_Ax2, radius, low, high)`. A whole sphere
    /// is one face without loops. Mass properties are the general certified
    /// ones.
    pub fn sphere_with(
        operation: OperationId,
        frame: Frame3,
        radius: f64,
        low: f64,
        high: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::sphere_in(
            &Context::new(operation),
            frame,
            radius,
            low,
            high,
            tolerance,
        )
    }

    /// [`Solid::sphere_with`] at a recorded algorithm level (H8).
    pub fn sphere_at(
        level: AlgorithmLevel,
        operation: OperationId,
        frame: Frame3,
        radius: f64,
        low: f64,
        high: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::sphere_in(
            &Context::new(operation).at(level),
            frame,
            radius,
            low,
            high,
            tolerance,
        )
    }

    /// [`Solid::sphere_with`] in an operation context. Every entity is
    /// generated from its meridian element; nothing carries an attribute.
    pub fn sphere_in(
        context: &Context,
        frame: Frame3,
        radius: f64,
        low: f64,
        high: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let solid = Self::build_sphere(operation, frame, radius, low, high, tolerance)?;
        let history = solid.revolve_history(level);
        solid.debug_check(&[], &history);
        Ok((solid, history))
    }

    /// A revolved primitive's construction history: every entity generated
    /// from its meridian element.
    fn revolve_history(&self, level: AlgorithmLevel) -> History {
        let t = &self.topology;
        let relations = t
            .ids()
            .map(|(id, _)| {
                let d = t.derivation(id).expect("every id has a derivation");
                Relation::Generated {
                    from: d.parents.clone(),
                    to: id,
                    role: d.role,
                }
            })
            .collect();
        History::new(
            self.operation,
            OperationKind::Revolve,
            Vec::new(),
            vec![t.body_id()],
            relations,
            Vec::new(),
        )
        .at_level(level)
    }

    /// A right circular cone or frustum on the frame's axis (S3 of
    /// REVIEW_NOTES.md): radius `bottom` at the frame origin and `top` at
    /// `height` along its normal, a zero radius being an apex, as
    /// `BRepPrimAPI_MakeCone(gp_Ax2, bottom, top, height)`. Its lateral face
    /// is a `Cone` surface whose seam OCCT would place along the frame's x.
    /// Mass properties are the general certified ones (REVIEW_NOTES.md U2).
    pub fn cone_with(
        operation: OperationId,
        frame: Frame3,
        bottom: f64,
        top: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::cone_in(
            &Context::new(operation),
            frame,
            bottom,
            top,
            height,
            tolerance,
        )
    }

    /// [`Solid::cone_with`] at a recorded algorithm level (H8).
    pub fn cone_at(
        level: AlgorithmLevel,
        operation: OperationId,
        frame: Frame3,
        bottom: f64,
        top: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::cone_in(
            &Context::new(operation).at(level),
            frame,
            bottom,
            top,
            height,
            tolerance,
        )
    }

    /// [`Solid::cone_with`] in an operation context. Every entity is
    /// generated from its meridian element; nothing carries an attribute.
    pub fn cone_in(
        context: &Context,
        frame: Frame3,
        bottom: f64,
        top: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let solid = Self::build_cone(operation, frame, bottom, top, height, tolerance)?;
        let history = solid.revolve_history(level);
        solid.debug_check(&[], &history);
        Ok((solid, history))
    }

    /// Split and fuse rebuild prisms of one profile.
    pub(crate) fn prism_profile(&self) -> Result<&Profile> {
        match &self.construction {
            Construction::Prism(profile) => Ok(profile),
            Construction::Cone { .. } => Err(Error::OutOfDomain(
                "split and fuse rebuild prisms; this solid is a cone",
            )),
            Construction::Sphere { .. } => Err(Error::OutOfDomain(
                "split and fuse rebuild prisms; this solid is a sphere",
            )),
            Construction::Clipped(_) | Construction::Half(_) => Err(Error::OutOfDomain(
                "a prism operation on a split piece (S8b)",
            )),
            Construction::Stack(_) => Err(Error::OutOfDomain(
                "a prism operation on a Boolean's stack (S9b)",
            )),
            Construction::Torus { .. } => Err(Error::OutOfDomain(
                "split and fuse rebuild prisms; this solid is a torus",
            )),
        }
    }

    /// Rectangle `[0, width] x [0, depth]` in the XY frame extruded to
    /// `height`. Convenience constructors derive ids from
    /// [`OperationId::UNSPECIFIED`] and discard the history; the `_with`
    /// forms take an operation id and return it.
    pub fn cuboid(width: f64, depth: f64, height: f64, tolerance: Tolerance) -> Result<Self> {
        Ok(Self::cuboid_with(OperationId::UNSPECIFIED, width, depth, height, tolerance)?.0)
    }

    pub fn cuboid_with(
        operation: OperationId,
        width: f64,
        depth: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        finite(height, "cuboid height")?;
        if height <= tolerance.linear() {
            return Err(Error::Degenerate("cuboid height"));
        }
        let profile = Profile::new(
            Boundary::rectangle(width, depth, tolerance)?,
            vec![],
            tolerance,
        )?;
        Self::extrude_with(operation, profile, Frame3::xy(), 0.0, height)
    }

    /// Axis-aligned box extending from `origin` by signed dimensions.
    /// Negative dimensions move the minimum corner, as in OCCT's
    /// `BRepPrimAPI_MakeBox(P, dx, dy, dz)`; volume remains positive.
    pub fn box_at(origin: Point3, size: crate::Vec3, tolerance: Tolerance) -> Result<Self> {
        Ok(Self::box_at_with(OperationId::UNSPECIFIED, origin, size, tolerance)?.0)
    }

    /// [`Solid::box_at`] with its history: the cuboid construction composed
    /// with the translation to the minimum corner.
    pub fn box_at_with(
        operation: OperationId,
        origin: Point3,
        size: crate::Vec3,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        origin.checked(tolerance)?;
        tolerance.resolve(&size.to_array())?;
        // Source: BRepPrimAPI_MakeBox.cxx, pmin and the point/size constructor.
        // See rust/SOURCE_MAP.md for the pinned revision and intentional limits.
        let minimum = origin + crate::Vec3::new(size.x.min(0.0), size.y.min(0.0), size.z.min(0.0));
        let (solid, built) = Self::cuboid_with(
            operation,
            size.x.abs(),
            size.y.abs(),
            size.z.abs(),
            tolerance,
        )?;
        let translation = RigidTransform::translation(minimum - Point3::ORIGIN)?;
        let (solid, moved) = solid.transform_with(operation, translation)?;
        let history = built
            .then(&moved)
            .map_err(|_| Error::InvalidTopology("history composition"))?;
        Ok((solid, history))
    }

    pub fn cylinder(
        frame: Frame3,
        radius: f64,
        start: f64,
        end: f64,
        tolerance: Tolerance,
    ) -> Result<Self> {
        Ok(Self::cylinder_with(
            OperationId::UNSPECIFIED,
            frame,
            radius,
            start,
            end,
            tolerance,
        )?
        .0)
    }

    pub fn cylinder_with(
        operation: OperationId,
        frame: Frame3,
        radius: f64,
        start: f64,
        end: f64,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        let profile = Profile::new(
            Boundary::circle(Point2::default(), radius, tolerance)?,
            vec![],
            tolerance,
        )?;
        Self::extrude_with(operation, profile, frame, start, end)
    }

    /// The extruded profile; `None` for a cone.
    pub fn profile(&self) -> Option<&Profile> {
        match &self.construction {
            Construction::Prism(profile) => Some(profile),
            Construction::Cone { .. }
            | Construction::Sphere { .. }
            | Construction::Torus { .. }
            | Construction::Clipped(_)
            | Construction::Half(_)
            | Construction::Stack(_) => None,
        }
    }
    /// The body's resolution.
    pub fn resolution(&self) -> Tolerance {
        match &self.construction {
            Construction::Prism(profile) => profile.tolerance(),
            Construction::Clipped(clipped) => clipped.tolerance(),
            Construction::Half(half) => half.tolerance(),
            Construction::Stack(stack) => stack.tolerance(),
            Construction::Cone { tolerance, .. }
            | Construction::Sphere { tolerance, .. }
            | Construction::Torus { tolerance, .. } => *tolerance,
        }
    }
    pub fn frame(&self) -> Frame3 {
        self.frame
    }
    pub fn start_offset(&self) -> f64 {
        self.start
    }
    pub fn end_offset(&self) -> f64 {
        self.end
    }
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    /// A watertight mesh of the solid within certified deflection and angle
    /// bounds (`tessellation::tessellate`).
    pub fn tessellate(
        &self,
        parameters: crate::tessellation::Parameters,
    ) -> Result<crate::tessellation::Mesh> {
        crate::tessellation::tessellate(&self.topology, parameters)
    }
    pub fn operation(&self) -> OperationId {
        self.operation
    }
    pub fn mass_properties(&self) -> MassProperties {
        self.mass
    }
    pub fn bounds(&self) -> Bounds3 {
        self.bounds
    }

    pub fn classify(&self, point: Point3) -> Result<Location> {
        let tolerance = self.resolution();
        point.checked(tolerance)?;
        let [x, y, z] = self.frame.coordinates(point);
        let profile = match &self.construction {
            Construction::Prism(profile) => profile,
            Construction::Clipped(clipped) => {
                let (low, high) = (self.start.min(self.end), self.start.max(self.end));
                return clipped.classify([x, y, z], low, high, tolerance);
            }
            Construction::Half(half) => return half.classify([x, y, z], tolerance),
            Construction::Stack(stack) => {
                return stack.classify([x, y, z], [self.start, self.end], tolerance)
            }
            Construction::Cone { bottom, top, .. } => {
                let z = finite(z, "axial coordinate")?;
                let local = [finite(x, "coordinate")?, finite(y, "coordinate")?, z];
                return Ok(
                    match decide::cone_location(local, *bottom, *top, self.end, tolerance.linear())
                    {
                        0 => Location::Inside,
                        1 => Location::Boundary,
                        _ => Location::Outside,
                    },
                );
            }
            Construction::Torus {
                major,
                minor,
                low,
                high,
                angle,
                ..
            } => {
                if high - low != std::f64::consts::TAU || *angle != std::f64::consts::TAU {
                    return Err(Error::OutOfDomain(
                        "classification of a torus segment or wedge",
                    ));
                }
                let local = [
                    finite(x, "coordinate")?,
                    finite(y, "coordinate")?,
                    finite(z, "axial coordinate")?,
                ];
                return Ok(
                    match decide::torus_location(local, *major, *minor, tolerance.linear()) {
                        0 => Location::Inside,
                        1 => Location::Boundary,
                        _ => Location::Outside,
                    },
                );
            }
            Construction::Sphere { radius, .. } => {
                let local = [
                    finite(x, "coordinate")?,
                    finite(y, "coordinate")?,
                    finite(z, "axial coordinate")?,
                ];
                return Ok(
                    match decide::sphere_location(
                        local,
                        *radius,
                        self.start,
                        self.end,
                        tolerance.linear(),
                    ) {
                        0 => Location::Inside,
                        1 => Location::Boundary,
                        _ => Location::Outside,
                    },
                );
            }
        };
        let z = finite(z, "axial coordinate")?;
        let (low, high) = (self.start.min(self.end), self.start.max(self.end));
        let tol = tolerance.linear();
        if decide::sum_gt(&[low, -tol], &[z]) || decide::sum_gt(&[z], &[high, tol]) {
            return Ok(Location::Outside);
        }
        match profile.classify(Point2::new(x, y))? {
            Location::Outside => Ok(Location::Outside),
            Location::Boundary => Ok(Location::Boundary),
            Location::Inside
                if [low, high].iter().any(|end| {
                    decide::sum_le(&[z, -end], &[tol]) && decide::sum_le(&[*end, -z], &[tol])
                }) =>
            {
                Ok(Location::Boundary)
            }
            Location::Inside => Ok(Location::Inside),
        }
    }

    /// [`Solid::transform_with`] without an operation id; the history is
    /// discarded.
    #[deprecated(note = "use Solid::transform_with, which also returns the history")]
    pub fn transformed(&self, transform: RigidTransform) -> Result<Self> {
        Ok(self.transform_with(OperationId::UNSPECIFIED, transform)?.0)
    }

    /// Rebuilds equivalent analytic geometry in the new frame, retaining
    /// body-local topology order, builder provenance, the body id and every
    /// entity id. The history reports every entity `Modified` with its id.
    pub fn transform_with(
        &self,
        operation: OperationId,
        transform: RigidTransform,
    ) -> Result<(Self, History)> {
        self.transform_at(AlgorithmLevel::CURRENT, operation, transform)
    }

    /// [`Solid::transform_with`] at a recorded algorithm level (H8).
    pub fn transform_at(
        &self,
        level: AlgorithmLevel,
        operation: OperationId,
        transform: RigidTransform,
    ) -> Result<(Self, History)> {
        self.transform_in(&Context::new(operation).at(level), transform)
    }

    /// [`Solid::transform_with`] in an operation context: its level, the
    /// attribute policies and the callback for `Recompute` policies.
    pub fn transform_in(
        &self,
        context: &Context,
        transform: RigidTransform,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let frame = self.frame.transformed(transform, self.resolution())?;
        let mut solid = match &self.construction {
            // A split half's certified mass is costly: its source's, moved.
            Construction::Half(half) => {
                half.rebuilt_with(self.operation, frame, Some(self.mass.moved(transform)))?
            }
            // A stack's too (S9a.2).
            Construction::Stack(stack) => {
                stack.rebuilt_with(self.operation, frame, Some(self.mass.moved(transform)))?
            }
            // An oblique piece's too (S8b's spline walls cost the most).
            Construction::Clipped(clipped) => clipped.rebuilt_with(
                self.operation,
                frame,
                self.start,
                self.end,
                Some(self.mass.moved(transform)),
            )?,
            _ => self.rebuilt(self.operation, frame)?,
        };
        // A rigid copy keeps every id, whatever operation named them.
        solid.topology = solid.topology.with_identity_of(&self.topology);
        let ids: Vec<EntityId> = self.topology.ids().map(|(id, _)| id).collect();
        let history = History::new(
            operation,
            OperationKind::Transform,
            vec![self.topology.body_id()],
            vec![solid.topology.body_id()],
            ids.into_iter()
                .map(|id| Relation::Modified { from: id, to: id })
                .collect(),
            Vec::new(),
        );
        let mut history = history.at_level(level);
        enclose::carry(&[self], &history, &mut [&mut solid]);
        let [moved] = attrs::propagate(context, &[self], &mut history, &[&solid])?
            .try_into()
            .expect("one output");
        solid.topology.set_attributes(moved);
        solid.debug_check(&[self.topology.entity_set(self.resolution())], &history);
        attrs::debug_check_attributes(context, &[self], &[&solid], &history);
        Ok((solid, history))
    }
}

pub(crate) fn replayable(level: AlgorithmLevel) -> Result<()> {
    if AlgorithmLevel::REPLAYABLE.contains(&level) {
        Ok(())
    } else {
        Err(Error::UnknownAlgorithmLevel(level.0))
    }
}

fn properties(profile: &Profile, frame: Frame3, low: f64, height: f64) -> Result<MassProperties> {
    let volume = finite(profile.area() * height, "volume")?;
    let surface_area = finite(
        2.0 * profile.area() + profile.perimeter() * height,
        "surface area",
    )?;
    let centroid = frame
        .point(profile.centroid(), low + height * 0.5)
        .checked(profile.tolerance())?;
    let [xx, xy, yy] = profile.moments.second;
    let axial = volume * height * height / 12.0;
    let local = [
        [height * yy + axial, -height * xy, 0.0],
        [-height * xy, height * xx + axial, 0.0],
        [0.0, 0.0, height * (xx + yy)],
    ];
    let columns = [
        frame.x().to_array(),
        frame.y().to_array(),
        frame.normal().to_array(),
    ];
    let mut inertia = [[0.0; 3]; 3];
    for (i, row) in inertia.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            for a in 0..3 {
                for b in 0..3 {
                    *entry += columns[a][i] * local[a][b] * columns[b][j];
                }
            }
            finite(*entry, "inertia")?;
        }
    }
    if volume <= 0.0 || surface_area <= 0.0 || (0..3).any(|i| inertia[i][i] <= 0.0) {
        return Err(Error::Degenerate("solid mass properties"));
    }
    Ok(MassProperties {
        volume,
        surface_area,
        centroid,
        inertia,
    })
}

fn bounds(profile: &Profile, frame: Frame3, low: f64, high: f64) -> Bounds3 {
    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    match &profile.outer().kind {
        BoundaryKind::Polygon(points) => {
            for height in [low, high] {
                for point in points {
                    let world = frame.point(*point, height).to_array();
                    for i in 0..3 {
                        min[i] = min[i].min(world[i]);
                        max[i] = max[i].max(world[i]);
                    }
                }
            }
        }
        BoundaryKind::Path { points, segments } => {
            // Points, and each arc's whole circle (conservative).
            for height in [low, high] {
                for point in points {
                    let world = frame.point(*point, height).to_array();
                    for i in 0..3 {
                        min[i] = min[i].min(world[i]);
                        max[i] = max[i].max(world[i]);
                    }
                }
            }
            let (x, y) = (frame.x().to_array(), frame.y().to_array());
            for segment in segments {
                // S8b: a spline lies in its control points' hull.
                if let crate::profile::Segment::Spline(span) = segment {
                    for height in [low, high] {
                        for p in span.curve().poles() {
                            let world = frame.point(p, height).to_array();
                            for i in 0..3 {
                                min[i] = min[i].min(world[i]);
                                max[i] = max[i].max(world[i]);
                            }
                        }
                    }
                }
                if let crate::profile::Segment::Arc { center, radius, .. } = segment {
                    for height in [low, high] {
                        let c = frame.point(*center, height).to_array();
                        for i in 0..3 {
                            let extent = radius * x[i].hypot(y[i]);
                            min[i] = min[i].min(c[i] - extent);
                            max[i] = max[i].max(c[i] + extent);
                        }
                    }
                }
            }
        }
        BoundaryKind::Circle { center, radius } => {
            let bottom = frame.point(*center, low).to_array();
            let top = frame.point(*center, high).to_array();
            let (x, y) = (frame.x().to_array(), frame.y().to_array());
            for i in 0..3 {
                let extent = radius * x[i].hypot(y[i]);
                min[i] = bottom[i].min(top[i]) - extent;
                max[i] = bottom[i].max(top[i]) + extent;
            }
        }
    }
    Bounds3 {
        min: Point3::new(min[0], min[1], min[2]),
        max: Point3::new(max[0], max[1], max[2]),
    }
}
