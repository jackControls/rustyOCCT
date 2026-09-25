use crate::history::{self, History, Relation};
use crate::identity::{EntityId, OperationId, OperationKind};
use crate::math::finite;
use crate::profile::BoundaryKind;
use crate::topology::Topology;
use crate::{
    Boundary, Bounds3, Error, Frame3, Location, Point2, Point3, Profile, Result, RigidTransform,
    Tolerance,
};

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

/// An immutable, validated normal extrusion of one planar material region.
/// The retained construction supports exact queries without a triangle mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Solid {
    profile: Profile,
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
        solid.debug_check(&[], &history);
        Ok((solid, history))
    }

    /// Every operation checks its history independently in debug builds.
    fn debug_check(&self, inputs: &[history::EntitySet], history: &History) {
        if cfg!(debug_assertions) {
            let outputs = [self.topology.entity_set(self.profile.tolerance())];
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
        if height <= tolerance.linear() {
            return Err(Error::Degenerate("extrusion height"));
        }
        let bounds = bounds(&profile, frame, low, high);
        bounds.min.checked(tolerance)?;
        bounds.max.checked(tolerance)?;
        let mass = properties(&profile, frame, low, height)?;
        let topology = Topology::prism(&profile, frame, low, high, start < end, operation)?;
        Ok(Self {
            profile,
            frame,
            start,
            end,
            topology,
            mass,
            bounds,
            operation,
        })
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

    pub fn profile(&self) -> &Profile {
        &self.profile
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
        let tolerance = self.profile.tolerance();
        point.checked(tolerance)?;
        let [x, y, z] = self.frame.coordinates(point);
        let (low, high) = (self.start.min(self.end), self.start.max(self.end));
        if z < low - tolerance.linear() || z > high + tolerance.linear() {
            return Ok(Location::Outside);
        }
        match self.profile.classify(Point2::new(x, y))? {
            Location::Outside => Ok(Location::Outside),
            Location::Boundary => Ok(Location::Boundary),
            Location::Inside
                if (z - low).abs() <= tolerance.linear()
                    || (z - high).abs() <= tolerance.linear() =>
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
        let solid = Self::build(
            self.operation,
            self.profile.clone(),
            self.frame
                .transformed(transform, self.profile.tolerance())?,
            self.start,
            self.end,
        )?;
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
        solid.debug_check(
            &[self.topology.entity_set(self.profile.tolerance())],
            &history,
        );
        Ok((solid, history))
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
