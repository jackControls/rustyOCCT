use crate::identity::OperationId;
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
    /// Either direction is allowed; zero or sub-tolerance thickness is rejected.
    /// Ids derive from [`OperationId::UNSPECIFIED`]; see [`Solid::extrude_with`].
    pub fn extrude(profile: Profile, frame: Frame3, start: f64, end: f64) -> Result<Self> {
        Self::extrude_with(OperationId::UNSPECIFIED, profile, frame, start, end)
    }

    /// [`Solid::extrude`] with a caller operation id, the root of every entity
    /// id together with the profile's labels (or indices).
    pub fn extrude_with(
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

    pub fn cuboid(width: f64, depth: f64, height: f64, tolerance: Tolerance) -> Result<Self> {
        finite(height, "cuboid height")?;
        if height <= tolerance.linear() {
            return Err(Error::Degenerate("cuboid height"));
        }
        let profile = Profile::new(
            Boundary::rectangle(width, depth, tolerance)?,
            vec![],
            tolerance,
        )?;
        Self::extrude(profile, Frame3::xy(), 0.0, height)
    }

    /// Axis-aligned box extending from `origin` by signed dimensions.
    /// Negative dimensions move the minimum corner, as in OCCT's
    /// `BRepPrimAPI_MakeBox(P, dx, dy, dz)`; volume remains positive.
    pub fn box_at(origin: Point3, size: crate::Vec3, tolerance: Tolerance) -> Result<Self> {
        origin.checked(tolerance)?;
        tolerance.resolve(&size.to_array())?;
        // Source: BRepPrimAPI_MakeBox.cxx, pmin and the point/size constructor.
        // See rust/SOURCE_MAP.md for the pinned revision and intentional limits.
        let minimum = origin + crate::Vec3::new(size.x.min(0.0), size.y.min(0.0), size.z.min(0.0));
        let solid = Self::cuboid(size.x.abs(), size.y.abs(), size.z.abs(), tolerance)?;
        solid.transformed(RigidTransform::translation(minimum - Point3::ORIGIN)?)
    }

    pub fn cylinder(
        frame: Frame3,
        radius: f64,
        start: f64,
        end: f64,
        tolerance: Tolerance,
    ) -> Result<Self> {
        let profile = Profile::new(
            Boundary::circle(Point2::default(), radius, tolerance)?,
            vec![],
            tolerance,
        )?;
        Self::extrude(profile, frame, start, end)
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

    /// Rebuilds equivalent analytic geometry in the new frame, retaining
    /// body-local topology order, builder provenance and every id.
    pub fn transformed(&self, transform: RigidTransform) -> Result<Self> {
        Self::extrude_with(
            self.operation,
            self.profile.clone(),
            self.frame
                .transformed(transform, self.profile.tolerance())?,
            self.start,
            self.end,
        )
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
