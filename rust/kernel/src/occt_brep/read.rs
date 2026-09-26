//! The reader: every record of OCCT's `.brep` text format, versions 1 to 3
//! (`dox/specification/brep_format.md`, mirrored against
//! `BRepTools_ShapeSet` and `TopTools_LocationSet` at the pinned revision).
//! Records the kernel cannot represent are read in full and kept by name, so
//! that nothing after them is misread.
use super::{BrepError, Transform};

/// A 3D curve record.
#[derive(Debug, Clone, PartialEq)]
pub enum Curve3 {
    /// `P + u D`.
    Line {
        p: [f64; 3],
        d: [f64; 3],
    },
    /// `P + r (cos u X + sin u Y)`; `n` is the record's main direction.
    Circle {
        p: [f64; 3],
        n: [f64; 3],
        x: [f64; 3],
        y: [f64; 3],
        r: f64,
    },
    Other(&'static str),
}

/// A 2D curve record.
#[derive(Debug, Clone, PartialEq)]
pub enum Curve2 {
    Line {
        p: [f64; 2],
        d: [f64; 2],
    },
    Circle {
        c: [f64; 2],
        x: [f64; 2],
        y: [f64; 2],
        r: f64,
    },
    Other(&'static str),
}

/// A surface record.
#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    Plane {
        p: [f64; 3],
        n: [f64; 3],
        x: [f64; 3],
        y: [f64; 3],
    },
    Cylinder {
        p: [f64; 3],
        n: [f64; 3],
        x: [f64; 3],
        y: [f64; 3],
        r: f64,
    },
    /// `Geom_ConicalSurface`: reference radius `r` in the plane through `p`,
    /// semi-angle `a`.
    Cone {
        p: [f64; 3],
        n: [f64; 3],
        x: [f64; 3],
        y: [f64; 3],
        r: f64,
        a: f64,
    },
    Other(&'static str),
}

/// One representation of an edge.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgeRep {
    /// 3D curve index (1-based), location index, parameter range.
    Curve {
        curve: usize,
        location: usize,
        range: [f64; 2],
    },
    /// A pcurve (two on a closed surface: the first for the edge's forward
    /// use) on surface `surface` with location `location`.
    OnSurface {
        pcurves: Vec<usize>,
        surface: usize,
        location: usize,
        range: [f64; 2],
    },
    /// Regularity, polygons and triangulation data: not geometry the
    /// kernel uses.
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
}

/// `+`, `-`, `i` or `e`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orient {
    Forward,
    Reversed,
    Internal,
    External,
}

/// A reference to a shape record (0-based index in file order) with its
/// orientation and location index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sub {
    pub orient: Orient,
    pub shape: usize,
    pub location: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Data {
    Vertex {
        tolerance: f64,
        point: [f64; 3],
    },
    Edge {
        tolerance: f64,
        degenerated: bool,
        reps: Vec<EdgeRep>,
    },
    Face {
        tolerance: f64,
        surface: usize,
        location: usize,
    },
    /// A face given only by a triangulation.
    MeshFace,
    None,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Shape {
    pub kind: Kind,
    pub data: Data,
    pub subs: Vec<Sub>,
}

/// A whole `.brep` document.
#[derive(Debug, Clone, PartialEq)]
pub struct Document {
    pub version: u8,
    /// Composed locations, index `i` at `locations[i - 1]`; index 0 is the
    /// identity.
    pub locations: Vec<Transform>,
    pub curves2d: Vec<Curve2>,
    pub curves: Vec<Curve3>,
    pub surfaces: Vec<Surface>,
    pub triangulations: usize,
    /// Shape records in file order; every `Sub` holds a 0-based index into
    /// this list (the file numbers records backwards from its count).
    pub shapes: Vec<Shape>,
    pub root: Sub,
}

struct Tokens<'a> {
    list: Vec<(usize, &'a str)>,
    at: usize,
}

impl<'a> Tokens<'a> {
    fn new(text: &'a str) -> Self {
        let list = text
            .lines()
            .enumerate()
            .flat_map(|(i, l)| l.split_whitespace().map(move |w| (i + 1, w)))
            .collect();
        Self { list, at: 0 }
    }
    fn line(&self) -> usize {
        self.list
            .get(self.at.min(self.list.len().saturating_sub(1)))
            .map_or(0, |t| t.0)
    }
    fn err(&self, what: &'static str) -> BrepError {
        BrepError::Syntax {
            line: self.line(),
            what,
        }
    }
    fn next(&mut self, what: &'static str) -> Result<&'a str, BrepError> {
        let t = self.list.get(self.at).ok_or(self.err(what))?.1;
        self.at += 1;
        Ok(t)
    }
    fn peek(&self) -> Option<&'a str> {
        self.list.get(self.at).map(|t| t.1)
    }
    fn word(&mut self, want: &'static str) -> Result<(), BrepError> {
        if self.next(want)? == want {
            Ok(())
        } else {
            self.at -= 1;
            Err(self.err(want))
        }
    }
    fn real(&mut self) -> Result<f64, BrepError> {
        let t = self.next("a real number")?;
        let v: f64 = t.parse().map_err(|_| {
            self.at -= 1;
            self.err("a real number")
        })?;
        if !v.is_finite() {
            self.at -= 1;
            return Err(self.err("a finite real number"));
        }
        Ok(v)
    }
    fn int(&mut self) -> Result<i64, BrepError> {
        let t = self.next("an integer")?;
        t.parse().map_err(|_| {
            self.at -= 1;
            self.err("an integer")
        })
    }
    /// A count or index, bounded so that malformed input cannot demand an
    /// unbounded allocation.
    fn count(&mut self, limit: usize) -> Result<usize, BrepError> {
        let v = self.int()?;
        if v < 0 || v as u64 > limit as u64 {
            self.at -= 1;
            return Err(self.err("a count within the supported limit"));
        }
        Ok(v as usize)
    }
    fn flag(&mut self) -> Result<bool, BrepError> {
        match self.int()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => {
                self.at -= 1;
                Err(self.err("a 0/1 flag"))
            }
        }
    }
    fn reals<const N: usize>(&mut self) -> Result<[f64; N], BrepError> {
        let mut out = [0.0; N];
        for v in &mut out {
            *v = self.real()?;
        }
        Ok(out)
    }
    fn skip_reals(&mut self, n: usize) -> Result<(), BrepError> {
        for _ in 0..n {
            self.real()?;
        }
        Ok(())
    }
}

/// Records and poles are bounded so malformed counts fail rather than
/// allocate.
const LIMIT: usize = 1 << 24;

fn curve3(t: &mut Tokens) -> Result<Curve3, BrepError> {
    Ok(match t.int()? {
        1 => {
            let [a, b, c, d, e, f] = t.reals::<6>()?;
            Curve3::Line {
                p: [a, b, c],
                d: [d, e, f],
            }
        }
        2 => {
            let v = t.reals::<13>()?;
            Curve3::Circle {
                p: [v[0], v[1], v[2]],
                n: [v[3], v[4], v[5]],
                x: [v[6], v[7], v[8]],
                y: [v[9], v[10], v[11]],
                r: v[12],
            }
        }
        3 => {
            t.skip_reals(14)?;
            Curve3::Other("Ellipse")
        }
        4 => {
            t.skip_reals(13)?;
            Curve3::Other("Parabola")
        }
        5 => {
            t.skip_reals(14)?;
            Curve3::Other("Hyperbola")
        }
        6 => {
            let rational = t.flag()?;
            let degree = t.count(64)?;
            t.skip_reals((degree + 1) * if rational { 4 } else { 3 })?;
            Curve3::Other("BezierCurve")
        }
        7 => {
            // Rational and periodic flags (the specification's example shows
            // a periodic flag of 0; it is a flag).
            let rational = t.flag()?;
            t.flag()?;
            t.count(64)?;
            let poles = t.count(LIMIT)?;
            let knots = t.count(LIMIT)?;
            t.skip_reals(poles * if rational { 4 } else { 3 })?;
            for _ in 0..knots {
                t.real()?;
                t.int()?;
            }
            Curve3::Other("BSplineCurve")
        }
        8 => {
            // A trimmed line or circle is its basis: an edge's own range
            // bounds the part it uses.
            t.skip_reals(2)?;
            match curve3(t)? {
                c @ (Curve3::Line { .. } | Curve3::Circle { .. }) => c,
                _ => Curve3::Other("TrimmedCurve"),
            }
        }
        9 => {
            t.skip_reals(4)?;
            curve3(t)?;
            Curve3::Other("OffsetCurve")
        }
        _ => {
            t.at -= 1;
            return Err(t.err("a 3D curve type 1-9"));
        }
    })
}

fn curve2(t: &mut Tokens) -> Result<Curve2, BrepError> {
    Ok(match t.int()? {
        1 => {
            let [a, b, c, d] = t.reals::<4>()?;
            Curve2::Line {
                p: [a, b],
                d: [c, d],
            }
        }
        2 => {
            let v = t.reals::<7>()?;
            Curve2::Circle {
                c: [v[0], v[1]],
                x: [v[2], v[3]],
                y: [v[4], v[5]],
                r: v[6],
            }
        }
        3 => {
            t.skip_reals(8)?;
            Curve2::Other("Ellipse2d")
        }
        4 => {
            t.skip_reals(7)?;
            Curve2::Other("Parabola2d")
        }
        5 => {
            t.skip_reals(8)?;
            Curve2::Other("Hyperbola2d")
        }
        6 => {
            let rational = t.flag()?;
            let degree = t.count(64)?;
            t.skip_reals((degree + 1) * if rational { 3 } else { 2 })?;
            Curve2::Other("BezierCurve2d")
        }
        7 => {
            let rational = t.flag()?;
            t.flag()?;
            t.count(64)?;
            let poles = t.count(LIMIT)?;
            let knots = t.count(LIMIT)?;
            t.skip_reals(poles * if rational { 3 } else { 2 })?;
            for _ in 0..knots {
                t.real()?;
                t.int()?;
            }
            Curve2::Other("BSplineCurve2d")
        }
        8 => {
            t.skip_reals(2)?;
            match curve2(t)? {
                c @ (Curve2::Line { .. } | Curve2::Circle { .. }) => c,
                _ => Curve2::Other("TrimmedCurve2d"),
            }
        }
        9 => {
            t.skip_reals(1)?;
            curve2(t)?;
            Curve2::Other("OffsetCurve2d")
        }
        _ => {
            t.at -= 1;
            return Err(t.err("a 2D curve type 1-9"));
        }
    })
}

fn surface(t: &mut Tokens) -> Result<Surface, BrepError> {
    Ok(match t.int()? {
        1 => {
            let v = t.reals::<12>()?;
            Surface::Plane {
                p: [v[0], v[1], v[2]],
                n: [v[3], v[4], v[5]],
                x: [v[6], v[7], v[8]],
                y: [v[9], v[10], v[11]],
            }
        }
        2 => {
            let v = t.reals::<13>()?;
            Surface::Cylinder {
                p: [v[0], v[1], v[2]],
                n: [v[3], v[4], v[5]],
                x: [v[6], v[7], v[8]],
                y: [v[9], v[10], v[11]],
                r: v[12],
            }
        }
        3 => {
            let v = t.reals::<14>()?;
            Surface::Cone {
                p: [v[0], v[1], v[2]],
                n: [v[3], v[4], v[5]],
                x: [v[6], v[7], v[8]],
                y: [v[9], v[10], v[11]],
                r: v[12],
                a: v[13],
            }
        }
        4 => {
            t.skip_reals(13)?;
            Surface::Other("SphericalSurface")
        }
        5 => {
            t.skip_reals(14)?;
            Surface::Other("ToroidalSurface")
        }
        6 => {
            t.skip_reals(3)?;
            curve3(t)?;
            Surface::Other("SurfaceOfLinearExtrusion")
        }
        7 => {
            t.skip_reals(6)?;
            curve3(t)?;
            Surface::Other("SurfaceOfRevolution")
        }
        8 => {
            let (ru, rv) = (t.flag()?, t.flag()?);
            let (du, dv) = (t.count(64)?, t.count(64)?);
            t.skip_reals((du + 1) * (dv + 1) * if ru || rv { 4 } else { 3 })?;
            Surface::Other("BezierSurface")
        }
        9 => {
            let (ru, rv) = (t.flag()?, t.flag()?);
            // Periodic flags in u and v.
            t.flag()?;
            t.flag()?;
            t.count(64)?;
            t.count(64)?;
            let (pu, pv) = (t.count(LIMIT)?, t.count(LIMIT)?);
            let (ku, kv) = (t.count(LIMIT)?, t.count(LIMIT)?);
            let poles = pu
                .checked_mul(pv)
                .filter(|n| *n <= LIMIT)
                .ok_or(t.err("a pole count within the supported limit"))?;
            t.skip_reals(poles * if ru || rv { 4 } else { 3 })?;
            for _ in 0..ku + kv {
                t.real()?;
                t.int()?;
            }
            Surface::Other("BSplineSurface")
        }
        10 => {
            t.skip_reals(4)?;
            surface(t)?;
            Surface::Other("RectangularTrimmedSurface")
        }
        11 => {
            t.skip_reals(1)?;
            surface(t)?;
            Surface::Other("OffsetSurface")
        }
        _ => {
            t.at -= 1;
            return Err(t.err("a surface type 1-11"));
        }
    })
}

fn continuity(word: &str) -> bool {
    matches!(word, "C0" | "C1" | "C2" | "C3" | "CN" | "G1" | "G2")
}

/// An index and the continuity after it. OCCT's stream reader stops the
/// number at the first letter, so version 1 files write them glued
/// together (`6CN`); a separate continuity token is read as well.
fn index_and_regularity(t: &mut Tokens) -> Result<usize, BrepError> {
    let w = t.next("an index and a continuity")?;
    let digits = w.bytes().take_while(u8::is_ascii_digit).count();
    let (number, rest) = w.split_at(digits);
    let index = number
        .parse::<usize>()
        .ok()
        .filter(|i| *i <= LIMIT)
        .ok_or(t.err("an index"))?;
    if rest.is_empty() {
        regularity(t)?;
    } else if !continuity(rest) {
        return Err(t.err("a continuity C0-C3, CN, G1 or G2"));
    }
    Ok(index)
}

fn regularity(t: &mut Tokens) -> Result<(), BrepError> {
    match t.next("a continuity")? {
        "C0" | "C1" | "C2" | "C3" | "CN" | "G1" | "G2" => Ok(()),
        _ => {
            t.at -= 1;
            Err(t.err("a continuity C0-C3, CN, G1 or G2"))
        }
    }
}

fn sub(t: &mut Tokens, word: &str) -> Result<Sub, BrepError> {
    let mut chars = word.chars();
    let orient = match chars.next() {
        Some('+') => Orient::Forward,
        Some('-') => Orient::Reversed,
        Some('i') => Orient::Internal,
        Some('e') => Orient::External,
        _ => return Err(t.err("a subshape orientation")),
    };
    let shape = chars
        .as_str()
        .parse::<usize>()
        .map_err(|_| t.err("a subshape number"))?;
    let location = t.count(LIMIT)?;
    Ok(Sub {
        orient,
        shape,
        location,
    })
}

/// The root reference. OCCT reads its location with `operator>>` and reads
/// nothing after it, so a location glued to trailing bytes (`0c`, seen in
/// upstream test data) counts by its leading digits, as there.
fn root_sub(t: &mut Tokens, word: &str) -> Result<Sub, BrepError> {
    let mut chars = word.chars();
    let orient = match chars.next() {
        Some('+') => Orient::Forward,
        Some('-') => Orient::Reversed,
        Some('i') => Orient::Internal,
        Some('e') => Orient::External,
        _ => return Err(t.err("a subshape orientation")),
    };
    let shape = chars
        .as_str()
        .parse::<usize>()
        .map_err(|_| t.err("a subshape number"))?;
    let location = t.next("a location index")?;
    let digits: String = location
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    let location = digits
        .parse::<usize>()
        .ok()
        .filter(|l| *l <= LIMIT)
        .ok_or_else(|| t.err("a location index"))?;
    Ok(Sub {
        orient,
        shape,
        location,
    })
}

fn section(t: &mut Tokens, name: &'static str) -> Result<usize, BrepError> {
    t.word(name)?;
    t.count(LIMIT)
}

/// Read a whole document. Every malformed input is a typed error with its
/// line; nothing panics and no count allocates beyond the supported limit.
pub fn read(text: &str) -> Result<Document, BrepError> {
    let mut t = Tokens::new(text);
    // Header: optional DBRep line, then "CASCADE Topology Vn, (c) Matra-Datavision".
    let mut version = 0;
    while let Some(w) = t.peek() {
        t.at += 1;
        if w == "Topology" {
            version = match t.next("a format version")? {
                "V1," => 1,
                "V2," => 2,
                "V3," => 3,
                _ => return Err(t.err("format version V1, V2 or V3")),
            };
            t.word("(c)")?;
            t.word("Matra-Datavision")?;
            break;
        }
        if t.at > 8 {
            return Err(t.err("a CASCADE Topology header"));
        }
    }
    if version == 0 {
        return Err(t.err("a CASCADE Topology header"));
    }
    // Locations (TopTools_LocationSet::Read): a matrix record always gets
    // the next number, even the identity matrix; a composite record whose
    // chain is empty is the identity location and gets none.
    let n = section(&mut t, "Locations")?;
    let mut locations: Vec<Transform> = Vec::new();
    for _ in 0..n {
        let (l, empty) = match t.int()? {
            1 => (Transform(t.reals::<12>()?), false),
            2 => {
                let mut l = Transform::IDENTITY;
                let mut empty = true;
                loop {
                    let index = t.count(LIMIT)?;
                    if index == 0 {
                        break;
                    }
                    empty = false;
                    let power = t.int()?;
                    let base = *locations
                        .get(index - 1)
                        .ok_or(t.err("a location defined earlier"))?;
                    l = base.powered(power).map_err(|what| t.err(what))?.times(&l);
                }
                (l, empty)
            }
            _ => {
                t.at -= 1;
                return Err(t.err("a location record type 1 or 2"));
            }
        };
        if !empty {
            locations.push(l);
        }
    }
    let n = section(&mut t, "Curve2ds")?;
    let curves2d: Vec<Curve2> = (0..n).map(|_| curve2(&mut t)).collect::<Result<_, _>>()?;
    let n = section(&mut t, "Curves")?;
    let curves: Vec<Curve3> = (0..n).map(|_| curve3(&mut t)).collect::<Result<_, _>>()?;
    let n = section(&mut t, "Polygon3D")?;
    for _ in 0..n {
        let nodes = t.count(LIMIT)?;
        let params = t.flag()?;
        t.real()?;
        t.skip_reals(nodes * 3)?;
        if params {
            t.skip_reals(nodes)?;
        }
    }
    let n = section(&mut t, "PolygonOnTriangulations")?;
    for _ in 0..n {
        let nodes = t.count(LIMIT)?;
        for _ in 0..nodes {
            t.int()?;
        }
        t.word("p")?;
        t.real()?;
        if t.flag()? {
            t.skip_reals(nodes)?;
        }
    }
    let n = section(&mut t, "Surfaces")?;
    let surfaces: Vec<Surface> = (0..n).map(|_| surface(&mut t)).collect::<Result<_, _>>()?;
    let triangulations = section(&mut t, "Triangulations")?;
    for _ in 0..triangulations {
        let (nodes, triangles) = (t.count(LIMIT)?, t.count(LIMIT)?);
        let uv = t.flag()?;
        let normals = version >= 3 && t.flag()?;
        t.real()?;
        t.skip_reals(nodes * 3)?;
        if uv {
            t.skip_reals(nodes * 2)?;
        }
        for _ in 0..triangles * 3 {
            t.int()?;
        }
        if normals {
            t.skip_reals(nodes * 3)?;
        }
    }
    let n = section(&mut t, "TShapes")?;
    let mut shapes = Vec::with_capacity(n.min(1 << 16));
    for _ in 0..n {
        let kind = match t.next("a shape type")? {
            "Ve" => Kind::Vertex,
            "Ed" => Kind::Edge,
            "Wi" => Kind::Wire,
            "Fa" => Kind::Face,
            "Sh" => Kind::Shell,
            "So" => Kind::Solid,
            "CS" => Kind::CompSolid,
            "Co" => Kind::Compound,
            _ => {
                t.at -= 1;
                return Err(t.err("a shape type Ve, Ed, Wi, Fa, Sh, So, CS or Co"));
            }
        };
        let data = match kind {
            Kind::Vertex => {
                let tolerance = t.real()?;
                let point = t.reals::<3>()?;
                loop {
                    t.real()?;
                    match t.int()? {
                        0 => break,
                        1 => {
                            t.count(LIMIT)?;
                        }
                        2 => {
                            t.count(LIMIT)?;
                            t.count(LIMIT)?;
                        }
                        3 => {
                            t.real()?;
                            t.count(LIMIT)?;
                        }
                        _ => {
                            t.at -= 1;
                            return Err(t.err("a vertex representation type 0-3"));
                        }
                    }
                    t.count(LIMIT)?;
                }
                Data::Vertex { tolerance, point }
            }
            Kind::Edge => {
                let tolerance = t.real()?;
                t.flag()?;
                t.flag()?;
                let degenerated = t.flag()?;
                let mut reps = Vec::new();
                loop {
                    let rep = match t.int()? {
                        0 => break,
                        1 => {
                            let (curve, location) = (t.count(LIMIT)?, t.count(LIMIT)?);
                            EdgeRep::Curve {
                                curve,
                                location,
                                range: t.reals::<2>()?,
                            }
                        }
                        k @ (2 | 3) => {
                            let mut pcurves = vec![t.count(LIMIT)?];
                            if k == 3 {
                                pcurves.push(index_and_regularity(&mut t)?);
                            }
                            let (surface, location) = (t.count(LIMIT)?, t.count(LIMIT)?);
                            let range = t.reals::<2>()?;
                            if version == 2 {
                                t.skip_reals(4)?;
                            }
                            EdgeRep::OnSurface {
                                pcurves,
                                surface,
                                location,
                                range,
                            }
                        }
                        4 => {
                            regularity(&mut t)?;
                            for _ in 0..4 {
                                t.count(LIMIT)?;
                            }
                            EdgeRep::Other
                        }
                        5 => {
                            t.count(LIMIT)?;
                            t.count(LIMIT)?;
                            EdgeRep::Other
                        }
                        k @ (6 | 7) => {
                            for _ in 0..if k == 7 { 4 } else { 3 } {
                                t.count(LIMIT)?;
                            }
                            EdgeRep::Other
                        }
                        _ => {
                            t.at -= 1;
                            return Err(t.err("an edge representation type 0-7"));
                        }
                    };
                    reps.push(rep);
                }
                Data::Edge {
                    tolerance,
                    degenerated,
                    reps,
                }
            }
            Kind::Face => match t.int()? {
                0 | 1 => {
                    let tolerance = t.real()?;
                    let (surface, location) = (t.count(LIMIT)?, t.count(LIMIT)?);
                    // An optional "2 <triangulation>" line precedes the flags.
                    if t.peek() == Some("2") {
                        t.at += 1;
                        t.count(LIMIT)?;
                    }
                    Data::Face {
                        tolerance,
                        surface,
                        location,
                    }
                }
                2 => {
                    t.count(LIMIT)?;
                    Data::MeshFace
                }
                _ => {
                    t.at -= 1;
                    return Err(t.err("a face natural-restriction flag"));
                }
            },
            _ => Data::None,
        };
        let flags = t.next("a shape flag word")?;
        if flags.len() != 7 || !flags.bytes().all(|b| b == b'0' || b == b'1') {
            t.at -= 1;
            return Err(t.err("a seven-flag word"));
        }
        let mut subs = Vec::new();
        loop {
            let w = t.next("a subshape or *")?;
            if w == "*" {
                break;
            }
            let s = sub(&mut t, w)?;
            // Records are numbered backwards: the first is `n`, the last 1.
            // A subshape refers to a record read earlier.
            if s.shape == 0 || s.shape > n || n - s.shape >= shapes.len() {
                return Err(t.err("a subshape read earlier"));
            }
            subs.push(Sub {
                shape: n - s.shape,
                ..s
            });
        }
        shapes.push(Shape { kind, data, subs });
    }
    let w = t.next("the root shape")?;
    let root = root_sub(&mut t, w)?;
    if root.shape == 0 || root.shape > n {
        return Err(t.err("a root shape in the list"));
    }
    let root = Sub {
        shape: n - root.shape,
        ..root
    };
    if root.location > locations.len() {
        return Err(BrepError::Reference {
            what: "a location index beyond the table",
        });
    }
    // Every table reference is in range (0 is "none" where the format
    // allows it); the converter indexes the tables without checks.
    let beyond = |what| Err(BrepError::Reference { what });
    for s in &shapes {
        if s.subs.iter().any(|sub| sub.location > locations.len()) {
            return beyond("a location index beyond the table");
        }
        match &s.data {
            Data::Edge { reps, .. } => {
                for rep in reps {
                    match rep {
                        EdgeRep::Curve {
                            curve, location, ..
                        } => {
                            if *curve > curves.len() {
                                return beyond("a curve index beyond the table");
                            }
                            if *location > locations.len() {
                                return beyond("a location index beyond the table");
                            }
                        }
                        EdgeRep::OnSurface {
                            pcurves,
                            surface,
                            location,
                            ..
                        } => {
                            if pcurves.iter().any(|c| *c > curves2d.len()) {
                                return beyond("a pcurve index beyond the table");
                            }
                            if *surface > surfaces.len() {
                                return beyond("a surface index beyond the table");
                            }
                            if *location > locations.len() {
                                return beyond("a location index beyond the table");
                            }
                        }
                        EdgeRep::Other => {}
                    }
                }
            }
            Data::Face {
                surface, location, ..
            } => {
                if *surface > surfaces.len() {
                    return beyond("a surface index beyond the table");
                }
                if *location > locations.len() {
                    return beyond("a location index beyond the table");
                }
            }
            _ => {}
        }
    }
    Ok(Document {
        version,
        locations,
        curves2d,
        curves,
        surfaces,
        triangulations,
        shapes,
        root,
    })
}
