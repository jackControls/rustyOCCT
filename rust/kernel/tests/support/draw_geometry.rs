//! S7 (the `lowalgos` group): the curves DRAW's `intersect` names, their
//! bounds and points, distances to the analytic surfaces for `xdistcs`, and
//! DRAW's numeric expressions (`dval`). Test-only, used by `draw_worker`.
use rusty_occt::intersection::{
    surface_surface, AnalyticItem, Branch, Component, ProceduralCurve, SurfaceIntersection,
    TracedCurve,
};
use rusty_occt::topology::Surface;
use rusty_occt::{Frame3, Point3, Tolerance, Vec3};
use std::f64::consts::TAU;
use std::rc::Rc;

/// OCCT's `Precision::Infinite()`: an unbounded curve's DRAW bounds.
pub const INFINITE: f64 = 2e100;

/// A DRAW 3D curve from `intersect`.
#[derive(Clone)]
pub enum Curve {
    /// Through a point along a unit direction.
    Line(Point3, Vec3),
    Circle(Frame3, f64),
    Ellipse(Frame3, f64, f64),
    /// One branch, `o + a cosh t x + b sinh t y`.
    Hyperbola(Frame3, f64, f64),
    /// A procedural loop: `tau` in `[0, 1]` runs the plus branch from `u0`
    /// to `u1`, `[1, 2]` the minus branch back.
    Loop {
        curve: Rc<ProceduralCurve>,
        u: [f64; 2],
    },
    /// A procedural ring (one branch) or figure-eight (both, one after the
    /// other) over whole turns from `start`: `tau` in `[0, 2 pi n]`.
    Turn {
        curve: Rc<ProceduralCurve>,
        start: f64,
        branches: Vec<Branch>,
    },
    /// One track of a traced curve, over its `phi` range.
    Track {
        curve: Rc<TracedCurve>,
        track: usize,
        phi: [f64; 2],
    },
}

fn mid([lo, hi]: [f64; 2]) -> f64 {
    0.5 * lo + 0.5 * hi
}

fn mid3(e: [[f64; 2]; 3]) -> Point3 {
    Point3::new(mid(e[0]), mid(e[1]), mid(e[2]))
}

fn frame(centre: Point3, normal: Vec3, x: Option<Vec3>) -> Result<Frame3, String> {
    let hint = x.unwrap_or_else(|| {
        if normal.x.abs() < 0.5 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        }
    });
    Frame3::new(centre, normal, hint, Tolerance::default()).map_err(|e| e.to_string())
}

impl Curve {
    /// DRAW's `bounds`.
    pub fn bounds(&self) -> [f64; 2] {
        match self {
            Curve::Line(..) | Curve::Hyperbola(..) => [-INFINITE, INFINITE],
            Curve::Circle(..) | Curve::Ellipse(..) => [0.0, TAU],
            Curve::Loop { .. } => [0.0, 2.0],
            Curve::Turn { branches, .. } => [0.0, TAU * branches.len() as f64],
            Curve::Track { phi, .. } => *phi,
        }
    }

    /// The curve's point at `t` (binary64, from the kernel's enclosures).
    pub fn point(&self, t: f64) -> Result<Point3, String> {
        let e = |r: rusty_occt::Result<[[f64; 2]; 3]>| r.map(mid3).map_err(|e| e.to_string());
        match self {
            Curve::Line(p, d) => Ok(*p + *d * t),
            Curve::Circle(f, r) => Ok(f.origin() + (f.x() * t.cos() + f.y() * t.sin()) * *r),
            Curve::Ellipse(f, a, b) => {
                Ok(f.origin() + f.x() * (a * t.cos()) + f.y() * (b * t.sin()))
            }
            Curve::Hyperbola(f, a, b) => {
                Ok(f.origin() + f.x() * (a * t.cosh()) + f.y() * (b * t.sinh()))
            }
            Curve::Loop { curve, u } => {
                let (branch, s) = if t <= 1.0 {
                    (Branch::Plus, t)
                } else {
                    (Branch::Minus, 2.0 - t)
                };
                let v = u[0] + (u[1] - u[0]) * s.clamp(0.0, 1.0);
                e(curve.point_at([v, v], branch))
            }
            Curve::Turn {
                curve,
                start,
                branches,
            } => {
                let k = ((t / TAU).floor().max(0.0) as usize).min(branches.len() - 1);
                let v = start + t - TAU * k as f64;
                e(curve.point_at([v, v], branches[k]))
            }
            Curve::Track { curve, track, phi } => {
                e(curve.point_at(*track, t.clamp(phi[0], phi[1])))
            }
        }
    }

    /// The first line of DRAW's `dump` of the curve.
    pub fn kind(&self) -> &'static str {
        match self {
            Curve::Line(..) => "Line",
            Curve::Circle(..) => "Circle",
            Curve::Ellipse(..) => "Ellipse",
            Curve::Hyperbola(..) => "Hyperbola",
            Curve::Loop { .. } | Curve::Turn { .. } | Curve::Track { .. } => {
                "Procedural curve (rusty-occt)"
            }
        }
    }
}

/// `intersect`'s curves and points of two surfaces; `None` for the same
/// surface (DRAW reports no curve).
pub fn intersect(a: &Surface, b: &Surface) -> Result<(Vec<Curve>, Vec<Point3>), String> {
    let result = surface_surface(a, b).map_err(|e| e.to_string())?;
    let (mut curves, mut points) = (Vec::new(), Vec::new());
    match result {
        SurfaceIntersection::Empty | SurfaceIntersection::Same => {}
        SurfaceIntersection::NotConic => return Err("a curve not yet parameterised".into()),
        SurfaceIntersection::Items(items) => {
            for item in items {
                match item {
                    AnalyticItem::Point(p) => points.push(mid3(p)),
                    AnalyticItem::Line { point, direction } => {
                        let d = mid3(direction) - Point3::ORIGIN;
                        curves.push(Curve::Line(mid3(point), d * (1.0 / d.length())));
                    }
                    AnalyticItem::Circle {
                        centre,
                        normal,
                        radius,
                    } => curves.push(Curve::Circle(
                        frame(mid3(centre), mid3(normal) - Point3::ORIGIN, None)?,
                        mid(radius),
                    )),
                    AnalyticItem::Ellipse {
                        centre,
                        normal,
                        major,
                        semi_major,
                        semi_minor,
                    } => curves.push(Curve::Ellipse(
                        frame(
                            mid3(centre),
                            mid3(normal) - Point3::ORIGIN,
                            Some(mid3(major) - Point3::ORIGIN),
                        )?,
                        mid(semi_major),
                        mid(semi_minor),
                    )),
                    AnalyticItem::Hyperbola {
                        centre,
                        normal,
                        transverse,
                        semi_transverse,
                        semi_conjugate,
                    } => {
                        let x = mid3(transverse) - Point3::ORIGIN;
                        for axis in [x, x * -1.0] {
                            curves.push(Curve::Hyperbola(
                                frame(mid3(centre), mid3(normal) - Point3::ORIGIN, Some(axis))?,
                                mid(semi_transverse),
                                mid(semi_conjugate),
                            ));
                        }
                    }
                }
            }
        }
        SurfaceIntersection::Procedural(c) => {
            let c: Rc<ProceduralCurve> = Rc::from(c);
            for component in c.components() {
                curves.push(match component {
                    Component::Loop { u } => Curve::Loop {
                        curve: c.clone(),
                        u: [u[0][1], u[1][0]],
                    },
                    Component::Ring { branch } => Curve::Turn {
                        curve: c.clone(),
                        start: 0.0,
                        branches: vec![*branch],
                    },
                    Component::FigureEight { node } => Curve::Turn {
                        curve: c.clone(),
                        start: mid(*node),
                        branches: vec![Branch::Plus, Branch::Minus],
                    },
                });
            }
        }
        SurfaceIntersection::Traced(c) => {
            let c: Rc<TracedCurve> = Rc::from(c);
            for (k, track) in c.tracks().iter().enumerate() {
                curves.push(Curve::Track {
                    curve: c.clone(),
                    track: k,
                    phi: track.phi,
                });
            }
        }
    }
    Ok((curves, points))
}

/// The distance from a point to a surface's exact point set, in binary64.
pub fn distance(s: &Surface, p: Point3) -> f64 {
    let f = match s {
        Surface::Plane(f)
        | Surface::Cylinder { frame: f, .. }
        | Surface::Cone { frame: f, .. }
        | Surface::Sphere { frame: f, .. }
        | Surface::Torus { frame: f, .. } => *f,
        Surface::BSpline(_) => return f64::NAN,
    };
    let rel = p - f.origin();
    let n = f.normal();
    let h = rel.dot(n);
    let radial = (rel - n * h).length();
    match s {
        Surface::Plane(_) => h.abs(),
        Surface::Cylinder { radius, .. } => (radial - radius).abs(),
        Surface::Sphere { radius, .. } => (rel.length() - radius).abs(),
        Surface::Cone {
            radius, half_angle, ..
        } => {
            let (sa, ca) = half_angle.sin_cos();
            let shift = radius * ca + h * sa;
            (radial * ca - shift).abs().min((radial * ca + shift).abs())
        }
        Surface::Torus { major, minor, .. } => ((radial - major).hypot(h) - minor).abs(),
        Surface::BSpline(_) => f64::NAN,
    }
}

/// C's `%.<digits>g`.
pub fn format_g(x: f64, digits: usize) -> String {
    if x == 0.0 || !x.is_finite() {
        return format!("{x}");
    }
    let e = format!("{:.*e}", digits - 1, x);
    // Round first: the exponent after rounding decides the style.
    let (mantissa, exp) = e.split_once('e').expect("an exponent");
    let exp: i32 = exp.parse().expect("an exponent");
    let trim = |s: &str| {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if exp < -4 || exp >= digits as i32 {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mantissa), exp.abs())
    } else {
        let decimals = (digits as i32 - 1 - exp).max(0) as usize;
        trim(&format!("{:.*}", decimals, x))
    }
}

/// DRAW's `Draw::Atof`: `+ - * /`, parentheses, numbers, `pi` and named
/// numbers.
pub fn evaluate(text: &str, var: &dyn Fn(&str) -> Option<f64>) -> Result<f64, String> {
    struct P<'a> {
        s: &'a [u8],
        i: usize,
        var: &'a dyn Fn(&str) -> Option<f64>,
    }
    impl P<'_> {
        fn skip(&mut self) {
            while self.i < self.s.len() && (self.s[self.i] as char).is_whitespace() {
                self.i += 1;
            }
        }
        fn expr(&mut self) -> Result<f64, String> {
            let mut v = self.term()?;
            loop {
                self.skip();
                match self.s.get(self.i) {
                    Some(b'+') => {
                        self.i += 1;
                        v += self.term()?;
                    }
                    Some(b'-') => {
                        self.i += 1;
                        v -= self.term()?;
                    }
                    _ => return Ok(v),
                }
            }
        }
        fn term(&mut self) -> Result<f64, String> {
            let mut v = self.factor()?;
            loop {
                self.skip();
                match self.s.get(self.i) {
                    Some(b'*') => {
                        self.i += 1;
                        v *= self.factor()?;
                    }
                    Some(b'/') => {
                        self.i += 1;
                        v /= self.factor()?;
                    }
                    _ => return Ok(v),
                }
            }
        }
        fn factor(&mut self) -> Result<f64, String> {
            self.skip();
            match self.s.get(self.i) {
                Some(b'-') => {
                    self.i += 1;
                    Ok(-self.factor()?)
                }
                Some(b'+') => {
                    self.i += 1;
                    self.factor()
                }
                Some(b'(') => {
                    self.i += 1;
                    let v = self.expr()?;
                    self.skip();
                    if self.s.get(self.i) != Some(&b')') {
                        return Err("unbalanced parentheses".into());
                    }
                    self.i += 1;
                    Ok(v)
                }
                Some(c) if c.is_ascii_digit() || *c == b'.' => {
                    let start = self.i;
                    while self.i < self.s.len() {
                        let c = self.s[self.i];
                        let exponent_sign = (c == b'+' || c == b'-')
                            && self.i > start
                            && matches!(self.s[self.i - 1], b'e' | b'E');
                        if c.is_ascii_digit()
                            || c == b'.'
                            || c == b'e'
                            || c == b'E'
                            || exponent_sign
                        {
                            self.i += 1;
                        } else {
                            break;
                        }
                    }
                    let t = std::str::from_utf8(&self.s[start..self.i]).expect("ASCII");
                    t.parse().map_err(|_| format!("not a number: {t}"))
                }
                Some(c) if c.is_ascii_alphabetic() || *c == b'_' => {
                    let start = self.i;
                    while self.i < self.s.len()
                        && (self.s[self.i].is_ascii_alphanumeric() || self.s[self.i] == b'_')
                    {
                        self.i += 1;
                    }
                    let name = std::str::from_utf8(&self.s[start..self.i]).expect("ASCII");
                    if name == "pi" {
                        return Ok(std::f64::consts::PI);
                    }
                    // Draw's `sqrt(...)` (the Boolean group's `dset r sqrt(2)`).
                    self.skip();
                    if name == "sqrt" && self.s.get(self.i) == Some(&b'(') {
                        return Ok(self.factor()?.sqrt());
                    }
                    (self.var)(name).ok_or_else(|| format!("{name} is not a number"))
                }
                _ => Err("a malformed expression".into()),
            }
        }
    }
    // Draw::Atof tolerates trailing separators such as `4.99,`.
    let text = text.trim().trim_end_matches(',');
    let mut p = P {
        s: text.as_bytes(),
        i: 0,
        var,
    };
    let v = p.expr()?;
    p.skip();
    if p.i != p.s.len() {
        return Err(format!("a malformed expression: {text}"));
    }
    Ok(v)
}
