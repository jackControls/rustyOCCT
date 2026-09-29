//! Pairs of lines, circles, ellipses and hyperbolas' branches (S7d.1 of
//! `REVIEW_NOTES.md`).
//!
//! A line is `p0 + s (p1 - p0)`; a circle lies in the plane through its
//! origin normal to its stored normal, an ellipse or a hyperbola in the plane
//! through its origin spanned by its stored axes (normal `x × y`), all
//! rational.
//!
//! * Two lines: exact linear algebra.
//! * A line and a plane conic: a line crossing the plane meets it in one
//!   rational point, on the conic exactly when the conic's equation vanishes
//!   there; a line in the plane substitutes into the equation, a quadratic
//!   with rational coefficients (a double root a tangency).
//! * Two conics in one plane: S7c's resultant (`curve_surface::section`) with
//!   the second conic's equation for the surface's function (a multiple root a
//!   tangency, `c0 = c1 = 0` coincidence). In two crossing planes both
//!   equations along the common line are quadratics with rational
//!   coefficients; the points are the real roots of their gcd, tangent where
//!   the root is double in both (the curves then share that line as tangent).
//!
//! A hyperbola keeps the points of its branch, decided by a certified sign.
use super::analytic::{Enclosure, Enclosure3};
use super::conic_surface::{turn_angle, Conic};
use super::curve_surface::{section, Section};
use super::procedural::{bounds, bounds3, cross, dot, e3, edot, esub, limit, q, sub, zero, E, X};
use super::tangency::{badd, bdot, bmul, bscale, cauchy, padd, pmul, pscale, B, P};
use crate::certified::Interval as I;
use crate::polynomial::real::{isolate, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::topology::Curve3;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A curve `curve_curve` takes: a topology edge's line or circle (a
/// segment's whole line, an arc's whole circle) or an S7c.2 conic.
#[derive(Debug, Clone, PartialEq)]
pub enum AnalyticCurve {
    Edge(Curve3),
    Conic(Conic),
}

/// A point where two curves meet.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveCurvePoint {
    /// Both curves' parameters: a line's `s`, a circle's or an ellipse's
    /// angle, a hyperbola's `t`.
    pub parameters: [Enclosure; 2],
    pub point: Enclosure3,
    /// The curves share their tangent line here, else they cross.
    pub tangent: bool,
}

/// The intersection of two curves.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveCurveIntersection {
    Empty,
    /// The same point sets.
    Coincident,
    /// Points sorted by the first curve's parameter.
    Points(Vec<CurveCurvePoint>),
}

/// The intersection of two analytic curves (S7d.1).
pub fn curve_curve(a: &AnalyticCurve, b: &AnalyticCurve) -> Result<CurveCurveIntersection> {
    let (ea, eb) = (Exact::of(a)?, Exact::of(b)?);
    let found = match (&ea, &eb) {
        (Exact::Line(l1), Exact::Line(l2)) => line_line(l1, l2),
        (Exact::Line(l), Exact::Plane(c)) => line_conic(l, c)?,
        (Exact::Plane(c), Exact::Line(l)) => line_conic(l, c)?,
        (Exact::Plane(c1), Exact::Plane(c2)) => conic_conic(c1, c2)?,
    };
    Ok(match found {
        Found::Coincident => CurveCurveIntersection::Coincident,
        Found::Points(points) => {
            let mut out = Vec::new();
            for (p, tangent) in points {
                out.push(CurveCurvePoint {
                    parameters: [ea.parameter(&p)?, eb.parameter(&p)?],
                    point: bounds3(&p),
                    tangent,
                });
            }
            out.sort_by(|x, y| {
                (x.parameters[0][0] + x.parameters[0][1])
                    .total_cmp(&(y.parameters[0][0] + y.parameters[0][1]))
            });
            if out.is_empty() {
                CurveCurveIntersection::Empty
            } else {
                CurveCurveIntersection::Points(out)
            }
        }
    })
}

enum Found {
    Coincident,
    /// Points (enclosed) and their contacts.
    Points(Vec<(E<I>, bool)>),
}

pub(super) struct Line {
    p0: X,
    d: X,
}

/// A plane conic: its plane (origin, normal) and its equation.
pub(super) struct PlaneConic {
    o: X,
    normal: X,
    kind: Kind,
}

/// A conic's equation in its frame coordinates `(xi, eta)`.
type ConicEquation<'a> = Box<dyn Fn(&B, &B) -> B + 'a>;

enum Kind {
    /// Radius; the stored axes (for the angle); a rational basis of the
    /// plane (for the resultant).
    Circle { r: R, x: X, y: X, u: X, v: X },
    /// Semi-axes, stored axes, and the rows of `[x y (x × y)]^-1` giving
    /// the frame coordinates `xi`, `eta`.
    Conic {
        a: R,
        b: R,
        hyperbola: bool,
        x: X,
        y: X,
        rows: [X; 2],
    },
}

// Two per call, never stored: their sizes do not matter.
#[allow(clippy::large_enum_variant)]
pub(super) enum Exact {
    Line(Line),
    Plane(PlaneConic),
}

impl Exact {
    pub(super) fn of(c: &AnalyticCurve) -> Result<Self> {
        Ok(match c {
            AnalyticCurve::Edge(Curve3::LineSegment { start, end }) => {
                let p0 = start.to_array().map(q);
                let d = sub(&end.to_array().map(q), &p0);
                if zero(&dot(&d, &d)) {
                    return Err(Error::Degenerate("a line through one point"));
                }
                Exact::Line(Line { p0, d })
            }
            AnalyticCurve::Edge(
                Curve3::Circle { frame, radius } | Curve3::CircularArc { frame, radius, .. },
            ) => {
                if !(radius.is_finite() && *radius > 0.0) {
                    return Err(Error::OutOfDomain("a circle's radius"));
                }
                let n = frame.normal().to_array().map(q);
                let u = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                    .iter()
                    .map(|e| cross(&n, &e.map(|k| R::from_integer(k.into()))))
                    .find(|u| !u.iter().all(zero))
                    .expect("a nonzero normal");
                let v = cross(&n, &u);
                Exact::Plane(PlaneConic {
                    o: frame.origin().to_array().map(q),
                    normal: n,
                    kind: Kind::Circle {
                        r: q(*radius),
                        x: frame.x().to_array().map(q),
                        y: frame.y().to_array().map(q),
                        u,
                        v,
                    },
                })
            }
            AnalyticCurve::Edge(Curve3::BSpline(_)) => {
                return Err(Error::OutOfDomain("a spline edge (S7d.2)"))
            }
            // A hyperbola edge's whole branch (S8d.2).
            AnalyticCurve::Edge(Curve3::HyperbolaArc {
                frame,
                major,
                minor,
                ..
            }) => {
                return Self::of(&AnalyticCurve::Conic(Conic::Hyperbola {
                    frame: *frame,
                    major: *major,
                    minor: *minor,
                }))
            }
            AnalyticCurve::Edge(Curve3::ParabolaArc { .. }) => {
                return Err(Error::OutOfDomain("a parabola edge"))
            }
            AnalyticCurve::Edge(Curve3::Section(_)) => {
                return Err(Error::OutOfDomain("a torus section edge"))
            }
            AnalyticCurve::Edge(Curve3::Meet(_) | Curve3::Rise(_)) => {
                return Err(Error::OutOfDomain("two cylinders' meeting edge"))
            }
            // An ellipse edge's whole ellipse (S8a.2).
            AnalyticCurve::Edge(Curve3::EllipseArc {
                frame,
                major,
                minor,
                ..
            }) => {
                return Self::of(&AnalyticCurve::Conic(Conic::Ellipse {
                    frame: *frame,
                    major: *major,
                    minor: *minor,
                }))
            }
            AnalyticCurve::Conic(conic) => {
                let (frame, major, minor, hyperbola) = match *conic {
                    Conic::Ellipse {
                        frame,
                        major,
                        minor,
                    } => (frame, major, minor, false),
                    Conic::Hyperbola {
                        frame,
                        major,
                        minor,
                    } => (frame, major, minor, true),
                };
                if !(major.is_finite() && minor.is_finite() && minor > 0.0 && major > 0.0)
                    || (!hyperbola && major < minor)
                {
                    return Err(Error::OutOfDomain(
                        "a conic's semi-axes (positive, an ellipse's major the larger)",
                    ));
                }
                let (x, y) = (frame.x().to_array().map(q), frame.y().to_array().map(q));
                let m = cross(&x, &y);
                let det = dot(&x, &cross(&y, &m));
                let rows = [
                    cross(&y, &m).map(|c| c / &det),
                    cross(&m, &x).map(|c| c / &det),
                ];
                Exact::Plane(PlaneConic {
                    o: frame.origin().to_array().map(q),
                    normal: m,
                    kind: Kind::Conic {
                        a: q(major),
                        b: q(minor),
                        hyperbola,
                        x,
                        y,
                        rows,
                    },
                })
            }
        })
    }

    /// The curve's parameter at one of its points.
    fn parameter(&self, p: &E<I>) -> Result<Enclosure> {
        match self {
            Exact::Line(l) => {
                let t = edot(&esub(p, &e3(&l.p0)), &e3(&l.d))
                    .div(&I::exact(dot(&l.d, &l.d)))
                    .ok_or(limit("a line's parameter"))?;
                Ok(bounds(&t))
            }
            Exact::Plane(c) => c.parameter(p),
        }
    }
}

impl PlaneConic {
    /// The conic's equation at a point given by polynomials (valid in its
    /// plane).
    fn equation(&self, p: &[B; 3]) -> B {
        let w: [B; 3] = std::array::from_fn(|i| badd(&p[i], &vec![vec![-self.o[i].clone()]]));
        match &self.kind {
            Kind::Circle { r, .. } => badd(&bdot(&w, &w), &vec![vec![-(r * r)]]),
            Kind::Conic {
                a,
                b,
                hyperbola,
                rows,
                ..
            } => {
                let coordinate = |row: &X| {
                    (0..3).fold(Vec::new(), |acc: B, i| badd(&acc, &bscale(&w[i], &row[i])))
                };
                let (xi, eta) = (coordinate(&rows[0]), coordinate(&rows[1]));
                let sign = if *hyperbola { -1 } else { 1 };
                badd(
                    &badd(
                        &bscale(&bmul(&xi, &xi), &(b * b)),
                        &bscale(&bmul(&eta, &eta), &(a * a * R::from_integer(sign.into()))),
                    ),
                    &vec![vec![-(a * a * b * b)]],
                )
            }
        }
    }

    /// Frame coordinates `(xi, eta)` of a point of the plane (a conic's).
    fn coordinates(&self, p: &E<I>) -> Option<(I, I)> {
        match &self.kind {
            Kind::Conic { rows, .. } => {
                let w = esub(p, &e3(&self.o));
                Some((edot(&w, &e3(&rows[0])), edot(&w, &e3(&rows[1]))))
            }
            Kind::Circle { .. } => None,
        }
    }

    /// Whether a point of the conic lies on its branch (always, but for a
    /// hyperbola: `xi > 0`, certainly).
    pub(super) fn on_branch(&self, p: &E<I>) -> Result<bool> {
        match &self.kind {
            Kind::Conic {
                hyperbola: true, ..
            } => match self.coordinates(p).expect("a conic").0.sign() {
                Some(Ordering::Greater) => Ok(true),
                Some(Ordering::Less) => Ok(false),
                _ => Err(limit("a hyperbola's branch")),
            },
            _ => Ok(true),
        }
    }

    pub(super) fn parameter(&self, p: &E<I>) -> Result<Enclosure> {
        let t = match &self.kind {
            Kind::Circle { x, y, .. } => {
                let w = esub(p, &e3(&self.o));
                turn_angle(&edot(&w, &e3(y)), &edot(&w, &e3(x)))?
            }
            Kind::Conic {
                a, b, hyperbola, ..
            } => {
                let (xi, eta) = self.coordinates(p).expect("a conic");
                let u = xi
                    .div(&I::exact(a.clone()))
                    .ok_or(limit("a conic's parameter"))?;
                let v = eta
                    .div(&I::exact(b.clone()))
                    .ok_or(limit("a conic's parameter"))?;
                if *hyperbola {
                    crate::certified::asinh(&v)
                } else {
                    turn_angle(&v, &u)?
                }
            }
        };
        Ok(bounds(&t))
    }

    /// `(origin, e1, e2, conic(xi, eta))` for S7c's resultant.
    fn section_basis(&self) -> (X, X, ConicEquation<'_>) {
        match &self.kind {
            Kind::Circle { r, u, v, .. } => (
                u.clone(),
                v.clone(),
                Box::new(move |xi: &B, eta: &B| {
                    let w: [B; 3] =
                        std::array::from_fn(|i| badd(&bscale(xi, &u[i]), &bscale(eta, &v[i])));
                    badd(&bdot(&w, &w), &vec![vec![-(r * r)]])
                }),
            ),
            Kind::Conic {
                a,
                b,
                hyperbola,
                x,
                y,
                ..
            } => (
                x.clone(),
                y.clone(),
                Box::new(move |xi: &B, eta: &B| {
                    let sign = if *hyperbola { -1 } else { 1 };
                    badd(
                        &badd(
                            &bscale(&bmul(xi, xi), &(b * b)),
                            &bscale(&bmul(eta, eta), &(a * a * R::from_integer(sign.into()))),
                        ),
                        &vec![vec![-(a * a * b * b)]],
                    )
                }),
            ),
        }
    }

    /// The plane's and the equation's homogeneous forms `n . (X - o W)`
    /// and `W^2 E(X / W)` for polynomials `(X, W)` (a spline's span).
    pub(super) fn homogeneous(&self, h: &[P; 4]) -> (P, P) {
        let w: [P; 3] = std::array::from_fn(|i| padd(&h[i], &pscale(&h[3], &-&self.o[i])));
        let dot_p = |a: &X| (0..3).fold(Vec::new(), |acc: P, i| padd(&acc, &pscale(&w[i], &a[i])));
        let w2 = pmul(&h[3], &h[3]);
        let equation = match &self.kind {
            Kind::Circle { r, .. } => {
                let ww = (0..3).fold(Vec::new(), |acc: P, i| padd(&acc, &pmul(&w[i], &w[i])));
                padd(&ww, &pscale(&w2, &-(r * r)))
            }
            Kind::Conic {
                a,
                b,
                hyperbola,
                rows,
                ..
            } => {
                let (xi, eta) = (dot_p(&rows[0]), dot_p(&rows[1]));
                let sign = if *hyperbola { -1 } else { 1 };
                padd(
                    &padd(
                        &pscale(&pmul(&xi, &xi), &(b * b)),
                        &pscale(&pmul(&eta, &eta), &(a * a * R::from_integer(sign.into()))),
                    ),
                    &pscale(&w2, &-(a * a * b * b)),
                )
            }
        };
        (dot_p(&self.normal), equation)
    }

    /// The equation along `p0 + t d` (a polynomial in `t`).
    fn along(&self, p0: &X, d: &X) -> P {
        let p: [B; 3] = std::array::from_fn(|i| vec![vec![p0[i].clone(), d[i].clone()]]);
        let mut f = self.equation(&p).into_iter().next().unwrap_or_default();
        while f.last().is_some_and(zero) {
            f.pop();
        }
        f
    }
}

fn line_line(l1: &Line, l2: &Line) -> Found {
    let n = cross(&l1.d, &l2.d);
    let w = sub(&l2.p0, &l1.p0);
    if n.iter().all(zero) {
        return if cross(&w, &l1.d).iter().all(zero) {
            Found::Coincident
        } else {
            Found::Points(Vec::new())
        };
    }
    if !zero(&dot(&w, &n)) {
        return Found::Points(Vec::new());
    }
    let s = dot(&cross(&w, &l2.d), &n) / dot(&n, &n);
    let p: X = std::array::from_fn(|i| &l1.p0[i] + &l1.d[i] * &s);
    Found::Points(vec![(e3::<I>(&p), false)])
}

/// The real roots of an exact polynomial with their multiplicities.
fn roots(f: &P) -> Result<Vec<(I, usize)>> {
    let int = IntPolynomial::from_rationals(f);
    if int.is_zero() {
        return Err(limit("a vanishing polynomial"));
    }
    let bound = cauchy(f);
    let mut budget = Budget::new(RootIsolationOptions::default());
    let mut out = Vec::new();
    for mut root in isolate(&int, -bound.clone(), bound, &mut budget)? {
        let m = root.multiplicity();
        root.refine_for_signs(300);
        let (lo, hi) = root.isolator();
        out.push((I::new(lo.clone(), hi.clone()), m));
    }
    Ok(out)
}

fn line_conic(l: &Line, c: &PlaneConic) -> Result<Found> {
    let dn = dot(&l.d, &c.normal);
    let off = dot(&sub(&l.p0, &c.o), &c.normal);
    let mut out = Vec::new();
    if !zero(&dn) {
        // One rational point of the plane: on the conic exactly or not.
        let s = -off / dn;
        let p: X = std::array::from_fn(|i| &l.p0[i] + &l.d[i] * &s);
        let at: [B; 3] = std::array::from_fn(|i| vec![vec![p[i].clone()]]);
        let value = c.equation(&at).first().and_then(|v| v.first().cloned());
        let point = e3::<I>(&p);
        if value.as_ref().is_none_or(zero) && c.on_branch(&point)? {
            out.push((point, false));
        }
        return Ok(Found::Points(out));
    }
    if !zero(&off) {
        return Ok(Found::Points(out));
    }
    // In the plane: the conic's equation along the line.
    for (t, m) in roots(&c.along(&l.p0, &l.d))? {
        let point: E<I> = std::array::from_fn(|i| {
            I::exact(l.p0[i].clone()).add(&t.mul(&I::exact(l.d[i].clone())))
        });
        if c.on_branch(&point)? {
            out.push((point, m > 1));
        }
    }
    Ok(Found::Points(out))
}

fn conic_conic(c1: &PlaneConic, c2: &PlaneConic) -> Result<Found> {
    let dl = cross(&c1.normal, &c2.normal);
    if dl.iter().all(zero) {
        if !zero(&dot(&sub(&c2.o, &c1.o), &c1.normal)) {
            return Ok(Found::Points(Vec::new()));
        }
        // One plane: S7c's resultant with the second conic's equation.
        let (e1, e2, conic) = c1.section_basis();
        let found = section(&c1.o, &e1, &e2, conic.as_ref(), &|p| c2.equation(p))?;
        return Ok(match found {
            Section::Contained => {
                // The same conic; two hyperbolas may take opposite branches.
                if let Kind::Conic {
                    a,
                    hyperbola: true,
                    x,
                    ..
                } = &c1.kind
                {
                    let vertex: X = std::array::from_fn(|i| &c1.o[i] + &x[i] * a);
                    if !c2.on_branch(&e3::<I>(&vertex))? {
                        return Ok(Found::Points(Vec::new()));
                    }
                }
                Found::Coincident
            }
            Section::Points(points) => {
                let mut out = Vec::new();
                for p in points {
                    if c1.on_branch(&p.point)? && c2.on_branch(&p.point)? {
                        out.push((p.point, p.tangent));
                    }
                }
                Found::Points(out)
            }
        });
    }
    // Crossing planes: their common line x0 + t dl.
    let rhs = [dot(&c1.normal, &c1.o), dot(&c2.normal, &c2.o)];
    // x0 = (rhs0 (n2 × dl) + rhs1 (dl × n1)) / |dl|^2 solves n1.x0 = rhs0,
    // n2.x0 = rhs1, dl.x0 = 0.
    let (a1, a2) = (cross(&c2.normal, &dl), cross(&dl, &c1.normal));
    let ll = dot(&dl, &dl);
    let x0: X = std::array::from_fn(|i| (&rhs[0] * &a1[i] + &rhs[1] * &a2[i]) / &ll);
    let (q1, q2) = (c1.along(&x0, &dl), c2.along(&x0, &dl));
    let (i1, i2) = (
        IntPolynomial::from_rationals(&q1),
        IntPolynomial::from_rationals(&q2),
    );
    let common = i1.gcd(&i2);
    let mut out = Vec::new();
    if common.is_constant() {
        return Ok(Found::Points(out));
    }
    let g: P = common.0.iter().cloned().map(R::from_integer).collect();
    // A root double in both quadratics: each touches the line there.
    let double =
        |f: &P| f.len() == 3 && zero(&(&f[1] * &f[1] - R::from_integer(4.into()) * &f[0] * &f[2]));
    let tangent = double(&q1) && double(&q2);
    for (t, _) in roots(&g)? {
        let point: E<I> =
            std::array::from_fn(|i| I::exact(x0[i].clone()).add(&t.mul(&I::exact(dl[i].clone()))));
        if c1.on_branch(&point)? && c2.on_branch(&point)? {
            out.push((point, tangent));
        }
    }
    Ok(Found::Points(out))
}
