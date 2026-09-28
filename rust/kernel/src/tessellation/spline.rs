//! Certified bounds of spline cells for tessellation (T-b of the tessellation
//! track, `MATHEMATICS.md`, tessellation bounds for splines).
//!
//! A spline's range is cut into its exact Bézier pieces (patches), each
//! halved by de Casteljau in the `Fast` tier into cells whose homogeneous
//! controls enclose the exact subdivision's. On a cell of parameter length
//! `L`, degree `p`, controls `(w_i P_i, w_i)` and `Q_i = P_i - c` for its
//! first control `c`: `R = max |Q_i|`, `ω = min w_i`, the numerator
//! `A = w (C - c)` and `w` bounded by their nets' first and second
//! differences, then `|C'| <= D1 = (|A'| + |w'| R) / ω` and
//! `|C''| <= D2 = (|A''| + 2 |w'| D1 + |w''| R) / ω`; a surface's `a`, `b`,
//! `c` bound `|S_uu|`, `|S_uv|`, `|S_vv|` the same way. Points and first
//! derivatives are evaluated by de Casteljau on a cell in the same tier.
use crate::certified::{Fast, Real};
use crate::{BSplineCurve3, BSplineSurface3, Error, Point2, Result, Vec3};
use num_rational::BigRational as R;
use poly::{Poly, Poly3};

mod poly;

type V3 = [Fast; 3];

/// Halvings of each Bézier piece of a curve, and at most this many cells.
const CURVE_DEPTH: usize = 3;
const MAX_CURVE_CELLS: usize = 4096;
/// Halvings of each patch in each direction, fewer when a surface has many
/// patches: at most this many cells, and this many controls in all.
const SURFACE_DEPTH: usize = 3;
const MAX_SURFACE_CELLS: usize = 1024;
const MAX_SURFACE_CONTROLS: usize = 1 << 18;
/// Degrees up to which a rational cell's derivatives are bounded through the
/// Bernstein coefficients of their numerators (above, the quotient rule's
/// triangle inequality), and up to which a surface cell carries its
/// normal's cone.
const MAX_RATIONAL_CURVE_DEGREE: usize = 6;
const MAX_RATIONAL_SURFACE_DEGREE: usize = 3;
const MAX_CONE_DEGREE: usize = 6;
/// Degrees up to which a surface cell bounds its normal's rates through
/// `M × M_u` and `M × M_v` (nonrational, rational).
const MAX_RATE_DEGREE: [usize; 2] = [4, 2];

fn c(x: f64) -> Fast {
    Fast::exact_f64(x)
}

fn upper(x: &Fast) -> f64 {
    x.bounds_f64().1
}

fn lower(x: &Fast) -> f64 {
    x.bounds_f64().0
}

/// An upper bound of a vector's length.
fn norm_upper(v: &V3) -> f64 {
    upper(&v[0].square().add(&v[1].square()).add(&v[2].square()).sqrt())
}

/// A lower bound of a vector's length (each component's least magnitude).
fn norm_lower(v: &V3) -> f64 {
    let least = |x: &Fast| {
        let (lo, hi) = x.bounds_f64();
        if lo > 0.0 {
            lo
        } else if hi < 0.0 {
            -hi
        } else {
            0.0
        }
    };
    let squared = c(least(&v[0]))
        .square()
        .add(&c(least(&v[1])).square())
        .add(&c(least(&v[2])).square());
    lower(&squared.sqrt()).max(0.0)
}

fn sub3(a: &V3, b: &V3) -> V3 {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

fn scale3(a: &V3, k: &Fast) -> V3 {
    [a[0].mul(k), a[1].mul(k), a[2].mul(k)]
}

fn cross3(a: &V3, b: &V3) -> V3 {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

/// A homogeneous control `(w x, w y, w z, w)`, enclosed.
#[derive(Debug, Clone, Copy)]
struct H([Fast; 4]);

impl H {
    fn exact(x: &[R; 4]) -> Self {
        H([
            Fast::from_r(&x[0]),
            Fast::from_r(&x[1]),
            Fast::from_r(&x[2]),
            Fast::from_r(&x[3]),
        ])
    }
    fn mix(&self, o: &Self, s: &Fast, t: &Fast) -> Self {
        H(std::array::from_fn(|k| {
            self.0[k].mul(s).add(&o.0[k].mul(t))
        }))
    }
}

/// De Casteljau's halves of Bernstein controls, enclosed.
fn halves(a: &[H]) -> (Vec<H>, Vec<H>) {
    let half = c(0.5);
    let mut row = a.to_vec();
    let (mut left, mut right) = (vec![row[0]], vec![row[row.len() - 1]]);
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|w| w[0].mix(&w[1], &half, &half))
            .collect();
        left.push(row[0]);
        right.push(row[row.len() - 1]);
    }
    right.reverse();
    (left, right)
}

/// The value of Bernstein controls at an enclosed `t`.
fn value<T: Copy>(a: &[T], t: &Fast, mix: impl Fn(&T, &T, &Fast, &Fast) -> T) -> T {
    let s = c(1.0).sub(t);
    let mut row = a.to_vec();
    while row.len() > 1 {
        row = row.windows(2).map(|w| mix(&w[0], &w[1], &s, t)).collect();
    }
    row[0]
}

fn mix3(a: &V3, b: &V3, s: &Fast, t: &Fast) -> V3 {
    std::array::from_fn(|k| a[k].mul(s).add(&b[k].mul(t)))
}

fn mix1(a: &Fast, b: &Fast, s: &Fast, t: &Fast) -> Fast {
    a.mul(s).add(&b.mul(t))
}

/// Translated controls `(w_i (P_i - c), w_i)` of a net about the point
/// nearest its first control, with `R = max |P_i - c|` and `ω = min w_i`
/// (`None` when a weight's sign is not certain).
#[allow(clippy::type_complexity)]
fn translated(net: &[H]) -> Option<([f64; 3], Vec<V3>, Vec<Fast>, f64, f64)> {
    let first = net[0].0;
    let mid = |x: &Fast| {
        let (lo, hi) = x.bounds_f64();
        0.5 * lo + 0.5 * hi
    };
    let w0 = mid(&first[3]);
    let origin = [
        mid(&first[0]) / w0,
        mid(&first[1]) / w0,
        mid(&first[2]) / w0,
    ];
    if !origin.iter().all(|x| x.is_finite()) {
        return None;
    }
    let o = [c(origin[0]), c(origin[1]), c(origin[2])];
    let mut wq = Vec::with_capacity(net.len());
    let mut w = Vec::with_capacity(net.len());
    let (mut reach, mut least) = (0.0f64, f64::INFINITY);
    for h in net {
        let weight = h.0[3];
        if weight.sign() != Some(std::cmp::Ordering::Greater) {
            return None;
        }
        let q: V3 = std::array::from_fn(|k| h.0[k].sub(&weight.mul(&o[k])));
        let p: V3 = [q[0].div(&weight)?, q[1].div(&weight)?, q[2].div(&weight)?];
        reach = reach.max(norm_upper(&p));
        least = least.min(lower(&weight));
        wq.push(q);
        w.push(weight);
    }
    Some((origin, wq, w, reach, least))
}

/// `max |Δ^k x_i|` over a row of vectors (`k` 1 or 2).
fn difference3(x: &[V3], k: usize) -> f64 {
    if x.len() <= k {
        return 0.0;
    }
    x.windows(k + 1)
        .map(|w| {
            if k == 1 {
                norm_upper(&sub3(&w[1], &w[0]))
            } else {
                norm_upper(&sub3(&sub3(&w[2], &w[1]), &sub3(&w[1], &w[0])))
            }
        })
        .fold(0.0, f64::max)
}

fn difference1(x: &[Fast], k: usize) -> f64 {
    if x.len() <= k {
        return 0.0;
    }
    x.windows(k + 1)
        .map(|w| {
            let d = if k == 1 {
                w[1].sub(&w[0])
            } else {
                w[2].sub(&w[1]).sub(&w[1].sub(&w[0]))
            };
            let (lo, hi) = d.bounds_f64();
            (-lo).max(hi)
        })
        .fold(0.0, f64::max)
}

/// `p / L` and `p (p - 1) / L²`, rounded up (infinite when `L` is not
/// certainly positive).
fn scales(p: usize, length: &Fast) -> (Fast, Fast) {
    let p = c(p as f64);
    let first = p.div(length).unwrap_or(c(f64::INFINITY));
    let second = p
        .mul(&p.sub(&c(1.0)))
        .div(&length.square())
        .unwrap_or(c(f64::INFINITY));
    (first, second)
}

fn up(x: Fast) -> f64 {
    let u = upper(&x);
    if u.is_nan() {
        f64::INFINITY
    } else {
        u.max(0.0)
    }
}

/// Increasing parameter breaks (cell ends), exact, with their binary64
/// values where exact (compared directly) and their enclosures.
#[derive(Debug, Clone)]
struct Breaks {
    exact: Vec<R>,
    value: Vec<Option<f64>>,
    enclosed: Vec<Fast>,
}

impl Breaks {
    fn new(exact: Vec<R>) -> Self {
        let value = exact
            .iter()
            .map(|r| {
                let (lo, hi) = Fast::from_r(r).bounds_f64();
                (lo == hi).then_some(lo)
            })
            .collect();
        let enclosed = exact.iter().map(Fast::from_r).collect();
        Self {
            exact,
            value,
            enclosed,
        }
    }

    fn cmp(&self, k: usize, x: f64) -> std::cmp::Ordering {
        match self.value[k] {
            Some(b) => x.total_cmp(&b),
            None => R::from_float(x).map_or(std::cmp::Ordering::Equal, |x| x.cmp(&self.exact[k])),
        }
    }

    /// The cell `k` with `x <= b[k + 1]` first (the last when beyond).
    fn find(&self, x: f64) -> usize {
        let n = self.exact.len() - 1;
        let (mut lo, mut hi) = (0, n - 1);
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.cmp(mid + 1, x) == std::cmp::Ordering::Greater {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// `x`'s local parameter in cell `k`, clamped to `[0, 1]`.
    fn local(&self, k: usize, x: f64, length: &Fast) -> Option<Fast> {
        use std::cmp::Ordering::{Greater, Less};
        if self.cmp(k, x) != Greater {
            Some(c(0.0))
        } else if self.cmp(k + 1, x) != Less {
            Some(c(1.0))
        } else {
            c(x).sub(&self.enclosed[k]).div(length)
        }
    }
}

// ------------------------------------------------------------------ curves

#[derive(Debug, Clone)]
struct CurveCell {
    lo: R,
    hi: R,
    /// Outward binary64 bounds of `lo` and `hi`.
    span: (f64, f64),
    length: Fast,
    origin: [f64; 3],
    wq: Vec<V3>,
    w: Vec<Fast>,
    /// Upper bounds of `|C'|` and `|C''|` over the cell.
    d1: f64,
    d2: f64,
    /// A positive multiple `N` of `C'` in the cell's own parameter, with
    /// upper bounds of `|N'|` and `|N × N'|`: the unit tangent turns at the
    /// rate `|N × N'| / |N|²`.
    tangent: Option<(Poly3, f64, f64)>,
}

/// The cells of a spline curve over a range, in increasing parameter.
#[derive(Debug, Clone)]
pub(super) struct CurveCells {
    cells: Vec<CurveCell>,
    breaks: Breaks,
    degree: usize,
}

impl CurveCells {
    pub(super) fn new(curve: &BSplineCurve3, range: [f64; 2]) -> Result<Self> {
        let arcs = curve.bezier_arcs_in(range[0], range[1])?;
        let degree = curve.degree();
        let mut depth = CURVE_DEPTH;
        while depth > 0 && arcs.len() << depth > MAX_CURVE_CELLS {
            depth -= 1;
        }
        let mut cells = Vec::with_capacity(arcs.len() << depth);
        for arc in &arcs {
            let net: Vec<H> = arc.homogeneous_poles().iter().map(H::exact).collect();
            let [a, b] = arc.domain().clone();
            let mut pieces = vec![(a, b, net)];
            for _ in 0..depth {
                let mut next = Vec::with_capacity(2 * pieces.len());
                for (a, b, net) in pieces {
                    let m = (&a + &b) / R::from_integer(2.into());
                    let (l, r) = halves(&net);
                    next.push((a, m.clone(), l));
                    next.push((m, b, r));
                }
                pieces = next;
            }
            for (lo, hi, net) in pieces {
                cells.push(curve_cell(lo, hi, &net, degree));
            }
        }
        if cells.is_empty() {
            return Err(Error::ComputationLimit("tessellation spline cells"));
        }
        let breaks = Breaks::new(
            cells
                .iter()
                .map(|x| x.lo.clone())
                .chain(std::iter::once(cells[cells.len() - 1].hi.clone()))
                .collect(),
        );
        Ok(Self {
            cells,
            breaks,
            degree,
        })
    }

    /// The largest `|C'|` bound: the curve's speed in its parameter.
    pub(super) fn speed(&self) -> f64 {
        self.cells.iter().map(|x| x.d1).fold(0.0, f64::max)
    }

    /// The largest `|C''|` bound over the whole range.
    pub(super) fn curvature(&self) -> f64 {
        self.cells.iter().map(|x| x.d2).fold(0.0, f64::max)
    }

    /// The cell whose exact range contains `u` (clamped to the range).
    fn locate(&self, u: f64) -> Option<(&CurveCell, Fast)> {
        if !u.is_finite() {
            return None;
        }
        let k = self.breaks.find(u);
        let cell = &self.cells[k];
        Some((cell, self.breaks.local(k, u, &cell.length)?))
    }

    /// The exact point at parameter `u` (the range's end nearest when `u`
    /// lies beyond it), enclosed.
    pub(super) fn point(&self, u: f64) -> Option<V3> {
        let (cell, t) = self.locate(u)?;
        if cell.wq.is_empty() {
            return None;
        }
        let a = value(&cell.wq, &t, mix3);
        let w = value(&cell.w, &t, mix1);
        Some([
            c(cell.origin[0]).add(&a[0].div(&w)?),
            c(cell.origin[1]).add(&a[1].div(&w)?),
            c(cell.origin[2]).add(&a[2].div(&w)?),
        ])
    }

    /// `C'(u)`, enclosed.
    fn derivative(&self, u: f64) -> Option<V3> {
        let (cell, t) = self.locate(u)?;
        if cell.wq.len() < 2 {
            return None;
        }
        let (k, _) = scales(self.degree, &cell.length);
        let dq: Vec<V3> = cell
            .wq
            .windows(2)
            .map(|x| scale3(&sub3(&x[1], &x[0]), &k))
            .collect();
        let dw: Vec<Fast> = cell.w.windows(2).map(|x| x[1].sub(&x[0]).mul(&k)).collect();
        let a = value(&cell.wq, &t, mix3);
        let w = value(&cell.w, &t, mix1);
        let a1 = value(&dq, &t, mix3);
        let w1 = value(&dw, &t, mix1);
        let p: V3 = [a[0].div(&w)?, a[1].div(&w)?, a[2].div(&w)?];
        Some([
            a1[0].sub(&w1.mul(&p[0])).div(&w)?,
            a1[1].sub(&w1.mul(&p[1])).div(&w)?,
            a1[2].sub(&w1.mul(&p[2])).div(&w)?,
        ])
    }

    /// The certified deviation (without node gaps) and tangent turn of the
    /// arc between parameters `a` and `b` against its chord: `D2 h² / 8`
    /// and `D2 h / s`, `s` a lower bound of the speed on the arc.
    pub(super) fn segment(&self, a: f64, b: f64) -> (f64, f64) {
        let (lo, hi) = (a.min(b), a.max(b));
        let d2 = self
            .cells
            .iter()
            .filter(|x| x.span.0 <= hi && x.span.1 >= lo)
            .map(|x| x.d2)
            .fold(0.0, f64::max);
        let h = c(hi).sub(&c(lo));
        let deviation = up(c(d2).mul(&h.square()).mul(&c(0.125)));
        let mid = 0.5 * lo + 0.5 * hi;
        let reach = c(mid).sub(&c(lo)).union(&c(hi).sub(&c(mid)));
        let speed = match self.derivative(mid) {
            Some(d) => lower(&c(norm_lower(&d)).sub(&c(d2).mul(&reach))),
            None => 0.0,
        };
        let turn = if speed > 0.0 {
            up(c(d2).mul(&h).div(&c(speed)).unwrap_or(c(f64::INFINITY)))
        } else {
            f64::INFINITY
        };
        (deviation, turn.min(self.tangent_turn(lo, hi)))
    }

    /// The tangent's turn over `[lo, hi]` from the cells' tangent
    /// numerators: over each cell's part `[τ_a, τ_b]` of its own parameter,
    /// at most `|N × N'| (τ_b - τ_a) / μ²`, `μ` the enclosed `|N|` at the
    /// part's middle less `|N'| (τ_b - τ_a)`.
    fn tangent_turn(&self, lo: f64, hi: f64) -> f64 {
        let mut total = c(0.0);
        for cell in self
            .cells
            .iter()
            .filter(|x| x.span.0 <= hi && x.span.1 >= lo)
        {
            let Some((n, slope, rate)) = &cell.tangent else {
                return f64::INFINITY;
            };
            let local = |x: f64| {
                c(x).sub(&Fast::from_r(&cell.lo))
                    .div(&cell.length)
                    .map(|t| t.bounds_f64())
            };
            let (Some(a), Some(b)) = (local(lo), local(hi)) else {
                return f64::INFINITY;
            };
            let (ta, tb) = (a.0.max(0.0), b.1.min(1.0));
            if ta > tb {
                continue;
            }
            let width = c(tb).sub(&c(ta));
            let mid = c(0.5 * ta + 0.5 * tb);
            let zero = c(0.0);
            let at = [
                n[0].value(&mid, &zero),
                n[1].value(&mid, &zero),
                n[2].value(&mid, &zero),
            ];
            let least = lower(&c(norm_lower(&at)).sub(&c(*slope).mul(&width)));
            if least.is_nan() || least <= 0.0 {
                return f64::INFINITY;
            }
            let Some(part) = c(*rate).mul(&width).div(&c(least).square()) else {
                return f64::INFINITY;
            };
            total = total.add(&part);
        }
        up(total)
    }
}

fn curve_cell(lo: R, hi: R, net: &[H], p: usize) -> CurveCell {
    let length = Fast::from_r(&(&hi - &lo));
    let span = (
        Fast::from_r(&lo).bounds_f64().0,
        Fast::from_r(&hi).bounds_f64().1,
    );
    let Some((origin, wq, w, reach, least)) = translated(net) else {
        return CurveCell {
            lo,
            hi,
            span,
            length,
            origin: [0.0; 3],
            wq: Vec::new(),
            w: Vec::new(),
            d1: f64::INFINITY,
            d2: f64::INFINITY,
            tangent: None,
        };
    };
    let rational = w.iter().any(|x| x != &w[0]);
    let a: Poly3 = std::array::from_fn(|k| Poly::new([p, 0], wq.iter().map(|x| x[k]).collect()));
    let ww = Poly::new([p, 0], w.clone());
    let tangent = (!rational || p <= MAX_RATIONAL_CURVE_DEGREE).then(|| {
        let n = if rational {
            poly::sub33(
                &poly::mul3(&poly::derivative3(&a, 0), &ww),
                &poly::mul3(&a, &ww.derivative(0)),
            )
        } else {
            poly::derivative3(&a, 0)
        };
        let dn = poly::derivative3(&n, 0);
        let rate = poly::largest3(&poly::cross3(&n, &dn));
        let slope = poly::largest3(&dn);
        (n, slope, rate)
    });
    let (d1, d2) = if rational && p <= MAX_RATIONAL_CURVE_DEGREE {
        // C' = N1 / w², C'' = N2 / w³ in the cell's own parameter, with
        // N1 = A' w - A w' and N2 = (A'' w - A w'') w - 2 w' N1.
        let (a1, w1) = (poly::derivative3(&a, 0), ww.derivative(0));
        let (a2, w2) = (poly::derivative3(&a1, 0), w1.derivative(0));
        let n1 = poly::sub33(&poly::mul3(&a1, &ww), &poly::mul3(&a, &w1));
        let n2 = poly::sub33(
            &poly::mul3(
                &poly::sub33(&poly::mul3(&a2, &ww), &poly::mul3(&a, &w2)),
                &ww,
            ),
            &poly::scale3(&poly::mul3(&n1, &w1), 2.0),
        );
        let least = c(least);
        let d1 = c(poly::largest3(&n1))
            .div(&least.square())
            .and_then(|x| x.div(&length));
        let d2 = c(poly::largest3(&n2))
            .div(&least.square().mul(&least))
            .and_then(|x| x.div(&length.square()));
        (
            d1.map(up).unwrap_or(f64::INFINITY),
            d2.map(up).unwrap_or(f64::INFINITY),
        )
    } else {
        let (k1, k2) = scales(p, &length);
        let a1 = k1.mul(&c(difference3(&wq, 1)));
        let a2 = k2.mul(&c(difference3(&wq, 2)));
        let w1 = k1.mul(&c(difference1(&w, 1)));
        let w2 = k2.mul(&c(difference1(&w, 2)));
        let (reach, least) = (c(reach), c(least));
        let d1 = a1.add(&w1.mul(&reach)).div(&least);
        let d1 = d1.map(up).unwrap_or(f64::INFINITY);
        let d2 = a2
            .add(&c(2.0).mul(&w1).mul(&c(d1)))
            .add(&w2.mul(&reach))
            .div(&least);
        (d1, d2.map(up).unwrap_or(f64::INFINITY))
    };
    CurveCell {
        lo,
        hi,
        span,
        length,
        origin,
        wq,
        w,
        d1,
        d2,
        tangent,
    }
}

// ------------------------------------------------------------------ surfaces

#[derive(Debug, Clone)]
struct SurfaceCell {
    lengths: [Fast; 2],
    origin: [f64; 3],
    /// Controls `(p + 1) × (q + 1)`, u-major.
    wq: Vec<V3>,
    w: Vec<Fast>,
    /// Upper bounds of `|S_uu|`, `|S_uv|`, `|S_vv|`, `|S_u|`, `|S_v|`.
    bounds: [f64; 5],
    /// The cone of the cell's normals, from the Bernstein coefficients of a
    /// positive multiple `M` of `S_u × S_v`.
    cone: Option<Cone>,
    /// Upper bounds of `|M × M_u|` and `|M × M_v|` per unit `u` and `v`:
    /// the unit normal turns at the rates `|M × M_u| / |M|²` and likewise.
    rates: Option<[f64; 2]>,
}

/// Every coefficient of `M` within `spread` of `axis`, so every normal of the
/// cell (a convex combination's direction); `|M| >= least` on the cell (its
/// component along `axis`); `|S_u × S_v| >= factor |M|`.
#[derive(Debug, Clone)]
struct Cone {
    axis: Vec3,
    spread: f64,
    least: f64,
    factor: f64,
}

/// The cone of a vector polynomial's coefficients about its value at the
/// middle of the cell; `None` unless every coefficient is certainly within
/// 90° of it.
fn cone_of(m: &[V3], centre: &V3, factor: f64) -> Option<Cone> {
    let mid = |x: &Fast| {
        let (lo, hi) = x.bounds_f64();
        0.5 * lo + 0.5 * hi
    };
    let axis = Vec3::new(mid(&centre[0]), mid(&centre[1]), mid(&centre[2]));
    let length = axis.length();
    if !(length.is_finite() && length > 0.0) {
        return None;
    }
    let d = [c(axis.x), c(axis.y), c(axis.z)];
    let (mut ratio, mut along) = (0.0f64, f64::INFINITY);
    for v in m {
        let a = lower(&v[0].mul(&d[0]).add(&v[1].mul(&d[1])).add(&v[2].mul(&d[2])));
        if a.is_nan() || a <= 0.0 {
            return None;
        }
        ratio = ratio.max(up(c(norm_upper(&cross3(v, &d))).div(&c(a))?));
        along = along.min(a);
    }
    let spread = upper(&Fast::atan2(&c(ratio), &c(1.0))?);
    let least = lower(&c(along).div(&c(length))?).max(0.0);
    Some(Cone {
        axis,
        spread,
        least,
        factor,
    })
}

/// The cells of a nonperiodic spline surface over its domain, a tensor grid.
#[derive(Debug, Clone)]
pub(super) struct SurfaceCells {
    ubreaks: Breaks,
    vbreaks: Breaks,
    /// Outward binary64 bounds of the breaks.
    uf: Vec<(f64, f64)>,
    vf: Vec<(f64, f64)>,
    cells: Vec<SurfaceCell>,
    degrees: [usize; 2],
    domain: [f64; 4],
    /// The domain corners where `S_u × S_v` vanishes exactly.
    corners: Vec<Corner>,
}

/// Halvings of a singular corner's patch toward the corner: its Taylor data
/// on the patch and on each nested box (deeper boxes bound the third
/// derivatives more tightly, until rounding takes over).
const CORNER_LEVELS: usize = 3;

/// A domain corner `c` where the parametric normal vanishes: its boundary
/// rows leave it in one direction (a cusp), and the normal's limit depends
/// on the direction of approach. On boxes `[c, c ± 2^-k L]` of its patch,
/// in each box's own parameters `y` from `c` (`M` the patch's positive
/// multiple of `S_u × S_v`, exactly zero at `c`),
/// `M(c + y) = G y + H[y, y] / 2 + E`, `G = [M_s M_t](c)`, `H` the second
/// derivatives at `c` and `|E| <= T(y) / 6`, `T` the third derivatives'
/// bounds on the box at `|y_s|`, `|y_t|`.
#[derive(Debug, Clone)]
struct Corner {
    at: Point2,
    levels: Vec<CornerBox>,
}

#[derive(Debug, Clone)]
struct CornerBox {
    /// The box, inner binary64 bounds `[u0, u1, v0, v1]`.
    inner: [f64; 4],
    lengths: [Fast; 2],
    /// `M_s`, `M_t` at the corner.
    g: [V3; 2],
    /// `M_ss`, `M_st`, `M_tt` at the corner.
    h: [V3; 3],
    /// Bounds of `|M_sss|`, `|M_sst|`, `|M_stt|`, `|M_ttt|` on the box.
    third: [f64; 4],
}

impl CornerBox {
    /// The normal's turn over a parameter triangle `p` inside the box, the
    /// corner excluded. With `Y` the vertices' offsets from `c` but a zero
    /// one and `h(y) = G y + H[y, y] / 2`, every `M(c + y)` for `y` in `Y`'s
    /// hull lies within `|E| <= T / 6` of `h(y)`, and `h` over the hull
    /// within the convex cone of its quadratic Bézier controls, the blossoms
    /// `h[y_i, y_j]`. A triangle at the corner is its rays `c + r y`,
    /// `0 < r <= 1`: `M / r` lies within `r² T / 6` of `G y + r H[y, y] / 2`,
    /// in the cone of the `G y_i` and the `h[y_i, y_j]`. So every normal is
    /// within `asin ε` of the cone, `ε = T / (6 μ)`, `T` at the largest
    /// `|y_s|`, `|y_t|` and `μ` the controls' least component along their
    /// mean direction; the controls pairwise within 90°, the cone's
    /// directions span their widest pairwise angle, and the normal turns by
    /// at most that plus `2 asin ε <= π ε` (`ε < 1`).
    fn turn(&self, at: Point2, p: [Point2; 3]) -> f64 {
        let mut offsets: Vec<[Fast; 2]> = Vec::with_capacity(3);
        let mut at_corner = false;
        let (mut reach_s, mut reach_t) = (0.0f64, 0.0f64);
        let abs = |x: &Fast| {
            let (lo, hi) = x.bounds_f64();
            (-lo).max(hi)
        };
        for q in p {
            if q.x == at.x && q.y == at.y {
                at_corner = true;
                continue;
            }
            let (Some(s), Some(t)) = (
                c(q.x).sub(&c(at.x)).div(&self.lengths[0]),
                c(q.y).sub(&c(at.y)).div(&self.lengths[1]),
            ) else {
                return f64::INFINITY;
            };
            reach_s = reach_s.max(abs(&s));
            reach_t = reach_t.max(abs(&t));
            offsets.push([s, t]);
        }
        let linear = |y: &[Fast; 2]| -> V3 {
            std::array::from_fn(|k| self.g[0][k].mul(&y[0]).add(&self.g[1][k].mul(&y[1])))
        };
        // H[y, z] and the blossom h[y, z] = G (y + z) / 2 + H[y, z] / 2.
        let form = |y: &[Fast; 2], z: &[Fast; 2]| -> V3 {
            let mixed = y[0].mul(&z[1]).add(&y[1].mul(&z[0]));
            std::array::from_fn(|k| {
                self.h[0][k]
                    .mul(&y[0].mul(&z[0]))
                    .add(&self.h[1][k].mul(&mixed))
                    .add(&self.h[2][k].mul(&y[1].mul(&z[1])))
            })
        };
        let half = c(0.5);
        let blossom = |y: &[Fast; 2], z: &[Fast; 2]| -> V3 {
            let (a, b, q) = (linear(y), linear(z), form(y, z));
            std::array::from_fn(|k| a[k].add(&b[k]).add(&q[k]).mul(&half))
        };
        let mut images: Vec<V3> = Vec::with_capacity(6);
        for (i, y) in offsets.iter().enumerate() {
            if at_corner {
                images.push(linear(y));
            }
            for z in &offsets[i..] {
                images.push(blossom(y, z));
            }
        }
        if images.is_empty() {
            return f64::INFINITY;
        }
        // The widest pairwise angle, each certainly below 90°: every
        // direction of their convex cone lies within it of every other. Its
        // tangent `|x × y| / x · y` bounded first, one arctangent after.
        let mut tangent = 0.0f64;
        for i in 0..images.len() {
            for j in i + 1..images.len() {
                let (x, y) = (&images[i], &images[j]);
                let dot = lower(&x[0].mul(&y[0]).add(&x[1].mul(&y[1])).add(&x[2].mul(&y[2])));
                if dot.is_nan() || dot <= 0.0 {
                    return f64::INFINITY;
                }
                let Some(ratio) = c(norm_upper(&cross3(x, y))).div(&c(dot)) else {
                    return f64::INFINITY;
                };
                tangent = tangent.max(up(ratio));
            }
        }
        let Some(spread) = Fast::atan2(&c(tangent), &c(1.0)).map(|x| upper(&x)) else {
            return f64::INFINITY;
        };
        // μ: the images' least component along their mean direction.
        let mid = |x: &Fast| {
            let (lo, hi) = x.bounds_f64();
            0.5 * lo + 0.5 * hi
        };
        let mut sum = Vec3::new(0.0, 0.0, 0.0);
        for x in &images {
            let v = Vec3::new(mid(&x[0]), mid(&x[1]), mid(&x[2]));
            let l = v.length();
            if !(l.is_finite() && l > 0.0) {
                return f64::INFINITY;
            }
            sum = sum + v * (1.0 / l);
        }
        let d = [c(sum.x), c(sum.y), c(sum.z)];
        let length = c(norm_upper(&d));
        let mut least = f64::INFINITY;
        for x in &images {
            let along = x[0].mul(&d[0]).add(&x[1].mul(&d[1])).add(&x[2].mul(&d[2]));
            least = least.min(c(lower(&along)).div(&length).map_or(0.0, |x| lower(&x)));
        }
        if least.is_nan() || least <= 0.0 {
            return f64::INFINITY;
        }
        let [a, b, e, f] = self.third;
        let (s, t) = (c(reach_s), c(reach_t));
        let three = c(3.0);
        let cubic = c(a)
            .mul(&s.square().mul(&s))
            .add(&three.mul(&c(b)).mul(&s.square().mul(&t)))
            .add(&three.mul(&c(e)).mul(&s.mul(&t.square())))
            .add(&c(f).mul(&t.square().mul(&t)));
        let Some(epsilon) = cubic.div(&c(6.0).mul(&c(least))) else {
            return f64::INFINITY;
        };
        if upper(&epsilon).is_nan() || upper(&epsilon) >= 1.0 {
            return f64::INFINITY;
        }
        // 2 asin ε <= π ε on [0, 1]; π rounded up.
        up(c(spread).add(&c(3.141_592_653_589_794).mul(&epsilon)))
    }
}

/// The corners of a patch grid's domain where the exact controls make
/// `S_u × S_v` vanish: there `S_u` is a nonzero multiple of the corner
/// control's difference to its neighbour along `u` (rational or not),
/// likewise `S_v`.
fn singular_corners(
    patches: &[crate::surface::ExactBezierSurface3],
    nu: usize,
    nv: usize,
) -> Vec<Corner> {
    let mut out = Vec::new();
    for (pi, a) in [(0, 0usize), (nu - 1, 1)] {
        for (pj, b) in [(0, 0usize), (nv - 1, 1)] {
            if let Some(corner) = singular_corner(&patches[pi * nv + pj], a, b) {
                out.push(corner);
            }
        }
    }
    out
}

fn singular_corner(
    patch: &crate::surface::ExactBezierSurface3,
    a: usize,
    b: usize,
) -> Option<Corner> {
    let [p, q] = patch.degrees();
    // Rational `M = N_u × N_v` has degrees `(4p - 1, 4q - 1)`: binomials
    // stay exact in binary64 up to 56.
    if p == 0 || q == 0 || 4 * p > 57 || 4 * q > 57 {
        return None;
    }
    let controls = patch.homogeneous_poles();
    let point = |i: usize, j: usize| -> [R; 3] {
        let h = &controls[i * (q + 1) + j];
        [&h[0] / &h[3], &h[1] / &h[3], &h[2] / &h[3]]
    };
    let (ic, jc) = (a * p, b * q);
    let centre = point(ic, jc);
    let along_u = point(if a == 0 { 1 } else { p - 1 }, jc);
    let along_v = point(ic, if b == 0 { 1 } else { q - 1 });
    let du: [R; 3] = std::array::from_fn(|k| &along_u[k] - &centre[k]);
    let dv: [R; 3] = std::array::from_fn(|k| &along_v[k] - &centre[k]);
    let zero = R::from_integer(0.into());
    let vanishes = (0..3).all(|k| {
        let (m, n) = ((k + 1) % 3, (k + 2) % 3);
        &du[m] * &dv[n] - &du[n] * &dv[m] == zero
    });
    if !vanishes {
        return None;
    }
    let exact = |x: &R| {
        let (l, h) = Fast::from_r(x).bounds_f64();
        (l == h).then_some(l)
    };
    let [[u0, u1], [v0, v1]] = patch.domain().clone();
    let at = Point2::new(
        exact(if a == 0 { &u0 } else { &u1 })?,
        exact(if b == 0 { &v0 } else { &v1 })?,
    );
    let two = R::from_integer(2.into());
    let (mut u, mut v) = ([u0, u1], [v0, v1]);
    let mut net: Vec<H> = controls.iter().map(H::exact).collect();
    let mut levels = Vec::with_capacity(CORNER_LEVELS + 1);
    for level in 0..=CORNER_LEVELS {
        if level > 0 {
            let (l, r) = halve(&net, p, q, 0);
            let (l, r) = halve(if a == 0 { &l } else { &r }, p, q, 1);
            net = if b == 0 { l } else { r };
            let (um, vm) = ((&u[0] + &u[1]) / &two, (&v[0] + &v[1]) / &two);
            u[1 - a] = um;
            v[1 - b] = vm;
        }
        match corner_box(&net, p, q, a, b, &u, &v) {
            Some(x) => levels.push(x),
            None => break,
        }
    }
    (!levels.is_empty()).then_some(Corner { at, levels })
}

fn corner_box(
    net: &[H],
    p: usize,
    q: usize,
    a: usize,
    b: usize,
    u: &[R; 2],
    v: &[R; 2],
) -> Option<CornerBox> {
    let (_, wq, w, _, _) = translated(net)?;
    let poly: Poly3 = std::array::from_fn(|k| Poly::new([p, q], wq.iter().map(|x| x[k]).collect()));
    let (ps, pt) = (poly::derivative3(&poly, 0), poly::derivative3(&poly, 1));
    let m = if w.iter().any(|x| x != &w[0]) {
        let ww = Poly::new([p, q], w);
        let (ws, wt) = (ww.derivative(0), ww.derivative(1));
        let ns = poly::sub33(&poly::mul3(&ps, &ww), &poly::mul3(&poly, &ws));
        let nt = poly::sub33(&poly::mul3(&pt, &ww), &poly::mul3(&poly, &wt));
        poly::cross3(&ns, &nt)
    } else {
        poly::cross3(&ps, &pt)
    };
    let d = |x: &Poly3, axis: usize| poly::derivative3(x, axis);
    let (ms, mt) = (d(&m, 0), d(&m, 1));
    let (mss, mst, mtt) = (d(&ms, 0), d(&ms, 1), d(&mt, 1));
    let third = [
        poly::largest3(&d(&mss, 0)),
        poly::largest3(&d(&mss, 1)),
        poly::largest3(&d(&mtt, 0)),
        poly::largest3(&d(&mtt, 1)),
    ];
    if !third.iter().all(|x| x.is_finite()) {
        return None;
    }
    let (s, t) = (c(a as f64), c(b as f64));
    let value = |x: &Poly3| -> V3 { std::array::from_fn(|k| x[k].value(&s, &t)) };
    // Inner bounds: the lower end's upper, the upper end's lower.
    let (lo, hi) = (
        |x: &R| Fast::from_r(x).bounds_f64().1,
        |x: &R| Fast::from_r(x).bounds_f64().0,
    );
    Some(CornerBox {
        inner: [lo(&u[0]), hi(&u[1]), lo(&v[0]), hi(&v[1])],
        lengths: [
            Fast::from_r(&(&u[1] - &u[0])),
            Fast::from_r(&(&v[1] - &v[0])),
        ],
        g: [value(&ms), value(&mt)],
        h: [value(&mss), value(&mst), value(&mtt)],
        third,
    })
}

/// Coefficients over a parameter box: `[a, b, c, |S_u|, |S_v|]` bounds.
pub(super) type Coefficients = [f64; 5];

impl SurfaceCells {
    pub(super) fn new(surface: &BSplineSurface3) -> Result<Self> {
        let patches = surface.bezier_patches()?;
        let first = patches
            .first()
            .ok_or(Error::ComputationLimit("tessellation spline cells"))?;
        let [p, q] = first.degrees();
        let nv = patches
            .iter()
            .take_while(|x| x.domain()[0] == first.domain()[0])
            .count();
        let nu = patches.len() / nv;
        let controls = (p + 1) * (q + 1);
        let mut depth = SURFACE_DEPTH;
        while depth > 0
            && (patches.len() << (2 * depth) > MAX_SURFACE_CELLS
                || (patches.len() << (2 * depth)) * controls > MAX_SURFACE_CONTROLS)
        {
            depth -= 1;
        }
        let split = 1usize << depth;
        let two = R::from_integer(2.into());
        let breaks = |ends: Vec<[R; 2]>| {
            let mut out = Vec::new();
            for [a, b] in ends {
                let mut row = vec![a, b];
                for _ in 0..depth {
                    let mut next = Vec::with_capacity(2 * row.len());
                    for w in row.windows(2) {
                        next.push(w[0].clone());
                        next.push((&w[0] + &w[1]) / &two);
                    }
                    next.push(row[row.len() - 1].clone());
                    row = next;
                }
                if out.is_empty() {
                    out = row;
                } else {
                    out.extend(row.into_iter().skip(1));
                }
            }
            out
        };
        let ubreaks = breaks(
            (0..nu)
                .map(|i| patches[i * nv].domain()[0].clone())
                .collect(),
        );
        let vbreaks = breaks((0..nv).map(|j| patches[j].domain()[1].clone()).collect());
        let bounds = |x: &R| {
            let f = Fast::from_r(x).bounds_f64();
            (f.0, f.1)
        };
        let uf: Vec<(f64, f64)> = ubreaks.iter().map(bounds).collect();
        let vf: Vec<(f64, f64)> = vbreaks.iter().map(bounds).collect();
        let (mu, mv) = (nu * split, nv * split);
        let mut cells: Vec<Option<SurfaceCell>> = vec![None; mu * mv];
        for (index, patch) in patches.iter().enumerate() {
            let (pi, pj) = (index / nv, index % nv);
            let net: Vec<H> = patch.homogeneous_poles().iter().map(H::exact).collect();
            // Halve along u, then each piece along v.
            let mut rows = vec![net];
            for _ in 0..depth {
                rows = rows
                    .iter()
                    .flat_map(|net| {
                        let (l, r) = halve(net, p, q, 0);
                        [l, r]
                    })
                    .collect();
            }
            for (a, net) in rows.into_iter().enumerate() {
                let mut pieces = vec![net];
                for _ in 0..depth {
                    pieces = pieces
                        .iter()
                        .flat_map(|net| {
                            let (l, r) = halve(net, p, q, 1);
                            [l, r]
                        })
                        .collect();
                }
                for (b, net) in pieces.into_iter().enumerate() {
                    let (i, j) = (pi * split + a, pj * split + b);
                    let lengths = [
                        Fast::from_r(&(&ubreaks[i + 1] - &ubreaks[i])),
                        Fast::from_r(&(&vbreaks[j + 1] - &vbreaks[j])),
                    ];
                    cells[i * mv + j] = Some(surface_cell(lengths, &net, p, q));
                }
            }
        }
        let cells: Vec<SurfaceCell> = cells.into_iter().map(|x| x.expect("every cell")).collect();
        let (ubreaks, vbreaks) = (Breaks::new(ubreaks), Breaks::new(vbreaks));
        let ((ua, ub), (va, vb)) = surface.domain();
        Ok(Self {
            ubreaks,
            vbreaks,
            uf,
            vf,
            cells,
            degrees: [p, q],
            domain: [ua, ub, va, vb],
            corners: singular_corners(&patches, nu, nv),
        })
    }

    /// A parameter point moved into the domain.
    pub(super) fn clamp(&self, uv: Point2) -> Point2 {
        let [a, b, c, d] = self.domain;
        Point2::new(uv.x.clamp(a, b), uv.y.clamp(c, d))
    }

    fn cell_at(&self, uv: Point2) -> Option<(&SurfaceCell, Fast, Fast)> {
        if !(uv.x.is_finite() && uv.y.is_finite()) {
            return None;
        }
        let (i, j) = (self.ubreaks.find(uv.x), self.vbreaks.find(uv.y));
        let cell = &self.cells[i * (self.vbreaks.exact.len() - 1) + j];
        let s = self.ubreaks.local(i, uv.x, &cell.lengths[0])?;
        let t = self.vbreaks.local(j, uv.y, &cell.lengths[1])?;
        Some((cell, s, t))
    }

    /// `S`, `S_u` and `S_v` at a parameter point in the domain, enclosed.
    pub(super) fn jet(&self, uv: Point2) -> Option<[V3; 3]> {
        let (cell, s, t) = self.cell_at(uv)?;
        if cell.wq.is_empty() {
            return None;
        }
        let [p, q] = self.degrees;
        let (ku, _) = scales(p, &cell.lengths[0]);
        let (kv, _) = scales(q, &cell.lengths[1]);
        // Each u-row along v: its value and v-derivative.
        let mut a = Vec::with_capacity(p + 1);
        let mut av = Vec::with_capacity(p + 1);
        let mut w = Vec::with_capacity(p + 1);
        let mut wv = Vec::with_capacity(p + 1);
        for i in 0..=p {
            let row = &cell.wq[i * (q + 1)..(i + 1) * (q + 1)];
            let wr = &cell.w[i * (q + 1)..(i + 1) * (q + 1)];
            a.push(value(row, &t, mix3));
            w.push(value(wr, &t, mix1));
            let d: Vec<V3> = row
                .windows(2)
                .map(|x| scale3(&sub3(&x[1], &x[0]), &kv))
                .collect();
            let dw: Vec<Fast> = wr.windows(2).map(|x| x[1].sub(&x[0]).mul(&kv)).collect();
            av.push(value(&d, &t, mix3));
            wv.push(value(&dw, &t, mix1));
        }
        let du: Vec<V3> = a
            .windows(2)
            .map(|x| scale3(&sub3(&x[1], &x[0]), &ku))
            .collect();
        let dwu: Vec<Fast> = w.windows(2).map(|x| x[1].sub(&x[0]).mul(&ku)).collect();
        let s0 = value(&a, &s, mix3);
        let w0 = value(&w, &s, mix1);
        let su = value(&du, &s, mix3);
        let wu = value(&dwu, &s, mix1);
        let sv = value(&av, &s, mix3);
        let wvv = value(&wv, &s, mix1);
        let rel: V3 = [s0[0].div(&w0)?, s0[1].div(&w0)?, s0[2].div(&w0)?];
        let point: V3 = std::array::from_fn(|k| c(cell.origin[k]).add(&rel[k]));
        let d_u = [
            su[0].sub(&wu.mul(&rel[0])).div(&w0)?,
            su[1].sub(&wu.mul(&rel[1])).div(&w0)?,
            su[2].sub(&wu.mul(&rel[2])).div(&w0)?,
        ];
        let d_v = [
            sv[0].sub(&wvv.mul(&rel[0])).div(&w0)?,
            sv[1].sub(&wvv.mul(&rel[1])).div(&w0)?,
            sv[2].sub(&wvv.mul(&rel[2])).div(&w0)?,
        ];
        Some([point, d_u, d_v])
    }

    /// The lengths of `S_u` and `S_v` at a parameter point, about (1 where
    /// they are not finite and positive).
    pub(super) fn speeds(&self, uv: Point2) -> (f64, f64) {
        let length = |v: &V3| {
            let m = v.map(|x| {
                let (lo, hi) = x.bounds_f64();
                0.5 * lo + 0.5 * hi
            });
            let l = Vec3::new(m[0], m[1], m[2]).length();
            if l.is_finite() && l > 0.0 {
                l
            } else {
                1.0
            }
        };
        self.jet(uv)
            .map_or((1.0, 1.0), |j| (length(&j[1]), length(&j[2])))
    }

    /// The cells meeting a parameter box: `(i0, i1, j0, j1)`, inclusive.
    fn range(&self, u0: f64, u1: f64, v0: f64, v1: f64) -> (usize, usize, usize, usize) {
        let range = |f: &[(f64, f64)], lo: f64, hi: f64| {
            let n = f.len() - 1;
            let first = (0..n).find(|k| f[k + 1].1 >= lo).unwrap_or(n - 1);
            let last = (0..n).rev().find(|k| f[*k].0 <= hi).unwrap_or(0);
            (first, last.max(first))
        };
        let (i0, i1) = range(&self.uf, u0, u1);
        let (j0, j1) = range(&self.vf, v0, v1);
        (i0, i1, j0, j1)
    }

    /// `S_u × S_v` at a parameter point, about.
    pub(super) fn normal(&self, uv: Point2) -> Option<Vec3> {
        let [_, su, sv] = self.jet(uv)?;
        let m = cross3(&su, &sv).map(|x| {
            let (lo, hi) = x.bounds_f64();
            0.5 * lo + 0.5 * hi
        });
        Some(Vec3::new(m[0], m[1], m[2]))
    }

    /// The largest coefficients over the cells meeting a parameter box.
    pub(super) fn coefficients(&self, u0: f64, u1: f64, v0: f64, v1: f64) -> Coefficients {
        let (i0, i1, j0, j1) = self.range(u0, u1, v0, v1);
        let nv = self.vbreaks.exact.len() - 1;
        let mut out = [0.0f64; 5];
        for i in i0..=i1 {
            for j in j0..=j1 {
                let b = &self.cells[i * nv + j].bounds;
                for k in 0..5 {
                    out[k] = out[k].max(b[k]);
                }
            }
        }
        out
    }

    /// The certified deviation (without node gaps) and normal turn of a
    /// parameter triangle, and the parametric normal `S_u × S_v` at its box's
    /// centre: T-a's `(a U² + 2 b U V + c V²) / 8` over the cells meeting the
    /// box, and `(M_u U + M_v V) / μ` with `M_u = a D_v + D_u b`,
    /// `M_v = b D_v + D_u c` and `μ` the least `|S_u × S_v|` over the box;
    /// inside a singular corner's patch also its corner turn.
    pub(super) fn triangle_bound(&self, p: [Point2; 3]) -> (f64, f64, Option<Vec3>) {
        let span = |f: fn(&Point2) -> f64| {
            p.iter()
                .map(f)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                    (l.min(x), h.max(x))
                })
        };
        let (u0, u1) = span(|q| q.x);
        let (v0, v1) = span(|q| q.y);
        let (deviation, turn, normal) = self.box_bound(u0, u1, v0, v1);
        (deviation, turn.min(self.corner_turn(p)), normal)
    }

    /// The normal's turn over a parameter triangle (or segment, a vertex
    /// repeated) inside a singular corner's boxes, the corner itself (where
    /// the normal is undefined) excluded: the least of `CornerBox::turn`
    /// over the boxes containing it; infinite outside them.
    fn corner_turn(&self, p: [Point2; 3]) -> f64 {
        let mut best = f64::INFINITY;
        for corner in &self.corners {
            for level in &corner.levels {
                let [u0, u1, v0, v1] = level.inner;
                if !p
                    .iter()
                    .all(|q| q.x >= u0 && q.x <= u1 && q.y >= v0 && q.y <= v1)
                {
                    break;
                }
                best = best.min(level.turn(corner.at, p));
            }
        }
        best
    }

    /// `triangle_bound` of a parameter box.
    pub(super) fn box_bound(&self, u0: f64, u1: f64, v0: f64, v1: f64) -> (f64, f64, Option<Vec3>) {
        let [a, b, cc, du, dv] = self.coefficients(u0, u1, v0, v1);
        let (fu, fv) = (c(u1).sub(&c(u0)), c(v1).sub(&c(v0)));
        let deviation = up(c(a)
            .mul(&fu.square())
            .add(&c(2.0 * b).mul(&fu.mul(&fv)))
            .add(&c(cc).mul(&fv.square()))
            .mul(&c(0.125)));
        let mu_u = c(a).mul(&c(dv)).add(&c(du).mul(&c(b)));
        let mu_v = c(b).mul(&c(dv)).add(&c(du).mul(&c(cc)));
        let centre = Point2::new(0.5 * u0 + 0.5 * u1, 0.5 * v0 + 0.5 * v1);
        let reach_u = c(centre.x).sub(&c(u0)).union(&c(u1).sub(&c(centre.x)));
        let reach_v = c(centre.y).sub(&c(v0)).union(&c(v1).sub(&c(centre.y)));
        let Some([_, su, sv]) = self.jet(centre) else {
            return (deviation, f64::INFINITY, None);
        };
        let m = cross3(&su, &sv);
        let mid = |x: &Fast| {
            let (lo, hi) = x.bounds_f64();
            0.5 * lo + 0.5 * hi
        };
        let normal = Vec3::new(mid(&m[0]), mid(&m[1]), mid(&m[2]));
        let (cone, cone_least, rated) = self.cone(u0, u1, v0, v1, normal, &fu, &fv);
        let least = lower(&c(norm_lower(&m)).sub(&mu_u.mul(&reach_u).add(&mu_v.mul(&reach_v))))
            .max(cone_least);
        let turn = if least > 0.0 {
            up(mu_u
                .mul(&fu)
                .add(&mu_v.mul(&fv))
                .div(&c(least))
                .unwrap_or(c(f64::INFINITY)))
        } else {
            f64::INFINITY
        };
        (deviation, turn.min(cone).min(rated), Some(normal))
    }

    /// Over the cells meeting a box of extents `U`, `V`: twice the largest
    /// angle between `d` and a normal of a cell (at most the angle between
    /// `d` and the cell's cone axis plus its spread), so any two normals'
    /// angle; the least `|S_u × S_v|` of a cell; and the largest turn
    /// `(|M × M_u| U + |M × M_v| V) / |M|²` of a cell (`|N_u| = |M × M_u| /
    /// |M|²` pointwise for the cell's own `M`); infinite and zero when a
    /// cell has no cone or leaves 90° of `d`.
    #[allow(clippy::too_many_arguments)]
    fn cone(
        &self,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        d: Vec3,
        fu: &Fast,
        fv: &Fast,
    ) -> (f64, f64, f64) {
        let none = (f64::INFINITY, 0.0, f64::INFINITY);
        let d = [c(d.x), c(d.y), c(d.z)];
        let (i0, i1, j0, j1) = self.range(u0, u1, v0, v1);
        let nv = self.vbreaks.exact.len() - 1;
        let (mut angle, mut least, mut rated) = (0.0f64, f64::INFINITY, 0.0f64);
        for i in i0..=i1 {
            for j in j0..=j1 {
                let cell = &self.cells[i * nv + j];
                let Some(cone) = &cell.cone else {
                    return none;
                };
                let a = [c(cone.axis.x), c(cone.axis.y), c(cone.axis.z)];
                let along = lower(&a[0].mul(&d[0]).add(&a[1].mul(&d[1])).add(&a[2].mul(&d[2])));
                let between = if along > 0.0 {
                    Fast::atan2(&c(norm_upper(&cross3(&a, &d))), &c(along))
                        .map_or(f64::INFINITY, |x| upper(&x))
                } else {
                    f64::INFINITY
                };
                angle = angle.max(up(c(between).add(&c(cone.spread))));
                least = least.min(lower(&c(cone.least).mul(&c(cone.factor))));
                rated = rated.max(match cell.rates {
                    Some([ru, rv]) if cone.least > 0.0 => c(ru)
                        .mul(fu)
                        .add(&c(rv).mul(fv))
                        .div(&c(cone.least).square())
                        .map_or(f64::INFINITY, up),
                    _ => f64::INFINITY,
                });
            }
        }
        let angle = if angle < std::f64::consts::FRAC_PI_2 {
            up(c(angle).mul(&c(2.0)))
        } else {
            f64::INFINITY
        };
        (angle, least.max(0.0), rated)
    }
}

/// De Casteljau's halves of a `(p + 1) × (q + 1)` net along u (axis 0) or
/// v (axis 1).
fn halve(net: &[H], p: usize, q: usize, axis: usize) -> (Vec<H>, Vec<H>) {
    let (mut left, mut right) = (net.to_vec(), net.to_vec());
    if axis == 0 {
        for j in 0..=q {
            let column: Vec<H> = (0..=p).map(|i| net[i * (q + 1) + j]).collect();
            let (l, r) = halves(&column);
            for i in 0..=p {
                left[i * (q + 1) + j] = l[i];
                right[i * (q + 1) + j] = r[i];
            }
        }
    } else {
        for i in 0..=p {
            let (l, r) = halves(&net[i * (q + 1)..(i + 1) * (q + 1)]);
            left[i * (q + 1)..(i + 1) * (q + 1)].copy_from_slice(&l);
            right[i * (q + 1)..(i + 1) * (q + 1)].copy_from_slice(&r);
        }
    }
    (left, right)
}

fn surface_cell(lengths: [Fast; 2], net: &[H], p: usize, q: usize) -> SurfaceCell {
    let Some((origin, wq, w, reach, least)) = translated(net) else {
        return SurfaceCell {
            lengths,
            origin: [0.0; 3],
            wq: Vec::new(),
            w: Vec::new(),
            bounds: [f64::INFINITY; 5],
            cone: None,
            rates: None,
        };
    };
    let rational = w.iter().any(|x| x != &w[0]);
    // M in the cell's own parameters is L_u L_v times the surface's
    // M = w^e S_u × S_v, the same function on every cell: its coefficients
    // and rates are divided by k = L_u L_v (and its derivative by L_u or
    // L_v), so that cells of different sizes compare.
    let k = lengths[0].mul(&lengths[1]);
    let rates = |m: &Poly3| -> Option<[f64; 2]> {
        let limit = MAX_RATE_DEGREE[usize::from(rational)];
        if p > limit || q > limit {
            return None;
        }
        let per = |axis: usize, length: &Fast| {
            let r = poly::largest3(&poly::cross3(m, &poly::derivative3(m, axis)));
            c(r).div(&length.mul(&k.square())).map_or(f64::INFINITY, up)
        };
        Some([per(0, &lengths[0]), per(1, &lengths[1])])
    };
    // The normals' cone, from M's coefficients about its middle value.
    let cone = |m: &Poly3, factor: f64| -> Option<Cone> {
        let [x, y, z] = m.each_ref().map(|p| p.coefficients());
        let vectors: Vec<V3> = (0..x.len())
            .map(|i| Some([x[i].div(&k)?, y[i].div(&k)?, z[i].div(&k)?]))
            .collect::<Option<_>>()?;
        let half = c(0.5);
        let centre = [
            m[0].value(&half, &half).div(&k)?,
            m[1].value(&half, &half).div(&k)?,
            m[2].value(&half, &half).div(&k)?,
        ];
        cone_of(&vectors, &centre, factor)
    };
    let a: Poly3 = std::array::from_fn(|k| Poly::new([p, q], wq.iter().map(|x| x[k]).collect()));
    let ww = Poly::new([p, q], w.clone());
    let w_max = w.iter().map(upper).fold(0.0, f64::max);
    // A lower bound of 1 / w^e: |S_u × S_v| / |M|.
    let factor = |e: i32| {
        let mut x = c(1.0);
        for _ in 0..e {
            x = x.mul(&c(w_max));
        }
        c(1.0).div(&x).map_or(0.0, |x| lower(&x).max(0.0))
    };
    if rational && p <= MAX_RATIONAL_SURFACE_DEGREE && q <= MAX_RATIONAL_SURFACE_DEGREE {
        // The quotient's partials through their numerators' coefficients:
        // S_u = N_u / w², S_uu = N_uu / w³, S_uv = N_uv / w³, S_vv = N_vv / w³.
        let (au, av) = (poly::derivative3(&a, 0), poly::derivative3(&a, 1));
        let (wu, wv) = (ww.derivative(0), ww.derivative(1));
        let (auu, avv, auv) = (
            poly::derivative3(&au, 0),
            poly::derivative3(&av, 1),
            poly::derivative3(&au, 1),
        );
        let (wuu, wvv, wuv) = (wu.derivative(0), wv.derivative(1), wu.derivative(1));
        let times = poly::mul3;
        let nu = poly::sub33(&times(&au, &ww), &times(&a, &wu));
        let nv = poly::sub33(&times(&av, &ww), &times(&a, &wv));
        let nuu = poly::sub33(
            &times(&poly::sub33(&times(&auu, &ww), &times(&a, &wuu)), &ww),
            &poly::scale3(&times(&nu, &wu), 2.0),
        );
        let nvv = poly::sub33(
            &times(&poly::sub33(&times(&avv, &ww), &times(&a, &wvv)), &ww),
            &poly::scale3(&times(&nv, &wv), 2.0),
        );
        let mixed = poly::sub33(
            &poly::add33(&times(&auv, &ww), &times(&au, &wv)),
            &poly::add33(&times(&av, &wu), &times(&a, &wuv)),
        );
        let nuv = poly::sub33(&times(&mixed, &ww), &poly::scale3(&times(&nu, &wv), 2.0));
        let least = c(least);
        let (w2, w3) = (least.square(), least.square().mul(&least));
        let bound = |n: &Poly3, w: &Fast, l: Fast| {
            c(poly::largest3(n))
                .div(w)
                .and_then(|x| x.div(&l))
                .map_or(f64::INFINITY, up)
        };
        let [lu, lv] = lengths;
        let bounds = [
            bound(&nuu, &w3, lu.square()),
            bound(&nuv, &w3, lu.mul(&lv)),
            bound(&nvv, &w3, lv.square()),
            bound(&nu, &w2, lu),
            bound(&nv, &w2, lv),
        ];
        let m = poly::cross3(&nu, &nv);
        return SurfaceCell {
            lengths,
            origin,
            wq,
            w,
            bounds,
            rates: rates(&m),
            cone: cone(&m, factor(4)),
        };
    }
    let (cone, rate) = if !rational && p <= MAX_CONE_DEGREE && q <= MAX_CONE_DEGREE {
        let m = poly::cross3(&poly::derivative3(&a, 0), &poly::derivative3(&a, 1));
        (cone(&m, factor(2)), rates(&m))
    } else {
        (None, None)
    };
    let row = |i: usize| i * (q + 1)..(i + 1) * (q + 1);
    let column = |j: usize| (0..=p).map(move |i| i * (q + 1) + j);
    // Differences along u (columns), along v (rows), and mixed.
    let along_u = |k: usize| {
        let mut a = 0.0f64;
        let mut b = 0.0f64;
        for j in 0..=q {
            let x: Vec<V3> = column(j).map(|i| wq[i]).collect();
            let y: Vec<Fast> = column(j).map(|i| w[i]).collect();
            a = a.max(difference3(&x, k));
            b = b.max(difference1(&y, k));
        }
        (a, b)
    };
    let along_v = |k: usize| {
        let mut a = 0.0f64;
        let mut b = 0.0f64;
        for i in 0..=p {
            a = a.max(difference3(&wq[row(i)], k));
            b = b.max(difference1(&w[row(i)], k));
        }
        (a, b)
    };
    let (mut mixed, mut mixed_w) = (0.0f64, 0.0f64);
    for i in 0..p {
        for j in 0..q {
            let at = |i: usize, j: usize| i * (q + 1) + j;
            let d = sub3(
                &sub3(&wq[at(i + 1, j + 1)], &wq[at(i, j + 1)]),
                &sub3(&wq[at(i + 1, j)], &wq[at(i, j)]),
            );
            mixed = mixed.max(norm_upper(&d));
            let e = w[at(i + 1, j + 1)]
                .sub(&w[at(i, j + 1)])
                .sub(&w[at(i + 1, j)].sub(&w[at(i, j)]));
            let (lo, hi) = e.bounds_f64();
            mixed_w = mixed_w.max((-lo).max(hi));
        }
    }
    let (ku1, ku2) = scales(p, &lengths[0]);
    let (kv1, kv2) = scales(q, &lengths[1]);
    let kuv = c(p as f64)
        .mul(&c(q as f64))
        .div(&lengths[0].mul(&lengths[1]))
        .unwrap_or(c(f64::INFINITY));
    let ((au1, wu1), (au2, wu2)) = (along_u(1), along_u(2));
    let ((av1, wv1), (av2, wv2)) = (along_v(1), along_v(2));
    let (reach, least) = (c(reach), c(least));
    let a_u = ku1.mul(&c(au1));
    let a_v = kv1.mul(&c(av1));
    let a_uu = ku2.mul(&c(au2));
    let a_vv = kv2.mul(&c(av2));
    let a_uv = kuv.mul(&c(mixed));
    let w_u = ku1.mul(&c(wu1));
    let w_v = kv1.mul(&c(wv1));
    let w_uu = ku2.mul(&c(wu2));
    let w_vv = kv2.mul(&c(wv2));
    let w_uv = kuv.mul(&c(mixed_w));
    let finish = |x: Option<Fast>| x.map(up).unwrap_or(f64::INFINITY);
    let d_u = finish(a_u.add(&w_u.mul(&reach)).div(&least));
    let d_v = finish(a_v.add(&w_v.mul(&reach)).div(&least));
    let (fd_u, fd_v) = (c(d_u), c(d_v));
    let s_uu = finish(
        a_uu.add(&c(2.0).mul(&w_u).mul(&fd_u))
            .add(&w_uu.mul(&reach))
            .div(&least),
    );
    let s_uv = finish(
        a_uv.add(&w_u.mul(&fd_v))
            .add(&w_v.mul(&fd_u))
            .add(&w_uv.mul(&reach))
            .div(&least),
    );
    let s_vv = finish(
        a_vv.add(&c(2.0).mul(&w_v).mul(&fd_v))
            .add(&w_vv.mul(&reach))
            .div(&least),
    );
    SurfaceCell {
        lengths,
        origin,
        wq,
        w,
        bounds: [s_uu, s_uv, s_vv, d_u, d_v],
        cone,
        rates: rate,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::{DerivativeOrder, KnotSide};

    fn contains(x: &Fast, y: f64, slack: f64) -> bool {
        let (lo, hi) = x.bounds_f64();
        lo - slack <= y && y <= hi + slack
    }

    /// Rational and nonrational surfaces of degrees 2 and 5, 3 and 8, with
    /// interior knots: the cells' points and partials enclose the exact
    /// evaluation's, at points on and between knots, and the cells' bounds
    /// hold there.
    #[test]
    fn surface_cells_enclose_the_exact_jets() {
        for (p, q, rational) in [(2, 5, true), (2, 5, false), (3, 8, true), (2, 8, true)] {
            let u = crate::KnotVector::new(p, vec![0.0, 0.25, 1.0], vec![p + 1, 1, p + 1]).unwrap();
            let v = crate::KnotVector::new(q, vec![0.0, 0.5, 1.0], vec![q + 1, 2, q + 1]).unwrap();
            let (nu, nv) = (u.pole_count(), v.pole_count());
            let mut poles = Vec::new();
            let mut weights = Vec::new();
            for i in 0..nu {
                for j in 0..nv {
                    let (x, y) = (i as f64, j as f64 * 0.5);
                    poles.push(crate::Point3::new(
                        x + 0.1 * y,
                        y,
                        (x * 1.7 + y * 0.9).sin(),
                    ));
                    weights.push(if rational {
                        1.0 + 0.25 * ((i + 2 * j) % 3) as f64
                    } else {
                        1.0
                    });
                }
            }
            let s = BSplineSurface3::new(u, v, poles, Some(weights)).unwrap();
            let cells = SurfaceCells::new(&s).unwrap();
            for a in 0..=12 {
                for b in 0..=12 {
                    let uv = Point2::new(a as f64 / 12.0, b as f64 / 12.0);
                    let side = |x: f64| {
                        if x >= 1.0 {
                            KnotSide::Left
                        } else {
                            KnotSide::Right
                        }
                    };
                    let e = s
                        .evaluate(
                            uv.x,
                            uv.y,
                            DerivativeOrder::Second,
                            [side(uv.x), side(uv.y)],
                        )
                        .unwrap();
                    let [point, su, sv] = cells.jet(uv).unwrap();
                    let exact = e.position().to_array();
                    let d = |i, j, k: usize| e.derivative_bounds(i, j).unwrap()[k].representative();
                    for k in 0..3 {
                        assert!(contains(&point[k], exact[k], 1e-12), "{p} {q} {uv:?}");
                        assert!(contains(&su[k], d(1, 0, k), 1e-9), "{p} {q} {uv:?}");
                        assert!(contains(&sv[k], d(0, 1, k), 1e-9), "{p} {q} {uv:?}");
                    }
                    let [a2, b2, c2, du, dv] = cells.coefficients(uv.x, uv.x, uv.y, uv.y);
                    let norm = |i, j| (0..3).map(|k| d(i, j, k).powi(2)).sum::<f64>().sqrt();
                    let within = |x: f64, bound: f64| x <= bound * (1.0 + 1e-9) + 1e-12;
                    assert!(within(norm(2, 0), a2), "{p} {q} {uv:?}");
                    assert!(within(norm(1, 1), b2), "{p} {q} {uv:?}");
                    assert!(within(norm(0, 2), c2), "{p} {q} {uv:?}");
                    assert!(within(norm(1, 0), du), "{p} {q} {uv:?}");
                    assert!(within(norm(0, 1), dv), "{p} {q} {uv:?}");
                }
            }
        }
    }

    /// A cusp corner (both boundary rows leave `(1, 0)` along `-x`, so
    /// `S_u × S_v` vanishes there and the normal's limit depends on the
    /// direction): on triangles and segments at and near the corner, the
    /// corner turn bounds the widest angle between sampled exact normals,
    /// and shrinks with the triangles (rational or not).
    #[test]
    fn corner_turns_enclose_the_normals_at_a_cusp() {
        for rational in [false, true] {
            let (p, q) = (2, 3);
            let u = crate::KnotVector::new(p, vec![0.0, 1.0], vec![p + 1, p + 1]).unwrap();
            let v = crate::KnotVector::new(q, vec![0.0, 0.5, 1.0], vec![q + 1, 1, q + 1]).unwrap();
            let (nu, nv) = (u.pole_count(), v.pole_count());
            let mut poles = Vec::new();
            let mut weights = Vec::new();
            for i in 0..nu {
                for j in 0..nv {
                    let (x, y) = ((nu - 1 - i) as f64, j as f64);
                    poles.push(if i == nu - 1 {
                        crate::Point3::new(-0.3 * y, 0.0, 0.0)
                    } else if i == nu - 2 && j == 0 {
                        crate::Point3::new(1.0, 0.0, 0.0)
                    } else {
                        crate::Point3::new(x - 0.3 * y, 0.6 * x + 0.1 * y, 0.4 * x * y)
                    });
                    weights.push(if rational {
                        1.0 + 0.25 * ((i + 2 * j) % 3) as f64
                    } else {
                        1.0
                    });
                }
            }
            let s = BSplineSurface3::new(u, v, poles, Some(weights)).unwrap();
            let cells = SurfaceCells::new(&s).unwrap();
            assert_eq!(cells.corners.len(), 1);
            assert_eq!(cells.corners[0].at, Point2::new(1.0, 0.0));
            let normal = |uv: Point2| {
                let side = |x: f64| {
                    if x >= 1.0 {
                        KnotSide::Left
                    } else {
                        KnotSide::Right
                    }
                };
                let e = s
                    .evaluate(uv.x, uv.y, DerivativeOrder::First, [side(uv.x), side(uv.y)])
                    .unwrap();
                let d = |i, j| {
                    let b = e.derivative_bounds(i, j).unwrap();
                    Vec3::new(
                        b[0].representative(),
                        b[1].representative(),
                        b[2].representative(),
                    )
                };
                let m = d(1, 0).cross(d(0, 1));
                m * (1.0 / m.length())
            };
            let c0 = Point2::new(1.0, 0.0);
            let mut previous = f64::INFINITY;
            for h in [0.2, 0.05, 0.01] {
                // A fan at the corner, a segment along each boundary row from
                // it and triangles beside it.
                let shapes = [
                    [c0, Point2::new(1.0, h), Point2::new(1.0 - h, 0.0)],
                    [c0, Point2::new(1.0, h), Point2::new(1.0 - h, h)],
                    [c0, Point2::new(1.0, h), Point2::new(1.0, h)],
                    [c0, Point2::new(1.0 - h, 0.0), Point2::new(1.0 - h, 0.0)],
                    [
                        Point2::new(1.0, h),
                        Point2::new(1.0, 2.0 * h),
                        Point2::new(1.0 - h, 2.0 * h),
                    ],
                    [
                        Point2::new(1.0 - h, h),
                        Point2::new(1.0 - 2.0 * h, h),
                        Point2::new(1.0 - h, 0.0),
                    ],
                ];
                for (k, t) in shapes.iter().enumerate() {
                    let bound = cells.corner_turn(*t);
                    if h == 0.01 {
                        // Along a row from the corner the normal barely turns;
                        // a fan's triangles turn by their sector's spread.
                        let most = if k == 2 || k == 3 { 2.0 * h } else { 0.15 };
                        assert!(bound < most, "{rational} {k}: {bound}");
                    }
                    if k == 0 {
                        assert!(bound <= previous, "{rational} {h}: {bound}");
                        previous = bound;
                    }
                    if !bound.is_finite() {
                        continue;
                    }
                    let mut normals = Vec::new();
                    for i in 0..=8 {
                        for j in 0..=8 - i {
                            let (a, b) = (i as f64 / 8.0, j as f64 / 8.0);
                            // Toward the first vertex geometrically too.
                            for r in [1.0, 1e-3, 1e-6] {
                                let (a, b) = (r * a, r * b);
                                let w = 1.0 - a - b;
                                let x = Point2::new(
                                    t[0].x * w + t[1].x * a + t[2].x * b,
                                    t[0].y * w + t[1].y * a + t[2].y * b,
                                );
                                if x != c0 {
                                    normals.push(normal(x));
                                }
                            }
                        }
                    }
                    let mut widest = 0.0f64;
                    for x in &normals {
                        for y in &normals {
                            widest = widest.max(x.cross(*y).length().atan2(x.dot(*y)));
                        }
                    }
                    assert!(
                        widest <= bound + 1e-9,
                        "{rational} {h} {k}: {widest} > {bound}"
                    );
                }
            }
            // Away from the corner's patch: no corner turn.
            let far = [
                Point2::new(0.5, 0.6),
                Point2::new(0.6, 0.6),
                Point2::new(0.5, 0.7),
            ];
            assert_eq!(cells.corner_turn(far), f64::INFINITY);
        }
    }
}
