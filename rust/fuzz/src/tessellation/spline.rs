//! Spline bodies for the tessellation target (T-b of REVIEW_NOTES.md):
//! prisms whose profile's right side (and maybe its top) is a clamped
//! quadratic or cubic B-spline, rational or not, with up to two spans,
//! extruded along z from a dyadic origin (its walls ruled spline surfaces, as
//! generate_brep_fixtures.prism builds them), and face bodies on whole
//! nonperiodic spline surfaces bounded by their boundary iso-curves. Both
//! are built as parts, measured and validated; an input whose body does not
//! validate is skipped. The distance to a spline surface is a projection
//! evaluated here in binary64 by Cox-de Boor, independent of the kernel's
//! evaluation and of the tessellation's cells.
use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::curve::{DerivativeOrder, KnotSide};
use rusty_occt::topology::{
    Curve2, Curve3, Edge, EdgeId, Enclosure, Face, FaceId, Fin, FinId, Loop, LoopId, Orientation,
    Region, RegionId, RegionKind, Shell, ShellId, Side, SplineSpan, Surface, Topology,
    TopologyParts, Vertex, VertexId,
};
use rusty_occt::{
    BSplineCurve2, BSplineCurve3, BSplineSurface3, Frame3, KnotVector, Point2, Point3, Tolerance,
    Vec3,
};

fn unit(u: &mut Unstructured) -> Result<f64> {
    Ok(f64::from(u.arbitrary::<u16>()?) / 65535.0)
}

/// A clamped basis of `spans` uniform spans on `[0, 1]`.
fn basis(degree: usize, spans: usize) -> (Vec<f64>, Vec<usize>) {
    let knots = (0..=spans).map(|k| k as f64 / spans as f64).collect();
    let mut mults = vec![1; spans + 1];
    mults[0] = degree + 1;
    mults[spans] = degree + 1;
    (knots, mults)
}

fn weight(u: &mut Unstructured, rational: bool) -> Result<f64> {
    Ok(if rational {
        [0.7, 0.85, 1.0, 1.2, 1.5][u.int_in_range(0..=4usize)?]
    } else {
        1.0
    })
}

/// A profile piece: a line, or a spline's interior poles and weights.
enum Piece {
    Line,
    Spline {
        degree: usize,
        spans: usize,
        interior: Vec<Point2>,
        weights: Vec<f64>,
    },
}

/// A frame: axis-aligned, or tilted about a random axis.
fn frame(u: &mut Unstructured, scale: f64) -> Result<(Point3, Vec3, Vec3, Vec3)> {
    if u.arbitrary::<bool>()? {
        return Ok((
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ));
    }
    let n = Vec3::new(2.0 * unit(u)? - 1.0, 2.0 * unit(u)? - 1.0, 0.3 + unit(u)?);
    let Ok(f) = Frame3::new(
        Point3::new(scale * unit(u)?, -scale * unit(u)?, scale),
        n,
        Vec3::new(1.0, 0.0, 0.0),
        Tolerance::default(),
    ) else {
        return Ok((
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ));
    };
    Ok((f.origin(), f.x(), f.y(), f.normal()))
}

/// The parts with the resolution declared as every gap's enclosure (the
/// check still certifies each gap within it), validated.
fn validated(mut parts: TopologyParts, scale: f64) -> Option<Topology> {
    let linear = 1e-7 * scale;
    let tolerance = Tolerance::new(linear, 1e-12).ok()?;
    let declared = Some(Enclosure::computed(linear));
    parts
        .vertices
        .iter_mut()
        .for_each(|v| v.enclosure = declared);
    parts.fins.iter_mut().for_each(|f| f.enclosure = declared);
    parts.faces.iter_mut().for_each(|f| f.enclosure = declared);
    Topology::from_parts(parts, tolerance).ok()
}

/// Gauss-Legendre nodes and weights of order 8 on `[-1, 1]`.
const GAUSS: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_2, 0.101_228_536_290_376_26),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_47),
    (-0.525_532_409_916_329_0, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329_0, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_47),
    (0.960_289_856_497_536_2, 0.101_228_536_290_376_26),
];

/// `(∫ (x dy - y dx), length)` of a planar spline, by Gauss-Legendre of
/// order 8 on each knot span with the kernel's exact evaluation (exact for
/// the polynomial integrand of a nonrational piece, far within the volume
/// band's slack for the fuzzed rational ones).
fn green(curve: &BSplineCurve2) -> Option<(f64, f64)> {
    let c = curve.as_curve3();
    let (mut twice, mut length) = (0.0, 0.0);
    for w in c.knots().windows(2) {
        let (a, b) = (w[0], w[1]);
        for (x, weight) in GAUSS {
            let t = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let e = c
                .evaluate(t, DerivativeOrder::First, KnotSide::Automatic)
                .ok()?;
            let p = e.position();
            let d = e.derivative_bounds(1)?;
            let (dx, dy) = (d[0].representative(), d[1].representative());
            let k = 0.5 * (b - a) * weight;
            twice += (p.x * dy - p.y * dx) * k;
            length += dx.hypot(dy) * k;
        }
    }
    Some((twice, length))
}

/// A prism of a rectangle whose right side (and maybe its top) is a spline,
/// with its volume and area by quadrature of its profile.
pub(super) fn prism(u: &mut Unstructured) -> Result<Option<(Topology, f64, f64)>> {
    let scale = 2f64.powi(u.int_in_range(-4..=4)?);
    let (w, d) = (
        scale * (1.0 + 3.0 * unit(u)?),
        scale * (1.0 + 3.0 * unit(u)?),
    );
    let h = scale * (0.25 + 2.0 * unit(u)?);
    let rational = u.arbitrary::<bool>()?;
    let spline = |u: &mut Unstructured, from: Point2, to: Point2, out: Vec3| {
        let degree = u.int_in_range(2..=3usize)?;
        let spans = u.int_in_range(1..=2usize)?;
        let n = degree + spans;
        let mut interior = Vec::new();
        let mut weights = vec![1.0];
        for k in 1..n - 1 {
            let f = k as f64 / (n - 1) as f64;
            // Along the chord, pushed out (or a little in) across it.
            let push = (unit(u)? - 0.3) * 0.8;
            interior.push(Point2::new(
                from.x + (to.x - from.x) * f + out.x * push,
                from.y + (to.y - from.y) * f + out.y * push,
            ));
            weights.push(weight(u, rational)?);
        }
        weights.push(1.0);
        Ok(Piece::Spline {
            degree,
            spans,
            interior,
            weights,
        })
    };
    let points = [
        Point2::new(0.0, 0.0),
        Point2::new(w, 0.0),
        Point2::new(w, d),
        Point2::new(0.0, d),
    ];
    let right = spline(u, points[1], points[2], Vec3::new(w, 0.0, 0.0))?;
    let top = if u.arbitrary::<bool>()? {
        spline(u, points[2], points[3], Vec3::new(0.0, d, 0.0))?
    } else {
        Piece::Line
    };
    let pieces = [Piece::Line, right, top, Piece::Line];
    // Axis-aligned, moved by a dyadic offset: a tilted frame's rounded walls
    // cost the validator's spline fluxes seconds under the sanitizer (the
    // sheets are tilted).
    let o = Point3::new(
        scale * f64::from(u.int_in_range(-8i8..=8)?),
        scale * f64::from(u.int_in_range(-8i8..=8)?),
        scale * f64::from(u.int_in_range(-8i8..=8)?),
    );
    let axes = (
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    // The profile's twice signed area and perimeter.
    let (mut twice, mut perimeter) = (0.0, 0.0);
    for (i, piece) in pieces.iter().enumerate() {
        let (p, q) = (points[i], points[(i + 1) % points.len()]);
        match piece {
            Piece::Line => {
                twice += p.x * q.y - q.x * p.y;
                perimeter += (q.x - p.x).hypot(q.y - p.y);
            }
            Piece::Spline {
                degree,
                spans,
                interior,
                weights,
            } => {
                let poles = std::iter::once(p)
                    .chain(interior.iter().copied())
                    .chain(std::iter::once(q))
                    .collect();
                let (knots, mults) = basis(*degree, *spans);
                let curve = BSplineCurve2::new(*degree, poles, Some(weights.clone()), knots, mults)
                    .expect("a clamped basis");
                let Some((t, l)) = green(&curve) else {
                    return Ok(None);
                };
                twice += t;
                perimeter += l;
            }
        }
    }
    let area = 0.5 * twice.abs();
    Ok(validated(
        build_prism(&points, &pieces, h, (o, axes.0, axes.1, axes.2)),
        scale,
    )
    .map(|t| (t, area * h, 2.0 * area + perimeter * h)))
}

fn build_prism(
    points: &[Point2],
    pieces: &[Piece],
    h: f64,
    (o, fx, fy, fz): (Point3, Vec3, Vec3, Vec3),
) -> TopologyParts {
    let n = points.len();
    let map = |p: Point2, z: f64| o + fx * p.x + fy * p.y + fz * z;
    let mut parts = TopologyParts::default();
    for z in [0.0, h] {
        for p in points {
            parts.vertices.push(Vertex {
                position: map(*p, z),
                enclosure: None,
            });
        }
    }
    let line = |a: Point3, b: Point3| Curve3::LineSegment { start: a, end: b };
    let edge = |parts: &mut TopologyParts, a: usize, b: usize, curve: Curve3| {
        parts.edges.push(Edge {
            start: Some(VertexId::new(a)),
            end: Some(VertexId::new(b)),
            curve,
            fins: Vec::new(),
        });
        parts.edges.len() - 1
    };
    let vertical: Vec<usize> = (0..n)
        .map(|i| {
            let (a, b) = (map(points[i], 0.0), map(points[i], h));
            edge(&mut parts, i, n + i, line(a, b))
        })
        .collect();
    // Each piece's bottom and top edges, its curve (for pcurves) and wall.
    let mut bottom_uses = Vec::new();
    let mut top_uses = Vec::new();
    let mut walls = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        let j = (i + 1) % n;
        let (p, q) = (points[i], points[j]);
        match piece {
            Piece::Line => {
                let eb = edge(&mut parts, i, j, line(map(p, 0.0), map(q, 0.0)));
                let et = edge(&mut parts, n + i, n + j, line(map(p, h), map(q, h)));
                bottom_uses.push((
                    eb,
                    Curve2::LineSegment {
                        start: Point2::new(q.x, -q.y),
                        end: Point2::new(p.x, -p.y),
                    },
                ));
                top_uses.push((et, Curve2::LineSegment { start: p, end: q }));
                let length = (q.x - p.x).hypot(q.y - p.y);
                let t = ((q.x - p.x) / length, (q.y - p.y) / length);
                let surface = Surface::Plane(
                    Frame3::new(
                        map(p, 0.0),
                        fx * t.1 - fy * t.0,
                        fx * t.0 + fy * t.1,
                        Tolerance::default(),
                    )
                    .expect("a unit direction"),
                );
                walls.push((surface, eb, et, i, j, [0.0, length]));
            }
            Piece::Spline {
                degree,
                spans,
                interior,
                weights,
            } => {
                let poles: Vec<Point2> = std::iter::once(p)
                    .chain(interior.iter().copied())
                    .chain(std::iter::once(q))
                    .collect();
                let (knots, mults) = basis(*degree, *spans);
                let at = |z: f64| poles.iter().map(|p| map(*p, z)).collect::<Vec<_>>();
                let curve = |z: f64| {
                    BSplineCurve3::new(
                        *degree,
                        at(z),
                        Some(weights.clone()),
                        knots.clone(),
                        mults.clone(),
                    )
                    .expect("a clamped basis")
                };
                let eb = edge(
                    &mut parts,
                    i,
                    j,
                    Curve3::BSpline(SplineSpan::whole(curve(0.0))),
                );
                let et = edge(
                    &mut parts,
                    n + i,
                    n + j,
                    Curve3::BSpline(SplineSpan::whole(curve(h))),
                );
                let plane = |f: &dyn Fn(Point2) -> Point2| {
                    BSplineCurve2::new(
                        *degree,
                        poles.iter().map(|p| f(*p)).collect(),
                        Some(weights.clone()),
                        knots.clone(),
                        mults.clone(),
                    )
                    .expect("a clamped basis")
                };
                bottom_uses.push((
                    eb,
                    Curve2::BSpline(
                        SplineSpan::whole(plane(&|p| Point2::new(p.x, -p.y))).reversed(),
                    ),
                ));
                top_uses.push((et, Curve2::BSpline(SplineSpan::whole(plane(&|p| p)))));
                let rows: Vec<Point3> = poles
                    .iter()
                    .flat_map(|p| [map(*p, 0.0), map(*p, h)])
                    .collect();
                let surface = BSplineSurface3::new(
                    KnotVector::new(*degree, knots.clone(), mults.clone())
                        .expect("a clamped basis"),
                    KnotVector::new(1, vec![0.0, h], vec![2, 2]).expect("a line"),
                    rows,
                    Some(weights.iter().flat_map(|w| [*w, *w]).collect()),
                )
                .expect("a ruled surface");
                walls.push((Surface::BSpline(surface), eb, et, i, j, [0.0, 1.0]));
            }
        }
    }
    let face = |parts: &mut TopologyParts, surface: Surface, uses: Vec<(usize, bool, Curve2)>| {
        let mut fins = Vec::new();
        for (e, forward, pcurve) in uses {
            let id = FinId::new(parts.fins.len());
            parts.fins.push(Fin {
                edge: EdgeId::new(e),
                sense: if forward {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                pcurve,
                enclosure: None,
            });
            parts.edges[e].fins.push(id);
            fins.push(id);
        }
        parts.loops.push(Loop::Edges {
            fins,
            winding: [0, 0],
        });
        parts.faces.push(Face {
            surface,
            sense: Orientation::Forward,
            loops: vec![LoopId::new(parts.loops.len() - 1)],
            front: ShellId::new(0),
            back: ShellId::new(1),
            enclosure: None,
        });
    };
    let bottom = Frame3::new(
        map(Point2::new(0.0, 0.0), 0.0),
        fz * -1.0,
        fx,
        Tolerance::default(),
    )
    .expect("the frame's axes");
    let top = Frame3::new(map(Point2::new(0.0, 0.0), h), fz, fx, Tolerance::default())
        .expect("the frame's axes");
    face(
        &mut parts,
        Surface::Plane(bottom),
        bottom_uses
            .into_iter()
            .rev()
            .map(|(e, c)| (e, false, c))
            .collect(),
    );
    face(
        &mut parts,
        Surface::Plane(top),
        top_uses.into_iter().map(|(e, c)| (e, true, c)).collect(),
    );
    for (surface, eb, et, i, j, [a, b]) in walls {
        let segment = |p: [f64; 2], q: [f64; 2]| Curve2::LineSegment {
            start: Point2::new(p[0], p[1]),
            end: Point2::new(q[0], q[1]),
        };
        face(
            &mut parts,
            surface,
            vec![
                (eb, true, segment([a, 0.0], [b, 0.0])),
                (vertical[j], true, segment([b, 0.0], [b, h])),
                (et, false, segment([b, h], [a, h])),
                (vertical[i], false, segment([a, h], [a, 0.0])),
            ],
        );
    }
    let faces = parts.faces.len();
    parts.shells = vec![
        Shell {
            region: RegionId::new(1),
            sides: (0..faces).map(|f| (FaceId::new(f), Side::Front)).collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        },
        Shell {
            region: RegionId::new(0),
            sides: (0..faces).map(|f| (FaceId::new(f), Side::Back)).collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        },
    ];
    parts.regions = vec![
        Region {
            kind: RegionKind::Void,
            shells: vec![ShellId::new(1)],
        },
        Region {
            kind: RegionKind::Solid,
            shells: vec![ShellId::new(0)],
        },
    ];
    parts
}

/// A face body on a whole spline surface: poles on a jittered grid with
/// random heights, maybe rational; its edges the boundary rows and columns.
pub(super) fn sheet(u: &mut Unstructured) -> Result<Option<Topology>> {
    let scale = 2f64.powi(u.int_in_range(-4..=4)?);
    let (p, q) = (u.int_in_range(2..=3usize)?, u.int_in_range(2..=3usize)?);
    let (su, sv) = (u.int_in_range(1..=2usize)?, u.int_in_range(1..=2usize)?);
    let rational = u.arbitrary::<bool>()?;
    let (nu, nv) = (p + su, q + sv);
    let (w, d, h) = (
        scale * (1.0 + 3.0 * unit(u)?),
        scale * (1.0 + 3.0 * unit(u)?),
        scale * 2.0 * unit(u)?,
    );
    let mut poles = Vec::with_capacity(nu * nv);
    let mut weights = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            let jitter = |u: &mut Unstructured| -> Result<f64> { Ok(0.2 * (unit(u)? - 0.5)) };
            let (fi, fj) = (i as f64 / (nu - 1) as f64, j as f64 / (nv - 1) as f64);
            let (ji, jj) = if 0 < i && i + 1 < nu && 0 < j && j + 1 < nv {
                (jitter(u)?, jitter(u)?)
            } else {
                (0.0, 0.0)
            };
            poles.push(Point3::new(
                w * (fi + ji / (nu - 1) as f64),
                d * (fj + jj / (nv - 1) as f64),
                h * (unit(u)? - 0.5),
            ));
            weights.push(weight(u, rational)?);
        }
    }
    let (uk, um) = basis(p, su);
    let (vk, vm) = basis(q, sv);
    let surface = BSplineSurface3::new(
        KnotVector::new(p, uk.clone(), um.clone()).expect("a clamped basis"),
        KnotVector::new(q, vk.clone(), vm.clone()).expect("a clamped basis"),
        poles.clone(),
        Some(weights.clone()),
    )
    .expect("a surface");
    let (o, x, y, z) = frame(u, scale)?;
    let map = |p: Point3| o + x * p.x + y * p.y + z * p.z;
    let surface = BSplineSurface3::new(
        surface.u_knots().clone(),
        surface.v_knots().clone(),
        poles.iter().map(|p| map(*p)).collect(),
        Some(weights.clone()),
    )
    .expect("a surface");
    let at = |i: usize, j: usize| i * nv + j;
    let row = |j: usize| -> Curve3 {
        Curve3::BSpline(SplineSpan::whole(
            BSplineCurve3::new(
                p,
                (0..nu).map(|i| map(poles[at(i, j)])).collect(),
                Some((0..nu).map(|i| weights[at(i, j)]).collect()),
                uk.clone(),
                um.clone(),
            )
            .expect("a boundary row"),
        ))
    };
    let column = |i: usize| -> Curve3 {
        Curve3::BSpline(SplineSpan::whole(
            BSplineCurve3::new(
                q,
                (0..nv).map(|j| map(poles[at(i, j)])).collect(),
                Some((0..nv).map(|j| weights[at(i, j)]).collect()),
                vk.clone(),
                vm.clone(),
            )
            .expect("a boundary column"),
        ))
    };
    let mut parts = TopologyParts::default();
    for (i, j) in [(0, 0), (nu - 1, 0), (nu - 1, nv - 1), (0, nv - 1)] {
        parts.vertices.push(Vertex {
            position: map(poles[at(i, j)]),
            enclosure: None,
        });
    }
    let edges = [
        (0, 1, row(0)),
        (1, 2, column(nu - 1)),
        (3, 2, row(nv - 1)),
        (0, 3, column(0)),
    ];
    for (a, b, curve) in edges {
        parts.edges.push(Edge {
            start: Some(VertexId::new(a)),
            end: Some(VertexId::new(b)),
            curve,
            fins: vec![FinId::new(parts.edges.len())],
        });
    }
    let corner = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)].map(|(a, b)| Point2::new(a, b));
    for (k, forward) in [(0, true), (1, true), (2, false), (3, false)] {
        parts.fins.push(Fin {
            edge: EdgeId::new(k),
            sense: if forward {
                Orientation::Forward
            } else {
                Orientation::Reversed
            },
            pcurve: Curve2::LineSegment {
                start: corner[k],
                end: corner[(k + 1) % 4],
            },
            enclosure: None,
        });
    }
    parts.loops.push(Loop::Edges {
        fins: (0..4).map(FinId::new).collect(),
        winding: [0, 0],
    });
    parts.faces.push(Face {
        surface: Surface::BSpline(surface),
        sense: Orientation::Forward,
        loops: vec![LoopId::new(0)],
        front: ShellId::new(0),
        back: ShellId::new(0),
        enclosure: None,
    });
    parts.shells.push(Shell {
        region: RegionId::new(0),
        sides: vec![(FaceId::new(0), Side::Front), (FaceId::new(0), Side::Back)],
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    });
    parts.regions.push(Region {
        kind: RegionKind::Void,
        shells: vec![ShellId::new(0)],
    });
    Ok(validated(parts, scale))
}

/// Distances to a spline surface over its domain, independent of the
/// kernel's evaluation and of the tessellation's cells: the surface's basis
/// functions and their derivatives by the Cox-de Boor recurrence in
/// binary64 (The NURBS Book, A2.2 and A2.3), the rational quotient rule;
/// from the four nearest points of a grid with four samples per knot span
/// (at least 9 per direction), Gauss-Newton iterations clamped to the
/// domain. A projection can only overstate a distance (by rounding, far
/// below the check's allowance).
pub(super) struct Projection<'a> {
    surface: &'a BSplineSurface3,
    flat: [Vec<f64>; 2],
    grid: Vec<(Point3, f64, f64)>,
}

/// The flat knot sequence of a nonperiodic basis.
fn flat(k: &KnotVector) -> Vec<f64> {
    k.knots()
        .iter()
        .zip(k.multiplicities())
        .flat_map(|(x, m)| std::iter::repeat_n(*x, *m))
        .collect()
}

/// The span `s` with `flat[s] <= x < flat[s + 1]` (the last nonempty one at
/// the domain's end), and the basis values and first derivatives there.
fn basis_at(flat: &[f64], p: usize, x: f64) -> (usize, Vec<f64>, Vec<f64>) {
    let n = flat.len() - p - 1;
    let mut s = p;
    while s + 1 < n && flat[s + 1] <= x {
        s += 1;
    }
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    ndu[0][0] = 1.0;
    let (mut left, mut right) = (vec![0.0; p + 1], vec![0.0; p + 1]);
    for j in 1..=p {
        left[j] = x - flat[s + 1 - j];
        right[j] = flat[s + j] - x;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let values: Vec<f64> = (0..=p).map(|j| ndu[j][p]).collect();
    // First derivatives: p (N_{r,p-1} / du - N_{r+1,p-1} / du').
    let derivatives = (0..=p)
        .map(|r| {
            let mut d = 0.0;
            if r >= 1 {
                d += ndu[r - 1][p - 1] / ndu[p][r - 1];
            }
            if r < p {
                d -= ndu[r][p - 1] / ndu[p][r];
            }
            p as f64 * d
        })
        .collect();
    (s, values, derivatives)
}

impl<'a> Projection<'a> {
    pub(super) fn new(surface: &'a BSplineSurface3) -> Self {
        let flat = [flat(surface.u_knots()), flat(surface.v_knots())];
        let mut out = Self {
            surface,
            flat,
            grid: Vec::new(),
        };
        let ((u0, u1), (v0, v1)) = surface.domain();
        let count = |k: &KnotVector| (4 * (k.knots().len() - 1)).max(8);
        let (nu, nv) = (count(surface.u_knots()), count(surface.v_knots()));
        for i in 0..=nu {
            for j in 0..=nv {
                let (u, v) = (
                    (u0 + (u1 - u0) * i as f64 / nu as f64).min(u1),
                    (v0 + (v1 - v0) * j as f64 / nv as f64).min(v1),
                );
                let point = out.jet(u, v).0;
                out.grid.push((point, u, v));
            }
        }
        out
    }

    /// `(S, S_u, S_v)` in binary64.
    fn jet(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        let s = self.surface;
        let (p, q) = (s.u_knots().degree(), s.v_knots().degree());
        let nv = s.v_knots().pole_count();
        let (su, nu_, du_) = basis_at(&self.flat[0], p, u);
        let (sv, nv_, dv_) = basis_at(&self.flat[1], q, v);
        let (mut a, mut au, mut av) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut w, mut wu, mut wv) = (0.0, 0.0, 0.0);
        for i in 0..=p {
            for j in 0..=q {
                let index = (su - p + i) * nv + (sv - q + j);
                let (pole, weight) = (s.poles()[index].to_array(), s.weights()[index]);
                let (b, bu, bv) = (nu_[i] * nv_[j], du_[i] * nv_[j], nu_[i] * dv_[j]);
                for k in 0..3 {
                    a[k] += b * weight * pole[k];
                    au[k] += bu * weight * pole[k];
                    av[k] += bv * weight * pole[k];
                }
                w += b * weight;
                wu += bu * weight;
                wv += bv * weight;
            }
        }
        let point = [a[0] / w, a[1] / w, a[2] / w];
        let d = |x: [f64; 3], dw: f64| {
            Vec3::new(
                (x[0] - dw * point[0]) / w,
                (x[1] - dw * point[1]) / w,
                (x[2] - dw * point[2]) / w,
            )
        };
        (
            Point3::new(point[0], point[1], point[2]),
            d(au, wu),
            d(av, wv),
        )
    }

    pub(super) fn distance(&self, p: Point3) -> f64 {
        let mut near: Vec<(f64, f64, f64)> = self
            .grid
            .iter()
            .map(|(q, u, v)| (q.distance(p), *u, *v))
            .collect();
        near.sort_by(|a, b| a.0.total_cmp(&b.0));
        near.iter()
            .take(4)
            .map(|&(d, u, v)| d.min(self.descend(p, u, v)))
            .fold(f64::INFINITY, f64::min)
    }

    fn descend(&self, p: Point3, mut u: f64, mut v: f64) -> f64 {
        let ((u0, u1), (v0, v1)) = self.surface.domain();
        let mut best = f64::INFINITY;
        for _ in 0..30 {
            let (point, su, sv) = self.jet(u, v);
            let r = point - p;
            best = best.min(r.length());
            let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
            let (g, h) = (-r.dot(su), -r.dot(sv));
            let det = a * c - b * b;
            if det.is_nan() || det <= 0.0 {
                break;
            }
            let (nu, nv) = (
                (u + (g * c - h * b) / det).clamp(u0, u1),
                (v + (a * h - b * g) / det).clamp(v0, v1),
            );
            if (nu - u).abs() <= 1e-15 * (u1 - u0) && (nv - v).abs() <= 1e-15 * (v1 - v0) {
                break;
            }
            (u, v) = (nu, nv);
        }
        best
    }
}
