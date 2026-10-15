//! S9e.4b.4d.2a: an imported turned body of cylinders, cones and planes
//! about one axis meeting at corners (REVIEW_NOTES.md, "S9e.4b.4d.2
//! refined").
//!
//! Its meridian profile is read off the stored faces: every curved face a
//! cylinder or a cone about one axis line within the resolution, every plane
//! normal to it, every edge a whole circle about it, the faces a path from
//! the axis to the axis (an apex or a disc's centre at each end) whose
//! heights never fall. The axis frame is the first curved face's stored
//! frame with its axes' components that move the body by less than a
//! thousandth of the resolution taken to zero (a stored axis within rounding
//! of the world's). Each curved face's generatrix is a line `r = r0 + k w`
//! in exact rationals of its stored data on the axis frame (`w` its height,
//! a cone's slope its half-angle's tangent rounded once), each plane its
//! height rounded once; each corner of the profile the two lines' (or a
//! line's and a plane's) common point, which must lie within the resolution
//! of its stored circle. The body is the lattice of its primitives
//! (Ovchinnikov's max-min form of a piecewise linear function): for each
//! curved face `i` the common of the primitives whose lines lie on or
//! outside `i`'s over its heights and whose ranges may cover them, each such
//! term then made as small as the form allows (a member dropped where the
//! form is still the profile), the fuse of those commons, terms holding
//! another's dropped and the common of the primitives every term holds taken
//! out. Each primitive's range: at a cap
//! (a step's or an end's plane on its wider side) that plane's height, at
//! a step's narrower side past the plane into the wider side by a quarter
//! of that one's height, at its apex its apex, and past a corner the
//! heights of the faces whose terms it bounds and a quarter of its own
//! height beyond. The form is checked exactly on the lines (between every
//! two breakpoints, at two points: the profile's radius, empty outside it);
//! the engine builds its tree as any chain's, every primitive on the axis
//! frame (a cone's exact model on the axis at its height, `build_cone_on`),
//! and its model's match to the stored topology is the arbiter.
//!
//! S9e.4b.4d.2b: a cone within the resolution of a plane (a step or an end
//! a few millionths from flat, `flat_kind`) is taken as the plane normal to
//! the axis through its wider rim (its stored circle of the larger radius),
//! its other rim the plane's meeting with its neighbour's line, a vertex
//! loop on it (its apex) the disc's centre within the resolution; the body
//! is then stored so (the face on that plane, its rims the model's circles,
//! their fins' pcurves again), its import, history and Booleans those of a
//! body of planes.
//!
//! The fix of S9e.4b.4d.2c.1a's consistency: every turned body is stored as
//! its model (`modelled`: every face on its model's surface on the axis
//! frame, every edge its model's circle at its node, every pcurve again), not
//! its flats alone. A Boolean's results keep the stored surfaces of the
//! faces and the curves of the edges they keep beside the model's new edges,
//! so a stored surface or edge off the model within the resolution (an axis
//! tilted within it, a corner's ring off its lines' meeting, a torus off its
//! snap) left faces bounded by edges off their surfaces, whose certified
//! measures were no set's: `OCC485a`'s Booleans with its rib disagreed with
//! one another beyond their enclosures, and boxes across it failed the
//! validator.
//!
//! S9e.4b.4d.2c.1a: a torus face between a plane and a cylinder is a fillet
//! (its material outside the tube) or a round, snapped tangent to both: the
//! plane's height, the cylinder's radius and the tube's radius on the
//! multiples of `2^-40` of the body's size, its centre where the two lines
//! offset by the tube toward its stored centre meet (within the resolution
//! of that), its ends the centre's feet on the lines (each within the
//! resolution of its stored circle), all of them binary64 numbers exactly.
//! The body is the sharp profile's lattice (each arc replaced by the corner
//! its lines meet at, `Meridian::sharp`) fused with each fillet's piece and
//! less each round's, a piece the box from the tube's centre past the corner
//! (an annular prism on the axis frame, `corner_boxes`) less the whole torus
//! on the axis frame (`build_torus_on`): the piece's own Boolean meets the
//! tube transversally, and the lattice's with it has the tube touch the
//! lattice's plane or cylinder along the piece's edge at the tangent point
//! (`graph.rs`'s touching circles, its tree's `arrange_touching`).
//!
//! S9e.4b.4d.2c.1b: a torus face between a cone and a plane or a cylinder
//! (`cone_corner`): the other line's number and the tube's radius on that
//! grid, the centre the tube from that line, the direction from the stored
//! centre toward the stored cone turned to the rational one whose half-angle
//! tangent is the nearest multiple of `2^-40`, the cone the line through the
//! tube's point in that direction normal to it: no binary64 data, so its
//! primitive carries its exact radii (`Solid::build_cone_exact_on`). The
//! piece's box beside a cone is the parallelogram in the meridian along the
//! two directions from the centre to its tangent points (`slant_box`): the
//! frusta of cones normal to the cone through the centre and past the
//! corner, between the plane through the centre and a far plane beside a
//! cylinder, common an annulus between the cylinder through the centre and
//! a far one beside a plane. A body of one curved face is its one primitive
//! with its pieces. A cone beside two tori, or a torus between two cones, is
//! no turned body's (S9e.4b.4d.2c.1b refined).
use super::*;
use crate::topology::{Curve2, LoopId, Projection};

/// The choice a turned body's piece records (`Imported::moved` reads it off
/// the moved stored topology again by it).
pub(crate) const TURNED: u32 = u32::MAX;

fn int(n: i64) -> R {
    R::from_integer(n.into())
}

/// A meridian line `r = r0 + k w` on the axis frame's heights.
#[derive(Debug, Clone, PartialEq)]
struct Line {
    r0: R,
    k: R,
}

impl Line {
    fn at(&self, w: &R) -> R {
        &self.r0 + &self.k * w
    }

    /// Its apex's height (a cone's).
    fn apex(&self) -> Option<R> {
        (self.k != zero()).then(|| -(&self.r0 / &self.k))
    }
}

/// A stored face's meridian: a curved face's line, a plane's height;
/// S9e.4b.4d.2c.1a: a torus face's tube.
#[derive(Debug, Clone)]
enum Kind {
    Curved(Line),
    Plane(R),
    Arc(Tube),
}

/// S9e.4b.4d.2c.1a: a torus face's tube in the meridian: its stored centre
/// `(r, w)` and radius, and whether its material lies inside the tube (a
/// round) or outside it (a fillet).
#[derive(Debug, Clone)]
struct Tube {
    centre: [R; 2],
    radius: R,
    round: bool,
}

/// S9e.4b.4d.2c.1a: a fillet or a round snapped tangent to its neighbouring
/// lines (a plane's and a cylinder's): its path position, the tube's centre
/// and radius, the sharp profile's corner where the lines meet, whether it
/// is a round.
#[derive(Debug, Clone)]
struct Corner {
    at: usize,
    centre: [R; 2],
    radius: R,
    corner: [R; 2],
    round: bool,
    /// S9e.4b.4d.2c.1b: beside a cone, the unit directions from the centre
    /// to its tangent points on the cone and on the other line (its box's
    /// axes).
    slant: Option<[[R; 2]; 2]>,
}

/// One end of a face's segment in the meridian.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum End {
    Ring(usize),
    Apex(usize),
    Centre,
}

/// The profile read off the stored faces: the axis frame, each face's
/// meridian, and the path from the bottom to the top (each face between
/// two nodes `(r, w)`).
#[derive(Debug, Clone)]
struct Meridian {
    axis: Frame3,
    kinds: Vec<Kind>,
    path: Vec<usize>,
    nodes: Vec<[R; 2]>,
    /// S9e.4b.4d.2b: each edge's node.
    rings: Vec<[R; 2]>,
    /// S9e.4b.4d.2c.1a: each torus face's snapped tube.
    corners: Vec<Corner>,
    /// S9e.4b.4d.2c.1b: each face a cone snapped tangent to a torus (its
    /// line no binary64 data: its primitive's exact radii).
    exact: Vec<bool>,
}

impl Meridian {
    /// S9e.4b.4d.2c.1a: the sharp profile, each fillet or round replaced by
    /// the corner its lines meet at: its path and nodes.
    fn sharp(&self) -> (Vec<usize>, Vec<[R; 2]>) {
        let mut path = Vec::new();
        let mut nodes = vec![self.nodes[0].clone()];
        for (s, &i) in self.path.iter().enumerate() {
            if let Kind::Arc(_) = self.kinds[i] {
                let c = self.corners.iter().find(|c| c.at == s).expect("a corner");
                *nodes.last_mut().expect("a node") = c.corner.clone();
                continue;
            }
            path.push(i);
            nodes.push(self.nodes[s + 1].clone());
        }
        (path, nodes)
    }
}

/// S9e.4b.4d.2c.1a: `x` on the multiples of `g` (a power of two).
fn on_grid(x: &R, g: f64) -> R {
    q((rational_f64(x) / g).round() * g)
}

/// Whether two directions agree within `tol` over `size`.
fn along(a: Vec3, b: Vec3, size: f64, tol: f64) -> bool {
    a.cross(b).length() * size <= tol * a.length() * b.length()
}

/// The meridian of a turned body of cylinders, cones and planes, or `None`.
#[allow(clippy::too_many_lines)]
fn meridian(t: &Topology, tol: f64) -> Option<Meridian> {
    if !one_shell(t) {
        return None;
    }
    let faces = t.faces();
    let first = faces
        .iter()
        .position(|f| matches!(f.surface, Surface::Cylinder { .. } | Surface::Cone { .. }))?;
    let bounds = super::super::split::edge_bounds(t);
    let stored = axis_of(t, first)?;
    let size = (bounds.max - bounds.min).length() + (stored.origin() - bounds.min).length();
    let axis = snapped(stored, size, tol)?;
    let local = Local::new(&axis).ok()?;
    let n = axis.normal();
    let on_axis_line = |frame: &Frame3| {
        let c = axis.coordinates(frame.origin());
        along(n, frame.normal(), size, tol) && c[0].hypot(c[1]) <= tol
    };
    let mut kinds = Vec::new();
    // S9e.4b.4d.2b: each cone within the resolution of a plane, taken as the
    // plane through its wider rim (below); one within the resolution of a
    // cylinder is no turned body's.
    let flats: Vec<bool> = (0..faces.len())
        .map(|i| flat_kind(t, i, tol) == Some(Flat::Plane))
        .collect();
    for (i, f) in faces.iter().enumerate() {
        if let (true, Surface::Cone { frame, .. }) = (flats[i], &f.surface) {
            if !on_axis_line(frame) {
                return None;
            }
            kinds.push(Kind::Plane(zero()));
            continue;
        }
        kinds.push(match &f.surface {
            Surface::Plane(p) => {
                if !along(p.normal(), n, size, tol) {
                    return None;
                }
                Kind::Plane(q(rational_f64(&local.of(p.origin())[2])))
            }
            Surface::Cylinder { frame, radius } => {
                if !on_axis_line(frame) || material_outside(t, i) {
                    return None;
                }
                Kind::Curved(Line {
                    r0: q(*radius),
                    k: zero(),
                })
            }
            // S9e.4b.4d.2c.1a: a torus's tube in the meridian.
            Surface::Torus {
                frame,
                major,
                minor,
            } => {
                if !on_axis_line(frame) {
                    return None;
                }
                Kind::Arc(Tube {
                    centre: [q(*major), local.of(frame.origin())[2].clone()],
                    radius: q(*minor),
                    round: !material_outside(t, i),
                })
            }
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => {
                if !on_axis_line(frame) || material_outside(t, i) {
                    return None;
                }
                let sense = if n.dot(frame.normal()) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                let k = q(sense * half_angle.tan());
                let c = local.of(frame.origin())[2].clone();
                Kind::Curved(Line {
                    r0: q(*radius) - &k * &c,
                    k,
                })
            }
            _ => return None,
        });
    }
    // Every edge a whole circle about the axis between two faces.
    let owner = fin_faces(t);
    let mut rings = Vec::new();
    for e in t.edges() {
        let (frame, radius) = match &e.curve {
            Curve3::Circle { frame, radius } if e.is_ring() => (frame, radius),
            Curve3::CircularArc {
                frame,
                radius,
                sweep_angle,
                ..
            } if sweep_angle.abs() == std::f64::consts::TAU && e.start == e.end => (frame, radius),
            _ => return None,
        };
        if !on_axis_line(frame) {
            return None;
        }
        let [a, b] = e.fins.as_slice() else {
            return None;
        };
        let (fa, fb) = (owner[a.index()], owner[b.index()]);
        if fa == usize::MAX || fb == usize::MAX || fa == fb {
            return None;
        }
        rings.push((fa, fb, *radius, local.of(frame.origin())[2].clone()));
    }
    // S9e.4b.4d.2b: a cone within the resolution of a plane taken as the
    // plane through its wider rim (its stored circle of the larger radius,
    // that circle's height rounded once).
    for (i, &flat) in flats.iter().enumerate() {
        if !flat {
            continue;
        }
        let wider = rings
            .iter()
            .filter(|(a, b, _, _)| *a == i || *b == i)
            .max_by(|x, y| x.2.total_cmp(&y.2))?;
        kinds[i] = Kind::Plane(q(rational_f64(&wider.3)));
    }
    // Each face's two ends: its rings, a cone's apex, a disc's centre
    // (S9e.4b.4d.2b: a vertex loop of a plane, a cone's apex where the cone
    // is taken as the plane, the stored vertex within the resolution of the
    // centre).
    let mut ends: Vec<Vec<End>> = vec![Vec::new(); faces.len()];
    for (e, (a, b, _, _)) in rings.iter().enumerate() {
        ends[*a].push(End::Ring(e));
        ends[*b].push(End::Ring(e));
    }
    for (i, f) in faces.iter().enumerate() {
        for l in &f.loops {
            if let Loop::Vertex(v) = &t.loops()[l.index()] {
                match &kinds[i] {
                    Kind::Arc(_) => return None,
                    Kind::Curved(_) => ends[i].push(End::Apex(v.index())),
                    Kind::Plane(h) => {
                        let c = local.of(t.vertices()[v.index()].position);
                        let off = rational_f64(&c[0]).hypot(rational_f64(&c[1]));
                        if off > tol || rational_f64(&(&c[2] - h)).abs() > tol {
                            return None;
                        }
                    }
                }
            }
        }
        if matches!(kinds[i], Kind::Plane(_)) && ends[i].len() == 1 {
            ends[i].push(End::Centre);
        }
        if ends[i].len() != 2 || ends[i][0] == ends[i][1] {
            return None;
        }
    }
    // The path from a face with an end on the axis to the other.
    let on_axis = |x: &End| !matches!(x, End::Ring(_));
    let start = (0..faces.len()).find(|&i| ends[i].iter().any(on_axis))?;
    let mut path = vec![start];
    let mut from = *ends[start].iter().find(|x| on_axis(x))?;
    let mut at = start;
    let mut path_ends = vec![from];
    loop {
        let next = *ends[at].iter().find(|&&x| x != from)?;
        path_ends.push(next);
        let End::Ring(e) = next else { break };
        let (a, b, _, _) = &rings[e];
        let face = if *a == at { *b } else { *a };
        if path.contains(&face) {
            return None;
        }
        path.push(face);
        from = next;
        at = face;
    }
    if path.len() != faces.len() {
        return None;
    }
    // S9e.4b.4d.2c.1a: each torus face a fillet or a round snapped tangent
    // to its two lines; between a plane and a cylinder, their height and
    // radius and its tube's radius on the multiples of `2^-40` of the body's
    // size (a power of two), so that its centre, its tangent points and its
    // piece's box are binary64 numbers exactly; its centre where the two
    // lines offset by the tube toward its stored centre meet, within the
    // resolution of that. S9e.4b.4d.2c.1b: beside a cone, `cone_corner`.
    let grid = 2f64.powi(size.log2().ceil() as i32 - 40);
    let mut corners = Vec::new();
    let mut exact = vec![false; faces.len()];
    for (s, &i) in path.iter().enumerate() {
        let Kind::Arc(tube) = kinds[i].clone() else {
            continue;
        };
        let (a, b) = (path[s.checked_sub(1)?], *path.get(s + 1)?);
        // S9e.4b.4d.2c.1b: beside a cone.
        let cone = match (&kinds[a], &kinds[b]) {
            (Kind::Curved(l), _) if l.k != zero() => Some((a, b)),
            (_, Kind::Curved(l)) if l.k != zero() => Some((b, a)),
            _ => None,
        };
        if let Some(pair) = cone {
            let c = cone_corner(&mut kinds, &mut exact, pair, &tube, grid, tol)?;
            corners.push(Corner { at: s, ..c });
            continue;
        }
        let (cyl, plane) = match (&kinds[a], &kinds[b]) {
            (Kind::Curved(l), Kind::Plane(_)) if l.k == zero() => (a, b),
            (Kind::Plane(_), Kind::Curved(l)) if l.k == zero() => (b, a),
            _ => return None,
        };
        let Kind::Curved(line) = &kinds[cyl] else {
            return None;
        };
        let radius = on_grid(&line.r0, grid);
        let Kind::Plane(h) = &kinds[plane] else {
            return None;
        };
        let h = on_grid(h, grid);
        kinds[cyl] = Kind::Curved(Line {
            r0: radius.clone(),
            k: zero(),
        });
        kinds[plane] = Kind::Plane(h.clone());
        let t = on_grid(&tube.radius, grid);
        let sign = |x: R| match x.cmp(&zero()) {
            std::cmp::Ordering::Greater => Some(int(1)),
            std::cmp::Ordering::Less => Some(int(-1)),
            std::cmp::Ordering::Equal => None,
        };
        let centre = [
            &radius + sign(&tube.centre[0] - &radius)? * &t,
            &h + sign(&tube.centre[1] - &h)? * &t,
        ];
        let moved = rational_f64(&(&centre[0] - &tube.centre[0]))
            .hypot(rational_f64(&(&centre[1] - &tube.centre[1])));
        if t <= zero() || moved > tol || centre[0] <= zero() {
            return None;
        }
        corners.push(Corner {
            at: s,
            centre,
            radius: t,
            corner: [radius, h],
            round: tube.round,
            slant: None,
        });
    }
    // Each node exactly: a corner the two meridians' common point, within
    // the resolution of its stored circle.
    let mut nodes: Vec<[R; 2]> = Vec::new();
    let mut at_ring: Vec<[R; 2]> = vec![[zero(), zero()]; rings.len()];
    for (s, end) in path_ends.iter().enumerate() {
        let own = &kinds[path[s.min(path.len() - 1)]];
        let node = match end {
            End::Ring(e) => {
                // S9e.4b.4d.2c.1a: a fillet's or a round's end its tangent
                // point on its neighbour's line (the centre's foot there).
                // S9e.4b.4d.2c.1b: on a cone, the tube's point along its
                // direction to it.
                let tangent = |at: usize, line: &Kind| -> Option<[R; 2]> {
                    let c = corners.iter().find(|c| c.at == at)?;
                    match line {
                        Kind::Plane(h) => Some([c.centre[0].clone(), h.clone()]),
                        Kind::Curved(l) if l.k == zero() => {
                            Some([l.r0.clone(), c.centre[1].clone()])
                        }
                        Kind::Curved(_) => {
                            let [e, _] = c.slant.as_ref()?;
                            Some([
                                &c.centre[0] + &c.radius * &e[0],
                                &c.centre[1] + &c.radius * &e[1],
                            ])
                        }
                        Kind::Arc(_) => None,
                    }
                };
                let node = match (&kinds[path[s - 1]], own) {
                    (Kind::Arc(_), line) => tangent(s - 1, line)?,
                    (line, Kind::Arc(_)) => tangent(s, line)?,
                    (Kind::Curved(a), Kind::Curved(b)) => {
                        if a.k == b.k {
                            return None;
                        }
                        let w = (&b.r0 - &a.r0) / (&a.k - &b.k);
                        [a.at(&w), w]
                    }
                    (Kind::Curved(a), Kind::Plane(h)) | (Kind::Plane(h), Kind::Curved(a)) => {
                        [a.at(h), h.clone()]
                    }
                    _ => return None,
                };
                let (_, _, radius, w) = &rings[*e];
                if (rational_f64(&node[0]) - radius).abs() > tol
                    || rational_f64(&(&node[1] - w)).abs() > tol
                {
                    return None;
                }
                at_ring[*e] = node.clone();
                node
            }
            End::Apex(v) => {
                let Kind::Curved(l) = own else { return None };
                let w = l.apex()?;
                let c = local.of(t.vertices()[*v].position);
                if rational_f64(&(&c[2] - &w)).abs() > tol {
                    return None;
                }
                [zero(), w]
            }
            End::Centre => {
                let Kind::Plane(h) = own else { return None };
                [zero(), h.clone()]
            }
        };
        nodes.push(node);
    }
    // From the bottom up, the heights never falling, a curved face rising
    // by more than the resolution, every corner off the axis.
    let reversed = nodes[0][1] > nodes[nodes.len() - 1][1];
    if reversed {
        nodes.reverse();
        path.reverse();
    }
    for (s, &i) in path.iter().enumerate() {
        let (a, b) = (&nodes[s][1], &nodes[s + 1][1]);
        match kinds[i] {
            Kind::Curved(_) | Kind::Arc(_) if rational_f64(&(b - a)) > tol => {}
            Kind::Plane(_) if a == b => {}
            _ => return None,
        }
    }
    if nodes[1..nodes.len() - 1]
        .iter()
        .any(|x| rational_f64(&x[0]) <= tol)
    {
        return None;
    }
    // Two curved faces on one line (a face split on its surface) are no
    // corner.
    let lines: Vec<&Line> = path
        .iter()
        .filter_map(|&i| match &kinds[i] {
            Kind::Curved(l) => Some(l),
            Kind::Plane(_) | Kind::Arc(_) => None,
        })
        .collect();
    for (a, x) in lines.iter().enumerate() {
        if lines[a + 1..].contains(x) {
            return None;
        }
    }
    // A path read from the top down: its corners' places from the bottom.
    if reversed {
        let n = path.len();
        for c in &mut corners {
            c.at = n - 1 - c.at;
        }
    }
    Some(Meridian {
        axis,
        kinds,
        path,
        nodes,
        rings: at_ring,
        corners,
        exact,
    })
}

/// S9e.4b.4d.2c.1b: a fillet or a round between a cone (`cone`, its face)
/// and a plane or a cylinder (`other`) snapped tangent to both: the other's
/// height or radius and the tube's radius on the multiples of `grid`, the
/// centre the tube from it (its other number its stored one's on the grid),
/// the direction from the stored centre toward the stored cone turned to the
/// rational one whose half-angle tangent lies on the multiples of `2^-40`,
/// the cone the line through the tube's point in that direction normal to
/// it (exact, no binary64 data); the centre within the resolution of its
/// stored one. A cone beside two tori is no turned body's.
fn cone_corner(
    kinds: &mut [Kind],
    exact: &mut [bool],
    (cone, other): (usize, usize),
    tube: &Tube,
    grid: f64,
    tol: f64,
) -> Option<Corner> {
    let Kind::Curved(stored) = kinds[cone].clone() else {
        return None;
    };
    if exact[cone] {
        return None;
    }
    let t = on_grid(&tube.radius, grid);
    let sign = |x: R| match x.cmp(&zero()) {
        std::cmp::Ordering::Greater => Some(int(1)),
        std::cmp::Ordering::Less => Some(int(-1)),
        std::cmp::Ordering::Equal => None,
    };
    // The other line snapped, the centre the tube from it, the direction to
    // it from the centre.
    let (centre, toward) = match kinds[other].clone() {
        Kind::Plane(h) => {
            let h = on_grid(&h, grid);
            kinds[other] = Kind::Plane(h.clone());
            let up = sign(&tube.centre[1] - &h)?;
            let centre = [on_grid(&tube.centre[0], grid), &h + &up * &t];
            (centre, [zero(), -up])
        }
        Kind::Curved(l) if l.k == zero() => {
            let radius = on_grid(&l.r0, grid);
            kinds[other] = Kind::Curved(Line {
                r0: radius.clone(),
                k: zero(),
            });
            let out = sign(&tube.centre[0] - &radius)?;
            let centre = [&radius + &out * &t, on_grid(&tube.centre[1], grid)];
            (centre, [-out, zero()])
        }
        _ => return None,
    };
    // The direction from the stored centre to the stored cone (binary64),
    // its half-angle tangent on the multiples of `2^-40`.
    let (k, r0) = (rational_f64(&stored.k), rational_f64(&stored.r0));
    let (cr, cw) = (rational_f64(&tube.centre[0]), rational_f64(&tube.centre[1]));
    let len = (1.0 + k * k).sqrt();
    let side = if cr - k * cw - r0 > 0.0 { -1.0 } else { 1.0 };
    let n = [side / len, -side * k / len];
    let half = if n[0] >= 0.0 {
        n[1] / (1.0 + n[0])
    } else {
        (1.0 - n[0]) / n[1]
    };
    let unit = 2f64.powi(-40);
    if !half.is_finite() || half.abs() > 4096.0 {
        return None;
    }
    let half = q((half / unit).round() * unit);
    let den = int(1) + &half * &half;
    let e = [(int(1) - &half * &half) / &den, (&half * int(2)) / &den];
    if e[0] == zero() || e[1] == zero() {
        return None;
    }
    let point = [&centre[0] + &t * &e[0], &centre[1] + &t * &e[1]];
    let slope = -(&e[1] / &e[0]);
    let line = Line {
        r0: &point[0] - &slope * &point[1],
        k: slope,
    };
    let moved = rational_f64(&(&centre[0] - &tube.centre[0]))
        .hypot(rational_f64(&(&centre[1] - &tube.centre[1])));
    if t <= zero() || moved > tol || centre[0] <= zero() {
        return None;
    }
    // The sharp profile's corner: the cone's and the other line's.
    let corner = match &kinds[other] {
        Kind::Plane(h) => [line.at(h), h.clone()],
        Kind::Curved(l) => {
            let w = (&l.r0 - &line.r0) / &line.k;
            [l.r0.clone(), w]
        }
        Kind::Arc(_) => return None,
    };
    kinds[cone] = Kind::Curved(line);
    exact[cone] = true;
    Some(Corner {
        at: 0,
        centre,
        radius: t,
        corner,
        round: tube.round,
        slant: Some([e, toward]),
    })
}

/// A turned body's stored topology as its model (the fix of S9e.4b.4d.2c.1a's
/// consistency, generalizing S9e.4b.4d.2b's flats): every face on its model's
/// surface on the axis frame, every edge its model's circle at its node, so
/// that a Boolean's results, which keep the stored surfaces and the curves of
/// the stored edges they keep, bound the model's sets. A plane (S9e.4b.4d.2b:
/// a cone within the resolution of a plane among them) the plane normal to
/// the axis at its height, its sense keeping the side the face leaves its
/// material on, its loops' windings zero and its outer loop first; a
/// cylinder its line's radius, a cone its line through its stored frame's
/// height rounded once, a torus (S9e.4b.4d.2c.1a) its snapped tube; each on
/// the axis frame (turned over where the stored frame's normal points down
/// the axis), its stored frame's `x` turned about the axis to the frame's.
/// Each edge the circle of its node's radius at its node's height on the
/// axis frame (rounded once), turned over as its stored circle; a vertex on
/// it at its curve's start; every fin's pcurve given again (a plane's arc of
/// the circle, a curved face's line at the circle's parameters, its stored
/// span and sheet kept); a vertex loop (a cone's apex, a disc's centre) at
/// its model's point on the axis. `None` where it is no turned body or a fin
/// is no such pcurve's.
#[allow(clippy::too_many_lines)]
pub(super) fn modelled(t: &Topology, tolerance: Tolerance) -> Option<Topology> {
    let m = meridian(t, tolerance.linear())?;
    let local = Local::new(&m.axis).ok()?;
    let n = m.axis.normal();
    let owner = fin_faces(t);
    let mut parts = t.to_parts();
    // A frame on the axis at a local height, turned over where `down`.
    let on_axis = |w: &R, down: bool| {
        let f = m.axis.at(local.axis_point(w));
        if down {
            f.flipped()
        } else {
            f
        }
    };
    let down = |f: &Frame3| f.normal().dot(n) < 0.0;
    // Each face's model surface.
    for (i, f) in t.faces().iter().enumerate() {
        parts.faces[i].surface = match (&m.kinds[i], &f.surface) {
            (Kind::Plane(h), _) => {
                let fin = f.loops.iter().find_map(|l| match &t.loops()[l.index()] {
                    Loop::Edges { fins, .. } => fins.first().copied(),
                    Loop::Vertex(_) => None,
                })?;
                let at = t.edges()[t.fins()[fin.index()].edge.index()]
                    .curve
                    .point(0.5);
                let out = oriented_normal(t, i, at)?;
                let frame = on_axis(h, false);
                parts.faces[i].sense = if frame.normal().dot(out) > 0.0 {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                // Its loops wind about no period on the plane, its outer loop
                // (its wider rim) first.
                let radius = |l: &LoopId| match &t.loops()[l.index()] {
                    Loop::Edges { fins, .. } => fins.first().map_or(0.0, |x| {
                        rational_f64(&m.rings[t.fins()[x.index()].edge.index()][0])
                    }),
                    Loop::Vertex(_) => 0.0,
                };
                parts.faces[i]
                    .loops
                    .sort_by(|a, b| radius(b).total_cmp(&radius(a)));
                for l in &f.loops {
                    if let Loop::Edges { winding, .. } = &mut parts.loops[l.index()] {
                        *winding = [0, 0];
                    }
                }
                Surface::Plane(frame)
            }
            (Kind::Curved(line), Surface::Cylinder { frame, .. }) if line.k == zero() => {
                Surface::Cylinder {
                    frame: on_axis(&zero(), down(frame)),
                    radius: rational_f64(&line.r0),
                }
            }
            (Kind::Curved(line), Surface::Cone { frame, .. }) if line.k != zero() => {
                let w = q(rational_f64(&local.of(frame.origin())[2]));
                let radius = rational_f64(&line.at(&w));
                if radius < 0.0 {
                    return None;
                }
                let sense = if down(frame) { -1.0 } else { 1.0 };
                Surface::Cone {
                    frame: on_axis(&w, down(frame)),
                    radius,
                    half_angle: (sense * rational_f64(&line.k)).atan(),
                }
            }
            (Kind::Arc(_), Surface::Torus { frame, .. }) => {
                let c = m.corners.iter().find(|c| m.path[c.at] == i)?;
                Surface::Torus {
                    frame: on_axis(&c.centre[1], down(frame)),
                    major: rational_f64(&c.centre[0]),
                    minor: rational_f64(&c.radius),
                }
            }
            _ => return None,
        };
    }
    // Each edge its model's circle, a vertex on it at its start.
    for (e, edge) in t.edges().iter().enumerate() {
        let [r, w] = &m.rings[e];
        let radius = rational_f64(r);
        parts.edges[e].curve = match &edge.curve {
            Curve3::Circle { frame, .. } => Curve3::Circle {
                frame: on_axis(w, down(frame)),
                radius,
            },
            Curve3::CircularArc {
                frame, sweep_angle, ..
            } => {
                let at = on_axis(w, down(frame));
                let c = at.coordinates(edge.curve.point(0.0));
                Curve3::CircularArc {
                    frame: at,
                    radius,
                    start_angle: c[1].atan2(c[0]),
                    sweep_angle: *sweep_angle,
                }
            }
            _ => return None,
        };
        let start = parts.edges[e].curve.point(0.0);
        for v in [edge.start, edge.end].into_iter().flatten() {
            parts.vertices[v.index()].position = start;
        }
    }
    // Each vertex loop at its model's point on the axis: a cone's apex, a
    // disc's centre.
    for (i, f) in t.faces().iter().enumerate() {
        for l in &f.loops {
            if let Loop::Vertex(v) = &t.loops()[l.index()] {
                let w = match &m.kinds[i] {
                    Kind::Curved(line) => line.apex()?,
                    Kind::Plane(h) => h.clone(),
                    Kind::Arc(_) => return None,
                };
                parts.vertices[v.index()].position = local.axis_point(&w);
            }
        }
    }
    // Every fin's pcurve again: a plane's arc of its circle, a curved face's
    // line at the circle's parameters (its stored span and sheet kept).
    for (x, fin) in t.fins().iter().enumerate() {
        let (fi, e) = (owner[x], fin.edge.index());
        if fi == usize::MAX {
            return None;
        }
        let curve = &parts.edges[e].curve;
        let surface = &parts.faces[fi].surface;
        parts.fins[x].pcurve = match (surface, &fin.pcurve) {
            (Surface::Plane(p), _) => plane_arc(p, curve, fin.sense)?,
            (
                Surface::Cylinder { .. } | Surface::Cone { .. } | Surface::Torus { .. },
                Curve2::LineSegment { start, end },
            ) if start.y == end.y => {
                let near = Projection::inverse(surface, curve.point(0.0))?;
                let tau = std::f64::consts::TAU;
                // On the stored pcurve's sheets.
                let u = near.x + ((start.x - near.x) / tau).round() * tau;
                let v = match surface {
                    Surface::Torus { .. } => near.y + ((start.y - near.y) / tau).round() * tau,
                    _ => near.y,
                };
                Curve2::LineSegment {
                    start: Point2::new(u, v),
                    end: Point2::new(u + (end.x - start.x), v),
                }
            }
            _ => return None,
        };
    }
    let model = Topology::from_parts(parts.with_measured_enclosures(), tolerance).ok()?;
    Some(model.with_identity_of(t))
}

/// A circle's pcurve on a plane: the arc of its image in the plane's
/// coordinates, at the circle's parameters the way the fin runs.
fn plane_arc(p: &Frame3, curve: &Curve3, sense: Orientation) -> Option<Curve2> {
    let (frame, radius, start, sweep) = match curve {
        Curve3::Circle { frame, radius } => (frame, *radius, 0.0, std::f64::consts::TAU),
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => (frame, *radius, *start_angle, *sweep_angle),
        _ => return None,
    };
    let c = p.coordinates(frame.origin());
    let (xx, xy) = (frame.x().dot(p.x()), frame.x().dot(p.y()));
    let (yx, yy) = (frame.y().dot(p.x()), frame.y().dot(p.y()));
    let turn = if xx * yy - xy * yx < 0.0 { -1.0 } else { 1.0 };
    let (start, sweep) = match sense {
        Orientation::Forward => (start, sweep),
        Orientation::Reversed => (start + sweep, -sweep),
    };
    Some(Curve2::CircularArc {
        center: Point2::new(c[0], c[1]),
        radius,
        start_angle: xy.atan2(xx) + turn * start,
        sweep_angle: turn * sweep,
    })
}

/// How far a part's end may reach: to a cap's height, a step's narrower
/// side's into the wider side, its apex, or past a corner.
#[derive(Debug, Clone, Copy)]
enum Reach {
    Cap(f64),
    Near(f64),
    Apex,
    Corner,
}

/// A primitive of the lattice: its line, heights and their ends'
/// reach, its range (exact; `None` at a cone's apex, its line's sign bounding
/// it) and its built heights (binary64).
#[derive(Debug, Clone)]
struct Part {
    line: Line,
    interval: [R; 2],
    reach: [Reach; 2],
    range: [Option<R>; 2],
    built: [f64; 2],
    /// S9e.4b.4d.2c.1b: a cone tangent to a torus, its primitive's radii
    /// exact.
    exact: bool,
}

impl Part {
    /// Its built heights `[lo, hi]`, a cone clamped at its apex.
    fn set(&mut self, built: [f64; 2]) {
        self.built = built;
        self.range = [Some(q(built[0])), Some(q(built[1]))];
        if let Some(apex) = self.line.apex() {
            let at = rational_f64(&apex);
            let up = self.line.k > zero();
            let end = usize::from(!up);
            let past = if up { built[0] <= at } else { built[1] >= at };
            if past || matches!(self.reach[end], Reach::Apex) {
                self.built[end] = at;
                self.range[end] = None;
            }
        }
    }

    /// Whether its range covers heights.
    fn covers(&self, [a, b]: &[R; 2]) -> bool {
        self.range[0].as_ref().is_none_or(|lo| lo <= a)
            && self.range[1].as_ref().is_none_or(|hi| hi >= b)
    }

    /// Its value at a height (`None`: empty there).
    fn value(&self, w: &R) -> Option<R> {
        if self.range[0].as_ref().is_some_and(|lo| w < lo)
            || self.range[1].as_ref().is_some_and(|hi| w > hi)
        {
            return None;
        }
        let r = self.line.at(w);
        (r > zero()).then_some(r)
    }
}

/// The max-min form at a height.
fn form_at(parts: &[Part], terms: &[Vec<usize>], w: &R) -> Option<R> {
    let mut best: Option<R> = None;
    for term in terms {
        let least = term
            .iter()
            .map(|&j| parts[j].value(w))
            .try_fold(None::<R>, |least, r| {
                let r = r?;
                Some(Some(match least {
                    Some(l) if l <= r => l,
                    _ => r,
                }))
            });
        if let Some(Some(r)) = least {
            if best.as_ref().is_none_or(|b| r > *b) {
                best = Some(r);
            }
        }
    }
    best
}

/// A turned body of corners as the lattice of its primitives, or
/// `not_their_chain` where it is no such body or no form of its primitives
/// builds it.
#[allow(clippy::too_many_lines)]
pub(super) fn turned_piece(t: &Topology, tolerance: Tolerance, op: OperationId) -> Result<Piece> {
    let tol = tolerance.linear();
    let m = meridian(t, tol).ok_or_else(not_their_chain)?;
    // S9e.4b.4d.2c.1a: the lattice of the sharp profile (each fillet and
    // round its lines' corner), its pieces added and taken away below.
    let (path, nodes) = m.sharp();
    let (nodes, path) = (&nodes, &path);
    let last = nodes.len() - 1;
    let (w_min, w_max) = (rational_f64(&nodes[0][1]), rational_f64(&nodes[last][1]));
    let r_max = nodes
        .iter()
        .map(|x| rational_f64(&x[0]))
        .fold(0.0f64, f64::max);
    let far = 0.25 * (w_max - w_min).max(r_max);
    let up = m.axis.normal();
    // Each curved face's part and how far each end may reach.
    let mut parts: Vec<Part> = Vec::new();
    for (s, &i) in path.iter().enumerate() {
        let Kind::Curved(line) = &m.kinds[i] else {
            continue;
        };
        let mut reach = [Reach::Corner, Reach::Corner];
        for (end, beyond) in [(0usize, s.checked_sub(1)), (1, Some(s + 1))] {
            if nodes[s + end][0] == zero() {
                reach[end] = Reach::Apex;
                continue;
            }
            let Some(&other) = beyond.and_then(|b| path.get(b)) else {
                return Err(not_their_chain());
            };
            let Kind::Plane(h) = &m.kinds[other] else {
                continue;
            };
            let sign = if end == 0 { -1.0 } else { 1.0 };
            let facing = outward(t, other).ok_or_else(general)?.dot(up) * sign;
            reach[end] = if facing > 0.0 {
                Reach::Cap(rational_f64(h))
            } else {
                // Into the wider side beyond the step.
                let wider = if end == 0 {
                    s.checked_sub(2)
                } else {
                    Some(s + 2)
                };
                let Some(wider) = wider.filter(|&x| x < path.len()) else {
                    return Err(not_their_chain());
                };
                let span = rational_f64(&(&nodes[wider + 1][1] - &nodes[wider][1]));
                Reach::Near(rational_f64(h) + sign * 0.25 * span)
            };
        }
        let mut part = Part {
            line: line.clone(),
            interval: [nodes[s][1].clone(), nodes[s + 1][1].clone()],
            reach,
            range: [None, None],
            built: [0.0; 2],
            exact: m.exact[i],
        };
        let potential = [0, 1].map(|end| match reach[end] {
            Reach::Cap(h) | Reach::Near(h) => h,
            Reach::Apex | Reach::Corner if end == 0 => w_min - far,
            Reach::Apex | Reach::Corner => w_max + far,
        });
        part.set(potential);
        parts.push(part);
    }
    // S9e.4b.4d.2c.1b: one part where the body has fillets or rounds (its
    // tree their pieces').
    if parts.len() < 2 && m.corners.is_empty() {
        return Err(not_their_chain());
    }
    // Each curved face's term: the parts on or outside its line over its
    // heights whose ranges may cover them.
    let mut every: Vec<Vec<usize>> = parts
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            let [a, b] = &pi.interval;
            (0..parts.len())
                .filter(|&j| {
                    let pj = &parts[j];
                    j == i
                        || (pj.covers(&pi.interval)
                            && pj.line.at(a) >= pi.line.at(a)
                            && pj.line.at(b) >= pi.line.at(b))
                })
                .collect()
        })
        .collect();
    // Each term then made as small as the form allows: a member dropped
    // where the form, its ranges taken again, is still the profile (in the
    // faces' order, the later members first): fewer Booleans and shorter
    // primitives.
    let span = (w_min, w_max);
    let mut terms = settle(&mut parts, &every, nodes, span);
    for i in 0..every.len() {
        for j in every[i].clone().into_iter().rev() {
            if j == i {
                continue;
            }
            let kept = every[i].clone();
            every[i].retain(|&x| x != j);
            let fewer = settle(&mut parts, &every, nodes, span);
            if check_form(&parts, &fewer, nodes, path, &m.kinds).is_ok() {
                terms = fewer;
            } else {
                every[i] = kept;
                terms = settle(&mut parts, &every, nodes, span);
            }
        }
    }
    check_form(&parts, &terms, nodes, path, &m.kinds)?;
    // The common of every term's parts taken out.
    let common: Vec<usize> = (0..parts.len())
        .filter(|j| terms.iter().all(|x| x.contains(j)))
        .collect();
    let chain = |op2: Op2, xs: &mut dyn Iterator<Item = Tree>| -> Option<Tree> {
        let first = xs.next()?;
        Some(xs.fold(first, |a, b| Tree::Op(op2, Box::new(a), Box::new(b))))
    };
    let mut rest = Vec::new();
    if terms.len() > 1 {
        for x in &terms {
            let mut own = x
                .iter()
                .filter(|j| !common.contains(j))
                .map(|&j| Tree::Primitive(j));
            rest.push(chain(Op2::Common, &mut own).ok_or_else(not_their_chain)?);
        }
    }
    let fused = chain(Op2::Fuse, &mut rest.into_iter());
    let mut leaves = fused
        .into_iter()
        .chain(common.iter().map(|&j| Tree::Primitive(j)));
    let mut tree = chain(Op2::Common, &mut leaves).ok_or_else(not_their_chain)?;
    if !matches!(tree, Tree::Op(..)) && m.corners.is_empty() {
        return Err(not_their_chain());
    }
    // Each part's primitive on the axis frame.
    let mut solids = Vec::new();
    for (j, p) in parts.iter().enumerate() {
        solids.push(primitive(&m.axis, p, tolerance, mine(op, j))?);
    }
    // S9e.4b.4d.2c.1a: each fillet's piece fused and each round's cut: the
    // box from its tube's centre past its corner (into the lattice for a
    // fillet, out of it for a round) less its whole torus; S9e.4b.4d.2c.1b:
    // beside a cone the slant box's solids' tree.
    let ends: Vec<f64> = parts.iter().flat_map(|p| p.built).collect();
    let boxes = corner_boxes(&m, nodes, &ends)?;
    let annulus = |radii: [f64; 2], heights: [f64; 2], op: OperationId| -> Result<Solid> {
        let profile = Profile::new(
            Boundary::circle(Point2::new(0.0, 0.0), radii[1], tolerance)?,
            vec![Boundary::circle(
                Point2::new(0.0, 0.0),
                radii[0],
                tolerance,
            )?],
            tolerance,
        )?;
        Solid::build(op, profile, m.axis, heights[0], heights[1])
    };
    for (c, b) in m.corners.iter().zip(&boxes) {
        let box_tree = match b {
            CornerBox::Prism(prism) => {
                let [lo, hi] = &**prism;
                let j = solids.len();
                solids.push(annulus(
                    [rational_f64(&lo[0]), rational_f64(&hi[0])],
                    [rational_f64(&lo[1]), rational_f64(&hi[1])],
                    mine(op, j),
                )?);
                Tree::Primitive(j)
            }
            // S9e.4b.4d.2c.1b: beside a cone, its slant box's solids.
            CornerBox::Slant(parts) => {
                let mut leaves = Vec::new();
                for part in parts {
                    let j = solids.len();
                    solids.push(match part {
                        SlantPart::Annulus { radii, heights } => {
                            annulus(*radii, *heights, mine(op, j))?
                        }
                        SlantPart::Frustum { heights, radii } => Solid::build_cone_exact_on(
                            mine(op, j),
                            m.axis,
                            *heights,
                            radii.clone(),
                            tolerance,
                        )?,
                    });
                    leaves.push(j);
                }
                slant_tree(&leaves)
            }
        };
        let j = solids.len();
        solids.push(Solid::build_torus_on(
            mine(op, j),
            m.axis,
            rational_f64(&c.centre[1]),
            rational_f64(&c.centre[0]),
            rational_f64(&c.radius),
            tolerance,
        )?);
        let piece = Tree::Op(Op2::Cut, Box::new(box_tree), Box::new(Tree::Primitive(j)));
        let op2 = if c.round { Op2::Cut } else { Op2::Fuse };
        tree = Tree::Op(op2, Box::new(tree), Box::new(piece));
    }
    let mut solids = solids.into_iter();
    let primitive = solids.next().ok_or_else(not_their_chain)?;
    Ok(Piece {
        primitive: Box::new(primitive),
        others: solids.collect(),
        planes: Vec::new(),
        form: Form::Tree,
        tree: Some(tree),
        choice: TURNED,
        chain: true,
        stack: false,
        nested: false,
    })
}

/// The terms of the faces' terms `every` (each once, a term holding
/// another's dropped: it lies below it everywhere) and each part's range:
/// past a corner the heights of the faces whose terms hold it and a quarter
/// of its own height beyond, at a step's narrower side its reach into the
/// wider side; such an end moved on by a thirty-second of the body's height
/// (`span`) while it lies nearer than that to a node, a cap, an apex or an
/// earlier part's end (two ends on one plane are refused by S9's rule, and a
/// thin band of a face between two nearby ones by the engine's grouping).
fn settle(
    parts: &mut [Part],
    every: &[Vec<usize>],
    nodes: &[[R; 2]],
    (w_min, w_max): (f64, f64),
) -> Vec<Vec<usize>> {
    let mut terms: Vec<Vec<usize>> = Vec::new();
    for x in every {
        if !terms.contains(x) {
            terms.push(x.clone());
        }
    }
    let all = terms.clone();
    terms.retain(|x| {
        !all.iter()
            .any(|y| y != x && y.iter().all(|j| x.contains(j)))
    });
    let intervals: Vec<[R; 2]> = parts.iter().map(|p| p.interval.clone()).collect();
    let apart = (w_max - w_min) / 32.0;
    let mut taken: Vec<f64> = nodes.iter().map(|x| rational_f64(&x[1])).collect();
    for p in parts.iter() {
        taken.extend(p.line.apex().map(|a| rational_f64(&a)));
        for how in p.reach {
            if let Reach::Cap(h) = how {
                taken.push(h);
            }
        }
    }
    for (j, p) in parts.iter_mut().enumerate() {
        let (mut lo, mut hi) = (p.interval[0].clone(), p.interval[1].clone());
        for (i, x) in every.iter().enumerate() {
            if x.contains(&j) {
                lo = lo.min(intervals[i][0].clone());
                hi = hi.max(intervals[i][1].clone());
            }
        }
        let margin = 0.25 * rational_f64(&(&p.interval[1] - &p.interval[0]));
        let want = [rational_f64(&lo) - margin, rational_f64(&hi) + margin];
        let mut built = [0.0; 2];
        for end in 0..2 {
            let sign = if end == 0 { -1.0 } else { 1.0 };
            built[end] = match p.reach[end] {
                Reach::Cap(h) | Reach::Near(h) => h,
                Reach::Apex | Reach::Corner => want[end],
            };
            if !matches!(p.reach[end], Reach::Cap(_)) {
                for _ in 0..64 {
                    if taken.iter().all(|&x| (x - built[end]).abs() >= apart) {
                        break;
                    }
                    built[end] += sign * apart;
                }
            }
        }
        p.set(built);
        taken.extend(p.built);
    }
    terms
}

/// The max-min form is the profile at every height: between every two of
/// its breakpoints (the nodes' heights, the ranges' ends, the apexes and the
/// lines' crossings) at two points, and beyond them, empty outside the body.
fn check_form(
    parts: &[Part],
    terms: &[Vec<usize>],
    nodes: &[[R; 2]],
    path: &[usize],
    kinds: &[Kind],
) -> Result<()> {
    let mut breaks: Vec<R> = nodes.iter().map(|x| x[1].clone()).collect();
    for p in parts {
        breaks.extend(p.range.iter().flatten().cloned());
        breaks.extend(p.line.apex());
    }
    for (a, x) in parts.iter().enumerate() {
        for y in &parts[a + 1..] {
            if x.line.k != y.line.k {
                breaks.push((&y.line.r0 - &x.line.r0) / (&x.line.k - &y.line.k));
            }
        }
    }
    breaks.sort();
    breaks.dedup();
    let int = |n: i64| R::from_integer(n.into());
    let mut points = vec![&breaks[0] - int(1), &breaks[breaks.len() - 1] + int(1)];
    for x in breaks.windows(2) {
        let d = &x[1] - &x[0];
        points.push(&x[0] + &d / int(3));
        points.push(&x[0] + &(&d * int(2)) / int(3));
    }
    let profile = |w: &R| -> Option<R> {
        path.iter().enumerate().find_map(|(s, &i)| match &kinds[i] {
            Kind::Curved(l) if nodes[s][1] < *w && *w < nodes[s + 1][1] => Some(l.at(w)),
            _ => None,
        })
    };
    if points
        .iter()
        .any(|w| form_at(parts, terms, w) != profile(w))
    {
        return Err(not_their_chain());
    }
    Ok(())
}

/// S9e.4b.4d.2c.1a: each corner's box in the meridian, `[(r, w) low, (r, w)
/// high]`: from its tube's centre past its corner by half the tube (halved
/// up to four times), on the grid of its numbers; the profile's other
/// segments and the other boxes clear of it, the two lines' segments at the
/// corner leaving it, its heights none of the lattice's (`ends`, the parts'
/// built heights) nor its nodes'. S9e.4b.4d.2c.1b: beside a cone, its slant
/// box (`slant_box`).
fn corner_boxes(m: &Meridian, nodes: &[[R; 2]], ends: &[f64]) -> Result<Vec<CornerBox>> {
    let mut out: Vec<[[R; 2]; 2]> = Vec::new();
    let mut boxes = Vec::new();
    for c in &m.corners {
        if let Some(e) = &c.slant {
            let (parts, rect) = slant_box(c, e, nodes, &out, ends).ok_or_else(not_their_chain)?;
            out.push(rect);
            boxes.push(CornerBox::Slant(parts));
            continue;
        }
        let toward = |k: usize| {
            if c.corner[k] > c.centre[k] {
                int(1)
            } else {
                int(-1)
            }
        };
        let mut found = None;
        let mut d = &c.radius / int(2);
        for _ in 0..5 {
            let far = [&c.corner[0] + toward(0) * &d, &c.corner[1] + toward(1) * &d];
            let lo = [
                c.centre[0].clone().min(far[0].clone()),
                c.centre[1].clone().min(far[1].clone()),
            ];
            let hi = [
                c.centre[0].clone().max(far[0].clone()),
                c.centre[1].clone().max(far[1].clone()),
            ];
            let inside =
                |x: &[R; 2]| lo[0] <= x[0] && x[0] <= hi[0] && lo[1] <= x[1] && x[1] <= hi[1];
            let clear = nodes.windows(2).all(|seg| {
                let at_corner = seg[0] == c.corner || seg[1] == c.corner;
                if at_corner {
                    // A line at the corner leaves the box.
                    let other = if seg[0] == c.corner { &seg[1] } else { &seg[0] };
                    !inside(other)
                } else {
                    !segment_meets(&seg[0], &seg[1], &lo, &hi)
                }
            }) && out
                .iter()
                .all(|[a, b]| b[0] < lo[0] || hi[0] < a[0] || b[1] < lo[1] || hi[1] < a[1])
                && [&c.centre[1], &far[1]].iter().all(|w| {
                    let w = rational_f64(w);
                    !ends.contains(&w) && !nodes.iter().any(|x| rational_f64(&x[1]) == w)
                })
                && lo[0] > zero();
            if clear {
                found = Some([lo, hi]);
                break;
            }
            d = &d / int(2);
        }
        let found = found.ok_or_else(not_their_chain)?;
        out.push(found.clone());
        boxes.push(CornerBox::Prism(Box::new(found)));
    }
    Ok(boxes)
}

/// A corner's box: S9e.4b.4d.2c.1a's annular prism `[(r, w) low, (r, w)
/// high]`, or S9e.4b.4d.2c.1b's slant box's solids beside a cone.
enum CornerBox {
    Prism(Box<[[R; 2]; 2]>),
    Slant(Vec<SlantPart>),
}

/// S9e.4b.4d.2c.1b: a slant box's solid on the axis frame: an annular prism
/// (its radii and heights), or a frustum of a rational line (its heights
/// and its exact radii there).
enum SlantPart {
    Annulus { radii: [f64; 2], heights: [f64; 2] },
    Frustum { heights: [f64; 2], radii: [R; 2] },
}

/// S9e.4b.4d.2c.1b: a slant box's tree over its solids' leaves: the outer
/// frustum less the inner (beside a cylinder), or the annulus common the
/// outer less the inner (beside a plane).
fn slant_tree(leaves: &[usize]) -> Tree {
    let p = |j: usize| Box::new(Tree::Primitive(leaves[j]));
    if leaves.len() == 2 {
        return Tree::Op(Op2::Cut, p(0), p(1));
    }
    Tree::Op(Op2::Cut, Box::new(Tree::Op(Op2::Common, p(0), p(1))), p(2))
}

/// S9e.4b.4d.2c.1b: a corner's slant box beside a cone: the parallelogram
/// in the meridian from the tube's centre along its directions `e` (to the
/// cone, `x`; to the other line, `y`) past the corner by half the tube
/// (halved up to four times), its far plane or cylinder and its far cone's
/// point on the near plane or cylinder on the multiples of `2^-24` of the
/// tube's power of two, outward; the profile's other segments and the
/// other boxes (`rects`, their bounding rectangles) clear of it, the two
/// lines at the corner leaving it, its solids' heights none of the
/// lattice's nor the nodes'. Its solids and its bounding rectangle.
fn slant_box(
    c: &Corner,
    e: &[[R; 2]; 2],
    nodes: &[[R; 2]],
    rects: &[[[R; 2]; 2]],
    ends: &[f64],
) -> Option<(Vec<SlantPart>, [[R; 2]; 2])> {
    let [ec, ef] = e;
    let det = &ec[0] * &ef[1] - &ec[1] * &ef[0];
    if det == zero() {
        return None;
    }
    let oblique = |p: &[R; 2]| -> [R; 2] {
        let v = [&p[0] - &c.centre[0], &p[1] - &c.centre[1]];
        [
            (&v[0] * &ef[1] - &v[1] * &ef[0]) / &det,
            (&ec[0] * &v[1] - &ec[1] * &v[0]) / &det,
        ]
    };
    let corner = oblique(&c.corner);
    if corner[0] <= zero() || corner[1] <= zero() {
        return None;
    }
    let g = q(2f64.powi(rational_f64(&c.radius).log2().ceil() as i32 - 24));
    // The coordinate the far side along `ef` fixes (a plane's height beside
    // a cylinder, a cylinder's radius beside a plane), and the other.
    let k = if ef[1] == zero() { 1 } else { 0 };
    let j = 1 - k;
    let mut d = &c.radius / int(2);
    for _ in 0..5 {
        let far = &c.centre[k] + (&corner[0] + &d) * &ec[k];
        let x = (grid_toward(&far, &g, ec[k] > zero()) - &c.centre[k]) / &ec[k];
        let near = &c.centre[j] + (&corner[1] + &d) * &ef[j];
        let y = (grid_toward(&near, &g, ef[j] > zero()) - &c.centre[j]) / &ef[j];
        let inside = |p: &[R; 2]| zero() <= p[0] && p[0] <= x && zero() <= p[1] && p[1] <= y;
        let segments = nodes.windows(2).all(|seg| {
            if seg[0] == c.corner || seg[1] == c.corner {
                let other = if seg[0] == c.corner { &seg[1] } else { &seg[0] };
                !inside(&oblique(other))
            } else {
                let box_hi = [x.clone(), y.clone()];
                !segment_meets(
                    &oblique(&seg[0]),
                    &oblique(&seg[1]),
                    &[zero(), zero()],
                    &box_hi,
                )
            }
        });
        let at = |a: &R, b: &R| {
            [
                &c.centre[0] + a * &ec[0] + b * &ef[0],
                &c.centre[1] + a * &ec[1] + b * &ef[1],
            ]
        };
        let vertices = [
            at(&zero(), &zero()),
            at(&x, &zero()),
            at(&zero(), &y),
            at(&x, &y),
        ];
        let least = |k: usize| vertices.iter().map(|p| p[k].clone()).min().expect("four");
        let most = |k: usize| vertices.iter().map(|p| p[k].clone()).max().expect("four");
        let (lo, hi) = ([least(0), least(1)], [most(0), most(1)]);
        let others = rects
            .iter()
            .all(|[a, b]| b[0] < lo[0] || hi[0] < a[0] || b[1] < lo[1] || hi[1] < a[1]);
        let parts = slant_parts(c, e, k, (&x, &y), [&lo[1], &hi[1]], &g);
        let heights = parts.as_ref().is_some_and(|ps| {
            ps.iter().all(|p| {
                let (SlantPart::Annulus { heights, .. } | SlantPart::Frustum { heights, .. }) = p;
                heights
                    .iter()
                    .all(|w| !ends.contains(w) && !nodes.iter().any(|x| rational_f64(&x[1]) == *w))
            })
        });
        if segments && others && heights && lo[0] > zero() {
            return Some((parts?, [lo, hi]));
        }
        d = &d / int(2);
    }
    None
}

/// `x` on the multiples of `g`, up or down.
fn grid_toward(x: &R, g: &R, up: bool) -> R {
    let n = x / g;
    (if up { n.ceil() } else { n.floor() }) * g
}

/// S9e.4b.4d.2c.1b: a slant box's solids: the cones normal to the cone (of
/// its tube's direction `e[0]`) through the centre and through the far
/// point `y` along `e[1]`, the outer's frustum less the inner's (past the
/// box's ends by an eighth of the tube); beside a cylinder (`k` one)
/// between the plane through the centre and the far plane `x` along
/// `e[0]`, beside a plane over the box's heights `w` on the grid `g`,
/// outward, common the annulus between the cylinder through the centre and
/// the far one (past them by an eighth of the tube, the frusta by two and
/// three). Every height a binary64 number, every radius positive.
fn slant_parts(
    c: &Corner,
    [ec, ef]: &[[R; 2]; 2],
    k: usize,
    (x, y): (&R, &R),
    [w_lo, w_hi]: [&R; 2],
    g: &R,
) -> Option<Vec<SlantPart>> {
    let eps = &c.radius / int(8);
    let slope = &ec[0] / &ec[1];
    let radius_at = |p: &[R; 2], h: &R| &p[0] + &slope * (h - &p[1]);
    let near = c.centre.clone();
    let far = [&c.centre[0] + y * &ef[0], &c.centre[1] + y * &ef[1]];
    let frustum = |p: &[R; 2], [a, b]: [R; 2]| -> Option<SlantPart> {
        let radii = [radius_at(p, &a), radius_at(p, &b)];
        let heights = [rational_f64(&a), rational_f64(&b)];
        let exact = q(heights[0]) == a && q(heights[1]) == b && a < b;
        (exact && radii.iter().all(|r| *r > zero()))
            .then_some(SlantPart::Frustum { heights, radii })
    };
    let (outer, inner) = if radius_at(&far, &c.centre[1]) > c.centre[0] {
        (&far, &near)
    } else {
        (&near, &far)
    };
    if k == 1 {
        let wf = &c.centre[1] + x * &ec[1];
        let (a, b) = if wf > c.centre[1] {
            (c.centre[1].clone(), wf)
        } else {
            (wf, c.centre[1].clone())
        };
        return Some(vec![
            frustum(outer, [a.clone(), b.clone()])?,
            frustum(inner, [&a - &eps, &b + &eps])?,
        ]);
    }
    let rf = &c.centre[0] + x * &ec[0];
    let radii = if rf > c.centre[0] {
        [c.centre[0].clone(), rf]
    } else {
        [rf, c.centre[0].clone()]
    };
    let (a, b) = (grid_toward(w_lo, g, false), grid_toward(w_hi, g, true));
    let two = &eps * int(2);
    let three = &eps * int(3);
    let annulus = [&a - &eps, &b + &eps];
    let binary = |x: &R| (q(rational_f64(x)) == *x).then(|| rational_f64(x));
    Some(vec![
        SlantPart::Annulus {
            radii: [binary(&radii[0])?, binary(&radii[1])?],
            heights: [binary(&annulus[0])?, binary(&annulus[1])?],
        },
        frustum(outer, [&a - &two, &b + &two])?,
        frustum(inner, [&a - &three, &b + &three])?,
    ])
}

/// Whether the segment from `a` to `b` meets the closed box `[lo, hi]`
/// (exactly: its parameter's range clipped by the box's four sides).
fn segment_meets(a: &[R; 2], b: &[R; 2], lo: &[R; 2], hi: &[R; 2]) -> bool {
    let (mut t0, mut t1) = (zero(), int(1));
    for k in 0..2 {
        let d = &b[k] - &a[k];
        if d == zero() {
            if a[k] < lo[k] || a[k] > hi[k] {
                return false;
            }
            continue;
        }
        let (mut e0, mut e1) = ((&lo[k] - &a[k]) / &d, (&hi[k] - &a[k]) / &d);
        if e0 > e1 {
            std::mem::swap(&mut e0, &mut e1);
        }
        t0 = t0.max(e0);
        t1 = t1.min(e1);
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// A part's primitive on the axis frame: a cylinder's prism of its circle,
/// a cone's frustum from its lower height, its exact model on the axis
/// (its radii its line's at its heights, rounded once).
fn primitive(axis: &Frame3, p: &Part, tolerance: Tolerance, op: OperationId) -> Result<Solid> {
    let [lo, hi] = p.built;
    if p.line.k == zero() {
        let profile = Profile::new(
            Boundary::circle(Point2::new(0.0, 0.0), rational_f64(&p.line.r0), tolerance)?,
            Vec::new(),
            tolerance,
        )?;
        return Solid::build(op, profile, *axis, lo, hi);
    }
    // S9e.4b.4d.2c.1b: a cone tangent to a torus, its line's radii exactly.
    if p.exact {
        let exact = |end: usize, w: f64| match &p.range[end] {
            None => zero(),
            Some(_) => p.line.at(&q(w)),
        };
        let radii = [exact(0, lo), exact(1, hi)];
        return Solid::build_cone_exact_on(op, *axis, [lo, hi], radii, tolerance);
    }
    let radius = |end: usize, w: f64| match &p.range[end] {
        None => 0.0,
        Some(_) => rational_f64(&p.line.at(&q(w))),
    };
    let (bottom, top) = (radius(0, lo), radius(1, hi));
    Solid::build_cone_on(op, *axis, [lo, hi], bottom, top, tolerance)
}

/// The axis frame: the stored frame with its axes' components that move the
/// body (of `size`) by less than a thousandth of the resolution taken to
/// zero (a stored axis within rounding of the world's: its exact data
/// short).
fn snapped(frame: Frame3, size: f64, tol: f64) -> Option<Frame3> {
    let small = tol / 1024.0 / size.max(1.0);
    let cut = |v: Vec3| {
        let a = v.to_array().map(|c| if c.abs() <= small { 0.0 } else { c });
        Vec3::new(a[0], a[1], a[2])
    };
    Frame3::new(
        frame.origin(),
        cut(frame.normal()),
        cut(frame.x()),
        Tolerance::default(),
    )
    .ok()
}
