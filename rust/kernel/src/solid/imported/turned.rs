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
use super::*;

/// The choice a turned body's piece records (`Imported::moved` reads it off
/// the moved stored topology again by it).
pub(super) const TURNED: u32 = u32::MAX;

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

/// A stored face's meridian: a curved face's line, a plane's height.
#[derive(Debug, Clone)]
enum Kind {
    Curved(Line),
    Plane(R),
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
    for (i, f) in faces.iter().enumerate() {
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
    // Each face's two ends: its rings, a cone's apex, a disc's centre.
    let mut ends: Vec<Vec<End>> = vec![Vec::new(); faces.len()];
    for (e, (a, b, _, _)) in rings.iter().enumerate() {
        ends[*a].push(End::Ring(e));
        ends[*b].push(End::Ring(e));
    }
    for (i, f) in faces.iter().enumerate() {
        for l in &f.loops {
            if let Loop::Vertex(v) = &t.loops()[l.index()] {
                if !matches!(kinds[i], Kind::Curved(_)) {
                    return None;
                }
                ends[i].push(End::Apex(v.index()));
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
    // Each node exactly: a corner the two meridians' common point, within
    // the resolution of its stored circle.
    let mut nodes: Vec<[R; 2]> = Vec::new();
    for (s, end) in path_ends.iter().enumerate() {
        let own = &kinds[path[s.min(path.len() - 1)]];
        let node = match end {
            End::Ring(e) => {
                let node = match (&kinds[path[s - 1]], own) {
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
    if nodes[0][1] > nodes[nodes.len() - 1][1] {
        nodes.reverse();
        path.reverse();
    }
    for (s, &i) in path.iter().enumerate() {
        let (a, b) = (&nodes[s][1], &nodes[s + 1][1]);
        match kinds[i] {
            Kind::Curved(_) if rational_f64(&(b - a)) > tol => {}
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
            Kind::Plane(_) => None,
        })
        .collect();
    for (a, x) in lines.iter().enumerate() {
        if lines[a + 1..].contains(x) {
            return None;
        }
    }
    Some(Meridian {
        axis,
        kinds,
        path,
        nodes,
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
    let (nodes, path) = (&m.nodes, &m.path);
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
        };
        let potential = [0, 1].map(|end| match reach[end] {
            Reach::Cap(h) | Reach::Near(h) => h,
            Reach::Apex | Reach::Corner if end == 0 => w_min - far,
            Reach::Apex | Reach::Corner => w_max + far,
        });
        part.set(potential);
        parts.push(part);
    }
    if parts.len() < 2 {
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
    let tree = chain(Op2::Common, &mut leaves).ok_or_else(not_their_chain)?;
    if !matches!(tree, Tree::Op(..)) {
        return Err(not_their_chain());
    }
    // Each part's primitive on the axis frame.
    let mut solids = Vec::new();
    for (j, p) in parts.iter().enumerate() {
        solids.push(primitive(&m.axis, p, tolerance, mine(op, j))?);
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
    let radius = |end: usize, w: f64| match &p.range[end] {
        None => 0.0,
        Some(_) => rational_f64(&p.line.at(&q(w))),
    };
    let (bottom, top) = (radius(0, lo), radius(1, hi));
    Solid::build_cone_on(op, *axis, lo, bottom, top, hi - lo, tolerance)
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
