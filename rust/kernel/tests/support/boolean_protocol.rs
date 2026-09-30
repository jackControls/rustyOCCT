//! The case protocol of `boolean-cases.txt` (S9a) and the kernel's rows for
//! it, shared by `boolean_probe` and `tests/booleans.rs`: a case is the
//! object's prism rows, `boolean fuse|cut|common ID`, and the tool's rows
//! (its `case` row the object's). Rows: `limit`, `unsupported` (a stack
//! before S9a.2), `refused` (a documented `Degenerate`), `empty`, or per
//! result solid `solid vol_lo vol_hi area_lo area_hi cx_lo cx_hi cy_lo cy_hi
//! cz_lo cz_hi faces edges vertices` (`Topology::occt_counts`).
//!
//! S9e.1: a case may chain a second Boolean on the first's result: after
//! the tool's rows a row `then fuse|cut|common ID` (`swapped` after it when
//! the first result is the tool) and a third solid's rows. The first
//! Boolean's result must be one solid; the rows are the second's. S9e.2: the
//! `then` row may end `solid X Y Z`: the first result may hold several
//! solids, and the one holding the point strictly inside (`Solid::classify`,
//! exactly one) is the second's argument.
#[path = "identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
pub use identity_protocol::{build, parse, CaseSpec};
use rusty_occt::history::History;
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Solid};

pub struct Case {
    pub name: String,
    pub object: CaseSpec,
    pub tool: CaseSpec,
    pub op: String,
    pub operation: OperationId,
    /// A second Boolean on the first's result (S9e.1).
    pub then: Option<Then>,
}

/// The second Boolean of a chained case: its operation, id, third solid
/// and whether the first result is its tool.
pub struct Then {
    pub op: String,
    pub operation: OperationId,
    pub third: CaseSpec,
    pub swapped: bool,
    /// The first result's solid holding this point (S9e.2), none for its
    /// one solid.
    pub pick: Option<[f64; 3]>,
}

pub fn cases(text: &str) -> Vec<Case> {
    text.split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|block| {
            let lines: Vec<&str> = block.trim().lines().collect();
            let at = lines
                .iter()
                .position(|l| l.starts_with("boolean "))
                .expect("a boolean row");
            let w: Vec<&str> = lines[at].split_whitespace().collect();
            let object = parse(&lines[..at].join("\n"));
            let then_at = lines
                .iter()
                .position(|l| l.starts_with("then "))
                .unwrap_or(lines.len());
            let tool_rows: Vec<&str> = std::iter::once(lines[0])
                .chain(lines[at + 1..then_at].iter().copied())
                .collect();
            let tool = parse(&tool_rows.join("\n"));
            let then = (then_at < lines.len()).then(|| {
                let t: Vec<&str> = lines[then_at].split_whitespace().collect();
                let rows: Vec<&str> = std::iter::once(lines[0])
                    .chain(lines[then_at + 1..].iter().copied())
                    .collect();
                Then {
                    op: t[1].to_string(),
                    operation: OperationId(t[2].parse().expect("an operation id")),
                    third: parse(&rows.join("\n")),
                    swapped: t.get(3) == Some(&"swapped"),
                    pick: t.iter().position(|w| *w == "solid").map(|k| {
                        [1, 2, 3].map(|i| t[k + i].parse::<f64>().expect("a pick point"))
                    }),
                }
            });
            Case {
                name: object.name.clone(),
                object,
                tool,
                op: w[1].to_string(),
                operation: OperationId(w[2].parse().expect("an operation id")),
                then,
            }
        })
        .collect()
}

/// The inputs, the results and the history of a case.
pub type Run = (Solid, Solid, Vec<Solid>, History);

pub fn run(case: &Case) -> Result<Run, Error> {
    let (a, b, out, history) = run_first(case)?;
    let Some(then) = &case.then else {
        return Ok((a, b, out, history));
    };
    // S9e.1: the second Boolean on the first's one solid (S9e.2: the one
    // holding the pick point).
    let first = given(case, out);
    let third = build(&then.third);
    let (a, b) = if then.swapped {
        (third, first)
    } else {
        (first, third)
    };
    let (out, history) = boolean(&a, &then.op, then.operation, &b)?;
    Ok((a, b, out, history))
}

/// The first result's solid the second Boolean takes: its one solid, or
/// the one holding the pick point strictly inside (S9e.2).
pub fn given(case: &Case, out: Vec<Solid>) -> Solid {
    let then = case.then.as_ref().expect("a chained case");
    let Some([x, y, z]) = then.pick else {
        let [first] = <[Solid; 1]>::try_from(out).unwrap_or_else(|_| {
            panic!("{}: the first Boolean's result is not one solid", case.name)
        });
        return first;
    };
    let point = rusty_occt::Point3::new(x, y, z);
    let mut inside: Vec<Solid> = out
        .into_iter()
        .filter(|s| {
            s.classify(point)
                .unwrap_or_else(|e| panic!("{}: the pick point: {e}", case.name))
                == rusty_occt::Location::Inside
        })
        .collect();
    assert_eq!(
        inside.len(),
        1,
        "{}: the pick point is not inside exactly one solid",
        case.name
    );
    inside.remove(0)
}

/// The first Boolean of a case (of a chained case, its first): its inputs,
/// result and history.
#[allow(dead_code)]
pub fn run_first(case: &Case) -> Result<Run, Error> {
    let (a, b) = (build(&case.object), build(&case.tool));
    let (out, history) = boolean(&a, &case.op, case.operation, &b)?;
    Ok((a, b, out, history))
}

fn boolean(
    a: &Solid,
    op: &str,
    operation: OperationId,
    b: &Solid,
) -> Result<(Vec<Solid>, History), Error> {
    match op {
        "fuse" => a.fuse(operation, b),
        "cut" => a.cut(operation, b),
        "common" => a.common(operation, b),
        other => panic!("unknown operation {other}"),
    }
}

/// The kernel's rows for a case; `Err` for an unexpected error.
pub fn rows(case: &Case) -> Result<Vec<String>, Error> {
    let out = match run(case) {
        Ok((_, _, out, _)) => out,
        Err(Error::ComputationLimit(_)) => return Ok(vec!["limit".into()]),
        Err(Error::OutOfDomain(_)) => return Ok(vec!["unsupported".into()]),
        Err(Error::Degenerate(_)) => return Ok(vec!["refused".into()]),
        Err(e) => return Err(e),
    };
    if out.is_empty() {
        return Ok(vec!["empty".into()]);
    }
    let mut rows = Vec::new();
    for solid in out {
        let t = solid.topology();
        let m = t
            .mass_enclosure()
            .ok_or(Error::ComputationLimit("a result's mass"))?;
        let c = t.occt_counts();
        let mut words = vec!["solid".to_string()];
        words.push(format!("{:?} {:?}", m.volume[0], m.volume[1]));
        words.push(format!("{:?} {:?}", m.surface_area[0], m.surface_area[1]));
        for [lo, hi] in m.centroid {
            words.push(format!("{lo:?} {hi:?}"));
        }
        words.push(format!("{} {} {}", c.faces, c.edges, c.vertices));
        rows.push(words.join(" "));
    }
    Ok(rows)
}
