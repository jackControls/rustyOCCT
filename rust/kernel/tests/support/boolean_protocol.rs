//! The case protocol of `boolean-cases.txt` (S9a) and the kernel's rows for
//! it, shared by `boolean_probe` and `tests/booleans.rs`: a case is the
//! object's prism rows, `boolean fuse|cut|common ID`, and the tool's rows
//! (its `case` row the object's). Rows: `limit`, `unsupported` (a stack
//! before S9a.2), `refused` (a documented `Degenerate`), `empty`, or per
//! result solid `solid vol_lo vol_hi area_lo area_hi cx_lo cx_hi cy_lo cy_hi
//! cz_lo cz_hi faces edges vertices` (`Topology::occt_counts`).
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
            let tool_rows: Vec<&str> = std::iter::once(lines[0])
                .chain(lines[at + 1..].iter().copied())
                .collect();
            let tool = parse(&tool_rows.join("\n"));
            Case {
                name: object.name.clone(),
                object,
                tool,
                op: w[1].to_string(),
                operation: OperationId(w[2].parse().expect("an operation id")),
            }
        })
        .collect()
}

/// The inputs, the results and the history of a case.
pub type Run = (Solid, Solid, Vec<Solid>, History);

pub fn run(case: &Case) -> Result<Run, Error> {
    let (a, b) = (build(&case.object), build(&case.tool));
    let (out, history) = match case.op.as_str() {
        "fuse" => a.fuse(case.operation, &b)?,
        "cut" => a.cut(case.operation, &b)?,
        "common" => a.common(case.operation, &b)?,
        other => panic!("unknown operation {other}"),
    };
    Ok((a, b, out, history))
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
