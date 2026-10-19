//! S10's budget fixtures (`BUDGETS.md`): the slowest inputs of the boolean,
//! split and brep_validation corpora (`rust/fixtures/budgets/`), each run
//! with its target's oracle and S10's invariant, and the work units each of
//! its metered operations counts, pinned in `rust/fixtures/budget-work.tsv`.
//! `--check` compares (work units are the same on every host); `--write`
//! records them again, keeping each row's measured seconds.
use std::{fs, path::Path, process::ExitCode, time::Instant};

fn main() -> ExitCode {
    let mode = std::env::args().nth(1).unwrap_or_default();
    if mode != "--check" && mode != "--write" {
        eprintln!("usage: replay_budgets --check | --write");
        return ExitCode::from(2);
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures");
    let table = root.join("budget-work.tsv");
    let text = fs::read_to_string(&table).expect("budget-work.tsv");
    let mut out = String::new();
    let mut failed = 0;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let (target, input, seconds, pinned) = (cols[0], cols[1], cols[2], cols.get(3));
        let data = fs::read(
            root.join("budgets")
                .join(target)
                .join(format!("{input}.bin")),
        )
        .expect("a fixture input");
        let start = Instant::now();
        let work = rusty_occt_fuzz::budget_work(|| match target {
            "boolean" => rusty_occt_fuzz::check_boolean(&data),
            "split" => rusty_occt_fuzz::check_split(&data),
            "brep_validation" => rusty_occt_fuzz::check_brep_validation(&data),
            other => panic!("no target {other}"),
        });
        let got = work
            .iter()
            .map(|(label, units)| format!("{label}={units}"))
            .collect::<Vec<_>>()
            .join(" ");
        println!(
            "{target} {input}: {:.3}s, {got}",
            start.elapsed().as_secs_f64()
        );
        if pinned.is_some_and(|p| *p != got) {
            failed += 1;
            eprintln!(
                "{target} {input}: pinned {}, counted {got}",
                pinned.unwrap()
            );
        }
        out.push_str(&format!("{target}\t{input}\t{seconds}\t{got}\n"));
    }
    if mode == "--write" {
        fs::write(&table, out).expect("budget-work.tsv written");
        return ExitCode::SUCCESS;
    }
    if failed > 0 {
        eprintln!("{failed} fixtures differ");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
