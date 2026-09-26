#!/usr/bin/env python3
"""Execute pinned, unmodified OCCT DRAW tests through Rust and/or native DRAW.

Python 3.9+, Tcl 8.5+. Unix process isolation; no third-party Python packages.
Exit zero means the explicit coverage contract holds, not that skipped tests pass.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "rust/fixtures/upstream-draw.json"
HOST = ROOT / "rust/tools/draw_bridge.tcl"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_worker():
    build = subprocess.run(
        ["cargo", "build", "--locked", "--example", "draw_worker", "--message-format=json"],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, check=True,
    )
    for line in build.stdout.splitlines():
        item = json.loads(line)
        if item.get("target", {}).get("name") == "draw_worker" and item.get("executable"):
            return Path(item["executable"])
    raise RuntimeError("Cargo did not report the DRAW worker executable")


def source_files(case, context=None):
    """DRAW group/grid begin scripts, the test, then reverse end scripts.

    A derived case (not an original OCCT test) lives outside tests/ and runs
    in the pinned group/grid `context` (e.g. tests/bugs/modalg_7)."""
    path = ROOT / case
    grid = ROOT / (context or str(Path(case).parent))
    directories = list(reversed(grid.relative_to(ROOT / "tests").parents))
    # Only ancestors below tests/ are groups/grids; tests/begin does not exist.
    directories = [ROOT / "tests" / p for p in directories if str(p) != "."]
    directories.append(grid)
    return ([p / "begin" for p in directories if (p / "begin").is_file()]
            + [path]
            + [p / "end" for p in reversed(directories) if (p / "end").is_file()])


def run_case(backend, sources, directory, worker=None, draw_exe=None,
             tclsh="tclsh", timeout=30.0, data_dirs=(), names=None, extra_env=None):
    directory.mkdir(parents=True, exist_ok=True)
    result_file = directory / "result.txt"
    result_file.unlink(missing_ok=True)  # A crashed rerun must not reuse an old pass.
    log_file = directory / "output.log"
    env = os.environ.copy()
    case_path = next(p for p in sources if p.name not in {"begin", "end"})
    try:
        group, grid, name = names or case_path.relative_to(ROOT / "tests").parts
    except ValueError:
        group, grid, name = "bugs", "bridge-self-test", case_path.name
    env.update({
        "RUSTY_DRAW_BACKEND": backend, "RUSTY_DRAW_ROOT": str(ROOT),
        "RUSTY_DRAW_WORKER": str(worker or ""),
        "RUSTY_DRAW_RESULT": str(result_file.resolve()),
        "RUSTY_DRAW_SOURCE_COUNT": str(len(sources)),
        "RUSTY_DRAW_DATA_COUNT": str(len(data_dirs)),
        "RUSTY_DRAW_GROUP": group, "RUSTY_DRAW_GRID": grid, "RUSTY_DRAW_CASE": name,
    })
    env.update(extra_env or {})
    env.update({f"RUSTY_DRAW_SOURCE_{i}": str(p.resolve()) for i, p in enumerate(sources)})
    env.update({f"RUSTY_DRAW_DATA_{i}": str(p.resolve()) for i, p in enumerate(data_dirs)})
    command = [str(draw_exe), "-b", "-f", str(HOST)] if backend == "occt" else [tclsh, str(HOST)]
    started = time.monotonic()
    timed_out = False
    with log_file.open("w") as log:
        try:
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as error:
            message = f"cannot start test runner: {error}"
            log.write(message + "\n")
            return {"backend": backend, "status": "failed", "queries": "0",
                    "seconds": round(time.monotonic() - started, 6),
                    "log": str(log_file.resolve()), "error": message}
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
    result = {"backend": backend, "status": "failed", "queries": "0",
              "seconds": round(time.monotonic() - started, 6), "log": str(log_file.resolve())}
    if timed_out:
        result.update(status="timeout", error=f"exceeded {timeout:g} seconds")
    elif process.returncode != 0:
        result["error"] = f"runner exited with code {process.returncode}"
    elif not result_file.is_file():
        result["error"] = "runner did not produce a result"
    else:
        try:
            fields = {}
            for line in result_file.read_text().splitlines():
                key, encoded = line.split(" ", 1)
                if key in fields:
                    raise ValueError(f"duplicate result field {key}")
                fields[key] = bytes.fromhex(encoded).decode("utf-8")
            required = {"status", "backend", "version", "queries", "unsupported", "missing", "error", "commands", "adjudication"}
            if fields.keys() != required or fields["backend"] != backend:
                raise ValueError("incomplete or inconsistent result record")
            if fields["status"] not in {"pass", "failed", "unsupported", "missing_fixture", "unverified", "known_failure", "unexpected_improvement", "skipped"}:
                raise ValueError("unknown result status")
            if fields["status"] == "pass" and int(fields["queries"]) <= 0:
                raise ValueError("pass without geometric observations")
            result.update(fields)
        except (ValueError, UnicodeError) as error:
            result.update(status="failed", error=f"invalid result: {error}")
    return result


# ---------------------------------------------------------------- coverage ledger
#
# Every assertion in every upstream case carries one status (TOPOLOGY_MODEL.md,
# T2). An assertion is model-independent when its value does not depend on
# how OCCT structures a body: validity verdicts, volumes, areas, points,
# images. It is structure-dependent otherwise; it counts as mapped-and-verified
# only when a mapping carries it to the cell model AND the case runs on both
# backends with its original assertions passing. Every other
# structure-dependent assertion is lost, with the reason.

# Assertion procedures and how their values depend on OCCT's structure.
MODEL_INDEPENDENT = {
    "checkshape", "checkreal", "checkview", "checktrinfo", "checkcolor", "checkdump",
    "checkpoint", "checkplatform", "checkarea", "checkgravitycenter", "checktime",
    "checkselfintersection", "checkMultilineStrings",
}
MAPPED = {
    "checknbshapes": "count synthesizer",
    "checklength": "per-use length synthesis",
}
UNMAPPED = {
    "checkfreebounds": "free boundaries count seams and degenerate edges",
    "checksection": "section wire and vertex counts follow OCCT's splitting",
    "checkmaxtol": "per-entity tolerances (enclosures are M5)",
    "checkfaults": "per-subshape fault statuses",
    "checkloc": "locations are OCCT structure",
    "checkoverlapedges": "per-edge overlap follows OCCT's edges",
    "checkcurveonsurf": "per-edge deviations follow OCCT's edges",
    "check_fsd": "OCCT persistence formats",
}
SELECTOR = "native selector"
LEDGER_SKIP = {"begin", "end", "parse.rules", "cases.list", "grids.list"}


def case_files():
    tests = ROOT / "tests"
    return sorted(p for p in tests.glob("*/*/**/*") if p.is_file() and p.name not in LEDGER_SKIP)


def assertions(text):
    """(kind, line) for every assertion call: a check* procedure at a command
    position, or a custom `puts "Error..."`/`puts "Faulty..."` report."""
    import re
    found = []
    call = re.compile(r'(?:^|[\[;{])\s*(check[A-Za-z_0-9]+)\b')
    custom = re.compile(r'(?:^|[\[;{])\s*puts\s+"?\s*(Error|Faulty|ERROR)')
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        for match in call.finditer(line):
            found.append((match.group(1), stripped))
        if custom.search(line):
            found.append(("custom", stripped))
    return found


def classify(kind, line, picks):
    """(structure dependence, mapping or loss reason) of one assertion."""
    import re
    if picks and (kind == "custom" or any(re.search(rf"\b{re.escape(b)}_\d", line) for b in picks)):
        return "structure", SELECTOR
    if kind == "checkprops":
        # Lengths sum per edge use; volumes, areas and centres do not.
        return ("structure", "per-use length synthesis") if re.search(r"\s-l\b", line) else ("independent", None)
    if kind in MAPPED:
        return "structure", MAPPED[kind]
    if kind in UNMAPPED:
        return "structure", None
    # Unknown procedures are the tests' own helpers; like custom reports they
    # compare values the test computed.
    return "independent", None


def ledger(manifest):
    import re
    verified = {c["path"] for c in manifest["cases"]
                if not c.get("derived") and c["expected_rust"] == "pass" and c["expected_occt"] == "pass"}
    statuses, kinds, lost = Counter(), {}, Counter()
    survey = Counter()
    rows = []
    for path in case_files():
        text = path.read_text(errors="replace")
        relative = str(path.relative_to(ROOT))
        survey["cases"] += 1
        survey["load_external_data"] += "locate_data_file" in text
        survey["touch_the_viewer"] += bool(re.search(r"\b(checkview|vinit|vdisplay|vdump|pload[^\n]*VISUALIZATION)\b", text))
        survey["property_assertions"] += bool(re.search(r"\b(checkprops|vprops|sprops|lprops)\b", text))
        survey["subshape_count_assertions"] += bool(re.search(r"\b(checknbshapes|nbshapes)\b", text))
        survey["validity_assertions"] += bool(re.search(r"\bcheckshape\b", text))
        survey["tolerance_maxima_assertions"] += bool(re.search(r"\bcheckmaxtol\b", text))
        picks = set(re.findall(r"(?:^|[\[;{])\s*explode\s+(\S+)", text, re.M))
        for kind, line in assertions(text):
            dependence, mapping = classify(kind, line, picks)
            if dependence == "independent":
                status = "model-independent"
            elif mapping and relative in verified:
                status = "mapped-and-verified"
            else:
                status = "lost"
                lost["mapping unverified: "+mapping if mapping else "no mapping: "+UNMAPPED.get(kind, kind)] += 1
            statuses[status] += 1
            kinds.setdefault(kind, Counter())[status] += 1
            rows.append({"case": relative, "assertion": kind, "status": status, "mapping": mapping})
    summary = {
        "survey": dict(survey), "assertions": sum(statuses.values()), "statuses": dict(sorted(statuses.items())),
        "lost_by_reason": dict(sorted(lost.items())),
        "by_assertion": {k: dict(sorted(v.items())) for k, v in sorted(kinds.items()) if sum(v.values()) >= 20},
    }
    return summary, rows


def junit(results, output):
    suite = ET.Element("testsuite", name="OCCT DRAW compatibility", tests=str(len(results)))
    failures = skipped = 0
    for result in results:
        test = ET.SubElement(suite, "testcase", name=result["case"],
                             classname=result["backend"], time=str(result["seconds"]))
        text = json.dumps(result, indent=2)
        if not result["contract_ok"] or result["status"] in {"failed", "timeout"}:
            failures += 1
            ET.SubElement(test, "failure", message=result["status"]).text = text
        elif result["status"] != "pass":
            skipped += 1
            ET.SubElement(test, "skipped", message=result["status"]).text = text
        ET.SubElement(test, "system-out").text = (
            Path(result["log"]).read_text(errors="replace") if result["log"] else result.get("error", ""))
    suite.set("failures", str(failures))
    suite.set("skipped", str(skipped))
    ET.ElementTree(suite).write(output, encoding="utf-8", xml_declaration=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["rust", "occt", "both"], default="rust")
    parser.add_argument("--case", action="append", help="manifest test path; repeat to select cases")
    parser.add_argument("--draw-exe", default=shutil.which("DRAWEXE") or shutil.which("occt-draw") or "/opt/homebrew/opt/opencascade/bin/DRAWEXE")
    parser.add_argument("--tclsh", default=shutil.which("tclsh") or "tclsh")
    parser.add_argument("--data-dir", action="append", type=Path, default=[])
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--output", type=Path, default=ROOT / "target/upstream-tests")
    parser.add_argument("--ledger", action="store_true",
                        help="only recompute the coverage ledger and check it against the manifest")
    parser.add_argument("--write-ledger", action="store_true",
                        help="record the recomputed ledger in the manifest (review the change)")
    args = parser.parse_args()
    if os.name != "posix":
        parser.error("the DRAW bridge currently requires Unix process-group supervision")
    if not 0 < args.timeout < 3600:
        parser.error("--timeout must be positive and less than one hour")
    manifest = json.loads(MANIFEST.read_text())
    for path, expected in manifest["sources"].items():
        if digest(ROOT / path) != expected:
            parser.error(f"pinned upstream source changed: {path}; review before updating the manifest")
    for path, expected in manifest.get("derived_sources", {}).items():
        if digest(ROOT / path) != expected:
            parser.error(f"derived case changed: {path}; review before updating the manifest")
    summary, rows = ledger(manifest)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "ledger.json").write_text(json.dumps({"summary": summary, "assertions": rows}, indent=1) + "\n")
    if args.write_ledger:
        manifest["ledger"] = summary
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    ledger_ok = manifest.get("ledger") == summary
    if not ledger_ok:
        print("coverage ledger differs from the manifest; review it, then --write-ledger", flush=True)
        print(json.dumps({"recomputed": summary["statuses"], "recorded": manifest.get("ledger", {}).get("statuses")}))
    if args.ledger:
        print(json.dumps(summary["statuses"], sort_keys=True))
        return 0 if ledger_ok else 1
    cases = manifest["cases"]
    if args.case:
        unknown = set(args.case) - {c["path"] for c in cases}
        if unknown:
            parser.error(f"unregistered tests: {sorted(unknown)}")
        cases = [c for c in cases if c["path"] in args.case]
    backends = ["rust", "occt"] if args.backend == "both" else [args.backend]
    if "occt" in backends and not shutil.which(args.draw_exe):
        parser.error("native DRAWEXE unavailable; supply --draw-exe or run --backend rust")
    if "rust" in backends and not shutil.which(args.tclsh):
        parser.error("Tcl interpreter unavailable; supply --tclsh")
    worker = build_worker() if "rust" in backends else None
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for case in cases:
        derived = case.get("derived", False)
        if derived and case["path"] not in manifest.get("derived_sources", {}):
            parser.error(f"unrecorded derived case: {case['path']}")
        sources = source_files(case["path"], case.get("context") if derived else None)
        grid_path = Path(case["context"]) if derived else Path(case["path"]).parent
        group, grid = grid_path.relative_to("tests").parts
        names = (group, grid, Path(case["path"]).name)
        rules = [p for p in [ROOT / "tests/parse.rules", ROOT / "tests" / group / "parse.rules", ROOT / "tests" / group / grid / "parse.rules"] if p.is_file()]
        for source in sources + rules:
            if derived and str(source.relative_to(ROOT)) == case["path"]:
                continue
            if str(source.relative_to(ROOT)) not in manifest["sources"]:
                parser.error(f"unrecorded upstream context: {source}")
        data_dirs = [ROOT / "data"] + args.data_dir
        for backend in backends:
            extra_env, note = {}, None
            if backend == "rust" and case.get("selector"):
                # The native selector: the same test in native DRAW records
                # each explode pick's geometry for the adapter to match.
                if not shutil.which(args.draw_exe):
                    note = "native selector unavailable"
                else:
                    picks = args.output / "selector" / case["path"] / "picks.txt"
                    picks.parent.mkdir(parents=True, exist_ok=True)
                    picks.unlink(missing_ok=True)
                    native = run_case("occt", sources, picks.parent, None, args.draw_exe, args.tclsh,
                                      args.timeout, data_dirs, names, {"RUSTY_DRAW_SELECTOR_OUT": str(picks)})
                    if native["status"] != "pass" or not picks.is_file():
                        note = "native selector run did not pass: " + native["status"]
                    else:
                        extra_env = {"RUSTY_DRAW_SELECTOR": str(picks.resolve())}
            if note and note != "native selector unavailable":
                result = {"backend": backend, "status": "failed", "queries": "0", "seconds": 0.0,
                          "log": "", "error": note}
            elif note:
                result = {"backend": backend, "status": "skipped", "queries": "0", "seconds": 0.0,
                          "log": "", "error": note}
            else:
                result = run_case(backend, sources, args.output / backend / case["path"],
                                  worker, args.draw_exe, args.tclsh, args.timeout, data_dirs, names, extra_env)
            result.update(case=case["path"], expected=case[f"expected_{backend}"], derived=derived)
            # A Rust-only run cannot select natively; it neither passes nor fails.
            result["contract_ok"] = result["status"] in {result["expected"], "skipped" if note else None}
            results.append(result)
            marker = "" if result["contract_ok"] else " [UNEXPECTED]"
            print(f"{backend:4} {result['status']:15} {case['path']}{marker}", flush=True)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--", "rust", "Cargo.toml", "Cargo.lock"], cwd=ROOT, text=True)
    paired_all = [case for case in cases if len(backends) == 2 and all(
        r["status"] == "pass" for r in results if r["case"] == case["path"])]
    # Derived cases are never counted as original upstream passes.
    paired = [c["path"] for c in paired_all if not c.get("derived")]
    paired_derived = [c["path"] for c in paired_all if c.get("derived")]
    report = {
        "upstream_revision": manifest["upstream_revision"], "rust_revision": revision,
        "working_tree_changes": dirty.splitlines(), "manifest_sha256": digest(MANIFEST),
        "host_sha256": digest(HOST), "runner_sha256": digest(Path(__file__)),
        "worker_sha256": digest(worker) if worker else None,
        "counts": {b: dict(Counter(r["status"] for r in results if r["backend"] == b)) for b in backends},
        "cases_passing_original_assertions_on_both_backends": paired,
        "derived_cases_passing_on_both_backends": paired_derived,
        # Mappings a derived case confirms against native DRAW; original
        # assertions count as mapped-and-verified only in their own cases.
        "mappings_confirmed_natively": sorted({m for c in paired_all for m in c.get("verifies", [])}),
        "coverage_ledger": summary,
        "results": results,
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    junit(results, args.output / "junit.xml")
    print(json.dumps(report["counts"], sort_keys=True))
    print(f"Report: {args.output / 'report.json'}")
    return 0 if ledger_ok and all(r["contract_ok"] and r["status"] not in {"failed", "timeout"} for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
