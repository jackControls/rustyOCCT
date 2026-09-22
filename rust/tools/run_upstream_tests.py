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


def source_files(case):
    """DRAW group/grid begin scripts, the original test, then reverse end scripts."""
    path = ROOT / case
    directories = list(reversed(path.parent.relative_to(ROOT / "tests").parents))
    # Only ancestors below tests/ are groups/grids; tests/begin does not exist.
    directories = [ROOT / "tests" / p for p in directories if str(p) != "."]
    directories.append(path.parent)
    return ([p / "begin" for p in directories if (p / "begin").is_file()]
            + [path]
            + [p / "end" for p in reversed(directories) if (p / "end").is_file()])


def run_case(backend, sources, directory, worker=None, draw_exe=None,
             tclsh="tclsh", timeout=30.0, data_dirs=()):
    directory.mkdir(parents=True, exist_ok=True)
    result_file = directory / "result.txt"
    result_file.unlink(missing_ok=True)  # A crashed rerun must not reuse an old pass.
    log_file = directory / "output.log"
    env = os.environ.copy()
    case_path = next(p for p in sources if p.name not in {"begin", "end"})
    try:
        group, grid, name = case_path.relative_to(ROOT / "tests").parts
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
        ET.SubElement(test, "system-out").text = Path(result["log"]).read_text(errors="replace")
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
    args = parser.parse_args()
    if os.name != "posix":
        parser.error("the DRAW bridge currently requires Unix process-group supervision")
    if not 0 < args.timeout < 3600:
        parser.error("--timeout must be positive and less than one hour")
    manifest = json.loads(MANIFEST.read_text())
    for path, expected in manifest["sources"].items():
        if digest(ROOT / path) != expected:
            parser.error(f"pinned upstream source changed: {path}; review before updating the manifest")
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
        sources = source_files(case["path"])
        group, grid, _ = Path(case["path"]).relative_to("tests").parts
        rules = [p for p in [ROOT / "tests/parse.rules", ROOT / "tests" / group / "parse.rules", ROOT / "tests" / group / grid / "parse.rules"] if p.is_file()]
        for source in sources + rules:
            if str(source.relative_to(ROOT)) not in manifest["sources"]:
                parser.error(f"unrecorded upstream context: {source}")
        for backend in backends:
            result = run_case(backend, sources, args.output / backend / case["path"],
                              worker, args.draw_exe, args.tclsh, args.timeout,
                              [ROOT / "data"] + args.data_dir)
            result.update(case=case["path"], expected=case[f"expected_{backend}"])
            result["contract_ok"] = result["status"] == result["expected"]
            results.append(result)
            marker = "" if result["contract_ok"] else " [UNEXPECTED]"
            print(f"{backend:4} {result['status']:15} {case['path']}{marker}", flush=True)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--", "rust", "Cargo.toml", "Cargo.lock"], cwd=ROOT, text=True)
    paired = [case["path"] for case in cases if len(backends) == 2 and all(
        r["status"] == "pass" for r in results if r["case"] == case["path"])]
    report = {
        "upstream_revision": manifest["upstream_revision"], "rust_revision": revision,
        "working_tree_changes": dirty.splitlines(), "manifest_sha256": digest(MANIFEST),
        "host_sha256": digest(HOST), "runner_sha256": digest(Path(__file__)),
        "worker_sha256": digest(worker) if worker else None,
        "counts": {b: dict(Counter(r["status"] for r in results if r["backend"] == b)) for b in backends},
        "cases_passing_original_assertions_on_both_backends": paired,
        "results": results,
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    junit(results, args.output / "junit.xml")
    print(json.dumps(report["counts"], sort_keys=True))
    print(f"Report: {args.output / 'report.json'}")
    return 0 if all(r["contract_ok"] and r["status"] not in {"failed", "timeout"} for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
