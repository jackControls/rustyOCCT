#!/usr/bin/env python3
"""Compare one shared fixture corpus against the pure Rust kernel and OCCT.

Requires an installed modeling-only OCCT SDK and a C++17 compiler. The normal
Cargo build/tests need neither. All scratch output stays in target/occt-oracle.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "rust/fixtures/prisms.txt"


def run(command, **kwargs):
    return subprocess.run(command, check=True, text=True, capture_output=True, **kwargs)


def observations(text):
    result = {}
    for line in text.splitlines():
        name, *values = line.split()
        if name in result:
            raise ValueError(f"duplicate case: {name}")
        result[name] = list(map(float, values))
    return result


def compare(rust, occt):
    if rust.keys() != occt.keys() or not rust:
        raise ValueError("oracle and Rust case sets differ or are empty")
    largest_fraction = 0.0
    query_count = 0
    for name, actual in rust.items():
        expected = occt[name]
        if len(actual) != len(expected) or len(actual) < 23:
            raise ValueError(f"{name}: malformed result")
        for index, (a, e) in enumerate(zip(actual, expected)):
            if not math.isfinite(a) or not math.isfinite(e):
                raise ValueError(f"{name}: nonfinite result")
            if index >= 20:
                if a != e:
                    raise AssertionError(f"{name}: topology/classification column {index}: Rust {a}, OCCT {e}")
            else:
                absolute = 1e-6 if index >= 11 else 1e-7
                budget = absolute + 2e-9 * abs(e)
                fraction = abs(a - e) / budget
                largest_fraction = max(largest_fraction, fraction)
                if fraction > 1:
                    raise AssertionError(f"{name}: property column {index}: Rust {a}, OCCT {e}, tolerance {budget}")
        query_count += len(actual) - 23
    return {"cases": len(rust), "point_classifications": query_count, "largest_error_fraction_of_tolerance": largest_fraction}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occt-root", type=Path, default=Path(os.environ.get("OCCT_ROOT", "/opt/homebrew/opt/opencascade" if sys.platform == "darwin" else "/usr")))
    parser.add_argument("--write-baseline", action="store_true", help="Replace the checked-in oracle results after a successful comparison")
    args = parser.parse_args()
    prefix = args.occt_root.resolve()
    include = next((p for p in [prefix / "include/opencascade", prefix / "inc", prefix / "include"] if (p / "Standard_Version.hxx").exists()), None)
    lib = next((p for p in [prefix / "lib", prefix / "lib64", prefix / "lib/x86_64-linux-gnu", prefix / "lib/aarch64-linux-gnu"] if any(p.glob("libTKPrim.*"))), None)
    if not include or not lib:
        parser.error("OCCT SDK not found; set --occt-root or OCCT_ROOT")
    build = ROOT / "target/occt-oracle"
    build.mkdir(parents=True, exist_ok=True)
    executable = build / "occt_oracle"
    libraries = ["TKPrim", "TKTopAlgo", "TKBRep", "TKGeomAlgo", "TKGeomBase", "TKG3d", "TKG2d", "TKMath", "TKernel"]
    command = shlex.split(os.environ.get("CXX", "c++")) + ["-std=c++17", "-O2", str(ROOT / "rust/tools/occt_oracle.cpp"), "-I" + str(include), "-L" + str(lib), "-Wl,-rpath," + str(lib), "-o", str(executable)] + ["-l" + name for name in libraries]
    run(command, cwd=ROOT)
    fixtures = FIXTURES.read_text()
    native = run([str(executable)], input=fixtures, cwd=ROOT)
    rust = run(["cargo", "run", "--quiet", "--locked", "--example", "prism_oracle"], input=fixtures, cwd=ROOT)
    report = compare(observations(rust.stdout), observations(native.stdout))
    report["oracle"] = native.stderr.strip()
    report["properties"] = ["volume", "surface_area", "centroid", "bounds", "central_inertia", "topology_counts", "point_classification"]
    report["tolerances"] = {"property_absolute": 1e-7, "inertia_absolute": 1e-6, "relative": 2e-9, "topology_and_classification": "exact"}
    (build / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    (build / "rust.tsv").write_text(rust.stdout)
    (build / "occt.tsv").write_text(native.stdout)
    if args.write_baseline:
        (ROOT / "rust/fixtures/occt-baseline.tsv").write_text(native.stdout)
        (ROOT / "rust/fixtures/occt-baseline.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(error.stderr, file=sys.stderr)
        raise SystemExit(error.returncode) from error
