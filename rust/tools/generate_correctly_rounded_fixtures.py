#!/usr/bin/env python3
"""Independent oracle for the kernel's correctly rounded elementary functions
(`rust/kernel/src/cr.rs`, R13 of the review of 18 October 2026).

Each value is mpmath's at 2400 bits (3600 for arguments past 2^60), rounded
once to the nearest binary64, ties to even. The arguments are seeded
samples per binade, the cases X1 found (macOS's `hypot` and `tan` one ulp
off), arguments next to multiples of pi/2 and to 1, huge arguments and
exact cases. Inputs and outputs are stored as bits. --check verifies
reproducibility without changing the fixture.
"""
import argparse
from pathlib import Path
import random
import struct

import mpmath as mp

OUTPUT = Path(__file__).resolve().parents[1] / "fixtures/correctly-rounded.tsv"


def bits(x):
    return struct.unpack(">Q", struct.pack(">d", x))[0]


def value(word):
    return struct.unpack(">d", struct.pack(">Q", word))[0]


def rn(f, *args):
    """f(args) rounded once to the nearest binary64 (mpmath's float()
    rounds to nearest, ties to even)."""
    huge = max(abs(a) for a in args) > 2.0**60
    with mp.workprec(3600 if huge else 2400):
        return float(f(*[mp.mpf(a) for a in args]))


FUNCTIONS = {
    "sin": mp.sin,
    "cos": mp.cos,
    "tan": mp.tan,
    "atan": mp.atan,
    "atan2": mp.atan2,
    "asin": mp.asin,
    "acos": mp.acos,
    "hypot": lambda x, y: mp.sqrt(x * x + y * y),
    "exp": mp.exp,
    "ln": mp.log,
    "log2": lambda x: mp.log(x, 2),
    "sinh": mp.sinh,
    "cosh": mp.cosh,
    "asinh": mp.asinh,
    "atanh": mp.atanh,
}


def sample(rng, lo, hi, signed=True):
    """A binary64 with exponent in [lo, hi), random significand and sign."""
    e = rng.randrange(lo, hi)
    m = 1.0 + rng.getrandbits(52) * 2.0**-52
    x = m * 2.0**e
    return -x if signed and rng.getrandbits(1) else x


def up(x, n=1):
    return value(bits(x) + n) if x > 0 else value(bits(x) - n)


def arguments():
    rng = random.Random(20261018)
    cases = []

    def add(name, label, *args):
        cases.append((name, label, args))

    # X1: the two quantities that differed between macOS and glibc.
    add("hypot", "x1_octant_tilt_sphere_bounds", value(0x3fdc71c71c71c71d), value(0xbfbc71c71c71c71e))
    add("tan", "x1_consistency_ribbed_cone_slope", value(0xbfe921fb54442d16))
    n = 120
    for _ in range(n):
        add("sin", "sample", sample(rng, -26, 22))
        add("cos", "sample", sample(rng, -26, 22))
        add("tan", "sample", sample(rng, -26, 22))
        add("atan", "sample", sample(rng, -27, 40))
        add("atan2", "sample", sample(rng, -30, 30), sample(rng, -30, 30))
        add("asin", "sample", sample(rng, -26, 0))
        add("acos", "sample", sample(rng, -26, 0))
        add("hypot", "sample", sample(rng, -40, 40), sample(rng, -40, 40))
        add("exp", "sample", sample(rng, -50, 9))
        add("ln", "sample", sample(rng, -1000, 1000, False))
        add("log2", "sample", sample(rng, -1000, 1000, False))
        add("sinh", "sample", sample(rng, -27, 5))
        add("cosh", "sample", sample(rng, -27, 5))
        add("asinh", "sample", sample(rng, -27, 40))
        add("atanh", "sample", sample(rng, -27, 0))
    # Next to multiples of pi/2: the sine's and cosine's smallest results.
    for k in (1, 2, 3, 4, 7, 100, 355, 1_000_003, 2**30):
        v = float(mp.pi / 2 * k)
        for d in (-2, -1, 0, 1, 2):
            x = up(v, d) if d >= 0 else value(bits(v) + d)
            add("sin", "near_half_pi", x)
            add("cos", "near_half_pi", x)
            add("tan", "near_half_pi", x)
    # Huge arguments (the second reduction).
    for x in (1e22, -1e22, 1e300, 2.0**1023, 5.319372648326541e255):
        add("sin", "huge", x)
        add("cos", "huge", x)
    # Next to 1 for the inverse sine and cosine, and atanh.
    for d in range(1, 6):
        x = value(bits(1.0) - d)
        add("asin", "near_one", x)
        add("acos", "near_one", x)
        add("acos", "near_one", -x)
        add("atanh", "near_one", x)
    # Exact cases and ties.
    add("hypot", "exact", 3.0, 4.0)
    add("hypot", "exact", -5.0, 12.0)
    add("hypot", "exact", 3.0 * 2.0**1000, 4.0 * 2.0**1000)
    add("log2", "exact", 1024.0)
    add("log2", "exact", 2.0**-1074)
    add("atan2", "exact", 1.0, 1.0)
    add("atan2", "exact", -1.0, -1.0)
    return cases


def generate():
    rows = []
    for name, label, args in arguments():
        y = rn(FUNCTIONS[name], *args)
        rows.append("\t".join([name, label] + [f"{bits(a):016x}" for a in args] + [f"{bits(y):016x}"]))
    header = ("# Correctly rounded elementary functions: mpmath at 2400 bits, rounded once.\n"
              "# function label argument_bits... result_bits\n"
              f"# cases: {len(rows)}; generator: rust/tools/generate_correctly_rounded_fixtures.py\n")
    return header + "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = generate()
    if args.check:
        if OUTPUT.read_text() != expected:
            parser.error("correctly rounded fixtures differ; investigate before updating them")
        print("Correctly rounded fixtures are reproducible")
    else:
        OUTPUT.write_text(expected)
        print(expected.splitlines()[2])


if __name__ == "__main__":
    main()
