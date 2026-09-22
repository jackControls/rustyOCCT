#!/usr/bin/env python3
"""Independent exact oracle: Python rational arithmetic on represented f64 inputs.

Expected signs never come from Rust or OCCT. Store the input bits, not rounded
decimal text. --check verifies reproducibility without changing the baseline.
"""
import argparse
from fractions import Fraction
import math
from pathlib import Path
import random
import struct

OUTPUT = Path(__file__).resolve().parents[1] / "fixtures/orient2d.tsv"


def bits(value):
    return struct.unpack(">Q", struct.pack(">d", value))[0]


def value(word):
    return struct.unpack(">d", struct.pack(">Q", word))[0]


def exact_sign(coordinates):
    ax, ay, bx, by, cx, cy = map(Fraction.from_float, coordinates)
    determinant = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    return (determinant > 0) - (determinant < 0)


def generate():
    rows = []

    def add(label, coordinates):
        if all(math.isfinite(x) for x in coordinates):
            rows.append(f"{label}\t" + "\t".join(f"{bits(x):016x}" for x in coordinates)
                        + f"\t{exact_sign(coordinates)}")

    tiny = value(1)
    maximum = value(0x7fefffffffffffff)
    for label, coordinates in [
        ("zero", [0., 0., 0., 0., 0., 0.]),
        ("signed_zero", [-0., 0., 0., -0., -0., -0.]),
        ("duplicate", [1., 2., 1., 2., 3., 4.]),
        ("unit", [0., 0., 1., 0., 0., 1.]),
        ("subnormal_product", [0., 0., tiny, 0., 0., tiny]),
        ("subnormal_difference", [tiny, tiny, 2*tiny, tiny, tiny, 2*tiny]),
        ("overflow_product", [0., 0., maximum, 0., 0., maximum]),
        ("overflow_difference", [-maximum, -maximum, maximum, -maximum, 0., maximum]),
        ("full_exponent_span", [tiny, maximum, maximum, tiny, -maximum, -tiny]),
        ("maximum_collinear", [-maximum, -maximum, maximum, maximum, tiny, tiny]),
        ("maximum_adjacent", [maximum, maximum, math.nextafter(maximum, 0.), maximum, maximum, math.nextafter(maximum, 0.)]),
    ]:
        add(label, coordinates)

    for power in [-1074, -1073, -1023, -1022, -1000, -800, -500, -401, -400,
                  -399, -100, -1, 0, 1, 399, 400, 401, 500, 800, 1000, 1023]:
        x = math.ldexp(1., power)
        for direction in [-math.inf, math.inf]:
            add(f"exponent_boundary_{power}_{direction}",
                [x, x, math.nextafter(x, direction), x, x, math.nextafter(x, direction)])

    for magnitude in [26, 27, 40, 52]:
        u = float(2**magnitude)
        for power in range(-1050, 901, 25):
            coordinates = [math.ldexp(x, power) for x in [0., 0., u, u-1., u+1., u]]
            add(f"cancel_{magnitude}_{power}", coordinates)

    rng = random.Random(0x0CC7_2D)
    for i in range(1024):
        coordinates = []
        while len(coordinates) < 6:
            x = value(rng.getrandbits(64))
            if math.isfinite(x):
                coordinates.append(x)
        add(f"finite_bits_{i}", coordinates)

    for i in range(1024):
        ax, ay = [rng.randint(-2**30, 2**30) for _ in range(2)]
        dx, dy = [rng.randint(-2**20, 2**20) for _ in range(2)]
        t = rng.randint(-8, 8)
        exponent = rng.randint(-500, 500)
        coordinates = [math.ldexp(float(x), exponent)
                       for x in [ax, ay, ax+dx, ay+dy, ax+t*dx, ay+t*dy]]
        if i % 4:
            axis = 4 + i % 2
            coordinates[axis] = math.nextafter(coordinates[axis], math.inf if i % 3 else -math.inf)
        add(f"near_collinear_{i}", coordinates)

    header = ("# Exact orientation oracle from Python fractions.Fraction.from_float.\n"
              "# label ax_bits ay_bits bx_bits by_bits cx_bits cy_bits sign\n"
              f"# cases: {len(rows)}; generator: rust/tools/generate_predicate_fixtures.py\n")
    return header + "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = generate()
    if args.check:
        if OUTPUT.read_text() != expected:
            parser.error("exact predicate fixtures differ; investigate before updating the baseline")
        print("Exact rational predicate fixtures are reproducible")
    else:
        OUTPUT.write_text(expected)
        print(expected.splitlines()[2])


if __name__ == "__main__":
    main()
