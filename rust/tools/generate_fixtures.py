#!/usr/bin/env python3
"""Regenerate deterministic prism inputs. Never generates expected results."""
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[2]


def polygon(*points):
    return ("P", list(points))


def circle(x, y, radius):
    return ("C", (x, y, radius))


TEMPLATES = [
    ("box", [polygon((0,0),(20,0),(20,12),(0,12))]),
    ("triangle", [polygon((0,0),(10,0),(2,8))]),
    ("concave", [polygon((0,0),(20,0),(20,6),(6,6),(6,20),(0,20))]),
    ("frame", [polygon((0,0),(20,0),(20,12),(0,12)), polygon((3,2),(8,2),(8,8),(3,8))]),
    ("cylinder", [circle(2,-3,10)]),
    ("tube", [circle(-2,3,10),circle(-2,3,6)]),
    ("plate", [polygon((0,0),(80,0),(80,40),(0,40))] + [circle(x,y,3) for x,y in [(10,10),(70,10),(70,30),(10,30)]]),
    ("mixed", [polygon((0,0),(40,0),(40,20),(0,20)), circle(10,10,2), polygon((28,8),(32,8),(32,12),(28,12))]),
    ("round_square_hole", [circle(0,0,10),polygon((-2,-2),(2,-2),(2,2),(-2,2))]),
    ("collinear", [polygon((0,0),(2,0),(4,0),(4,3),(0,3))]),
    ("concave_hole", [polygon((0,0),(20,0),(20,6),(6,6),(6,20),(0,20)),circle(3,12,1)]),
]


def main():
    rng = random.Random(20260922)
    cases = []
    for template_name, boundaries in TEMPLATES:
        for variant in range(6):
            scale = [1.0, 0.02, 3.5, 1.0, 0.5, 1.0][variant]
            tx, ty = (1e4, -1e4) if variant == 5 else (0.0, 0.0)
            origin = [rng.uniform(-100,100) for _ in range(3)] if variant in (2,3) else [0,0,0]
            normal = [1,2,3] if variant in (2,3,5) else ([0,1,0] if variant == 4 else [0,0,1])
            hint = [1,0,0]
            start, end = (-3.0*scale, 7.0*scale)
            if variant % 2: start, end = end, start
            lines = [f"{template_name}_{variant} {len(boundaries)}", " ".join(map(str,origin + normal + hint + [start,end]))]
            boundary_samples = []
            extent = []
            for index, (kind, geometry) in enumerate(boundaries):
                if kind == "P":
                    points = [(x*scale+tx,y*scale+ty) for x,y in geometry]
                    if variant % 2: points = points[:1] + points[:0:-1]
                    lines.append("P " + str(len(points)) + " " + " ".join(str(value) for point in points for value in point))
                    boundary_samples.append(points[0])
                    if index == 0: extent = points
                else:
                    x,y,r = geometry
                    x,y,r = x*scale+tx,y*scale+ty,r*scale
                    lines.append(f"C {x} {y} {r}")
                    boundary_samples.append((x+r,y))
                    if index == 0: extent = [(x-r,y-r),(x+r,y+r)]
            xmin,xmax = min(p[0] for p in extent),max(p[0] for p in extent)
            ymin,ymax = min(p[1] for p in extent),max(p[1] for p in extent)
            low,high = min(start,end),max(start,end)
            mid = (low+high)/2
            center = ((xmin+xmax)/2,(ymin+ymax)/2)
            queries = [(x,y,z) for x,y in boundary_samples for z in (low,mid,high)]
            queries += [(*center,z) for z in (low,mid,high,low-scale,high+scale)]
            queries += [(rng.uniform(xmin-scale,xmax+scale),rng.uniform(ymin-scale,ymax+scale),rng.uniform(low-scale,high+scale)) for _ in range(24)]
            lines.append(str(len(queries)))
            lines.extend(" ".join(map(str,query)) for query in queries)
            cases.append("\n".join(lines))
    destination = ROOT / "rust/fixtures/prisms.txt"
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(str(len(cases)) + "\n" + "\n".join(cases) + "\n")
    print(f"Wrote {len(cases)} deterministic cases to {destination}")


if __name__ == "__main__":
    main()
