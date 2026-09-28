// Source-pinned observations of GeomAPI_IntCS (S7c of REVIEW_NOTES.md): each
// case block on stdin (`case NAME`, a `curve line x y z dx dy dz` row, the
// line through a point along a direction, or a `curve circle origin normal x
// radius` row with the kernel's stored frame axes, a `surface KIND origin
// normal x [radius [half-angle]]` row (planes, cylinders, cones, spheres and
// tori as `major minor`), `end`) gives `NAME done|not_done POINTS SEGMENTS`,
// then per point `P x y z w` (w the curve's parameter: a line's distance
// along its unit direction from the point, a circle's angle) and per segment
// `S w1 w2`. Nothing is computed from the kernel's output.
#include <GeomAPI_IntCS.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Row {
  std::string kind;
  double v[11] = {0};
};

Handle(Geom_Surface) surface(const Row& s) {
  const gp_Ax3 ax(gp_Pnt(s.v[0], s.v[1], s.v[2]), gp_Dir(s.v[3], s.v[4], s.v[5]),
                  gp_Dir(s.v[6], s.v[7], s.v[8]));
  if (s.kind == "plane") return new Geom_Plane(ax);
  if (s.kind == "cylinder") return new Geom_CylindricalSurface(ax, s.v[9]);
  if (s.kind == "cone") return new Geom_ConicalSurface(ax, s.v[10], s.v[9]);
  if (s.kind == "torus") return new Geom_ToroidalSurface(ax, s.v[9], s.v[10]);
  return new Geom_SphericalSurface(ax, s.v[9]);
}

Handle(Geom_Curve) curve(const Row& c) {
  if (c.kind == "line")
    return new Geom_Line(gp_Pnt(c.v[0], c.v[1], c.v[2]), gp_Dir(c.v[3], c.v[4], c.v[5]));
  const gp_Ax2 ax(gp_Pnt(c.v[0], c.v[1], c.v[2]), gp_Dir(c.v[3], c.v[4], c.v[5]),
                  gp_Dir(c.v[6], c.v[7], c.v[8]));
  return new Geom_Circle(ax, c.v[9]);
}

void run(const std::string& name, const Row& c, const Row& s) {
  GeomAPI_IntCS inter(curve(c), surface(s));
  std::ostringstream out;
  out << std::setprecision(17);
  if (!inter.IsDone()) {
    out << name << " not_done 0 0\n";
    std::cout << out.str();
    return;
  }
  out << name << " done " << inter.NbPoints() << " " << inter.NbSegments() << "\n";
  for (int i = 1; i <= inter.NbPoints(); ++i) {
    const gp_Pnt p = inter.Point(i);
    double u, v, w;
    inter.Parameters(i, u, v, w);
    out << "P " << p.X() << " " << p.Y() << " " << p.Z() << " " << w << "\n";
  }
  for (int i = 1; i <= inter.NbSegments(); ++i) {
    const Handle(Geom_Curve) seg = inter.Segment(i);
    out << "S " << seg->FirstParameter() << " " << seg->LastParameter() << "\n";
  }
  std::cout << out.str();
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " GeomAPI_IntCS" << std::endl;
  std::string line, name;
  Row c, s;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string word;
    if (!(in >> word)) continue;
    if (word == "case") {
      in >> name;
      c = Row();
      s = Row();
    } else if (word == "curve" || word == "surface") {
      Row& r = word == "curve" ? c : s;
      in >> r.kind;
      for (double& x : r.v) in >> x;
    } else if (word == "end") {
      try {
        run(name, c, s);
      } catch (const Standard_Failure& e) {
        std::cout << name << " failure 0 0\n";
      }
      std::cout << std::flush;
    }
  }
  return 0;
}
