// Source-pinned observations of GeomAPI_IntCS (S7c of REVIEW_NOTES.md): each
// case block on stdin (`case NAME`, a `curve line x y z dx dy dz` row, the
// line through a point along a direction, or a `curve circle origin normal x
// radius` row with the kernel's stored frame axes (S7c.2: `curve ellipse` and
// `curve hyperbola origin normal x major minor`, and `curve spline degree n
// x y z w ... k knot mult ...`, a clamped rational B-spline), a `surface KIND origin
// normal x [radius [half-angle]]` row (planes, cylinders, cones, spheres and
// tori as `major minor`), `end`) gives `NAME done|not_done POINTS SEGMENTS`,
// then per point `P x y z w` (w the curve's parameter: a line's distance
// along its unit direction from the point, a circle's angle) and per segment
// `S w1 w2`. Nothing is computed from the kernel's output.
#include <GeomAPI_IntCS.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
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
  std::vector<double> v = std::vector<double>(11, 0.0);
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
  if (c.kind == "spline") {
    const int degree = int(c.v[0]), n = int(c.v[1]);
    TColgp_Array1OfPnt poles(1, n);
    TColStd_Array1OfReal weights(1, n);
    for (int i = 0; i < n; ++i) {
      const double* q = &c.v[2 + 4 * i];
      poles.SetValue(i + 1, gp_Pnt(q[0], q[1], q[2]));
      weights.SetValue(i + 1, q[3]);
    }
    const int at = 2 + 4 * n, k = int(c.v[at]);
    TColStd_Array1OfReal knots(1, k);
    TColStd_Array1OfInteger mults(1, k);
    for (int i = 0; i < k; ++i) {
      knots.SetValue(i + 1, c.v[at + 1 + 2 * i]);
      mults.SetValue(i + 1, int(c.v[at + 2 + 2 * i]));
    }
    return new Geom_BSplineCurve(poles, weights, knots, mults, degree);
  }
  const gp_Ax2 ax(gp_Pnt(c.v[0], c.v[1], c.v[2]), gp_Dir(c.v[3], c.v[4], c.v[5]),
                  gp_Dir(c.v[6], c.v[7], c.v[8]));
  if (c.kind == "ellipse") return new Geom_Ellipse(ax, c.v[9], c.v[10]);
  if (c.kind == "hyperbola") return new Geom_Hyperbola(ax, c.v[9], c.v[10]);
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
      r.v.clear();
      double x;
      while (in >> x) r.v.push_back(x);
      if (r.v.size() < 11) r.v.resize(11, 0.0);
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
