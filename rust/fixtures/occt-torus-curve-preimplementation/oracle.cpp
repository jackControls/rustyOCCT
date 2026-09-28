// Source-pinned observations of GeomInt_IntSS, the engine under GeomAPI_IntSS
// that also reports isolated points (S7b of REVIEW_NOTES.md):
// each case block on stdin (`case NAME`, two `surface KIND origin normal x
// [radius [half-angle]]` rows with the kernel's stored frame axes, planes,
// cylinders, cones, spheres and tori (`major minor`), `end`) gives `NAME done|not_done LINES POINTS`, then per intersection line
// `L closed|open first last` and 17 rows `p x y z` at equally spaced
// parameters, and per isolated point `P x y z`. Tolerance 1e-7. Nothing is
// computed from the kernel's output.
#include <GeomInt_IntSS.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <gp_Ax3.hxx>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Surf {
  std::string kind;
  double v[11] = {0};
};

Handle(Geom_Surface) make(const Surf& s) {
  const gp_Ax3 ax(gp_Pnt(s.v[0], s.v[1], s.v[2]), gp_Dir(s.v[3], s.v[4], s.v[5]),
                  gp_Dir(s.v[6], s.v[7], s.v[8]));
  if (s.kind == "cylinder") return new Geom_CylindricalSurface(ax, s.v[9]);
  if (s.kind == "cone") return new Geom_ConicalSurface(ax, s.v[10], s.v[9]);
  if (s.kind == "plane") return new Geom_Plane(ax);
  if (s.kind == "torus") return new Geom_ToroidalSurface(ax, s.v[9], s.v[10]);
  return new Geom_SphericalSurface(ax, s.v[9]);
}

void run(const std::string& name, const Surf& a, const Surf& b) {
  GeomInt_IntSS inter(make(a), make(b), 1.0e-7, true, false, false);
  std::ostringstream out;
  out << std::setprecision(17);
  if (!inter.IsDone()) {
    out << name << " not_done 0 0\n";
    std::cout << out.str();
    return;
  }
  out << name << " done " << inter.NbLines() << " " << inter.NbPoints() << "\n";
  for (int i = 1; i <= inter.NbLines(); ++i) {
    const Handle(Geom_Curve) c = inter.Line(i);
    const double f = c->FirstParameter(), l = c->LastParameter();
    out << "L " << (c->IsClosed() ? "closed" : "open") << " " << f << " " << l << "\n";
    for (int k = 0; k <= 16; ++k) {
      const gp_Pnt p = c->Value(f + (l - f) * k / 16.0);
      out << "p " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
    }
  }
  for (int i = 1; i <= inter.NbPoints(); ++i) {
    const gp_Pnt p = inter.Point(i);
    out << "P " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
  }
  std::cout << out.str();
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " GeomInt_IntSS" << std::endl;
  std::string line, name;
  std::vector<Surf> surfaces;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string word;
    if (!(in >> word)) continue;
    if (word == "case") {
      in >> name;
      surfaces.clear();
    } else if (word == "surface") {
      Surf s;
      in >> s.kind;
      for (double& x : s.v) in >> x;
      surfaces.push_back(s);
    } else if (word == "end") {
      try {
        run(name, surfaces.at(0), surfaces.at(1));
      } catch (const Standard_Failure& e) {
        std::cout << name << " failure 0 0\n";
      }
      std::cout << std::flush;
    }
  }
  return 0;
}
