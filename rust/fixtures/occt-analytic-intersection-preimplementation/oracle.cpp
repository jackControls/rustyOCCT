// Source-pinned observations of IntAna_QuadQuadGeo (S7a of REVIEW_NOTES.md):
// each case block on stdin (`case NAME`, two `surface KIND origin normal x
// [radius [half-angle]]` rows with the kernel's stored frame axes, `end`)
// gives a row `NAME TYPE N`, the result type and number of solutions, then
// one row per solution: `P x y z`, `L point direction`, `C centre normal
// radius`, `E centre normal major-direction major minor`, `H centre normal
// transverse-direction major minor` or `B centre normal axis focal`
// (a parabola). Tolerances are Precision::Angular() and
// Precision::Confusion(), as IntPatch passes them. Nothing is computed from
// the kernel's output.
#include <IntAna_QuadQuadGeo.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Elips.hxx>
#include <gp_Hypr.hxx>
#include <gp_Lin.hxx>
#include <gp_Parab.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>

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

gp_Ax3 frame(const Surf& s) {
  return gp_Ax3(gp_Pnt(s.v[0], s.v[1], s.v[2]), gp_Dir(s.v[3], s.v[4], s.v[5]),
                gp_Dir(s.v[6], s.v[7], s.v[8]));
}

int rank(const std::string& kind) {
  return kind == "plane" ? 0 : kind == "cylinder" ? 1 : kind == "cone" ? 2 : 3;
}

void xyz(std::ostream& out, const gp_XYZ& p) { out << " " << p.X() << " " << p.Y() << " " << p.Z(); }

const char* name(IntAna_ResultType t) {
  switch (t) {
    case IntAna_Point: return "point";
    case IntAna_Line: return "line";
    case IntAna_Circle: return "circle";
    case IntAna_PointAndCircle: return "point_and_circle";
    case IntAna_Ellipse: return "ellipse";
    case IntAna_Parabola: return "parabola";
    case IntAna_Hyperbola: return "hyperbola";
    case IntAna_Empty: return "empty";
    case IntAna_Same: return "same";
    case IntAna_NoGeometricSolution: return "no_geometric_solution";
  }
  return "unknown";
}

void run(const std::string& case_name, Surf a, Surf b) {
  if (rank(a.kind) > rank(b.kind)) std::swap(a, b);
  const double ang = Precision::Angular(), tol = Precision::Confusion();
  IntAna_QuadQuadGeo q;
  auto pln = [](const Surf& s) { return gp_Pln(frame(s)); };
  auto cyl = [](const Surf& s) { return gp_Cylinder(frame(s), s.v[9]); };
  auto con = [](const Surf& s) { return gp_Cone(frame(s), s.v[10], s.v[9]); };
  auto sph = [](const Surf& s) { return gp_Sphere(frame(s), s.v[9]); };
  const std::string k = a.kind + "/" + b.kind;
  if (k == "plane/plane") q.Perform(pln(a), pln(b), ang, tol);
  else if (k == "plane/cylinder") q.Perform(pln(a), cyl(b), ang, tol);
  else if (k == "plane/sphere") q.Perform(pln(a), sph(b));
  else if (k == "plane/cone") q.Perform(pln(a), con(b), ang, tol);
  else if (k == "cylinder/cylinder") q.Perform(cyl(a), cyl(b), tol);
  else if (k == "cylinder/sphere") q.Perform(cyl(a), sph(b), tol);
  else if (k == "cylinder/cone") q.Perform(cyl(a), con(b), tol);
  else if (k == "sphere/sphere") q.Perform(sph(a), sph(b), tol);
  else if (k == "cone/sphere") q.Perform(sph(b), con(a), tol);
  else if (k == "cone/cone") q.Perform(con(a), con(b), tol);
  else {
    std::cout << case_name << " unsupported_pair 0\n";
    return;
  }
  if (!q.IsDone()) {
    std::cout << case_name << " not_done 0\n";
    return;
  }
  const IntAna_ResultType t = q.TypeInter();
  const int n = (t == IntAna_Empty || t == IntAna_Same || t == IntAna_NoGeometricSolution) ? 0 : q.NbSolutions();
  std::ostringstream out;
  out << std::setprecision(17) << case_name << " " << name(t) << " " << n << "\n";
  for (int i = 1; i <= n; ++i) {
    switch (t) {
      case IntAna_Point: out << "P"; xyz(out, q.Point(i).XYZ()); break;
      case IntAna_Line: {
        const gp_Lin l = q.Line(i);
        out << "L"; xyz(out, l.Location().XYZ()); xyz(out, l.Direction().XYZ());
        break;
      }
      case IntAna_Circle: {
        const gp_Circ c = q.Circle(i);
        out << "C"; xyz(out, c.Location().XYZ()); xyz(out, c.Axis().Direction().XYZ());
        out << " " << c.Radius();
        break;
      }
      case IntAna_PointAndCircle:
        if (i == 1) {
          out << "P"; xyz(out, q.Point(1).XYZ());
        } else {
          const gp_Circ c = q.Circle(1);
          out << "C"; xyz(out, c.Location().XYZ()); xyz(out, c.Axis().Direction().XYZ());
          out << " " << c.Radius();
        }
        break;
      case IntAna_Ellipse: {
        const gp_Elips e = q.Ellipse(i);
        out << "E"; xyz(out, e.Location().XYZ()); xyz(out, e.Axis().Direction().XYZ());
        xyz(out, e.XAxis().Direction().XYZ());
        out << " " << e.MajorRadius() << " " << e.MinorRadius();
        break;
      }
      case IntAna_Hyperbola: {
        const gp_Hypr h = q.Hyperbola(i);
        out << "H"; xyz(out, h.Location().XYZ()); xyz(out, h.Axis().Direction().XYZ());
        xyz(out, h.XAxis().Direction().XYZ());
        out << " " << h.MajorRadius() << " " << h.MinorRadius();
        break;
      }
      case IntAna_Parabola: {
        const gp_Parab p = q.Parabola(i);
        out << "B"; xyz(out, p.Location().XYZ()); xyz(out, p.Axis().Direction().XYZ());
        xyz(out, p.XAxis().Direction().XYZ());
        out << " " << p.Focal();
        break;
      }
      default: break;
    }
    out << "\n";
  }
  std::cout << out.str();
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " IntAna_QuadQuadGeo" << std::endl;
  std::string line, case_name;
  std::vector<Surf> surfaces;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string word;
    if (!(in >> word)) continue;
    if (word == "case") {
      in >> case_name;
      surfaces.clear();
    } else if (word == "surface") {
      Surf s;
      in >> s.kind;
      for (double& x : s.v) in >> x;
      surfaces.push_back(s);
    } else if (word == "end") {
      try {
        run(case_name, surfaces.at(0), surfaces.at(1));
      } catch (const Standard_Failure& e) {
        std::cout << case_name << " failure " << (e.GetMessageString() ? e.GetMessageString() : "") << "\n";
      }
      std::cout << std::flush;
    }
  }
  return 0;
}
