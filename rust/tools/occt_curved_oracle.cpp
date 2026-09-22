// Test-only OCCT observations, independent of Rust results.
#include <IntAna_IntConicQuad.hxx>
#include <IntAna_Quadric.hxx>
#include <IntAna2d_AnaIntersection.hxx>
#include <IntAna2d_IntPoint.hxx>
#include <Standard_Version.hxx>
#include <gp_Circ2d.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Lin.hxx>
#include <gp_Lin2d.hxx>
#include <gp_Sphere.hxx>
#include <gp_Vec.hxx>
#include <math_DirectPolynomialRoots.hxx>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

gp_Pnt read_point() {
  double x,y,z;
  if (!(std::cin >> x >> y >> z)) throw std::runtime_error("incomplete point");
  return gp_Pnt(x,y,z);
}

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::cout << std::setprecision(17);
  std::string kind,name;
  while (std::cin >> kind >> name) {
    std::cout << name;
    if (kind == "q") {
      double a,b,c;
      if (!(std::cin >> a >> b >> c)) return 2;
      const math_DirectPolynomialRoots roots(a,b,c);
      if (!roots.IsDone()) return 3;
      if (roots.InfiniteRoots()) { std::cout << " A\n"; continue; }
      std::vector<double> values;
      for (int i=1;i<=roots.NbSolutions();++i) values.push_back(roots.Value(i));
      std::sort(values.begin(),values.end());
      std::cout << ' ' << values.size();
      for (double value: values) std::cout << ' ' << value;
      std::cout << '\n';
      continue;
    }
    if (kind != "s" && kind != "y" && kind != "c") return 2;
    const gp_Pnt center = read_point();
    const gp_Pnt axis = read_point();
    double radius;
    if (!(std::cin >> radius)) return 2;
    const gp_Pnt p = read_point();
    const gp_Pnt q = read_point();
    const gp_Vec direction(p,q);
    std::vector<std::array<double,4>> points;
    if (kind == "c") {
      // Native circle comparison is explicitly the common XY-plane domain.
      const gp_Circ2d circle(gp_Ax2d(gp_Pnt2d(center.X(),center.Y()),gp_Dir2d(1,0)),radius);
      const gp_Lin2d line(gp_Pnt2d(p.X(),p.Y()),gp_Dir2d(q.X()-p.X(),q.Y()-p.Y()));
      const IntAna2d_AnaIntersection intersection(line,circle);
      if (!intersection.IsDone()) return 3;
      for (int i=1;i<=intersection.NbPoints();++i) {
        const auto& hit=intersection.Point(i);
        points.push_back({hit.ParamOnFirst()/direction.Magnitude(),hit.Value().X(),hit.Value().Y(),center.Z()});
      }
    } else {
      const gp_Ax3 position(center,gp_Dir(axis.X(),axis.Y(),axis.Z()));
      const IntAna_Quadric quadric = kind == "s"
        ? IntAna_Quadric(gp_Sphere(position,radius))
        : IntAna_Quadric(gp_Cylinder(position,radius));
      const IntAna_IntConicQuad intersection(gp_Lin(p,gp_Dir(direction)),quadric);
      if (!intersection.IsDone()) return 3;
      if (intersection.IsInQuadric()) { std::cout << " A\n"; continue; }
      for (int i=1;i<=intersection.NbPoints();++i) {
        const gp_Pnt hit = intersection.Point(i);
        points.push_back({intersection.ParamOnConic(i)/direction.Magnitude(),hit.X(),hit.Y(),hit.Z()});
      }
    }
    std::sort(points.begin(),points.end());
    // OCCT's polynomial intersector may emit a double root twice. Compare
    // geometric hit count after merging exactly equal native parameters only.
    points.erase(std::unique(points.begin(),points.end(),[](const auto& a,const auto& b) {return a[0]==b[0];}),points.end());
    std::cout << ' ' << points.size();
    for (const auto& point: points) for (double value:point) std::cout << ' ' << value;
    std::cout << '\n';
  }
}
