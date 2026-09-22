// Test-only reference. Never linked into the Rust kernel.
#include <IntAna_IntConicQuad.hxx>
#include <Standard_Version.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Vec.hxx>
#include <iomanip>
#include <iostream>
#include <string>

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::cout << std::setprecision(17);
  std::string name;
  while (std::cin >> name) {
    gp_Pnt points[5];
    for (auto &point : points) {
      double x, y, z;
      if (!(std::cin >> x >> y >> z)) return 2;
      point.SetCoord(x,y,z);
    }
    const gp_Vec direction(points[3],points[4]);
    const gp_Vec normal = gp_Vec(points[0],points[1]).Crossed(gp_Vec(points[0],points[2]));
    const gp_Pln plane(points[0],gp_Dir(normal));
    const gp_Lin line(points[3],gp_Dir(direction));
    const IntAna_IntConicQuad intersection(line,plane,1e-12,1e-12);
    if (!intersection.IsDone()) return 3;
    std::cout << name;
    if (intersection.IsInQuadric()) std::cout << " C\n";
    else if (intersection.IsParallel()) std::cout << " N\n";
    else {
      if (intersection.NbPoints() != 1) return 4;
      const gp_Pnt p = intersection.Point(1);
      std::cout << " P " << intersection.ParamOnConic(1)/direction.Magnitude()
                << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z() << '\n';
    }
  }
}
