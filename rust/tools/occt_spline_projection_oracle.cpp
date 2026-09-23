// Research-only native observations, captured before any Rust projection API.
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string row;
  while (std::getline(std::cin, row)) {
    std::istringstream in(row);
    std::string name, extra;
    int degree, periodic, np, nk;
    double first, last, qx, qy, qz;
    if (!(in >> name >> degree >> periodic >> np >> nk >> first >> last >> qx >> qy >> qz)) return 2;
    try {
      TColgp_Array1OfPnt poles(1,np);
      TColStd_Array1OfReal weights(1,np), knots(1,nk);
      TColStd_Array1OfInteger mults(1,nk);
      for (int i=1; i<=np; ++i) {
        double x,y,z,w;
        if (!(in >> x >> y >> z >> w)) return 2;
        poles(i)=gp_Pnt(x,y,z); weights(i)=w;
      }
      for (int i=1; i<=nk; ++i) if (!(in >> knots(i) >> mults(i))) return 2;
      if (in >> extra) return 2;
      Handle(Geom_BSplineCurve) curve=new Geom_BSplineCurve(poles,weights,knots,mults,degree,periodic!=0,false);
      GeomAPI_ProjectPointOnCurve projector(gp_Pnt(qx,qy,qz),curve,first,last);
      const auto& extrema=projector.Extrema();
      int count=projector.NbPoints();
      std::ostringstream out;
      out << std::setprecision(17) << name << " R " << extrema.IsDone() << ' ' << count;
      for (int i=1; i<=count; ++i) {
        const auto point=projector.Point(i);
        out << ' ' << projector.Parameter(i) << ' ' << extrema.IsMin(i) << ' ' << extrema.SquareDistance(i)
            << ' ' << point.X() << ' ' << point.Y() << ' ' << point.Z();
      }
      double a,b; gp_Pnt pa,pb;
      extrema.TrimmedSquareDistances(a,b,pa,pb);
      out << " E " << a << ' ' << b << ' ' << pa.X() << ' ' << pa.Y() << ' ' << pa.Z()
          << ' ' << pb.X() << ' ' << pb.Y() << ' ' << pb.Z();
      if (count) out << " M " << projector.LowerDistanceParameter() << ' ' << projector.LowerDistance();
      std::cout << out.str() << '\n';
    } catch (const Standard_Failure& failure) {
      std::cout << name << " X\n";
      std::cerr << name << ": " << failure.GetMessageString() << '\n';
    }
  }
}
