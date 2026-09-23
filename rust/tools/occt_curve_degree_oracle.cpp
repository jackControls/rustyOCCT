// Test-only complete degree-elevation observations; never linked into the Rust kernel.
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
#include <vector>

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string row;
  while (std::getline(std::cin, row)) {
    std::istringstream in(row); std::string name, kind, extra;
    int degree, np, nk, operations;
    if (!(in >> name >> kind >> degree >> np >> nk >> operations)) return 2;
    try {
      TColgp_Array1OfPnt poles(1, np);
      TColStd_Array1OfReal weights(1, np), knots(1, nk);
      TColStd_Array1OfInteger mults(1, nk);
      for (int i=1; i<=np; ++i) {
        double x, y, z, w; if (!(in >> x >> y >> z >> w)) return 2;
        poles(i)=gp_Pnt(x,y,z); weights(i)=w;
      }
      for (int i=1; i<=nk; ++i) if (!(in >> knots(i) >> mults(i))) return 2;
      Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(
        poles, weights, knots, mults, degree, kind=="P", false);
      std::vector<bool> results;
      for (int i=0; i<operations; ++i) {
        char op; double u; int target;
        if (!(in >> op >> u >> target)) return 2;
        if (op=='D') { curve->IncreaseDegree(target); results.push_back(true); }
        else if (op=='I') { curve->InsertKnot(u, target, 0., false); results.push_back(true); }
        else if (op=='R') {
          int index=0;
          for (int k=1; k<=curve->NbKnots(); ++k) if (curve->Knot(k)==u) index=k;
          if (!index) return 3;
          results.push_back(curve->RemoveKnot(index, target, 1.e-9));
        } else return 2;
      }
      if (in >> extra) return 2;
      std::ostringstream out;
      out << std::setprecision(17) << name << " R " << results.size();
      for (bool result: results) out << ' ' << result;
      out << ' ' << curve->Degree() << ' ' << curve->IsPeriodic()
          << ' ' << curve->NbPoles() << ' ' << curve->NbKnots()
          << ' ' << curve->FirstParameter() << ' ' << curve->LastParameter();
      for (int i=1; i<=curve->NbPoles(); ++i) {
        const gp_Pnt& p=curve->Pole(i);
        out << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z() << ' ' << curve->Weight(i);
      }
      for (int i=1; i<=curve->NbKnots(); ++i) out << ' ' << curve->Knot(i) << ' ' << curve->Multiplicity(i);
      std::cout << out.str() << '\n';
    } catch (const Standard_Failure& e) {
      std::cout << name << " E\n";
      std::cerr << name << ": " << e.GetMessageString() << '\n';
    }
  }
}
