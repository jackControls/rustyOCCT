// Research-only observations of the newer original OCCT ExtremaPC family.
#include <ExtremaPC_Curve.hxx>
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
  std::cout << std::setprecision(17);
  std::string row;
  while (std::getline(std::cin,row)) {
    std::istringstream in(row);
    std::string name,extra;
    int degree,periodic,np,nk;
    double first,last,qx,qy,qz;
    if (!(in>>name>>degree>>periodic>>np>>nk>>first>>last>>qx>>qy>>qz)) return 2;
    try {
      TColgp_Array1OfPnt poles(1,np);
      TColStd_Array1OfReal weights(1,np),knots(1,nk);
      TColStd_Array1OfInteger mults(1,nk);
      for (int i=1;i<=np;++i) {
        double x,y,z,w;if (!(in>>x>>y>>z>>w)) return 2;
        poles(i)=gp_Pnt(x,y,z);weights(i)=w;
      }
      for (int i=1;i<=nk;++i) if (!(in>>knots(i)>>mults(i))) return 2;
      if (in>>extra) return 2;
      Handle(Geom_Curve) curve=new Geom_BSplineCurve(poles,weights,knots,mults,degree,periodic!=0,false);
      ExtremaPC_Curve solver(curve,first,last);
      for (int mode=0;mode<2;++mode) {
        const auto& result=mode?solver.PerformWithEndpoints(gp_Pnt(qx,qy,qz),1e-6):solver.Perform(gp_Pnt(qx,qy,qz),1e-6);
        std::ostringstream text;text<<std::setprecision(17)<<name<<" R "<<mode<<' '<<static_cast<int>(result.Status)
            <<' '<<result.IsDone()<<' '<<result.IsInfinite()<<' '<<result.NbExt();
        for (size_t i=0;i<result.NbExt();++i) {
          const auto& e=result[i];text<<' '<<e.Parameter<<' '<<e.IsMinimum<<' '<<e.IsMaximum<<' '<<e.SquareDistance
              <<' '<<e.Point.X()<<' '<<e.Point.Y()<<' '<<e.Point.Z();
        }
        if (result.IsInfinite()) text<<" I "<<result.InfiniteSquareDistance;
        if (result.NbExt()) text<<" M "<<result.MinIndex()<<' '<<result.MinSquareDistance();
        std::cout<<text.str()<<'\n';
      }
    } catch (const Standard_Failure& failure) {
      std::cout<<name<<" X\n";std::cerr<<name<<": "<<failure.GetMessageString()<<'\n';
    }
  }
}
