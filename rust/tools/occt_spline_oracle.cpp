// Independent test-only OCCT spline observations; never linked by the kernel.
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <iomanip>
#include <iostream>
#include <string>

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::cout << std::setprecision(17);
  std::string label,kind,side;
  int degree,np,nk,order;
  double u;
  while (std::cin >> label >> kind >> degree >> np >> nk >> u >> side >> order) {
    TColgp_Array1OfPnt poles(1,np);
    TColStd_Array1OfReal weights(1,np),knots(1,nk);
    TColStd_Array1OfInteger mults(1,nk);
    for (int i=1;i<=np;++i) {
      double x,y,z,w;
      if (!(std::cin>>x>>y>>z>>w)) return 2;
      poles(i)=gp_Pnt(x,y,z); weights(i)=w;
    }
    for (int i=1;i<=nk;++i) if (!(std::cin>>knots(i)>>mults(i))) return 2;
    gp_Pnt p;
    gp_Vec d1,d2;
    if (kind=="B") {
      Geom_BezierCurve curve(poles,weights);
      if (order==0) curve.D0(u,p);
      if (order==1) curve.D1(u,p,d1);
      if (order==2) curve.D2(u,p,d1,d2);
    } else if (kind=="S") {
      Geom_BSplineCurve curve(poles,weights,knots,mults,degree,false,false);
      int span=1;
      while (span<nk-1 && (knots(span+1)<u || (knots(span+1)==u && side!="L"))) ++span;
      if (order==0) curve.LocalD0(u,span,span+1,p);
      if (order==1) curve.LocalD1(u,span,span+1,p,d1);
      if (order==2) curve.LocalD2(u,span,span+1,p,d1,d2);
    } else return 2;
    std::cout<<label<<' '<<p.X()<<' '<<p.Y()<<' '<<p.Z();
    if (order>=1) std::cout<<' '<<d1.X()<<' '<<d1.Y()<<' '<<d1.Z();
    if (order>=2) std::cout<<' '<<d2.X()<<' '<<d2.Y()<<' '<<d2.Z();
    std::cout<<'\n';
  }
}
