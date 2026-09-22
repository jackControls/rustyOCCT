// Test-only GeomAPI_IntCS observations; never linked into the kernel.
#include <GeomAPI_IntCS.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <gp_Pln.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::cerr<<"OCCT "<<OCC_VERSION_COMPLETE<<'\n';
  std::cout<<std::setprecision(17);
  std::string name,kind;
  int degree,np,nk;
  while (std::cin>>name>>kind>>degree>>np>>nk) {
    gp_Pnt plane_points[3];
    for (auto& p:plane_points) {
      double x,y,z; if (!(std::cin>>x>>y>>z)) return 2; p=gp_Pnt(x,y,z);
    }
    TColgp_Array1OfPnt poles(1,np);
    TColStd_Array1OfReal weights(1,np),knots(1,nk);
    TColStd_Array1OfInteger mults(1,nk);
    for (int i=1;i<=np;++i) {
      double x,y,z,w; if (!(std::cin>>x>>y>>z>>w)) return 2;
      poles(i)=gp_Pnt(x,y,z);weights(i)=w;
    }
    for (int i=1;i<=nk;++i) if (!(std::cin>>knots(i)>>mults(i))) return 2;
    const bool trimmed=kind.size()==2 && kind.back()=='T';
    double first=0,last=0;
    if (trimmed && !(std::cin>>first>>last)) return 2;
    try {
      Handle(Geom_Curve) curve;
      if (kind.front()=='B') curve=new Geom_BezierCurve(poles,weights);
      else curve=new Geom_BSplineCurve(poles,weights,knots,mults,degree,kind.front()=='P',false);
      if (trimmed) curve=new Geom_TrimmedCurve(curve,first,last,true,false);
      const gp_Vec normal=gp_Vec(plane_points[0],plane_points[1]).Crossed(gp_Vec(plane_points[0],plane_points[2]));
      Handle(Geom_Surface) plane=new Geom_Plane(gp_Pln(plane_points[0],gp_Dir(normal)));
      GeomAPI_IntCS inter(curve,plane);
      if (!inter.IsDone()) {std::cout<<name<<" FAIL\n";continue;}
      std::vector<std::array<double,4>> points;
      std::vector<std::array<double,2>> spans;
      for (int i=1;i<=inter.NbPoints();++i) {
        double u,v,t;inter.Parameters(i,u,v,t);
        const auto& p=inter.Point(i); points.push_back({t,p.X(),p.Y(),p.Z()});
      }
      for (int i=1;i<=inter.NbSegments();++i) {
        const auto c=inter.Segment(i); spans.push_back({c->FirstParameter(),c->LastParameter()});
      }
      std::sort(points.begin(),points.end());std::sort(spans.begin(),spans.end());
      std::cout<<name<<' '<<points.size()<<' '<<spans.size();
      for (auto p:points) for (auto x:p) std::cout<<' '<<x;
      for (auto s:spans) for (auto x:s) std::cout<<' '<<x;
      std::cout<<'\n';
    } catch (const Standard_Failure&) {std::cout<<name<<" FAIL\n";}
  }
}
