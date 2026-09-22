// Test-only native Bézier extraction/editing. Never linked to the Rust kernel.
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
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

struct Arc { Handle(Geom_BezierCurve) curve; double first,last; };
static Arc trim(const Arc& a,double low,double high) {
  Handle(Geom_BezierCurve) curve=Handle(Geom_BezierCurve)::DownCast(a.curve->Copy());
  curve->Segment(low,high);
  return {curve,a.first+(a.last-a.first)*low,a.first+(a.last-a.first)*high};
}
int main() {
  std::cerr<<"OCCT "<<OCC_VERSION_COMPLETE<<'\n';
  std::string row;
  while (std::getline(std::cin,row)) {
    std::istringstream input(row);std::string name,kind,extra;
    int degree,np,nk,op,elevation;double first,last;
    if (!(input>>name>>kind>>degree>>np>>nk>>first>>last>>op>>elevation)) return 2;
    try {
      TColgp_Array1OfPnt poles(1,np);TColStd_Array1OfReal weights(1,np),knots(1,nk);TColStd_Array1OfInteger mults(1,nk);
      for (int i=1;i<=np;++i) {
        double x,y,z,w;if (!(input>>x>>y>>z>>w)) return 2;
        poles(i)=gp_Pnt(x,y,z);weights(i)=w;
      }
      for (int i=1;i<=nk;++i) if (!(input>>knots(i)>>mults(i))) return 2;
      if (input>>extra) return 2;
      Handle(Geom_BSplineCurve) spline=new Geom_BSplineCurve(poles,weights,knots,mults,degree,kind=="P",false);
      GeomConvert_BSplineCurveToBezierCurve decomposition(spline,first,last,0.);
      TColStd_Array1OfReal boundaries(1,decomposition.NbArcs()+1);decomposition.Knots(boundaries);
      std::vector<Arc> arcs;
      for (int i=1;i<=decomposition.NbArcs();++i) {
        Arc a{decomposition.Arc(i),boundaries(i),boundaries(i+1)};
        if (op==1 || op==5) a=trim(a,.25,.75);
        if (op==4 || op==5) a.curve->Increase(elevation);
        if (op==3 || op==5) a.curve->Reverse();
        if (op==2 || op==5) {
          arcs.push_back(trim(a,0.,.375));arcs.push_back(trim(a,.375,1.));
        } else arcs.push_back(a);
      }
      std::ostringstream out;out<<std::setprecision(17)<<name<<" R "<<arcs.size();
      for (const Arc& a:arcs) {
        out<<' '<<a.curve->Degree()<<' '<<a.first<<' '<<a.last;
        for (int i=1;i<=a.curve->NbPoles();++i) {
          const gp_Pnt& p=a.curve->Pole(i);
          out<<' '<<p.X()<<' '<<p.Y()<<' '<<p.Z()<<' '<<a.curve->Weight(i);
        }
      }
      std::cout<<out.str()<<'\n';
    } catch (const Standard_Failure& e) {
      std::cout<<name<<" E\n";std::cerr<<name<<": "<<e.GetMessageString()<<'\n';
    }
  }
}
