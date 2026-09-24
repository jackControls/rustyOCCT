// Source-pinned IntTools_EdgeEdge observations for spline/line and segment queries.
// Rows use the pre-implementation protocol. IntTools_EdgeEdge needs finite
// edges, so a line ("L") uses the finite range [-2e,2e] of the SAME gp_Lin,
// where e bounds every pole's L1 distance from A. Positive weights keep every curve point in the pole hull, and each
// hull point's unit-direction projection has magnitude below e, so no
// intersection is excluded. No tolerance, sampling or result is adjusted.
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Geom_BSplineCurve.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <gp_Lin.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

static void point(std::ostream& out, const gp_Pnt& p) {
  out << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z();
}
int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string row;
  while (std::getline(std::cin,row)) {
    std::istringstream input(row);
    std::string name,kind,extra;
    int degree,periodic,np,nk;
    double first,last;
    if (!(input>>name>>kind>>degree>>periodic>>np>>nk>>first>>last)) return 2;
    try {
      gp_Pnt endpoints[2];
      for (auto& p:endpoints) {
        double x,y,z;if (!(input>>x>>y>>z)) return 2;p=gp_Pnt(x,y,z);
      }
      TColgp_Array1OfPnt poles(1,np);
      TColStd_Array1OfReal weights(1,np),knots(1,nk);
      TColStd_Array1OfInteger mults(1,nk);
      for (int i=1;i<=np;++i) {
        double x,y,z,w;if (!(input>>x>>y>>z>>w)) return 2;
        poles(i)=gp_Pnt(x,y,z);weights(i)=w;
      }
      for (int i=1;i<=nk;++i) if (!(input>>knots(i)>>mults(i))) return 2;
      if (input>>extra || (kind!="L" && kind!="S")) return 2;
      Handle(Geom_BSplineCurve) curve=new Geom_BSplineCurve(poles,weights,knots,mults,degree,periodic!=0,false);
      BRepBuilderAPI_MakeEdge makeCurve(curve,curve->FirstParameter(),curve->LastParameter());
      if (!makeCurve.IsDone()) {
        std::cout<<name<<" E curve_edge "<<static_cast<int>(makeCurve.Error())<<'\n';continue;
      }
      const gp_Lin line(endpoints[0],gp_Dir(gp_Vec(endpoints[0],endpoints[1])));
      const double length=endpoints[0].Distance(endpoints[1]);
      BRepBuilderAPI_MakeEdge makeLine;
      if (kind=="L") {
        double extent=1.;
        for (int i=1;i<=np;++i) {
          const gp_Pnt& p=poles(i);
          extent=std::max(extent,1.+std::abs(p.X()-endpoints[0].X())+std::abs(p.Y()-endpoints[0].Y())+std::abs(p.Z()-endpoints[0].Z()));
        }
        makeLine=BRepBuilderAPI_MakeEdge(line,-2.*extent,2.*extent);
      } else {
        makeLine=BRepBuilderAPI_MakeEdge(line,0.,length);
      }
      if (!makeLine.IsDone()) {
        std::cout<<name<<" E line_edge "<<static_cast<int>(makeLine.Error())<<'\n';continue;
      }
      const auto cedge=makeCurve.Edge(),ledge=makeLine.Edge();
      BRepAdaptor_Curve adaptedLine(ledge);
      IntTools_EdgeEdge inter(cedge,first,last,ledge,adaptedLine.FirstParameter(),adaptedLine.LastParameter());
      inter.SetFuzzyValue(0.);
      inter.UseQuickCoincidenceCheck(false);
      inter.Perform();
      std::ostringstream out;out<<std::setprecision(17)<<name<<" R "<<inter.IsDone()<<' '<<inter.CommonParts().Length();
      for (const auto& part:inter.CommonParts()) {
        const bool curveFirst=part.Edge1().IsSame(cedge);
        if (part.Type()==TopAbs_VERTEX) {
          const double u=part.VertexParameter1(),v=part.VertexParameter2();
          out<<" P "<<curveFirst<<' '<<u<<' '<<v;
          point(out,BRepAdaptor_Curve(part.Edge1()).Value(u));
          point(out,BRepAdaptor_Curve(part.Edge2()).Value(v));
        } else if (part.Type()==TopAbs_EDGE) {
          out<<" I "<<curveFirst<<' '<<part.Range1().First()<<' '<<part.Range1().Last()<<' '<<part.Ranges2().Length();
          for (const auto& range:part.Ranges2()) out<<' '<<range.First()<<' '<<range.Last();
        } else {
          throw Standard_Failure("unexpected common-part type");
        }
      }
      std::cout<<out.str()<<'\n';
    } catch (const Standard_Failure& failure) {
      std::cout<<name<<" E native_exception\n";
      std::cerr<<name<<": "<<failure.what()<<'\n';
    }
  }
}
