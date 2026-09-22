// Test-only OCCT surface extraction/editing. Never linked to the Rust kernel.
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BezierCurve.hxx>
#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>

struct Patch { Handle(Geom_BezierSurface) surface; double ua,ub,va,vb; };
static Patch trim(const Patch& p,double ua,double ub,double va,double vb) {
  Handle(Geom_BezierSurface) s=Handle(Geom_BezierSurface)::DownCast(p.surface->Copy());
  s->Segment(ua,ub,va,vb);
  return {s,p.ua+(p.ub-p.ua)*ua,p.ua+(p.ub-p.ua)*ub,p.va+(p.vb-p.va)*va,p.va+(p.vb-p.va)*vb};
}
static void emit(std::ostream& out,const Patch& p,int op) {
  if (op==8 || op==9) {
    Handle(Geom_BezierCurve) c=Handle(Geom_BezierCurve)::DownCast(op==8 ? p.surface->UIso(.375) : p.surface->VIso(.625));
    out<<" C "<<c->Degree()<<" 0 "<<(op==8?p.va:p.ua)<<' '<<(op==8?p.vb:p.ub);
    for(int i=1;i<=c->NbPoles();++i) {
      const gp_Pnt& q=c->Pole(i); out<<' '<<q.X()<<' '<<q.Y()<<' '<<q.Z()<<' '<<c->Weight(i);
    }
  } else {
    out<<" P "<<p.surface->UDegree()<<' '<<p.surface->VDegree()<<' '<<p.ua<<' '<<p.ub<<' '<<p.va<<' '<<p.vb;
    for(int i=1;i<=p.surface->NbUPoles();++i) for(int j=1;j<=p.surface->NbVPoles();++j) {
      const gp_Pnt& q=p.surface->Pole(i,j); out<<' '<<q.X()<<' '<<q.Y()<<' '<<q.Z()<<' '<<p.surface->Weight(i,j);
    }
  }
}
int main() {
  std::cerr<<"OCCT "<<OCC_VERSION_COMPLETE<<'\n'; std::string row;
  while(std::getline(std::cin,row)) {
    std::istringstream in(row); std::string name,kind,extra;
    int du,dv,nu,nv,nuk,nvk,pu,pv,op,eu,ev; double ua,ub,va,vb;
    if(!(in>>name>>kind>>du>>dv>>nu>>nv>>nuk>>nvk>>pu>>pv>>ua>>ub>>va>>vb>>op>>eu>>ev)) return 2;
    try {
      TColgp_Array2OfPnt poles(1,nu,1,nv); TColStd_Array2OfReal weights(1,nu,1,nv);
      TColStd_Array1OfReal uk(1,nuk),vk(1,nvk); TColStd_Array1OfInteger um(1,nuk),vm(1,nvk);
      for(int i=1;i<=nu;++i) for(int j=1;j<=nv;++j) {
        double x,y,z,w; if(!(in>>x>>y>>z>>w)) return 2;
        poles(i,j)=gp_Pnt(x,y,z); weights(i,j)=w;
      }
      for(int i=1;i<=nuk;++i) if(!(in>>uk(i)>>um(i))) return 2;
      for(int i=1;i<=nvk;++i) if(!(in>>vk(i)>>vm(i))) return 2;
      if(in>>extra) return 2;
      Handle(Geom_BSplineSurface) spline=new Geom_BSplineSurface(poles,weights,uk,vk,um,vm,du,dv,pu!=0,pv!=0);
      GeomConvert_BSplineSurfaceToBezierSurface decomposition(spline,ua,ub,va,vb,0.);
      TColStd_Array1OfReal us(1,decomposition.NbUPatches()+1),vs(1,decomposition.NbVPatches()+1);
      decomposition.UKnots(us); decomposition.VKnots(vs);
      // OCCT may return equivalent periodic knots in its fundamental turn.
      // Preserve the caller's original parameter units using the actual offset.
      double uperiod=uk(nuk)-uk(1), vperiod=vk(nvk)-vk(1);
      double uoffset=pu ? std::round((ua-us(1))/uperiod)*uperiod : 0.;
      double voffset=pv ? std::round((va-vs(1))/vperiod)*vperiod : 0.;
      std::vector<Patch> patches;
      for(int i=1;i<=decomposition.NbUPatches();++i) for(int j=1;j<=decomposition.NbVPatches();++j) {
        Patch p{decomposition.Patch(i,j),us(i)+uoffset,us(i+1)+uoffset,vs(j)+voffset,vs(j+1)+voffset};
        if(op==1 || op==10) p=trim(p,.25,.75,.25,.75);
        if(op==7 || op==10) p.surface->Increase(eu,ev);
        if(op==4 || op==10) p.surface->UReverse();
        if(op==5) p.surface->VReverse();
        if(op==6 || op==10) { p.surface->ExchangeUV(); std::swap(p.ua,p.va); std::swap(p.ub,p.vb); }
        if(op==2) { patches.push_back(trim(p,0.,.375,0.,1.)); patches.push_back(trim(p,.375,1.,0.,1.)); }
        else if(op==3) { patches.push_back(trim(p,0.,1.,0.,.625)); patches.push_back(trim(p,0.,1.,.625,1.)); }
        else if(op==10) {
          patches.push_back(trim(p,0.,.375,0.,.625)); patches.push_back(trim(p,0.,.375,.625,1.));
          patches.push_back(trim(p,.375,1.,0.,.625)); patches.push_back(trim(p,.375,1.,.625,1.));
        } else patches.push_back(p);
      }
      std::ostringstream out; out<<std::setprecision(17)<<name<<" R "<<patches.size();
      for(const auto& p:patches) emit(out,p,op);
      std::cout<<out.str()<<'\n';
    } catch(const Standard_Failure& e) { std::cout<<name<<" E\n"; std::cerr<<name<<": "<<e.GetMessageString()<<'\n'; }
  }
}
