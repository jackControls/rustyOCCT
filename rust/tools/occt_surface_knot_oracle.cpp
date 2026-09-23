// Test-only native surface knot observations; never linked to the Rust kernel.
#include <Geom_BSplineSurface.hxx>
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

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string row;
  while (std::getline(std::cin,row)) {
    std::istringstream in(row); std::string name,extra;
    int du,dv,nu,nv,nuk,nvk,pu,pv,nop;
    if (!(in>>name>>du>>dv>>nu>>nv>>nuk>>nvk>>pu>>pv>>nop)) return 2;
    try {
      TColgp_Array2OfPnt poles(1,nu,1,nv); TColStd_Array2OfReal weights(1,nu,1,nv);
      TColStd_Array1OfReal uk(1,nuk),vk(1,nvk); TColStd_Array1OfInteger um(1,nuk),vm(1,nvk);
      for (int i=1;i<=nu;++i) for (int j=1;j<=nv;++j) {
        double x,y,z,w; if (!(in>>x>>y>>z>>w)) return 2;
        poles(i,j)=gp_Pnt(x,y,z); weights(i,j)=w;
      }
      for (int i=1;i<=nuk;++i) if (!(in>>uk(i)>>um(i))) return 2;
      for (int i=1;i<=nvk;++i) if (!(in>>vk(i)>>vm(i))) return 2;
      Handle(Geom_BSplineSurface) s=new Geom_BSplineSurface(poles,weights,uk,vk,um,vm,du,dv,pu!=0,pv!=0);
      std::vector<bool> flags;
      for (int i=0;i<nop;++i) {
        char op,axis; double u; int target;
        if (!(in>>op>>axis>>u>>target) || (axis!='U' && axis!='V')) return 2;
        if (op=='I') {
          if (axis=='U') s->InsertUKnot(u,target,0.,false);
          else s->InsertVKnot(u,target,0.,false);
          flags.push_back(true);
        } else if (op=='R') {
          int index=0, n=axis=='U'?s->NbUKnots():s->NbVKnots();
          for (int k=1;k<=n;++k) if ((axis=='U'?s->UKnot(k):s->VKnot(k))==u) index=k;
          if (!index) return 3;
          flags.push_back(axis=='U'?s->RemoveUKnot(index,target,1.e-9):s->RemoveVKnot(index,target,1.e-9));
        } else return 2;
      }
      if (in>>extra) return 2;
      double ua,ub,va,vb; s->Bounds(ua,ub,va,vb);
      std::ostringstream out; out<<std::setprecision(17)<<name<<" R "<<flags.size();
      for (bool f:flags) out<<' '<<f;
      out<<' '<<s->UDegree()<<' '<<s->VDegree()<<' '<<s->IsUPeriodic()<<' '<<s->IsVPeriodic()
         <<' '<<s->NbUPoles()<<' '<<s->NbVPoles()<<' '<<s->NbUKnots()<<' '<<s->NbVKnots()
         <<' '<<ua<<' '<<ub<<' '<<va<<' '<<vb;
      for (int i=1;i<=s->NbUPoles();++i) for (int j=1;j<=s->NbVPoles();++j) {
        const gp_Pnt& p=s->Pole(i,j); out<<' '<<p.X()<<' '<<p.Y()<<' '<<p.Z()<<' '<<s->Weight(i,j);
      }
      for (int i=1;i<=s->NbUKnots();++i) out<<' '<<s->UKnot(i)<<' '<<s->UMultiplicity(i);
      for (int i=1;i<=s->NbVKnots();++i) out<<' '<<s->VKnot(i)<<' '<<s->VMultiplicity(i);
      std::cout<<out.str()<<'\n';
    } catch (const Standard_Failure& e) {
      std::cout<<name<<" E\n";
      std::cerr<<name<<": "<<e.GetMessageString()<<'\n';
    }
  }
}
