// Independent test-only OCCT surface observations; never linked by the kernel.
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <iomanip>
#include <iostream>
#include <string>

int span(double u, const std::string& side, const TColStd_Array1OfReal& knots) {
  int i=1;
  while (i<knots.Upper()-1 && (knots(i+1)<u || (knots(i+1)==u && side!="L"))) ++i;
  return i;
}
void seam(double& u, const std::string& side, const TColStd_Array1OfReal& knots) {
  if (u==knots(1) && side=="L") u=knots(knots.Upper());
  if (u==knots(knots.Upper()) && side!="L") u=knots(1);
}
int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::cout << std::setprecision(17);
  std::string label,kind,su,sv;
  int du,dv,nu,nv,ku,kv,pu,pv,order;
  double u,v;
  while (std::cin>>label>>kind>>du>>dv>>nu>>nv>>ku>>kv>>pu>>pv>>u>>v>>su>>sv>>order) {
    TColgp_Array2OfPnt poles(1,nu,1,nv);
    TColStd_Array2OfReal weights(1,nu,1,nv);
    TColStd_Array1OfReal uk(1,ku),vk(1,kv);
    TColStd_Array1OfInteger um(1,ku),vm(1,kv);
    for (int i=1;i<=nu;++i) for (int j=1;j<=nv;++j) {
      double x,y,z,w;
      if (!(std::cin>>x>>y>>z>>w)) return 2;
      poles(i,j)=gp_Pnt(x,y,z); weights(i,j)=w;
    }
    for (int i=1;i<=ku;++i) if (!(std::cin>>uk(i)>>um(i))) return 2;
    for (int i=1;i<=kv;++i) if (!(std::cin>>vk(i)>>vm(i))) return 2;
    gp_Pnt p;
    gp_Vec d1u,d1v,d2u,d2v,d2uv;
    if (kind=="B") {
      Geom_BezierSurface surface(poles,weights);
      if (order==0) surface.D0(u,v,p);
      if (order==1) surface.D1(u,v,p,d1u,d1v);
      if (order==2) surface.D2(u,v,p,d1u,d1v,d2u,d2v,d2uv);
    } else if (kind=="S") {
      Geom_BSplineSurface surface(poles,weights,uk,vk,um,vm,du,dv,pu!=0,pv!=0);
      surface.PeriodicNormalization(u,v);
      if (pu) seam(u,su,uk);
      if (pv) seam(v,sv,vk);
      double uf,ul,vf,vl;
      surface.Bounds(uf,ul,vf,vl);
      if (!pu && u==ul && su=="A") su="L";
      if (!pv && v==vl && sv=="A") sv="L";
      int i=span(u,su,uk),j=span(v,sv,vk);
      if (order==0) surface.LocalD0(u,v,i,i+1,j,j+1,p);
      if (order==1) surface.LocalD1(u,v,i,i+1,j,j+1,p,d1u,d1v);
      if (order==2) surface.LocalD2(u,v,i,i+1,j,j+1,p,d1u,d1v,d2u,d2v,d2uv);
    } else return 2;
    std::cout<<label<<' '<<p.X()<<' '<<p.Y()<<' '<<p.Z();
    if (order>=1) for (auto d: {d1u,d1v}) std::cout<<' '<<d.X()<<' '<<d.Y()<<' '<<d.Z();
    if (order>=2) for (auto d: {d2u,d2v,d2uv}) std::cout<<' '<<d.X()<<' '<<d.Y()<<' '<<d.Z();
    std::cout<<'\n';
  }
}
