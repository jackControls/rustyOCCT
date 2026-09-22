// Test-only OCCT distance observations. Never linked into the Rust kernel.
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Vec.hxx>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

struct Primitive {
  char kind;
  gp_Pnt p[3];
  explicit Primitive(std::istream& input) {
    if (!(input >> kind)) throw Standard_Failure("missing primitive kind");
    int count = kind == 'P' ? 1 : (kind == 'L' || kind == 'S') ? 2 : 3;
    for (int i = 0; i < count; ++i) {
      double x,y,z;
      if (!(input >> x >> y >> z)) throw Standard_Failure("missing primitive point");
      p[i].SetCoord(x,y,z);
    }
  }
  gp_Lin line() const { return gp_Lin(p[0],gp_Dir(gp_Vec(p[0],p[1]))); }
  gp_Pln plane() const { return gp_Pln(p[0],gp_Dir(gp_Vec(p[0],p[1]).Crossed(gp_Vec(p[0],p[2])))); }
};

static TopoDS_Shape shape(const Primitive& input) {
  const char kind=input.kind;
  const gp_Pnt* p=input.p;
  if (kind == 'P') return BRepBuilderAPI_MakeVertex(p[0]).Vertex();
  // A collapsed closed segment is the singleton point as a geometric set.
  if (kind == 'S' && p[0].IsEqual(p[1],0.)) return BRepBuilderAPI_MakeVertex(p[0]).Vertex();
  if (kind == 'S') return BRepBuilderAPI_MakeEdge(p[0],p[1]).Edge();
  if (kind == 'L') return BRepBuilderAPI_MakeEdge(input.line()).Edge();
  const gp_Pln plane=input.plane();
  if (kind == 'F') return BRepBuilderAPI_MakeFace(plane).Face();
  if (kind == 'T') {
    TopoDS_Wire wire = BRepBuilderAPI_MakePolygon(p[0],p[1],p[2],true).Wire();
    return BRepBuilderAPI_MakeFace(plane,wire,true).Face();
  }
  throw Standard_Failure("invalid primitive kind");
}

static std::optional<double> affine_distance(const Primitive& a,const Primitive& b) {
  if (a.kind=='P' && b.kind=='P') return a.p[0].Distance(b.p[0]);
  if (a.kind=='L' && b.kind=='P') return a.line().Distance(b.p[0]);
  if (a.kind=='P' && b.kind=='L') return b.line().Distance(a.p[0]);
  if (a.kind=='F' && b.kind=='P') return a.plane().Distance(b.p[0]);
  if (a.kind=='P' && b.kind=='F') return b.plane().Distance(a.p[0]);
  if (a.kind=='L' && b.kind=='L') return a.line().Distance(b.line());
  if (a.kind=='L' && b.kind=='F') return b.plane().Distance(a.line());
  if (a.kind=='F' && b.kind=='L') return a.plane().Distance(b.line());
  if (a.kind=='F' && b.kind=='F') return a.plane().Distance(b.plane());
  return std::nullopt;
}

int main(int argc,char** argv) {
  bool affine=argc==2 && std::string(argv[1])=="--affine";
  if (argc!=1 && !affine) return 2;
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::cout << std::setprecision(17);
  std::string line;
  while (std::getline(std::cin,line)) {
    std::istringstream input(line);
    std::string name,extra;
    if (!(input >> name)) return 2;
    try {
      const Primitive first(input),second(input);
      if (input >> extra) return 2;
      if (affine) {
        auto distance=affine_distance(first,second);
        if (first.kind==second.kind && (first.kind=='L' || first.kind=='F')) {
          const gp_Dir a=first.kind=='L'?first.line().Direction():first.plane().Axis().Direction();
          const gp_Dir b=second.kind=='L'?second.line().Direction():second.plane().Axis().Direction();
          std::cerr << std::setprecision(17) << name << " resolution=" << gp::Resolution()
                    << " angle=" << a.Angle(b) << " parallel=" << a.IsParallel(b,gp::Resolution());
          for (const gp_Dir* p : {&a,&b})
            std::cerr << " direction=" << p->X() << ',' << p->Y() << ',' << p->Z();
          std::cerr << '\n';
        }
        std::cout << name;
        if (distance) std::cout << " P " << *distance << '\n';
        else std::cout << " N\n";
        continue;
      }
      TopoDS_Shape a=shape(first), b=shape(second);
      BRepExtrema_DistShapeShape result(a,b);
      if (!result.IsDone()) { std::cout << name << " F\n"; continue; }
      std::cout << name << " P " << result.Value() << ' ' << result.NbSolution();
      for (int i=1;i<=result.NbSolution();++i) {
        for (gp_Pnt p : {result.PointOnShape1(i),result.PointOnShape2(i)})
          std::cout << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z();
      }
      std::cout << '\n';
    } catch (const Standard_Failure& failure) {
      std::cout << name << " E\n";
      std::cerr << name << ": " << failure.GetMessageString() << '\n';
    }
  }
  return 0;
}
