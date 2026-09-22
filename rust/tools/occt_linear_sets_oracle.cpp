// Test-only observations of the native COMMON and SECTION results.
// No native code is linked to the Rust kernel.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRep_TVertex.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Vec.hxx>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static TopoDS_Shape primitive(std::istream& input) {
  char kind;
  if (!(input >> kind) || std::string("PLSFT").find(kind)==std::string::npos)
    throw Standard_Failure("invalid primitive kind");
  const int count=kind=='P'?1:(kind=='L'||kind=='S')?2:3;
  gp_Pnt p[3];
  for (int i=0;i<count;++i) {
    double x,y,z;
    if (!(input>>x>>y>>z)) throw Standard_Failure("invalid coordinates");
    p[i]=gp_Pnt(x,y,z);
  }
  if (kind=='P' || (kind=='S' && p[0].IsEqual(p[1],0.)))
    return BRepBuilderAPI_MakeVertex(p[0]).Vertex();
  if (kind=='S') return BRepBuilderAPI_MakeEdge(p[0],p[1]).Edge();
  if (kind=='L') return BRepBuilderAPI_MakeEdge(gp_Lin(p[0],gp_Dir(gp_Vec(p[0],p[1])))).Edge();
  const gp_Pln plane(p[0],gp_Dir(gp_Vec(p[0],p[1]).Crossed(gp_Vec(p[0],p[2]))));
  if (kind=='F') return BRepBuilderAPI_MakeFace(plane).Face();
  const TopoDS_Wire wire=BRepBuilderAPI_MakePolygon(p[0],p[1],p[2],true).Wire();
  return BRepBuilderAPI_MakeFace(plane,wire,true).Face();
}

static void point(std::ostream& out,const gp_Pnt& p) {
  out << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z();
}
static void direction(std::ostream& out,const gp_Dir& p) {
  out << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z();
}
// Same stored point + location operation as BRep_Tool::Pnt. Using the
// public BRep_TVertex accessors avoids pulling Poly_Triangulation headers into
// this modeling-only oracle (the packaged Ubuntu 7.6 modeling SDK omits
// NCollection_AliasedArray.hxx, a dependency of those triangulation headers).
static gp_Pnt vertex_point(const TopoDS_Vertex& vertex) {
  const Handle(BRep_TVertex) storage=Handle(BRep_TVertex)::DownCast(vertex.TShape());
  if (storage.IsNull()) throw Standard_Failure("non-BRep vertex");
  return vertex.Location().IsIdentity() ? storage->Pnt()
      : storage->Pnt().Transformed(vertex.Location().Transformation());
}

static std::string describe(const TopoDS_Shape& shape) {
  std::ostringstream out;
  out << std::setprecision(17);
  std::vector<std::string> pieces;
  for (TopExp_Explorer it(shape,TopAbs_FACE);it.More();it.Next()) {
    const TopoDS_Face face=TopoDS::Face(it.Current());
    const BRepAdaptor_Surface surface(face);
    if (surface.GetType()!=GeomAbs_Plane) throw Standard_Failure("nonlinear face result");
    TopTools_IndexedMapOfShape vertices;
    TopExp::MapShapes(face,TopAbs_VERTEX,vertices);
    out.str("");
    if (vertices.IsEmpty()) {
      out << "F"; point(out,surface.Plane().Location()); direction(out,surface.Plane().Axis().Direction());
    } else {
      out << "G " << vertices.Extent();
      for (int i=1;i<=vertices.Extent();++i) point(out,vertex_point(TopoDS::Vertex(vertices(i))));
    }
    pieces.push_back(out.str());
  }
  for (TopExp_Explorer it(shape,TopAbs_EDGE,TopAbs_FACE);it.More();it.Next()) {
    const BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
    if (curve.GetType()!=GeomAbs_Line) throw Standard_Failure("nonlinear edge result");
    const double first=curve.FirstParameter(),last=curve.LastParameter();
    out.str("");
    if (Precision::IsInfinite(first) && Precision::IsInfinite(last)) {
      out << "L"; point(out,curve.Line().Location()); direction(out,curve.Line().Direction());
    } else if (Precision::IsInfinite(first) || Precision::IsInfinite(last)) {
      throw Standard_Failure("unexpected ray result");
    } else {
      out << "S"; point(out,curve.Value(first)); point(out,curve.Value(last));
    }
    pieces.push_back(out.str());
  }
  for (TopExp_Explorer it(shape,TopAbs_VERTEX,TopAbs_EDGE);it.More();it.Next()) {
    out.str(""); out << "P"; point(out,vertex_point(TopoDS::Vertex(it.Current())));
    pieces.push_back(out.str());
  }
  out.str(""); out << "R " << pieces.size();
  for (const auto& piece:pieces) out << ' ' << piece;
  return out.str();
}

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string row;
  while (std::getline(std::cin,row)) {
    std::istringstream input(row);
    std::string name,extra;
    if (!(input>>name)) return 2;
    try {
      const TopoDS_Shape a=primitive(input),b=primitive(input);
      if (input>>extra) return 2;
      for (const char op : {'C','S'}) {
        try {
          std::string result;
          if (op=='C') {
            BRepAlgoAPI_Common common(a,b);
            result=common.IsDone()?describe(common.Shape()):"F";
          } else {
            BRepAlgoAPI_Section section(a,b);
            result=section.IsDone()?describe(section.Shape()):"F";
          }
          std::cout << name << ' ' << op << ' ' << result << '\n';
        } catch (const Standard_Failure& error) {
          std::cout << name << ' ' << op << " E\n";
          std::cerr << name << ' ' << op << ": " << error.GetMessageString() << '\n';
        }
      }
    } catch (const Standard_Failure& error) {
      std::cout << name << " C E\n" << name << " S E\n";
      std::cerr << name << ": " << error.GetMessageString() << '\n';
    }
  }
}
