// Source-pinned observations of OCCT's primitives (S3 of REVIEW_NOTES.md).
// Each input line is `cone NAME ox oy oz nx ny nz xx xy xz r1 r2 h`, built
// with BRepPrimAPI_MakeCone(gp_Ax2(origin, normal, x), r1, r2, h). For each
// solid it prints the BRepCheck_Analyzer verdict and the distinct subshape
// counts, the BRepGProp volume, surface area, centre of mass and matrix of
// inertia about the centre, and every face (surface type, area, centre),
// edge (degenerated, closed, length, point at mid parameter) and vertex, in
// TopExp::MapShapes order. Nothing is computed from any kernel output.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <gp_Ax2.hxx>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {
const char* surface_type(GeomAbs_SurfaceType t) {
  switch (t) {
    case GeomAbs_Plane: return "plane";
    case GeomAbs_Cylinder: return "cylinder";
    case GeomAbs_Cone: return "cone";
    case GeomAbs_Sphere: return "sphere";
    case GeomAbs_Torus: return "torus";
    default: return "other";
  }
}

void point(std::ostream& out, const gp_Pnt& p) { out << ' ' << p.X() << ' ' << p.Y() << ' ' << p.Z(); }
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepPrimAPI BRepCheck_Analyzer BRepGProp" << std::endl;
  std::cout << std::setprecision(17);
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string kind, name;
    double v[12];
    if (!(in >> kind >> name) || kind != "cone") return 2;
    for (double& x : v)
      if (!(in >> x)) return 2;
    std::ostringstream out;
    out << std::setprecision(17);
    try {
      gp_Ax2 axes(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]), gp_Dir(v[6], v[7], v[8]));
      BRepPrimAPI_MakeCone make(axes, v[9], v[10], v[11]);
      const TopoDS_Shape& solid = make.Shape();
      BRepCheck_Analyzer analyzer(solid);
      out << "case " << name << ' ' << (analyzer.IsValid() ? "valid" : "invalid");
      TopTools_IndexedMapOfShape vertices, edges, faces;
      TopExp::MapShapes(solid, TopAbs_VERTEX, vertices);
      TopExp::MapShapes(solid, TopAbs_EDGE, edges);
      TopExp::MapShapes(solid, TopAbs_FACE, faces);
      for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID}) {
        TopTools_IndexedMapOfShape map;
        TopExp::MapShapes(solid, type, map);
        out << ' ' << map.Extent();
      }
      GProp_GProps volume, area;
      BRepGProp::VolumeProperties(solid, volume);
      BRepGProp::SurfaceProperties(solid, area);
      const gp_Mat m = volume.MatrixOfInertia();
      out << "\nprops " << volume.Mass() << ' ' << area.Mass();
      point(out, volume.CentreOfMass());
      for (int r = 1; r <= 3; ++r)
        for (int c = r; c <= 3; ++c) out << ' ' << m.Value(r, c);
      for (int i = 1; i <= faces.Extent(); ++i) {
        const TopoDS_Face& f = TopoDS::Face(faces(i));
        GProp_GProps g;
        BRepGProp::SurfaceProperties(f, g);
        out << "\nface " << i << ' ' << surface_type(BRepAdaptor_Surface(f).GetType()) << ' ' << g.Mass();
        point(out, g.CentreOfMass());
      }
      for (int i = 1; i <= edges.Extent(); ++i) {
        const TopoDS_Edge& e = TopoDS::Edge(edges(i));
        const bool degenerated = BRep_Tool::Degenerated(e);
        GProp_GProps g;
        BRepGProp::LinearProperties(e, g);
        gp_Pnt mid;
        if (degenerated) {
          mid = BRep_Tool::Pnt(TopExp::FirstVertex(e));
        } else {
          BRepAdaptor_Curve c(e);
          mid = c.Value(0.5 * (c.FirstParameter() + c.LastParameter()));
        }
        TopoDS_Vertex a, b;
        TopExp::Vertices(e, a, b);
        out << "\nedge " << i << ' ' << (degenerated ? "degenerated" : "regular") << ' '
            << (a.IsSame(b) ? "closed" : "open") << ' ' << g.Mass();
        point(out, mid);
      }
      for (int i = 1; i <= vertices.Extent(); ++i) {
        out << "\nvertex " << i;
        point(out, BRep_Tool::Pnt(TopoDS::Vertex(vertices(i))));
      }
      std::cout << out.str() << "\nend\n";
    } catch (const Standard_Failure& failure) {
      std::cout << "case " << name << " exception\nend\n";
      std::cerr << name << ": " << failure.what() << '\n';
    }
  }
  return 0;
}
