// Source-pinned OCCT history observations for revolved primitives (S3 of
// REVIEW_NOTES.md). Each input line is `cone NAME ox oy oz nx ny nz xx xy xz
// r1 r2 h`, the arguments of BRepPrimAPI_MakeCone with a gp_Ax2. The probe
// builds the cone's meridian instead: the polygon (0, 0), (r1, 0), (r2, h),
// (0, h) in the half-plane of the frame's x and axis (a zero radius point is
// left out, as it is the axis point), as a planar face, and revolves it a
// full turn about the axis with BRepPrimAPI_MakeRevol. No history or
// identity is computed here.
//
// For each meridian subshape (vertex j and edge j from vertex j to j+1 in
// that order, then the face) BRepTools_History reports Generated, Modified
// and IsRemoved, and MakeRevol reports FirstShape and LastShape. Every output
// subshape is listed with whether any query reaches it, directly or through a
// reached shape. N rows count distinct subshapes; S rows list the seams and
// the vertices only they use, and D rows the degenerated edges. Shapes are
// printed as geometric signatures; no expected value lives here.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepTools_History.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string num(double x) {
  std::ostringstream s;
  s << std::setprecision(17) << x;
  return s.str();
}

std::string point(const gp_Pnt& p) { return num(p.X()) + " " + num(p.Y()) + " " + num(p.Z()); }

const char* surface_type(GeomAbs_SurfaceType t) {
  switch (t) {
    case GeomAbs_Plane: return "plane";
    case GeomAbs_Cylinder: return "cylinder";
    case GeomAbs_Cone: return "cone";
    case GeomAbs_Sphere: return "sphere";
    case GeomAbs_Torus: return "torus";
    case GeomAbs_SurfaceOfRevolution: return "revolution";
    default: return "other";
  }
}

std::string signature(const TopoDS_Shape& s) {
  switch (s.ShapeType()) {
    case TopAbs_VERTEX:
      return "V " + point(BRep_Tool::Pnt(TopoDS::Vertex(s)));
    case TopAbs_EDGE: {
      const TopoDS_Edge& e = TopoDS::Edge(s);
      if (BRep_Tool::Degenerated(e)) return "E degenerated " + point(BRep_Tool::Pnt(TopExp::FirstVertex(e)));
      BRepAdaptor_Curve c(TopoDS::Edge(s.Oriented(TopAbs_FORWARD)));
      const char* type = c.GetType() == GeomAbs_Line     ? "line"
                         : c.GetType() == GeomAbs_Circle ? "circle"
                                                         : "other";
      double f = c.FirstParameter(), l = c.LastParameter();
      return std::string("E ") + type + " " + point(c.Value(f)) + " " + point(c.Value(0.5 * (f + l))) +
             " " + point(c.Value(l));
    }
    case TopAbs_FACE: {
      GProp_GProps g;
      BRepGProp::SurfaceProperties(s, g);
      return std::string("F ") + surface_type(BRepAdaptor_Surface(TopoDS::Face(s)).GetType()) + " " +
             num(g.Mass()) + " " + point(g.CentreOfMass());
    }
    default: {
      GProp_GProps g;
      BRepGProp::VolumeProperties(s, g);
      return "S " + num(g.Mass()) + " " + point(g.CentreOfMass());
    }
  }
}

using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

// Seam edges (closed on one of their faces) and vertices used only by seams,
// by degenerated edges and by edges closed on them.
ShapeMap structure_only(const TopoDS_Shape& shape) {
  ShapeMap seams, out;
  for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
    for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
      if (BRep_Tool::IsClosed(TopoDS::Edge(e.Current()), TopoDS::Face(f.Current()))) {
        seams.Add(e.Current());
        out.Add(e.Current());
      }
  TopTools_IndexedDataMapOfShapeListOfShape users;
  TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, users);
  for (int i = 1; i <= users.Extent(); ++i) {
    bool seam = false, only = true;
    for (const auto& e : users(i)) {
      const TopoDS_Edge& edge = TopoDS::Edge(e);
      seam = seam || seams.Contains(e);
      only = only && (seams.Contains(e) || BRep_Tool::Degenerated(edge) ||
                      TopExp::FirstVertex(edge).IsSame(TopExp::LastVertex(edge)));
    }
    if (seam && only) out.Add(users.FindKey(i));
  }
  return out;
}

std::string counts(const TopoDS_Shape& shape) {
  std::string out = "N";
  for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID}) {
    ShapeMap map;
    TopExp::MapShapes(shape, type, map);
    out += " " + std::to_string(map.Extent());
  }
  return out;
}

struct Input {
  std::string label;  // "vertex j", "edge j" or "face 0"
  TopoDS_Shape shape;
};
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepPrimAPI_MakeRevol BRepTools_History" << std::endl;
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
      const gp_Pnt o = axes.Location();
      const gp_Vec x(axes.XDirection()), n(axes.Direction());
      auto at = [&](double r, double z) { return o.Translated(x * r + n * z); };
      std::vector<gp_Pnt> polygon{at(0.0, 0.0)};
      if (v[9] > 0.0) polygon.push_back(at(v[9], 0.0));
      if (v[10] > 0.0) polygon.push_back(at(v[10], v[11]));
      polygon.push_back(at(0.0, v[11]));
      const int count = static_cast<int>(polygon.size());
      std::vector<TopoDS_Vertex> vertices;
      for (const auto& p : polygon) vertices.push_back(BRepBuilderAPI_MakeVertex(p));
      std::vector<Input> inputs;
      BRepBuilderAPI_MakeWire wire;
      std::vector<TopoDS_Edge> edges;
      for (int j = 0; j < count; ++j) {
        edges.push_back(BRepBuilderAPI_MakeEdge(vertices[j], vertices[(j + 1) % count]));
        wire.Add(edges.back());
      }
      for (int j = 0; j < count; ++j) inputs.push_back({"vertex " + std::to_string(j), vertices[j]});
      for (int j = 0; j < count; ++j) inputs.push_back({"edge " + std::to_string(j), edges[j]});
      TopoDS_Face face = BRepBuilderAPI_MakeFace(wire.Wire(), true);
      inputs.push_back({"face 0", face});
      out << "R " << name << ' ' << BRepCheck_Analyzer(face).IsValid() << '\n';
      for (const auto& input : inputs) out << "I " << input.label << ' ' << signature(input.shape) << '\n';
      BRepPrimAPI_MakeRevol maker(face, gp_Ax1(o, axes.Direction()));
      TopoDS_Shape result = maker.Shape();
      NCollection_List<TopoDS_Shape> arguments;
      arguments.Append(face);
      BRepTools_History history(arguments, maker);
      ShapeMap reached;
      auto report = [&](const std::string& label, const char* query, const TopoDS_Shape& s) {
        if (s.IsNull()) {
          out << "Q " << label << ' ' << query << " null\n";
          return;
        }
        out << "Q " << label << ' ' << query << ' ' << signature(s) << '\n';
        reached.Add(s);
      };
      for (const auto& input : inputs) {
        for (const auto& s : history.Generated(input.shape)) report(input.label, "gen", s);
        for (const auto& s : history.Modified(input.shape)) report(input.label, "mod", s);
        if (history.IsRemoved(input.shape)) out << "Q " << input.label << " deleted\n";
        if (input.label == "face 0") {
          report(input.label, "first", maker.FirstShape());
          report(input.label, "last", maker.LastShape());
        } else {
          report(input.label, "first", maker.FirstShape(input.shape));
          report(input.label, "last", maker.LastShape(input.shape));
        }
      }
      ShapeMap covered;
      for (int i = 1; i <= reached.Extent(); ++i) TopExp::MapShapes(reached(i), covered);
      out << "B " << signature(result) << ' ' << BRepCheck_Analyzer(result).IsValid() << '\n';
      for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE}) {
        ShapeMap map;
        TopExp::MapShapes(result, type, map);
        for (int i = 1; i <= map.Extent(); ++i)
          out << "O " << signature(map(i)) << ' '
              << (reached.Contains(map(i)) ? "direct" : covered.Contains(map(i)) ? "through" : "none") << '\n';
      }
      out << counts(result) << '\n';
      ShapeMap structure = structure_only(result);
      for (int i = 1; i <= structure.Extent(); ++i) out << "S " << signature(structure(i)) << '\n';
      ShapeMap all_edges;
      TopExp::MapShapes(result, TopAbs_EDGE, all_edges);
      for (int i = 1; i <= all_edges.Extent(); ++i)
        if (BRep_Tool::Degenerated(TopoDS::Edge(all_edges(i)))) out << "D " << signature(all_edges(i)) << '\n';
      std::cout << out.str() << "end\n";
    } catch (const Standard_Failure& failure) {
      std::cout << "R " << name << " E native_exception\nend\n";
      std::cerr << name << ": " << failure.what() << '\n';
    }
  }
  return 0;
}
