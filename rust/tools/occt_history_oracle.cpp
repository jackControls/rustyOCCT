// Source-pinned OCCT history observations for extrusions and rigid transforms.
// Input blocks are explicit constructions produced by identity_reference.py:
// the profile face at the extrusion's start offset (outer wire and holes in
// counter-clockwise stored order; holes are reversed here), the prism vector
// and zero or more rigid transforms. No history or identity is computed here.
//
// For each profile subshape (per boundary: vertex j and edge j from vertex j
// to j+1; then the face), BRepTools_History wrapping BRepPrimAPI_MakePrism
// reports Generated, Modified and IsRemoved, and MakePrism reports FirstShape
// and LastShape. Every output subshape is listed with whether any of those
// queries reaches it. Each transform step reports Modified for every subshape.
// N rows count distinct subshapes; S rows list the structure-only subshapes
// of the result (step -) and of each transform step's input.
// Shapes are printed as geometric signatures; no expected value lives here.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools_History.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<double> numbers(std::istringstream& in, int n) {
  std::vector<double> v(n);
  for (auto& x : v)
    if (!(in >> x)) throw Standard_Failure("truncated numbers");
  return v;
}

std::string num(double x) {
  std::ostringstream s;
  s << std::setprecision(17) << x;
  return s.str();
}

std::string point(const gp_Pnt& p) { return num(p.X()) + " " + num(p.Y()) + " " + num(p.Z()); }

std::string signature(const TopoDS_Shape& s) {
  switch (s.ShapeType()) {
    case TopAbs_VERTEX:
      return "V " + point(BRep_Tool::Pnt(TopoDS::Vertex(s)));
    case TopAbs_EDGE: {
      BRepAdaptor_Curve c(TopoDS::Edge(s.Oriented(TopAbs_FORWARD)));
      const char* type = c.GetType() == GeomAbs_Line     ? "line"
                         : c.GetType() == GeomAbs_Circle ? "circle"
                                                         : "other";
      double f = c.FirstParameter(), l = c.LastParameter();
      return std::string("E ") + type + " " + point(c.Value(f)) + " " + point(c.Value(0.5 * (f + l))) +
             " " + point(c.Value(l));
    }
    case TopAbs_FACE: {
      BRepAdaptor_Surface a(TopoDS::Face(s));
      const char* type = a.GetType() == GeomAbs_Plane      ? "plane"
                         : a.GetType() == GeomAbs_Cylinder ? "cylinder"
                                                           : "other";
      GProp_GProps g;
      BRepGProp::SurfaceProperties(s, g);
      return std::string("F ") + type + " " + num(g.Mass()) + " " + point(g.CentreOfMass());
    }
    default: {
      GProp_GProps g;
      BRepGProp::VolumeProperties(s, g);
      return "S " + num(g.Mass()) + " " + point(g.CentreOfMass());
    }
  }
}

using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

// OCCT structure a seamless model does not have: seam edges (closed on one of
// their faces) and vertices used only by seams and by edges closed on them.
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
      only = only && (seams.Contains(e) || TopExp::FirstVertex(edge).IsSame(TopExp::LastVertex(edge)));
    }
    if (seam && only) out.Add(users.FindKey(i));
  }
  return out;
}

// Distinct subshapes, as DRAW's nbshapes counts them.
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
  std::string label;  // "b vertex j", "b edge j" or "face 0 0"
  TopoDS_Shape shape;
};
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream head(line);
    std::string tag, name;
    if (!(head >> tag >> name) || tag != "case") return 2;
    std::ostringstream out;
    try {
      gp_Ax3 frame;
      std::vector<TopoDS_Wire> wires;
      std::vector<Input> inputs;
      gp_Vec prism;
      std::vector<gp_Trsf> transforms;
      while (std::getline(std::cin, line) && line != "end") {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "plane") {
          auto v = numbers(in, 9);
          frame = gp_Ax3(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]), gp_Dir(v[6], v[7], v[8]));
        } else if (kind == "wire") {
          std::string shape;
          in >> shape;
          int b = static_cast<int>(wires.size());
          BRepBuilderAPI_MakeWire mw;
          if (shape == "P") {
            int n;
            in >> n;
            std::vector<TopoDS_Vertex> vs;
            for (int j = 0; j < n; ++j) {
              auto p = numbers(in, 3);
              vs.push_back(BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])));
            }
            std::vector<TopoDS_Edge> es;
            for (int j = 0; j < n; ++j) {
              es.push_back(BRepBuilderAPI_MakeEdge(vs[j], vs[(j + 1) % n]));
              mw.Add(es.back());
            }
            for (int j = 0; j < n; ++j)
              inputs.push_back({std::to_string(b) + " vertex " + std::to_string(j), vs[j]});
            for (int j = 0; j < n; ++j)
              inputs.push_back({std::to_string(b) + " edge " + std::to_string(j), es[j]});
          } else {
            auto c = numbers(in, 4);
            gp_Ax2 axis(gp_Pnt(c[0], c[1], c[2]), frame.Direction(), frame.XDirection());
            TopoDS_Edge e = BRepBuilderAPI_MakeEdge(gp_Circ(axis, c[3]));
            mw.Add(e);
            TopoDS_Vertex v = TopExp::FirstVertex(e);
            inputs.push_back({std::to_string(b) + " vertex 0", v});
            inputs.push_back({std::to_string(b) + " edge 0", e});
          }
          TopoDS_Wire w = mw.Wire();
          if (b > 0) w.Reverse();
          wires.push_back(w);
        } else if (kind == "prism") {
          auto v = numbers(in, 3);
          prism = gp_Vec(v[0], v[1], v[2]);
        } else if (kind == "transform") {
          auto m = numbers(in, 12);
          gp_Trsf t;
          t.SetValues(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
          transforms.push_back(t);
        } else {
          throw Standard_Failure("unknown row");
        }
      }
      BRepBuilderAPI_MakeFace mf(gp_Pln(frame), wires.at(0), true);
      for (size_t i = 1; i < wires.size(); ++i) mf.Add(wires[i]);
      TopoDS_Face face = mf.Face();
      inputs.push_back({"face 0 0", face});
      out << "R " << name << ' ' << BRepCheck_Analyzer(face).IsValid() << '\n';
      BRepPrimAPI_MakePrism maker(face, prism);
      TopoDS_Shape result = maker.Shape();
      NCollection_List<TopoDS_Shape> arguments;
      arguments.Append(face);
      BRepTools_History history(arguments, maker);
      ShapeMap reached;
      auto report = [&](const std::string& label, const char* query, const TopoDS_Shape& s) {
        out << "Q " << label << ' ' << query << ' ' << signature(s) << '\n';
        reached.Add(s);
      };
      for (const auto& input : inputs) {
        for (const auto& s : history.Generated(input.shape)) report(input.label, "gen", s);
        for (const auto& s : history.Modified(input.shape)) report(input.label, "mod", s);
        if (history.IsRemoved(input.shape)) out << "Q " << input.label << " deleted\n";
        if (input.label == "face 0 0") {
          report(input.label, "first", maker.FirstShape());
          report(input.label, "last", maker.LastShape());
        } else {
          report(input.label, "first", maker.FirstShape(input.shape));
          report(input.label, "last", maker.LastShape(input.shape));
        }
      }
      // Faces and edges reached through a cap count as covered subshapes.
      ShapeMap covered;
      for (int i = 1; i <= reached.Extent(); ++i) TopExp::MapShapes(reached(i), covered);
      out << "B " << signature(result) << ' ' << BRepCheck_Analyzer(result).IsValid() << '\n';
      for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE}) {
        ShapeMap map;
        TopExp::MapShapes(result, type, map);
        for (int i = 1; i <= map.Extent(); ++i)
          out << "O " << signature(map(i)) << ' ' << (reached.Contains(map(i)) ? "direct" : covered.Contains(map(i)) ? "through" : "none")
              << '\n';
      }
      out << counts(result) << '\n';
      ShapeMap structure = structure_only(result);
      for (int i = 1; i <= structure.Extent(); ++i) out << "S - " << signature(structure(i)) << '\n';
      TopoDS_Shape current = result;
      for (size_t k = 0; k < transforms.size(); ++k) {
        BRepBuilderAPI_Transform moved(current, transforms[k], true);
        TopoDS_Shape next = moved.Shape();
        ShapeMap before = structure_only(current);
        for (int i = 1; i <= before.Extent(); ++i) out << "S " << k << ' ' << signature(before(i)) << '\n';
        for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE}) {
          ShapeMap map;
          TopExp::MapShapes(current, type, map);
          for (int i = 1; i <= map.Extent(); ++i) {
            const auto& images = moved.Modified(map(i));
            out << "T " << k << ' ' << signature(map(i)) << " -> " << images.Extent() << '\n';
            for (const auto& s : images) out << "  " << signature(s) << '\n';
          }
        }
        out << "B " << signature(next) << ' ' << BRepCheck_Analyzer(next).IsValid() << '\n';
        current = next;
      }
      std::cout << out.str() << "end\n";
    } catch (const Standard_Failure& failure) {
      std::cout << "R " << name << " E native_exception\nend\n";
      std::cerr << name << ": " << failure.what() << '\n';
      while (line != "end" && std::getline(std::cin, line)) {
      }
    }
  }
}
