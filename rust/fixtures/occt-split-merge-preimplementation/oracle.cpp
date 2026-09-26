// Source-pinned OCCT observations for the height split and the stacked fuse.
// Input blocks come from generate_split_merge_fixtures.native_rows: named
// prisms (the profile face at the start offset, the prism vector and rigid
// transforms, as in occt_history_oracle.cpp), then one operation:
//   split P point normal half      BRepAlgoAPI_Splitter of P by a square
//                                  planar face tool
//   fuse A B                       BRepAlgoAPI_Fuse, then
//                                  ShapeUpgrade_UnifySameDomain
//   splitfuse P point normal half  the split, then the fuse and unification
//                                  of its two solids
// For every stage, each argument subshape (vertex, edge, face, solid in
// TopExp::MapShapes order) is queried in the stage's BRepTools_History
// (merged across the algorithms of that stage): Modified, Generated and
// IsRemoved. Each result solid is listed with its validity and distinct
// subshape counts, and every result vertex, edge and face with whether it is
// an argument subshape kept as is, an image of a query, or neither. Shapes are
// geometric signatures; no expected value or history is computed here.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
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
#include <NCollection_IndexedMap.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
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
#include <map>
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
using Shapes = NCollection_List<TopoDS_Shape>;

// One prism from its construction rows.
struct Spec {
  gp_Ax3 frame;
  std::vector<TopoDS_Wire> wires;
  gp_Vec prism;
  std::vector<gp_Trsf> transforms;
};

TopoDS_Shape build(const Spec& s) {
  BRepBuilderAPI_MakeFace mf(gp_Pln(s.frame), s.wires.at(0), true);
  for (size_t i = 1; i < s.wires.size(); ++i) mf.Add(s.wires[i]);
  TopoDS_Shape shape = BRepPrimAPI_MakePrism(mf.Face(), s.prism).Shape();
  for (const auto& t : s.transforms) shape = BRepBuilderAPI_Transform(shape, t, true).Shape();
  return shape;
}

void row(Spec& s, const std::string& kind, std::istringstream& in) {
  if (kind == "plane") {
    auto v = numbers(in, 9);
    s.frame = gp_Ax3(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]), gp_Dir(v[6], v[7], v[8]));
  } else if (kind == "wire") {
    std::string shape;
    in >> shape;
    BRepBuilderAPI_MakeWire mw;
    if (shape == "P") {
      int n;
      in >> n;
      std::vector<TopoDS_Vertex> vs;
      for (int j = 0; j < n; ++j) {
        auto p = numbers(in, 3);
        vs.push_back(BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])));
      }
      for (int j = 0; j < n; ++j) mw.Add(BRepBuilderAPI_MakeEdge(vs[j], vs[(j + 1) % n]));
    } else {
      auto c = numbers(in, 4);
      gp_Ax2 axis(gp_Pnt(c[0], c[1], c[2]), s.frame.Direction(), s.frame.XDirection());
      mw.Add(BRepBuilderAPI_MakeEdge(gp_Circ(axis, c[3])));
    }
    TopoDS_Wire w = mw.Wire();
    if (!s.wires.empty()) w.Reverse();
    s.wires.push_back(w);
  } else if (kind == "prism") {
    auto v = numbers(in, 3);
    s.prism = gp_Vec(v[0], v[1], v[2]);
  } else if (kind == "transform") {
    auto m = numbers(in, 12);
    gp_Trsf t;
    t.SetValues(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
    s.transforms.push_back(t);
  } else {
    throw Standard_Failure("unknown row");
  }
}

std::string counts(const TopoDS_Shape& shape) {
  std::string out;
  for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID}) {
    ShapeMap map;
    TopExp::MapShapes(shape, type, map);
    out += " " + std::to_string(map.Extent());
  }
  return out;
}

// Queries of every argument subshape in one stage's history, then the result.
void report(std::ostream& out, const std::string& stage, const std::vector<std::pair<std::string, TopoDS_Shape>>& args,
            const BRepTools_History& history, const TopoDS_Shape& result) {
  ShapeMap images, inputs;
  for (const auto& [name, arg] : args) {
    const char* tags[] = {"v", "e", "f", "s"};
    int t = 0;
    for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID}) {
      ShapeMap map;
      TopExp::MapShapes(arg, type, map);
      for (int i = 1; i <= map.Extent(); ++i) {
        const TopoDS_Shape& s = map(i);
        inputs.Add(s);
        const Shapes& mod = history.Modified(s);
        const Shapes& gen = history.Generated(s);
        out << "Q " << stage << ' ' << name << ' ' << tags[t] << ' ' << i - 1 << " mod " << mod.Extent()
            << " gen " << gen.Extent() << " del " << history.IsRemoved(s) << " | " << signature(s) << '\n';
        for (const auto& m : mod) {
          out << "  M " << signature(m) << '\n';
          images.Add(m);
        }
        for (const auto& g : gen) {
          out << "  G " << signature(g) << '\n';
          images.Add(g);
        }
      }
      ++t;
    }
  }
  int k = 0;
  for (TopExp_Explorer e(result, TopAbs_SOLID); e.More(); e.Next(), ++k) {
    out << "B " << stage << ' ' << k << ' ' << signature(e.Current()) << ' '
        << BRepCheck_Analyzer(e.Current()).IsValid() << '\n';
    out << "N " << stage << ' ' << k << counts(e.Current()) << '\n';
  }
  for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE}) {
    ShapeMap map;
    TopExp::MapShapes(result, type, map);
    for (int i = 1; i <= map.Extent(); ++i)
      out << "O " << stage << ' ' << signature(map(i)) << ' '
          << (inputs.Contains(map(i)) ? "kept" : images.Contains(map(i)) ? "image" : "none") << '\n';
  }
}

// Split one shape by a square planar face.
TopoDS_Shape split(const TopoDS_Shape& shape, std::istringstream& in, BRepTools_History& history) {
  auto v = numbers(in, 7);
  gp_Pln plane(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]));
  TopoDS_Face tool = BRepBuilderAPI_MakeFace(plane, -v[6], v[6], -v[6], v[6]).Face();
  BRepAlgoAPI_Splitter splitter;
  Shapes arguments, tools;
  arguments.Append(shape);
  tools.Append(tool);
  splitter.SetArguments(arguments);
  splitter.SetTools(tools);
  splitter.Build();
  if (!splitter.IsDone() || splitter.HasErrors()) throw Standard_Failure("splitter failed");
  history.Merge(splitter.History());
  return splitter.Shape();
}

TopoDS_Shape fuse(const TopoDS_Shape& a, const TopoDS_Shape& b, BRepTools_History& history) {
  BRepAlgoAPI_Fuse fused(a, b);
  if (!fused.IsDone() || fused.HasErrors()) throw Standard_Failure("fuse failed");
  history.Merge(fused.History());
  ShapeUpgrade_UnifySameDomain unify(fused.Shape(), true, true, true);
  unify.Build();
  history.Merge(unify.History());
  return unify.Shape();
}
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
      std::map<std::string, Spec> specs;
      std::map<std::string, TopoDS_Shape> bodies;
      std::string current;
      out << "R " << name << '\n';
      while (std::getline(std::cin, line) && line != "end") {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "body") {
          in >> current;
          specs[current];
        } else if (kind == "split" || kind == "splitfuse" || kind == "fuse") {
          for (const auto& [n, s] : specs) {
            bodies[n] = build(s);
            out << "I " << n << ' ' << signature(bodies[n]) << ' ' << BRepCheck_Analyzer(bodies[n]).IsValid()
                << '\n';
          }
          if (kind == "fuse") {
            std::string a, b;
            in >> a >> b;
            BRepTools_History history;
            TopoDS_Shape result = fuse(bodies.at(a), bodies.at(b), history);
            report(out, "fuse", {{a, bodies.at(a)}, {b, bodies.at(b)}}, history, result);
          } else {
            std::string p;
            in >> p;
            BRepTools_History history;
            TopoDS_Shape pieces = split(bodies.at(p), in, history);
            report(out, "split", {{p, bodies.at(p)}}, history, pieces);
            if (kind == "splitfuse") {
              std::vector<TopoDS_Shape> solids;
              for (TopExp_Explorer e(pieces, TopAbs_SOLID); e.More(); e.Next()) solids.push_back(e.Current());
              if (solids.size() != 2) throw Standard_Failure("split did not give two solids");
              BRepTools_History fused;
              TopoDS_Shape result = fuse(solids[0], solids[1], fused);
              report(out, "fuse", {{"L", solids[0]}, {"U", solids[1]}}, fused, result);
              history.Merge(fused);
              report(out, "composed", {{p, bodies.at(p)}}, history, result);
            }
          }
        } else {
          row(specs.at(current), kind, in);
        }
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
