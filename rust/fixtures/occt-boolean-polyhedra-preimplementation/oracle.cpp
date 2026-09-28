// Source-pinned observations of BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and
// BRepAlgoAPI_Common on two prisms (S9a of REVIEW_NOTES.md). Each input block
// is identity_reference.native_boolean_case: the object's explicit prism
// construction (a `plane` row with the profile's frame at the extrusion's
// start offset, `wire` rows: P polygons, S paths of line and arc segments, C
// circles, holes reversed here; a `prism` vector, as occt_split_oracle.cpp
// reads them), a `boolean fuse|cut|common` row, then the tool's construction
// rows, and `end`. Each prism is BRepPrimAPI_MakePrism of its profile face;
// the Boolean is the object (argument) with the tool.
//
// Output: `NAME done N valid warnings` (N solids in the result, the result
// checked by BRepCheck_Analyzer, 1 if the operation reported warnings), then
// per solid, ordered by volume and centre, `S volume area cx cy cz faces
// edges vertices ufaces uedges uvertices valid` (BRepGProp's volume, surface
// area and centre of mass; counts of the solid as OCCT builds it, then after
// ShapeUpgrade_UnifySameDomain merges its coplanar faces and collinear edges;
// the solid checked alone), or `NAME not_done 0 0 0` or `NAME failure 0 0 0`.
// Nothing is computed from the kernel's output.
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

std::vector<double> numbers(std::istringstream& in, int n) {
  std::vector<double> v(n);
  for (double& x : v)
    if (!(in >> x)) throw Standard_Failure("short row");
  return v;
}

int count(const TopoDS_Shape& s, TopAbs_ShapeEnum type) {
  ShapeMap map;
  TopExp::MapShapes(s, type, map);
  return map.Extent();
}

// One prism's construction rows, as occt_split_oracle.cpp reads them.
struct Prism {
  gp_Ax3 frame;
  std::vector<TopoDS_Wire> wires;
  gp_Vec vector;
  bool has_vector = false;

  void row(const std::string& kind, std::istringstream& in) {
    if (kind == "plane") {
      auto v = numbers(in, 9);
      frame = gp_Ax3(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]), gp_Dir(v[6], v[7], v[8]));
    } else if (kind == "wire") {
      std::string shape;
      in >> shape;
      BRepBuilderAPI_MakeWire mw;
      if (shape == "P" || shape == "S") {
        int n;
        in >> n;
        std::vector<TopoDS_Vertex> vs;
        std::vector<std::vector<double>> arcs;
        for (int j = 0; j < n; ++j) {
          auto p = numbers(in, 3);
          vs.push_back(BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])));
          std::vector<double> arc;
          if (shape == "S") {
            std::string seg;
            in >> seg;
            if (seg == "A")
              arc = numbers(in, 5);
            else if (seg != "L")
              throw Standard_Failure("unknown segment");
          }
          arcs.push_back(arc);
        }
        for (int j = 0; j < n; ++j) {
          const auto& a = arcs[j];
          if (a.empty()) {
            mw.Add(BRepBuilderAPI_MakeEdge(vs[j], vs[(j + 1) % n]));
          } else {
            gp_Dir normal = a[4] != 0 ? frame.Direction() : frame.Direction().Reversed();
            gp_Ax2 axis(gp_Pnt(a[0], a[1], a[2]), normal, frame.XDirection());
            BRepBuilderAPI_MakeEdge me(gp_Circ(axis, a[3]), vs[j], vs[(j + 1) % n]);
            if (!me.IsDone()) throw Standard_Failure("arc edge");
            mw.Add(me.Edge());
          }
        }
      } else {
        auto c = numbers(in, 4);
        gp_Ax2 axis(gp_Pnt(c[0], c[1], c[2]), frame.Direction(), frame.XDirection());
        mw.Add(BRepBuilderAPI_MakeEdge(gp_Circ(axis, c[3])));
      }
      TopoDS_Wire w = mw.Wire();
      if (!wires.empty()) w.Reverse();
      wires.push_back(w);
    } else if (kind == "prism") {
      auto v = numbers(in, 3);
      vector = gp_Vec(v[0], v[1], v[2]);
      has_vector = true;
    } else {
      throw Standard_Failure("unknown row");
    }
  }

  TopoDS_Shape shape() const {
    if (wires.empty() || !has_vector) throw Standard_Failure("incomplete prism");
    BRepBuilderAPI_MakeFace mf(gp_Pln(frame), wires.at(0), true);
    for (size_t i = 1; i < wires.size(); ++i) mf.Add(wires[i]);
    return BRepPrimAPI_MakePrism(mf.Face(), vector).Shape();
  }
};

struct Solid {
  double volume, area;
  gp_Pnt centre;
  int counts[6];
  bool valid;
};
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut, BRepAlgoAPI_Common"
            << std::endl;
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::istringstream head(line);
    std::string tag, name;
    if (!(head >> tag >> name) || tag != "case") return 2;
    std::ostringstream out;
    out << std::setprecision(17);
    try {
      Prism prisms[2];
      int current = 0;
      std::string operation;
      while (std::getline(std::cin, line) && line != "end") {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "boolean") {
          if (current != 0 || !(in >> operation)) throw Standard_Failure("boolean row");
          current = 1;
        } else {
          prisms[current].row(kind, in);
        }
      }
      if (current != 1) throw Standard_Failure("no tool");
      const TopoDS_Shape object = prisms[0].shape(), tool = prisms[1].shape();
      std::unique_ptr<BRepAlgoAPI_BooleanOperation> op;
      if (operation == "fuse")
        op.reset(new BRepAlgoAPI_Fuse(object, tool));
      else if (operation == "cut")
        op.reset(new BRepAlgoAPI_Cut(object, tool));
      else if (operation == "common")
        op.reset(new BRepAlgoAPI_Common(object, tool));
      else
        throw Standard_Failure("unknown operation");
      if (!op->IsDone() || op->HasErrors()) {
        std::cout << name << " not_done 0 0 0\n" << std::flush;
        continue;
      }
      const TopoDS_Shape result = op->Shape();
      std::vector<Solid> solids;
      for (TopExp_Explorer e(result, TopAbs_SOLID); e.More(); e.Next()) {
        GProp_GProps v, a;
        BRepGProp::VolumeProperties(e.Current(), v);
        BRepGProp::SurfaceProperties(e.Current(), a);
        ShapeUpgrade_UnifySameDomain unify(e.Current(), true, true, false);
        unify.Build();
        const TopoDS_Shape merged = unify.Shape();
        solids.push_back({v.Mass(),
                          a.Mass(),
                          v.CentreOfMass(),
                          {count(e.Current(), TopAbs_FACE), count(e.Current(), TopAbs_EDGE),
                           count(e.Current(), TopAbs_VERTEX), count(merged, TopAbs_FACE),
                           count(merged, TopAbs_EDGE), count(merged, TopAbs_VERTEX)},
                          BRepCheck_Analyzer(e.Current()).IsValid() == Standard_True});
      }
      std::sort(solids.begin(), solids.end(), [](const Solid& x, const Solid& y) {
        if (x.volume != y.volume) return x.volume < y.volume;
        if (x.centre.X() != y.centre.X()) return x.centre.X() < y.centre.X();
        if (x.centre.Y() != y.centre.Y()) return x.centre.Y() < y.centre.Y();
        return x.centre.Z() < y.centre.Z();
      });
      out << name << " done " << solids.size() << ' ' << (BRepCheck_Analyzer(result).IsValid() ? 1 : 0) << ' '
          << (op->HasWarnings() ? 1 : 0) << "\n";
      for (const auto& s : solids) {
        out << "S " << s.volume << ' ' << s.area << ' ' << s.centre.X() << ' ' << s.centre.Y() << ' '
            << s.centre.Z();
        for (int c : s.counts) out << ' ' << c;
        out << ' ' << s.valid << "\n";
      }
      std::cout << out.str() << std::flush;
    } catch (const Standard_Failure& failure) {
      std::cout << name << " failure 0 0 0\n" << std::flush;
      while (line != "end" && std::getline(std::cin, line)) {
      }
    }
  }
  return 0;
}
