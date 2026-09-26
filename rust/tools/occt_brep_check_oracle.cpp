// Source-pinned BRepCheck_Analyzer observations for generic B-rep cases.
// Input blocks are explicit OCCT constructions produced by brep_reference.py:
// every curve parameter already matches its pcurves (SameParameter), wires are
// built in the underlying FORWARD face, and a reversed face is reversed last.
// A seam uses UpdateEdge(E, C_forward, C_reversed, F). No geometry, tolerance
// or orientation is computed here. A second row counts the distinct subshapes.
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax22d.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Use {
  int edge;
  bool forward;
  Handle(Geom2d_Curve) curve;
};

std::vector<double> numbers(std::istringstream& in, int n) {
  std::vector<double> v(n);
  for (auto& x : v)
    if (!(in >> x)) throw Standard_Failure("truncated numbers");
  return v;
}

gp_Ax3 frame(const std::vector<double>& v, int i) {
  return gp_Ax3(gp_Pnt(v[i], v[i + 1], v[i + 2]), gp_Dir(v[i + 3], v[i + 4], v[i + 5]),
                gp_Dir(v[i + 6], v[i + 7], v[i + 8]));
}

void statuses(std::ostream& out, const BRepCheck_Analyzer& ana, const TopoDS_Shape& s,
              const std::string& label) {
  // Results are keyed by oriented shape: collect both orientations. A shape
  // outside the analyzed solid (e.g. an unused vertex) has neither.
  std::set<int> found;
  bool present = false;
  for (const TopoDS_Shape& key : {s.Oriented(TopAbs_FORWARD), s.Oriented(TopAbs_REVERSED)}) {
    Handle(BRepCheck_Result) result;
    try {
      result = ana.Result(key);
    } catch (const Standard_Failure&) {
      continue;
    }
    if (result.IsNull()) continue;
    present = true;
    for (const auto& status : result->Status())
      if (status != BRepCheck_NoError) found.insert(status);
    for (result->InitContextIterator(); result->MoreShapeInContext(); result->NextShapeInContext())
      for (const auto& status : result->StatusOnShape())
        if (status != BRepCheck_NoError) found.insert(status);
  }
  if (!present) {
    out << ' ' << label << ":absent";
    return;
  }
  if (found.empty()) return;
  out << ' ' << label;
  for (int status : found) out << ':' << status;
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream head(line);
    std::string tag, name;
    double tol;
    if (!(head >> tag >> name >> tol) || tag != "case") return 2;
    std::vector<TopoDS_Vertex> vertices;
    std::vector<TopoDS_Edge> edges;
    std::vector<std::pair<double, double>> ranges;
    std::vector<TopoDS_Face> faces;
    std::vector<std::vector<TopoDS_Wire>> wires;
    std::vector<TopoDS_Shell> shells;
    std::ostringstream out;
    out << std::setprecision(17);
    try {
      BRep_Builder b;
      std::map<int, std::vector<Use>> pending;  // pcurves of the current face
      bool reversed = false, open = false;
      auto finish_face = [&]() {
        if (!open) return;
        open = false;
        TopoDS_Face& f = faces.back();
        // A shape is frozen once added, so complete wires join the face here.
        for (const auto& w : wires.back()) b.Add(f, w);
        for (auto& [edge, uses] : pending) {
          const Use* fw = nullptr;
          const Use* rv = nullptr;
          for (const auto& u : uses) {
            if (u.forward && !fw) fw = &u;
            if (!u.forward && !rv) rv = &u;
          }
          if (fw && rv)
            b.UpdateEdge(edges[edge], fw->curve, rv->curve, f, tol);
          else
            b.UpdateEdge(edges[edge], uses.front().curve, f, tol);
        }
        pending.clear();
        if (reversed) f.Orientation(TopAbs_REVERSED);
      };
      while (std::getline(std::cin, line) && line != "end") {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "v") {
          auto p = numbers(in, 3);
          TopoDS_Vertex v;
          b.MakeVertex(v, gp_Pnt(p[0], p[1], p[2]), tol);
          vertices.push_back(v);
        } else if (kind == "e") {
          int s, e;
          std::string curve;
          in >> s >> e >> curve;
          Handle(Geom_Curve) c;
          std::vector<double> v;
          if (curve == "line") {
            v = numbers(in, 8);
            c = new Geom_Line(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]));
          } else {
            v = numbers(in, 12);
            c = new Geom_Circle(frame(v, 0).Ax2(), v[9]);
          }
          double first = v[v.size() - 2], last = v[v.size() - 1];
          TopoDS_Edge edge;
          b.MakeEdge(edge, c, tol);
          b.Add(edge, vertices.at(s).Oriented(TopAbs_FORWARD));
          b.Add(edge, vertices.at(e).Oriented(TopAbs_REVERSED));
          edges.push_back(edge);
          ranges.emplace_back(first, last);
        } else if (kind == "f") {
          finish_face();
          std::string surface, orient;
          in >> surface;
          Handle(Geom_Surface) g;
          if (surface == "plane") {
            auto v = numbers(in, 9);
            g = new Geom_Plane(frame(v, 0));
          } else {
            auto v = numbers(in, 10);
            g = new Geom_CylindricalSurface(frame(v, 0), v[9]);
          }
          in >> orient;
          reversed = orient == "R";
          TopoDS_Face f;
          b.MakeFace(f, g, tol);
          faces.push_back(f);
          wires.emplace_back();
          open = true;
        } else if (kind == "w") {
          TopoDS_Wire w;
          b.MakeWire(w);
          wires.back().push_back(w);
        } else if (kind == "u") {
          int edge;
          std::string orient, curve;
          in >> edge >> orient >> curve;
          Handle(Geom2d_Curve) c;
          if (curve == "line") {
            auto v = numbers(in, 4);
            c = new Geom2d_Line(gp_Pnt2d(v[0], v[1]), gp_Dir2d(v[2], v[3]));
          } else {
            auto v = numbers(in, 6);
            gp_Ax22d axis(gp_Pnt2d(v[0], v[1]), gp_Dir2d(v[2], v[3]), v[5] > 0);
            c = new Geom2d_Circle(axis, v[4]);
          }
          bool forward = orient == "F";
          b.Add(wires.back().back(), edges.at(edge).Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED));
          pending[edge].push_back({edge, forward, c});
        } else if (kind == "s") {
          finish_face();
          TopoDS_Shell s;
          b.MakeShell(s);
          int face;
          while (in >> face) b.Add(s, faces.at(face));
          s.Closed(BRep_Tool::IsClosed(s));
          shells.push_back(s);
        } else {
          throw Standard_Failure("unknown row");
        }
      }
      finish_face();
      for (size_t i = 0; i < edges.size(); ++i) b.Range(edges[i], ranges[i].first, ranges[i].second);
      TopoDS_Solid solid;
      b.MakeSolid(solid);
      for (const auto& s : shells) b.Add(solid, s);
      BRepCheck_Analyzer ana(solid);
      out << name << " R " << ana.IsValid();
      for (size_t i = 0; i < vertices.size(); ++i) statuses(out, ana, vertices[i], "v" + std::to_string(i));
      for (size_t i = 0; i < edges.size(); ++i) statuses(out, ana, edges[i], "e" + std::to_string(i));
      for (size_t i = 0; i < faces.size(); ++i) {
        statuses(out, ana, faces[i], "f" + std::to_string(i));
        for (size_t j = 0; j < wires[i].size(); ++j)
          statuses(out, ana, wires[i][j], "w" + std::to_string(i) + "." + std::to_string(j));
      }
      for (size_t i = 0; i < shells.size(); ++i) statuses(out, ana, shells[i], "s" + std::to_string(i));
      statuses(out, ana, solid, "solid");
      // Distinct subshapes, as DRAW's nbshapes counts them.
      out << '\n' << name << " N";
      for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID}) {
        TopTools_IndexedMapOfShape map;
        TopExp::MapShapes(solid, type, map);
        out << ' ' << map.Extent();
      }
      std::cout << out.str() << '\n';
    } catch (const Standard_Failure& failure) {
      std::cout << name << " E native_exception\n";
      std::cerr << name << ": " << failure.what() << '\n';
      while (line != "end" && std::getline(std::cin, line)) {
      }
    }
  }
}
