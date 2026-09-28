// Source-pinned BRepCheck_Analyzer observations for generic B-rep cases.
// Input blocks are explicit OCCT constructions produced by brep_reference.py:
// every curve parameter already matches its pcurves (SameParameter), wires are
// built in the underlying FORWARD face, and a reversed face is reversed last.
// A seam uses UpdateEdge(E, C_forward, C_reversed, F). No geometry, tolerance
// or orientation is computed here. A second row counts the distinct subshapes.
// A third row (M5 observations) gives each vertex's BRep_Tool::Tolerance and
// its largest distance to the curve ends of its edges, each edge's
// tolerance, each use's curve-on-surface deviation (BRepLib_ValidateEdge,
// exact method) and each face's tolerance. B-spline curves, 2D curves and
// surfaces (S4) are given by degree, periodicity, knots, multiplicities,
// poles and weights, in OCCT's own conventions. With BREP_ORACLE_PROPERTIES
// set, a fourth row gives BRepGProp's volume and surface area with their
// error estimates at a requested relative precision of 1e-12, the centre of
// mass and the matrix of inertia about it. A `result` row (S6) checks a free
// face, a shell, a wire (`W` edges), an edge or a vertex instead of a solid of the
// shells; a `G` row then gives its area or length and centre.
#include <Adaptor3d_CurveOnSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepLib_ValidateEdge.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColgp_Array2OfPnt.hxx>
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
#include <gp_Mat.hxx>
#include <gp_Pnt.hxx>
#include <cmath>
#include <cstdlib>
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

// A use as built: face, wire and position in the wire, edge, orientation.
struct Placed {
  size_t face, wire, index;
  int edge;
  bool forward;
};

// The exact-method maximum distance between an edge's 3D curve and its
// curve on the face for the use's orientation; NaN if not computed.
double use_deviation(const TopoDS_Edge& edge, bool forward, const TopoDS_Face& face) {
  const TopoDS_Edge oriented =
      TopoDS::Edge(edge.Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED));
  const TopoDS_Face unoriented = TopoDS::Face(face.Oriented(TopAbs_FORWARD));
  double f, l;
  Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface(oriented, unoriented, f, l);
  if (pc.IsNull()) return std::nan("");
  Handle(Geom2dAdaptor_Curve) pca = new Geom2dAdaptor_Curve(pc, f, l);
  Handle(GeomAdaptor_Surface) sa = new GeomAdaptor_Surface(BRep_Tool::Surface(unoriented));
  Handle(Adaptor3d_CurveOnSurface) cos = new Adaptor3d_CurveOnSurface(pca, sa);
  Handle(BRepAdaptor_Curve) c3 = new BRepAdaptor_Curve(edge);
  BRepLib_ValidateEdge validate(c3, cos, true);
  validate.SetExactMethod(true);
  validate.Process();
  return validate.IsDone() ? validate.GetMaxDistance() : std::nan("");
}

std::vector<double> numbers(std::istringstream& in, int n) {
  std::vector<double> v(n);
  for (auto& x : v)
    if (!(in >> x)) throw Standard_Failure("truncated numbers");
  return v;
}

// `DEG PERIODIC NK knots... mults...`.
struct Basis {
  int degree;
  bool periodic;
  TColStd_Array1OfReal knots;
  TColStd_Array1OfInteger mults;
};

Basis basis(std::istringstream& in) {
  int degree, periodic, n;
  if (!(in >> degree >> periodic >> n) || n < 2) throw Standard_Failure("bad basis");
  Basis b{degree, periodic == 1, TColStd_Array1OfReal(1, n), TColStd_Array1OfInteger(1, n)};
  for (int i = 1; i <= n; ++i)
    if (!(in >> b.knots(i))) throw Standard_Failure("truncated knots");
  for (int i = 1; i <= n; ++i)
    if (!(in >> b.mults(i))) throw Standard_Failure("truncated multiplicities");
  return b;
}

int count(std::istringstream& in) {
  int n;
  if (!(in >> n) || n < 1) throw Standard_Failure("bad pole count");
  return n;
}

Handle(Geom_Curve) bspline3(std::istringstream& in) {
  Basis b = basis(in);
  int n = count(in);
  TColgp_Array1OfPnt poles(1, n);
  for (int i = 1; i <= n; ++i) {
    auto p = numbers(in, 3);
    poles(i) = gp_Pnt(p[0], p[1], p[2]);
  }
  TColStd_Array1OfReal weights(1, n);
  for (int i = 1; i <= n; ++i) weights(i) = numbers(in, 1)[0];
  return new Geom_BSplineCurve(poles, weights, b.knots, b.mults, b.degree, b.periodic);
}

Handle(Geom2d_Curve) bspline2(std::istringstream& in) {
  Basis b = basis(in);
  int n = count(in);
  TColgp_Array1OfPnt2d poles(1, n);
  for (int i = 1; i <= n; ++i) {
    auto p = numbers(in, 2);
    poles(i) = gp_Pnt2d(p[0], p[1]);
  }
  TColStd_Array1OfReal weights(1, n);
  for (int i = 1; i <= n; ++i) weights(i) = numbers(in, 1)[0];
  return new Geom2d_BSplineCurve(poles, weights, b.knots, b.mults, b.degree, b.periodic);
}

// Poles U-major (index u * v_count + v), as the kernel stores them.
Handle(Geom_Surface) bspline_surface(std::istringstream& in) {
  Basis u = basis(in), v = basis(in);
  int n = count(in);
  auto poles_in = [](const Basis& b) {
    int total = 0;
    for (int i = b.mults.Lower(); i <= b.mults.Upper(); ++i) total += b.mults(i);
    return b.periodic ? total - b.mults(b.mults.Lower()) : total - b.degree - 1;
  };
  int nu = poles_in(u), nv = poles_in(v);
  if (nu * nv != n) throw Standard_Failure("pole grid");
  TColgp_Array2OfPnt poles(1, nu, 1, nv);
  TColStd_Array2OfReal weights(1, nu, 1, nv);
  for (int i = 1; i <= nu; ++i)
    for (int j = 1; j <= nv; ++j) {
      auto p = numbers(in, 3);
      poles(i, j) = gp_Pnt(p[0], p[1], p[2]);
    }
  for (int i = 1; i <= nu; ++i)
    for (int j = 1; j <= nv; ++j) weights(i, j) = numbers(in, 1)[0];
  return new Geom_BSplineSurface(poles, weights, u.knots, v.knots, u.mults, v.mults, u.degree,
                                 v.degree, u.periodic, v.periodic);
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
    // S6: what is checked, when not a solid of the shells: a free face,
    // a shell, a wire of listed edges or a vertex.
    std::string result_kind = "solid";
    TopoDS_Wire free_wire;
    std::vector<std::pair<int, int>> ends;  // start and end vertex per edge
    std::vector<Placed> placed;
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
          if (curve == "bspline") {
            c = bspline3(in);
            v = numbers(in, 2);
          } else if (curve == "line") {
            v = numbers(in, 8);
            c = new Geom_Line(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]));
          } else {
            v = numbers(in, 12);
            c = new Geom_Circle(frame(v, 0).Ax2(), v[9]);
          }
          double first = v[v.size() - 2], last = v[v.size() - 1];
          TopoDS_Edge edge;
          b.MakeEdge(edge, c, tol);
          ends.emplace_back(s, e);
          b.Add(edge, vertices.at(s).Oriented(TopAbs_FORWARD));
          b.Add(edge, vertices.at(e).Oriented(TopAbs_REVERSED));
          edges.push_back(edge);
          ranges.emplace_back(first, last);
        } else if (kind == "f") {
          finish_face();
          std::string surface, orient;
          in >> surface;
          Handle(Geom_Surface) g;
          if (surface == "bspline") {
            g = bspline_surface(in);
          } else if (surface == "plane") {
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
          if (curve == "bspline") {
            c = bspline2(in);
          } else if (curve == "line") {
            auto v = numbers(in, 4);
            c = new Geom2d_Line(gp_Pnt2d(v[0], v[1]), gp_Dir2d(v[2], v[3]));
          } else {
            auto v = numbers(in, 6);
            gp_Ax22d axis(gp_Pnt2d(v[0], v[1]), gp_Dir2d(v[2], v[3]), v[5] > 0);
            c = new Geom2d_Circle(axis, v[4]);
          }
          bool forward = orient == "F";
          size_t index = 0;
          for (const auto& p : placed)
            if (p.face == faces.size() - 1 && p.wire == wires.back().size() - 1) ++index;
          placed.push_back({faces.size() - 1, wires.back().size() - 1, index, edge, forward});
          b.Add(wires.back().back(), edges.at(edge).Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED));
          pending[edge].push_back({edge, forward, c});
        } else if (kind == "W") {
          b.MakeWire(free_wire);
          int edge;
          while (in >> edge) b.Add(free_wire, edges.at(edge));
        } else if (kind == "result") {
          in >> result_kind;
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
      TopoDS_Shape solid;
      if (result_kind == "solid") {
        TopoDS_Solid made;
        b.MakeSolid(made);
        for (const auto& s : shells) b.Add(made, s);
        solid = made;
      } else if (result_kind == "face") {
        solid = faces.at(0);
      } else if (result_kind == "shell") {
        solid = shells.at(0);
      } else if (result_kind == "wire") {
        solid = free_wire;
      } else if (result_kind == "edge") {
        solid = edges.at(0);
      } else if (result_kind == "vertex") {
        solid = vertices.at(0);
      } else {
        throw Standard_Failure("unknown result");
      }
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
      if (result_kind == "wire") statuses(out, ana, free_wire, "wire");
      if (result_kind == "solid") statuses(out, ana, solid, "solid");
      // Distinct subshapes, as DRAW's nbshapes counts them.
      out << '\n' << name << " N";
      for (TopAbs_ShapeEnum type : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID}) {
        TopTools_IndexedMapOfShape map;
        TopExp::MapShapes(solid, type, map);
        out << ' ' << map.Extent();
      }
      // Tolerances and measured deviations, as observations (T6).
      out << '\n' << name << " T";
      for (size_t i = 0; i < vertices.size(); ++i) {
        const gp_Pnt p = BRep_Tool::Pnt(vertices[i]);
        double gap = 0.0;
        for (size_t k = 0; k < edges.size(); ++k) {
          double f, l;
          Handle(Geom_Curve) c = BRep_Tool::Curve(edges[k], f, l);
          if (c.IsNull()) continue;
          if (ends[k].first == int(i)) gap = std::max(gap, p.Distance(c->Value(f)));
          if (ends[k].second == int(i)) gap = std::max(gap, p.Distance(c->Value(l)));
        }
        out << " v" << i << ':' << BRep_Tool::Tolerance(vertices[i]) << ':' << gap;
      }
      for (size_t i = 0; i < edges.size(); ++i) out << " e" << i << ':' << BRep_Tool::Tolerance(edges[i]);
      for (const auto& p : placed)
        out << " u" << p.face << '.' << p.wire << '.' << p.index << ':'
            << use_deviation(edges[p.edge], p.forward, faces[p.face]);
      for (size_t i = 0; i < faces.size(); ++i) out << " f" << i << ':' << BRep_Tool::Tolerance(faces[i]);
      if (result_kind != "solid") {
        // S6: the area (faces, shells) or length (wires) and its centre.
        GProp_GProps g;
        if (result_kind == "wire" || result_kind == "edge")
          BRepGProp::LinearProperties(solid, g);
        else if (result_kind != "vertex")
          BRepGProp::SurfaceProperties(solid, g);
        const gp_Pnt c = result_kind == "vertex" ? BRep_Tool::Pnt(vertices.at(0)) : g.CentreOfMass();
        out << '\n' << name << " G " << (result_kind == "vertex" ? 0.0 : g.Mass()) << ' ' << c.X() << ' '
            << c.Y() << ' ' << c.Z();
      }
      if (std::getenv("BREP_ORACLE_PROPERTIES")) {
        GProp_GProps volume, surface;
        const double volume_error = BRepGProp::VolumeProperties(solid, volume, 1e-12);
        const double area_error = BRepGProp::SurfaceProperties(solid, surface, 1e-12);
        const gp_Pnt c = volume.CentreOfMass();
        const gp_Mat m = volume.MatrixOfInertia();
        out << '\n' << name << " P " << volume.Mass() << ' ' << volume_error << ' ' << surface.Mass() << ' '
            << area_error << ' ' << c.X() << ' ' << c.Y() << ' ' << c.Z();
        for (int a = 1; a <= 3; ++a)
          for (int b = 1; b <= 3; ++b) out << ' ' << m.Value(a, b);
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
