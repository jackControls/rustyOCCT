// Source-pinned observations of BRepAlgoAPI_Splitter (S8 of REVIEW_NOTES.md):
// each input block is an explicit prism construction produced by
// identity_reference.native_case (a `plane` row with the profile's frame at
// the extrusion's start offset, `wire` rows: P polygons, S paths of line and
// arc segments (S8b: and nonrational B-spline segments `B p n` with n 3D
// poles on the profile plane, `k` knots and multiplicities, each an edge of
// a Geom_BSplineCurve between the path's vertices), C circles, holes reversed
// here; a `prism` vector), or (S8c)
// a `cone ox oy oz nx ny nz xx xy xz r1 r2 h` or `sphere ox oy oz nx ny nz xx
// xy xz R a1 a2` row built by BRepPrimAPI_MakeCone or MakeSphere on that
// gp_Ax2 (S8d: `torus ... R r angle`, MakeTorus), with a `split ox oy oz nx
// ny nz` row: the solid is split by a
// planar face through
// the point with the normal, far larger than the solid. Output: `NAME done N`
// (N result solids), then per solid `S side volume area cx cy cz faces edges
// vertices valid` (side: below or above the plane, by its centre), or
// `NAME failure`. Nothing is computed from the kernel's output.
//
// S8e: a `make face` row in place of the `prism` vector makes the planar face
// of the wires (BRepBuilderAPI_MakeFace on the frame's plane) the argument,
// and each face of the result is a piece: `S side area perimeter cx cy cz 1
// edges vertices valid` (the perimeter BRepGProp::LinearProperties over its
// edges). A `make wire` row makes the one wire the argument. The splitter
// keeps a wire whole with its edges split, so the pieces are grouped here:
// the input wire's edges in their stored order, each replaced by its images
// (`Modified`) in order along it; each image takes the side of its midpoint,
// or, lying on the plane (1e-9), the side of the image before it; each
// maximal cyclic run of images on one side is a piece (every image on the
// plane: one piece below), `S side length 0 cx cy cz 0 edges vertices
// valid` (validity: a wire made of the run's edges in order, checked). The
// images must cover the result's edges.
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <iomanip>
#include <iostream>
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

// A path segment: a line (no data), an arc (centre, radius, ccw) or a
// B-spline curve (S8b).
struct Segment {
  std::vector<double> arc;
  Handle(Geom_BSplineCurve) spline;
};

Handle(Geom_BSplineCurve) spline(std::istringstream& in) {
  int degree, n, k;
  if (!(in >> degree >> n)) throw Standard_Failure("spline header");
  NCollection_Array1<gp_Pnt> poles(1, n);
  for (int i = 1; i <= n; ++i) {
    auto p = numbers(in, 3);
    poles.SetValue(i, gp_Pnt(p[0], p[1], p[2]));
  }
  if (!(in >> k)) throw Standard_Failure("spline knots");
  auto u = numbers(in, k);
  NCollection_Array1<double> knots(1, k);
  NCollection_Array1<int> mults(1, k);
  for (int i = 1; i <= k; ++i) knots.SetValue(i, u[i - 1]);
  for (int i = 1; i <= k; ++i) {
    int m;
    if (!(in >> m)) throw Standard_Failure("spline multiplicities");
    mults.SetValue(i, m);
  }
  return new Geom_BSplineCurve(poles, knots, mults, degree);
}

struct Piece {
  std::string side;
  double volume, area;
  gp_Pnt centre;
  int faces, edges, vertices;
  bool valid;
};

// S8e: a split wire's pieces, the maximal runs of its edges' images on one
// side in the wire's stored order.
std::vector<Piece> wire_runs(const TopoDS_Wire& wire, BRepAlgoAPI_Splitter& splitter, const gp_Pnt& q,
                             const gp_Dir& m) {
  struct Part {
    TopoDS_Edge edge;
    int side;
    double length;
    gp_Pnt centre;
  };
  std::vector<Part> parts;
  ShapeMap covered;
  for (TopoDS_Iterator it(wire); it.More(); it.Next()) {
    const TopoDS_Edge edge = TopoDS::Edge(it.Value());
    double f, l;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, f, l);
    std::vector<TopoDS_Edge> list;
    const TopTools_ListOfShape& modified = splitter.Modified(edge);
    if (modified.IsEmpty()) {
      if (splitter.IsDeleted(edge)) throw Standard_Failure("edge deleted");
      list.push_back(edge);
    } else {
      for (const TopoDS_Shape& s : modified) list.push_back(TopoDS::Edge(s));
    }
    std::vector<std::pair<double, Part>> images;
    for (const TopoDS_Edge& image : list) {
      double a, b;
      Handle(Geom_Curve) c = BRep_Tool::Curve(image, a, b);
      const gp_Pnt mid = c->Value(0.5 * (a + b));
      GeomAPI_ProjectPointOnCurve projection(mid, curve, f, l);
      if (projection.NbPoints() == 0) throw Standard_Failure("image off its edge");
      GProp_GProps g;
      BRepGProp::LinearProperties(image, g);
      const double s = gp_Vec(q, mid).Dot(gp_Vec(m));
      images.push_back({projection.LowerDistanceParameter(),
                        {image, std::abs(s) <= 1e-9 ? 0 : (s < 0 ? -1 : 1), g.Mass(), g.CentreOfMass()}});
      covered.Add(image);
    }
    std::sort(images.begin(), images.end(),
              [](const auto& x, const auto& y) { return x.first < y.first; });
    if (edge.Orientation() == TopAbs_REVERSED) std::reverse(images.begin(), images.end());
    for (const auto& image : images) parts.push_back(image.second);
  }
  if (covered.Extent() != count(splitter.Shape(), TopAbs_EDGE))
    throw Standard_Failure("images do not cover the result");
  const int n = static_cast<int>(parts.size());
  std::vector<int> side(n, -1);
  const bool sided = std::any_of(parts.begin(), parts.end(), [](const Part& p) { return p.side != 0; });
  for (int i = 0; sided && i < n; ++i) {
    int j = i;
    while (parts[j].side == 0) j = (j + n - 1) % n;
    side[i] = parts[j].side;
  }
  int start = 0;
  for (int i = 0; i < n; ++i)
    if (side[i] != side[(i + n - 1) % n]) {
      start = i;
      break;
    }
  std::vector<Piece> pieces;
  for (int k = 0; k < n;) {
    const int s = side[(start + k) % n];
    double length = 0;
    gp_XYZ moment(0, 0, 0);
    ShapeMap edges, vertices;
    BRepBuilderAPI_MakeWire mw;
    for (; k < n && side[(start + k) % n] == s; ++k) {
      const Part& p = parts[(start + k) % n];
      length += p.length;
      moment += p.length * p.centre.XYZ();
      edges.Add(p.edge);
      TopExp::MapShapes(p.edge, TopAbs_VERTEX, vertices);
      mw.Add(p.edge);
    }
    const bool valid = mw.IsDone() && BRepCheck_Analyzer(mw.Wire()).IsValid();
    pieces.push_back({s < 0 ? "below" : "above", length, 0.0, gp_Pnt(moment / length), 0, edges.Extent(),
                      vertices.Extent(), valid});
  }
  return pieces;
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepAlgoAPI_Splitter" << std::endl;
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::istringstream head(line);
    std::string tag, name;
    if (!(head >> tag >> name) || tag != "case") return 2;
    std::ostringstream out;
    out << std::setprecision(17);
    try {
      gp_Ax3 frame;
      std::vector<TopoDS_Wire> wires;
      gp_Vec prism;
      gp_Pnt q;
      gp_Dir m;
      TopoDS_Shape primitive;
      std::string make;
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
          BRepBuilderAPI_MakeWire mw;
          if (shape == "P" || shape == "S") {
            int n;
            in >> n;
            std::vector<TopoDS_Vertex> vs;
            std::vector<Segment> segments;
            for (int j = 0; j < n; ++j) {
              auto p = numbers(in, 3);
              vs.push_back(BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])));
              Segment segment;
              if (shape == "S") {
                std::string seg;
                in >> seg;
                if (seg == "A")
                  segment.arc = numbers(in, 5);
                else if (seg == "B")
                  segment.spline = spline(in);
                else if (seg != "L")
                  throw Standard_Failure("unknown segment");
              }
              segments.push_back(segment);
            }
            for (int j = 0; j < n; ++j) {
              const auto& a = segments[j].arc;
              if (!segments[j].spline.IsNull()) {
                BRepBuilderAPI_MakeEdge me(segments[j].spline, vs[j], vs[(j + 1) % n]);
                if (!me.IsDone()) throw Standard_Failure("spline edge");
                mw.Add(me.Edge());
              } else if (a.empty()) {
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
        } else if (kind == "cone" || kind == "sphere" || kind == "torus") {
          auto v = numbers(in, 12);
          gp_Ax2 axis(gp_Pnt(v[0], v[1], v[2]), gp_Dir(v[3], v[4], v[5]), gp_Dir(v[6], v[7], v[8]));
          if (kind == "cone")
            primitive = BRepPrimAPI_MakeCone(axis, v[9], v[10], v[11]).Shape();
          else if (kind == "sphere")
            primitive = BRepPrimAPI_MakeSphere(axis, v[9], v[10], v[11]).Shape();
          else
            primitive = BRepPrimAPI_MakeTorus(axis, v[9], v[10], v[11]).Shape();
        } else if (kind == "prism") {
          auto v = numbers(in, 3);
          prism = gp_Vec(v[0], v[1], v[2]);
        } else if (kind == "make") {
          in >> make;
          if (make != "face" && make != "wire") throw Standard_Failure("unknown body");
        } else if (kind == "split") {
          auto v = numbers(in, 6);
          q = gp_Pnt(v[0], v[1], v[2]);
          m = gp_Dir(v[3], v[4], v[5]);
        } else {
          throw Standard_Failure("unknown row");
        }
      }
      TopoDS_Shape solid = primitive;
      if (make == "wire") {
        solid = wires.at(0);
      } else if (solid.IsNull()) {
        BRepBuilderAPI_MakeFace mf(gp_Pln(frame), wires.at(0), true);
        for (size_t i = 1; i < wires.size(); ++i) mf.Add(wires[i]);
        if (make == "face")
          solid = mf.Face();
        else
          solid = BRepPrimAPI_MakePrism(mf.Face(), prism).Shape();
      }
      TopoDS_Face tool = BRepBuilderAPI_MakeFace(gp_Pln(q, m), -1e3, 1e3, -1e3, 1e3);
      TopTools_ListOfShape arguments, tools;
      arguments.Append(solid);
      tools.Append(tool);
      BRepAlgoAPI_Splitter splitter;
      splitter.SetArguments(arguments);
      splitter.SetTools(tools);
      splitter.Build();
      if (!splitter.IsDone()) {
        std::cout << name << " not_done 0\n";
        continue;
      }
      std::vector<Piece> pieces;
      if (make == "wire") pieces = wire_runs(TopoDS::Wire(solid), splitter, q, m);
      for (TopExp_Explorer e(splitter.Shape(), TopAbs_FACE); make == "face" && e.More(); e.Next()) {
        GProp_GProps a, l;
        BRepGProp::SurfaceProperties(e.Current(), a);
        BRepGProp::LinearProperties(e.Current(), l, true);
        const gp_Pnt c = a.CentreOfMass();
        const double s = gp_Vec(q, c).Dot(gp_Vec(m));
        pieces.push_back({s < 0 ? "below" : "above", a.Mass(), l.Mass(), c, 1, count(e.Current(), TopAbs_EDGE),
                          count(e.Current(), TopAbs_VERTEX), BRepCheck_Analyzer(e.Current()).IsValid()});
      }
      for (TopExp_Explorer e(splitter.Shape(), TopAbs_SOLID); e.More(); e.Next()) {
        GProp_GProps v, a;
        BRepGProp::VolumeProperties(e.Current(), v);
        BRepGProp::SurfaceProperties(e.Current(), a);
        const gp_Pnt c = v.CentreOfMass();
        const double s = gp_Vec(q, c).Dot(gp_Vec(m));
        pieces.push_back({s < 0 ? "below" : "above", v.Mass(), a.Mass(), c,
                          count(e.Current(), TopAbs_FACE), count(e.Current(), TopAbs_EDGE),
                          count(e.Current(), TopAbs_VERTEX), BRepCheck_Analyzer(e.Current()).IsValid()});
      }
      std::sort(pieces.begin(), pieces.end(), [](const Piece& x, const Piece& y) {
        return x.side != y.side ? x.side < y.side : x.volume < y.volume;
      });
      out << name << " done " << pieces.size() << "\n";
      for (const auto& p : pieces)
        out << "S " << p.side << ' ' << p.volume << ' ' << p.area << ' ' << p.centre.X() << ' '
            << p.centre.Y() << ' ' << p.centre.Z() << ' ' << p.faces << ' ' << p.edges << ' '
            << p.vertices << ' ' << p.valid << "\n";
      std::cout << out.str() << std::flush;
    } catch (const Standard_Failure& failure) {
      std::cout << name << " failure 0\n" << std::flush;
      while (line != "end" && std::getline(std::cin, line)) {
      }
    }
  }
  return 0;
}
