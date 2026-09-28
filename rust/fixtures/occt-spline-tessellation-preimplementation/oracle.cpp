// Source-pinned observations of BRepMesh_IncrementalMesh (tessellation T-a of
// REVIEW_NOTES.md). Stdin: `body NAME`, the lines of a .brep text (the
// kernel's writer's), `endbody`; then `mesh NAME SETTING DEFLECTION ANGLE`
// rows, each meshing a clean copy of that body with an absolute linear
// deflection and an angle, sequentially (isRelative false, isInParallel
// false). Per row it prints
//
//   S NAME SETTING status faces unmeshed nodes welded triangles degenerate
//     free nonmanifold misoriented weld_gap face_deflection edge_deflection
//     face_distance edge_distance area volume
//
// nodes: the sum over faces of Poly_Triangulation::NbNodes; welded: the
// nodes left after joining, through each edge's Poly_PolygonOnTriangulation
// on every face using it, the nodes at the same edge parameter index, every
// polygon end with its topological vertex, and every node of a degenerated
// edge's polygon with its vertex (weld_gap: the largest distance moved);
// degenerate: triangles with two equal welded nodes (left out of the rest);
// free, nonmanifold: mesh edges in one or more than two triangles;
// misoriented: directed mesh edges used twice (a face's triangles are taken
// reversed when the face is); face_deflection, edge_deflection: OCCT's own
// Poly_Triangulation and Poly_PolygonOnTriangulation deflections, largest;
// face_distance: the largest distance from 12 barycentric samples of every
// triangle to its face's surface (ElSLib projections on planes, cylinders,
// cones, spheres and tori); edge_distance: from the quarter points of every
// edge polygon segment to its curve; area and volume of the welded mesh.
// With the argument `dump` it also prints each welded mesh as `mesh NAME
// SETTING deflection angle`, `v x y z`, `t a b c face bound 0` (bound: the
// face's Poly_Triangulation deflection), `e edge bound 0 nodes...`, `end`.
// Nothing is computed from the kernel's output beyond reading its text.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <ElSLib.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Circ.hxx>
#include <gp_Lin.hxx>

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Find {
  std::vector<int> up;
  int root(int a) {
    while (up[a] != a) a = up[a] = up[up[a]];
    return a;
  }
  // The smaller index stays the root: deterministic.
  void join(int a, int b) {
    a = root(a);
    b = root(b);
    if (a != b) up[std::max(a, b)] = std::min(a, b);
  }
};

double surface_distance(const BRepAdaptor_Surface& s, const TopoDS_Face& face, const gp_Pnt& p) {
  double u = 0, v = 0;
  switch (s.GetType()) {
    case GeomAbs_Plane:
      return s.Plane().Distance(p);
    case GeomAbs_Cylinder:
      ElSLib::Parameters(s.Cylinder(), p, u, v);
      return ElSLib::Value(u, v, s.Cylinder()).Distance(p);
    case GeomAbs_Cone:
      ElSLib::Parameters(s.Cone(), p, u, v);
      return ElSLib::Value(u, v, s.Cone()).Distance(p);
    case GeomAbs_Sphere:
      ElSLib::Parameters(s.Sphere(), p, u, v);
      return ElSLib::Value(u, v, s.Sphere()).Distance(p);
    case GeomAbs_Torus:
      ElSLib::Parameters(s.Torus(), p, u, v);
      return ElSLib::Value(u, v, s.Torus()).Distance(p);
    default: {
      GeomAPI_ProjectPointOnSurf proj(p, BRep_Tool::Surface(face));
      return proj.NbPoints() > 0 ? proj.LowerDistance() : -1.0;
    }
  }
}

double curve_distance(const BRepAdaptor_Curve& c, const TopoDS_Edge& edge, const gp_Pnt& p) {
  switch (c.GetType()) {
    case GeomAbs_Line:
      return c.Line().Distance(p);
    case GeomAbs_Circle:
      return c.Circle().Distance(p);
    default: {
      double f = 0, l = 0;
      occ::handle<Geom_Curve> g = BRep_Tool::Curve(edge, f, l);
      if (g.IsNull()) return 0.0;
      GeomAPI_ProjectPointOnCurve proj(p, g);
      return proj.NbPoints() > 0 ? proj.LowerDistance() : -1.0;
    }
  }
}

struct Body {
  TopoDS_Shape shape;
};

void mesh(const std::string& name, const std::string& setting, const TopoDS_Shape& shape,
          double deflection, double angle, bool dump) {
  BRepTools::Clean(shape, true);
  BRepMesh_IncrementalMesh mesher(shape, deflection, false, angle, false);
  TopTools_IndexedMapOfShape faces, edges, vertices;
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  TopExp::MapShapes(shape, TopAbs_EDGE, edges);
  TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
  // Global nodes: every face's triangulation nodes, then one per vertex.
  std::vector<gp_Pnt> nodes;
  std::vector<int> base(faces.Extent() + 1, 0);
  std::vector<occ::handle<Poly_Triangulation>> tris(faces.Extent() + 1);
  std::vector<TopLoc_Location> locs(faces.Extent() + 1);
  int unmeshed = 0, triangles = 0;
  double face_deflection = 0, edge_deflection = 0;
  for (int i = 1; i <= faces.Extent(); ++i) {
    const TopoDS_Face face = TopoDS::Face(faces(i));
    tris[i] = BRep_Tool::Triangulation(face, locs[i]);
    base[i] = static_cast<int>(nodes.size());
    if (tris[i].IsNull()) {
      ++unmeshed;
      continue;
    }
    face_deflection = std::max(face_deflection, tris[i]->Deflection());
    for (int k = 1; k <= tris[i]->NbNodes(); ++k)
      nodes.push_back(tris[i]->Node(k).Transformed(locs[i].Transformation()));
    triangles += tris[i]->NbTriangles();
  }
  const int raw = static_cast<int>(nodes.size());
  for (int v = 1; v <= vertices.Extent(); ++v)
    nodes.push_back(BRep_Tool::Pnt(TopoDS::Vertex(vertices(v))));
  Find find;
  find.up.resize(nodes.size());
  std::iota(find.up.begin(), find.up.end(), 0);
  // Each edge's polygons, by face, in the order faces list them.
  std::map<int, std::vector<std::vector<int>>> polygons;
  for (int i = 1; i <= faces.Extent(); ++i) {
    if (tris[i].IsNull()) continue;
    for (TopExp_Explorer ex(faces(i), TopAbs_EDGE); ex.More(); ex.Next()) {
      const TopoDS_Edge edge = TopoDS::Edge(ex.Current());
      const occ::handle<Poly_PolygonOnTriangulation>& poly =
          BRep_Tool::PolygonOnTriangulation(edge, tris[i], locs[i]);
      if (poly.IsNull()) continue;
      edge_deflection = std::max(edge_deflection, poly->Deflection());
      std::vector<int> ids;
      for (int k = 1; k <= poly->NbNodes(); ++k) ids.push_back(base[i] + poly->Node(k) - 1);
      const int e = edges.FindIndex(edge);
      TopoDS_Vertex first, last;
      TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), first, last);
      if (!first.IsNull()) find.join(ids.front(), raw + vertices.FindIndex(first) - 1);
      if (!last.IsNull()) find.join(ids.back(), raw + vertices.FindIndex(last) - 1);
      if (BRep_Tool::Degenerated(edge) && !first.IsNull())
        for (int id : ids) find.join(id, raw + vertices.FindIndex(first) - 1);
      polygons[e].push_back(ids);
    }
  }
  int length_mismatch = 0;
  for (auto& [e, list] : polygons)
    for (const auto& ids : list) {
      if (ids.size() != list.front().size()) {
        ++length_mismatch;
        continue;
      }
      for (size_t k = 0; k < ids.size(); ++k) find.join(ids[k], list.front()[k]);
    }
  double weld_gap = 0;
  for (size_t k = 0; k < nodes.size(); ++k)
    weld_gap = std::max(weld_gap, nodes[k].Distance(nodes[find.root(static_cast<int>(k))]));
  // Welded ids in first-seen order over the triangles.
  std::map<int, int> welded;
  std::vector<std::array<int, 4>> tri_list;
  int degenerate = 0;
  double face_distance = 0;
  for (int i = 1; i <= faces.Extent(); ++i) {
    if (tris[i].IsNull()) continue;
    const TopoDS_Face face = TopoDS::Face(faces(i));
    const bool reversed = face.Orientation() == TopAbs_REVERSED;
    const BRepAdaptor_Surface surface(face);
    for (int k = 1; k <= tris[i]->NbTriangles(); ++k) {
      int n[3];
      tris[i]->Triangle(k).Get(n[0], n[1], n[2]);
      if (reversed) std::swap(n[1], n[2]);
      int w[3];
      for (int j = 0; j < 3; ++j) w[j] = find.root(base[i] + n[j] - 1);
      const gp_Pnt a = nodes[base[i] + n[0] - 1], b = nodes[base[i] + n[1] - 1],
                   c = nodes[base[i] + n[2] - 1];
      for (int p = 0; p <= 4; ++p)
        for (int q = 0; q <= 4 - p; ++q) {
          const int r = 4 - p - q;
          if (p == 4 || q == 4 || r == 4) continue;
          const gp_Pnt s((p * a.X() + q * b.X() + r * c.X()) / 4,
                         (p * a.Y() + q * b.Y() + r * c.Y()) / 4,
                         (p * a.Z() + q * b.Z() + r * c.Z()) / 4);
          face_distance = std::max(face_distance, surface_distance(surface, face, s));
        }
      if (w[0] == w[1] || w[1] == w[2] || w[0] == w[2]) {
        ++degenerate;
        continue;
      }
      for (int j = 0; j < 3; ++j)
        if (!welded.count(w[j])) welded.emplace(w[j], static_cast<int>(welded.size()));
      tri_list.push_back({welded[w[0]], welded[w[1]], welded[w[2]], i});
    }
  }
  std::vector<gp_Pnt> points(welded.size());
  for (auto& [root, id] : welded) points[id] = nodes[root];
  std::map<std::pair<int, int>, int> directed, undirected;
  double area = 0, volume = 0;
  for (const auto& t : tri_list) {
    for (int j = 0; j < 3; ++j) {
      const int a = t[j], b = t[(j + 1) % 3];
      ++directed[{a, b}];
      ++undirected[{std::min(a, b), std::max(a, b)}];
    }
    const gp_XYZ p = points[t[0]].XYZ(), q = points[t[1]].XYZ(), r = points[t[2]].XYZ();
    area += (q - p).Crossed(r - p).Modulus() / 2;
    volume += p.Dot(q.Crossed(r)) / 6;
  }
  int free_edges = length_mismatch, nonmanifold = 0, misoriented = 0;
  for (const auto& [e, k] : undirected) {
    if (k == 1) ++free_edges;
    if (k > 2) ++nonmanifold;
  }
  for (const auto& [e, k] : directed)
    if (k > 1) ++misoriented;
  double edge_distance = 0;
  for (const auto& [e, list] : polygons) {
    const TopoDS_Edge edge = TopoDS::Edge(edges(e));
    if (BRep_Tool::Degenerated(edge)) continue;
    const BRepAdaptor_Curve curve(edge);
    const auto& ids = list.front();
    for (size_t k = 0; k + 1 < ids.size(); ++k)
      for (double t : {0.25, 0.5, 0.75}) {
        const gp_Pnt s(nodes[ids[k]].XYZ() * (1 - t) + nodes[ids[k + 1]].XYZ() * t);
        edge_distance = std::max(edge_distance, curve_distance(curve, edge, s));
      }
  }
  std::cout << "S " << name << " " << setting << " " << (mesher.IsDone() ? "done" : "not_done") << " "
            << faces.Extent() << " " << unmeshed << " " << raw << " " << welded.size() << " "
            << triangles << " " << degenerate << " " << free_edges << " " << nonmanifold << " "
            << misoriented << " " << weld_gap << " " << face_deflection << " " << edge_deflection
            << " " << face_distance << " " << edge_distance << " " << area << " " << volume
            << "\n";
  if (!dump) return;
  std::cout << "mesh " << name << " " << setting << " " << face_deflection << " 0\n";
  for (const gp_Pnt& p : points) std::cout << "v " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
  for (const auto& t : tri_list)
    std::cout << "t " << t[0] << " " << t[1] << " " << t[2] << " " << t[3] - 1 << " "
              << tris[t[3]]->Deflection() << " 0\n";
  for (const auto& [e, list] : polygons) {
    if (BRep_Tool::Degenerated(TopoDS::Edge(edges(e)))) continue;
    std::cout << "e " << e - 1 << " " << edge_deflection << " 0";
    for (int id : list.front()) {
      auto it = welded.find(find.root(id));
      std::cout << " " << (it == welded.end() ? -1 : it->second);
    }
    std::cout << "\n";
  }
  std::cout << "end\n";
}
}  // namespace

int main(int argc, char** argv) {
  const bool dump = argc > 1 && std::string(argv[1]) == "dump";
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepMesh_IncrementalMesh" << std::endl;
  std::cout << std::setprecision(17);
  std::map<std::string, TopoDS_Shape> bodies;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string word;
    if (!(in >> word)) continue;
    if (word == "body") {
      std::string name, text, row;
      in >> name;
      while (std::getline(std::cin, row) && row != "endbody") text += row + "\n";
      std::istringstream brep(text);
      TopoDS_Shape shape;
      BRep_Builder builder;
      try {
        BRepTools::Read(shape, brep, builder);
      } catch (const Standard_Failure&) {
        shape.Nullify();
      }
      bodies[name] = shape;
    } else if (word == "mesh") {
      std::string name, setting;
      double deflection = 0, angle = 0;
      in >> name >> setting >> deflection >> angle;
      const auto it = bodies.find(name);
      if (it == bodies.end() || it->second.IsNull()) {
        std::cout << "S " << name << " " << setting << " unreadable\n";
        continue;
      }
      try {
        mesh(name, setting, it->second, deflection, angle, dump);
      } catch (const Standard_Failure&) {
        std::cout << "S " << name << " " << setting << " failure\n";
      }
      std::cout << std::flush;
    }
  }
  return 0;
}
