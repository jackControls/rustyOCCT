// Source-pinned observations of IntTools_EdgeEdge (S7d of REVIEW_NOTES.md):
// each case block on stdin (`case NAME`, two `curve` rows, `end`) builds one
// edge per curve and gives `NAME done|not_done PARTS`, then per common part
// `P x y z t1 t2` (a vertex: the point on the first edge and both curves'
// parameters) or `S a b` (an edge part: its range on the first edge).
// Curves: `line x y z dx dy dz` (the line through a point along a unit
// direction, the edge over [-10, 10]), `circle|ellipse|hyperbola origin
// normal x radius [minor]` with the kernel's stored frame axes (a circle or
// an ellipse over a whole turn, a hyperbola over [-3, 3]). Nothing is
// computed from the kernel's output.
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_Line.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax2.hxx>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Row {
  std::string kind;
  std::vector<double> v;
};

TopoDS_Edge edge(const Row& c) {
  if (c.kind == "line") {
    Handle(Geom_Line) l =
        new Geom_Line(gp_Pnt(c.v[0], c.v[1], c.v[2]), gp_Dir(c.v[3], c.v[4], c.v[5]));
    return BRepBuilderAPI_MakeEdge(l, -10.0, 10.0);
  }
  const gp_Ax2 ax(gp_Pnt(c.v[0], c.v[1], c.v[2]), gp_Dir(c.v[3], c.v[4], c.v[5]),
                  gp_Dir(c.v[6], c.v[7], c.v[8]));
  if (c.kind == "ellipse") return BRepBuilderAPI_MakeEdge(new Geom_Ellipse(ax, c.v[9], c.v[10]));
  if (c.kind == "hyperbola")
    return BRepBuilderAPI_MakeEdge(new Geom_Hyperbola(ax, c.v[9], c.v[10]), -3.0, 3.0);
  return BRepBuilderAPI_MakeEdge(new Geom_Circle(ax, c.v[9]));
}

void run(const std::string& name, const Row& a, const Row& b) {
  const TopoDS_Edge e1 = edge(a), e2 = edge(b);
  IntTools_EdgeEdge inter(e1, e2);
  inter.Perform();
  std::ostringstream out;
  out << std::setprecision(17);
  if (!inter.IsDone()) {
    out << name << " not_done 0\n";
    std::cout << out.str();
    return;
  }
  const auto& parts = inter.CommonParts();
  out << name << " done " << parts.Length() << "\n";
  const BRepAdaptor_Curve curve(e1);
  for (int i = 1; i <= parts.Length(); ++i) {
    const IntTools_CommonPrt& part = parts(i);
    if (part.Type() == TopAbs_VERTEX) {
      const double t1 = part.VertexParameter1(), t2 = part.VertexParameter2();
      const gp_Pnt p = curve.Value(t1);
      out << "P " << p.X() << " " << p.Y() << " " << p.Z() << " " << t1 << " " << t2 << "\n";
    } else {
      double f, l;
      part.Range1(f, l);
      out << "S " << f << " " << l << "\n";
    }
  }
  std::cout << out.str();
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " IntTools_EdgeEdge" << std::endl;
  std::string line, name;
  std::vector<Row> rows;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string word;
    if (!(in >> word)) continue;
    if (word == "case") {
      in >> name;
      rows.clear();
    } else if (word == "curve") {
      Row r;
      in >> r.kind;
      double x;
      while (in >> x) r.v.push_back(x);
      r.v.resize(std::max<size_t>(r.v.size(), 11), 0.0);
      rows.push_back(r);
    } else if (word == "end") {
      try {
        run(name, rows.at(0), rows.at(1));
      } catch (const Standard_Failure&) {
        std::cout << name << " failure 0\n";
      }
      std::cout << std::flush;
    }
  }
  return 0;
}
