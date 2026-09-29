// The free-shape probe of occt_free_shape_oracle.cpp with one change: faces
// and shells are integrated adaptively, BRepGProp::SurfaceProperties with
// Eps 1e-12, instead of by the default fixed Gauss rule, whose error on the
// corpus's spline faces reaches 7e-6 relative. compare_brep_io.py runs it
// only on the files holding a free face whose certified enclosure excludes
// the default measure, as reviewed evidence; the S6 capture pins the
// default probe's text and observations, which stay unchanged. Every other
// line is the default probe's (test_brep_io_oracle.py checks the pair).
#include <BRepCheck_Analyzer.hxx>
#include <cstdlib>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shape.hxx>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
void frees(const TopoDS_Shape& s, std::vector<TopoDS_Shape>& out) {
  switch (s.ShapeType()) {
    case TopAbs_COMPOUND:
      for (TopoDS_Iterator it(s); it.More(); it.Next()) frees(it.Value(), out);
      break;
    case TopAbs_SHELL:
    case TopAbs_FACE:
    case TopAbs_WIRE:
    case TopAbs_EDGE:
    case TopAbs_VERTEX:
      out.push_back(s);
      break;
    default:
      break;
  }
}

int count(const TopoDS_Shape& s, TopAbs_ShapeEnum kind) {
  TopTools_IndexedMapOfShape map;
  TopExp::MapShapes(s, kind, map);
  return map.Extent();
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepTools::Read BRepCheck_Analyzer BRepGProp adaptive 1e-12"
            << std::endl;
  std::cout << std::setprecision(17);
  std::string path;
  while (std::getline(std::cin, path)) {
    if (path.empty()) continue;
    try {
      TopoDS_Shape shape;
      BRep_Builder builder;
      if (!BRepTools::Read(shape, path.c_str(), builder) || shape.IsNull()) {
        std::cout << "F " << path << " unreadable\n";
        continue;
      }
      {
        std::vector<TopoDS_Shape> loose;
        frees(shape, loose);
        std::cout << "F " << path << " " << loose.size() << "\n";
        for (const TopoDS_Shape& s : loose) {
          BRepCheck_Analyzer analyzer(s);
          GProp_GProps g;
          const bool linear = s.ShapeType() == TopAbs_WIRE || s.ShapeType() == TopAbs_EDGE;
          if (linear)
            BRepGProp::LinearProperties(s, g);
          else if (s.ShapeType() != TopAbs_VERTEX)
            BRepGProp::SurfaceProperties(s, g, 1e-12);
          const char* type = s.ShapeType() == TopAbs_SHELL  ? "Sh"
                             : s.ShapeType() == TopAbs_FACE ? "Fa"
                             : linear && s.ShapeType() == TopAbs_WIRE ? "Wi"
                             : s.ShapeType() == TopAbs_EDGE ? "Ed"
                                                            : "Ve";
          std::cout << "X " << type << " " << (analyzer.IsValid() ? "valid" : "invalid");
          for (TopAbs_ShapeEnum kind : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE,
                                        TopAbs_SHELL, TopAbs_SOLID})
            std::cout << " " << count(s, kind);
          const gp_Pnt c = s.ShapeType() == TopAbs_VERTEX ? BRep_Tool::Pnt(TopoDS::Vertex(s)) : g.CentreOfMass();
          std::cout << " " << (s.ShapeType() == TopAbs_VERTEX ? 0.0 : g.Mass()) << " " << c.X() << " "
                    << c.Y() << " " << c.Z() << "\n";
        }
      }
    } catch (const Standard_Failure& e) {
      std::cout << "F " << path << " failure " << e.GetMessageString() << "\n";
    }
    std::cout << std::flush;
  }
  return 0;
}
