// Source-pinned observations of .brep files (T2): each input line is a path
// read with BRepTools::Read. Solids are the ones reached from the root through
// compounds (TopoDS_Iterator composes locations and orientations), in order.
// Per solid: the BRepCheck_Analyzer verdict, the distinct subshape counts
// (TopExp::MapShapes, same TShape and location), and BRepGProp volume,
// surface area and centroid. Nothing is computed from the kernel's output
// beyond reading it.
#include <BRepCheck_Analyzer.hxx>
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
void solids(const TopoDS_Shape& s, std::vector<TopoDS_Shape>& out) {
  if (s.ShapeType() == TopAbs_SOLID) {
    out.push_back(s);
  } else if (s.ShapeType() == TopAbs_COMPOUND) {
    for (TopoDS_Iterator it(s); it.More(); it.Next()) solids(it.Value(), out);
  }
}

int count(const TopoDS_Shape& s, TopAbs_ShapeEnum kind) {
  TopTools_IndexedMapOfShape map;
  TopExp::MapShapes(s, kind, map);
  return map.Extent();
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " BRepTools::Read BRepCheck_Analyzer BRepGProp"
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
      std::vector<TopoDS_Shape> found;
      solids(shape, found);
      std::cout << "F " << path << " " << found.size() << "\n";
      for (const TopoDS_Shape& s : found) {
        BRepCheck_Analyzer analyzer(s);
        GProp_GProps volume, surface;
        BRepGProp::VolumeProperties(s, volume);
        BRepGProp::SurfaceProperties(s, surface);
        const gp_Pnt c = volume.CentreOfMass();
        std::cout << "S " << (analyzer.IsValid() ? "valid" : "invalid");
        for (TopAbs_ShapeEnum kind : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE,
                                      TopAbs_SHELL, TopAbs_SOLID})
          std::cout << " " << count(s, kind);
        std::cout << " " << volume.Mass() << " " << surface.Mass() << " " << c.X() << " "
                  << c.Y() << " " << c.Z() << "\n";
      }
    } catch (const Standard_Failure& e) {
      std::cout << "F " << path << " failure " << e.GetMessageString() << "\n";
    }
    std::cout << std::flush;
  }
  return 0;
}
