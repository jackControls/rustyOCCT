// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
// Optional validation executable. Never linked into the Rust kernel.
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T> T token() {
  T value{};
  if (!(std::cin >> value)) throw std::runtime_error("invalid fixture token");
  return value;
}
static gp_Pnt point() {
  const double x = token<double>(), y = token<double>(), z = token<double>();
  return gp_Pnt(x, y, z);
}
static gp_Pnt placed(const gp_Ax2& frame, double x, double y, double z) {
  return frame.Location().Translated(gp_Vec(frame.XDirection()) * x +
                                    gp_Vec(frame.YDirection()) * y +
                                    gp_Vec(frame.Direction()) * z);
}
static TopoDS_Wire boundary(const gp_Ax2& frame, double height, bool inner) {
  const std::string kind = token<std::string>();
  TopoDS_Wire wire;
  if (kind == "P") {
    const int count = token<int>();
    std::vector<std::array<double, 2>> points;
    for (int i = 0; i < count; ++i) points.push_back({token<double>(), token<double>()});
    if (points.size() > 1 && points.front() == points.back()) points.pop_back();
    if (points.size() < 3) throw std::runtime_error("invalid polygon fixture");
    double area = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
      const auto& a = points[i];
      const auto& b = points[(i + 1) % points.size()];
      area += (a[0] - points[0][0]) * (b[1] - points[0][1]) -
              (a[1] - points[0][1]) * (b[0] - points[0][0]);
    }
    if (area < 0.0) std::reverse(points.begin() + 1, points.end());
    BRepBuilderAPI_MakePolygon polygon;
    for (const auto& p : points) polygon.Add(placed(frame, p[0], p[1], height));
    polygon.Close();
    if (!polygon.IsDone()) throw std::runtime_error("OCCT wire failed");
    wire = polygon.Wire();
  } else if (kind == "C") {
    const double x = token<double>(), y = token<double>(), radius = token<double>();
    const gp_Circ circle(gp_Ax2(placed(frame, x, y, height), frame.Direction(), frame.XDirection()), radius);
    wire = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circle).Edge()).Wire();
  } else {
    throw std::runtime_error("unknown fixture boundary kind");
  }
  if (inner) wire.Reverse();
  return wire;
}

int main() {
  try {
    std::cerr << "OCCT " << OCC_VERSION_COMPLETE << '\n';
    std::cout << std::setprecision(17) << std::scientific;
    const int count = token<int>();
    for (int c = 0; c < count; ++c) {
      const std::string name = token<std::string>();
      const int boundaries = token<int>();
      const gp_Pnt origin = point();
      const gp_Pnt normal = point();
      const gp_Pnt hint = point();
      const gp_Ax2 frame(origin, gp_Dir(normal.X(), normal.Y(), normal.Z()), gp_Dir(hint.X(), hint.Y(), hint.Z()));
      const double start = token<double>(), end = token<double>();
      const double low = std::min(start, end), high = std::max(start, end);
      const gp_Pln plane(placed(frame, 0, 0, low), frame.Direction());
      BRepBuilderAPI_MakeFace face(plane, boundary(frame, low, false), true);
      for (int i = 1; i < boundaries; ++i) face.Add(boundary(frame, low, true));
      if (!face.IsDone()) throw std::runtime_error("OCCT face failed: " + name);
      BRepPrimAPI_MakePrism prism(face.Face(), gp_Vec(frame.Direction()) * (high - low));
      if (!prism.IsDone()) throw std::runtime_error("OCCT prism failed: " + name);
      const TopoDS_Shape solid = prism.Shape();
      if (!BRepCheck_Analyzer(solid).IsValid()) throw std::runtime_error("OCCT invalid solid: " + name);
      GProp_GProps volume, area;
      BRepGProp::VolumeProperties(solid, volume);
      BRepGProp::SurfaceProperties(solid, area);
      const gp_Pnt center = volume.CentreOfMass();
      Bnd_Box bounds;
      BRepBndLib::AddOptimal(solid, bounds, false, false);
      double xmin, ymin, zmin, xmax, ymax, zmax;
      bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
      std::cout << name << ' ' << volume.Mass() << ' ' << area.Mass();
      for (double value : {center.X(), center.Y(), center.Z(), xmin, ymin, zmin, xmax, ymax, zmax}) std::cout << ' ' << value;
      // OCCT's global-origin inertia accumulation loses low bits for small
      // parts far from the origin. Recenter a geometry copy for the central
      // tensor; keep the original shape for centroid, bounds and queries.
      gp_Trsf centering;
      centering.SetTranslation(gp_Vec(center, gp_Pnt(0, 0, 0)));
      const TopoDS_Shape centered = BRepBuilderAPI_Transform(solid, centering, true).Shape();
      GProp_GProps centeredVolume;
      BRepGProp::VolumeProperties(centered, centeredVolume);
      const gp_Mat inertia = centeredVolume.MatrixOfInertia();
      for (int i = 1; i <= 3; ++i) for (int j = 1; j <= 3; ++j) std::cout << ' ' << inertia.Value(i, j);
      for (auto kind : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE}) {
        TopTools_IndexedMapOfShape shapes;
        TopExp::MapShapes(solid, kind, shapes);
        std::cout << ' ' << shapes.Extent();
      }
      const int queries = token<int>();
      for (int i = 0; i < queries; ++i) {
        const gp_Pnt p = point();
        BRepClass3d_SolidClassifier classifier(solid, placed(frame, p.X(), p.Y(), p.Z()), 1e-7);
        const auto state = classifier.State();
        if (state != TopAbs_OUT && state != TopAbs_ON && state != TopAbs_IN) throw std::runtime_error("OCCT unknown point state");
        std::cout << ' ' << (state == TopAbs_OUT ? 0 : state == TopAbs_ON ? 1 : 2);
      }
      std::cout << '\n';
    }
    std::string trailing;
    if (std::cin >> trailing) throw std::runtime_error("unexpected fixture tokens");
  } catch (const Standard_Failure& failure) {
    std::cerr << failure.GetMessageString() << '\n';
    return 1;
  } catch (const std::exception& failure) {
    std::cerr << failure.what() << '\n';
    return 1;
  }
  return 0;
}
