// Source-pinned observations of OCCT's STEP reader (the STEP import track of
// REVIEW_NOTES.md): each input block is `case NAME LENGTH` followed by LENGTH
// bytes of a Part 21 file, read by STEPControl_Reader::ReadStream with the
// default parameters (lengths in millimetres, the file's uncertainty as
// precision, shape healing on), every root transferred. Output: `NAME
// status roots bodies` (status `done`, or the read status), then per body in
// the result, solids first and then the shells outside a solid, `B kind V E W
// F SH SO valid volume area cx cy cz tolerance`: kind `solid` or `sheet`,
// OCCT's counts of distinct sub-shapes, BRepCheck's verdict, the volume
// (`-` for a sheet), the area, the centre of the volume (of the area for a
// sheet) and the largest vertex, edge or face tolerance. Nothing is computed
// from the kernel's output.
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <NCollection_IndexedMap.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

int count(const TopoDS_Shape& s, TopAbs_ShapeEnum type) {
  ShapeMap map;
  TopExp::MapShapes(s, type, map);
  return map.Extent();
}

double tolerance(const TopoDS_Shape& s) {
  double t = 0.0;
  for (TopExp_Explorer e(s, TopAbs_VERTEX); e.More(); e.Next())
    t = std::max(t, BRep_Tool::Tolerance(TopoDS::Vertex(e.Current())));
  for (TopExp_Explorer e(s, TopAbs_EDGE); e.More(); e.Next())
    t = std::max(t, BRep_Tool::Tolerance(TopoDS::Edge(e.Current())));
  for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
    t = std::max(t, BRep_Tool::Tolerance(TopoDS::Face(e.Current())));
  return t;
}

void body(std::ostream& out, const TopoDS_Shape& s, bool solid) {
  GProp_GProps volume, area;
  BRepGProp::SurfaceProperties(s, area);
  if (solid) BRepGProp::VolumeProperties(s, volume);
  const gp_Pnt c = solid ? volume.CentreOfMass() : area.CentreOfMass();
  out << "B " << (solid ? "solid" : "sheet") << ' ' << count(s, TopAbs_VERTEX) << ' '
      << count(s, TopAbs_EDGE) << ' ' << count(s, TopAbs_WIRE) << ' ' << count(s, TopAbs_FACE)
      << ' ' << count(s, TopAbs_SHELL) << ' ' << count(s, TopAbs_SOLID) << ' '
      << BRepCheck_Analyzer(s).IsValid() << ' ';
  if (solid)
    out << volume.Mass();
  else
    out << '-';
  out << ' ' << area.Mass() << ' ' << c.X() << ' ' << c.Y() << ' ' << c.Z() << ' ' << tolerance(s)
      << "\n";
}
}  // namespace

int main() {
  std::cerr << "OCCT " << OCC_VERSION_COMPLETE << " STEPControl_Reader" << std::endl;
  // The reader's messages would interleave with the observations.
  Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::istringstream head(line);
    std::string tag, name;
    size_t length = 0;
    if (!(head >> tag >> name >> length) || tag != "case") return 2;
    std::string text(length, '\0');
    if (!std::cin.read(&text[0], static_cast<std::streamsize>(length))) return 2;
    std::ostringstream out;
    out << std::setprecision(17);
    try {
      STEPControl_Reader reader;
      std::istringstream stream(text);
      const IFSelect_ReturnStatus status = reader.ReadStream(name.c_str(), stream);
      if (status != IFSelect_RetDone) {
        std::cout << name << " read_status_" << static_cast<int>(status) << " 0 0\n" << std::flush;
        continue;
      }
      const int roots = reader.TransferRoots();
      const TopoDS_Shape shape = reader.OneShape();
      std::vector<TopoDS_Shape> solids, sheets;
      for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next()) solids.push_back(e.Current());
      for (TopExp_Explorer e(shape, TopAbs_SHELL, TopAbs_SOLID); e.More(); e.Next())
        sheets.push_back(e.Current());
      out << name << " done " << roots << ' ' << solids.size() + sheets.size() << "\n";
      for (const auto& s : solids) body(out, s, true);
      for (const auto& s : sheets) body(out, s, false);
      std::cout << out.str() << std::flush;
    } catch (const Standard_Failure& failure) {
      std::cout << name << " failure 0 0\n" << std::flush;
    }
  }
  return 0;
}
