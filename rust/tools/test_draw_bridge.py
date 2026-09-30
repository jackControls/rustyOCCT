#!/usr/bin/env python3
"""Regression tests for the test bridge itself, not additional upstream coverage."""
import os
from pathlib import Path
import tempfile
import unittest

import shutil

from run_upstream_tests import ROOT, build_worker, classify_missing, release_expectation, run_case, source_files


BOX_CHECKS = """
foreach sx {-1 1} {
    foreach sy {-1 1} {
        foreach sz {-1 1} {
            box b 5 7 9 [expr $sx*2] [expr $sy*3] [expr $sz*4]
            checkshape b
            checkprops b -v 24 -deps 1e-10
            checknbshapes b -vertex 8 -edge 12 -wire 6 -face 6 -shell 1 -solid 1 -shape 34
            checkgravitycenter b -v [expr 5+$sx] [expr 7+$sy*1.5] [expr 9+$sz*2] 1e-9
        }
    }
}
box a 1 1 1
foreach gap {0 1.9e-7 2.1e-7 1} {
    box b 1+$gap 0 0 1 1 1
    set separate [regexp {NOT interfered} [isbbinterf a b]]
    checkreal separation $separate [expr $gap>2e-7] 0 0
}
copy a c
ttranslate c 3 4 5
trotate c 3 4 5 0 0 1 90
checkprops c -v 1 -deps 1e-10
checkgravitycenter c -v 2.5 4.5 5.5 1e-9
"""


class BridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = build_worker()
        cls.temp = tempfile.TemporaryDirectory(prefix="bridge-self-tests-", dir=ROOT / "target")
        cls.root = Path(cls.temp.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def execute(self, script, backend="rust", completion=True, **kwargs):
        directory = self.root / self.id().split(".")[-1] / backend
        directory.mkdir(parents=True, exist_ok=True)
        source = directory / "case.tcl"
        source.write_text(script + ('\nputs "TEST COMPLETED"\n' if completion else "\n"))
        return run_case(backend, [source], directory, worker=self.worker, **kwargs)

    def expect(self, script, status, **kwargs):
        result = self.execute(script, **kwargs)
        self.assertEqual(result["status"], status, Path(result["log"]).read_text())
        return result

    def test_original_helpers_with_real_tcl_control_flow(self):
        self.expect(BOX_CHECKS, "pass")

    @unittest.skipUnless(os.environ.get("RUSTY_TEST_DRAW_EXE"), "optional native DRAW cross-check")
    def test_box_semantics_against_native_occt(self):
        self.expect(BOX_CHECKS, "pass", backend="occt", draw_exe=os.environ["RUSTY_TEST_DRAW_EXE"])

    def test_print_only_assertion_failures(self):
        for assertion in ["checkreal wrong 2 1 0 0", "checkprops b -v 999", "checknbshapes b -face 7"]:
            with self.subTest(assertion=assertion):
                self.expect("box b 1 2 3\ncheckshape b\n" + assertion, "failed")

    def test_caught_unknown_command_cannot_pass(self):
        result = self.expect("box b 1 2 3\ncheckshape b\ncatch {not_a_command b}", "unsupported")
        self.assertIn("not_a_command", result["unsupported"])

    def test_caught_unsupported_signature_cannot_pass(self):
        self.expect("box b 1 2 3\ncheckshape b\ncatch {isbbinterf b b -o}", "unsupported")

    def test_caught_missing_fixture_cannot_pass(self):
        self.expect("box b 1 2 3\ncheckshape b\ncatch {locate_data_file absent.brep}", "missing_fixture")

    def test_geometry_error_is_a_failure(self):
        self.expect("box b 1 0 3", "failed")

    def test_constructor_without_observations_is_unverified(self):
        self.expect("box b 1 2 3", "unverified")

    def test_an_older_native_release_fails_only_as_stated(self):
        case = {"expected_occt": "viewer_skipped", "expected_rust": "unsupported",
                "expected_occt_by_version": {"7.6.3": "known_failure"}}
        old = "Open CASCADE Technology 7.6.3\nOS: Linux\n"
        result = {"status": "failed", "version": old}
        self.assertEqual(release_expectation(case, "occt", result), "known_failure")
        self.assertEqual((result["status"], result["release"]), ("known_failure", "7.6.3"))
        # The pinned release keeps its own expectation, and a failure stays one.
        result = {"status": "failed", "version": "Open CASCADE Technology 8.1.0\n"}
        self.assertEqual(release_expectation(case, "occt", result), "viewer_skipped")
        self.assertEqual(result["status"], "failed")
        # The old release passing is an outcome compared with its statement.
        result = {"status": "viewer_skipped", "version": old}
        self.assertEqual(release_expectation(case, "occt", result), "known_failure")
        self.assertEqual(result["status"], "viewer_skipped")
        # The Rust backend is never keyed by a native release.
        result = {"status": "failed", "version": old}
        self.assertEqual(release_expectation(case, "rust", result), "unsupported")
        self.assertEqual(result["status"], "failed")

    def test_original_known_failure_is_not_a_pass(self):
        self.expect('puts "TODO All: Error: wrong"\nbox b 1 2 3\ncheckshape b\ncheckreal wrong 2 1 0 0', "known_failure")

    def test_original_required_message_is_enforced(self):
        self.expect('puts "REQUIRED All: this message never arrives"\nbox b 1 2 3\ncheckshape b', "failed")

    def test_original_ignored_diagnostic_is_not_an_error(self):
        self.expect('box b 1 2 3\ncheckshape b\nputs "Relative error of mass computation : 0"', "pass")

    def test_missing_completion_marker_fails(self):
        self.expect("box b 1 2 3\ncheckshape b", "failed", completion=False)

    def test_shape_names_use_framed_utf8(self):
        self.expect('set name {box μ with spaces}\nbox $name 1 2 3\nset props [vprops $name]\nregexp {Mass : ([^\\n]+)} $props all mass\ncheckreal mass $mass 6 1e-9 0', "pass")

    def test_timeout_is_a_failure(self):
        self.expect("while {1} {}", "timeout", timeout=0.3)

    def test_missing_executable_is_a_reported_failure(self):
        self.expect("box b 1 2 3", "failed", tclsh=str(self.root / "absent-tclsh"))

    def test_crash_cannot_reuse_old_success(self):
        source = self.root / "crash.tcl"
        source.write_text("box b 1 2 3\ncheckshape b\n")
        directory = self.root / "crash"
        directory.mkdir()
        (directory / "result.txt").write_text("status 70617373\n")
        result = run_case("rust", [source], directory, worker=Path("/usr/bin/false"))
        self.assertEqual(result["status"], "failed")

    def test_history_assertions_fail_when_wrong(self):
        prism = "polyline w 0 0 0 10 0 0 10 5 0 0 5 0 0 0 0\nmkplane f w\nprism p f 0 0 3 Copy\nsavehistory h\ncheckshape p\nexplode w e\n"
        self.expect(prism + "generated g h w_1\ncheckprops g -s 30", "pass")
        for wrong in ["generated g h w_1\ncheckprops g -s 15",
                      "generated g h w_1\nchecknbshapes g -face 2",
                      'if {[isdeleted h w_1] != "Deleted."} {puts "Error: expected deleted"}',
                      "generated g h f\ncheckprops g -v 151 -deps 1e-9"]:
            with self.subTest(wrong=wrong):
                self.expect(prism + wrong, "failed")
        # History of a wire is outside the subset; so are a prism without
        # Copy's counts and history (OCCT reuses its start shapes), not its
        # properties.
        self.expect(prism + "catch {generated g h w}", "unsupported")
        self.expect(prism + "prism q f 0 0 3\ncheckprops q -v 150 -s 190 -l 144", "pass")
        self.expect(prism + "prism q f 0 0 3\ncatch {nbshapes q}", "unsupported")
        self.expect(prism + "prism q f 0 0 3\ncatch {savehistory hq}", "unsupported")

    def test_split_by_a_plane_face(self):
        # A 5 x 2 x 3 prism split at x = 2 (S8e): OCCT's two solids share the
        # cut face, its edges and vertices; areas and lengths count per use.
        prism = ("polyline w 0 0 0 5 0 0 5 2 0 0 2 0 0 0 0\nmkplane s w\nprism p s 0 0 3 Copy\n"
                 "checkshape p\nplane t 2 0 0 1 0 0\nmkface f t -10 10 -10 10\n")
        split = prism + "bclearobjects\nbcleartools\nbaddobjects p\nbaddtools f\nbfillds\n"
        self.expect(split + "bsplit r\ncheckshape r\n"
                    "checknbshapes r -vertex 12 -edge 20 -wire 11 -face 11 -shell 2 -solid 2 -compound 1\n"
                    "checkprops r -v 30 -s 74 -l 120\n", "pass")
        # The kernel's pieces have a cut face each; the result must not say so.
        for wrong in ["checknbshapes r -face 12", "checkprops r -s 68", "checkprops r -l 100"]:
            with self.subTest(wrong=wrong):
                self.expect(split + "bsplit r\n" + wrong, "failed")
        # Outside the subset: several tools, a tool face not reaching across
        # the object, a plane through a vertex, objects that may interfere,
        # a seam (ring edges), arguments changed after bfillds, and the
        # split's history.
        for gap, why in [
                (prism + "baddobjects p\nbaddtools f f\nbfillds\ncatch {bsplit r}", "one planar tool"),
                (prism + "mkface g t -1 1 -1 1\nbaddobjects p\nbaddtools g\ncatch {bapisplit r}", "does not reach"),
                (prism + "plane u 5 0 0 1 0 0\nmkface g u -9 9 -9 9\nbaddobjects p\nbaddtools g\n"
                 "catch {bapisplit r}", "through a vertex"),
                (prism + "baddobjects p p\nbaddtools f\ncatch {bapisplit r}", "interfere"),
                (prism + "pcylinder c 1 2\nbaddobjects c\nbaddtools f\ncatch {bapisplit r}", "closed edge"),
                (split + "bclearobjects\nbaddobjects p\ncatch {bsplit r}", "bsplit r"),
                (split + "bsplit r\ncatch {savehistory h}", "savehistory h")]:
            with self.subTest(why=why):
                self.assertIn(why, self.expect(gap, "unsupported")["unsupported"])

    def test_boolean_of_prisms(self):
        # Two 4 x 4 x 2 boxes overlapping in a 2 x 2 square (S9a): OCCT's
        # result is a compound; the kernel's is merged, as after unifysamedom.
        boxes = "box a 0 0 0 4 4 2\nbox b 2 2 0 4 4 2\n"
        self.expect(boxes + "bfuse r a b\ncheckshape r\ncheckprops r -v 56 -s 104\nunifysamedom u r\n"
                    "checkprops u -l 128\n"
                    "checknbshapes u -vertex 16 -edge 24 -wire 10 -face 10 -shell 1 -solid 1 -compound 1\n"
                    "bop a b\nboptuc t\ncheckprops t -v 24 -s 56\nbopcommon m\ncheckprops m -v 8\n"
                    "bclearobjects\nbcleartools\nbaddobjects a\nbaddtools b\nbfillds\nbbop c 2\n"
                    "checkprops c -v 24\nbapibop f 1\ncheckprops f -v 56\n"
                    # A result is an argument again; DRAW variables in numbers.
                    "dset h sqrt(4)\nbox e 2 2 h 2 2 h/2\nbfuse s m e\ncheckprops s -v 12\n"
                    # A whole circle of `profile`, cut from a box as a hole.
                    "profile o O 0 0 0 F 2 1 C 1 360\nprism p o 0 0 2 Copy\nbcut k a p\n"
                    "checkprops k -v [expr {32 - 2 * acos(-1)}] -deps 1e-5\n", "pass")
        for wrong in ["bfuse r a b\ncheckprops r -v 64", "bcut r a b\ncheckprops r -s 64",
                      "bcommon r a b\nchecknbshapes r -compound 0", "bfuse r a b\nchecknbshapes r -solid 2"]:
            with self.subTest(wrong=wrong):
                self.expect(boxes + wrong, "failed")
        # Outside S9a, never an answer: arcs in frames with other axes
        # (S9c.2's: two cylinders, one turned, meeting in quartics), a stack
        # as an argument whose result's top face touches
        # itself at a vertex (Degenerate: a tower's corner on the corner of
        # a notch), prisms touching along an edge (Degenerate), a torus's
        # half turn against a cylinder (S9d.4b), a section,
        # several objects, a result of two solids as an argument,
        # unifysamedom of other shapes, counts of an uncopied prism's result
        # and a Boolean's history.
        # Two cylinders, one turned, meeting in a quartic that crosses the
        # turned one's cap circle: decided since S9c.2b.2 (its values are
        # the fixture comparisons' to check).
        self.expect(boxes + "pcylinder k 1 2\ntrotate k 0 0 0 1 0 0 30\npcylinder m 1.5 3\n"
                    "bfuse r m k\n", "unverified")
        # A torus inside the box (S9d.4a): their common is the torus, 2 pi^2
        # R r^2.
        self.expect(boxes + "ptorus t 1.25 0.5\nttranslate t 2 2 1\nbcommon r a t\n"
                    "checkprops r -v [expr {0.625 * acos(-1) ** 2}] -deps 1e-6\n", "pass")
        # A half turn of it inside the box (S9d.4b.1): their common is the
        # wedge, pi^2 R r^2.
        self.expect(boxes + "ptorus t 1.25 0.5 180\nttranslate t 2 2 1\nbcommon r a t\n"
                    "checkprops r -v [expr {0.3125 * acos(-1) ** 2}] -deps 1e-6\n", "pass")
        # A coaxial pipe through a torus's tube (S9d.4b.2a): their common is
        # the tube's inner part, pi^2 / 3 - 5 sqrt(3) pi / 16.
        self.expect(boxes + "pcylinder c 1.75 2\nptorus t 2 0.5\nttranslate t 0 0 1\nbcommon r c t\n"
                    "checkprops r -v [expr {acos(-1) ** 2 / 3 - 5 * sqrt(3) * acos(-1) / 16}] "
                    "-deps 1e-6\n", "pass")
        # A coaxial pipe cut from a frustum (S9d.3b.1): the frustum less the
        # core, 14 pi / 3 - pi / 2.
        self.expect(boxes + "pcone k 2 1 2\npcylinder c 0.5 4\nttranslate c 0 0 -1\nbcut r k c\n"
                    "checkprops r -v [expr {25 * acos(-1) / 6}] -deps 1e-6\n", "pass")
        # A cone standing on the box's floor inside it (S9d.3a): its common
        # is the cone, its cut the box less it.
        self.expect(boxes + "pcone k 1 0 1.5\nttranslate k 2 2 0\nbcommon r a k\n"
                    "checkprops r -v [expr {acos(-1) / 2}] -deps 1e-6\nbcut s a k\n"
                    "checkprops s -v [expr {32 - acos(-1) / 2}] -deps 1e-6\n", "pass")
        for gap, why in [
                (boxes + "box t 1 1 2 1 1 1\nbfuse r a t\ncheckprops r -v 33\ncatch {bcut s r b}", "egenerate"),
                (boxes + "box t 4 4 0 1 1 2\ncatch {bfuse r a t}", "egenerate"),
                (boxes + "pcylinder c 1 2\nptorus t 2 0.5 180\nttranslate t 1 0 1\n"
                 "catch {bcommon r c t}", "segment or wedge"),
                (boxes + "baddobjects a\nbaddtools b\ncatch {bapibop r 4}", "bapibop r 4"),
                (boxes + "baddobjects a b\nbaddtools b\ncatch {bapibop r 1}", "one object"),
                (boxes + "box d 9 9 9 1 1 1\nbfuse r a d\ncatch {bcut s r b}", "several solids"),
                (boxes + "catch {unifysamedom u a}", "unifysamedom u a"),
                (boxes + "polyline w 0 0 0 1 0 0 1 1 0 0 0 0\nmkplane p w\nprism q p 0 0 2\nbfuse r a q\n"
                 "checkprops r -v 32\ncatch {nbshapes r}", "uncopied"),
                (boxes + "bfuse r a b\ncatch {savehistory h}", "savehistory h"),
                (boxes + "catch {dset x atan2(1,2)}", "atan2")]:
            with self.subTest(why=why):
                self.assertIn(why, self.expect(gap, "unsupported")["unsupported"])

    def test_boolean_stacks(self):
        # S9a.2's stacks: a 4 x 4 x 4 box less a 2 x 2 x 2 box inside it is
        # one solid of two shells; a 2 x 4 x 2 box on half of a 4 x 4 x 2
        # box is an L-shaped step. The kernel's stacks are unified already.
        stacks = "box o 0 0 0 4 4 4\nbox c 1 1 1 2 2 2\nbox a 0 0 0 4 4 2\nbox b 0 0 2 2 4 2\n"
        self.expect(stacks + "bcut h o c\ncheckshape h\ncheckprops h -v 56 -s 120\nunifysamedom u h\n"
                    "checkprops u -l 144\n"
                    "checknbshapes u -vertex 16 -edge 24 -wire 12 -face 12 -shell 2 -solid 1 -compound 1\n"
                    "bop a b\nbopfuse s\ncheckshape s\ncheckprops s -v 48 -s 88\nunifysamedom us s\n"
                    "checkprops us -l 112\n"
                    "checknbshapes us -vertex 12 -edge 18 -wire 8 -face 8 -shell 1 -solid 1 -compound 1\n",
                    "pass")
        # A stack's measures and counts are the kernel's, never echoed.
        for wrong in ["bcut h o c\ncheckprops h -v 64", "bcut h o c\ncheckprops h -s 96",
                      "bcut h o c\nunifysamedom u h\nchecknbshapes u -shell 1",
                      "bcut h o c\nunifysamedom u h\nchecknbshapes u -solid 2",
                      "bfuse s a b\nunifysamedom u s\ncheckprops u -l 128",
                      "bfuse s a b\nunifysamedom u s\nchecknbshapes u -face 10 -vertex 16"]:
            with self.subTest(wrong=wrong):
                self.expect(stacks + wrong, "failed")
        # A stack is a Boolean argument again (S9b.2, on its stored
        # geometry): the box's lower half fills the cavity's lower half, the
        # step lies inside the cube, and fused with it is the cube.
        for again in self.STACKS_AGAIN:
            with self.subTest(again=again):
                self.expect(self.STACKS + again, "pass")

    @unittest.skipUnless(os.environ.get("RUSTY_TEST_DRAW_EXE"), "optional native DRAW cross-check")
    def test_boolean_stacks_against_native_occt(self):
        # The same stacks' values on native DRAW, after unifysamedom, and
        # the stacks given to a Boolean again.
        self.expect("box o 0 0 0 4 4 4\nbox c 1 1 1 2 2 2\nbcut h o c\ncheckshape h\n"
                    "checkprops h -v 56 -s 120\nunifysamedom u h\ncheckprops u -l 144\n"
                    "checknbshapes u -vertex 16 -edge 24 -wire 12 -face 12 -shell 2 -solid 1 -compound 1\n",
                    "pass", backend="occt", draw_exe=os.environ["RUSTY_TEST_DRAW_EXE"])
        for again in self.STACKS_AGAIN:
            with self.subTest(again=again):
                self.expect(self.STACKS + again, "pass", backend="occt",
                            draw_exe=os.environ["RUSTY_TEST_DRAW_EXE"])

    # S9b.2: stacks given to a Boolean again.
    STACKS = "box o 0 0 0 4 4 4\nbox c 1 1 1 2 2 2\nbox a 0 0 0 4 4 2\nbox b 0 0 2 2 4 2\n"
    STACKS_AGAIN = ["bcut h o c\ncheckprops h -v 56\nbfuse r h a\ncheckprops r -v 60",
                    "bfuse s a b\ncheckprops s -v 48\nbcommon r o s\ncheckprops r -v 48",
                    "bfuse s a b\ncheckprops s -v 48\nbop s o\nbopfuse r\ncheckprops r -v 64"]
    # S9b.2: polyhedra given to a Boolean again.
    POLYHEDRA_AGAIN = ["bcut h o c\ncheckprops h -v 63\nbfuse r h a\ncheckprops r -v 63",
                       "bfuse l a b\ncheckprops l -v 4\nbop l o\nbopcut r\ncheckprops r -v 2"]

    # S9b.1's polyhedra: a 2 x 1 x 1 box and its copy turned a quarter turn
    # about z fuse into an L (the walls at y = 0 one face); a unit box
    # turned 30 degrees about (1, 1, 1) inside a 4 x 4 x 4 box leaves a
    # cavity, one solid of two shells.
    POLYHEDRA = ("box a 0 0 0 2 1 1\ntcopy a b\ntrotate b 0 0 0 0 0 1 90\n"
                 "box o 0 0 0 4 4 4\nbox c 1.5 1.5 1.5 1 1 1\ntrotate c 2 2 2 1 1 1 30\n")
    POLYHEDRA_CHECKS = ("bfuse l a b\ncheckshape l\ncheckprops l -v 4 -s 18\nunifysamedom ul l\n"
                        "checkprops ul -l 52\n"
                        "checknbshapes ul -vertex 12 -edge 18 -wire 8 -face 8 -shell 1 -solid 1 -compound 1\n"
                        "bcut h o c\ncheckshape h\ncheckprops h -v 63 -s 102\nunifysamedom uh h\n"
                        "checkprops uh -l 120\n"
                        "checknbshapes uh -vertex 16 -edge 24 -wire 12 -face 12 -shell 2 -solid 1 -compound 1\n")

    def test_boolean_polyhedra(self):
        self.expect(self.POLYHEDRA + self.POLYHEDRA_CHECKS, "pass")
        # A quarter turn about any coordinate axis, either way, turns the
        # frame exactly: each L's coplanar walls are one face.
        for box, turn in [("0 0 0 1 2 1", "1 0 0 90"), ("0 0 0 2 1 1", "0 1 0 90"),
                          ("0 0 0 2 1 1", "0 0 -1 270"), ("0 0 0 2 1 1", "0 0 1 -270")]:
            with self.subTest(turn=turn):
                self.expect(f"box p {box}\ncopy p q\ntrotate q 0 0 0 {turn}\nbfuse l p q\n"
                            "checkprops l -v 4 -s 18\nunifysamedom u l\ncheckprops u -l 52\n"
                            "checknbshapes u -vertex 12 -edge 18 -wire 8 -face 8 -shell 1 -solid 1\n", "pass")
        # A polyhedron's measures and counts are the kernel's, never echoed:
        # among them the counts of the L with its walls 6.1e-17 apart, as
        # a quarter turn rounded in radians would leave them.
        for wrong in ["bfuse l a b\ncheckprops l -v 4.5", "bfuse l a b\ncheckprops l -s 20",
                      "bfuse l a b\nunifysamedom u l\nchecknbshapes u -vertex 14 -edge 21 -face 9",
                      "bfuse l a b\nunifysamedom u l\ncheckprops u -l 54",
                      "bcut h o c\ncheckprops h -v 64", "bcut h o c\nunifysamedom u h\ncheckprops u -l 60",
                      "bcut h o c\nunifysamedom u h\nchecknbshapes u -shell 1",
                      "bcut h o c\nunifysamedom u h\nchecknbshapes u -solid 2"]:
            with self.subTest(wrong=wrong):
                self.expect(self.POLYHEDRA + wrong, "failed")
        # A polyhedron is a Boolean argument again (S9b.2): the box at the
        # cube's corner leaves the cube with its cavity, and the L less the
        # cube is its arm outside, touching the cube's face.
        for again in self.POLYHEDRA_AGAIN:
            with self.subTest(again=again):
                self.expect(self.POLYHEDRA + again, "pass")
        # A cylinder turned inside a cube cut from it (S9c.1).
        self.expect(self.POLYHEDRA + "pcylinder k 1 1\ntrotate k 0 0 0 1 1 1 30\nttranslate k 2 2 2\n"
                    "bcut r o k\ncheckprops r -v [expr {64 - acos(-1)}] -deps 1e-9\n", "pass")
        # Two cylinders apart, one turned: decided since S9c.2b.1 (refused
        # as quartics before), the cut the other cylinder whole.
        self.expect(self.POLYHEDRA + "pcylinder k 1 1\ntrotate k 2 2 2 1 1 1 30\npcylinder m 1.5 4\n"
                    "ttranslate m 2 2 0\nbcut r m k\ncheckprops r -v [expr {9 * acos(-1)}] -deps 1e-9\n", "pass")
        # Parallel cylinders, one in a turned frame, meeting in generatrices
        # at an ellipse's and a circle's crossings, and a turned pipe leaving
        # a cylinder through its top cap's rim: decided since S9c.2b.2 (the
        # fixture comparisons check their values).
        for script in ["pcylinder k 1 4\nttranslate k 1 0 0\ntrotate k 1 0 0 0 0 1 -120\n"
                       "trotate k 0 0 0 0 0 1 60\npcylinder m 1 4\nbcut r m k\n",
                       "pcylinder k 0.5 6\ntrotate k 0 0 0 1 1 1 30\nttranslate k 2.4 2.8 1\n"
                       "pcylinder m 1.5 4\nttranslate m 2 2 0\nbcut r m k\n"]:
            # Each takes 24 to 34 seconds on the debug worker of a loaded
            # machine, about the default 30 second limit.
            self.expect(self.POLYHEDRA + script, "unverified", timeout=120)
        # A turned box's corner on a face is degenerate.
        for gap, why in [("box d 0 0 0 2 2 1\ntrotate d 0 0 0 0 0 1 45\nttranslate d 2 0 0\n"
                          "box e 0 0 0 4 4 1\ncatch {bcut r e d}", "egenerate")]:
            with self.subTest(why=why):
                self.assertIn(why, self.expect(self.POLYHEDRA + gap, "unsupported")["unsupported"])

    @unittest.skipUnless(os.environ.get("RUSTY_TEST_DRAW_EXE"), "optional native DRAW cross-check")
    def test_boolean_polyhedra_against_native_occt(self):
        # The same polyhedra's values on native DRAW, after unifysamedom.
        self.expect(self.POLYHEDRA + self.POLYHEDRA_CHECKS, "pass", backend="occt",
                    draw_exe=os.environ["RUSTY_TEST_DRAW_EXE"])
        for again in self.POLYHEDRA_AGAIN:
            with self.subTest(again=again):
                self.expect(self.POLYHEDRA + again, "pass", backend="occt",
                            draw_exe=os.environ["RUSTY_TEST_DRAW_EXE"])

    def test_native_selector_picks_by_geometry(self):
        picks = self.root / "picks.txt"
        # A synthetic native record: explode 1 of a box's faces, deliberately
        # out of OCCT order, with one pick no face has; explode 2 claims edges.
        picks.write_text("pick 1 1 face 200 6 12 33\npick 1 2 face 600 1 12 18\n"
                         "pick 1 3 face 600 1 12 99\npick 2 1 edge 30 1 2 18\n")
        env = {"RUSTY_DRAW_SELECTOR": str(picks)}
        box = "box b 1 2 3 10 20 30\nexplode b f\n"
        self.expect(box + "checkprops b_1 -s 200 -deps 1e-9\ncheckgravitycenter b_1 -s 6 12 33 1e-9\n"
                    "checkprops b_2 -s 600 -deps 1e-9\ncheckgravitycenter b_2 -s 1 12 18 1e-9\n",
                    "pass", extra_env=env)
        # Wrong expectations still fail: selection is not self-confirming.
        self.expect(box + "checkgravitycenter b_2 -s 11 12 18 1e-9\n", "failed", extra_env=env)
        # A pick without an entity is lost; using it is a capability gap.
        self.expect(box + "checkprops b_1 -s 200\nsprops b_3\n", "unsupported", extra_env=env)
        # A kind that disagrees with the native record, or no record at all.
        self.expect(box + "explode b f\ncheckprops b_1 -s 200\n", "unsupported", extra_env=env)
        self.expect(box + "checkprops b_1 -s 200\n", "unsupported")

    def test_viewer_commands_are_recorded_not_run(self):
        script = ("box b 1 2 3\ncheckshape b\nvinit\nvdisplay b\nsmallview\n"
                  "checkview -display b -2d -path ${imagedir}/${test_image}.png\n")
        result = self.expect(script, "viewer_skipped")
        self.assertIn("smallview", result["viewer"])
        # A viewer command does not rescue a failed assertion.
        self.expect("box b 1 2 3\ncheckprops b -v 7\nvinit\n", "failed")

    def test_data_lookup_searches_every_subdirectory(self):
        data = self.root / "data-tree"
        nested = data / "brep" / "deeper"
        nested.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "data/occ/wedge_ok.brep", nested / "wedge_ok.brep")
        script = ("restore [locate_data_file wedge_ok.brep] w\ncheckshape w\n"
                  "checknbshapes w -vertex 8 -edge 12 -face 6 -shell 1 -solid 1\n")
        self.expect(script, "pass", data_dirs=[data])
        # A restored body has certified mass properties (S3); moving it is
        # still a capability gap.
        self.expect(script + "vprops w\n", "pass", data_dirs=[data])
        self.expect(script + "ttranslate w 1 2 3\n", "unsupported", data_dirs=[data])

    def test_restore_reports_constructs_the_kernel_cannot_represent(self):
        result = self.expect("restore [locate_data_file bottle.brep] b\ncheckshape b\n", "unsupported",
                             data_dirs=[ROOT / "data"])
        # B-splines are representable since S4; trimmed surfaces and
        # ellipses are not.
        self.assertIn("RectangularTrimmedSurface", result["unsupported"])
        self.assertIn("Ellipse2d", result["unsupported"])
        self.assertNotIn("BSpline", result["unsupported"])

    def test_missing_data_is_not_fetched_or_private(self):
        missing = {"status": "missing_fixture", "missing": "secret.brep"}
        self.assertEqual(classify_missing(dict(missing), None)["status"], "not_fetched")
        self.assertEqual(classify_missing(dict(missing), {"public.brep"})["status"], "private_data")
        self.assertEqual(classify_missing(dict(missing), {"secret.brep"})["status"], "missing_fixture")
        self.assertEqual(classify_missing({"status": "pass"}, None)["status"], "pass")

    def test_upstream_context_order(self):
        self.assertEqual([str(p.relative_to(ROOT)) for p in source_files("tests/bugs/modalg_7/bug29311_5")], [
            "tests/bugs/begin", "tests/bugs/modalg_7/begin",
            "tests/bugs/modalg_7/bug29311_5", "tests/bugs/end",
        ])
        derived = source_files("rust/fixtures/draw-derived/prism_history_rectangle", "tests/bugs/modalg_7")
        self.assertEqual([str(p.relative_to(ROOT)) for p in derived], [
            "tests/bugs/begin", "tests/bugs/modalg_7/begin",
            "rust/fixtures/draw-derived/prism_history_rectangle", "tests/bugs/end",
        ])


if __name__ == "__main__":
    unittest.main()
