#!/usr/bin/env python3
"""Regression tests for the test bridge itself, not additional upstream coverage."""
import os
from pathlib import Path
import tempfile
import unittest

import shutil

from run_upstream_tests import ROOT, build_worker, classify_missing, run_case, source_files


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
                    "dset h sqrt(4)\nbox e 2 2 h 2 2 h/2\nbfuse s m e\ncheckprops s -v 12\n", "pass")
        for wrong in ["bfuse r a b\ncheckprops r -v 64", "bcut r a b\ncheckprops r -s 64",
                      "bcommon r a b\nchecknbshapes r -compound 0", "bfuse r a b\nchecknbshapes r -solid 2"]:
            with self.subTest(wrong=wrong):
                self.expect(boxes + wrong, "failed")
        # Outside S9a, never an answer: frames with other axes (S9b), a cut
        # leaving a stack (S9a.2), prisms touching along an edge
        # (Degenerate), a cone, a section, several objects, a result of two
        # solids as an argument, unifysamedom of other shapes, counts of an
        # uncopied prism's result and a Boolean's history.
        for gap, why in [
                (boxes + "trotate b 0 0 0 0 0 1 30\ncatch {bfuse r a b}", "different axes"),
                (boxes + "box t 1 1 1 1 1 1\ncatch {bcut r a t}", "S9a.2"),
                (boxes + "box t 4 4 0 1 1 2\ncatch {bfuse r a t}", "egenerate"),
                (boxes + "pcone k 1 0 2\ncatch {bcommon r a k}", "prisms"),
                (boxes + "baddobjects a\nbaddtools b\ncatch {bapibop r 4}", "bapibop r 4"),
                (boxes + "baddobjects a b\nbaddtools b\ncatch {bapibop r 1}", "one object"),
                (boxes + "box d 9 9 9 1 1 1\nbfuse r a d\ncatch {bcut s r b}", "several solids"),
                (boxes + "catch {unifysamedom u a}", "unifysamedom u a"),
                (boxes + "polyline w 0 0 0 1 0 0 1 1 0 0 0 0\nmkplane p w\nprism q p 0 0 2\nbfuse r a q\n"
                 "checkprops r -v 32\ncatch {nbshapes r}", "uncopied"),
                (boxes + "bfuse r a b\ncatch {savehistory h}", "savehistory h")]:
            with self.subTest(why=why):
                self.assertIn(why, self.expect(gap, "unsupported")["unsupported"])

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
