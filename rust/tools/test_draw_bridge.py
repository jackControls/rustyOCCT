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
        # History of a wire, and a prism without Copy, are outside the subset.
        self.expect(prism + "catch {generated g h w}", "unsupported")
        self.expect(prism + "catch {prism q f 0 0 3}", "unsupported")

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
        # Anything else a restored body cannot answer yet is a capability gap.
        self.expect(script + "vprops w\n", "unsupported", data_dirs=[data])

    def test_restore_reports_constructs_the_kernel_cannot_represent(self):
        result = self.expect("restore [locate_data_file Ball.brep] b\ncheckshape b\n", "unsupported",
                             data_dirs=[ROOT / "data"])
        self.assertIn("SphericalSurface", result["unsupported"])

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
