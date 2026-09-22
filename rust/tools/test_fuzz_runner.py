"""Exercise campaign accounting and process containment, without a nightly build."""
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
import run_fuzz


class FuzzRunnerTests(unittest.TestCase):
    def test_requires_completed_mutation_after_corpus_replay(self):
        text = '#99\tINITED cov: 12 ft: 50\n#102\tDONE cov: 14\nstat::number_of_executed_units: 102\n'
        self.assertEqual(run_fuzz.statistics(text),{
            'executions':102,'initial_executions':99,'mutation_executions':3,'coverage_edges':14})
        self.assertIsNone(run_fuzz.statistics('#99 INITED cov: 12')['mutation_executions'])
        self.assertEqual(run_fuzz.statistics('#99 INITED cov: 12\nstat::number_of_executed_units: 99')['mutation_executions'],0)

    def test_propagates_failure_and_preserves_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log_path = Path(directory)/'run.log'
            with log_path.open('w') as log:
                code = run_fuzz.run_process([sys.executable,'-c','print("failure evidence",flush=True); raise SystemExit(7)'],log,3,os.environ)
            self.assertEqual(code,7)
            self.assertIn('failure evidence',log_path.read_text())

    def test_mutation_budget_starts_after_delayed_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            script=('import time,pathlib,sys; time.sleep(.25); '
                    'print("#99 INITED cov: 12",flush=True); '
                    '\nwhile not pathlib.Path(sys.argv[1]).exists(): time.sleep(.01)\n'
                    'print("stat::number_of_executed_units: 102",flush=True)')
            timer=run_fuzz.MutationBudget(log_path,stop,.2)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c',script,str(stop)],log,3,os.environ,timer.tick)
            self.assertEqual(code,0)
            self.assertGreaterEqual(timer.initialized-timer.started,.25)
            self.assertGreaterEqual(timer.stopped-timer.initialized,.2)
            self.assertTrue(stop.read_bytes(),'libFuzzer ignores empty stop files')
            self.assertTrue(timer.evidence()['mutation_budget_completed'])

    def test_early_exit_does_not_complete_mutation_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            timer=run_fuzz.MutationBudget(log_path,stop,1.)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c','print("#99 INITED cov: 12",flush=True)'],log,3,os.environ,timer.tick)
            self.assertEqual(code,0)
            self.assertFalse(timer.evidence()['mutation_budget_completed'])
            self.assertFalse(stop.exists())

    def test_outer_deadline_still_applies_after_stop_request(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            timer=run_fuzz.MutationBudget(log_path,stop,.1)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c','import time; print("#99 INITED cov: 12",flush=True); time.sleep(10)'],log,.5,os.environ,timer.tick)
            self.assertEqual(code,124)
            self.assertTrue(timer.evidence()['mutation_budget_completed'])

    @unittest.skipUnless(os.name == 'posix','campaign runner requires Unix process groups')
    def test_deadline_kills_fuzzer_descendants(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            marker = directory/'orphan-survived'
            ready = directory/'started'
            child = 'import time,pathlib,sys; time.sleep(1); pathlib.Path(sys.argv[1]).touch()'
            parent = ('import subprocess,sys,time,pathlib; '
                      'subprocess.Popen([sys.executable,"-c",sys.argv[1],sys.argv[2]]); '
                      'pathlib.Path(sys.argv[3]).touch(); time.sleep(10)')
            with (directory/'log').open('w') as log:
                code = run_fuzz.run_process([sys.executable,'-c',parent,child,str(marker),str(ready)],log,0.4,os.environ)
            self.assertEqual(code,124)
            self.assertTrue(ready.exists())
            time.sleep(1.)
            self.assertFalse(marker.exists(),'timed-out campaign left a live descendant')


if __name__ == '__main__':
    unittest.main()
