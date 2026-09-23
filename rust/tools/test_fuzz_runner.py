"""Exercise campaign accounting and process containment, without a nightly build."""
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import run_fuzz


class FuzzRunnerTests(unittest.TestCase):
    def test_retained_corpus_startup_scales_and_remains_capped(self):
        self.assertEqual(run_fuzz.startup_budget(90),600)
        self.assertEqual(run_fuzz.startup_budget(648),1416)
        self.assertEqual(run_fuzz.startup_budget(100000),3600)
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); log_path=directory/'log'
            log_path.write_text('#649 INITED cov: 12\n')
            with patch('run_fuzz.time.monotonic',return_value=0.):
                timer=run_fuzz.MutationBudget(log_path,directory/'stop',60,startup_seconds=run_fuzz.startup_budget(648))
            with patch('run_fuzz.time.monotonic',return_value=1416.01):
                timer.tick()
                self.assertFalse(timer.evidence()['startup_budget_completed'])
            self.assertEqual(timer.deadline(),1416)

    def test_tensor_quarantine_budget_is_explicit_and_target_local(self):
        base={'ASAN_OPTIONS':'detect_stack_use_after_return=1','CARGO_NET_OFFLINE':'true'}
        env=run_fuzz.campaign_environment('surface_knots',base)
        self.assertEqual(env['ASAN_OPTIONS'],'detect_stack_use_after_return=1:quarantine_size_mb=64')
        self.assertEqual(env['CARGO_NET_OFFLINE'],'true')
        self.assertEqual(base['ASAN_OPTIONS'],'detect_stack_use_after_return=1')
        for target in run_fuzz.TARGETS:
            if target!='surface_knots':
                self.assertEqual(run_fuzz.campaign_environment(target,base),base)
        self.assertEqual(run_fuzz.campaign_environment('surface_knots',{})['ASAN_OPTIONS'],'quarantine_size_mb=64')

    def test_allocator_hook_is_linked_only_with_address_sanitizer(self):
        for target in run_fuzz.TARGETS:
            args=run_fuzz.sanitizer_build_args(target)
            self.assertEqual(args[:2],['--sanitizer','address'])
            if target=='surface_knots':
                self.assertEqual(args[2:],['--features','asan-allocator'])
            else:
                self.assertEqual(args[2:],[])

    def test_requires_completed_mutation_after_corpus_replay(self):
        text = '#99\tINITED cov: 12 ft: 50\n#102\tDONE cov: 14\nstat::number_of_executed_units: 102\n'
        self.assertEqual(run_fuzz.statistics(text),{
            'executions':102,'initial_executions':99,'mutation_executions':3,'coverage_edges':14,
            'slowest_input_seconds':None,'peak_rss_mb':None})
        self.assertIsNone(run_fuzz.statistics('#99 INITED cov: 12')['mutation_executions'])
        self.assertEqual(run_fuzz.statistics('#99 INITED cov: 12\nstat::number_of_executed_units: 99')['mutation_executions'],0)

    def test_successful_exit_does_not_exempt_reported_resource_overrun(self):
        # The first Linux surface-knot campaign exited zero but reported a 64s
        # input against a 60s alarm. Such a run must fail the resource gate.
        stats=run_fuzz.statistics('stat::slowest_unit_time_sec: 64\nstat::peak_rss_mb: 1400\n')
        report={**stats,'input_limit_seconds':60,'exit_code':0}
        self.assertFalse(run_fuzz.resource_limits_satisfied(report))
        report['slowest_input_seconds']=60
        self.assertTrue(run_fuzz.resource_limits_satisfied(report))
        report['peak_rss_mb']=2049
        self.assertFalse(run_fuzz.resource_limits_satisfied(report))
        for key in ['slowest_input_seconds','peak_rss_mb']:
            self.assertFalse(run_fuzz.resource_limits_satisfied({**report,key:None}))

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
            timer=run_fuzz.MutationBudget(log_path,stop,.1,startup_seconds=1,shutdown_seconds=.3)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c','import time; print("#99 INITED cov: 12",flush=True); time.sleep(10)'],log,3,os.environ,timer.tick,timer.deadline)
            self.assertEqual(code,124)
            self.assertTrue(timer.evidence()['mutation_budget_completed'])

    def test_late_replay_cannot_borrow_mutation_or_shutdown_time(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            log_path.write_text('')
            with patch('run_fuzz.time.monotonic',return_value=0.):
                timer=run_fuzz.MutationBudget(log_path,stop,60)
            self.assertEqual(timer.deadline(),600)
            log_path.write_text('#99 INITED cov: 12\n')
            with patch('run_fuzz.time.monotonic',return_value=600.01):
                timer.tick()
                self.assertFalse(timer.evidence()['startup_budget_completed'])
            self.assertEqual(timer.deadline(),600)

    def test_linux_boundary_keeps_startup_cap_and_final_input_grace(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            log_path.write_text('#332 INITED cov: 5689\n')
            with patch('run_fuzz.time.monotonic',return_value=0.):
                timer=run_fuzz.MutationBudget(log_path,stop,60)
            with patch('run_fuzz.time.monotonic',return_value=593.19): timer.tick()
            with patch('run_fuzz.time.monotonic',return_value=653.28):
                timer.tick()
                evidence=timer.evidence()
            self.assertTrue(evidence['startup_budget_completed'])
            self.assertTrue(evidence['mutation_budget_completed'])
            self.assertAlmostEqual(timer.deadline(),678.19)
            self.assertEqual(evidence['shutdown_grace_seconds'],25)

    def test_tensor_budget_allows_complete_final_input_without_extending_startup(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); log_path,stop=directory/'log',directory/'stop'
            log_path.write_text('#99 INITED cov: 12\n')
            input_seconds=run_fuzz.TARGET_INPUT_SECONDS['surface_knots']
            with patch('run_fuzz.time.monotonic',return_value=0.):
                timer=run_fuzz.MutationBudget(log_path,stop,60,shutdown_seconds=input_seconds+5)
            with patch('run_fuzz.time.monotonic',return_value=590.): timer.tick()
            with patch('run_fuzz.time.monotonic',return_value=650.): timer.tick()
            self.assertEqual(timer.deadline(),715.)
            self.assertEqual(timer.evidence()['startup_limit_seconds'],600)
            self.assertEqual(timer.evidence()['shutdown_grace_seconds'],65)
            self.assertTrue(stop.read_bytes())
            with patch('run_fuzz.time.monotonic',return_value=0.):
                late=run_fuzz.MutationBudget(log_path,directory/'late-stop',60,shutdown_seconds=input_seconds+5)
            with patch('run_fuzz.time.monotonic',return_value=600.01): late.tick()
            self.assertEqual(late.deadline(),600.)
            self.assertFalse(late.evidence()['startup_budget_completed'])

    def test_missing_initialization_ends_at_the_startup_deadline(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            timer=run_fuzz.MutationBudget(log_path,stop,1.,startup_seconds=.3,shutdown_seconds=1.)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c','import time; time.sleep(10)'],log,3,os.environ,timer.tick,timer.deadline)
            self.assertEqual(code,124)
            self.assertLess(time.monotonic()-timer.started,2.)
            self.assertFalse(timer.evidence()['startup_budget_completed'])
            self.assertFalse(timer.evidence()['mutation_budget_completed'])

    def test_inflight_input_finishes_after_the_mutation_stop(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            log_path,stop=directory/'log',directory/'stop'
            script=('import time,pathlib,sys; time.sleep(.65); '
                    'print("#99 INITED cov: 12",flush=True); '
                    '\nwhile not pathlib.Path(sys.argv[1]).exists(): time.sleep(.01)\n'
                    'time.sleep(.65); print("stat::number_of_executed_units: 102",flush=True)')
            timer=run_fuzz.MutationBudget(log_path,stop,.2,startup_seconds=1.,shutdown_seconds=1.)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c',script,str(stop)],log,2.2,os.environ,timer.tick,timer.deadline)
            self.assertEqual(code,0)
            self.assertTrue(timer.evidence()['startup_budget_completed'])
            self.assertTrue(timer.evidence()['mutation_budget_completed'])
            self.assertEqual(run_fuzz.statistics(log_path.read_text())['mutation_executions'],3)

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
