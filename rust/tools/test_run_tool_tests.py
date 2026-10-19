"""The rust/tools unit tests' shards: every test in exactly one."""
from pathlib import Path
import re
import sys
import tempfile
import unittest

import run_tool_tests as rt


class ToolTestShardTests(unittest.TestCase):
    def test_plan_places_every_module_once_and_balances_the_measured_seconds(self):
        modules = sorted(rt.WEIGHTS) + ['test_new_a', 'test_new_b']
        parts = rt.plan(modules)
        self.assertEqual(parts, rt.plan(list(reversed(modules))))
        placed = [m for part in parts for m in part]
        self.assertEqual(sorted(placed), sorted(modules))
        self.assertEqual(len(placed), len(set(placed)))
        loads = [sum(rt.WEIGHTS.get(m, rt.DEFAULT_SECONDS) for m in part) for part in parts]
        self.assertLessEqual(max(loads), max(rt.WEIGHTS.values()) + sum(loads)/len(loads))
        self.assertEqual(rt.plan(['b', 'a'], 1), [['a', 'b']])

    def test_shards_together_run_exactly_the_discovered_suite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for k in range(7):
                (root/f'test_mod{k}.py').write_text(
                    'import unittest\n'
                    f'class A(unittest.TestCase):\n    def test_one(self): pass\n    def test_two(self): pass\n'
                    f'class B(unittest.TestCase):\n    def test_three(self): pass\n')
            (root/'test_broken.py').write_text('import no_such_module_anywhere\n')
            sys.path.insert(0, str(root))
            try:
                suite = unittest.TestLoader().discover(str(root), pattern='test_*.py', top_level_dir=str(root))
                every = sorted(t.id() for t in rt.cases(suite))
                self.assertEqual(len(every), 22)
                self.assertIn('test_broken', rt.modules_of(suite))
                shards = [sorted(t.id() for t in rt.cases(rt.shard_suite(suite, k, 3))) for k in range(3)]
                union = [i for s in shards for i in s]
                self.assertEqual(sorted(union), every)
                self.assertEqual(len(union), len(set(union)))
                self.assertTrue(all(shards))
            finally:
                sys.path.remove(str(root))
                for k in list(sys.modules):
                    if k.startswith('test_mod') or k == 'test_broken':
                        del sys.modules[k]

    def test_workflow_runs_every_shard_in_place_of_the_upstream_discover(self):
        text = (rt.TOOLS.parents[1]/'.github/workflows/rust-kernel.yml').read_text()
        self.assertNotIn('-m unittest discover', text)
        job = text[text.index('\n  tool-tests:\n'):]
        job = job[:re.search(r'\n  [\w-]+:\n', job[1:]).start()+1]
        shards = re.search(r'\n\s+shard: \[([^\]]*)\]', job).group(1)
        self.assertEqual([int(k) for k in shards.split(',')], list(range(rt.SHARDS)))
        self.assertIn('target/math-oracle-venv/bin/python rust/tools/run_tool_tests.py --shard ${{ matrix.shard }}', job)
        self.assertIn('RUSTY_TEST_DRAW_EXE: /usr/bin/occt-draw', job)
        upstream = text[text.index('\n  upstream-draw:\n'):text.index('\n  tool-tests:\n')]
        # The same environment as the upstream job's.
        for line in re.findall(r'\n( {6}[A-Z_]+: .*)', upstream):
            self.assertIn(line, job)
        self.assertIn(re.search(r'sudo apt-get install [^\n]*', upstream).group(0), job)


if __name__ == '__main__':
    unittest.main()
