"""The fixture-generator groups' input audit and the skip it allows."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import generator_groups as gg


class GeneratorGroupTests(unittest.TestCase):
    def test_workflow_runs_every_group_through_the_audited_runner(self):
        # Text, not YAML: the oracle's venv that runs these tests has no YAML parser.
        text = (gg.ROOT/'.github/workflows/rust-kernel.yml').read_text()
        job = text[text.index('\n  fixture-generators:\n'):]
        job = job[:re.search(r'\n  [\w-]+:\n', job[1:]).start()+1]
        groups = re.search(r'\n\s+group: \[([^\]]*)\]', job).group(1)
        self.assertEqual([g.strip() for g in groups.split(',')], list(gg.GROUPS))
        self.assertNotIn('generate_', job, 'a generator runs outside generator_groups.py')
        self.assertIn('python3 rust/tools/generator_groups.py check "${{ matrix.group }}" --marker target/generator-pass/marker.json', job)
        self.assertIn('python3 rust/tools/generator_groups.py run "${{ matrix.group }}" --marker target/generator-pass/marker.json', job)
        # One group's cache prefix is never another's (given / given-curved).
        self.assertEqual(set(re.findall(r'\n\s+path: (\S+)', job)), {'target/generator-pass'})
        keys = re.findall(r'\n\s+key: (.+)', job)
        self.assertEqual(len(keys), 2)
        self.assertTrue(all(k.startswith('generator-pass-v1-${{ matrix.group }}--') for k in keys))
        self.assertEqual(re.findall(r'\n\s+restore-keys: (.+)', job), ['generator-pass-v1-${{ matrix.group }}--'])
        self.assertFalse(any('--' in g for g in gg.GROUPS))
        # Only pushes and pull requests may skip.
        self.assertIn("MAY_SKIP: ${{ github.event_name == 'push' || github.event_name == 'pull_request' }}", job)

    def test_every_command_is_a_check_of_a_tools_script(self):
        for group, commands in gg.GROUPS.items():
            self.assertTrue(commands, group)
            for interpreter, script, *args in commands:
                self.assertIn(interpreter, ('python3', 'venv'))
                self.assertTrue((gg.ROOT/script).is_file(), script)
                self.assertEqual(Path(script).parent, Path('rust/tools'))
                self.assertIn('--check', args)

    def test_every_root_a_group_names_is_keyed(self):
        # Static: the roots each group's rust/tools closure names in its path
        # constants are keyed as whole trees, or recorded per file (the
        # fixtures), or the job's own output (target/).
        known = set(gg.TREE_ROOTS) | gg.RECORDED_ROOTS
        for group, commands in gg.GROUPS.items():
            modules = sorted({m for c in commands for m in gg.tool_modules(c[1])})
            roots, _ = gg.static_audit(modules)
            self.assertLessEqual(set(roots), known, group)
            self.assertEqual(gg.key_trees(group), sorted(gg.TREE_ROOTS[r] for r in roots if r in gg.TREE_ROOTS), group)
        self.assertIn('rust/kernel', gg.key_trees('algebra'))
        self.assertIn('data', gg.key_trees('brep'))
        for rel in gg.FIXED_INPUTS:
            self.assertTrue((gg.ROOT/rel).is_file(), rel)

    def test_recorder_sees_reads_listings_existence_imports_workers_and_programs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            (root/'tools').mkdir(); (root/'fixtures').mkdir(); (root/'audit').mkdir(); (root/'listed').mkdir()
            (root/'fixtures/read.txt').write_text('a')
            (root/'fixtures/worker.txt').write_text('b')
            (root/'listed/x').write_text('')
            (root/'tools/helper.py').write_text('VALUE = 1\n')
            (root/'tools/main.py').write_text(
                'import multiprocessing, os, subprocess, sys\n'
                'from pathlib import Path\n'
                'import helper\n'
                'ROOT = Path(__file__).resolve().parents[1]\n'
                'def work(_):\n'
                '    return (ROOT/"fixtures/worker.txt").read_text()\n'
                'if __name__ == "__main__":\n'
                '    (ROOT/"fixtures/read.txt").read_text()\n'
                '    os.listdir(ROOT/"listed")\n'
                '    (ROOT/"fixtures/absent.txt").exists()\n'
                '    with open(ROOT/"fixtures/written.txt", "w") as f: f.write("w")\n'
                '    method = "fork" if "fork" in multiprocessing.get_all_start_methods() else "spawn"\n'
                '    for m in sorted({method, "spawn"}):\n'
                '        with multiprocessing.get_context(m).Pool(1) as pool: pool.map(work, [0])\n'
                '    subprocess.run(["/bin/echo"], stdout=subprocess.DEVNULL)\n')
            env = dict(os.environ, GENERATOR_AUDIT_DIR=str(root/'audit'), PYTHONPATH=str(gg.HOOK))
            subprocess.run([sys.executable, str(root/'tools/main.py')], env=env, check=True, cwd=root)
            inputs, programs = gg.read_audit(root/'audit', root)
            for rel in ['tools/main.py', 'tools/helper.py', 'fixtures/read.txt', 'fixtures/worker.txt', 'fixtures/absent.txt']:
                self.assertIn(rel, inputs)
            self.assertTrue(inputs['listed'])
            self.assertNotIn('fixtures/written.txt', inputs)
            self.assertFalse([r for r in inputs if r.startswith('audit')])
            self.assertIn('/bin/echo', programs)
            self.assertFalse(all(gg.is_python(p) for p in programs))
            self.assertTrue(gg.is_python('/usr/bin/python3.12') and gg.is_python('python'))

    def test_a_marker_allows_a_skip_only_on_identical_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            for rel in gg.FIXED_INPUTS:
                (root/rel).parent.mkdir(parents=True, exist_ok=True)
                (root/rel).write_text(rel)
            (root/'rust/fixtures').mkdir(parents=True)
            (root/'rust/fixtures/a.tsv').write_text('1')
            (root/'rust/fixtures/listed').mkdir()
            (root/'rust/fixtures/listed/one').write_text('')
            (root/'rust/kernel').mkdir(parents=True)
            (root/'rust/kernel/lib.rs').write_text('fn f() {}')
            audit = root/'audit'; audit.mkdir()
            (audit/'1.log').write_text(f'open\t{root}/rust/fixtures/a.tsv\nlist\t{root}/rust/fixtures/listed\n'
                                       f'stat\t{root}/rust/fixtures/absent.tsv\nopen\t{root}/target/venv.py\n'
                                       f'open\t{root}/rust/tools/__pycache__/generator_groups.cpython-312.pyc\n'
                                       'open\t/usr/lib/python3/os.py\nexec\t/usr/bin/python3\n')
            groups = {'demo': [('venv', 'rust/tools/generate_demo.py', '--check')]}
            with patch.object(gg, 'GROUPS', groups), patch.object(gg, 'key_trees', return_value=['rust/kernel']):
                marker = gg.build_marker('demo', audit, root, passed_in={'run_id': '42'})
                marker = json.loads(json.dumps(marker))
                self.assertEqual(sorted(marker['inputs']), sorted(gg.FIXED_INPUTS + [
                    'rust/fixtures/a.tsv', 'rust/fixtures/listed', 'rust/fixtures/absent.tsv']))
                self.assertTrue(marker['memoisable'])
                self.assertEqual(gg.check_marker('demo', marker, root), [])

                def reasons(change, undo):
                    change()
                    try:
                        return gg.check_marker('demo', marker, root)
                    finally:
                        undo()
                path = root/'rust/fixtures/a.tsv'
                self.assertTrue(reasons(lambda: path.write_text('2'), lambda: path.write_text('1')))
                new = root/'rust/fixtures/listed/two'
                self.assertTrue(reasons(lambda: new.write_text(''), new.unlink))
                absent = root/'rust/fixtures/absent.tsv'
                self.assertTrue(reasons(lambda: absent.write_text(''), absent.unlink))
                workflow = root/gg.FIXED_INPUTS[0]
                self.assertTrue(reasons(lambda: workflow.write_text('changed'), lambda: workflow.write_text(gg.FIXED_INPUTS[0])))
                kernel = root/'rust/kernel/lib.rs'
                self.assertTrue(reasons(lambda: kernel.write_text('fn g() {}'), lambda: kernel.write_text('fn f() {}')))
                # A file the group never read, outside its keyed trees, changes nothing.
                other = root/'rust/fixtures/b.tsv'
                self.assertEqual(reasons(lambda: other.write_text('x'), other.unlink), [])
                self.assertTrue(gg.check_marker('demo', {**marker, 'environment': {**marker['environment'], 'image': 'x'}}, root))
                self.assertTrue(gg.check_marker('demo', {**marker, 'memoisable': False}, root))
                self.assertTrue(gg.check_marker('demo', {**marker, 'commands': []}, root))
                self.assertTrue(gg.check_marker('other', marker, root))
            with patch.object(gg, 'GROUPS', groups), patch.object(gg, 'key_trees', return_value=['rust/kernel']):
                (audit/'2.log').write_text('exec\tcargo\n')
                self.assertFalse(gg.build_marker('demo', audit, root)['memoisable'])


if __name__ == '__main__':
    unittest.main()
