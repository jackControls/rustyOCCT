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
        self.assertEqual(run_fuzz.startup_budget(0),620)
        self.assertEqual(run_fuzz.startup_budget(90),2420)
        self.assertEqual(run_fuzz.startup_budget(648),3600)
        self.assertEqual(run_fuzz.startup_budget(100000),3600)
        self.assertEqual(run_fuzz.startup_budget(10,60),1260)
        self.assertEqual(run_fuzz.startup_budget(385,60),3600)
        # surface_knots' larger cap (REVIEW_NOTES R12).
        cap=run_fuzz.TARGET_MAX_STARTUP_SECONDS['surface_knots']
        self.assertEqual(run_fuzz.startup_budget(577,60,cap),7200)
        self.assertEqual(run_fuzz.startup_budget(100,60,cap),6660)
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); log_path=directory/'log'
            log_path.write_text('#649 INITED cov: 12\n')
            with patch('run_fuzz.time.monotonic',return_value=0.):
                timer=run_fuzz.MutationBudget(log_path,directory/'stop',60,startup_seconds=run_fuzz.startup_budget(648))
            with patch('run_fuzz.time.monotonic',return_value=3600.01):
                timer.tick()
                self.assertFalse(timer.evidence()['startup_budget_completed'])
            self.assertEqual(timer.deadline(),3600)

    def test_tensor_quarantine_budget_is_explicit_and_target_local(self):
        base={'ASAN_OPTIONS':'detect_stack_use_after_return=1','CARGO_NET_OFFLINE':'true'}
        env=run_fuzz.campaign_environment('surface_knots',base)
        self.assertEqual(env['ASAN_OPTIONS'],'detect_stack_use_after_return=1:quarantine_size_mb=64')
        self.assertEqual(env['CARGO_NET_OFFLINE'],'true')
        self.assertEqual(base['ASAN_OPTIONS'],'detect_stack_use_after_return=1')
        for target in run_fuzz.TARGETS:
            if target not in run_fuzz.ALLOCATOR_TARGETS:
                self.assertEqual(run_fuzz.campaign_environment(target,base),base)
        self.assertEqual(run_fuzz.campaign_environment('surface_knots',{})['ASAN_OPTIONS'],'quarantine_size_mb=64')
        # The split target's spiric sections fill the stack depot (37008675181).
        self.assertEqual(run_fuzz.campaign_environment('split',{})['ASAN_OPTIONS'],
                         'quarantine_size_mb=64:malloc_context_size=5')

    def test_allocator_hook_is_linked_only_with_address_sanitizer(self):
        for target in run_fuzz.TARGETS:
            args=run_fuzz.sanitizer_build_args(target)
            self.assertEqual(args[:2],['--sanitizer','address'])
            if target in run_fuzz.ALLOCATOR_TARGETS:
                self.assertEqual(args[2:],['--features','asan-allocator'])
            else:
                self.assertEqual(args[2:],[])

    def test_minimization_merges_with_each_targets_limits(self):
        command=run_fuzz.merge_command('surface_knots','nightly',Path('new'),Path('old'),Path('art'))
        self.assertEqual(command[:7],['cargo','+nightly','fuzz','run','surface_knots','new','old'])
        tail=command[command.index('--')+1:]
        self.assertIn('-merge=1',tail)
        self.assertIn('-timeout=60',tail)
        self.assertIn('-rss_limit_mb=2048',tail)
        self.assertIn('-max_len=4096',tail)
        self.assertIn('--features',command)
        self.assertIn('-timeout=20',run_fuzz.merge_command('brep_io','n',Path('a'),Path('b'),Path('c')))

    def test_failed_minimization_keeps_the_corpus(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); corpus=directory/'corpus'/'brep_io'; corpus.mkdir(parents=True)
            for k in range(3): (corpus/f'input{k}').write_bytes(bytes([k]))
            with patch('run_fuzz.FUZZ',directory), patch('run_fuzz.run_process',return_value=1):
                report=run_fuzz.minimize('brep_io','n',corpus,directory,{})
            self.assertFalse(report['replaced'])
            self.assertEqual(sorted(p.name for p in corpus.iterdir()),['input0','input1','input2'])
            self.assertFalse((corpus.parent/'.brep_io-minimized').exists())

            def merged(command,log,timeout,env):
                (Path(command[5])/'input1').write_bytes(b'\x01')
                return 0
            with patch('run_fuzz.FUZZ',directory), patch('run_fuzz.run_process',side_effect=merged):
                report=run_fuzz.minimize('brep_io','n',corpus,directory,{})
            self.assertTrue(report['replaced'])
            self.assertEqual((report['corpus_files_before'],report['corpus_files_after']),(3,1))
            self.assertEqual([p.name for p in corpus.iterdir()],['input1'])

    def test_push_replay_takes_regressions_new_inputs_and_a_seeded_sample(self):
        import hashlib
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            corpus=directory/'corpus'/'brep_io'; corpus.mkdir(parents=True)
            regressions=directory/'regressions'/'brep_io'; regressions.mkdir(parents=True)
            (regressions/'r.bin').write_bytes(b'regression')
            names=[f'{k:03}' for k in range(200)]
            for name in names: (corpus/name).write_bytes(name.encode())
            regression=hashlib.sha256(b'regression').hexdigest()
            (corpus/regression).write_bytes(b'regression')
            with patch('run_fuzz.FUZZ',directory):
                # No manifest: everything is new, so the whole corpus replays.
                chosen,evidence=run_fuzz.replay_plan('brep_io',corpus,7)
                self.assertEqual(len(chosen),201)
                self.assertFalse(evidence['manifest'])
                run_fuzz.write_manifest('brep_io',corpus)
                for name in ['new1','new2']: (corpus/name).write_bytes(name.encode())
                chosen,evidence=run_fuzz.replay_plan('brep_io',corpus,7)
                again,_=run_fuzz.replay_plan('brep_io',corpus,7)
                other,_=run_fuzz.replay_plan('brep_io',corpus,8)
            self.assertEqual(chosen,again)
            self.assertNotEqual(chosen,other)
            self.assertIn(regression,chosen)
            self.assertTrue({'new1','new2'} <= set(chosen))
            self.assertEqual(len(chosen),1+2+run_fuzz.SAMPLE_SIZE)
            self.assertEqual((evidence['regressions'],evidence['new_since_full_replay'],evidence['sampled']),
                             (1,2,run_fuzz.SAMPLE_SIZE))

    def test_replay_shards_cover_the_corpus_exactly_once_by_contents(self):
        with tempfile.TemporaryDirectory() as directory:
            corpus=Path(directory)/'boolean'; corpus.mkdir()
            for k in range(300): (corpus/f'{k:03}').write_bytes(k.to_bytes(2,'little')*3)
            (corpus/'empty').write_bytes(b'')
            plan=run_fuzz.shard_names(corpus,4)
            self.assertEqual(plan,run_fuzz.shard_names(corpus,4))
            names=[n for part in plan for n in part]
            self.assertEqual(sorted(names),sorted(p.name for p in corpus.iterdir()))
            self.assertEqual(len(names),len(set(names)))
            self.assertTrue(all(part==sorted(part) and len(part)>40 for part in plan))
            # The same bytes under another name (a seed and libFuzzer's copy) share a shard.
            (corpus/'copy').write_bytes((corpus/'123').read_bytes())
            again=run_fuzz.shard_names(corpus,4)
            self.assertEqual([k for k,part in enumerate(again) if 'copy' in part],
                             [k for k,part in enumerate(again) if '123' in part])
            self.assertEqual(run_fuzz.shard_names(corpus,1),[sorted(p.name for p in corpus.iterdir())])

    def test_shard_replay_runs_its_inputs_once_with_the_targets_limits(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            corpus=directory/'corpus'/'boolean'; corpus.mkdir(parents=True)
            for k in range(40): (corpus/f'input{k}').write_bytes(bytes([k,1]))
            n=run_fuzz.REPLAY_SHARDS['boolean']
            plan=run_fuzz.shard_names(corpus,n)
            seen={}
            def replay(command,log,timeout,env):
                run_corpus=Path(command[5])
                seen.update(files=sorted(p.name for p in run_corpus.iterdir()),timeout=timeout,env=env,command=command)
                log.write(f'#{len(seen["files"])+3}\tINITED cov: 9\nstat::number_of_executed_units: {len(seen["files"])+3}\n'
                          'stat::slowest_unit_time_sec: 42\nstat::peak_rss_mb: 543\n')
                return 0
            with patch('run_fuzz.FUZZ',directory), patch('run_fuzz.run_process',side_effect=replay):
                report=run_fuzz.replay_shard('boolean','n',corpus,2,directory,{})
            self.assertEqual(seen['files'],plan[2])
            self.assertEqual(report['names'],plan[2])
            self.assertFalse((corpus.parent/'.boolean-shard-2').exists())
            tail=seen['command'][seen['command'].index('--')+1:]
            for flag in ['-runs=0','-timeout=60','-rss_limit_mb=2048','-max_len=256','-print_final_stats=1']:
                self.assertIn(flag,tail)
            self.assertIn('asan-allocator',seen['command'])
            self.assertIn('malloc_context_size=5',seen['env']['ASAN_OPTIONS'])
            self.assertEqual(seen['timeout'],run_fuzz.startup_budget(len(plan[2]),60))
            self.assertEqual(report['budget_seconds'],seen['timeout'])
            self.assertEqual((report['shard'],report['shards'],report['corpus_files']),(2,n,40))
            self.assertTrue(report['passed'])
            # A replay killed before INITED, or with an input over its limit, fails.
            for change in [{'initial_executions':None},{'initial_executions':len(plan[2])},
                           {'slowest_input_seconds':61},{'peak_rss_mb':None},{'exit_code':124}]:
                self.assertFalse(run_fuzz.replay_shard_passed({**report,**change}),change)

    def test_shard_check_accepts_only_the_whole_snapshot_once(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            corpus=directory/'corpus'/'boolean'; corpus.mkdir(parents=True)
            for k in range(60): (corpus/f'input{k}').write_bytes(bytes([k,7]))
            n=run_fuzz.REPLAY_SHARDS['boolean']
            plan=run_fuzz.shard_names(corpus,n)
            digest=run_fuzz.listing_digest(p.name for p in corpus.iterdir())
            reports=[{'target':'boolean','shard':k,'shards':n,'corpus_digest':digest,'names':part,'passed':True,
                      'elapsed_seconds':100.+k,'slowest_input_seconds':k,'peak_rss_mb':500+k} for k,part in enumerate(plan)]
            with patch('run_fuzz.FUZZ',directory):
                def check(changed):
                    run_fuzz.manifest_path('boolean').unlink(missing_ok=True)
                    result=run_fuzz.check_replay_shards('boolean',corpus,changed)
                    self.assertEqual(result['manifest_written'],run_fuzz.manifest_path('boolean').exists())
                    return result
                result=check(reports)
                self.assertTrue(result['passed'])
                self.assertEqual((result['replayed_files'],result['slowest_input_seconds'],result['peak_rss_mb']),(60,n-1,500+n-1))
                self.assertEqual(run_fuzz.manifest_path('boolean').read_text().split(),sorted(p.name for p in corpus.iterdir()))
                missing=check(reports[:-1])
                self.assertFalse(missing['passed'])
                self.assertEqual(missing['unreplayed'],plan[-1])
                self.assertIn(f'shard {n-1} did not report',missing['problems'])
                self.assertFalse(check(reports+[reports[0]])['passed'])
                failed=check([{**r,'passed':k!=1} for k,r in enumerate(reports)])
                self.assertTrue(failed['complete'])
                self.assertEqual(failed['failed_shards'],[1])
                self.assertFalse(failed['passed'])
                short=[dict(r) for r in reports]; short[0]['names']=short[0]['names'][1:]
                self.assertEqual(check(short)['unreplayed'],[plan[0][0]])
                moved=[dict(r) for r in reports]
                moved[1]['names']=sorted(moved[1]['names']+[moved[0]['names'][0]])
                self.assertEqual(check(moved)['duplicated'],[plan[0][0]])
                self.assertFalse(check([{**r,'corpus_digest':'other'} for r in reports])['passed'])
                self.assertFalse(check([{**r,'shards':n-1} for r in reports])['passed'])
                # An input added after the snapshot was given to the shards.
                (corpus/'late').write_bytes(b'late')
                late=check(reports)
                self.assertEqual(late['unreplayed'],['late'])
                self.assertFalse(late['passed'])

    def test_scheduled_campaign_of_a_sharded_target_samples_whatever_the_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            corpus=directory/'corpus'/'boolean'; corpus.mkdir(parents=True)
            regressions=directory/'regressions'/'boolean'; regressions.mkdir(parents=True)
            (regressions/'r.bin').write_bytes(b'regression')
            import hashlib
            regression=hashlib.sha256(b'regression').hexdigest()
            (corpus/regression).write_bytes(b'regression')
            for k in range(300): (corpus/f'{k:03}').write_bytes(str(k).encode())
            with patch('run_fuzz.FUZZ',directory):
                chosen,evidence=run_fuzz.sharded_campaign_plan('boolean',corpus,11)
                self.assertEqual(chosen,run_fuzz.sharded_campaign_plan('boolean',corpus,11)[0])
            self.assertIn(regression,chosen)
            self.assertEqual(len(chosen),1+run_fuzz.SAMPLE_SIZE)
            self.assertEqual((evidence['replay'],evidence['shards'],evidence['regressions'],evidence['sampled']),
                             ('sharded',run_fuzz.REPLAY_SHARDS['boolean'],1,run_fuzz.SAMPLE_SIZE))

    def test_workflow_replay_jobs_match_replay_shards(self):
        import re
        text=(run_fuzz.ROOT/'.github/workflows/rust-fuzz.yml').read_text()
        jobs=text[text.index('\n  replay-snapshot:'):]
        targets={tuple(t.strip() for t in m.split(',')) for m in re.findall(r'\n\s+target: \[([^\]]*)\]',jobs)}
        self.assertEqual(targets,{tuple(sorted(run_fuzz.REPLAY_SHARDS))})
        shards=re.findall(r'\n\s+shard: \[([^\]]*)\]',jobs)
        self.assertEqual(len(shards),1)
        listed=[int(k) for k in shards[0].split(',')]
        self.assertEqual(listed,list(range(len(listed))))
        # The replay job's matrix less its exclusions: each target's shards.
        excluded=re.findall(r'\n\s+- \{target: (\w+), shard: (\d+)\}',jobs)
        self.assertEqual(len(excluded),len(set(excluded)))
        jobs_of={t:[k for k in listed if (t,str(k)) not in excluded] for t in run_fuzz.REPLAY_SHARDS}
        self.assertEqual(jobs_of,{t:list(range(n)) for t,n in run_fuzz.REPLAY_SHARDS.items()})
        self.assertLessEqual({t for t,_ in excluded},set(run_fuzz.REPLAY_SHARDS))
        self.assertEqual(len(listed),max(run_fuzz.REPLAY_SHARDS.values()))
        self.assertLessEqual(set(run_fuzz.REPLAY_SHARDS),set(run_fuzz.TARGETS))

    def test_workflow_minimize_jobs_match_minimize_shards(self):
        import re
        text=(run_fuzz.ROOT/'.github/workflows/rust-fuzz.yml').read_text()
        jobs=text[text.index('\n  minimize-snapshot:'):text.index('\n  replay-snapshot:')]
        targets=[tuple(t.strip() for t in m.split(',')) for m in re.findall(r'\n\s+target: \[([^\]]*)\]',jobs)]
        self.assertEqual(len(targets),3)
        self.assertEqual(set(targets),{tuple(sorted(run_fuzz.MINIMIZE_SHARDS))})
        shards=re.findall(r'\n\s+shard: \[([^\]]*)\]',jobs)
        self.assertEqual(len(shards),1)
        listed=[int(k) for k in shards[0].split(',')]
        self.assertEqual(listed,list(range(max(run_fuzz.MINIMIZE_SHARDS.values()))))
        excluded=re.findall(r'\n\s+- \{target: (\w+), shard: (\d+)\}',jobs)
        self.assertEqual(len(excluded),len(set(excluded)))
        jobs_of={t:[k for k in listed if (t,str(k)) not in excluded] for t in run_fuzz.MINIMIZE_SHARDS}
        self.assertEqual(jobs_of,{t:list(range(n)) for t,n in run_fuzz.MINIMIZE_SHARDS.items()})
        self.assertLessEqual({t for t,_ in excluded},set(run_fuzz.MINIMIZE_SHARDS))
        self.assertLessEqual(set(run_fuzz.MINIMIZE_SHARDS),set(run_fuzz.TARGETS))
        # The campaign job leaves these targets to the minimize jobs on a
        # minimisation run, and keeps its cache key free for their combine.
        fuzz=text[text.index('\n  fuzz:'):text.index('\n  minimize-snapshot:')]
        left=re.search(r"SHARD_MINIMIZED: \$\{\{ contains\(fromJSON\('\[([^\]]*)\]'\)",fuzz)
        self.assertEqual(sorted(t.strip().strip('"') for t in left.group(1).split(',')),sorted(run_fuzz.MINIMIZE_SHARDS))

    @staticmethod
    def fake_merge(features_of, crash=()):
        """A stand-in for the set-cover merge: it lists the run corpus's
        nonempty files smallest first, as libFuzzer does, and records each
        one's complete features (none for an input in `crash`)."""
        def merge(command,log,timeout,env):
            control=Path(next(a for a in command if a.startswith('-merge_control_file=')).split('=',1)[1])
            run_corpus=Path(command[-1])
            files=sorted((p for p in run_corpus.iterdir() if p.stat().st_size),key=lambda p:(p.stat().st_size,str(p)))
            lines=[str(len(files)),'0',*map(str,files)]
            for i,path in enumerate(files):
                lines.append(f'STARTED {i} {path.stat().st_size}')
                if path.name in crash: continue
                lines.append('FT '+' '.join([str(i)]+[str(f) for f in sorted(features_of(path.read_bytes()))]))
                lines.append(f'COV {i} {i}')
            control.write_text('\n'.join(lines)+'\n')
            return 0
        return merge

    def test_control_file_records_complete_crashed_and_unrun_inputs(self):
        text=('4\n0\n/a/x\n/a/y\n/b/z\n/a/w\nSTARTED 0 3\nFT 0 9 4 7\nCOV 0 1 2\n'
              'STARTED 1 5\nSTARTED 2 6\nFT 2\nCOV 2\n')
        inputs=run_fuzz.parse_control_file(text)
        self.assertEqual({**inputs['x'],'features':list(inputs['x']['features'])},
                         {'size':3,'features':[9,4,7],'coverage':[1,2]})
        self.assertIsNone(inputs['y']['features'])
        self.assertEqual(list(inputs['z']['features']),[])
        self.assertEqual(inputs['w'],{})
        self.assertEqual(run_fuzz.parse_control_file('0\n0\n'),{})
        with self.assertRaises(ValueError):
            run_fuzz.parse_control_file('1\n0\n/a/x\nDONE 0\n')

    def test_set_cover_keeps_every_feature_whatever_the_listing_order(self):
        import random
        rng=random.Random(5)
        inputs={f'{k:03}':{'size':rng.randrange(1,9),'features':rng.sample(range(60),rng.randrange(0,12)),'coverage':[]}
                for k in range(120)}
        chosen=run_fuzz.set_cover(inputs)
        every={f for e in inputs.values() for f in e['features']}
        self.assertEqual({f for n in chosen for f in inputs[n]['features']},every)
        self.assertLess(len(chosen),len([n for n,e in inputs.items() if e['features']]))
        shuffled=list(inputs.items()); rng.shuffle(shuffled)
        self.assertEqual(run_fuzz.set_cover(dict(shuffled)),chosen)
        # The lazy evaluation chooses what recomputing every gain would.
        covered,eager=set(),[]
        while True:
            gain,size,name=min((-len(set(e['features'])-covered),e['size'],n) for n,e in inputs.items())
            if not gain: break
            covered.update(inputs[name]['features']); eager.append(name)
        self.assertEqual(sorted(eager),chosen)
        # An input with no feature of its own is dropped; of two with the same, the smaller stays.
        self.assertEqual(run_fuzz.set_cover({'a':{'size':4,'features':[1,2]},'b':{'size':3,'features':[2,1]},
                                             'c':{'size':1,'features':[]}}),['b'])

    def test_minimization_shard_runs_its_inputs_through_the_set_cover_merge(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            corpus=directory/'corpus'/'boolean'; corpus.mkdir(parents=True)
            for k in range(40): (corpus/f'input{k}').write_bytes(bytes([k,1,k%3]))
            (corpus/'empty').write_bytes(b'')
            binary=directory/'boolean'; binary.write_bytes(b'binary')
            plan=run_fuzz.shard_names(corpus,8)
            seen={}
            merge=self.fake_merge(set)
            def run(command,log,timeout,env):
                seen.update(command=command,timeout=timeout,env=env,files=sorted(p.name for p in Path(command[-1]).iterdir()))
                return merge(command,log,timeout,env)
            with patch('run_fuzz.FUZZ',directory), patch('run_fuzz.run_process',side_effect=run):
                report=run_fuzz.minimize_shard('boolean',corpus,3,directory,{},binary)
            self.assertEqual(seen['files'],plan[3])
            self.assertEqual(seen['command'][0],str(binary))
            for flag in ['-set_cover_merge=1','-timeout=60','-rss_limit_mb=2048','-max_len=256']:
                self.assertIn(flag,seen['command'])
            self.assertIn('malloc_context_size=5',seen['env']['ASAN_OPTIONS'])
            self.assertTrue(seen['env']['ASAN_OPTIONS'].endswith('detect_odr_violation=0'))
            self.assertEqual(seen['timeout'],run_fuzz.startup_budget(len(plan[3]),60))
            self.assertEqual((report['shard'],report['shards'],report['corpus_files']),(3,8,41))
            self.assertEqual(report['listed_names'],sorted(n for n in plan[3] if n != 'empty'))
            self.assertTrue(report['passed'])
            self.assertEqual(sorted(p.name for p in corpus.parent.iterdir()),['boolean'])
            self.assertEqual(len(list(corpus.iterdir())),41)
            # An input that crashed, timed out or ran out of memory, or a merge
            # stopped at its deadline, fails the shard.
            victim=plan[3][0]
            with patch('run_fuzz.FUZZ',directory), patch('run_fuzz.run_process',side_effect=self.fake_merge(set,{victim})):
                crashed=run_fuzz.minimize_shard('boolean',corpus,3,directory,{},binary)
            self.assertEqual(crashed['failed_inputs'],[victim])
            self.assertFalse(crashed['passed'])
            for change in [{'exit_code':124},{'listed_names':report['listed_names'][1:]},{'unrun_inputs':['x']}]:
                self.assertFalse(run_fuzz.minimize_shard_passed({**report,**change}),change)

    def test_combined_minimization_keeps_every_feature_or_replaces_nothing(self):
        def features(data): return {data[0]%23,(data[0]*7+data[1])%31+100}
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            corpus=directory/'corpus'/'boolean'; corpus.mkdir(parents=True)
            for k in range(90): (corpus/f'input{k:02}').write_bytes(bytes([k,k%5,1+k%4])[:1+k%3])
            (corpus/'pad').write_bytes(b'\x05\x05')
            binary=directory/'boolean'; binary.write_bytes(b'binary')
            shards=run_fuzz.MINIMIZE_SHARDS['boolean']
            reports=[]
            with patch('run_fuzz.FUZZ',directory), patch('run_fuzz.run_process',side_effect=self.fake_merge(lambda d:features(d+b'\0'))):
                for k in range(shards):
                    out=directory/f'shard-{k}'; out.mkdir()
                    report=run_fuzz.minimize_shard('boolean',corpus,k,out,{},binary)
                    reports.append(dict(report,control_path=str(out/report['control_file'])))
            before=sorted(p.name for p in corpus.iterdir())
            every={f for n in before for f in features((corpus/n).read_bytes()+b'\0')}
            with patch('run_fuzz.FUZZ',directory):
                def refused(changed):
                    listing=sorted(p.name for p in corpus.iterdir())
                    result=run_fuzz.combine_minimized('boolean',corpus,changed)
                    self.assertFalse(result['passed'])
                    self.assertFalse(result['replaced'])
                    self.assertEqual(sorted(p.name for p in corpus.iterdir()),listing)
                    self.assertFalse(run_fuzz.manifest_path('boolean').exists())
                    return result
                self.assertIn('shard 2 did not report',refused(reports[:2]+reports[3:])['problems'])
                self.assertEqual(refused([{**r,'passed':k!=4} for k,r in enumerate(reports)])['failed_shards'],[4])
                refused(reports+[reports[0]])
                refused([{**r,'corpus_digest':'other'} for r in reports])
                refused([{**r,'binary_sha256':'other' if k==1 else r['binary_sha256']} for k,r in enumerate(reports)])
                refused([{**r,'control_path':str(directory/'missing')} if k==0 else r for k,r in enumerate(reports)])
                # A control file that lost an input's features.
                lost=directory/'lost.control'
                lines=Path(reports[5]['control_path']).read_text().splitlines()
                lost.write_text('\n'.join(l for l in lines if not l.startswith('FT 0 '))+'\n')
                refused([{**r,'control_path':str(lost)} if k==5 else r for k,r in enumerate(reports)])
                # An input added after the snapshot was given to the shards.
                (corpus/'late').write_bytes(b'\x63\x01')
                refused(reports)
                (corpus/'late').unlink()
                result=run_fuzz.combine_minimized('boolean',corpus,reports)
                self.assertTrue(result['passed'],result['problems'])
                after=sorted(p.name for p in corpus.iterdir())
                self.assertEqual(result['corpus_files_after'],len(after))
                self.assertLess(len(after),len(before))
                self.assertEqual({f for n in after for f in features((corpus/n).read_bytes()+b'\0')},every)
                self.assertEqual((result['features'],result['features_kept']),(len(every),len(every)))
                self.assertEqual(run_fuzz.manifest_path('boolean').read_text().split(),after)
                self.assertEqual(sorted(p.name for p in corpus.parent.iterdir()),['boolean','boolean.replayed'])
            # The choice depends on the features only, not on how the shards fell.
            inputs={}
            for r in reports: inputs.update(run_fuzz.parse_control_file(Path(r['control_path']).read_text()))
            self.assertEqual(after,run_fuzz.set_cover(inputs))

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

    def test_linux_boundary_keeps_startup_cap_and_final_batch_grace(self):
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
            self.assertAlmostEqual(timer.deadline(),758.19)
            self.assertEqual(evidence['shutdown_grace_seconds'],105)

    def test_tensor_budget_allows_complete_final_batch_without_extending_startup(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); log_path,stop=directory/'log',directory/'stop'
            log_path.write_text('#99 INITED cov: 12\n')
            input_seconds=run_fuzz.TARGET_INPUT_SECONDS['surface_knots']
            with patch('run_fuzz.time.monotonic',return_value=0.):
                timer=run_fuzz.MutationBudget(log_path,stop,60,shutdown_seconds=run_fuzz.shutdown_budget(input_seconds))
            with patch('run_fuzz.time.monotonic',return_value=590.): timer.tick()
            with patch('run_fuzz.time.monotonic',return_value=650.): timer.tick()
            self.assertEqual(timer.deadline(),955.)
            self.assertEqual(timer.evidence()['startup_limit_seconds'],600)
            self.assertEqual(timer.evidence()['shutdown_grace_seconds'],305)
            self.assertTrue(stop.read_bytes())
            with patch('run_fuzz.time.monotonic',return_value=0.):
                late=run_fuzz.MutationBudget(log_path,directory/'late-stop',60,shutdown_seconds=run_fuzz.shutdown_budget(input_seconds))
            with patch('run_fuzz.time.monotonic',return_value=600.01): late.tick()
            self.assertEqual(late.deadline(),600.)
            self.assertFalse(late.evidence()['startup_budget_completed'])

    def test_stop_file_is_observed_after_complete_mutation_batch(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); log_path,stop=directory/'log',directory/'stop'
            # Mirrors the pinned runtime's outer stop-file check. Five legal
            # callbacks can outlive a grace period sized for only one input.
            script=('import time,pathlib,sys; print("#99 INITED cov: 12",flush=True); '
                    '\nwhile not pathlib.Path(sys.argv[1]).exists():\n'
                    ' for _ in range(5): time.sleep(.1)\n'
                    'print("stat::number_of_executed_units: 104",flush=True)')
            timer=run_fuzz.MutationBudget(log_path,stop,.1,startup_seconds=1.,shutdown_seconds=.6)
            with log_path.open('w') as log:
                code=run_fuzz.run_process([sys.executable,'-c',script,str(stop)],log,2.,os.environ,timer.tick,timer.deadline)
            self.assertEqual(code,0)
            self.assertTrue(timer.evidence()['mutation_budget_completed'])
            self.assertEqual(run_fuzz.statistics(log_path.read_text())['mutation_executions'],5)

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
