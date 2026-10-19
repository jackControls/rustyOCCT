"""Record what a fixture generator reads, for generator_groups.py.

Loaded through PYTHONPATH when GENERATOR_AUDIT_DIR is set, in the generator
and in every Python process it starts (multiprocessing workers included).
Each process appends, to its own file in that directory, one line per new
event: 'open' for a file opened for reading, 'list' for a directory listed
(os.listdir, os.scandir, so glob, iterdir and walk too), 'stat' for a path
whose existence or kind was asked (os.stat, os.lstat, so Path.exists and
os.path.exists), and 'exec' for a program started. Paths are absolute;
generator_groups.py keeps those inside the checkout.
"""
import os
import sys

_DIR = os.environ.get('GENERATOR_AUDIT_DIR')


def _install():
    state = {'pid': None, 'fd': None, 'seen': set(), 'busy': False}
    write, getpid, fspath, abspath = os.write, os.getpid, os.fspath, os.path.abspath
    os_open, flags = os.open, os.O_WRONLY | os.O_CREAT | os.O_APPEND
    audit_dir = abspath(_DIR)

    def emit(kind, path):
        if state['busy']:
            return
        state['busy'] = True
        try:
            if isinstance(path, int) or path is None and kind != 'list':
                return
            path = abspath(fspath(path if path is not None else '.'))
            if isinstance(path, bytes):
                path = os.fsdecode(path)
            if path.startswith(audit_dir):
                return
            key = (kind, path)
            if state['pid'] != getpid():
                # A forked child: its own file, the parent's events already written.
                state['pid'], state['fd'] = getpid(), os_open(os.path.join(audit_dir, f'{getpid()}.log'), flags)
            if key in state['seen']:
                return
            state['seen'].add(key)
            write(state['fd'], f'{kind}\t{path}\n'.encode('utf-8', 'surrogateescape'))
        except Exception:
            pass
        finally:
            state['busy'] = False

    def hook(event, args):
        if event == 'open':
            path, mode, oflags = args
            if mode is None:
                reading = oflags is None or not oflags & (os.O_WRONLY | os.O_RDWR) or oflags & os.O_RDWR
            else:
                reading = 'r' in mode or '+' in mode
            if reading:
                emit('open', path)
        elif event in ('os.listdir', 'os.scandir'):
            emit('list', args[0])
        elif event in ('subprocess.Popen', 'os.posix_spawn', 'os.exec', 'os.spawn'):
            # Popen: (executable, args, cwd, env); posix_spawn and exec:
            # (path, argv, env); spawn: (mode, path, args, env).
            program, argv = (args[1], args[2]) if event == 'os.spawn' else (args[0], args[1])
            if program is None:
                program = argv if isinstance(argv, (str, bytes)) else (list(argv)[0] if argv else '')
            emit_exec(program)
        elif event == 'os.system':
            emit_exec('/bin/sh')

    def emit_exec(program):
        # Not a path to resolve: a bare name runs whatever PATH finds.
        if state['busy']:
            return
        name = os.fsdecode(os.fspath(program)) if isinstance(program, (str, bytes, os.PathLike)) else repr(program)
        state['busy'] = True
        try:
            if state['pid'] != getpid():
                state['pid'], state['fd'] = getpid(), os_open(os.path.join(audit_dir, f'{getpid()}.log'), flags)
            if ('exec', name) not in state['seen']:
                state['seen'].add(('exec', name))
                write(state['fd'], f'exec\t{name}\n'.encode('utf-8', 'surrogateescape'))
        except Exception:
            pass
        finally:
            state['busy'] = False

    class StatWrapper:
        # A callable object, not a function: Python 3.9's pathlib stores
        # os.stat as a class attribute, where a function would bind self.
        def __init__(self, original):
            self.__wrapped__ = original

        def __call__(self, path, *args, **kwargs):
            if not isinstance(path, int):
                emit('stat', path)
            return self.__wrapped__(path, *args, **kwargs)

    os.stat, os.lstat = StatWrapper(os.stat), StatWrapper(os.lstat)
    sys.addaudithook(hook)


if _DIR:
    _install()
