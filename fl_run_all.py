#!/usr/bin/env python3
"""Run the whole pipeline: preprocess, then the statewide analysis.

    python fl_run_all.py

Exactly equivalent to running the two scripts in order. It exists so that a
scheduled or unattended run is one command, and so that a failure in
preprocessing stops before the statewide script starts on a cache that is not
there.

Each script runs in its own namespace, so nothing leaks between them - which
matters, because both define names like `counties`.

Everything either script does is controlled by `fl_config.py`. To skip
preprocessing on a re-run, set `RUN_PREPROCESS = False` there rather than
editing this file.
"""
import os
import runpy
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

STEPS = [
    ('fl_01_preprocess.py', 'PREPROCESSING — build the cache'),
    ('fl_02_statewide.py',  'STATEWIDE — passes A, B and C'),
]

t_all = time.time()
for i, (script, title) in enumerate(STEPS, 1):
    path = os.path.join(HERE, script)
    print(f'\n\n{"#" * 74}\n#  STEP {i} of {len(STEPS)}  ·  {title}\n'
          f'#  {path}\n{"#" * 74}\n', flush=True)
    t0 = time.time()
    try:
        runpy.run_path(path, run_name='__main__')
    except SystemExit as e:
        if e.code:
            raise
    except Exception as e:
        print(f'\n{"!" * 74}\n'
              f'STEP {i} FAILED after {time.time() - t0:,.0f}s: '
              f'{type(e).__name__}: {e}\n'
              f'Nothing downstream ran. Fix the cause and re-run - both steps '
              f'resume\nfrom what is already on disk, so you do not repeat '
              f'work that succeeded.\n{"!" * 74}', flush=True)
        raise
    print(f'\n[step {i} finished in {time.time() - t0:,.0f}s]', flush=True)

print(f'\n{"#" * 74}\n#  ALL STEPS COMPLETE in {time.time() - t_all:,.0f}s\n'
      f'{"#" * 74}')
