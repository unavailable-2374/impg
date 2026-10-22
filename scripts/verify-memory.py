#!/usr/bin/env python3
"""Real six-genome regression runner. Never modifies original journals or inputs."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

p = argparse.ArgumentParser()
p.add_argument('binary', type=Path)
p.add_argument('output', type=Path)
p.add_argument('window', choices=['w5k', 'w20k'])
p.add_argument('--journal', type=Path)
p.add_argument('--keep-trees', action='store_true')
p.add_argument('--threads', type=int, default=12)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
d = Path('/scratch/10779/shuocao2374/VGP/ILS/partition')
env = os.environ.copy()
env['LD_LIBRARY_PATH'] = '/scratch/10779/shuocao2374/tool/htslib-shim/lib:/scratch/10779/shuocao2374/tool/miniconda3/lib'
env['IMPG_PARTITION_MEMORY_STATS'] = '1'
if a.keep_trees:
    env['IMPG_PARTITION_KEEP_TREES'] = '1'
if a.journal:
    shutil.copyfile(a.journal, a.output / 'resume.journal')
    env['IMPG_PARTITION_JOURNAL'] = str(a.output.resolve() / 'resume.journal')
    env['IMPG_PARTITION_JOURNAL_TAG'] = 'verify-perf'
w, dist, floor = (5000, 500, 2500) if a.window == 'w5k' else (20000, 2000, 10000)
cmd = ['/usr/bin/time', '-v', str(a.binary.resolve()), 'partition', '--alignment-list', str(d / 'verify/fresh/list'), '--index-mode', 'per-file', '-w', str(w), '-d', str(dist), '--min-missing-size', str(floor), '--min-boundary-distance', str(floor), '--max-depth', '1', '-t', str(a.threads), '-v', '1', '--output-folder', str(a.output.resolve())]
with (a.output / 'run.log').open('w') as log:
    log.write(repr(cmd) + '\n'); log.flush()
    subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
reference = d / 'verify_perf' / a.window / a.window / 'ref/partitions.bed'
subprocess.run(['cmp', str(reference), str(a.output / 'partitions.bed')], check=True)
print(f'PASS byte identical: {a.output}', flush=True)
