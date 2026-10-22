#!/usr/bin/env python3
"""Compare real PAF queries across disk-backed single/multi index and BFS/DFS."""
import os
from pathlib import Path
import subprocess
root=Path(__file__).resolve().parents[1]
out=root/'validation/query-regression-checked';out.mkdir(exist_ok=False)
d=Path('/scratch/10779/shuocao2374/VGP/ILS/partition');lst=d/'verify/fresh/list';paf=lst.read_text().splitlines()[0]
env=os.environ.copy();env['LD_LIBRARY_PATH']='/scratch/10779/shuocao2374/tool/htslib-shim/lib:/scratch/10779/shuocao2374/tool/miniconda3/lib';env['IMPG_TREE_LOAD_STATS']='1'
for mode in ['single','multi']:
    inp=['-a',paf,'--index',paf+'.impg'] if mode=='single' else ['--alignment-list',str(lst),'--index-mode','per-file']
    for walk in ['bfs','dfs']:
        args=['query',*inp,'-r','GCA_965226575.1#0#OZ245636.1:10000000-10005000','--no-merge','--transitive','--max-depth','1','-t','12','-v','1']
        if walk=='dfs':args+=['--transitive-dfs']
        for name,binary in [('old','/scratch/10779/shuocao2374/tool/bin-frozen/impg-3351b8fd'),('new',str(root/'validation/bin/final'))]:
            with (out/f'{mode}-{walk}-{name}.bed').open('wb') as bed,(out/f'{mode}-{walk}-{name}.log').open('wb') as log:
                subprocess.run([binary,*args],env=env,stdout=bed,stderr=log,check=True)
        old=out/f'{mode}-{walk}-old.bed';new=out/f'{mode}-{walk}-new.bed'
        assert old.stat().st_size>0
        subprocess.run(['cmp',str(old),str(new)],check=True)
        print(f'PASS real {mode} {walk} max-depth=1 query',flush=True)
