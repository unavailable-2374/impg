#!/usr/bin/env python3
"""Resume a validated real journal prefix, SIGKILL after new queries, resume again."""
import argparse
import os
from pathlib import Path
import struct
import subprocess
import time
p=argparse.ArgumentParser()
p.add_argument('binary',type=Path);p.add_argument('output',type=Path);p.add_argument('window',choices=['w5k','w20k'])
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
d=Path('/scratch/10779/shuocao2374/VGP/ILS/partition')
source=d/'verify_perf'/a.window/a.window/'p3.journal'
b=source.read_bytes();assert b[:8]==b'IMPGPJ01'
pos=12+struct.unpack_from('<I',b,8)[0];offsets=[pos]
while pos<len(b):
    magic,seq,start,end,n=struct.unpack_from('<IIiiI',b,pos);assert magic==0x50574A31
    stop=pos+20+n*12;assert stop+8<=len(b)
    h=0xcbf29ce484222325
    for v in b[pos:stop]: h=((h^v)*0x100000001b3)&((1<<64)-1)
    assert h==struct.unpack_from('<Q',b,stop)[0]
    pos=stop+8;offsets.append(pos)
cut=offsets[-1001];journal=a.output.resolve()/'resume.journal';journal.write_bytes(b[:cut])
# Include a genuine torn record from the old file. Recovery must discard it.
with journal.open('ab') as out: out.write(b[cut:cut+13])
env=os.environ.copy();env.update(LD_LIBRARY_PATH='/scratch/10779/shuocao2374/tool/htslib-shim/lib:/scratch/10779/shuocao2374/tool/miniconda3/lib',IMPG_PARTITION_MEMORY_STATS='1',IMPG_PARTITION_JOURNAL=str(journal),IMPG_PARTITION_JOURNAL_TAG='verify-perf',IMPG_PARTITION_JOURNAL_SYNC_SECS='0')
w,dist,floor=(5000,500,2500) if a.window=='w5k' else (20000,2000,10000)
cmd=[str(a.binary.resolve()),'partition','--alignment-list',str(d/'verify/fresh/list'),'--index-mode','per-file','-w',str(w),'-d',str(dist),'--min-missing-size',str(floor),'--min-boundary-distance',str(floor),'--max-depth','1','-t','12','-v','1','--output-folder',str(a.output.resolve())]
with (a.output/'interrupted.log').open('w') as log:
    proc=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
    while proc.poll() is None and journal.stat().st_size<offsets[-901]: time.sleep(.02)
    assert proc.poll() is None, 'did not interrupt new queries'
    proc.kill(); assert proc.wait()==-9
assert journal.stat().st_size>cut
with (a.output/'resumed.log').open('w') as log: subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
subprocess.run(['cmp',str(d/'verify_perf'/a.window/a.window/'ref/partitions.bed'),str(a.output/'partitions.bed')],check=True)
subprocess.run(['cmp',str(source),str(journal)],check=True)
print(f'PASS old journal + torn tail + SIGKILL + resume: BED and journal byte identical: {a.output}',flush=True)
