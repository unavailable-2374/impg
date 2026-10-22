#!/usr/bin/env python3
"""Read-only source journal; bounded continuation uses an independent copy."""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import struct
import subprocess
import time
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);p.add_argument('--queries',type=int,default=100);p.add_argument('--threads',type=int,default=12);p.add_argument('--source-journal',type=Path)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
d=Path('/scratch/10779/shuocao2374/VGP/ILS/partition');source=a.source_journal or d/'run_w5k/partition.journal';journal=a.output.resolve()/'partition.journal'
with source.open('rb') as f:
    assert f.read(8)==b'IMPGPJ01';n=struct.unpack('<I',f.read(4))[0];header=f.read(n).decode()
tag=header.split(' tag=',1)[1]
print(f'copy source={source} bytes={source.stat().st_size}',flush=True)
shutil.copyfile(source,journal)
env=os.environ.copy();env.update(LD_LIBRARY_PATH='/scratch/10779/shuocao2374/tool/htslib-shim/lib:/scratch/10779/shuocao2374/tool/miniconda3/lib',IMPG_PARTITION_MEMORY_STATS='1',IMPG_PARTITION_JOURNAL=str(journal),IMPG_PARTITION_JOURNAL_TAG=tag,IMPG_PARTITION_JOURNAL_SYNC_SECS='0')
cmd=[str(a.binary.resolve()),'partition','--alignment-list',str(d/'alignments.list'),'-w','5000','-d','500','--min-missing-size','2500','--min-boundary-distance','2500','--max-depth','1','--index-mode','per-file','-t',str(a.threads),'-v','1','--output-folder',str(a.output.resolve())]
with (a.output/'probe.log').open('w',buffering=1) as log:
    log.write(repr(cmd)+'\n')
    proc=subprocess.Popen(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    count=0;queried=0;querying=False;digest=hashlib.sha256();start=time.monotonic();stopped=False
    try:
        for raw in proc.stdout:
            if b'now querying' in raw:
                querying=True;log.write(raw.decode());print(raw.decode().strip(),flush=True)
            if b'Computed partition' in raw:
                line=raw[raw.index(b'Computed partition'):];digest.update(line);count+=1
                if querying:
                    queried+=1
                    log.write(line.decode())
                if count%100000==0:
                    text=f'computed={count} queried={queried} elapsed_s={time.monotonic()-start:.1f} digest={digest.hexdigest()} last={line.decode().strip()}'
                    log.write(text+'\n');print(text,flush=True)
                if queried>=a.queries and not stopped:
                    proc.terminate();stopped=True
            else:
                log.write(raw.decode(errors='replace'))
        rc=proc.wait();log.write(f'exit={rc} computed={count} queried={queried} elapsed_s={time.monotonic()-start:.1f} digest={digest.hexdigest()}\n')
    finally:
        if proc.poll() is None: proc.kill();proc.wait()
    if not stopped or rc not in (-15,0): raise RuntimeError(f'production probe incomplete: exit={rc} queried={queried}')
print(f'PASS production replay + {queried} continuing queries; stopped own process intentionally',flush=True)
