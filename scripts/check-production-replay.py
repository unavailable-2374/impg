#!/usr/bin/env python3
"""Compare normalized partition logs against the original production run."""
import hashlib
from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
source=Path('/scratch/10779/shuocao2374/VGP/ILS/partition/logs/partition.3566336.out')
expected={};count=0;digest=hashlib.sha256()
with source.open('rb') as f:
    for raw in f:
        if b'Computed partition' not in raw: continue
        if not raw.endswith(b'\n'): break
        line=raw[raw.index(b'Computed partition'):];digest.update(line);count+=1
        if count%100000==0: expected[count]=digest.hexdigest()
checked=0
for line in (root/'logs/production-probe.log').read_text().splitlines():
    m=re.match(r'computed=(\d+).*digest=([0-9a-f]+)',line)
    if m:
        n=int(m[1])
        if n in expected:
            assert expected[n]==m[2], f'replay divergence at partition count {n}'
            checked=n
assert checked>0
print(f'PASS normalized original production log: {checked} consecutive partitions; original complete lines={count}',flush=True)
