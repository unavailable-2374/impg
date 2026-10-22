#!/usr/bin/env python3
"""Check real production replay digests and durable continuation records."""
import hashlib
from pathlib import Path
import re
import struct
import sys
new_label,old_label = sys.argv[1:3] if len(sys.argv)==3 else ("production","production-original")
root=Path(__file__).resolve().parents[1]
def digests(name):
    result={}
    for line in (root/'logs'/name).read_text().splitlines():
        m=re.match(r'computed=(\d+).*digest=([0-9a-f]+)',line)
        if m: result[int(m[1])]=m[2]
    return result
new=digests(new_label+'-probe.log');old=digests(old_label+'-probe.log')
common=set(new)&set(old);assert common
for n in common: assert new[n]==old[n], f'different normalized partition log through {n}'
print(f'PASS production replay digests through {max(common)} partitions')
def replay_base(label):
    text=(root/'logs'/(label+'-probe.log')).read_text()
    size=int(re.search(r'^copy source=.* bytes=(\d+)$',text,re.M)[1])
    probe=(root/'validation'/label/'probe.log').read_text()
    torn=re.search(r'dropping (\d+) bytes of torn tail',probe)
    return size-(int(torn[1]) if torn else 0)
base=replay_base(new_label)
assert base==replay_base(old_label),'source snapshots differed'
paths=[root/'validation'/name/'partition.journal' for name in [new_label,old_label]]
size=min(p.stat().st_size for p in paths);assert size>base,'no durable continuation yet'
with paths[0].open('rb') as a,paths[1].open('rb') as b:
    pos=0
    while pos<size:
        length=min(8<<20,size-pos);assert a.read(length)==b.read(length),f'journal differs at chunk {pos}'
        pos+=length
with paths[0].open('rb') as f:
    f.seek(base);data=f.read(size-base)
pos=0;count=0
while pos+20<=len(data):
    magic,seq,start,end,n=struct.unpack_from('<IIiiI',data,pos);assert magic==0x50574A31
    stop=pos+20+n*12
    if stop+8>len(data): break
    h=0xcbf29ce484222325
    for v in data[pos:stop]: h=((h^v)*0x100000001b3)&((1<<64)-1)
    assert h==struct.unpack_from('<Q',data,stop)[0]
    count+=1;pos=stop+8
assert count>0
print(f'PASS common journal prefix {size} bytes; {count} durable newly queried records byte identical with valid checksums')
