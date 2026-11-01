#!/usr/bin/env python3
"""Identify duplicate chr:pos:allele markers; keep lowest missingness.
Usage: python3 get_duprem_var.py PREFIX (expects PREFIX.bim and PREFIX.lmiss)
"""
import argparse
from collections import defaultdict
p=argparse.ArgumentParser();p.add_argument('prefix');a=p.parse_args()
missing={}
with open(a.prefix+'.lmiss') as fh:
    header=fh.readline().split()
    for line in fh:
        row=dict(zip(header,line.split()))
        missing[row['SNP']]=float(row['F_MISS'])
groups=defaultdict(list)
with open(a.prefix+'.bim') as fh:
    for line in fh:
        f=line.split()
        groups[(f[0],f[3],tuple(sorted(f[4:6])))].append(f[1])
with open(a.prefix+'.remdup.snp','w') as out:
    for ids in groups.values():
        if len(ids)>1:
            winner=min(ids,key=lambda x:(missing.get(x,float('inf')),x))
            for vid in ids:
                if vid!=winner:out.write(vid+'\n')
