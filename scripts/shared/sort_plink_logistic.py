#!/usr/bin/env python3
"""Sort PLINK2 logistic results by named P column, retaining header and all rows."""
import csv, math, sys
from pathlib import Path
if len(sys.argv) != 3: sys.exit("Usage: sort_plink_logistic.py INPUT OUTPUT")
with open(sys.argv[1], newline="") as f:
    reader=csv.DictReader(f, delimiter="\t")
    if not reader.fieldnames or "P" not in reader.fieldnames: sys.exit("Missing P column")
    fields=reader.fieldnames
    rows=list(reader)
def key(row):
    try:
        p=float(row["P"])
        return p if math.isfinite(p) else float("inf")
    except (TypeError, ValueError): return float("inf")
rows.sort(key=key)
Path(sys.argv[2]).parent.mkdir(parents=True, exist_ok=True)
with open(sys.argv[2],"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    w.writeheader();w.writerows(rows)
