#!/usr/bin/env python3
"""Sort additive PLINK2 logistic rows by numeric P value.
Usage: python3 sort_plink_logistic.py INPUT.tsv OUTPUT.tsv
"""
import csv
import math
import sys

if len(sys.argv) != 3:
    sys.exit("Usage: sort_plink_logistic.py INPUT.tsv OUTPUT.tsv")
with open(sys.argv[1], newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    if not reader.fieldnames or not {"ID", "P", "TEST"}.issubset(reader.fieldnames):
        sys.exit("Input requires ID, P and TEST columns")
    fields = reader.fieldnames
    rows = []
    for row in reader:
        if row["TEST"] != "ADD":
            continue
        try:
            p = float(row["P"])
        except (TypeError, ValueError):
            continue
        if math.isfinite(p):
            rows.append((p, row))
rows.sort(key=lambda item: item[0])
with open(sys.argv[2], "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(row for _, row in rows)
