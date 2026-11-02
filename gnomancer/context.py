"""Short, population-level context for Y and mtDNA haplogroups (data/haplogroup_context.tsv)."""
from __future__ import annotations

import csv
from pathlib import Path

from .paths import DATA as _DATA

DATA = _DATA / "haplogroup_context.tsv"


def _rows():
    with open(DATA, newline="") as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def lookup(kind: str, *labels: str) -> list[dict]:
    """Return context rows for a haplogroup, most specific first.

    kind: 'Y' or 'MT'. labels: any names for the call (e.g. short 'R-M269' and long 'R1b1a2a1a2').
    Matching is by prefix on each label; results are the chain of matching prefixes, longest first,
    so a report can show the lineage from broad to specific.
    """
    hits = {}
    for r in _rows():
        if r["type"] != kind:
            continue
        for lab in labels:
            if lab and (lab == r["prefix"] or lab.startswith(r["prefix"])):
                hits.setdefault(r["name"], r)  # one entry per lineage name, whichever prefix matched first
    return sorted(hits.values(), key=lambda r: -int(r.get("depth") or len(r["prefix"])))
