"""Fetch 1000 Genomes phase 3 super-population frequencies for every panel site (Ensembl GRCh37 REST).

Writes gnomancer/data/marker_frequencies.tsv with the EFFECT allele frequency per super-population.
Run occasionally by a maintainer; the output is committed so users never need the network for this.
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MARKERS = ROOT / "gnomancer" / "data" / "markers.tsv"
OUT = ROOT / "gnomancer" / "data" / "marker_frequencies.tsv"
POPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def fetch(rsid: str) -> dict | None:
    url = f"https://grch37.rest.ensembl.org/variation/human/{rsid}?pops=1"
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:  # noqa: S310
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            code = getattr(e, "code", None)
            if code == 400:
                return None
            time.sleep(2 * (attempt + 1))
    return None


def effect_freqs(j: dict, ref: str, effect: str) -> tuple[dict, str]:
    """Return {pop: effect-allele frequency} and a note on how the allele was matched."""
    alleles = (j.get("mappings") or [{}])[0].get("allele_string", "")
    pops = {}
    for p in j.get("populations", []):
        name = p.get("population", "")
        if not name.startswith("1000GENOMES:phase_3:"):
            continue
        pop = name.split(":")[-1]
        if pop in POPS:
            pops.setdefault(pop, {})[p["allele"]] = float(p["frequency"])
    if not pops:
        return {}, "no 1000G data"
    present = set(a for d in pops.values() for a in d)
    target = effect
    note = "direct"
    if effect in ("I", "D"):
        if "-" in present:
            target = "-" if effect == "D" else next((a for a in present if a != "-"), None)
            note = "indel"
        else:
            return {}, "indel alleles not comparable"
    elif effect not in present:
        comp = COMP.get(effect)
        if comp in present and {ref, effect} not in ({"A", "T"}, {"C", "G"}):
            target = comp
            note = "complement"
        else:
            return {}, f"effect allele {effect} not in {alleles}"
    out = {pop: d.get(target, 0.0) for pop, d in pops.items()}
    return out, note


def main():
    with open(MARKERS, newline="") as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    markers = list(csv.DictReader(lines, delimiter="\t"))
    rows = []
    for i, m in enumerate(markers, 1):
        j = fetch(m["rsid"])
        if j is None:
            rows.append(dict(rsid=m["rsid"], note="lookup failed", **{p: "" for p in POPS}))
        else:
            freqs, note = effect_freqs(j, m["ref"], m["effect"])
            rows.append(dict(rsid=m["rsid"], note=note, **{p: (f"{freqs[p]:.4f}" if p in freqs else "") for p in POPS}))
        print(f"{i}/{len(markers)} {m['rsid']} {rows[-1]['note']} {rows[-1]['EUR']}", file=sys.stderr)
        time.sleep(0.15)
    with open(OUT, "w", newline="") as f:
        f.write("# Effect-allele frequency per 1000 Genomes phase 3 super-population (Ensembl GRCh37 REST). Regenerate with examples/fetch_frequencies.py\n")
        w = csv.DictWriter(f, fieldnames=["rsid", *POPS, "note"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} rows, {sum(1 for r in rows if r['EUR'])} with data)")


if __name__ == "__main__":
    main()
