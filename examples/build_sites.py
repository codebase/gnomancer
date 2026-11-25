"""Rebuild gnomancer/data/sites.tsv: where every site gnomancer genotypes by name sits on GRCh37 and GRCh38.

Array files name their sites (rsIDs), so the array path never needed coordinates. A sequencing VCF is the other
way round: positions and alleles, often no rsIDs, and nothing at all where the sample matches the reference. To
read the marker panel and the ancestry SNPs from a whole genome, gnomancer needs, for each rsID, its position on
both builds, the reference allele there and the alternate alleles, in VCF form.

Sources:
- SNVs: Ensembl REST variation endpoint (GRCh38 from rest.ensembl.org, GRCh37 from grch37.rest.ensembl.org).
- Indels: the ClinVar VCFs (already left-aligned, VCF-style), looked up by rsID; dbSNP places indels at the
  right-hand end of a repeat, so Ensembl's coordinates are not what a variant caller writes.
- INDEL_OVERRIDES: indel sites that are not in ClinVar, entered by hand from the reference sequence.

Usage: python examples/build_sites.py   (needs network, and `gnomancer setup` for the two ClinVar files)
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from gnomancer.paths import DATA, tool  # noqa: E402

SERVERS = {37: "https://grch37.rest.ensembl.org", 38: "https://rest.ensembl.org"}
CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
COMP = str.maketrans("ACGT", "TGCA")

# rsid -> {build: (chrom, pos, ref, alt)}; only for indel markers ClinVar does not carry.
# ABO: the reference genome carries the O allele (the 261delG deletion), so a caller writes the A/B allele as an
# insertion. Checked against the reference sequence on both builds.
INDEL_OVERRIDES: dict = {"rs8176719": {37: ("9", 136132908, "T", "TC"), 38: ("9", 133257521, "T", "TC")}}


def post(server: str, ids: list[str]) -> dict:
    req = urllib.request.Request(server + "/variation/homo_sapiens", data=json.dumps({"ids": ids}).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            print(f"  retry after error: {e}", file=sys.stderr)
            time.sleep(5)
    return {}


def ensembl(ids: list[str], build: int) -> dict:
    """rsid -> (chrom, pos, ref, [alts]) for SNVs on the primary assembly, forward strand."""
    out = {}
    for i in range(0, len(ids), 150):
        res = post(SERVERS[build], ids[i:i + 150])
        for rsid, v in res.items():
            names = {rsid, *v.get("synonyms", [])}
            for m in v.get("mappings", []):
                if m.get("seq_region_name") not in CHROMS or m.get("start") != m.get("end"):
                    continue
                alleles = m["allele_string"].split("/")
                if m.get("strand") == -1:
                    alleles = [a.translate(COMP)[::-1] for a in alleles]
                if any(len(a) != 1 or a not in "ACGT" for a in alleles):
                    continue
                for n in names & set(ids):
                    out[n] = (m["seq_region_name"], int(m["start"]), alleles[0], alleles[1:])
                break
    return out


def clinvar_indels(kinds: dict, build: int) -> dict:
    """rsid -> (chrom, pos, ref, [alt]) from the ClinVar VCF of that build. `kinds` says, per rsID, whether the
    panel's effect allele is the deletion ('D') or the insertion ('I'); an rsID can cover both at one site
    (BRCA1 c.68_69del and c.68_69dup share rs80357914), and only the matching record is the marker."""
    ids = set(kinds)
    path = tool("clinvar", f"clinvar_GRCh{build}.vcf.gz")
    want = {i[2:] for i in ids}
    out: dict = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t", 8)
            j = p[7].find("RS=")
            if j < 0:
                continue
            for rs in p[7][j + 3:].split(";", 1)[0].split("|"):
                kind = "D" if len(p[4]) < len(p[3]) else "I" if len(p[4]) > len(p[3]) else ""
                if rs in want and kind == kinds["rs" + rs] and "rs" + rs not in out:
                    out["rs" + rs] = (p[0], int(p[1]), p[3], [p[4]])
    return out


def main() -> None:
    markers = [r for r in csv.DictReader((l for l in open(DATA / "markers.tsv") if not l.startswith("#")), delimiter="\t")]
    indel_ids = {m["rsid"]: m["effect"] for m in markers if m["ref"] in "ID" or m["effect"] in "ID"}
    snv_ids = [m["rsid"] for m in markers if m["rsid"] not in indel_ids]
    sets = {m["rsid"]: "panel" for m in markers}
    ais = DATA / "aisnps.txt"
    if ais.exists():
        for rs in ais.read_text().split():
            if rs not in sets:
                sets[rs] = "ancestry"
                snv_ids.append(rs)
    found = {b: {**ensembl(snv_ids, b), **clinvar_indels(indel_ids, b)} for b in (37, 38)}
    for rs, per_build in INDEL_OVERRIDES.items():
        for b, (chrom, pos, ref, alt) in per_build.items():
            found[b][rs] = (chrom, pos, ref, [alt])
    rows, missing = [], []
    for rs, which in sets.items():
        a, b = found[37].get(rs), found[38].get(rs)
        if not a and not b:
            missing.append(rs)
            continue
        chrom = (b or a)[0]
        # The reference allele can differ between builds; keep each build's own.
        rows.append([rs, which, chrom, a[1] if a else "", a[2] if a else "", ",".join(a[3]) if a else "",
                     b[1] if b else "", b[2] if b else "", ",".join(b[3]) if b else ""])
    with open(DATA / "sites.tsv", "w") as f:
        f.write("# Positions and VCF-style alleles of every site gnomancer reads by rsID, on both builds.\n"
                "# Generated by examples/build_sites.py from Ensembl (SNVs) and ClinVar (indels); do not edit by hand.\n")
        f.write("rsid\tset\tchrom\tpos37\tref37\talt37\tpos38\tref38\talt38\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} sites written; missing: {missing or 'none'}")


if __name__ == "__main__":
    main()
