"""Synthetic family (two parents, two children) in AncestryDNA layout, with real inheritance structure:
each parent has two haplotypes per chromosome; each child inherits one recombined haplotype from each parent.
Used by tests and demos for `gnomancer family`. No real person is represented."""
from __future__ import annotations

import random
import sys
from pathlib import Path

CHR_LEN = {str(i): L for i, L in zip(range(1, 23), [249e6, 243e6, 198e6, 191e6, 181e6, 171e6, 159e6, 146e6, 141e6, 135e6, 135e6, 134e6, 115e6, 107e6, 102e6, 90e6, 81e6, 78e6, 59e6, 63e6, 48e6, 51e6])}
ANCHORS = [("rs3094315", "1", 752566), ("rs3131972", "1", 752721), ("rs11928389", "3", 50927009)]


def main(outdir: str, seed: int = 11, per_chrom: int = 4000):
    rng = random.Random(seed)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    sites = []  # (rsid, chrom, pos, ref, alt)
    for ch in CHR_LEN:
        n = max(800, int(per_chrom * CHR_LEN[ch] / 249e6))
        for pos in sorted(rng.sample(range(100_000, int(CHR_LEN[ch]) - 100_000), n)):
            ref, alt = rng.sample("ACGT", 2)
            sites.append((f"rs{len(sites) + 10_000_000}", ch, pos, ref, alt))
    for rsid, ch, pos in ANCHORS:
        sites.append((rsid, ch, pos, "A", "G"))
    # the curated panel sites, so family reports show side-by-side genotypes and inheritance
    import csv
    with open(Path(__file__).resolve().parent.parent / "gnomancer" / "data" / "markers.tsv", newline="") as f:
        for m in csv.DictReader([l for l in f if l.strip() and not l.startswith("#")], delimiter="\t"):
            ch = rng.choice(list(CHR_LEN))
            sites.append((m["rsid"], ch, rng.randint(1_000_000, int(CHR_LEN[ch]) - 1_000_000), m["ref"], m["effect"]))
    sites.sort(key=lambda s: (int(s[1]), s[2]))
    by_chr = {}
    for i, s in enumerate(sites):
        by_chr.setdefault(s[1], []).append(i)

    def haplotype():
        return [s[4] if rng.random() < 0.4 else s[3] for s in sites]

    def recombine(h1, h2):
        child = [None] * len(sites)
        for ch, idx in by_chr.items():
            cut = rng.choice(idx)
            src = rng.random() < 0.5
            for i in idx:
                if i == cut:
                    src = not src
                child[i] = h1[i] if src else h2[i]
        return child

    people = {}
    dad = (haplotype(), haplotype())
    mom = (haplotype(), haplotype())
    people["dad"] = dad
    people["mom"] = mom
    people["child1"] = (recombine(*dad), recombine(*mom))
    people["child2"] = (recombine(*dad), recombine(*mom))
    people["stranger"] = (haplotype(), haplotype())
    sexes = dict(dad="M", mom="F", child1="M", child2="F", stranger="F")
    for name, (h1, h2) in people.items():
        p = out / f"{name}_AncestryDNA.txt"
        with open(p, "w", newline="") as f:
            f.write("#AncestryDNA raw data download\r\n#synthetic family demo\r\n#Data was collected using AncestryDNA array version: V2.0\r\n")
            f.write("rsid\tchromosome\tposition\tallele1\tallele2\r\n")
            for s, a, b in zip(sites, h1, h2):
                f.write(f"{s[0]}\t{s[1]}\t{s[2]}\t{a}\t{b}\r\n")
            # a few X sites so sex inference works
            for k in range(200):
                pos = 3_000_000 + k * 700_000
                a = rng.choice("ACGT")
                b = a if sexes[name] == "M" else rng.choice("ACGT")
                f.write(f"rs{90_000_000 + k}\t23\t{pos}\t{a}\t{b}\r\n")
        print(p)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "family")
