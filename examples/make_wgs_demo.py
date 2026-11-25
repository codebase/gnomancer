"""Write a synthetic whole-genome VCF (GRCh38) for demos and tests. No real person's data is involved.

The background is random SNVs and indels with realistic proportions (transition/transversion about 2.05, more
heterozygous than homozygous calls, a few percent failing quality). On top of that a fixed list of real,
well-known variants is planted so every part of the pipeline has something to find, including the awkward
cases: an indel written without trimming, a low-quality call at a pathogenic site, two variants in one gene
phased onto opposite copies, a heteroplasmic mitochondrial variant, a long run of homozygosity.

Usage: python examples/make_wgs_demo.py OUT.vcf.gz [--snvs 3500000] [--gvcf] [--female] [--seed 7]
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RCRS = ROOT / "gnomancer" / "data" / "rCRS_NC_012920.1.fasta"

LENGTHS = {"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973,
           "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718,
           "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167, "21": 46709983,
           "22": 50818468, "X": 156040895, "Y": 57227415}
ORDER = list(LENGTHS) + ["M"]
ROH = ("5", 60_000_000, 67_000_000)          # no heterozygous calls in here
THIN = ("19", 44_900_000, 44_915_000)        # in the gVCF variant, reference blocks here are too shallow to call
TI = {"A": "G", "G": "A", "C": "T", "T": "C"}

# (chrom, pos, ref, alt, GT, AD, GQ, FILTER, PS): real GRCh38 coordinates and alleles.
PLANTED = [
    ("13", 32340300, "GT", "G", "0/1", (16, 15), 99, "PASS", ""),            # BRCA2 6174delT, heterozygous
    ("7", 117559590, "ATCTT", "AT", "0/1", (15, 14), 99, "PASS", ""),        # CFTR F508del, written with a shared trailing base
    ("6", 26092913, "G", "A", "1/1", (0, 31), 93, "PASS", ""),               # HFE C282Y, two copies
    ("13", 51944145, "G", "T", "0|1", (14, 16), 99, "PASS", "51940000"),     # ATP7B H1069Q ...
    ("13", 51958333, "C", "T", "1|0", (17, 13), 99, "PASS", "51940000"),     # ... and R778L on the other copy
    ("12", 102840493, "G", "A", "0/1", (27, 3), 40, "PASS", ""),             # PAH R408W, only 10% of reads: weak
    ("3", 38551513, "G", "A", "0/1", (14, 13), 99, "LowQual", ""),           # SCN5A, failed the caller's filter
    ("19", 44908684, "T", "C", "0/1", (15, 15), 99, "PASS", ""),             # APOE rs429358: e3/e4
    ("10", 94775367, "A", "G", "0|1", (15, 15), 99, "PASS", "94700000"),     # CYP2C19*2 is defined by this site ...
    ("10", 94781859, "G", "A", "0|1", (16, 14), 99, "PASS", "94700000"),     # ... together with this one
    ("10", 94761900, "C", "T", "1|0", (15, 15), 99, "PASS", "94700000"),     # CYP2C19*17, opposite copy
    ("10", 94842866, "A", "G", "1/1", (0, 30), 90, "PASS", ""),              # the reference base here is the rare *38 allele; *1, *2 and *17 all carry G
    ("2", 233760233, "C", "CAT", "1/1", (0, 28), 84, "PASS", ""),            # UGT1A1*28, two copies
    ("7", 99672916, "T", "C", "1/1", (0, 30), 90, "PASS", ""),               # CYP3A5*3/*3
    ("16", 31096368, "C", "T", "0/1", (15, 16), 99, "PASS", ""),             # VKORC1 -1639G>A
    ("9", 133257521, "T", "TC", "0/1", (14, 15), 99, "PASS", ""),            # ABO: one functional (A/B) allele, one O
    ("2", 135851076, "G", "A", "0/1", (16, 15), 99, "PASS", ""),             # LCT -13910: lactase persistent
    ("1", 11796321, "G", "A", "0/1", (2, 2), 12, "PASS", ""),                # MTHFR C677T with 4 reads: should be a no-call
]
# Mitochondrial differences from rCRS typical of haplogroup H1, an insertion, and a heteroplasmic m.3243A>G.
MT_FIXED = [(263, "G"), (750, "G"), (1438, "G"), (3010, "A"), (4769, "G"), (8860, "G"), (15326, "G")]
MT_HETEROPLASMIC = [(3243, "G", 0.2)]


def reserved_positions() -> dict:
    """Positions the random background must not touch: a random base there would look like a finding. Covers the
    sites gnomancer reads by name and, when `gnomancer setup` has run, PharmCAT's positions and all of ClinVar."""
    import sys
    sys.path.insert(0, str(ROOT))
    from gnomancer import clinvar, pharmcat
    from gnomancer.seqsites import load_sites
    from gnomancer.seqvcf import CHROMS
    out: dict = {}
    for s in load_sites(38).values():
        out.setdefault(s.chrom, set()).update(range(s.pos - 1, s.pos + len(s.ref) + 1))
    if pharmcat.POSITIONS.exists():
        for chrom, ps in pharmcat.named_keys(pharmcat.positions()).items():
            for p in ps:
                out.setdefault(chrom, set()).update(range(p - 1, p + 8))
    if clinvar.vcf_for(38).exists():
        for k in clinvar.positions(38):
            out.setdefault(CHROMS[k >> 28], set()).add(k & ((1 << 28) - 1))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--snvs", type=int, default=3_500_000)
    ap.add_argument("--gvcf", action="store_true", help="add reference blocks between variants")
    ap.add_argument("--female", action="store_true")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rnd = random.Random(a.seed)
    total = sum(LENGTHS[c] for c in LENGTHS if c != "Y")
    rcrs = "".join(l.strip() for l in open(RCRS) if not l.startswith(">"))
    planted: dict = {}
    for chrom, pos, *rest in PLANTED:
        planted.setdefault(chrom, {})[pos] = rest
    avoid = reserved_positions()
    opener = gzip.open if a.out.endswith(".gz") else open
    with opener(a.out, "wt") as f:
        f.write("##fileformat=VCFv4.2\n##source=gnomancer synthetic whole genome (examples/make_wgs_demo.py)\n##reference=GRCh38\n")
        for c in ORDER:
            f.write(f"##contig=<ID=chr{c},length={LENGTHS.get(c, 16569)}>\n")
        f.write('##FILTER=<ID=LowQual,Description="Low quality">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End of block">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in block">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n")
        male = not a.female
        for chrom in LENGTHS:
            if chrom == "Y" and not male:
                continue
            n = int(a.snvs * LENGTHS[chrom] / total * (0.25 if chrom == "X" else 0.004 if chrom == "Y" else 1))
            positions = set(rnd.sample(range(100_000, LENGTHS[chrom] - 100_000), n))
            positions -= set(planted.get(chrom, {}))
            positions -= avoid.get(chrom, set())
            rows = sorted(positions | set(planted.get(chrom, {})))
            haploid = male and chrom == "Y"
            last_end = 0
            for pos in rows:
                if a.gvcf and pos - last_end > 1 and not haploid:
                    thin = chrom == THIN[0] and last_end < THIN[2] and pos > THIN[1]
                    dp = 3 if thin else 30
                    f.write(f"chr{chrom}\t{last_end + 1}\t.\tN\t<NON_REF>\t.\t.\tEND={pos - 1}\tGT:DP:GQ:MIN_DP\t0/0:{dp}:{9 if thin else 90}:{dp}\n")
                if pos in planted.get(chrom, {}):
                    ref, alt, gt, ad, gq, filt, ps = planted[chrom][pos]
                    fmt, val = "GT:AD:DP:GQ", f"{gt}:{ad[0]},{ad[1]}:{sum(ad)}:{gq}"
                    if ps:
                        fmt, val = fmt + ":PS", val + ":" + ps
                    f.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t.\t{fmt}\t{val}\n")
                    last_end = pos + len(ref) - 1
                    continue
                ref = rnd.choice("ACGT")
                kind = rnd.random()
                if kind < 0.10:                                   # indel
                    extra = "".join(rnd.choice("ACGT") for _ in range(rnd.choice((1, 1, 1, 2, 3, 4))))
                    ref, alt = (ref + extra, ref) if rnd.random() < 0.5 else (ref, ref + extra)
                else:
                    alt = TI[ref] if rnd.random() < 0.672 else rnd.choice([b for b in "ACGT" if b not in (ref, TI[ref])])
                in_roh = chrom == ROH[0] and ROH[1] <= pos <= ROH[2]
                x_hemi = male and chrom == "X" and 2_781_479 < pos < 155_701_383
                hom = haploid or in_roh or x_hemi or rnd.random() < 0.38
                dp = max(10, int(rnd.gauss(30, 5)))
                bad = rnd.random() < 0.03
                if haploid:
                    gt, ad = "1", f"0,{dp // 2}"
                    dp //= 2
                elif hom:
                    gt, ad = "1/1", f"0,{dp}"
                else:
                    k = int(dp * rnd.uniform(0.38, 0.62))
                    gt, ad = "0/1", f"{dp - k},{k}"
                gq = rnd.randint(3, 15) if bad else 99
                f.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{'LowQual' if bad and rnd.random() < 0.5 else 'PASS'}\t.\tGT:AD:DP:GQ\t{gt}:{ad}:{dp}:{gq}\n")
                last_end = pos + len(ref) - 1
        mt = {p: (rcrs[p - 1], alt, "1/1", "0,2000", 2000) for p, alt in MT_FIXED}
        mt[310] = (rcrs[309], rcrs[309] + "C", "1/1", "0,1900", 1900)       # the common 315.1C insertion, left-aligned
        for p, alt, share in MT_HETEROPLASMIC:
            mt[p] = (rcrs[p - 1], alt, "0/1", f"{int(2000 * (1 - share))},{int(2000 * share)}", 2000)
        for p in sorted(mt):
            ref, alt, gt, ad, dp = mt[p]
            f.write(f"chrM\t{p}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:AD:DP:GQ\t{gt}:{ad}:{dp}:99\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
