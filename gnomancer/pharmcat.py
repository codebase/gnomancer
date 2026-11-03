"""PharmCAT (PharmGKB) star-allele calling from array data.

PharmCAT wants a GRCh38 VCF covering its allele-definition positions. Array files are GRCh37 and cover only a
fraction of those positions, so we build a VCF by matching rsIDs against PharmCAT's own positions file (which
carries GRCh38 coordinates and REF/ALT), and let PharmCAT report which positions were missing. Calls made
from partial data are provisional: PharmCAT assumes reference at untyped positions.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

from .paths import java_env, tool

JAR = tool("pharmcat", "pharmcat.jar")
POSITIONS = tool("pharmcat", "pharmcat_positions.vcf")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def build_vcf(df: pd.DataFrame, out: Path, sample: str = "sample") -> dict:
    stats = dict(positions=0, typed=0, written=0, mismatched=0, indels_skipped=0, genes={})
    gt_by_rs = df.genotype.to_dict()
    header, rows = [], []
    with open(POSITIONS) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line.rstrip("\n"))
                continue
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            chrom, pos, rsid, ref, alt, info = p[0], p[1], p[2], p[3], p[4], p[7]
            stats["positions"] += 1
            gene = next((kv[3:] for kv in info.split(";") if kv.startswith("PX=")), "?")
            stats["genes"].setdefault(gene, [0, 0])[0] += 1
            g = gt_by_rs.get(rsid)
            if not isinstance(g, str):
                continue
            stats["typed"] += 1
            alts = alt.split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts) or any(c in "ID" for c in g):
                stats["indels_skipped"] += 1
                continue
            alleles = {ref, *alts}
            if not set(g) <= alleles:
                flipped = "".join(COMP.get(c, "?") for c in g)
                if set(flipped) <= alleles and not any({ref, alt} in ({"A", "T"}, {"C", "G"}) for alt in alts):
                    g = flipped
                else:
                    stats["mismatched"] += 1
                    continue
            idx = {ref: "0", **{a: str(i + 1) for i, a in enumerate(alts)}}
            gt = "/".join(idx[c] for c in g) if len(g) == 2 else idx[g[0]]
            rows.append("\t".join([chrom, pos, rsid, ref, alt, ".", "PASS", f"PX={gene}", "GT", gt]))
            stats["genes"][gene][1] += 1
            stats["written"] += 1
    with open(out, "w") as f:
        f.write("\n".join(header) + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
        f.write("\n".join(rows) + "\n")
    return stats


def run(df: pd.DataFrame, outdir: Path, sample: str = "sample") -> dict:
    res = dict(available=JAR.exists() and POSITIONS.exists(), calls=pd.DataFrame(), stats={}, error=None, html=None)
    if not res["available"]:
        res["error"] = "PharmCAT not installed; run `gnomancer setup`"
        return res
    pdir = (outdir / "pharmcat").resolve()
    pdir.mkdir(parents=True, exist_ok=True)
    vcf = pdir / f"{sample}.vcf"
    res["stats"] = build_vcf(df, vcf, sample)
    env = java_env()
    p = subprocess.run([f"{env.get('JAVA_HOME', '')}/bin/java" if env.get("JAVA_HOME") else "java", "-jar", str(JAR),
                        "-vcf", str(vcf), "-o", str(pdir), "-reporterJson", "-reporterCallsOnlyTsv", "-del"],
                       capture_output=True, text=True, env=env, cwd=str(pdir))
    if p.returncode != 0:
        res["error"] = (p.stderr or p.stdout)[-1200:]
        return res
    tsv = next(pdir.glob(f"{sample}*.report.tsv"), None)
    html = next(pdir.glob(f"{sample}*.report.html"), None)
    res["html"] = str(html) if html else None
    if not tsv:
        res["error"] = "PharmCAT produced no calls TSV: " + (p.stdout or "")[-600:]
        return res
    t = pd.read_csv(tsv, sep="\t", skiprows=1, dtype=str).fillna("")
    rows = []
    for _, r in t.iterrows():
        gene = r["Gene"]
        typed = res["stats"]["genes"].get(gene, [0, 0])
        rows.append(dict(gene=gene, diplotype=r["Source Diplotype"], phenotype=r["Phenotype"].strip(),
                         activity=r["Activity Score"], h1=r["Haplotype 1 Function"], h2=r["Haplotype 2 Function"],
                         missing_positions=(r["Missing positions"] == "yes"), typed=f"{typed[1]}/{typed[0]}"))
    res["calls"] = pd.DataFrame(rows)
    return res


def interpret(res: dict) -> str:
    if not res["available"]:
        return res["error"]
    s = res["stats"]
    base = (f"PharmCAT v3 ran on {s['written']} of its {s['positions']} allele-defining positions (the chip types {s['typed']}; "
            f"{s['indels_skipped']} indels and {s['mismatched']} allele mismatches were dropped). ")
    if res.get("error"):
        return base + "PharmCAT failed: " + res["error"]
    return base + ("Untyped positions are assumed reference, so a called diplotype can miss rare alleles; the 'typed' column "
                   "shows how much of each gene's definition was actually observed. Treat results as provisional and confirm "
                   "with a clinical pharmacogenomic test before changing any prescription.")
