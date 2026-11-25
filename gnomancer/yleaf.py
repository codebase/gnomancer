"""Second-opinion Y haplogroup via Yleaf (Erasmus MC) on current trees (YFull v14, FTDNA, ISOGG).

yhaplo is fast and dependable but its tree is ISOGG 2016. Yleaf accepts PLINK array input and predicts on
up-to-date trees, which often resolves one or two branches deeper. It needs its own venv (Python <=3.11) and
downloads the hg19 reference (~900 MB) once on first use.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

from .paths import ancestry_bin, local_bin

TREES = ("yfull", "ftdna", "isogg")


def _bin() -> str | None:
    c = ancestry_bin("Yleaf")
    if c:
        return str(c)
    return local_bin("Yleaf")


def write_plink(df: pd.DataFrame, prefix: Path, sample: str = "sample") -> int:
    y = df[(df.chrom == "Y") & df.genotype.notna()].sort_values("pos")
    rows = []
    for rsid, r in y.iterrows():
        g = r.genotype
        if g[0] in "ACGT" and (len(g) == 1 or g[0] == g[1]):
            rows.append((rsid, int(r.pos), g[0]))
    with open(prefix.with_suffix(".map"), "w") as f:
        for rs, pos, _ in rows:
            f.write(f"Y\t{rs}\t0\t{pos}\n")
    with open(prefix.with_suffix(".ped"), "w") as f:
        f.write(f"FAM1\t{sample}\t0\t0\t1\t-9\t" + "\t".join(f"{a}\t{a}" for _, _, a in rows) + "\n")
    return len(rows)


def parse_combined(path: Path) -> dict:
    t = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    r = t.iloc[0]
    return {tree: dict(haplogroup=r.get(f"{tree}_Hg", ""), marker=r.get(f"{tree}_Hg_marker", ""),
                       valid_markers=r.get(f"{tree}_Valid_markers", ""), qc=r.get(f"{tree}_QC", ""))
            for tree in TREES if f"{tree}_Hg" in t.columns}


def run(df: pd.DataFrame, outdir: Path, sample: str = "sample") -> dict:
    res = dict(available=bool(_bin()), n_y_snps=0, trees={}, error=None)
    exe = _bin()
    if not exe:
        res["error"] = "Yleaf not installed (optional; `gnomancer setup` installs it when python3.11 exists)"
        return res
    ydir = (outdir / "yleaf").resolve()
    ydir.mkdir(parents=True, exist_ok=True)
    res["n_y_snps"] = write_plink(df, ydir / sample, sample)
    if res["n_y_snps"] < 100:
        res["error"] = f"only {res['n_y_snps']} usable Y SNPs"
        return res
    p = subprocess.run([exe, "-plink", str(ydir / f"{sample}.ped"), "-rg", "hg19", "-o", str(ydir / "out"), "-tree", *TREES,
                        "--no-update-check", "-force"], capture_output=True, text=True, cwd=str(ydir))
    comb = ydir / "out" / "hg_prediction_combined.hg"
    if not comb.exists():
        res["error"] = (p.stderr or p.stdout)[-800:]
        return res
    res["trees"] = parse_combined(comb)
    return res


def run_seq(d, outdir: Path, whole_genome: bool = True, sample: str = "sample") -> dict:
    """Y haplogroup from a sequencing VCF, with Yleaf's VCF mode on current trees.

    A VCF lists only the positions where this Y differs from the reference Y, and many tree markers are positions
    where the reference itself carries the derived allele. Yleaf handles that by reading the reference base at
    every marker the VCF does not list (right for a whole genome, where absence means 'read and matched'); for
    partial data that inference is switched off and only the listed positions are used.

    Needs bcftools (from `gnomancer setup --wgs`) and the reference genome for the file's build.
    """
    import gzip
    from . import wgstools
    res = dict(haplogroup=None, ycc="", n_y_snps=0, path=None, tool="Yleaf", error=None, yleaf=dict(trees={}, n_y_snps=0), basis="")
    exe = _bin()
    if not exe:
        res["error"] = "Yleaf not installed (`gnomancer setup` installs it when python3.11 exists)"
        return res
    if not wgstools.have("bcftools"):
        res["error"] = "the Y call from sequencing data needs bcftools; run `gnomancer setup --wgs`"
        return res
    calls = sorted((c for k, cs in d.calls.items() if (k >> 28) == 23 for c in cs
                    if c.ok and len(c.ref) == 1 and len(c.alt) == 1 and c.copies == c.ploidy), key=lambda c: c.pos)
    res["n_y_snps"] = res["yleaf"]["n_y_snps"] = len(calls)
    if len(calls) < 20:
        res["error"] = f"only {len(calls)} usable Y variants in the file"
        return res
    ydir = (outdir / "yleaf").resolve()
    ydir.mkdir(parents=True, exist_ok=True)
    plain = ydir / f"{sample}.chrY.vcf"
    with open(plain, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chrY>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for c in calls:
            f.write(f"chrY\t{c.pos}\t.\t{c.ref}\t{c.alt}\t.\tPASS\t.\tGT\t1\n")
    vcf = ydir / f"{sample}.chrY.vcf.gz"
    bgzip = wgstools.tool_bin("bgzip")
    if bgzip:
        subprocess.run([bgzip, "-f", str(plain)], check=True)
    else:
        with open(plain, "rb") as src, gzip.open(vcf, "wb") as dst:
            dst.write(src.read())
        plain.unlink()
    rg = "hg38" if d.build == 38 else "hg19"
    cmd = [exe, "-vcf", str(vcf), "-rg", rg, "-o", str(ydir / "out"), "-tree", *TREES, "--no-update-check", "-force"]
    ref = wgstools.reference(d.build)
    if ref:
        cmd += ["--ref-fasta", str(ref)]
    if not whole_genome:
        cmd.append("--no-ref-inference")
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ydir), env=wgstools.env())
    comb = next(iter(sorted((ydir / "out").rglob("hg_prediction_combined.hg"))), None) if (ydir / "out").exists() else None
    if comb is None:
        res["error"] = (p.stderr or p.stdout)[-800:] or "Yleaf produced no prediction"
        return res
    trees = parse_combined(comb)
    res["yleaf"]["trees"] = trees
    best = next((trees[t] for t in ("yfull", "ftdna", "isogg") if trees.get(t, {}).get("haplogroup") not in (None, "", "NA")), None)
    if not best:
        res["error"] = "Yleaf could not place the sample on any tree"
        return res
    res["haplogroup"] = best["haplogroup"].split("*")[0] or best["haplogroup"]
    res["ycc"] = f"via {best['marker']}" if best.get("marker") else ""
    res["basis"] = (f"Yleaf, from {len(calls):,} Y variants in the VCF" +
                    ("; marker positions the VCF does not list were read from the reference, the standard treatment of a whole genome. "
                     "Running Yleaf on the BAM or CRAM reads every marker directly." if whole_genome else "; only positions listed in the file were used."))
    return res
