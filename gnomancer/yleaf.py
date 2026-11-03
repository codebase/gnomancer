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
