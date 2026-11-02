"""Y-chromosome (yhaplo) and mitochondrial (Haplogrep3) haplogroup calling."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from .paths import DATA, tool

RCRS = DATA / "rCRS_NC_012920.1.fasta"


def _bin(name: str) -> str | None:
    cand = Path(sys.executable).parent / name
    if cand.exists():
        return str(cand)
    return shutil.which(name)


def call_y(df: pd.DataFrame, outdir: Path, sample: str = "sample") -> dict:
    """Write a sample-major genos.txt and run yhaplo. Returns dict with call and path."""
    y = df[(df.chrom == "Y") & df.genotype.notna()].sort_values("pos")
    rows = []
    for rsid, r in y.iterrows():
        g = r.genotype
        if len(g) == 2:
            if g[0] != g[1]:
                continue
            g = g[0]
        if g in "ACGT":
            rows.append((int(r.pos), g))
    res = dict(n_y_snps=len(rows), haplogroup=None, ycc=None, path=None, tool="yhaplo", error=None)
    if len(rows) < 100:
        res["error"] = f"only {len(rows)} usable Y SNPs; skipped"
        return res
    exe = _bin("yhaplo")
    if not exe:
        res["error"] = "yhaplo not installed (pip install git+https://github.com/23andMe/yhaplo)"
        return res
    ydir = (outdir / "yhaplo").resolve()
    ydir.mkdir(parents=True, exist_ok=True)
    genos = ydir / f"{sample}.genos.txt"
    with open(genos, "w") as f:
        f.write("ID\t" + "\t".join(str(p) for p, _ in rows) + "\n")
        f.write(sample + "\t" + "\t".join(g for _, g in rows) + "\n")
    p = subprocess.run([exe, "-i", str(genos), "-o", str(ydir), "--haplogroup_paths", "--der_snps_detail",
                        "--anc_der_counts"], capture_output=True, text=True)
    hg = ydir / f"haplogroups.{sample}.txt"
    if p.returncode != 0 or not hg.exists():
        res["error"] = (p.stderr or p.stdout)[-800:]
        return res
    parts = hg.read_text().split()
    # columns: id, short hg, hg-snp form, YCC long form
    res["haplogroup"], res["ycc"] = parts[2], parts[3]
    path_file = ydir / f"paths.{sample}.txt"
    res["path"] = path_file.read_text().split("|", 1)[1].strip() if path_file.exists() else None
    if res["haplogroup"] == "A":
        res["error"] = "assigned to root (no derived SNPs found): likely a female sample or no usable Y data"
    return res


def _haplogrep() -> str | None:
    env = os.environ.get("HAPLOGREP3")
    if env and Path(env).exists():
        return env
    cand = tool("haplogrep3", "haplogrep3")
    if cand.exists():
        return str(cand)
    return shutil.which("haplogrep3")


def _java_env() -> dict:
    env = dict(os.environ)
    for jh in ("/opt/homebrew/opt/openjdk", "/usr/local/opt/openjdk", "/usr/lib/jvm/default-java"):
        if Path(jh, "bin", "java").exists() and "JAVA_HOME" not in env:
            env["JAVA_HOME"] = jh
            env["PATH"] = f"{jh}/bin:" + env.get("PATH", "")
    return env


def call_mt(df: pd.DataFrame, outdir: Path, build_original: int, sample: str = "sample",
            tree: str = "phylotree-rcrs@17.3") -> dict:
    res = dict(n_mt_snps=0, haplogroup=None, quality=None, diffs=[], indels=[], tool="haplogrep3", error=None)
    if build_original == 36:
        res["error"] = "build 36 MT coordinates are not rCRS; skipped"
        return res
    mt = df[(df.chrom == "MT") & df.genotype.notna()].sort_values("pos")
    res["n_mt_snps"] = len(mt)
    if mt.empty:
        res["error"] = "no mitochondrial SNPs in file"
        return res
    ref = "".join(l.strip() for l in open(RCRS) if not l.startswith(">"))
    diffs, indels, positions = [], [], []
    for _, r in mt.iterrows():
        g = r.genotype
        a = g[0]
        if len(g) == 2 and g[0] != g[1]:
            continue  # heteroplasmic/ambiguous: skip
        pos = int(r.pos)
        if a in "ID":
            indels.append(f"{pos}{a}")
            continue
        positions.append(pos)
        if pos <= len(ref) and ref[pos - 1] != a:
            diffs.append(f"{pos}{a}")
    res["diffs"], res["indels"] = diffs, indels
    exe = _haplogrep()
    if not exe:
        res["error"] = "haplogrep3 not found; run `gnomancer setup` or set HAPLOGREP3"
        return res
    mdir = (outdir / "haplogrep3").resolve()
    mdir.mkdir(parents=True, exist_ok=True)
    hsd = mdir / f"{sample}.hsd"
    with open(hsd, "w") as f:
        f.write("SampleID\tRange\tHaplogroup\tPolymorphisms\n")
        f.write(f"{sample}\t{';'.join(map(str, positions))}\t?\t" + "\t".join(diffs) + "\n")
    out = mdir / f"{sample}.haplogrep3.tsv"
    p = subprocess.run([exe, "classify", "--tree", tree, "--in", str(hsd), "--out", str(out), "--extend-report"],
                       capture_output=True, text=True, env=_java_env(), cwd=str(Path(exe).parent))
    if p.returncode != 0 or not out.exists():
        res["error"] = (p.stderr or p.stdout)[-800:]
        return res
    t = pd.read_csv(out, sep="\t")
    res["haplogroup"] = str(t.iloc[0]["Haplogroup"])
    res["quality"] = float(t.iloc[0]["Quality"])
    return res
