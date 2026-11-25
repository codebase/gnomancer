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


def run(df: pd.DataFrame, outdir: Path, sample: str = "sample", seq=None, whole_genome: bool = True, outside: dict | None = None) -> dict:
    """Array input: `df`. Sequencing input: pass `seq` (a SeqData read with PharmCAT's positions tracked).

    `outside` is {gene: diplotype} from a dedicated caller (CYP2D6 from read depth); PharmCAT takes it as an
    'outside call' and assigns the phenotype, which it will not do for CYP2D6 from a VCF."""
    res = dict(available=JAR.exists() and POSITIONS.exists(), calls=pd.DataFrame(), stats={}, error=None, html=None)
    if not res["available"]:
        res["error"] = "PharmCAT not installed; run `gnomancer setup`"
        return res
    if seq is not None and seq.build != 38:
        res["error"] = "PharmCAT needs GRCh38 coordinates and this file is on GRCh37; pharmacogenomic calls were skipped."
        return res
    pdir = (outdir / "pharmcat").resolve()
    pdir.mkdir(parents=True, exist_ok=True)
    vcf = pdir / f"{sample}.vcf"
    res["stats"] = build_vcf(df, vcf, sample) if seq is None else build_vcf_seq(seq, vcf, whole_genome, sample)
    env = java_env()
    cmd = [f"{env.get('JAVA_HOME', '')}/bin/java" if env.get("JAVA_HOME") else "java", "-jar", str(JAR),
           "-vcf", str(vcf), "-o", str(pdir), "-reporterJson", "-reporterCallsOnlyTsv", "-del"]
    extra = []
    if outside:
        po = pdir / f"{sample}.outside.tsv"
        po.write_text("".join(f"{g}\t{d}\n" for g, d in outside.items()))
        extra = ["-po", str(po)]
    p = subprocess.run(cmd + extra, capture_output=True, text=True, env=env, cwd=str(pdir))
    if p.returncode != 0 and extra:       # an outside diplotype PharmCAT does not recognise must not cost the other genes
        res["outside_error"] = (p.stderr or p.stdout).strip()[-300:]
        p = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(pdir))
    elif extra:
        res["outside"] = dict(outside)
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
        if gene in (res.get("outside") or {}):
            typed = [0, 0]
        rows.append(dict(gene=gene, diplotype=r["Source Diplotype"], phenotype=r["Phenotype"].strip(),
                         activity=r["Activity Score"], h1=r["Haplotype 1 Function"], h2=r["Haplotype 2 Function"],
                         missing_positions=(r["Missing positions"] == "yes") and gene not in (res.get("outside") or {}),
                         typed=("read depth" if gene in (res.get("outside") or {}) else f"{typed[1]}/{typed[0]}")))
    res["calls"] = pd.DataFrame(rows)
    if seq is not None and not res["calls"].empty:
        res["calls"]["seq"] = True          # wording downstream: 'file' and 'read', not 'chip' and 'typed'
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


# ---------------------------------------------------------------------------------------------------------------
# Sequencing VCFs
#
# PharmCAT treats a position that is missing from its input as a no-call, and a plain VCF leaves out every
# position where the sample matches the reference. PharmCAT's own preprocessor offers `--absent-to-ref` for
# that case; this does the same thing with the evidence gnomancer has: a reference block (gVCF) is an observed
# reference call, absence from a whole-genome VCF is taken as reference, and a call that failed the quality
# gates, or anything in a partial file, stays missing. Indels are matched as alleles, not skipped, and phased
# genotypes keep their phase set so PharmCAT can use it.

def positions(path: Path | None = None) -> list[dict]:
    """PharmCAT's allele-defining positions (GRCh38): chrom, pos, rsid, ref, alts, gene."""
    out = []
    with open(path or POSITIONS) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            gene = next((kv[3:] for kv in p[7].split(";") if kv.startswith("PX=")), "?")
            out.append(dict(raw_chrom=p[0], chrom=p[0][3:] if p[0].startswith("chr") else p[0], pos=int(p[1]), rsid=p[2],
                            ref=p[3], alts=p[4].split(","), gene=gene))
    return out


def named_keys(pos_list: list[dict]) -> dict:
    """chrom -> positions to track, covering every trimmed allele start so calls and reference blocks are kept."""
    from .seqvcf import CODE, trim
    out: dict = {}
    for r in pos_list:
        if r["chrom"] not in CODE:
            continue
        s = out.setdefault(r["chrom"], set())
        s.add(r["pos"])
        for a in r["alts"]:
            s.add(trim(r["pos"], r["ref"], a)[0])
    return {c: sorted(v) for c, v in out.items()}


def build_vcf_seq(d, out: Path, whole_genome: bool, sample: str = "sample") -> dict:
    from .seqsites import expected_ploidy
    from .seqvcf import CODE, key, trim
    stats = dict(positions=0, typed=0, written=0, mismatched=0, indels_skipped=0, genes={}, called=0, reference_read=0,
                 reference_assumed=0, low_quality=0, not_covered=0, novel_alleles=0)
    header = [l.rstrip("\n") for l in open(POSITIONS) if l.startswith("##")]
    rows = []
    for r in positions():
        stats["positions"] += 1
        g = stats["genes"].setdefault(r["gene"], [0, 0])
        g[0] += 1
        if r["chrom"] not in CODE:
            continue
        ref, alts = r["ref"], list(r["alts"])
        targets = {trim(r["pos"], ref, a): i + 1 for i, a in enumerate(alts)}
        found, bad, novel = [], False, False
        for start in {t[0] for t in targets} | {r["pos"]}:
            for c in d.calls.get(key(r["chrom"], start), []):
                i = targets.get((c.pos, c.ref, c.alt))
                if i is None:
                    if not (len(ref) == 1 and c.pos == r["pos"] and len(c.ref) == 1 and len(c.alt) == 1):
                        continue
                    if c.ref != ref:                        # the file's reference base disagrees with GRCh38 here
                        stats["mismatched"] += 1
                        bad = True
                        continue
                    alts.append(c.alt)                      # an SNV PharmCAT does not list: pass it through as observed
                    i = len(alts)
                    novel = True
                if not c.ok:
                    bad = True
                found.append((i, c))
        k = key(r["chrom"], r["pos"])
        ploidy = expected_ploidy(r["chrom"], r["pos"], d.sex, d.build)
        if bad or any(c.alt == "." for c in d.spanning.get(k, [])):
            stats["low_quality"] += 1
            continue
        ps = ""
        if found:
            n = found[0][1].ploidy
            haps = [0] * n
            phased = all(c.phased for _, c in found) and n == 2
            slot = 0
            for i, c in found:
                if phased and c.hap >= 0:
                    haps[c.hap] = i
                else:
                    for _ in range(c.copies):
                        while slot < n and haps[slot]:
                            slot += 1
                        if slot < n:
                            haps[slot] = i
            if phased and len({c.ps for _, c in found}) == 1:
                ps = found[0][1].ps
            if ploidy == 1 and len(set(haps)) == 1:
                haps = haps[:1]
            if not phased:
                haps.sort()
            gt = ("|" if phased else "/").join(map(str, haps))
            stats["called"] += 1
            stats["novel_alleles"] += novel
        elif d.observed:
            cov = d.covered.get(k)
            if not cov or not cov[2]:
                stats["not_covered"] += 1
                continue
            gt = "/".join(["0"] * ploidy)
            stats["reference_read"] += 1
        elif whole_genome:
            gt = "/".join(["0"] * ploidy)
            stats["reference_assumed"] += 1
        else:
            stats["not_covered"] += 1
            continue
        fmt, val = ("GT:PS", f"{gt}:{ps}") if ps else ("GT", gt)
        rows.append("\t".join([r["raw_chrom"], str(r["pos"]), r["rsid"], ref, ",".join(alts), ".", "PASS", f"PX={r['gene']}", fmt, val]))
        g[1] += 1
        stats["written"] += 1
    stats["typed"] = stats["written"]
    with open(out, "w") as f:
        f.write("\n".join(header) + "\n")
        if not any(h.startswith("##FORMAT=<ID=PS") for h in header):
            f.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
        f.write("\n".join(rows) + "\n")
    return stats


def interpret_seq(res: dict, d) -> str:
    if not res["available"]:
        return res["error"]
    s = res["stats"]
    if d.build != 38:
        return res.get("error") or "PharmCAT needs GRCh38 coordinates; this file is on another build."
    how = (f"{s['reference_read']} were read as reference from the {'gVCF' if d.gvcf else 'BAM/CRAM'}" if d.observed else
           f"{s['reference_assumed']} are absent from the VCF and were taken as reference")
    base = (f"PharmCAT ran on {s['written']} of its {s['positions']} allele-defining positions: {s['called']} carry a variant, {how}; "
            f"{s['low_quality'] + s['not_covered']} were left out as low quality or not covered. ")
    if res.get("error"):
        return base + "PharmCAT failed: " + res["error"]
    return base + ("Insertions and deletions are included, so these calls are far more complete than an array's. Short reads still cannot "
                   "resolve CYP2D6 (deletions, duplications and hybrids with its pseudogene) from a VCF, and HLA types need a dedicated caller; "
                   "see the separate sections if a BAM or CRAM was analyzed. Confirm with a clinical pharmacogenomic test before changing any prescription.")
