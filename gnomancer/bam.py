"""Analyses that need the reads themselves: a BAM or CRAM alongside the sequencing VCF.

A VCF is the variant caller's summary. Several things a whole genome can tell you are not in it:

- whether a position with no variant was actually read (the VCF's silence is ambiguous; the reads are not)
- how many reads cover each chromosome, which gives depth and an independent check of sex
- mitochondrial variants carried by a small share of mitochondria, below what a germline caller reports
- the Y lineage read marker by marker, without leaning on the reference
- everything in later modules that depends on read depth or read structure: CYP2D6, repeat expansions, HLA, SMN1

All of it runs through samtools and bcftools from the tool environment (`gnomancer setup --wgs`). A CRAM needs the
reference it was compressed against; the bundled GRCh38 analysis set is tried, and samtools falls back to
fetching the exact sequences by checksum from the ENA reference registry if it differs.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from . import wgstools
from .seqvcf import CHR1_LENGTH, CODE, MIN_DP, norm_chrom

# Mitochondrial positions where short reads misalign around homopolymer runs, and the placeholder base at 3107:
# apparent low-level variants here are artifacts (the same sites GATK's mitochondrial pipeline blacklists).
MT_ARTIFACT = set(range(300, 318)) | set(range(513, 527)) | set(range(3105, 3110)) | set(range(16180, 16196))
MT_MIN_SHARE = 0.03         # report mitochondrial variants carried by at least this share of reads ...
MT_MIN_READS = 10           # ... supported by this many reads on each strand
MT_MIN_DEPTH = 200


class BamError(RuntimeError):
    pass


@dataclass
class Bam:
    path: Path
    cram: bool
    build: int
    contigs: dict                       # normalized name -> (name in the file, length)
    reference: Path | None
    sample: str = ""
    notes: list = field(default_factory=list)
    cov: dict = field(default_factory=dict)         # coverage(), once run
    confirm: dict = field(default_factory=dict)     # confirm_reference(), once run

    def name(self, chrom: str) -> str:
        return self.contigs[chrom][0]

    def ref_args(self, flag: str = "--reference") -> list:
        return [flag, str(self.reference)] if self.reference else []


def _run(cmd: list, group: str | None = None, **kw) -> subprocess.CompletedProcess:
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, env=wgstools.env(group), **kw)


def open_bam(path: str | Path) -> Bam:
    """Read the header: build, contig naming, sample. Raises BamError with a fixable message."""
    path = Path(path).expanduser()
    samtools = wgstools.tool_bin("samtools")
    if not samtools:
        raise BamError("reading a BAM or CRAM needs samtools; run `gnomancer setup --wgs`")
    if not path.exists():
        raise BamError(f"file not found: {path}")
    cram = path.suffix.lower() == ".cram"
    p = _run([samtools, "view", "-H", path])
    if p.returncode != 0:
        raise BamError(f"samtools could not read {path.name}: {p.stderr.strip()[-300:]}")
    contigs, build, sample = {}, 0, ""
    for line in p.stdout.splitlines():
        f = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
        if line.startswith("@SQ"):
            c = norm_chrom(f.get("SN", ""))
            if c and c not in contigs:
                contigs[c] = (f["SN"], int(f.get("LN", 0)))
                if c == "1":
                    build = CHR1_LENGTH.get(int(f.get("LN", 0)), 0)
        elif line.startswith("@RG") and not sample:
            sample = f.get("SM", "")
    if build not in (37, 38):
        raise BamError(f"{path.name} is not aligned to GRCh37 or GRCh38 (chromosome 1 length is not recognised)")
    index = [path.with_suffix(path.suffix + e) for e in (".bai", ".crai", ".csi")] + [path.with_suffix(".bai"), path.with_suffix(".crai")]
    if not any(i.exists() for i in index):
        print(f"      indexing {path.name} (once) ...", flush=True)
        q = _run([samtools, "index", "-@", "4", path])
        if q.returncode != 0:
            raise BamError(f"could not index {path.name}: {q.stderr.strip()[-300:]}")
    b = Bam(path=path, cram=cram, build=build, contigs=contigs, reference=wgstools.reference(build), sample=sample)
    if not b.reference:
        b.notes.append(f"No GRCh{build} reference genome is set up, so steps that need one were skipped" +
                       ("; run `gnomancer setup --wgs`." if build == 38 else "; set GNOMANCER_REF_37 to an indexed FASTA."))
    return b


def coverage(b: Bam) -> dict:
    """Mean depth per chromosome from the index (reads x read length / chromosome length), and sex from X and Y.

    This costs nothing: the index already counts reads per contig. It slightly overstates depth (the denominator
    includes unsequenceable stretches), which is fine for 'is this 30x or 3x' and for the X:autosome ratio.
    """
    samtools = wgstools.tool_bin("samtools")
    p = _run([samtools, "idxstats", b.path])
    if p.returncode != 0:
        raise BamError(p.stderr.strip()[-300:])
    q = _run(["sh", "-c", f'"{samtools}" view {" ".join(b.ref_args("-T"))} "{b.path}" "{b.name("1")}:10000000-10100000" | head -2000 | cut -f10'])
    lens = [len(x) for x in q.stdout.split()]
    read_len = round(sum(lens) / len(lens)) if lens else 150
    depth, mapped_total, unmapped = {}, 0, 0
    for line in p.stdout.splitlines():
        name, length, mapped, un = line.split("\t")
        c = norm_chrom(name)
        mapped_total += int(mapped)
        unmapped += int(un)
        if c and int(length):
            depth[c] = int(mapped) * read_len / int(length)
    auto = [depth[c] for c in map(str, range(1, 23)) if c in depth]
    mean = sum(auto) / len(auto) if auto else 0.0
    x, y = depth.get("X", 0.0) / mean if mean else 0, depth.get("Y", 0.0) / mean if mean else 0
    sex = "Male" if (x < 0.75 and y > 0.10) else "Female" if (x > 0.85 and y < 0.05) else "ambiguous"
    return dict(mean_depth=mean, read_length=read_len, per_chrom=depth, x_ratio=x, y_ratio=y, sex=sex, mt_depth=depth.get("MT", 0.0),
                mapped_reads=mapped_total, unmapped_reads=unmapped)


def pileup_sites(b: Bam, sites: list) -> dict:
    """Allele counts straight from the reads at a list of (chrom, pos): {(chrom, pos): (depth, [ref, alt...], [counts])}.

    bcftools mpileup with -R jumps through the index, so a thousand sites take seconds. Reads with mapping quality
    under 20 and bases under quality 20 are ignored, the same bar most callers use.
    """
    bcftools = wgstools.tool_bin("bcftools")
    if not bcftools or not b.reference:
        return {}
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        for chrom, pos in sorted(set(sites), key=lambda s: (CODE[s[0]], s[1])):
            if chrom in b.contigs:
                f.write(f"{b.name(chrom)}\t{pos}\n")
        regions = f.name
    try:
        p = _run([bcftools, "mpileup", "-f", b.reference, "-R", regions, "-a", "AD", "-q", "20", "-Q", "20", "-d", "1000", "-I", "--no-BAQ", b.path])
    finally:
        Path(regions).unlink(missing_ok=True)
    if p.returncode != 0:
        raise BamError(f"bcftools mpileup failed: {p.stderr.strip()[-300:]}")
    out = {}
    for line in p.stdout.splitlines():
        if line[0] == "#":
            continue
        c = line.split("\t")
        chrom = norm_chrom(c[0])
        fmt = c[8].split(":")
        ad = [int(x) if x.isdigit() else 0 for x in c[9].split(":")[fmt.index("AD")].split(",")] if "AD" in fmt else []
        alleles = [c[3].upper()] + [a for a in c[4].split(",") if a != "<*>"]
        out[(chrom, int(c[1]))] = (sum(ad), alleles, ad[:len(alleles)])
    return out


def confirm_reference(b: Bam, d, positions: dict) -> dict:
    """Check, in the reads, every named position the VCF is silent about, and record the answer on the SeqData.

    After this the genotyper treats the file like a gVCF at those positions: 'reference' means reads were seen
    and matched, a thinly covered position becomes a no-call, and a position where the reads plainly carry a
    variant the VCF does not list is a no-call too (and counted, because it should not happen).
    """
    from .seqvcf import key
    todo = [(c, p) for c, ps in positions.items() for p in ps
            if key(c, p) not in d.calls and key(c, p) not in d.spanning and not (c == "Y" and d.sex == "Female")]
    piles = pileup_sites(b, todo)
    stats = dict(checked=len(todo), confirmed=0, thin=0, disagree=0, examples=[])
    if not piles and todo:
        return stats
    for chrom, pos in todo:
        depth, alleles, counts = piles.get((chrom, pos), (0, [], []))
        ref = counts[0] if counts else 0
        need = MIN_DP if not (chrom in ("X", "Y") and d.sex == "Male") else max(MIN_DP // 2, 3)
        if depth < need:
            d.covered[key(chrom, pos)] = (depth, None, False)
            stats["thin"] += 1
        elif depth and (depth - ref) / depth >= 0.20:
            d.covered[key(chrom, pos)] = (depth, None, False)
            stats["disagree"] += 1
            if len(stats["examples"]) < 5:
                stats["examples"].append(f"chr{chrom}:{pos}")
        else:
            d.covered[key(chrom, pos)] = (depth, None, True)
            stats["confirmed"] += 1
    d.reads_checked = True
    return stats


def flags(b: Bam, d) -> list:
    """Quality flags that only the reads can raise."""
    out, cov, conf = [], b.cov, b.confirm
    if cov["mean_depth"] < 15:
        out.append(f"Mean depth in the BAM/CRAM is {cov['mean_depth']:.0f}x; below about 15x heterozygous variants are often missed.")
    if cov["sex"] in ("Male", "Female") and d.sex in ("Male", "Female") and cov["sex"] != d.sex:
        out.append(f"Sex from read depth ({cov['sex']}) disagrees with sex from the VCF ({d.sex}): the two files may not be from the same person.")
    if conf["checked"] and conf["disagree"] / conf["checked"] > 0.02:
        out.append(f"At {conf['disagree']} of {conf['checked']} positions the VCF is silent about, the reads show a variant (for example {', '.join(conf['examples'])}): "
                   "the BAM/CRAM and the VCF may not be from the same sample.")
    return out


def merge_mito(d, het: dict) -> int:
    """Add read-level mitochondrial variants the VCF does not list to the call set, so the ClinVar scan and the
    haplogroup step see them with their read share. Returns how many were added."""
    from .seqvcf import Call, key
    n = 0
    for v in het.get("variants") or []:
        k = key("MT", v["pos"])
        if any(c.ref == v["ref"] and c.alt == v["alt"] for c in d.calls.get(k, [])):
            continue
        d.calls.setdefault(k, []).append(Call("MT", v["pos"], v["ref"], v["alt"], 1, 2, v["depth"], None, v["share"], "PASS", True, ""))
        n += 1
    return n


def mito_heteroplasmy(b: Bam, known: set | None = None) -> dict:
    """Mitochondrial variants by read share, from a pileup of chrM.

    A germline caller treats chrM as diploid and reports roughly 20-100% shares. The mitochondrion is covered
    thousands of times in a whole genome, so smaller shares are measurable. To keep sequencing errors and reads
    from nuclear copies of mitochondrial DNA out, a variant needs at least 3% of reads, 10 supporting reads on
    each strand, and 200 reads of depth; indels are left to the main VCF.
    """
    bcftools = wgstools.tool_bin("bcftools")
    res = dict(available=bool(bcftools and b.reference and "MT" in b.contigs), variants=[], depth=None, error=None)
    if not res["available"]:
        res["error"] = "needs bcftools and the reference genome"
        return res
    if b.contigs["MT"][1] != 16569:
        res["error"] = f"the mitochondrial contig is {b.contigs['MT'][1]} bases, not rCRS"
        return res
    p = _run([bcftools, "mpileup", "-f", b.reference, "-r", b.name("MT"), "-a", "AD,ADF,ADR", "-q", "20", "-Q", "20", "-d", "1000000", "-I", "--no-BAQ", b.path])
    if p.returncode != 0:
        res["error"] = p.stderr.strip()[-300:]
        return res
    depths = []
    for line in p.stdout.splitlines():
        if line[0] == "#":
            continue
        c = line.split("\t")
        fmt, smp = c[8].split(":"), c[9].split(":")
        get = lambda k: [int(x) if x.isdigit() else 0 for x in smp[fmt.index(k)].split(",")]  # noqa: E731
        ad, adf, adr = get("AD"), get("ADF"), get("ADR")
        total = sum(ad)
        depths.append(total)
        alts = c[4].split(",")
        if int(c[1]) in MT_ARTIFACT or c[3].upper() not in "ACGT":
            continue
        for i, alt in enumerate(alts, start=1):
            if alt == "<*>" or i >= len(ad) or total < MT_MIN_DEPTH:
                continue
            share = ad[i] / total
            if share >= MT_MIN_SHARE and adf[i] >= MT_MIN_READS and adr[i] >= MT_MIN_READS:
                res["variants"].append(dict(pos=int(c[1]), ref=c[3].upper(), alt=alt, share=round(share, 4), depth=total,
                                            forward=adf[i], reverse=adr[i], homoplasmic=share >= 0.95))
    res["depth"] = sorted(depths)[len(depths) // 2] if depths else 0
    return res


def yleaf(b: Bam, outdir: Path, sample: str = "sample") -> dict:
    """Y haplogroup with Yleaf reading the BAM or CRAM directly: every tree marker is genotyped from the reads."""
    from .yleaf import TREES, _bin, parse_combined
    res = dict(haplogroup=None, ycc="", n_y_snps=0, path=None, tool="Yleaf", error=None, yleaf=dict(trees={}, n_y_snps=0), basis="")
    exe = _bin()
    if not exe:
        res["error"] = "Yleaf not installed (`gnomancer setup` installs it when python3.11 exists)"
        return res
    ydir = (outdir / "yleaf").resolve()
    ydir.mkdir(parents=True, exist_ok=True)
    rg = "hg38" if b.build == 38 else "hg19"
    cmd = [exe, "-cram" if b.cram else "-bam", str(b.path), "-rg", rg, "-o", str(ydir / "out_reads"), "-tree", *TREES, "--no-update-check", "-force", "-t", "4"]
    if b.reference:
        cmd += ["--ref-fasta", str(b.reference)] + (["-cr", str(b.reference)] if b.cram else [])
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ydir), env=wgstools.env())
    comb = next(iter(sorted((ydir / "out_reads").rglob("hg_prediction_combined.hg"))), None) if (ydir / "out_reads").exists() else None
    if comb is None:
        res["error"] = (p.stderr or p.stdout)[-600:] or "Yleaf produced no prediction"
        return res
    trees = parse_combined(comb)
    res["yleaf"]["trees"] = trees
    best = next((trees[t] for t in ("yfull", "ftdna", "isogg") if trees.get(t, {}).get("haplogroup") not in (None, "", "NA")), None)
    if not best:
        res["error"] = "Yleaf could not place the sample on any tree"
        return res
    res["haplogroup"] = best["haplogroup"].split("*")[0] or best["haplogroup"]
    res["ycc"] = f"via {best['marker']}" if best.get("marker") else ""
    try:
        res["n_y_snps"] = res["yleaf"]["n_y_snps"] = int(float(best.get("valid_markers") or 0))
    except ValueError:
        pass
    res["basis"] = "Yleaf, reading every tree marker directly from the Y-chromosome reads in the BAM/CRAM."
    return res
