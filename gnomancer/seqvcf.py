"""Read a sequencing VCF (whole genome or exome) in one streaming pass.

An array export names a fixed list of sites and gives a genotype at each. A sequencing VCF is a different kind
of file, and reading it as if it were an array throws most of it away:

- sites are identified by position and alleles; the ID column is often '.'
- insertions and deletions are written out as sequence, and matter (most severe pathogenic variants are indels)
- a plain VCF lists only the places where the sample differs from the reference, so a site that is absent is
  usually homozygous reference, not untested; a gVCF says so explicitly with reference blocks
- every call carries evidence (FILTER, depth, genotype quality, allele depths) that should gate what is reported
- coordinates are usually GRCh38, with 'chr' contig names and 'chrM' for the mitochondrion

The file can hold five million records, so nothing here keeps all of them. One pass collects: the calls at
positions the caller asked for (ClinVar positions, the sites gnomancer genotypes by name), every call on Y and
MT, reference-block coverage for the named sites, and running quality statistics.

Alleles are normalized per alternate allele by trimming shared bases (suffix first, then prefix), which is what
`bcftools norm -m-` does short of left-alignment. Variant callers already left-align indels, and so does the
ClinVar VCF, so trimmed alleles compare exactly.
"""
from __future__ import annotations

import gzip
import io
import zipfile
from array import array
from dataclasses import dataclass, field
from pathlib import Path

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
CODE = {c: i for i, c in enumerate(CHROMS)}
AUTOSOMES = set(CHROMS[:22])
CHR1_LENGTH = {248956422: 38, 249250621: 37, 247249719: 36}
# Pseudo-autosomal regions of X (1-based, inclusive): diploid in males, so excluded from sex inference.
PAR_X = {38: ((10001, 2781479), (155701383, 156030895)), 37: ((60001, 2699520), (154931044, 155260560))}
# Outside these a male Y call is haploid; inside, callers place the reads on X.
REF_BLOCK_ALTS = (".", "<NON_REF>", "<*>")

MIN_DP = 8          # reads covering the site
MIN_GQ = 20         # phred genotype quality: 1 in 100 chance the genotype is wrong
HET_AB = (0.20, 0.80)   # share of reads carrying the alternate allele at a heterozygous call
HOM_AB = 0.85           # minimum share at a homozygous-alternate call


def key(chrom: str, pos: int) -> int:
    """One integer per genomic position (chromosome index in the high bits), cheap to keep millions of in a set."""
    return (CODE[chrom] << 28) | pos


def trim(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Drop bases shared by both alleles: from the right first (keeps the left anchor), then from the left."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def norm_chrom(name: str) -> str | None:
    """'chr7' -> '7', 'chrM'/'M' -> 'MT'; anything off the primary assembly (alts, decoys, HLA contigs) -> None."""
    c = name[3:] if name[:3].lower() == "chr" else name
    c = "MT" if c in ("M", "MT", "m") else c.upper() if c in ("x", "y") else c
    return c if c in CODE else None


def open_text(path: str | Path):
    """Text handle for .vcf, .vcf.gz/.bgz (BGZF is multi-member gzip) or a .zip holding one of those."""
    path = Path(path)
    with open(path, "rb") as f:
        magic = f.read(4)
    if magic[:2] == b"PK":
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if n.lower().endswith((".vcf", ".vcf.gz", ".vcf.bgz")) and not n.startswith("__MACOSX")]
        if not names:
            raise SystemExit(f"No VCF inside {path.name}.")
        raw = z.open(sorted(names, key=lambda n: -z.getinfo(n).file_size)[0])
        if names[0].lower().endswith(("gz", "bgz")):
            return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8", errors="replace")
        return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    if magic[:2] == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


@dataclass
class Sniff:
    """What the header and the first records say about a VCF."""
    is_sequencing: bool = False
    build: int = 0
    build_from: str = ""
    caller: str = ""
    sample: str = ""
    n_samples: int = 0
    gvcf: bool = False              # has reference blocks, so absent-vs-reference is known, not assumed
    rsid_fill: float = 0.0          # share of records with an rsID
    formats: tuple = ()
    mt_length: int = 0              # 16569 = rCRS; 16571 = the old hg19 chrM (positions are not rCRS)
    records_seen: int = 0
    why: str = ""


CALLERS = (("DRAGEN", "dragen"), ("DeepVariant", "deepvariant"), ("GATK HaplotypeCaller", "haplotypecaller"),
           ("GATK", "gatk"), ("Strelka", "strelka"), ("bcftools", "bcftools"), ("FreeBayes", "freebayes"),
           ("Sentieon", "sentieon"), ("Clair3", "clair3"), ("Octopus", "octopus"))


def sniff(path: str | Path, limit: int = 20000) -> Sniff:
    """Decide whether a VCF came from sequencing, and read build, caller and layout from it.

    Array exports saved as VCF have an rsID on every row and nothing but GT; they belong on the array path."""
    s = Sniff()
    head_low = []
    n = ids = indels = with_depth = 0
    try:
        f = open_text(path)
    except (OSError, zipfile.BadZipFile, SystemExit):
        return s
    with f:
        first = f.readline()
        if not first.startswith("##fileformat=VCF"):
            return s
        for line in f:
            if line.startswith("##"):
                low = line.lower()
                head_low.append(low)
                if low.startswith("##contig=<"):
                    body = line[10:].rstrip(">\n")
                    kv = dict(p.split("=", 1) for p in body.split(",") if "=" in p)
                    cid, length = norm_chrom(kv.get("ID", "")), kv.get("length", "")
                    if length.isdigit():
                        if cid == "1" and int(length) in CHR1_LENGTH:
                            s.build, s.build_from = CHR1_LENGTH[int(length)], "contig lengths in the header"
                        elif cid == "MT":
                            s.mt_length = int(length)
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                s.n_samples = max(len(cols) - 9, 0)
                s.sample = cols[9] if len(cols) > 9 else ""
                continue
            c = line.split("\t", 10)
            if len(c) < 10:
                continue
            n += 1
            ids += c[2].startswith("rs")
            if c[4] in REF_BLOCK_ALTS or c[4].endswith((",<NON_REF>", ",<*>")):
                s.gvcf = True
            elif len(c[3]) != 1 or any(len(a) != 1 for a in c[4].split(",")):
                indels += 1
            fmt = c[8].split(":")
            if n == 1:
                s.formats = tuple(fmt)
            with_depth += ("DP" in fmt or "AD" in fmt or "GQ" in fmt)
            if n >= limit:
                break
    s.records_seen = n
    if not n:
        return s
    s.rsid_fill = ids / n
    header = "".join(head_low)
    for name, token in CALLERS:
        if token in header:
            s.caller = name
            break
    if not s.build:
        for line in head_low:
            if line.startswith(("##reference=", "##assembly=", "##dragencommandline", "##gatkcommandline")):
                for build, tokens in ((38, ("grch38", "hg38", "hs38", "gca_000001405.15")), (37, ("grch37", "hg19", "hs37", "b37", "human_g1k_v37"))):
                    if any(t in line for t in tokens):
                        s.build, s.build_from = build, "reference named in the header"
                        break
            if s.build:
                break
    # Sequencing output carries per-call evidence or spelled-out indels; an array export has neither.
    s.is_sequencing = bool(s.gvcf or with_depth / n > 0.5 or indels / n > 0.01 or (s.caller and s.caller != "bcftools"))
    s.why = (f"{with_depth / n:.0%} of records carry depth or quality fields, {indels / n:.1%} are indels, "
             f"{s.rsid_fill:.0%} have an rsID" + (f", caller {s.caller}" if s.caller else "") + (", reference blocks present" if s.gvcf else ""))
    return s


@dataclass(slots=True)
class Call:
    """One alternate allele the sample carries at one normalized site."""
    chrom: str
    pos: int
    ref: str
    alt: str
    copies: int             # chromosomes carrying this allele
    ploidy: int             # alleles in the genotype (2, or 1 on a haploid call)
    dp: int | None
    gq: int | None
    ab: float | None        # share of reads supporting this allele
    filt: str
    ok: bool
    why: str                # reason a call did not pass, '' when ok
    phased: bool = False
    ps: str = ""            # phase set, when the caller gives one
    other_alt: bool = False     # the other chromosome carries a different alternate allele
    hap: int = -1               # which haplotype (0/1) carries the allele, for a phased heterozygous call

    @property
    def zygosity(self) -> str:
        if self.ploidy == 1:
            return "hemizygous"
        return "homozygous" if self.copies == self.ploidy else "heterozygous"


@dataclass
class SeqData:
    path: str
    sniff: Sniff
    build: int
    calls: dict = field(default_factory=dict)          # key -> [Call] at requested positions, and all of Y and MT
    spanning: dict = field(default_factory=dict)       # key -> [Call]: a deletion or no-call covers this named site
    covered: dict = field(default_factory=dict)        # key -> (dp, gq, ok): reference block over a named site
    het_pos: dict = field(default_factory=dict)        # autosome -> array of passing heterozygous SNV positions
    hom_pos: dict = field(default_factory=dict)        # autosome -> array of passing homozygous SNV positions
    stats: dict = field(default_factory=dict)
    per_chrom: dict = field(default_factory=dict)      # chrom -> counters
    sex: str = "unknown"
    notes: list = field(default_factory=list)
    sites: dict = field(default_factory=dict)          # rsid -> Site, the named sites this pass tracked
    named: dict = field(default_factory=dict)          # chrom -> positions tracked (named sites and PharmCAT's)
    reads_checked: bool = False                        # `covered` was filled in from a BAM/CRAM for the silent positions

    @property
    def gvcf(self) -> bool:
        return self.sniff.gvcf

    @property
    def observed(self) -> bool:
        """True when a silent position has been looked at (gVCF blocks or the reads), so reference is seen, not assumed."""
        return self.sniff.gvcf or self.reads_checked


def _judge(filt_ok: bool, filt: str, dp, gq, ab, hom: bool, chrom: str) -> tuple[bool, str]:
    if not filt_ok:
        return False, f"caller filter: {filt}"
    if dp is not None and dp < MIN_DP:
        return False, f"only {dp} reads"
    if gq is not None and gq < MIN_GQ:
        return False, f"genotype quality {gq}"
    if ab is not None and chrom != "MT":     # on MT a low share of reads is heteroplasmy, not an error
        if hom and ab < HOM_AB:
            return False, f"called homozygous but only {ab:.0%} of reads carry the allele"
        if not hom and not (HET_AB[0] <= ab <= HET_AB[1]):
            return False, f"called heterozygous but {ab:.0%} of reads carry the allele"
    return True, ""


def _int(s: str):
    try:
        return int(s)
    except ValueError:
        try:
            return int(float(s))
        except ValueError:
            return None


def read(path: str | Path, want: set, named: dict, sniffed: Sniff | None = None, build: int | None = None,
         progress=None) -> SeqData:
    """Stream the VCF once.

    want:  set of key(chrom, pos) whose calls to keep (ClinVar positions plus the named sites)
    named: chrom -> sorted list of positions gnomancer genotypes by name; for these, reference blocks and
           overlapping deletions are also recorded so 'reference' and 'no call' can be told apart.
    """
    sn = sniffed or sniff(path)
    b = build or sn.build
    if b not in (37, 38):
        raise SystemExit("Could not tell which genome build this VCF uses (no contig lengths or reference name in the header). "
                         "Re-run with --build 38 or --build 37.")
    d = SeqData(path=str(path), sniff=sn, build=b)
    named_set = {key(c, p) for c, ps in named.items() for p in ps}
    ptr = {c: 0 for c in named}                      # per-chromosome cursor into `named` for reference blocks
    calls, spanning, covered = d.calls, d.spanning, d.covered
    het_pos = {c: array("I") for c in AUTOSOMES}
    hom_pos = {c: array("I") for c in AUTOSOMES}
    pc = {c: dict(snv_het=0, snv_hom=0, indel_het=0, indel_hom=0, failed=0, hap=0) for c in CHROMS}
    par = PAR_X[b]
    x_het = x_hom = x_hap = 0
    ti = tv = n_rec = n_var = n_failed = n_symbolic = n_missing = n_offprimary = n_blocks = 0
    dp_sum = dp_n = 0
    multi_warned = False
    contig: dict = {}
    fmt_cache: dict = {}
    gt_cache: dict = {}
    TI = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    with open_text(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 10)
            if len(c) < 10:
                continue
            n_rec += 1
            if progress and n_rec % 1_000_000 == 0:
                progress(n_rec)
            raw = c[0]
            chrom = contig.get(raw, 0)
            if chrom == 0:
                chrom = contig[raw] = norm_chrom(raw)
            if chrom is None:
                n_offprimary += 1
                continue
            if len(c) > 10 and not multi_warned:
                multi_warned = True
                d.notes.append(f"The VCF holds more than one sample; only the first ({sn.sample}) is analyzed.")
            pos = int(c[1])
            ref, alt, filt = c[3].upper(), c[4], c[6]
            filt_ok = filt == "PASS" or filt == "."
            fmt = c[8]
            idx = fmt_cache.get(fmt)
            if idx is None:
                keys = fmt.split(":")
                idx = fmt_cache[fmt] = tuple(keys.index(k) if k in keys else -1 for k in ("GT", "DP", "GQ", "AD", "PS", "MIN_DP"))
            i_gt, i_dp, i_gq, i_ad, i_ps, i_mindp = idx
            smp = c[9].rstrip("\n").split(":")
            gts = smp[i_gt] if 0 <= i_gt < len(smp) else "."
            g = gt_cache.get(gts)
            if g is None:
                phased = "|" in gts
                parts = gts.replace("|", "/").split("/")
                g = gt_cache[gts] = (None if any(not p.isdigit() for p in parts) else tuple(int(p) for p in parts), phased)
            gt, phased = g
            dp = _int(smp[i_dp]) if 0 <= i_dp < len(smp) else None
            gq = _int(smp[i_gq]) if 0 <= i_gq < len(smp) else None

            # Reference block (gVCF) or a homozygous-reference record: the sample was read here and matches.
            is_block = alt in REF_BLOCK_ALTS
            if not is_block and alt.endswith((",<NON_REF>", ",<*>")):
                alt = alt.rsplit(",", 1)[0]
            if is_block or (gt is not None and not any(gt)):
                n_blocks += is_block
                ps_named = named.get(chrom)
                if ps_named:
                    end = pos + len(ref) - 1
                    if is_block:
                        j = c[7].find("END=")
                        if j >= 0 and (j == 0 or c[7][j - 1] == ";"):
                            end = int(c[7][j + 4:].split(";", 1)[0])
                        if dp is None and 0 <= i_mindp < len(smp):
                            dp = _int(smp[i_mindp])
                    i = ptr[chrom]
                    while i < len(ps_named) and ps_named[i] < pos:
                        i += 1
                    ptr[chrom] = i
                    if i < len(ps_named) and ps_named[i] <= end:
                        # DeepVariant marks a candidate it judged to be reference with FILTER=RefCall: for a 0/0 record
                        # that is the caller saying 'reference', not a failed call.
                        ok = (filt_ok or filt == "RefCall") and gt is not None and (dp is None or dp >= MIN_DP) and (gq is None or gq >= MIN_GQ)
                        while i < len(ps_named) and ps_named[i] <= end:
                            covered[(CODE[chrom] << 28) | ps_named[i]] = (dp, gq, ok)
                            i += 1
                continue

            n_var += 1
            if gt is None:
                n_missing += 1
                for p in range(pos, pos + min(len(ref), 2000)):
                    k = (CODE[chrom] << 28) | p
                    if k in named_set:
                        spanning.setdefault(k, []).append(Call(chrom, pos, ref, ".", 0, 2, dp, gq, None, filt, False, "no genotype called"))
                continue
            alts = alt.upper().split(",")
            ad = None
            if 0 <= i_ad < len(smp) and smp[i_ad] not in (".", ""):
                ad = [_int(x) or 0 for x in smp[i_ad].split(",")]
                if dp is None:
                    dp = sum(ad)
            if dp is not None:
                dp_sum += dp
                dp_n += 1
            ploidy = len(gt)
            carried = sorted({a for a in gt if a})
            keep_all = chrom == "Y" or chrom == "MT"
            for a in carried:
                if a > len(alts):
                    continue
                A = alts[a - 1]
                if A == "*" or A[0] == "<" or "[" in A or "]" in A:
                    n_symbolic += 1
                    continue
                copies = gt.count(a)
                hom = copies == ploidy
                ab = None
                if ad and a < len(ad):
                    tot = sum(ad)
                    ab = ad[a] / tot if tot else None
                ok, why = _judge(filt_ok, filt, dp, gq, ab, hom, chrom)
                p2, r2, a2 = (pos, ref, A) if len(ref) == 1 or len(A) == 1 else trim(pos, ref, A)
                snv = len(r2) == 1 and len(a2) == 1
                st = pc[chrom]
                if not ok:
                    st["failed"] += 1
                    n_failed += 1
                else:
                    if ploidy == 1:
                        st["hap"] += 1
                    st[("snv_" if snv else "indel_") + ("hom" if hom else "het")] += 1
                    if snv:
                        if (r2, a2) in TI:
                            ti += 1
                        else:
                            tv += 1
                        if chrom in AUTOSOMES:
                            (hom_pos if hom else het_pos)[chrom].append(p2)
                        elif chrom == "X" and not (par[0][0] <= p2 <= par[0][1] or par[1][0] <= p2 <= par[1][1]):
                            if ploidy == 1:
                                x_hap += 1
                            elif hom:
                                x_hom += 1
                            else:
                                x_het += 1
                k = (CODE[chrom] << 28) | p2
                is_named = k in named_set
                call = None
                if keep_all or is_named or k in want:
                    call = Call(chrom, p2, r2, a2, copies, ploidy, dp, gq, ab, filt, ok, why, phased,
                                smp[i_ps] if 0 <= i_ps < len(smp) and smp[i_ps] != "." else "", len(carried) > 1,
                                gt.index(a) if phased and copies == 1 and ploidy == 2 else -1)
                    calls.setdefault(k, []).append(call)
                # A deletion removes the bases after its anchor; a named site inside it has no base to report.
                if len(r2) > 1 and named_set:
                    for p in range(p2 + 1, p2 + min(len(r2), 2000)):
                        kk = (CODE[chrom] << 28) | p
                        if kk in named_set:
                            spanning.setdefault(kk, []).append(call or Call(chrom, p2, r2, a2, copies, ploidy, dp, gq, ab, filt, ok, why, phased, "", len(carried) > 1))

    d.het_pos, d.hom_pos, d.per_chrom = het_pos, hom_pos, pc
    y_ok = sum(v for k_, v in pc["Y"].items() if k_ not in ("failed", "hap"))
    x_dip = x_het + x_hom
    if x_hap > max(50, x_dip):
        d.sex = "Male"                               # the caller wrote X as haploid, which it only does for a male
    elif x_dip >= 500:
        frac = x_het / x_dip
        d.sex = "Male" if frac < 0.25 else "Female" if frac > 0.40 else "ambiguous"
    passing = sum(v for st in pc.values() for k_, v in st.items() if k_ not in ("failed", "hap"))
    auto = [pc[c] for c in AUTOSOMES]
    het = sum(s["snv_het"] for s in auto)
    hom = sum(s["snv_hom"] for s in auto)
    d.stats = dict(records=n_rec, variant_records=n_var, reference_blocks=n_blocks, passing_alleles=passing, failed_alleles=n_failed,
                   symbolic_alleles=n_symbolic, missing_genotypes=n_missing, off_primary_records=n_offprimary,
                   snvs=sum(s["snv_het"] + s["snv_hom"] for s in pc.values()), indels=sum(s["indel_het"] + s["indel_hom"] for s in pc.values()),
                   ti_tv=(ti / tv if tv else None), het_hom_ratio=(het / hom if hom else None),
                   mean_depth_at_variants=(dp_sum / dp_n if dp_n else None),
                   x_het=x_het, x_hom=x_hom, x_haploid=x_hap, y_calls=y_ok)
    return d


def fetch(path: str | Path, chrom: str, start: int, end: int) -> list[dict]:
    """Records overlapping chrom:start-end (1-based, inclusive), for single-site lookups after the main pass.

    Uses the tabix index when one sits next to the file and tabix is available; otherwise streams the file and
    stops once it has passed the region (VCFs are sorted). Returns plain dicts: pos, ref, alts, filter, genotype
    (allele indices), depth, gq, allele depths, and whether the record is a reference block.
    """
    import subprocess
    from . import wgstools
    lines = None
    tabix = wgstools.tool_bin("tabix")
    if tabix and any(Path(str(path) + ext).exists() for ext in (".tbi", ".csi")):
        for name in (f"chr{'M' if chrom == 'MT' else chrom}", chrom):
            p = subprocess.run([tabix, str(path), f"{name}:{start}-{end}"], capture_output=True, text=True)
            if p.returncode == 0 and p.stdout:
                lines = p.stdout.splitlines()
                break
        if lines is None:
            lines = []
    if lines is None:
        lines, seen = [], False
        with open_text(path) as f:
            for line in f:
                if line[0] == "#":
                    continue
                c0, _, rest = line.partition("\t")
                if norm_chrom(c0) != chrom:
                    if seen:
                        break
                    continue
                seen = True
                pos = int(rest.split("\t", 1)[0])
                if pos > end:
                    break
                if pos + 20000 < start and "END=" not in line:
                    continue
                lines.append(line.rstrip("\n"))
    out = []
    for line in lines:
        c = line.split("\t")
        if len(c) < 10:
            continue
        pos, ref, alt = int(c[1]), c[3], c[4]
        block = alt in REF_BLOCK_ALTS
        stop = pos + len(ref) - 1
        if block:
            j = c[7].find("END=")
            if j >= 0:
                stop = int(c[7][j + 4:].split(";", 1)[0])
        if stop < start or pos > end:
            continue
        fmt, smp = c[8].split(":"), c[9].split(":")
        get = lambda k: smp[fmt.index(k)] if k in fmt and fmt.index(k) < len(smp) else None  # noqa: E731
        out.append(dict(pos=pos, end=stop, ref=ref, alts=[a for a in alt.split(",") if a not in ("<NON_REF>", "<*>")] if not block else [],
                        filter=c[6], genotype=get("GT"), depth=get("DP") or get("MIN_DP"), gq=get("GQ"), allele_depths=get("AD"), reference_block=block))
    return out
