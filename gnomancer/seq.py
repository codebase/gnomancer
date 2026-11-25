"""Load a sequencing VCF into the shapes the rest of gnomancer works with.

`load_seq` makes one pass over the file (see seqvcf) and returns a `Loaded` whose table holds the named sites
(marker panel, ancestry SNPs) genotyped by position, with the full per-call data kept on `Loaded.seq` for the
ClinVar scan, PharmCAT, the mitochondrial and Y calls and the quality summary.
"""
from __future__ import annotations

import sys
from pathlib import Path

from . import clinvar as clinvar_mod
from . import pharmcat as pharmcat_mod
from .load import Loaded
from .seqqc import is_whole_genome
from .seqsites import ASSUMED, load_sites, named_positions, view
from .seqvcf import Sniff, key, read, sniff


def detect(path: str) -> Sniff | None:
    """The sniff result if `path` is a sequencing VCF, else None (array exports and non-VCF files)."""
    name = str(path).lower()
    if not any(name.endswith(e) for e in (".vcf", ".vcf.gz", ".vcf.bgz", ".gvcf", ".gvcf.gz", ".g.vcf.gz", ".zip", ".gz")):
        return None
    s = sniff(path)
    return s if s.is_sequencing else None


def load_seq(path: str, sn: Sniff, build: int | None = None, use_clinvar: bool = True, use_pharmcat: bool = True,
             quiet: bool = False) -> Loaded:
    b = build or sn.build
    if b not in (37, 38):
        raise SystemExit("Could not tell which genome build this VCF uses (no contig lengths or reference name in the header). "
                         "Re-run with --build 38 or --build 37.")
    sites = load_sites(b)
    named = named_positions(sites.values())
    if use_pharmcat and b == 38 and pharmcat_mod.POSITIONS.exists():
        for chrom, ps in pharmcat_mod.named_keys(pharmcat_mod.positions()).items():
            named[chrom] = sorted(set(named.get(chrom, [])) | set(ps))
    if b == 38:
        from . import pgs
        for chrom, ps in pgs.positions().items():
            named[chrom] = sorted(set(named.get(chrom, [])) | set(ps))
    want = {key(c, p) for c, ps in named.items() for p in ps}
    if use_clinvar and clinvar_mod.vcf_for(b).exists():
        want |= clinvar_mod.positions(b)
    progress = None if quiet else (lambda n: print(f"      {n / 1e6:.0f} million records ...", file=sys.stderr))
    d = read(path, want, named, sniffed=sn, build=b, progress=progress)
    d.sites, d.named = sites, named
    whole = is_whole_genome(d)
    df = view(d, sites, whole)
    notes = list(d.notes)
    if build and sn.build and build != sn.build:
        notes.append(f"--build {build} overrides the header, which says GRCh{sn.build}.")
    if d.gvcf:
        notes.append("The file is a gVCF: positions matching the reference are read from its reference blocks, not assumed.")
    elif whole:
        n = int((df.basis == ASSUMED).sum())
        notes.append(f"This VCF lists only differences from the reference. {n} of the {len(df)} named sites are absent from it and were "
                     "taken as homozygous reference, the standard reading of a whole genome; a gVCF would confirm each one.")
    kind = "wgs" if whole else "partial"
    source = ("Whole-genome sequencing" if whole else "Sequencing (exome or panel)") + (f", {sn.caller}" if sn.caller else "")
    return Loaded(df=df, source=source, build=b, build_detected=bool(sn.build) and not build, build_original=b, remapped=False,
                  chip="", chip_version="", sex=d.sex if d.sex in ("Male", "Female") else "", count=int(d.stats["passing_alleles"]),
                  notes=notes, kind=kind, seq=d)


def refresh(loaded: Loaded) -> None:
    """Rebuild the named-site table after the reads were consulted (bam.confirm_reference), and fix the note about inference."""
    d = loaded.seq
    loaded.df = view(d, d.sites, is_whole_genome(d))
    loaded.notes = [n for n in loaded.notes if not n.startswith("This VCF lists only differences")]


def ancestry_file(loaded: Loaded, outdir) -> Path:
    """The ancestry-informative SNPs as a small array-style file, which is what ezancestry reads.

    Positions are written on GRCh37 (the build that reader expects); genotypes are forward-strand letters and do
    not depend on the build. Sites without a confident genotype are left out."""
    from .paths import DATA
    wanted = set((DATA / "aisnps.txt").read_text().split())
    pos37 = load_sites(37)
    path = Path(outdir) / "ancestry-sites.txt"
    with open(path, "w") as f:
        f.write("# rsid\tchromosome\tposition\tgenotype\n")     # the generic layout the reader recognises; no other comment lines
        for rsid, s in loaded.seq.sites.items():
            g = loaded.df.genotype.get(rsid)
            if rsid in wanted and isinstance(g, str) and rsid in pos37 and not s.indel:
                f.write(f"{rsid}\t{pos37[rsid].chrom}\t{pos37[rsid].pos}\t{g if len(g) == 2 else g * 2}\n")
    return path


def merge_mito_vcf(d, path: str) -> int:
    """Add variants from a mitochondrial heteroplasmy VCF (mity, Mutect2 mitochondrial mode, DRAGEN) to the call set.

    These callers report the share of reads per variant, as FORMAT AF / HF / VAF or as allele depths; the main
    VCF's diploid caller usually misses anything under about 20%. Only records the caller passed are used."""
    from .seqvcf import Call, norm_chrom, open_text, trim
    added = 0
    with open_text(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10 or norm_chrom(c[0]) != "MT" or c[6] not in ("PASS", "."):
                continue
            fmt, smp = c[8].split(":"), c[9].split(":")
            get = lambda k: smp[fmt.index(k)] if k in fmt and fmt.index(k) < len(smp) else ""  # noqa: E731
            ad = [int(x) for x in get("AD").split(",") if x.isdigit()]
            dp = int(get("DP")) if get("DP").isdigit() else (sum(ad) or None)
            for i, alt in enumerate(c[4].upper().split(","), start=1):
                if alt in ("<NON_REF>", "<*>", "*", "."):
                    continue
                share = None
                for k in ("AF", "HF", "VAF"):
                    vals = get(k).split(",")
                    if len(vals) >= i and vals[i - 1].replace(".", "", 1).isdigit():
                        share = float(vals[i - 1])
                        break
                if share is None and len(ad) > i and sum(ad):
                    share = ad[i] / sum(ad)
                if share is None or share <= 0:
                    continue
                pos, ref, a2 = trim(int(c[1]), c[3].upper(), alt)
                k = key("MT", pos)
                if any(x.ref == ref and x.alt == a2 for x in d.calls.get(k, [])):
                    continue
                d.calls.setdefault(k, []).append(Call("MT", pos, ref, a2, 1, 2, dp, None, share, c[6], True, ""))
                added += 1
    return added
