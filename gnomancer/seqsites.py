"""Genotypes at the sites gnomancer reads by name (marker panel, ancestry SNPs), from a sequencing VCF.

The rest of gnomancer works on a table indexed by rsID with a genotype such as 'AG', 'TT' or 'ID'. An array file
is that table. For a sequencing VCF this module builds the same table from positions: `data/sites.tsv` says
where each rsID sits on GRCh37 and GRCh38 and what the reference and alternate alleles are there.

What a missing record means depends on the file:
- gVCF: a reference block covers the site, so homozygous reference is observed; no block means not read.
- whole-genome VCF without reference blocks: an absent site is taken as homozygous reference. At 30x nearly all
  of the genome that can be mapped is covered, so this is the standard reading, but it is an inference and the
  table says so in its `basis` column.
- partial data (an exome, a gene panel): absent means unknown, and the site is left as a no-call.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

import pandas as pd

from .paths import DATA
from .seqvcf import PAR_X, SeqData, key, trim

SITES = DATA / "sites.tsv"

READ = "read as reference"
ASSUMED = "not in the VCF: taken as reference"
CALLED = "called"


@dataclass(frozen=True)
class Site:
    rsid: str
    set: str
    chrom: str
    pos: int
    ref: str
    alts: tuple

    @property
    def indel(self) -> bool:
        return len(self.ref) != 1 or any(len(a) != 1 for a in self.alts)


def load_sites(build: int, path=SITES) -> dict:
    """rsid -> Site on the given build (sites with no position on that build are left out)."""
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
            pos, ref, alt = r[f"pos{build}"], r[f"ref{build}"], r[f"alt{build}"]
            if pos:
                out[r["rsid"]] = Site(r["rsid"], r["set"], r["chrom"], int(pos), ref, tuple(alt.split(",")) if alt else ())
    return out


def named_positions(sites) -> dict:
    """chrom -> sorted positions, the shape seqvcf.read wants for reference-block tracking."""
    out: dict = {}
    for s in sites:
        out.setdefault(s.chrom, set()).add(s.pos)
    return {c: sorted(p) for c, p in out.items()}


def expected_ploidy(chrom: str, pos: int, sex: str, build: int) -> int:
    if chrom in ("MT", "Y"):
        return 1
    if chrom == "X" and sex == "Male":
        par = PAR_X[build]
        return 2 if (par[0][0] <= pos <= par[0][1] or par[1][0] <= pos <= par[1][1]) else 1
    return 2


def genotype(d: SeqData, site: Site, whole_genome: bool) -> tuple[str | None, str]:
    """(genotype, basis). Genotype uses bases for SNVs and I/D for indels (I = the longer allele), like an array."""
    if site.chrom == "Y" and d.sex == "Female":
        return None, "no Y chromosome"
    k = key(site.chrom, site.pos)
    ploidy = expected_ploidy(site.chrom, site.pos, d.sex, d.build)
    here = d.calls.get(k, [])
    if site.indel:
        target = {trim(site.pos, site.ref, a) for a in site.alts}
        longer = len(site.alts[0]) > len(site.ref)
        ref_letter, alt_letter = ("D", "I") if longer else ("I", "D")
        ref_out = ref_letter
        match = [c for c in here if (c.pos, c.ref, c.alt) in target]
        other = [c for c in here if (c.pos, c.ref, c.alt) not in target and len(c.ref) != len(c.alt)]
    else:
        ref_out = site.ref
        match = [c for c in here if len(c.ref) == 1 and len(c.alt) == 1]
        other = []
    if not site.indel and any(c.ref != site.ref for c in match):
        return None, "reference base differs from the genome build (wrong build?)"
    span = d.spanning.get(k, [])
    if any(c.alt == "." for c in span):
        return None, "no genotype called here"
    if match:
        bad = next((c for c in match if not c.ok), None)
        if bad:
            return None, f"low-quality call ({bad.why})"
        if other or span:
            return None, "another indel overlaps this site"
        n_ploidy = match[0].ploidy
        letters = [(alt_letter if site.indel else c.alt) for c in match for _ in range(c.copies)]
        letters = [ref_out] * max(n_ploidy - len(letters), 0) + letters
        if ploidy == 1 and len(set(letters)) == 1:
            letters = letters[:1]            # a caller that writes male X or MT as diploid 1/1
        return "".join(letters), CALLED
    if other or span:
        return None, "another indel overlaps this site"
    if d.observed:
        cov = d.covered.get(k)
        if cov is None:
            return None, "not covered"
        return (ref_out * ploidy, READ) if cov[2] else (None, "covered too thinly to call" if d.gvcf or (cov[0] or 0) < 8
                                                        else "the reads show a variant here that the VCF does not list")
    if not whole_genome:
        return None, "not in the VCF (partial data, so unknown)"
    return ref_out * ploidy, ASSUMED


def view(d: SeqData, sites: dict, whole_genome: bool) -> pd.DataFrame:
    """The array-shaped table: index rsid; columns chrom, pos, genotype (NaN = no-call), basis."""
    rows = []
    for s in sites.values():
        g, basis = genotype(d, s, whole_genome)
        rows.append((s.rsid, s.chrom, s.pos, g if g is not None else float("nan"), basis))
    df = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "genotype", "basis"]).set_index("rsid")
    return df
