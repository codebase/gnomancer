"""Structural and copy-number variants from the vendor's own SV and CNV VCFs.

A small-variant VCF cannot show a deletion that removes an exon or a duplicated gene; whole-genome vendors
deliver those calls as separate VCFs (Sequencing.com ships an SV file and a CNV file). This module reads them
and answers one narrow question well: does any deletion or duplication overlap a gene where copy number is known
to matter?

The lens is ClinGen's dosage-sensitivity curation (data/dosage_genes.tsv): a deletion matters when the gene is
haploinsufficient (score 3), a duplication when it is triplosensitive (score 3); in a gene tied to a recessive
condition (score 30) a one-copy deletion is a carrier state. The ACMG actionable genes and a few genes where
copy number changes drug handling (CYP2D6, CYP2A6, GSTM1, UGT2B17) or carrier status (SMN1, HBA1/HBA2) are
listed whenever something overlaps them.

Limits worth stating wherever this is shown: short-read SV calling produces many false positives and misses
events in repetitive regions; a person carries several thousand real structural variants, nearly all benign;
overlap is computed against the gene's span, so a deletion wholly inside an intron is listed too; and callers
are unreliable exactly in the paralogous genes that matter most (CYP2D6, SMN1, HBA1/HBA2, CYP21A2), which need a
dedicated caller on the BAM. Anything listed is a lead for a clinical test (MLPA or a chromosomal microarray),
not a result.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .paths import DATA
from .seqvcf import CHR1_LENGTH, norm_chrom, open_text

DOSAGE = DATA / "dosage_genes.tsv"
MIN_SIZE = 50
LARGE = 1_000_000

# Genes not in ClinGen's list that are still worth naming: (gene, chrom, start37, end37, start38, end38).
EXTRA_GENES = [
    ("DES", "2", 220283099, 220291461, 219417474, 219426740),
    ("BTD", "3", 15642848, 15687329, 15601341, 15722522),
    ("GAA", "17", 78075355, 78093678, 80101526, 80119881),
    ("RPE65", "1", 68894505, 68915642, 68428822, 68451103),
    ("CYP2D6", "22", 42522501, 42526908, 42120962, 42131727),
    ("CYP2A6", "19", 41349443, 41356352, 40841915, 40850493),
    ("GSTM1", "1", 110230436, 110251661, 109687309, 109709552),
    ("UGT2B17", "4", 69402902, 69434245, 68516301, 68627039),
    ("SMN1", "5", 70220768, 70249769, 70925030, 70953942),
    ("SMN2", "5", 69345350, 69374349, 70049635, 70078522),
    ("HBB", "11", 5246694, 5250625, 5225464, 5229395),
]
# Why a copy-number change in these is worth a line even without a ClinGen score.
SPECIAL = {
    "CYP2D6": "deletion (*5) or duplication changes how many drugs are metabolized; short-read callers are unreliable here, so confirm with a CYP2D6-specific caller",
    "CYP2A6": "deletion (*4) slows nicotine and some drug metabolism",
    "GSTM1": "the whole-gene deletion is common (about half of people lack both copies) and is not a disease finding",
    "UGT2B17": "the whole-gene deletion is common and affects testosterone glucuronidation; not a disease finding",
    "SMN1": "loss of one copy is spinal muscular atrophy carrier status; SMN1 and SMN2 are near-identical, so this needs a dedicated caller",
    "SMN2": "copy number modifies spinal muscular atrophy severity; needs a dedicated caller",
    "HBA1": "deletion of alpha-globin genes is alpha-thalassemia trait; the two genes are near-identical, so confirm clinically",
    "HBA2": "deletion of alpha-globin genes is alpha-thalassemia trait; the two genes are near-identical, so confirm clinically",
    "HBB": "deletions here cause beta-thalassemia or hereditary persistence of fetal hemoglobin",
    "CYP21A2": "deletions and conversions with its pseudogene cause congenital adrenal hyperplasia; needs a dedicated caller",
}


@dataclass
class SV:
    chrom: str
    start: int
    end: int
    kind: str           # DEL, DUP, INS, INV, BND, CNV
    length: int
    passed: bool
    filt: str
    genotype: str       # 'heterozygous', 'homozygous', 'hemizygous', or '' when not given
    copy_number: int | None
    source: str


def load_genes(build: int) -> dict:
    """chrom -> [(start, end, name, hi, ts, kind)] on the given build: ClinGen genes and recurrent regions, plus the extras above."""
    out: dict = {}
    seen = set()
    with open(DOSAGE, newline="") as f:
        for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
            if not r[f"start{build}"]:
                continue
            seen.add(r["gene"])
            out.setdefault(r["chrom"], []).append((int(r[f"start{build}"]), int(r[f"end{build}"]), r["gene"],
                                                   int(r["hi"]) if r["hi"] else None, int(r["ts"]) if r["ts"] else None, r.get("kind") or "gene"))
    for gene, chrom, s37, e37, s38, e38 in EXTRA_GENES:
        if gene not in seen:
            s, e = (s37, e37) if build == 37 else (s38, e38)
            out.setdefault(chrom, []).append((s, e, gene, None, None, "gene"))
    return {c: sorted(v) for c, v in out.items()}


def _info(s: str) -> dict:
    d = {}
    for kv in s.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v
    return d


def _int(s):
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return None


def parse(path: str | Path, sex: str = "", label: str = "") -> tuple[list, dict]:
    """Read one SV or CNV VCF. Returns (events the sample carries, header facts)."""
    events, meta = [], dict(build=0, records=0, path=str(path))
    label = label or Path(path).name
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith("##contig=<") and not meta["build"]:
                    kv = dict(p.split("=", 1) for p in line[10:].rstrip(">\n").split(",") if "=" in p)
                    if norm_chrom(kv.get("ID", "")) == "1" and kv.get("length", "").isdigit():
                        meta["build"] = CHR1_LENGTH.get(int(kv["length"]), 0)
                continue
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            meta["records"] += 1
            chrom = norm_chrom(c[0])
            if chrom is None:
                continue
            pos, ref, alt, filt = int(c[1]), c[3], c[4], c[6]
            info = _info(c[7])
            kind = info.get("SVTYPE", "")
            if not kind:
                kind = ("BND" if "[" in alt or "]" in alt else alt.strip("<>").split(":")[0] if alt.startswith("<")
                        else "DEL" if len(ref) > len(alt) else "INS" if len(alt) > len(ref) else "")
            kind = kind.upper().split(":")[0]
            fmt = c[8].split(":") if len(c) > 9 else []
            smp = c[9].split(":") if len(c) > 9 else []
            field = lambda k: smp[fmt.index(k)] if k in fmt and fmt.index(k) < len(smp) else None  # noqa: E731
            gt = field("GT") or ""
            cn = _int(field("CN"))
            alleles = [a for a in gt.replace("|", "/").split("/") if a not in ("", ".")]
            expected = 1 if (chrom == "Y" or (chrom == "X" and sex == "Male")) else 2
            if kind in ("CNV", "") and cn is not None:
                if cn == expected:
                    continue                       # a segment at the normal copy number is not a variant
                kind = "DEL" if cn < expected else "DUP"
            if alleles and all(a == "0" for a in alleles) and (cn is None or cn == expected):
                continue                           # homozygous reference: the sample does not carry it
            end = _int(info.get("END"))
            svlen = _int(info.get("SVLEN", "").split(",")[0]) if info.get("SVLEN") else None
            if end is None:
                end = pos + abs(svlen) if svlen and kind != "INS" else pos + max(len(ref) - 1, 0)
            length = abs(svlen) if svlen else (max(len(alt) - len(ref), 0) if kind == "INS" and not alt.startswith("<") else end - pos)
            if kind not in ("BND", "INV") and length < MIN_SIZE:
                continue
            if cn is not None and kind in ("DEL", "DUP"):
                zyg = ("homozygous" if cn == 0 and expected == 2 else "hemizygous" if expected == 1 and kind == "DEL"
                       else "heterozygous" if kind == "DEL" else "")
            elif alleles:
                zyg = "hemizygous" if len(alleles) == 1 else "homozygous" if len(set(alleles)) == 1 else "heterozygous"
            else:
                zyg = ""
            events.append(SV(chrom, pos, end, kind, length, filt in ("PASS", "."), filt, zyg, cn, label))
    return events, meta


def _size(n: int) -> str:
    return f"{n / 1e6:.2f} Mb" if n >= 1e6 else f"{n / 1e3:.1f} kb" if n >= 1000 else f"{n} bp"


def analyze(paths: list, build: int, sex: str = "", outdir: Path | None = None) -> dict:
    """Parse the files, find gene overlaps, and return {'stats', 'findings' (DataFrame), 'section', 'notes'}."""
    from .clinvar import acmg_genes
    events, notes = [], []
    for p in paths:
        ev, meta = parse(p, sex)
        if meta["build"] and meta["build"] != build:
            notes.append(f"{Path(p).name} is on GRCh{meta['build']} but the main file is on GRCh{build}; it was skipped.")
            continue
        if not meta["records"]:
            notes.append(f"{Path(p).name} holds no records.")
        events += ev
    genes = load_genes(build)
    acmg = acmg_genes()
    passed = [e for e in events if e.passed]
    by_kind: dict = {}
    for e in passed:
        by_kind[e.kind] = by_kind.get(e.kind, 0) + 1
    rows = []
    for e in passed:
        hits = [g for g in genes.get(e.chrom, []) if g[0] <= e.end and g[1] >= e.start] if e.kind in ("DEL", "DUP") else \
               [g for g in genes.get(e.chrom, []) if g[0] <= e.start <= g[1] or (e.kind == "INV" and g[0] <= e.end <= g[1])]
        for gs, ge, gene, hi, ts, what in hits:
            whole = e.start <= gs and e.end >= ge
            share = (min(e.end, ge) - max(e.start, gs) + 1) / (ge - gs + 1)
            why, rank = "", 9
            copies = {"homozygous": "both copies", "hemizygous": "the only copy"}.get(e.genotype, "one copy")
            if what == "region":
                # A recurrent region is a syndrome only when most of it is gained or lost, not for a small event inside it.
                if share >= 0.5 and e.kind == "DEL" and hi == 3:
                    why, rank = f"loss of {copies} of a recurrent region whose deletion causes a known syndrome (ClinGen: sufficient evidence)", 0
                elif share >= 0.5 and e.kind == "DUP" and ts == 3:
                    why, rank = "a duplication of a recurrent region whose duplication causes a known syndrome (ClinGen: sufficient evidence)", 0
            elif e.kind == "DEL" and hi == 3:
                why, rank = f"loss of {copies} of a haploinsufficient gene (ClinGen: sufficient evidence)", 0
            elif e.kind == "DUP" and ts == 3 and whole:
                why, rank = "an extra copy of a triplosensitive gene (ClinGen: sufficient evidence)", 0
            elif e.kind == "DEL" and hi == 30:
                why, rank = ("both copies deleted in a gene tied to a recessive condition" if e.genotype == "homozygous"
                             else "one copy deleted in a gene tied to a recessive condition: a carrier state"), (1 if e.genotype == "homozygous" else 3)
            elif e.kind in ("INV", "BND") and hi == 3:
                why, rank = "a rearrangement breakpoint inside a haploinsufficient gene, which may disrupt it", 2
            elif gene in SPECIAL and e.kind in ("DEL", "DUP"):
                why, rank = SPECIAL[gene], 4
            elif gene in acmg and e.kind in ("DEL", "DUP"):
                why, rank = f"overlaps an ACMG actionable gene ({acmg[gene]['condition']}); the effect of a {'deletion' if e.kind == 'DEL' else 'duplication'} here is not established", 3
            if not why:
                continue
            if gene in SPECIAL and gene not in why and rank < 4:
                why += "; " + SPECIAL[gene]
            rows.append(dict(rank=rank, gene=gene, type=e.kind, zygosity=e.genotype or "not given", chrom=e.chrom, start=e.start, end=e.end,
                             size=_size(e.length), gene_share=(f"{share:.0%} of the region" if what == "region" else "whole gene" if whole else f"{share:.0%} of the gene span"),
                             copy_number=e.copy_number if e.copy_number is not None else "", meaning=why, acmg=gene in acmg, source=e.source))
    large = [e for e in passed if e.kind in ("DEL", "DUP") and e.length >= LARGE]
    for e in large:
        if not any(r["chrom"] == e.chrom and r["start"] == e.start and r["end"] == e.end for r in rows):
            rows.append(dict(rank=5, gene="(none curated)", type=e.kind, zygosity=e.genotype or "not given", chrom=e.chrom, start=e.start, end=e.end,
                             size=_size(e.length), gene_share="", copy_number=e.copy_number if e.copy_number is not None else "",
                             meaning="larger than 1 Mb: uncommon in healthy people, and a frequent caller artifact near centromeres and segmental duplications",
                             acmg=False, source=e.source))
    findings = pd.DataFrame(rows)
    if not findings.empty:
        findings = findings.sort_values(["rank", "gene", "start"]).drop_duplicates(["gene", "type", "start", "end"]).reset_index(drop=True)
    stats = dict(files=[Path(p).name for p in paths], events=len(events), passed=len(passed), failed=len(events) - len(passed), by_type=by_kind,
                 large=len(large), findings=len(findings), build=build)
    if outdir is not None and not findings.empty:
        findings.drop(columns=["rank"]).to_csv(Path(outdir) / "structural-variants.tsv", sep="\t", index=False)
    names = dict(DEL="deletion", DUP="duplication", INS="insertion", INV="inversion", BND="breakend")
    kinds = ", ".join(f"{n:,} {names.get(k, k)}{'s' if n != 1 and k in names else ''}" for k, n in sorted(by_kind.items(), key=lambda kv: -kv[1])) or "none"
    strong = int((findings["rank"] <= 2).sum()) if not findings.empty else 0
    text = (f"The vendor's files list {len(passed):,} structural or copy-number variants of 50 bp or more that pass the caller's filters ({kinds}); "
            f"{len(events) - len(passed):,} more failed the caller's filters and were ignored. A person carries several thousand such variants, nearly all harmless. "
            + (f"{len(findings)} overlap a gene where copy number is known or suspected to matter" + (f", {strong} of them in a way ClinGen rates as sufficient evidence for disease" if strong else "")
               + ". Each is a lead for a clinical test, not a result: short-read structural calls are often wrong, and the overlap is with the gene's span, not specifically its coding sequence."
               if len(findings) else "None overlaps a gene where copy number is known to matter (ClinGen dosage-sensitive genes, the ACMG actionable genes, and a few drug and carrier genes)."))
    section = dict(id="sv", title="Structural and copy-number variants", nav="Structural variants", sensitive=bool(len(findings)),
                   lead="Deletions, duplications and rearrangements of 50 bases or more, from the vendor's own structural-variant and copy-number files, "
                        "checked against the genes where gaining or losing a copy is known to matter.",
                   text=text, columns=["Gene", "Type", "Zygosity", "Location", "Size", "Overlap", "What it could mean"],
                   rows=[[r.gene, r.type, r.zygosity, f"chr{r.chrom}:{r.start:,}-{r.end:,}", r["size"], r.gene_share, r.meaning] for _, r in findings.iterrows()],
                   notes=notes + ["Genes with a near-identical copy (CYP2D6, SMN1, HBA1/HBA2, CYP21A2, PMS2) are where structural callers are least reliable; "
                                  "read any row there as a prompt to run a dedicated caller or a clinical test."])
    return dict(stats=stats, findings=findings, section=section, notes=notes)


def glance_bullets(res: dict) -> tuple[list, list]:
    """(notable, reassuring) lines for the at-a-glance summary."""
    f = res["findings"]
    if f.empty:
        return [], ["No structural or copy-number variant overlaps a gene where copy number is known to matter."]
    notable = []
    for _, r in f[f["rank"] <= 2].iterrows():
        kind = dict(DEL="deletion", DUP="duplication", INV="inversion", BND="rearrangement").get(r.type, r.type)
        notable.append(f"Structural variant: {kind} of {r['size']} at {r.gene}: {r.meaning}. Short-read structural calls need clinical confirmation.")
    rest = f[f["rank"] > 2]
    if len(rest):
        genes = [g for g in rest.gene.drop_duplicates() if not g.startswith("(")][:8]
        notable.append(f"{len(rest)} other structural variant(s) are listed" + (f" ({', '.join(genes)})" if genes else "") + ": carrier-type, drug-gene or unconfirmed; see the table.")
    return notable, []
