"""Specialist callers for the parts of a genome a small-variant VCF cannot describe.

Each function wraps one established open-source tool, run on the BAM or CRAM, and returns a dict with the
parsed result and a report section. Each is optional: if its tool is not installed (`gnomancer setup --wgs`
installs them where a build exists for the platform) the function says so and the report lists the topic under
its limits instead. None of these is a clinical test, and each section says what would confirm it.

- repeats:  ExpansionHunter sizes known disease-causing repeats; stranger (Clinical Genomics Stockholm) supplies
            the catalog and the normal / pathogenic thresholds.
- cyp2d6:   PyPGx calls CYP2D6 star alleles including whole-gene deletion, duplication and CYP2D7 hybrids, from
            read depth plus variants. PharmCAT, which will not call CYP2D6 from a VCF, turns the result into a
            phenotype.
- hla:      T1K types the classical HLA genes against the IPD-IMGT/HLA database.
- smn:      Illumina's SMNCopyNumberCaller counts SMN1 and SMN2 copies (spinal muscular atrophy carrier status).
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

from . import wgstools
from .bam import Bam
from .paths import HOME

THREADS = "4"

# Names for the abbreviations in stranger's catalog.
REPEAT_DISEASE = {
    "HTT": "Huntington disease", "FMR1": "Fragile X syndrome; premutation: tremor/ataxia syndrome (FXTAS) and primary ovarian insufficiency",
    "DMPK": "Myotonic dystrophy type 1", "CNBP": "Myotonic dystrophy type 2", "C9ORF72": "ALS / frontotemporal dementia",
    "ATXN1": "Spinocerebellar ataxia 1", "ATXN2": "Spinocerebellar ataxia 2", "ATXN3": "Spinocerebellar ataxia 3 (Machado-Joseph)",
    "CACNA1A": "Spinocerebellar ataxia 6", "ATXN7": "Spinocerebellar ataxia 7", "ATXN8OS": "Spinocerebellar ataxia 8",
    "ATXN10": "Spinocerebellar ataxia 10", "PPP2R2B": "Spinocerebellar ataxia 12", "TBP": "Spinocerebellar ataxia 17",
    "FGF14": "Spinocerebellar ataxia 27B", "BEAN1": "Spinocerebellar ataxia 31", "NOP56": "Spinocerebellar ataxia 36",
    "DAB1": "Spinocerebellar ataxia 37", "ZFHX3": "Spinocerebellar ataxia 4", "THAP11": "Spinocerebellar ataxia 51",
    "ATN1": "Dentatorubral-pallidoluysian atrophy", "AR": "Spinal and bulbar muscular atrophy (Kennedy disease)",
    "FXN": "Friedreich ataxia", "RFC1": "CANVAS (cerebellar ataxia, neuropathy, vestibular areflexia)", "CSTB": "Progressive myoclonic epilepsy (Unverricht-Lundborg)",
    "PABPN1": "Oculopharyngeal muscular dystrophy", "PHOX2B": "Congenital central hypoventilation syndrome", "TCF4": "Fuchs endothelial corneal dystrophy",
    "NOTCH2NLC": "Neuronal intranuclear inclusion disease", "JPH3": "Huntington disease-like 2", "AFF2": "Fragile XE intellectual disability",
    "GLS": "Glutaminase deficiency", "HOXD13": "Synpolydactyly", "FOXL2": "Blepharophimosis syndrome (BPES)", "RUNX2": "Cleidocranial dysplasia",
    "ZIC2": "Holoprosencephaly 5", "ARX_EIEE": "ARX-related developmental and epileptic encephalopathy", "SOX3": "X-linked hypopituitarism",
    "LRP12": "Oculopharyngodistal myopathy 1", "GIPC1": "Oculopharyngodistal myopathy 2", "RILPL1": "Oculopharyngodistal myopathy 4",
    "ABCD3": "Oculopharyngodistal myopathy", "SAMD12": "Familial adult myoclonic epilepsy 1", "RAPGEF2": "Familial adult myoclonic epilepsy 7",
    "TNRC6A": "Familial adult myoclonic epilepsy 6", "XYLT1": "Desbuquois dysplasia 2", "DIP2B": "FRA12A intellectual disability",
    "CBL": "FRA11B / Jacobsen syndrome region", "EIF4A3": "Richieri-Costa-Pereira syndrome", "TBX1": "TBX1-related conotruncal anomalies",
    "ZIC3": "X-linked heterotaxy (VACTERL-X)",
}

# Loci where repeat length alone does not decide anything, because only certain motifs are harmful and this
# method sizes the repeat without telling the motifs apart.
MOTIF_DEPENDENT = {
    "RFC1": "an expanded allele is present, but only certain motifs (chiefly AAGGG) cause disease, the harmless AAAAG expansion is common "
            "and looks the same here, and the condition needs the pathogenic motif on both copies; not interpretable without a motif-specific test",
}

# HLA alleles with an established drug or disease link: (gene, allele prefix) -> what it means.
HLA_NOTES = [
    ("HLA-B", "B*57:01", "abacavir hypersensitivity: the drug is contraindicated"),
    ("HLA-B", "B*15:02", "carbamazepine, oxcarbazepine and phenytoin: high risk of Stevens-Johnson syndrome"),
    ("HLA-B", "B*58:01", "allopurinol: high risk of severe skin reactions"),
    ("HLA-A", "A*31:01", "carbamazepine: raised risk of hypersensitivity reactions"),
    ("HLA-B", "B*27", "associated with ankylosing spondylitis and related arthritis (most carriers never develop it)"),
    ("HLA-B", "B*51", "associated with Behcet disease"),
    ("HLA-C", "C*06:02", "the main psoriasis risk allele"),
    ("HLA-A", "A*29:02", "associated with birdshot chorioretinopathy"),
    ("HLA-DQB1", "DQB1*06:02", "strongly associated with narcolepsy type 1 (common in the population; most carriers are unaffected)"),
    ("HLA-DRB1", "DRB1*15:01", "the main multiple sclerosis risk allele"),
    ("HLA-DRB1", "DRB1*04:01", "the main rheumatoid arthritis risk allele; also raises type 1 diabetes risk"),
    ("HLA-DRB1", "DRB1*03:01", "DR3: type 1 diabetes and other autoimmune risk"),
]
HLA_GENES = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")


def _run(cmd: list, group: str, cwd: Path | None = None, extra_env: dict | None = None, timeout: int = 6 * 3600) -> subprocess.CompletedProcess:
    env = wgstools.env(group)
    env.update(extra_env or {})
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, env=env, cwd=str(cwd) if cwd else None, timeout=timeout)


def _err(p: subprocess.CompletedProcess) -> str:
    """The line of a failed tool's output most likely to say why, and its exit status. Tools log progress to the same
    stream, so the last line is often not the error."""
    lines = [l.strip() for l in ((p.stderr or "") + "\n" + (p.stdout or "")).splitlines() if l.strip()]
    telling = [l for l in lines if re.search(r"error|exception|fail|abort|cannot|could not|not found|no such|segmentation|killed", l, re.I)]
    last = telling[-1] if telling else (lines[-1] if lines else "no output")
    status = f"killed by signal {-p.returncode}" if p.returncode < 0 else f"exit status {p.returncode}"
    return f"{last[:240]} ({status})"


def _unavailable(what: str, why: str) -> dict:
    return dict(available=False, error=why, section=None, what=what)


# ------------------------------------------------------------------------------------------- repeat expansions

def repeats(b: Bam, outdir: Path, sex: str) -> dict:
    eh, stranger = wgstools.tool_bin("ExpansionHunter", "repeats"), wgstools.tool_bin("stranger", "repeats")
    if not eh or not stranger:
        return _unavailable("repeats", "ExpansionHunter and stranger are not installed (`gnomancer setup --wgs`)")
    if not b.reference:
        return _unavailable("repeats", "needs the reference genome")
    chr_named = b.name("1").startswith("chr")
    tag = {(38, True): "hg38", (38, False): "grch38", (37, True): "hg19", (37, False): "grch37"}[(b.build, chr_named)]
    res_dir = next(iter(Path(stranger).resolve().parent.parent.glob("lib/python*/site-packages/stranger/resources")), None)
    catalog = res_dir / f"variant_catalog_{tag}.json" if res_dir else None
    if not catalog or not catalog.exists():
        return _unavailable("repeats", "stranger's repeat catalog was not found")
    wd = (outdir / "repeats").resolve()
    wd.mkdir(parents=True, exist_ok=True)
    prefix = wd / "eh"
    reads = b.path
    if b.cram:
        # ExpansionHunter 5.0.0 crashes reading CRAM on some builds. It only looks at reads near each locus, its
        # off-target regions and their mates, so those are pulled into a small BAM (whole pairs, wherever the mate maps).
        reads = wd / "repeat_regions.bam"
        bed = wd / "regions.bed"
        rows = set()
        for x in json.loads(catalog.read_text()):
            regs = x["ReferenceRegion"] if isinstance(x["ReferenceRegion"], list) else [x["ReferenceRegion"]]
            for r in regs + list(x.get("OfftargetRegions") or []):
                chrom, span = r.rsplit(":", 1)
                start, end = span.split("-")
                rows.add((chrom, max(int(start) - 2000, 0), int(end) + 2000))
        bed.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in sorted(rows)))
        samtools = wgstools.tool_bin("samtools")
        q = _run([samtools, "view", "-@", THREADS, "-b", "--fetch-pairs", "-L", bed, "-T", b.reference, "-o", wd / "unsorted.bam", b.path], "repeats")
        if q.returncode == 0:
            q = _run([samtools, "sort", "-@", THREADS, "-o", reads, wd / "unsorted.bam"], "repeats")
        (wd / "unsorted.bam").unlink(missing_ok=True)
        if q.returncode != 0:
            return _unavailable("repeats", "could not extract the repeat regions from the CRAM: " + _err(q))
        _run([samtools, "index", reads], "repeats")
    p = _run([eh, "--reads", reads, "--reference", b.reference, "--variant-catalog", catalog, "--output-prefix", prefix,
              "--sex", "male" if sex == "Male" else "female", "--threads", THREADS], "repeats")
    if b.cram:
        reads.unlink(missing_ok=True)
        Path(str(reads) + ".bai").unlink(missing_ok=True)
    if p.returncode != 0 or not Path(str(prefix) + ".vcf").exists():
        return _unavailable("repeats", "ExpansionHunter failed: " + _err(p))
    q = _run([stranger, "-f", catalog, str(prefix) + ".vcf"], "repeats")
    if q.returncode != 0:
        return _unavailable("repeats", "stranger failed: " + _err(q))
    (wd / "repeats.annotated.vcf").write_text(q.stdout)
    rows = parse_repeats(q.stdout, sex)
    return dict(available=True, error=None, loci=rows, section=repeats_section(rows), what="repeats")


def parse_repeats(vcf_text: str, sex: str = "") -> list:
    """One row per repeat locus from a stranger-annotated ExpansionHunter VCF."""
    by_locus: dict = {}
    for line in vcf_text.splitlines():
        if not line or line[0] == "#":
            continue
        c = line.split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        fmt, smp = c[8].split(":"), c[9].split(":")
        get = lambda k: smp[fmt.index(k)] if k in fmt and fmt.index(k) < len(smp) else ""  # noqa: E731
        locus = info.get("REPID") or info.get("VARID") or c[2]
        if "STR_STATUS" not in info:
            continue        # a flanking repeat unit of a multi-part locus (HTT's CCG tract and the like): sized, but not what the thresholds refer to
        counts = [int(x) for x in get("REPCN").split("/") if x.isdigit()]
        ci = get("REPCI")
        row = dict(locus=locus, variant=info.get("VARID", locus), chrom=c[0], pos=int(c[1]), unit=info.get("DisplayRU") or info.get("RU", ""),
                   counts=counts, interval=ci, support=get("SO"), status=info.get("STR_STATUS", ""),
                   normal_max=_num(info.get("STR_NORMAL_MAX")), pathologic_min=_num(info.get("STR_PATHOLOGIC_MIN")),
                   inheritance=info.get("InheritanceMode", ""), disease=REPEAT_DISEASE.get(locus, info.get("Disease", "")), filter=c[6])
        prev = by_locus.get(locus)
        # A locus can have several repeat units; keep the record stranger scored (it carries the thresholds), worst status first.
        rank = {"full_mutation": 0, "pre_mutation": 1, "normal": 2, "": 3}
        if prev is None or (rank.get(row["status"], 3), -max(row["counts"] or [0])) < (rank.get(prev["status"], 3), -max(prev["counts"] or [0])):
            by_locus[locus] = row
    out = []
    for r in by_locus.values():
        longest = max(r["counts"]) if r["counts"] else None
        r["longest"] = longest
        r["estimated"] = "INREPEAT" in r["support"] or "FLANKING" in r["support"]
        expanded = [n for n in r["counts"] if r["pathologic_min"] and n >= r["pathologic_min"]]
        if r["status"] != "normal" and r["locus"] in MOTIF_DEPENDENT:
            r["reading"], r["level"] = MOTIF_DEPENDENT[r["locus"]], "uncertain"
        elif r["status"] == "full_mutation":
            if r["inheritance"] == "AR" and len(expanded) < 2 and len(r["counts"]) == 2:
                r["reading"] = "one allele in the disease range; the condition is recessive, so this is a carrier state"
                r["level"] = "carrier"
            else:
                r["reading"] = "in the disease range" + (" on both alleles" if len(expanded) == 2 else "")
                r["level"] = "expanded"
        elif r["status"] == "pre_mutation":
            r["reading"] = "above the normal range but below the disease threshold (intermediate or premutation)"
            r["level"] = "intermediate"
        elif r["status"] == "normal":
            r["reading"], r["level"] = "normal", "normal"
        else:
            r["reading"], r["level"] = "not scored", "unknown"
        out.append(r)
    order = {"expanded": 0, "carrier": 1, "intermediate": 2, "uncertain": 3, "unknown": 4, "normal": 5}
    return sorted(out, key=lambda r: (order[r["level"]], r["locus"]))


def _num(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def repeats_section(rows: list) -> dict:
    n = len(rows)
    flagged = [r for r in rows if r["level"] in ("expanded", "carrier", "intermediate")]
    unsure = ", ".join(r["locus"] for r in rows if r["level"] == "uncertain")
    text = (f"{n} repeat loci tied to known conditions were sized from the reads. "
            + (f"{len(flagged)} {'is' if len(flagged) == 1 else 'are'} outside the normal range and listed first. " if flagged else
               "None is in a disease or intermediate range" + (f"; {unsure} could not be interpreted and is explained in its row. " if unsure else ". "))
            + "Short reads size a repeat exactly only while it fits inside one read (about 150 bases); longer alleles are estimated from partial reads, "
              "and the very large expansions of some conditions can only be flagged as 'long', not measured. A result outside the normal range needs a "
              "clinical repeat test (repeat-primed PCR or Southern blot) before it means anything.")
    table = []
    for r in rows:
        counts = "/".join(map(str, r["counts"])) or "no call"
        rng = (f"normal up to {r['normal_max']}" if r["normal_max"] is not None else "") + (f", disease from {r['pathologic_min']}" if r["pathologic_min"] is not None else "")
        table.append([r["locus"], r["disease"], r["unit"], counts + (" (estimated)" if r["estimated"] and r["level"] != "normal" else ""), rng, r["inheritance"], r["reading"]])
    return dict(id="repeats", title="Repeat expansions", nav="Repeats", sensitive=True,
                lead="Some neurological and developmental conditions are caused by a short DNA motif repeated too many times. These are invisible "
                     "in a variant VCF and are sized here directly from the reads with ExpansionHunter, against published normal and disease ranges.",
                text=text, columns=["Gene", "Condition", "Motif", "Your repeat counts", "Reference ranges", "Inheritance", "Reading"], rows=table,
                notes=["Ranges are from the stranger catalog (Clinical Genomics Stockholm), which follows the published literature per locus. "
                       "Inheritance: AD dominant, AR recessive (both alleles must be expanded), XR/XD X-linked."])


def repeats_bullets(res: dict) -> tuple[list, list]:
    if not res.get("available"):
        return [], []
    rows = res["loci"]
    notable = [f"Repeat expansion: {r['locus']} ({r['disease']}) has {'/'.join(map(str, r['counts']))} repeats, {r['reading']}. Needs a clinical repeat test to confirm."
               for r in rows if r["level"] in ("expanded", "carrier")]
    mid = [r for r in rows if r["level"] == "intermediate"]
    if mid:
        notable.append("Repeat lengths in the intermediate range at: " + ", ".join(f"{r['locus']} ({max(r['counts'])})" for r in mid) + ". Usually without symptoms; some can expand in the next generation.")
    clear = sum(1 for r in rows if r["level"] == "normal")
    good = [] if notable else [f"{clear} of {len(rows)} disease-associated repeat loci checked are in the normal range (Huntington, fragile X, myotonic dystrophy and the ataxias among them)"
                               + ("; the rest could not be interpreted from short reads." if clear < len(rows) else ".")]
    return notable, good


# ----------------------------------------------------------------------------------------------------- CYP2D6

def pypgx_bundle() -> Path | None:
    for p in (HOME / "tools" / "pypgx-bundle", Path.home() / "pypgx-bundle"):
        if (p / "cnv").exists():
            return p
    return None


def cyp2d6(b: Bam, outdir: Path, genes=("CYP2D6",)) -> dict:
    pypgx = wgstools.tool_bin("pypgx", "pgx")
    bundle = pypgx_bundle()
    if not pypgx or not bundle:
        return _unavailable("cyp2d6", "PyPGx and its data bundle are not installed (`gnomancer setup --wgs`)")
    if not b.reference:
        return _unavailable("cyp2d6", "needs the reference genome")
    wd = (outdir / "pypgx").resolve()
    wd.mkdir(parents=True, exist_ok=True)
    env = {"PYPGX_BUNDLE": str(bundle)}
    asm = f"GRCh{b.build}"
    reads = b.path
    if b.cram:
        # PyPGx reads BAM only. It needs just the pharmacogene and control-gene regions, so slice those out of the CRAM.
        bed = wd / "regions.bed"
        q = _run([pypgx, "create-regions-bed", "--assembly", asm, "--merge"] + (["--add-chr-prefix"] if b.name("1").startswith("chr") else []), "pgx", extra_env=env)
        if q.returncode != 0 or not q.stdout.strip():
            return _unavailable("cyp2d6", "PyPGx could not list its regions: " + _err(q))
        bed.write_text(q.stdout)
        reads = wd / "pgx_regions.bam"
        samtools = wgstools.tool_bin("samtools")
        q = _run([samtools, "view", "-@", THREADS, "-b", "-M", "-L", bed, "-T", b.reference, "-o", reads, b.path], "pgx")
        if q.returncode != 0:
            return _unavailable("cyp2d6", "could not slice the CRAM for PyPGx: " + _err(q))
        _run([samtools, "index", reads], "pgx")
    steps = [
        ("variants", [pypgx, "create-input-vcf", wd / "variants.vcf.gz", b.reference, reads, "--assembly", asm, "--genes", *genes]),
        ("depth", [pypgx, "prepare-depth-of-coverage", wd / "depth.zip", reads, "--assembly", asm, "--genes", *genes]),
        ("control", [pypgx, "compute-control-statistics", "VDR", wd / "control.zip", reads, "--assembly", asm]),
    ]
    for name, cmd in steps:
        p = _run(cmd, "pgx", extra_env=env)
        if p.returncode != 0:
            return _unavailable("cyp2d6", f"PyPGx ({name}) failed: " + _err(p))
    calls = []
    for gene in genes:
        gdir = wd / gene
        if gdir.exists():
            shutil.rmtree(gdir)
        p = _run([pypgx, "run-ngs-pipeline", gene, gdir, "--variants", wd / "variants.vcf.gz", "--depth-of-coverage", wd / "depth.zip",
                  "--control-statistics", wd / "control.zip", "--assembly", asm, "--do-not-plot-copy-number", "--do-not-plot-allele-fraction"], "pgx", extra_env=env)
        if p.returncode != 0 or not (gdir / "results.zip").exists():
            return _unavailable("cyp2d6", f"PyPGx ({gene}) failed: " + _err(p))
        q = _run([pypgx, "print-data", gdir / "results.zip"], "pgx", extra_env=env)
        row = parse_pypgx(q.stdout)
        if row:
            calls.append(dict(gene=gene, **row))
    if b.cram:
        reads.unlink(missing_ok=True)
        Path(str(reads) + ".bai").unlink(missing_ok=True)
    if not calls:
        return _unavailable("cyp2d6", "PyPGx returned no genotype")
    return dict(available=True, error=None, calls=calls, section=cyp2d6_section(calls), what="cyp2d6")


def parse_pypgx(text: str) -> dict | None:
    """First sample row of `pypgx print-data results.zip` (tab-separated, header line first)."""
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    head, vals = lines[0].split("\t"), lines[1].split("\t")
    if len(vals) == len(head) + 1:          # the sample name column has no header
        head = ["Sample"] + head
    d = dict(zip(head, vals))
    return dict(genotype=d.get("Genotype", ""), phenotype=d.get("Phenotype", ""), cnv=d.get("CNV", ""),
                haplotype1=d.get("Haplotype1", ""), haplotype2=d.get("Haplotype2", ""), alternative=d.get("AlternativePhase", ""))


CNV_TEXT = {"Normal": "two copies, no structural change", "WholeDel1": "one copy deleted (*5)", "WholeDel2": "both copies deleted (*5/*5)",
            "WholeDup1": "one extra copy (duplication)", "WholeMultip1": "several extra copies (multiplication)", "Unknown": "copy number could not be classified",
            "Indeterminate": "copy number could not be determined"}


def cyp2d6_section(calls: list) -> dict:
    rows = []
    for c in calls:
        cnv = CNV_TEXT.get(c["cnv"], c["cnv"].replace(",", ", ") or "not assessed")
        rows.append([c["gene"], c["genotype"] or "not determined", c["phenotype"] or "not assigned", cnv])
    first = calls[0]
    text = (f"CYP2D6: {first['genotype'] or 'not determined'}" + (f", {first['phenotype'].lower()}" if first["phenotype"] else "") + ". "
            "CYP2D6 handles about a quarter of prescribed drugs (codeine, tramadol, many antidepressants and antipsychotics, tamoxifen, metoprolol, ondansetron). "
            "It sits beside a near-identical pseudogene and is deleted, duplicated or fused with it in many people, which is why it is called here from read depth "
            "and not from the VCF. Short-read calls are right for most people but can miss rare hybrid arrangements; confirm with a clinical pharmacogenomic test before a prescription depends on it.")
    return dict(id="cyp2d6", title="CYP2D6 from read depth", nav="CYP2D6", sensitive=False,
                lead="Star alleles with gene deletions, duplications and CYP2D7 hybrids, called by PyPGx from read depth and variants in the BAM.",
                text=text, columns=["Gene", "Diplotype", "Phenotype", "Structure"], rows=rows,
                notes=["Method: PyPGx (Lee et al., 2022) with VDR as the copy-number control; phasing against the 1000 Genomes panel."])


def cyp2d6_bullets(res: dict) -> tuple[list, list]:
    if not res.get("available"):
        return [], []
    c = res["calls"][0]
    line = f"CYP2D6 {c['genotype']}: {c['phenotype'] or 'phenotype not assigned'} (from read depth; affects codeine, tramadol, many antidepressants, tamoxifen)."
    return ([line], []) if c["phenotype"] and "normal" not in c["phenotype"].lower() else ([], [line])


# -------------------------------------------------------------------------------------------------------- HLA

def t1k_index() -> Path | None:
    d = HOME / "tools" / "t1k" / "hlaidx"
    return d if (d / "hlaidx_dna_seq.fa").exists() and (d / "hlaidx_dna_coord.fa").exists() else None


def hla(b: Bam, outdir: Path) -> dict:
    t1k, idx = wgstools.tool_bin("run-t1k", "hla"), t1k_index()
    if not t1k or not idx:
        return _unavailable("hla", "T1K and its HLA index are not installed (`gnomancer setup --wgs`)")
    if b.build != 38:
        return _unavailable("hla", "the HLA index is built for GRCh38")
    wd = (outdir / "hla").resolve()
    wd.mkdir(parents=True, exist_ok=True)
    # T1K's own extraction reads the whole file once (15 minutes or more for 30x). Pulling the extended MHC region and
    # the unmapped reads through the index first gives it the same reads in seconds.
    if b.cram and not b.reference:
        return _unavailable("hla", "a CRAM needs the reference genome")
    reads = wd / "hla_region.bam"
    samtools = wgstools.tool_bin("samtools")
    p = _run([samtools, "view", "-@", THREADS, "-b", *b.ref_args("-T"), "-o", reads, b.path, f"{b.name('6')}:28000000-34000000", "*"], "hla")
    if p.returncode != 0:
        return _unavailable("hla", "could not extract the HLA region: " + _err(p))
    _run([samtools, "index", reads], "hla")
    p = _run([t1k, "-b", reads, "-f", idx / "hlaidx_dna_seq.fa", "-c", idx / "hlaidx_dna_coord.fa", "--preset", "hla-wgs", "-t", THREADS,
              "--od", wd, "-o", "sample"], "hla")
    out = wd / "sample_genotype.tsv"
    if p.returncode != 0 or not out.exists():
        return _unavailable("hla", "T1K failed: " + _err(p))
    types = parse_t1k(out.read_text())
    reads.unlink(missing_ok=True)
    Path(str(reads) + ".bai").unlink(missing_ok=True)
    for f in wd.glob("sample_*.fq"):
        f.unlink()
    if not types:
        return _unavailable("hla", "T1K typed none of the classical HLA genes (too few reads over the region?)")
    return dict(available=True, error=None, types=types, section=hla_section(types), what="hla")


def _two_field(allele: str) -> str:
    """'HLA-A*02:01:01:01' -> 'A*02:01'."""
    a = allele.replace("HLA-", "")
    parts = a.split(":")
    return ":".join(parts[:2])


def parse_t1k(text: str) -> dict:
    """gene -> [(two-field allele, quality)], from T1K's genotype table (gene, n, allele1, abundance, quality, allele2, ...)."""
    out = {}
    for line in text.splitlines():
        c = line.split("\t")
        if len(c) < 8 or c[0] not in HLA_GENES:
            continue
        alleles = []
        for name, qual in ((c[2], c[4]), (c[5], c[7])):
            if name and name != ".":
                alleles.append((_two_field(name.split(",")[0]), _num(qual) or 0))
        if len(alleles) == 1:
            alleles.append(alleles[0])              # T1K lists a homozygous gene once
        out[c[0]] = alleles
    return out


def hla_notes(types: dict) -> list:
    """Clinically linked alleles the person carries: [(allele, copies, meaning)]."""
    found = []
    for gene, prefix, meaning in HLA_NOTES:
        mine = [a for a, _ in types.get(gene, [])]
        n = sum(1 for a in mine if a == prefix or a.startswith(prefix + ":"))
        if n:
            found.append((prefix, n, meaning))
    return found


def celiac(types: dict) -> str:
    """The HLA-DQ types celiac disease depends on. Nearly every patient carries DQ2.5, DQ8 or DQ2.2, so their absence
    makes the disease very unlikely; their presence is common (about a third of Europeans) and predicts little."""
    dqa = [a for a, _ in types.get("HLA-DQA1", [])]
    dqb = [a for a, _ in types.get("HLA-DQB1", [])]
    if not dqa or not dqb:
        return ""
    found = []
    if any(a.startswith("DQA1*05") for a in dqa) and any(b.startswith("DQB1*02") for b in dqb):
        found.append("DQ2.5")
    if any(b.startswith("DQB1*03:02") for b in dqb):
        found.append("DQ8")
    if any(a.startswith("DQA1*02:01") for a in dqa) and any(b.startswith("DQB1*02:02") for b in dqb):
        found.append("DQ2.2")
    if not found:
        return ("Celiac disease: none of the permissive types (DQ2.5, DQ8, DQ2.2) is present, which makes celiac disease very unlikely; "
                "about 99% of people with it carry one of them.")
    return (f"Celiac disease: {' and '.join(found)} present. These types are required for celiac disease but are carried by roughly a third of people, "
            "nearly all of whom never develop it; the result matters mainly for ruling it out, which it does not do here.")


def hla_section(types: dict) -> dict:
    rows = []
    for g in HLA_GENES:
        if g in types:
            (a1, q1), (a2, q2) = types[g][0], types[g][1]
            rows.append([g, a1, a2, "homozygous" if a1 == a2 else "", "ok" if min(q1, q2) > 0 else "low confidence"])
    notes = hla_notes(types)
    listed = "; ".join(f"{a} ({'two copies' if n > 1 else 'one copy'}): {m}" for a, n, m in notes)
    text = ("Classical HLA genes typed from the reads at two-field resolution. " +
            (f"Alleles with an established drug or disease link: {listed}. " if notes else
             "None of the alleles with a well-established drug-reaction link (B*57:01, B*15:02, B*58:01, A*31:01) is present. ") +
            (celiac(types) + " " if celiac(types) else "") +
            "HLA disease associations are risk factors carried by many healthy people, not diagnoses; the drug-reaction alleles are the ones that change prescribing. "
            "Typing from 30x short reads is accurate for most people but not at transplant-matching standard.")
    return dict(id="hla", title="HLA types", nav="HLA", sensitive=False,
                lead="The HLA genes are the most variable in the genome and cannot be read from a variant VCF. T1K assembles each person's alleles from the reads "
                     "against the IPD-IMGT/HLA database.",
                text=text, columns=["Gene", "Allele 1", "Allele 2", "", "Call"], rows=rows,
                notes=["The single-site HLA rows in the marker panel are tag SNPs standing in for these alleles; where they disagree, this table is the better evidence."])


def hla_bullets(res: dict) -> tuple[list, list]:
    if not res.get("available"):
        return [], []
    notes = hla_notes(res["types"])
    drug = [x for x in notes if any(k in x[2] for k in ("abacavir", "carbamazepine", "allopurinol"))]
    notable = [f"HLA-{a}: {m}. Tell any prescriber." for a, _, m in drug]
    other = [a for a, _, m in notes if (a, _, m) not in drug]
    if other:
        notable.append("HLA alleles with known disease associations: " + ", ".join(other) + ". Risk factors that many healthy people carry; see the HLA section.")
    good = [] if drug else ["No HLA allele linked to severe drug reactions (B*57:01, B*15:02, B*58:01, A*31:01) was found."]
    return notable, good


# -------------------------------------------------------------------------------------------------------- SMN

def smn(b: Bam, outdir: Path) -> dict:
    caller = wgstools.tool_bin("smn_caller.py", "smn")
    if not caller:
        return _unavailable("smn", "SMNCopyNumberCaller is not installed (`gnomancer setup --wgs`)")
    if b.cram and not b.reference:
        return _unavailable("smn", "a CRAM needs the reference genome")
    wd = (outdir / "smn").resolve()
    wd.mkdir(parents=True, exist_ok=True)
    manifest = wd / "manifest.txt"
    manifest.write_text(str(b.path) + "\n")
    genome = "38" if b.build == 38 else "19" if b.name("1").startswith("chr") else "37"
    cmd = [caller, "--manifest", manifest, "--genome", genome, "--prefix", "smn", "--outDir", wd, "--threads", THREADS]
    if b.cram:
        cmd += ["--reference", b.reference]
    p = _run(cmd, "smn")
    out = wd / "smn.tsv"
    if p.returncode != 0 or not out.exists():
        return _unavailable("smn", "SMNCopyNumberCaller failed: " + _err(p))
    r = parse_smn(out.read_text())
    if not r:
        return _unavailable("smn", "SMNCopyNumberCaller returned no call")
    return dict(available=True, error=None, call=r, section=smn_section(r), what="smn")


def parse_smn(text: str) -> dict | None:
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    d = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    flag = lambda v: {"True": True, "False": False}.get(v)  # noqa: E731
    return dict(smn1=_num(d.get("SMN1_CN")), smn2=_num(d.get("SMN2_CN")), is_sma=flag(d.get("isSMA", "")), is_carrier=flag(d.get("isCarrier", "")),
                silent_carrier_marker=_num(d.get("g.27134T>G_CN")), smn2_delta78=_num(d.get("SMN2delta7-8_CN")))


def smn_section(r: dict) -> dict:
    s1 = r["smn1"]
    if s1 is None:
        reading = "SMN1 copy number could not be called with confidence."
    elif s1 == 0:
        reading = "No functional SMN1 copy was found: the genotype that causes spinal muscular atrophy. Confirm clinically at once if this is unexpected."
    elif s1 == 1:
        reading = "One SMN1 copy: a carrier of spinal muscular atrophy. Carriers are healthy; it matters when a partner is also a carrier."
    else:
        reading = (f"{s1} SMN1 copies: not a carrier by copy number. " +
                   ("The g.27134T>G marker is present, which in some populations tags two SMN1 copies on one chromosome and none on the other (a 'silent' carrier)."
                    if (r.get("silent_carrier_marker") or 0) > 0 and s1 == 2 else
                    "About 3-4% of carriers have two copies on one chromosome and none on the other, which copy counting cannot see."))
    return dict(id="smn", title="SMN1 and SMN2 copy number", nav="SMN1", sensitive=True,
                lead="Spinal muscular atrophy is caused by losing both copies of SMN1, which has a near-identical neighbour, SMN2. Variant callers cannot tell "
                     "them apart; this counts copies of each from read depth at the bases that differ.",
                text=reading, columns=["SMN1 copies", "SMN2 copies", "Carrier call", "g.27134T>G copies"],
                rows=[[s1 if s1 is not None else "no call", r["smn2"] if r["smn2"] is not None else "no call",
                       {True: "carrier", False: "not a carrier", None: "no call"}[r["is_carrier"]], r.get("silent_carrier_marker") if r.get("silent_carrier_marker") is not None else ""]],
                notes=["Method: SMNCopyNumberCaller (Chen et al., Genetics in Medicine 2020). Carrier screening before pregnancy should use a clinical test."])


def smn_bullets(res: dict) -> tuple[list, list]:
    if not res.get("available"):
        return [], []
    s1 = res["call"]["smn1"]
    if s1 == 1:
        return ["SMN1: one copy, a carrier of spinal muscular atrophy (healthy; relevant if a partner is also a carrier). Confirm with a clinical carrier test."], []
    if s1 == 0:
        return ["SMN1: no functional copy found, the genotype of spinal muscular atrophy. Confirm clinically."], []
    if s1 is None:
        return [], []
    return [], [f"SMN1: {s1} copies, not a spinal muscular atrophy carrier by copy number."]


ALL = (("repeats", repeats, repeats_bullets), ("cyp2d6", cyp2d6, cyp2d6_bullets), ("hla", hla, hla_bullets), ("smn", smn, smn_bullets))
