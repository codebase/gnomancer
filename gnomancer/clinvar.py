"""Scan every typed site against ClinVar (GRCh37 VCF) and report pathogenic / likely pathogenic alleles carried.

Why this matters: the curated panel covers ~110 sites; consumer chips type several thousand sites that ClinVar
classifies. This finds the rest. Why to be careful: arrays have a high false-positive rate at rare sites
(Tandy-Connor et al. 2018, Genetics in Medicine, doi:10.1038/gim.2018.38: 40% of variants from consumer raw data
were false positives on clinical confirmation), so every hit here is "worth confirming with a clinical test",
never a diagnosis.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from .paths import tool

DEFAULT_VCF = tool("clinvar", "clinvar_GRCh37.vcf.gz")
PATHOGENIC = ("Pathogenic", "Likely_pathogenic")
STAR = {"practice_guideline": 4, "reviewed_by_expert_panel": 3, "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_single_submitter": 1, "criteria_provided,_conflicting_classifications": 1,
        "criteria_provided,_conflicting_interpretations": 1, "no_assertion_criteria_provided": 0, "no_classification_provided": 0}


VAGUE = ("inborn genetic diseases", "fetal anomalies with a likely genetic cause", "rare genetic disease", "abnormality of", "disorder")


def _conditions(raw: str, limit: int = 4) -> str:
    """'Disease_A|not_provided|Disease_B' -> 'Disease A; Disease B' (deduplicated, empty entries dropped, catch-all
    labels moved behind named conditions, at most `limit`)."""
    seen, out = set(), []
    for c in raw.replace("_", " ").split("|"):
        c = c.strip()
        if not c or c.lower() in ("not provided", "not specified", "see cases") or c in seen:
            continue
        seen.add(c)
        out.append(c)
    # Named diseases first; catch-all labels, susceptibility entries and OMIM-style all-caps trait names after them.
    out.sort(key=lambda c: (c.lower().startswith(VAGUE[:4]) or c.lower() == "disorder" or c.lower().endswith("-related disorder")
                            or "susceptibility" in c.lower() or c.isupper()))
    return "; ".join(out[:limit]) + (f" (+{len(out) - limit} more)" if len(out) > limit else "")


def _info(s: str) -> dict:
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def scan(df: pd.DataFrame, vcf: Path = DEFAULT_VCF) -> dict:
    """df: normalized SNP table (chrom, pos, genotype). Returns dict with hits table and coverage stats."""
    res = dict(available=vcf.exists(), vcf=str(vcf), sites_in_clinvar=0, pathogenic_sites_typed=0, hits=pd.DataFrame(), note="",
               carried={}, release=None)
    if not vcf.exists():
        res["note"] = "ClinVar VCF not found; run `gnomancer setup` (downloads ~200 MB from NCBI)."
        return res
    called = df[df.genotype.notna()]
    key = {}
    for chrom, pos, gt in zip(called.chrom, called.pos, called.genotype):
        key[(chrom, int(pos))] = gt
    hits = []
    carried_sites = {}
    n_in = n_path = 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                if line.startswith("##fileDate="):
                    res["release"] = line.strip().split("=", 1)[1]
                continue
            p = line.rstrip("\n").split("\t", 8)
            chrom = p[0]
            k = (chrom, int(p[1]))
            gt = key.get(k)
            if gt is None:
                continue
            n_in += 1
            info = _info(p[7])
            sig = info.get("CLNSIG", "")
            ref_, alts_ = p[3], p[4].split(",")
            if len(ref_) == 1 and all(len(a) == 1 for a in alts_) and any(a in alts_ for a in gt) and all(a in alts_ or a == ref_ for a in gt):
                carried_sites[f"{chrom}:{p[1]}:{p[4]}"] = dict(rsid=("rs" + info["RS"]) if info.get("RS") else p[2], gene=info.get("GENEINFO", "").split(":")[0],
                                                       sig=sig, rev=info.get("CLNREVSTAT", ""), gt=gt)
            if not any(s in sig for s in PATHOGENIC) or sig.startswith("Benign") or sig.startswith("Likely_benign"):
                continue
            n_path += 1
            ref, alts = p[3], p[4].split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts):
                continue  # indels are not reliably comparable to array I/D calls
            carried = [a for a in gt if a in alts]
            if not carried:
                continue
            if any(a not in alts and a != ref for a in gt):
                continue  # allele mismatch (likely strand issue); do not report
            zyg = "homozygous" if len(gt) == 2 and gt[0] == gt[1] else ("hemizygous" if len(gt) == 1 else "heterozygous")
            rev = info.get("CLNREVSTAT", "")
            # Allele frequency of the ALT allele, from whichever population resource ClinVar carries.
            af = None
            for k_ in ("AF_EXAC", "AF_TGP", "AF_ESP"):
                if k_ in info:
                    try:
                        af = max(float(x) for x in info[k_].split(","))
                    except ValueError:
                        pass
                    break
            flags = []
            ambiguous = {ref, alts[0]} in ({"A", "T"}, {"C", "G"})
            if ambiguous:
                flags.append("strand-ambiguous SNP (A/T or C/G): a homozygous-reference call reported on the other strand looks exactly like this")
            if af is not None and af >= 0.05:
                flags.append(f"ALT allele is common (AF {af:.1%}): classification is likely outdated or low-penetrance")
            if zyg == "homozygous" and (af is None or af < 0.01):
                flags.append("homozygous for a rare allele: statistically far more likely a genotyping/strand artifact than a true homozygote")
            if "conflicting" in rev:
                flags.append("submitters disagree on classification")
            confidence = "low" if (ambiguous or (zyg == "homozygous" and (af is None or af < 0.01)) or STAR.get(rev, 0) == 0) else ("medium" if (af or 0) >= 0.05 else "high")
            hits.append(dict(rsid=("rs" + info["RS"]) if info.get("RS") else p[2], chrom=chrom, pos=k[1], genotype=gt,
                             ref=ref, alt=",".join(alts), zygosity=zyg, gene=info.get("GENEINFO", "").split(":")[0],
                             significance=sig.replace("_", " "), review=rev.replace("_", " "), stars=STAR.get(rev, 0),
                             condition=_conditions(info.get("CLNDN", "")),
                             variant=info.get("CLNHGVS", ""), alt_af=af, confidence=confidence, flags="; ".join(flags)))
    res["sites_in_clinvar"], res["pathogenic_sites_typed"] = n_in, n_path
    res["carried"] = carried_sites
    out = pd.DataFrame(hits)
    if not out.empty:
        out["conf_rank"] = out.confidence.map({"high": 0, "medium": 1, "low": 2})
        out = out.sort_values(["conf_rank", "stars"], ascending=[True, False]).drop(columns="conf_rank")
    res["hits"] = out
    return res


def interpret(res: dict) -> str:
    if not res["available"]:
        return res["note"]
    h = res["hits"]
    base = (f"{res['sites_in_clinvar']:,} typed sites have a ClinVar record; {res['pathogenic_sites_typed']:,} of those are "
            f"classified pathogenic or likely pathogenic for at least one submitter. ")
    if h.empty:
        return base + "No pathogenic or likely-pathogenic alternate allele was called at any of them."
    hi = h[h.confidence == "high"]
    return base + (f"{len(h)} site(s) carry a pathogenic/likely-pathogenic allele; {len(hi)} pass the plausibility checks "
                   "(heterozygous or common, unambiguous strand, reviewed). The rest are flagged as probable artifacts. "
                   "Array calls at rare sites are often wrong: when a clinical lab re-tested variants reported in consumer raw data, "
                   "40% were false positives (Tandy-Connor et al., Genetics in Medicine 2018, doi:10.1038/gim.2018.38). Treat every row as "
                   "'confirm with a clinical test', not as a diagnosis. Heterozygous hits for recessive conditions are carrier states.")


# ---------------------------------------------------------------------------------------------------------------
# Sequencing VCFs
#
# The array scan above compares letters at typed sites and has to skip indels. A sequencing VCF gives exact
# alleles, so here a ClinVar record matches when chromosome, position, reference and alternate allele are all
# equal after trimming. ClinVar's VCF is left-aligned (checked: none of its ~318,000 indels can be shifted left)
# and so is caller output, so equal variants are written the same way.

ACMG_FILE = Path(__file__).resolve().parent / "data" / "acmg_sf.tsv"
LOF = ("frameshift_variant", "nonsense", "splice_donor_variant", "splice_acceptor_variant", "initiator_codon_variant", "stop_lost")
BENIGN = ("Benign", "Likely_benign", "Benign/Likely_benign")


def vcf_for(build: int) -> Path:
    return tool("clinvar", f"clinvar_GRCh{build}.vcf.gz")


def acmg_genes(path: Path = ACMG_FILE) -> dict:
    import csv
    with open(path, newline="") as f:
        return {r["gene"]: r for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")}


def positions(build: int, vcf: Path | None = None) -> set:
    """Every position ClinVar has a record at, as seqvcf keys. Cached next to the VCF (rebuilt when it changes)."""
    import numpy as np
    from .seqvcf import CODE, trim
    vcf = vcf or vcf_for(build)
    cache = vcf.with_suffix(".positions.npy")
    if cache.exists() and cache.stat().st_mtime >= vcf.stat().st_mtime:
        return set(np.load(cache).tolist())
    keys = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 5)
            code = CODE.get(p[0])
            if code is None:
                continue
            pos = int(p[1])
            if len(p[3]) > 1 and len(p[4]) > 1:
                pos = trim(pos, p[3], p[4])[0]
            keys.append((code << 28) | pos)
    arr = np.unique(np.array(keys, dtype=np.int64))
    try:
        np.save(cache, arr)
    except OSError:
        pass
    return set(arr.tolist())


def _allele(ref: str, alt: str) -> str:
    """Short description of a change for a table cell."""
    if len(ref) == 1 and len(alt) == 1:
        return f"{ref}>{alt}"
    if len(alt) == 1 and ref[0] == alt:
        return f"del{ref[1:]}" if len(ref) <= 13 else f"del {len(ref) - 1} bp"
    if len(ref) == 1 and alt[0] == ref:
        return f"ins{alt[1:]}" if len(alt) <= 13 else f"ins {len(alt) - 1} bp"
    return f"{ref}>{alt}" if len(ref) + len(alt) <= 16 else f"{len(ref)} bp replaced by {len(alt)} bp"


def scan_seq(d, vcf: Path | None = None) -> dict:
    """Match the calls kept in a SeqData against ClinVar on the file's own build."""
    from .seqsites import expected_ploidy
    from .seqvcf import CODE, trim
    vcf = vcf or vcf_for(d.build)
    res = dict(available=vcf.exists(), vcf=str(vcf), sites_in_clinvar=0, pathogenic_sites_typed=0, hits=pd.DataFrame(), note="",
               carried={}, release=None, kind="wgs", pathogenic_records=0, failed_quality=0)
    if not vcf.exists():
        res["note"] = f"ClinVar GRCh{d.build} VCF not found; run `gnomancer setup` (downloads ~200 MB from NCBI)."
        return res
    acmg = acmg_genes()
    calls = d.calls
    hits, carried = [], {}
    n_match = n_path_total = 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                if line.startswith("##fileDate="):
                    res["release"] = line.strip().split("=", 1)[1]
                continue
            p = line.split("\t", 8)
            code = CODE.get(p[0])
            if code is None:
                continue
            pos, ref, alt = int(p[1]), p[3], p[4]
            if len(ref) > 1 and len(alt) > 1:
                pos, ref, alt = trim(pos, ref, alt)
            here = calls.get((code << 28) | pos)
            path = "athogenic" in p[7] and ("CLNSIG=Pathogenic" in p[7] or "CLNSIG=Likely_pathogenic" in p[7])
            n_path_total += path
            if not here:
                continue
            c = next((c for c in here if c.ref == ref and c.alt == alt), None)
            if c is None:
                continue
            n_match += 1
            info = _info(p[7].rstrip("\n"))
            sig, rev = info.get("CLNSIG", ""), info.get("CLNREVSTAT", "")
            gene = info.get("GENEINFO", "").split(":")[0] or "(no gene)"
            rsid = ("rs" + info["RS"].split("|")[0]) if info.get("RS") else p[2]
            if c.ok and sig not in BENIGN:
                carried[f"{p[0]}:{p[1]}:{p[3]}:{p[4]}"] = dict(rsid=rsid, gene=gene, sig=sig, rev=rev, gt=c.zygosity)
            if not path:
                continue
            if not c.ok:
                res["failed_quality"] += 1
            zyg = c.zygosity
            if c.chrom == "MT":
                share = c.ab if c.ab is not None else (1.0 if c.copies == c.ploidy else 0.5)
                zyg = "homoplasmic" if share >= 0.9 else f"heteroplasmic ({share:.0%} of reads)"
            elif zyg == "homozygous" and expected_ploidy(c.chrom, c.pos, d.sex, d.build) == 1:
                zyg = "hemizygous"
            af = None
            for k_ in ("AF_EXAC", "AF_TGP", "AF_ESP"):
                if k_ in info:
                    try:
                        af = max(float(x) for x in info[k_].split(","))
                    except ValueError:
                        pass
                    break
            cons = ", ".join(sorted({m.split("|", 1)[1].replace("_", " ") for m in info.get("MC", "").split(",") if "|" in m}))
            stars = STAR.get(rev, 0)
            a = acmg.get(gene)
            flags = []
            if not c.ok:
                flags.append(f"weak call: {c.why}")
            if af is not None and af >= 0.05:
                flags.append(f"allele is common (AF {af:.1%}): classification is likely outdated or low-penetrance")
            if "low_penetrance" in sig or "risk_allele" in sig:
                flags.append("classified as low-penetrance or a risk allele, not a high-risk variant")
            if stars == 0:
                flags.append("no review criteria provided by the submitter")
            if c.other_alt:
                flags.append("the other copy carries a different alternate allele at this position")
            if a and a["note"]:
                flags.append(a["note"])
            confidence = "low" if (not c.ok or stars == 0) else ("medium" if (af or 0) >= 0.05 or "low_penetrance" in sig or "risk_allele" in sig else "high")
            role = ""
            if a:
                inh = a["inheritance"]
                if inh == "AR":
                    role = "actionable gene, recessive: two copies" if zyg == "homozygous" else "actionable gene, recessive: carrier unless a second variant is on the other copy"
                elif inh == "XL":
                    role = "actionable gene, X-linked" + (": hemizygous" if zyg == "hemizygous" else ": one copy")
                else:
                    role = "actionable gene, dominant: one copy is enough"
            hits.append(dict(rsid=rsid, chrom=c.chrom, pos=c.pos, genotype=_allele(c.ref, c.alt), ref=c.ref, alt=c.alt, zygosity=zyg,
                             gene=gene, significance=sig.replace("_", " "), review=rev.replace("_", " "), stars=stars,
                             condition=_conditions(info.get("CLNDN", "")), variant=info.get("CLNHGVS", ""), alt_af=af,
                             confidence=confidence, flags="; ".join(flags), consequence=cons, lof=any(t in info.get("MC", "") for t in LOF),
                             depth=c.dp, gq=c.gq, allele_share=(round(c.ab, 2) if c.ab is not None else None),
                             acmg=bool(a), acmg_condition=(a["condition"] if a else ""), role=role, phase_set=c.ps, hap=c.hap, clinvar_id=p[2]))
    res["sites_in_clinvar"], res["pathogenic_sites_typed"], res["pathogenic_records"] = n_match, n_path_total, n_path_total
    res["carried"] = carried
    out = pd.DataFrame(hits)
    if not out.empty:
        # Two pathogenic variants in one gene matter for recessive conditions; say whether phasing places them together.
        for gene, grp in out[(out.confidence != "low") & (out.zygosity == "heterozygous") & (out.chrom != "MT")].groupby("gene"):
            if len(grp) < 2 or not gene:
                continue
            ps, hp = set(grp.phase_set), set(grp.hap)
            if len(ps) == 1 and "" not in ps and -1 not in hp:
                msg = ("on opposite copies (phased): both copies of the gene are affected" if len(hp) > 1
                       else "on the same copy (phased): the other copy is unaffected")
            else:
                msg = "not phased: on opposite copies they would affect both copies of the gene; testing a parent settles it"
            for i in grp.index:
                out.at[i, "flags"] = "; ".join(x for x in (out.at[i, "flags"], f"{len(grp)} pathogenic variants in {gene}, {msg}") if x)
        out["rank"] = out.confidence.map({"high": 0, "medium": 1, "low": 2}) * 2 + (~out.acmg).astype(int)
        out = out.sort_values(["rank", "stars"], ascending=[True, False]).drop(columns=["rank", "hap"])
    res["hits"] = out
    return res


def interpret_seq(res: dict, whole_genome: bool = True) -> str:
    if not res["available"]:
        return res["note"]
    h = res["hits"]
    base = (f"Every variant in the file was compared with ClinVar's {res['pathogenic_records']:,} pathogenic or likely pathogenic records, "
            f"by exact position and alleles, insertions and deletions included. {res['sites_in_clinvar']:,} of your variants have a ClinVar record of any kind. ")
    if h.empty:
        return base + ("None is classified pathogenic or likely pathogenic. " +
                       ("" if whole_genome else "The file is not a whole genome, so genes outside it were not examined. ") +
                       "This rules out known small variants at covered positions; it does not rule out large deletions, repeat expansions or variants ClinVar has not classified.")
    groups = {k: len(df) for k, _, _, df in seq_groups(h)}
    return base + (f"{len(h)} carry a pathogenic or likely pathogenic classification: {groups.get('actionable', 0)} in one of the 84 genes the ACMG considers "
                   f"medically actionable, {groups.get('both', 0)} elsewhere on both copies, {groups.get('single', 0)} as single copies, "
                   + (f"{groups['mito_low']} low-level mitochondrial, " if groups.get("mito_low") else "") +
                   f"{groups.get('risk', 0)} common or low-penetrance, and {groups.get('aside', 0)} set aside as unreliable. "
                   "Most healthy people carry a few: usually single copies of recessive variants, which make you a carrier and do not cause disease. "
                   "A classification can also be wrong or apply only in some families, so confirm anything that matters with a clinical test before acting on it.")


BOTH_COPIES = "homozygous|hemizygous|homoplasmic|heteroplasmic"
MT_LOW_SHARE = 0.10


def seq_groups(hits: pd.DataFrame) -> list:
    """Split sequencing hits by what they mean: [(key, title, lead, rows)]. One place, so the summary, the HTML
    report and the Markdown report agree. Empty groups are left out."""
    if hits is None or hits.empty:
        return []
    # A mitochondrial variant in a small share of reads is its own case: disease needs a high share of mutant copies.
    share = pd.to_numeric(hits.allele_share, errors="coerce")
    low_mt = (hits.chrom.astype(str) == "MT") & (share < MT_LOW_SHARE) & (hits.confidence != "low")
    rest = hits[~low_mt]
    usable = rest[rest.confidence != "low"]
    good = rest[rest.confidence == "high"]
    both = good.zygosity.astype(str).str.contains(BOTH_COPIES)
    groups = [
        ("actionable", "In a medically actionable gene",
         "One of the 84 genes on the ACMG secondary-findings list, where a pathogenic variant changes medical care. The notes say whether one copy is enough.",
         usable[usable.acmg]),
        ("both", "On both copies, or the only copy",
         "Homozygous, hemizygous or mitochondrial. For a recessive condition this is the genotype that causes it, if the call and the classification are both right.",
         good[~good.acmg & both]),
        ("single", "Single copies",
         "One copy of a pathogenic variant. For recessive conditions this is a carrier state: it does not cause disease, and matters mainly when both parents "
         "carry a variant in the same gene. A few conditions here may be dominant; the condition name is the place to start.",
         good[~good.acmg & ~both]),
        ("mito_low", "Low-level mitochondrial variants",
         f"Carried by under {MT_LOW_SHARE:.0%} of mitochondrial reads. Mitochondrial variants generally cause disease only at high shares, well over half; "
         "low-level ones are common, differ between tissues, and can be sequencing artifacts. Worth raising only alongside matching symptoms.",
         hits[low_mt]),
        ("risk", "Common or low-penetrance alleles",
         "Real calls, but the allele is common in the population or is classified as low-penetrance or a risk factor: it shifts odds and does not, on its own, predict disease.",
         rest[(rest.confidence == "medium") & ~rest.acmg]),
        ("aside", "Set aside",
         "Matches that should not be relied on: a weak call, or a classification no one has reviewed. Listed so nothing is hidden.",
         hits[hits.confidence == "low"]),
    ]
    return [g for g in groups if not g[3].empty]
