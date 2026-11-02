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


def _conditions(raw: str, limit: int = 4) -> str:
    """'Disease_A|not_provided|Disease_B' -> 'Disease A; Disease B' (deduplicated, generic entries dropped, at most `limit`)."""
    seen, out = set(), []
    for c in raw.replace("_", " ").split("|"):
        c = c.strip()
        if not c or c.lower() in ("not provided", "not specified", "see cases") or c in seen:
            continue
        seen.add(c)
        out.append(c)
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
            if chrom == "MT":
                chrom = "MT"
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
            if zyg != "heterozygous" and (af is None or af < 0.01):
                flags.append("homozygous for a rare allele: statistically far more likely a genotyping/strand artifact than a true homozygote")
            if "conflicting" in rev:
                flags.append("submitters disagree on classification")
            confidence = "low" if (ambiguous or (zyg != "heterozygous" and (af is None or af < 0.01)) or STAR.get(rev, 0) == 0) else ("medium" if (af or 0) >= 0.05 else "high")
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
