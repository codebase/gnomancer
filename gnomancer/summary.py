"""Digestible summaries: a deterministic "at a glance" built from the results, a prompt any AI agent can
turn into a narrative, and an optional built-in Claude call to do that.

The narrative is the only part of gnomancer that may send data off the machine, and only when the user
passes --ai. What leaves is the structured result summary (haplogroups, panel readings, ClinVar and
PharmCAT tables), never the raw genotype file.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

POP_NAMES = {
    "ACB": "African Caribbean in Barbados", "ASW": "African ancestry in the south-west USA", "ESN": "Esan in Nigeria",
    "GWD": "Gambian in Western Division", "LWK": "Luhya in Kenya", "MSL": "Mende in Sierra Leone", "YRI": "Yoruba in Nigeria",
    "CLM": "Colombian in Medellín", "MXL": "Mexican ancestry in Los Angeles", "PEL": "Peruvian in Lima", "PUR": "Puerto Rican",
    "CDX": "Dai Chinese in Xishuangbanna", "CHB": "Han Chinese in Beijing", "CHS": "Southern Han Chinese", "JPT": "Japanese in Tokyo",
    "KHV": "Kinh in Ho Chi Minh City", "CEU": "Northern and Western European ancestry in Utah", "FIN": "Finnish", "GBR": "British",
    "IBS": "Iberian populations in Spain", "TSI": "Toscani in Italy", "BEB": "Bengali in Bangladesh", "GIH": "Gujarati Indian in Houston",
    "ITU": "Indian Telugu in the UK", "PJL": "Punjabi in Lahore", "STU": "Sri Lankan Tamil in the UK",
    "AFR": "African", "AMR": "Admixed American", "EAS": "East Asian", "EUR": "European", "SAS": "South Asian",
}

# Drugs per PharmCAT gene, from CPIC/DPWG guideline scope. Short, for orientation only.
PGX_DRUGS = {
    "CYP2C19": "clopidogrel, PPIs (omeprazole etc.), citalopram/escitalopram/sertraline, voriconazole",
    "CYP2C9": "warfarin, NSAIDs (celecoxib, ibuprofen…), phenytoin, siponimod",
    "VKORC1": "warfarin", "CYP4F2": "warfarin", "CYP3A5": "tacrolimus", "CYP3A4": "tacrolimus, statins (minor)",
    "CYP2B6": "efavirenz, bupropion, methadone", "CYP2D6": "codeine, tramadol, many antidepressants and antipsychotics, tamoxifen",
    "SLCO1B1": "statins (simvastatin, atorvastatin, rosuvastatin)", "ABCG2": "rosuvastatin, allopurinol",
    "DPYD": "5-fluorouracil, capecitabine", "TPMT": "azathioprine, mercaptopurine, thioguanine", "NUDT15": "azathioprine, mercaptopurine, thioguanine",
    "UGT1A1": "irinotecan, atazanavir, belinostat", "NAT2": "isoniazid, hydralazine, sulfasalazine",
    "G6PD": "rasburicase, primaquine, dapsone, sulfonamides, nitrofurantoin", "RYR1": "volatile anaesthetics and succinylcholine (malignant hyperthermia)",
    "CACNA1S": "volatile anaesthetics and succinylcholine (malignant hyperthermia)", "CFTR": "ivacaftor (CF modulator)",
    "F5": "estrogen-containing contraceptives and hormone therapy (clotting risk)", "F2": "estrogen-containing contraceptives and hormone therapy (clotting risk)",
    "IFNL3": "peginterferon alfa", "MT-RNR1": "aminoglycoside antibiotics", "HLA-B": "abacavir, carbamazepine, allopurinol",
}

NORMAL_PHENOS = {"normal metabolizer", "normal", "normal function", "normal risk", "b (reference)", ""}


def _clean_pheno(p: str) -> str:
    p = (p or "").strip().strip("/").strip()
    return p


def tidy_result(d: str, limit: int = 3) -> str:
    """'*1/*5 OR *1/*15 OR *1/*40 OR ...' -> first few alternatives plus a count."""
    alts = [a.strip() for a in str(d).split(" OR ")]
    if len(alts) <= limit:
        return " or ".join(alts)
    return " or ".join(alts[:limit]) + f" (+{len(alts) - limit} more possible)"


def tidy_pheno(p: str) -> str:
    """'Decreased Function / Decreased Function / Poor Function' -> 'Decreased Function or Poor Function'.
    'Indeterminate' is dropped when a definite phenotype is also possible."""
    seen, out = set(), []
    for x in [x.strip() for x in str(p).split("/")]:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    if len(out) > 1:
        out = [x for x in out if x.lower() != "indeterminate"]
    return " or ".join(out) if out else "—"


def at_a_glance(loaded, qc, roh, roh_text, combos, y, mt, ancestry, panel: pd.DataFrame, clinvar: dict, pharmcat: dict) -> dict:
    """Return {'who': [...], 'notable': [...], 'reassuring': [...], 'next': [...]} of short plain-language bullets."""
    who, notable, reassuring, nxt = [], [], [], []
    combo = {k: (v, n) for k, v, n in combos}
    # Who
    if y.get("haplogroup") and not y.get("error"):
        yl = (y.get("yleaf") or {}).get("trees") or {}
        best = yl.get("yfull", {}).get("haplogroup", "").split("*")[0] if yl else ""
        who.append(f"Father's line: Y haplogroup {y['haplogroup']}" + (f" (YFull: {best})" if best else "") + ".")
    if mt.get("haplogroup") and not mt.get("error"):
        who.append(f"Mother's line: mitochondrial haplogroup {mt['haplogroup']}.")
    if ancestry:
        sp = ancestry.get("superpopulation_code")
        pp = ancestry.get("population_code")
        spp = ancestry.get("superpopulation_percent")
        who.append(f"Closest 1000 Genomes group: {POP_NAMES.get(pp, pp)} ({POP_NAMES.get(sp, sp)}"
                   + (f", {spp:.0%} of neighbours" if isinstance(spp, (int, float)) else "") + "). A classification, not an admixture breakdown.")
    if "ABO blood group" in combo:
        who.append(f"Blood group {combo['ABO blood group'][0].split(' (rs')[0]} (Rh not typed).")
    if "Lactase persistence" in combo:
        who.append("Lactose: " + combo["Lactase persistence"][0].split(" (")[0] + ".")
    if "TAS2R38 bitter taste (PTC/PROP)" in combo:
        who.append("Bitter taste: " + combo["TAS2R38 bitter taste (PTC/PROP)"][0] + ".")
    # Notable: homozygous risk, carriers, PGx non-normal, ClinVar high-confidence
    if not panel.empty:
        p = panel[~panel.genotype.isin(["not on array", "no call"])]
        for _, r in p[p.category.isin(["health", "carrier"])].iterrows():
            n = str(r.effect_count).split("/")[0]
            if n == "2" and r.category == "health":
                if "protective" in str(r.interpretation).lower():
                    reassuring.append(f"{r.gene}: {r.interpretation}.")
                else:
                    notable.append(f"{r.gene} {r.label}: {r.interpretation}.")
            elif n == "1" and r.category == "carrier":
                notable.append(f"{r.gene} {r.label}: carrier ({r.interpretation}).")
            elif n == "2" and r.category == "carrier":
                notable.append(f"{r.gene} {r.label}: two copies ({r.interpretation}).")
        hets = p[(p.category == "health") & (p.effect_count.astype(str).str.startswith("1")) & ~p.interpretation.str.contains("benign", case=False)]
        if not hets.empty:
            notable.append("One copy of a risk allele at: " + ", ".join(f"{r.gene} ({r.label})" for _, r in hets.iterrows()) + ". Population-average territory, listed for completeness.")
    if "APOE" in combo:
        v = combo["APOE"][0].split(" (")[0]
        (reassuring if v == "e3/e3" else notable).append(f"APOE {v}: {combo['APOE'][1]}.")
    hits = clinvar.get("hits")
    if hits is not None and not hits.empty:
        hi = hits[hits.confidence == "high"]
        lo = len(hits) - len(hi)
        for _, r in hi.iterrows():
            conds = [c.strip() for c in str(r.condition).split(";")]
            cond = next((c for c in conds if "related disorder" not in c.lower()), conds[0])
            notable.append(f"ClinVar: {r.zygosity} for a reported pathogenic {r.gene} variant ({cond}). Needs clinical confirmation.")
        if lo:
            reassuring.append(f"{lo} other ClinVar 'pathogenic' hit(s) were flagged as probable array artifacts and should be ignored unless confirmed.")
    elif clinvar.get("available"):
        reassuring.append("No pathogenic ClinVar allele was found at any typed site.")
    calls = pharmcat.get("calls")
    if calls is not None and not calls.empty:
        odd = [r for _, r in calls.iterrows() if _clean_pheno(r.phenotype).lower() not in NORMAL_PHENOS
               and "reference" not in str(r.diplotype).lower() and "non-responsive" not in str(r.phenotype).lower()]
        for r in odd:
            note = " (the common European genotype; standard dosing)" if r.gene == "CYP3A5" and "*3/*3" in str(r.diplotype) else ""
            n_alt = len(str(r.diplotype).split(" OR "))
            result = tidy_result(r.diplotype, limit=2) if n_alt <= 2 else "several possible diplotypes (the chip cannot separate them)"
            notable.append(f"{r.gene} {result}: {tidy_pheno(r.phenotype)}{note}" + (f"; affects {PGX_DRUGS[r.gene]}" if r.gene in PGX_DRUGS else "") + ".")
        odd_genes = {r.gene for r in odd}
        normal = [r.gene for _, r in calls.iterrows() if r.gene not in odd_genes and _clean_pheno(r.phenotype).lower() in NORMAL_PHENOS]
        if normal:
            reassuring.append("Normal pharmacogenomic result for " + ", ".join(sorted(normal)) + " at the sites typed.")
    # Reassuring: negatives in carrier/health
    if not panel.empty:
        neg = panel[(panel.category.isin(["health", "carrier"])) & (panel.effect_count.astype(str).str.startswith("0"))
                    & ~panel.interpretation.astype(str).str.lower().str.contains("protective")]
        if len(neg):
            reassuring.append(f"Negative at {len(neg)} checked health and carrier sites, including " + ", ".join(neg.gene.drop_duplicates().head(8)) + ".")
    reassuring.append(roh_text.split(":")[0] + "." if roh_text else "")
    if not qc["flags"]:
        reassuring.append("The file itself passed every quality check.")
    # Next steps
    if hits is not None and not hits.empty and (hits.confidence == "high").any():
        nxt.append("If any ClinVar carrier result matters to you (for example before having children), ask a clinician for a confirmatory test of that specific variant.")
    nxt.append("Keep the pharmacogenomic table handy for any prescriber; a clinical PGx panel confirms it and covers the sites the chip misses.")
    if y.get("haplogroup") and not y.get("error"):
        nxt.append("For a deeper paternal-line result, Y sequencing (e.g. FTDNA Big Y) resolves far below what an array can.")
    return dict(who=who, notable=[n for n in notable if n], reassuring=[r for r in reassuring if r], next=nxt)


def build_prompt(glance: dict, combos, y, mt, ancestry, panel: pd.DataFrame, clinvar: dict, pharmcat: dict, qc: dict) -> str:
    """A self-contained prompt an AI agent can answer with a narrative summary."""
    facts = dict(
        quality=dict(no_call_rate=qc["no_call_rate"], sex=qc["sex_inferred"], flags=qc["flags"]),
        at_a_glance=glance,
        readings=[dict(name=k, value=v, note=n) for k, v, n in combos],
        y=dict(haplogroup=y.get("haplogroup"), ycc=y.get("ycc"), yleaf=(y.get("yleaf") or {}).get("trees")),
        mt=dict(haplogroup=mt.get("haplogroup"), quality=mt.get("quality")),
        ancestry=ancestry,
        panel=[dict(gene=r.gene, label=r.label, effect_count=r.effect_count, reading=r.interpretation, category=r.category)
               for _, r in panel.iterrows() if r.genotype not in ("not on array", "no call")],
        clinvar=(clinvar.get("hits").to_dict("records") if clinvar.get("hits") is not None and not clinvar["hits"].empty else []),
        pharmcat=(pharmcat.get("calls").to_dict("records") if pharmcat.get("calls") is not None and not pharmcat["calls"].empty else []),
    )
    return (
        "You are summarizing a consumer DNA report for the person it belongs to. Write 250-400 words of plain English, "
        "second person, warm but factual, no headings, no bullet lists, no medical advice. Cover, in this order: where the "
        "two lineages and the population estimate point; what is genuinely notable (two-copy risk alleles, carrier states, "
        "non-normal drug-metabolism results, high-confidence ClinVar hits) with the odds stated the way the data states them; "
        "what is reassuring; and what, if anything, is worth a conversation with a clinician. Say explicitly that a genotyping "
        "array is probabilistic and that rare-site calls need confirmation. Do not invent facts that are not in the data. "
        "Do not mention rsIDs. Do not use the words 'journey' or 'tapestry'.\n\n"
        "DATA (JSON):\n" + json.dumps(facts, indent=1, default=str)
    )


def ai_summary(prompt: str, model: str = "claude-opus-5-5") -> str:
    """Call Claude to write the narrative. Requires the anthropic SDK and credentials (ANTHROPIC_API_KEY or `ant auth login`)."""
    import anthropic  # imported lazily so the rest of gnomancer has no dependency on it

    client = anthropic.Anthropic()
    response = client.beta.messages.create(
        model=model,
        max_tokens=4000,
        output_config={"effort": "medium"},
        betas=["server-side-fallback-2026-07-01"],
        fallbacks="default",
        system="You write clear, careful summaries of personal genetic reports. You never diagnose.",
        messages=[{"role": "user", "content": prompt}],
    )
    if response.stop_reason == "refusal":
        raise RuntimeError("the model declined to summarize this report")
    return "".join(b.text for b in response.content if b.type == "text").strip()


def write_prompt_file(out: Path, prompt: str) -> Path:
    p = out / "summary-prompt.md"
    p.write_text("# Prompt for an AI agent to write the narrative summary\n\n"
                 "Paste the whole of this file to your assistant, save its answer as `summary.md`, then re-run gnomancer with "
                 "`--summary summary.md` to embed it at the top of the report. Only the structured results below are shared, "
                 "not the raw genotype file.\n\n```\n" + prompt + "\n```\n")
    return p
