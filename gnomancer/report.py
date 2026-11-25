"""Markdown report writer."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .panel import CATEGORY_ORDER, POP_LABEL, fmt_share

DISCLAIMER = (
    "Computed locally from a consumer genotyping array. The array reports a fixed list of common sites, "
    "not your whole genome; absence of a variant here never means a gene is normal, and single-marker "
    "readings are probabilistic. This is not medical advice. Discuss anything concerning with a clinician "
    "or genetic counselor and confirm with a clinical-grade test before acting on it."
)


def _md_table(df: pd.DataFrame, cols: list[str]) -> str:
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |")
    return "\n".join(lines)


def write(out: Path, loaded, qc: dict, chr_table: pd.DataFrame, roh: pd.DataFrame, roh_text: str,
          panel: pd.DataFrame, combos: list, y: dict, mt: dict, ancestry: dict | None, extra: dict | None = None) -> Path:
    extra = extra or {}
    seq = getattr(loaded, "kind", "array") != "array"
    if seq:
        from .clinvar import seq_groups as clinvar_groups
        from .report_html import SEQ_DISCLAIMER, seq_limits
    L = []
    L.append(f"# {'Genome' if seq else 'Raw DNA'} analysis report\n\n{SEQ_DISCLAIMER if seq else DISCLAIMER}\n")
    if extra.get("narrative"):
        L.append("## In plain language\n\n" + extra["narrative"] + "\n")
    g = extra.get("glance") or {}
    if g:
        L.append("## At a glance\n")
        for title, key in (("Who the data says you are", "who"), ("Worth knowing", "notable"), ("Reassuring", "reassuring"), ("Sensible next steps", "next")):
            if g.get(key):
                L.append(f"**{title}**\n" + "\n".join(f"- {b}" for b in g[key]) + "\n")
    L.append("Sections: input · quality · homozygosity · paternal line · maternal line · ancestry · combined readings · marker panel · ClinVar · drug-gene results"
             + "".join(f" · {x.get('nav') or x['title']}".lower() for x in extra.get("sections") or []) + " · limits\n")
    L.append("## Input\n")
    L.append(f"- Vendor/format detected: **{loaded.source}**" + (f" (chip {loaded.chip} {loaded.chip_version})" if loaded.chip else ""))
    L.append(f"- Genome build: GRCh{loaded.build}" + (f" (remapped from build {loaded.build_original})" if loaded.remapped else "") +
             ("" if loaded.build_detected else " (assumed, not detected)"))
    L.append(f"- {'Passing variants' if seq else 'Sites loaded'}: {loaded.count:,}")
    for n in loaded.notes:
        L.append(f"- Note: {n}")
    L.append("\n## Quality\n")
    L.append("| Metric | Value |\n|---|---|")
    if seq:
        f2 = lambda v, fmt: format(v, fmt) if isinstance(v, (int, float)) and v == v else "n/a"
        L.append(f"| Layout | {'gVCF with reference blocks' if qc.get('gvcf') else 'variant-only VCF'}{' · caller ' + qc['caller'] if qc.get('caller') else ''} |")
        L.append(f"| Passing variants | {qc['snvs']:,} SNVs, {qc['indel_calls']:,} insertions/deletions |")
        L.append(f"| Failed quality gates | {qc['no_call']:,} ({qc['no_call_rate']:.1%}), not used in this report |")
        L.append(f"| Transition / transversion | {f2(qc['ti_tv'], '.2f')} (about 2.0-2.1 expected for a whole genome) |")
        L.append(f"| Het / hom ratio (autosomal SNVs) | {f2(qc['het_hom_ratio'], '.2f')} |")
        L.append(f"| Mean depth at variants | {f2(qc['mean_depth'], '.0f')}x |")
        L.append(f"| X heterozygosity | {f2(qc['x_het_rate'], '.3f')} |")
        L.append(f"| Genetic sex | {qc['sex_inferred']} |")
        rd = qc.get("reads")
        if rd:
            L.append(f"| Mean depth in the reads ({rd['file']}) | {rd['mean_depth']:.0f}x, {rd['read_length']}-base reads; mitochondrion about {rd['mt_depth']:,.0f}x |")
            L.append(f"| X and Y depth relative to autosomes | {rd['x_ratio']:.2f} and {rd['y_ratio']:.2f}: {rd['sex']} |")
            L.append(f"| Silent positions checked in the reads | {rd['sites_checked']:,}: {rd['sites_confirmed']:,} read as reference, {rd['sites_thin']:,} too thin, "
                     f"{rd['sites_disagree']:,} where the reads show a variant the VCF lacks |")
    else:
        L.append(f"| No-calls | {qc['no_call']:,} ({qc['no_call_rate']:.2%}) |")
        L.append(f"| Indel calls (I/D) | {qc['indel_calls']:,} |")
        L.append(f"| Autosomal heterozygosity | {qc['autosomal_het']:.3f} |")
        L.append(f"| X heterozygosity | {qc['x_het_rate']:.4f} |")
        L.append(f"| Y sites called | {qc['y_called']:,} / {qc['y_total']:,} |")
        L.append(f"| Genetic sex | {qc['sex_inferred']} (file header: {qc['sex_library'] or 'n/a'}) |")
    for f in qc["flags"]:
        L.append(f"\n⚠ {f}")
    L.append("\nPer-chromosome table: `chromosome-qc.tsv`.\n")
    L.append("## Runs of homozygosity\n")
    L.append(roh_text)
    if not roh.empty:
        L.append("\n" + _md_table(roh.head(15), ["chrom", "start", "end", "length_mb", "snps"]))
    L.append("\n## Paternal line (Y chromosome)\n")
    if seq and y.get("haplogroup") and not y.get("error"):
        L.append(f"**Haplogroup {y['haplogroup']}**. {y.get('basis') or 'Yleaf, from the sequenced Y chromosome.'}")
        names = dict(yfull="YFull", ftdna="FamilyTreeDNA", isogg="ISOGG")
        for tree, t in ((y.get("yleaf") or {}).get("trees") or {}).items():
            L.append(f"- {names.get(tree, tree)}: **{t['haplogroup']}** (defining marker {t['marker']}, {t['valid_markers']} valid markers, QC {t['qc']})")
    elif y.get("haplogroup") and not y.get("error"):
        L.append(f"**Haplogroup {y['haplogroup']}** (ISOGG long form {y['ycc']}), called by yhaplo from {y['n_y_snps']:,} Y sites.")
        L.append(f"\nDerived-SNP path: `{y['path']}`" if y.get("path") else "")
        yl = y.get("yleaf") or {}
        if yl.get("trees"):
            L.append("\nSecond opinion from Yleaf on current trees (" + f"{yl['n_y_snps']:,} Y sites). Each tree names branches its own way, so the labels differ even when they point at the same lineage:")
            names = dict(yfull="YFull v14", ftdna="FamilyTreeDNA", isogg="ISOGG")
            for tree, t in yl["trees"].items():
                L.append(f"- {names.get(tree, tree)}: **{t['haplogroup']}** (defining marker {t['marker']}, {t['valid_markers']} valid markers, QC {t['qc']})")
        elif yl.get("error"):
            L.append(f"\nYleaf second opinion not available: {yl['error']}")
        L.append("\nArrays type only a few thousand Y sites, so this is the deepest branch the chip can see; "
                 "Y sequencing (e.g. Big Y-700) is needed for a terminal subclade.")
    else:
        L.append(f"Not called: {y.get('error', 'unknown')}")
    L.append("\n## Maternal line (mitochondrial DNA)\n")
    if mt.get("sequenced") and mt.get("haplogroup") and not mt.get("error"):
        L.append(f"**Haplogroup {mt['haplogroup']}** (Haplogrep3 quality {mt['quality']:.2f}) from the whole mitochondrial sequence.")
        L.append(f"\n{mt['n_mt_snps']} differences from rCRS: {', '.join(mt['diffs'] + mt['indels']) or 'none'}.")
        if mt.get("heteroplasmies"):
            L.append("Carried by only part of the reads (heteroplasmy; left out of the haplogroup call): "
                     + ", ".join(f"{h['pos']}{h['ref']}>{h['alt']} at {h['share']:.0%}" for h in mt["heteroplasmies"][:20]) + ".")
    elif mt.get("haplogroup") and not mt.get("error"):
        L.append(f"**Haplogroup {mt['haplogroup']}** (Haplogrep3 quality {mt['quality']:.2f}) from {mt['n_mt_snps']} typed MT positions.")
        L.append(f"\nDifferences from rCRS at typed sites: {', '.join(mt['diffs']) or 'none'}.")
        if mt["indels"]:
            L.append(f"Indel probes ignored: {', '.join(mt['indels'])}.")
        L.append("\nOnly a few hundred MT sites are typed, so sub-branches usually cannot be resolved; full mtDNA sequencing would.")
    else:
        L.append(f"Not called: {mt.get('error', 'unknown')}")
    L.append("\n## Ancestry estimate (ezancestry)\n")
    if not ancestry:
        failed = [n for n in loaded.notes if n.startswith("ezancestry failed")]
        if failed:
            L.append(f"Requested but failed: {failed[0]}")
        else:
            L.append("Not run. Add `--ancestry` for a 1000 Genomes population classification "
                     "(`gnomancer setup` installs ezancestry into a side venv when Python 3.11 is available).")
    if ancestry:
        from .summary import POP_NAMES
        pp, sp = ancestry.get("population_code"), ancestry.get("superpopulation_code")
        L.append(f"- Nearest reference population: **{POP_NAMES.get(pp, pp)}** ({pp}); super-population **{POP_NAMES.get(sp, sp)}** ({sp})")
        if isinstance(ancestry.get("superpopulation_percent"), (int, float)):
            L.append(f"- Share of nearest neighbours in that super-population: {ancestry['superpopulation_percent']:.0%}; in the single nearest population: {ancestry.get('population_percent', 0):.0%}")
        probs = ancestry.get("superpopulation_probabilities") or {}
        if probs:
            L.append("- Super-population probabilities: " + ", ".join(f"{k} {v:.0%}" for k, v in sorted(probs.items(), key=lambda kv: -kv[1])))
        L.append("- Method: k-nearest-neighbours on 1000 Genomes PCA; a classifier into reference groups, not an admixture breakdown.")
    L.append("\n## Combined readings\n")
    for name, val, note in combos:
        L.append(f"- **{name}:** {val}" + (f" — {note}" if note else ""))
    L.append("\n## Marker panel\n")
    pop = extra.get("pop") or "EUR"
    L.append(("Effect-allele counts are on the forward strand. `basis` says how each genotype was arrived at: called from a variant record, or taken as reference because the VCF has no variant there. "
              if seq else "Effect-allele counts are on the GRCh37 forward strand. `not on array` means the vendor's chip does not type that site. ") +
             f"`share` is the fraction of the {POP_LABEL.get(pop, pop)} 1000 Genomes reference population expected to carry the same genotype "
             f"(Hardy-Weinberg from the effect-allele frequency); the TSV has all five super-populations.\n")
    titles = dict(ancestry="Ancestry-informative markers", trait="Traits", health="Health-associated common variants",
                  carrier="Carrier / rare pathogenic sites typed on arrays", pgx="Pharmacogenomic sites (single markers)")
    for cat in CATEGORY_ORDER:
        sub = panel[panel.category == cat]
        if sub.empty:
            continue
        L.append(f"\n### {titles.get(cat, cat).replace(' typed on arrays', '') if seq else titles.get(cat, cat)}\n")
        sub = sub.copy()
        col = f"gf_{pop}"
        sub["share"] = sub[col].map(lambda v: fmt_share(v) if isinstance(v, float) and v == v else "") if col in sub else ""
        L.append(_md_table(sub.astype(object).fillna(""), ["rsid", "gene", "label", "genotype", "effect_count", "share", "basis", "interpretation"] if seq and "basis" in sub
                           else ["rsid", "gene", "label", "genotype", "effect_count", "share", "interpretation"]))
    cv = extra.get("clinvar") or {}
    L.append(f"\n## ClinVar scan of every {'variant' if seq else 'typed site'}\n")
    L.append(extra.get("clinvar_text", "not run"))
    hits = cv.get("hits")
    if seq and hits is not None and not hits.empty:
        for _, title, lead, df in clinvar_groups(hits):
            L.append(f"\n### {title}\n\n{lead}\n\n" + _md_table(df.astype(object).fillna(""), ["gene", "genotype", "zygosity", "significance", "stars", "depth", "allele_share", "condition", "role", "flags"]))
    elif hits is not None and not hits.empty:
        L.append("\n" + _md_table(hits, ["rsid", "gene", "genotype", "zygosity", "significance", "stars", "confidence", "condition", "flags"]))
    pc = extra.get("pharmcat") or {}
    L.append("\n## Drug-gene results (PharmCAT)\n")
    L.append(extra.get("pharmcat_text", "not run"))
    card = extra.get("card") or {}
    if card.get("rows"):
        L.append(f"\nOne-page card for a pharmacist: `{Path(card['pdf'] or card['html']).name}`")
    calls = pc.get("calls")
    if calls is not None and not calls.empty:
        from .summary import PGX_DRUGS, pgx_cells
        c = calls.copy()
        cells = [pgx_cells(r) for _, r in c.iterrows()]
        c["result"] = [x[0] for x in cells]
        c["phenotype"] = [x[1] for x in cells]
        c["drugs"] = c.gene.map(PGX_DRUGS).fillna("")
        c["coverage"] = c.missing_positions.map({True: "partial", False: "complete"})
        L.append("\n" + _md_table(c, ["gene", "result", "phenotype", "drugs", "typed", "coverage"]))
        L.append("\nRYR1/CACNA1S 'Uncertain Susceptibility' means a variant of uncertain significance for malignant hyperthermia was seen, not a positive result. "
                 "The CFTR row describes ivacaftor response, not cystic fibrosis status.")
    for x in extra.get("sections") or []:
        L.append(f"\n## {x['title']}\n")
        if x.get("lead"):
            L.append(x["lead"] + "\n")
        if x.get("text"):
            L.append(x["text"])
        if x.get("rows"):
            L.append("\n" + _md_table(pd.DataFrame(x["rows"], columns=x["columns"]), list(x["columns"])))
        for n in x.get("notes") or []:
            L.append(f"\n{n}")
    L.append("\n## What this file cannot tell you\n")
    if seq:
        L.extend(f"- {t}" for t in seq_limits(loaded, qc, extra, extra.get("pop") or "EUR"))
        out.write_text("\n".join(L) + "\n")
        return out
    L.append("- Rare pathogenic variants (the chip targets common variation; BRCA, CFTR, etc. are covered at a handful of sites only).")
    L.append("- Copy-number genes such as CYP2D6, or HLA types beyond tag SNPs.")
    L.append("- Admixture percentages at the resolution of the vendor's own report (that needs their reference panel).")
    L.append("- Anything on sites the vendor did not type; the panel marks those `not on array`. PharmCAT assumes reference at untyped positions.")
    L.append("- For a male sample, X-chromosome sites are hemizygous (one copy), so carrier/homozygous language does not apply there.")
    L.append(f"- `share` uses the {POP_LABEL.get(pop, pop)} reference population unless `--ancestry` placed the sample in another 1000 Genomes group.")
    out.write_text("\n".join(L) + "\n")
    return out
