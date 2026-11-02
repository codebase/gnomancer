"""Markdown report writer."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .panel import CATEGORY_ORDER, POPS, POP_LABEL

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
    L = []
    L.append(f"# Raw DNA analysis report\n\n{DISCLAIMER}\n")
    if extra.get("narrative"):
        L.append("## In plain language\n\n" + extra["narrative"] + "\n")
    g = extra.get("glance") or {}
    if g:
        L.append("## At a glance\n")
        for title, key in (("Who the data says you are", "who"), ("Worth knowing", "notable"), ("Reassuring", "reassuring"), ("Sensible next steps", "next")):
            if g.get(key):
                L.append(f"**{title}**\n" + "\n".join(f"- {b}" for b in g[key]) + "\n")
    L.append("## 1. Input\n")
    L.append(f"- Vendor/format detected: **{loaded.source}**" + (f" (chip {loaded.chip} {loaded.chip_version})" if loaded.chip else ""))
    L.append(f"- Genome build: GRCh{loaded.build}" + (f" (remapped from build {loaded.build_original})" if loaded.remapped else "") +
             ("" if loaded.build_detected else " (assumed, not detected)"))
    L.append(f"- SNPs loaded: {loaded.count:,}")
    for n in loaded.notes:
        L.append(f"- Note: {n}")
    L.append("\n## 2. Quality\n")
    L.append(f"| Metric | Value |\n|---|---|")
    L.append(f"| No-calls | {qc['no_call']:,} ({qc['no_call_rate']:.2%}) |")
    L.append(f"| Indel calls (I/D) | {qc['indel_calls']:,} |")
    L.append(f"| Autosomal heterozygosity | {qc['autosomal_het']:.3f} |")
    L.append(f"| X heterozygosity | {qc['x_het_rate']:.4f} |")
    L.append(f"| Y SNPs called | {qc['y_called']:,} / {qc['y_total']:,} |")
    L.append(f"| Sex inferred (X/Y) | {qc['sex_inferred']} (library: {qc['sex_library'] or 'n/a'}) |")
    for f in qc["flags"]:
        L.append(f"\n⚠ {f}")
    L.append("\nPer-chromosome table: `chromosome-qc.tsv`.\n")
    L.append("## 3. Runs of homozygosity\n")
    L.append(roh_text)
    if not roh.empty:
        L.append("\n" + _md_table(roh.head(15), ["chrom", "start", "end", "length_mb", "snps"]))
    L.append("\n## 4. Paternal line (Y chromosome)\n")
    if y.get("haplogroup") and not y.get("error"):
        L.append(f"**Haplogroup {y['haplogroup']}** (ISOGG long form {y['ycc']}), called by yhaplo from {y['n_y_snps']:,} Y SNPs.")
        L.append(f"\nDerived-SNP path: `{y['path']}`" if y.get("path") else "")
        yl = y.get("yleaf") or {}
        if yl.get("trees"):
            L.append("\nSecond opinion from Yleaf on current trees (" + f"{yl['n_y_snps']:,} Y SNPs):")
            names = dict(yfull="YFull v14", ftdna="FamilyTreeDNA", isogg="ISOGG")
            for tree, t in yl["trees"].items():
                L.append(f"- {names.get(tree, tree)}: **{t['haplogroup']}** (defining marker {t['marker']}, {t['valid_markers']} valid markers, QC {t['qc']})")
        elif yl.get("error"):
            L.append(f"\nYleaf second opinion not available: {yl['error']}")
        L.append("\nArrays type only a few thousand Y sites, so this is the deepest branch the chip can see; "
                 "Y sequencing (e.g. Big Y-700) is needed for a terminal subclade.")
    else:
        L.append(f"Not called: {y.get('error', 'unknown')}")
    L.append("\n## 5. Maternal line (mitochondrial DNA)\n")
    if mt.get("haplogroup") and not mt.get("error"):
        L.append(f"**Haplogroup {mt['haplogroup']}** (Haplogrep3 quality {mt['quality']:.2f}) from {mt['n_mt_snps']} typed MT positions.")
        L.append(f"\nDifferences from rCRS at typed sites: {', '.join(mt['diffs']) or 'none'}.")
        if mt["indels"]:
            L.append(f"Indel probes ignored: {', '.join(mt['indels'])}.")
        L.append("\nOnly a few hundred MT sites are typed, so sub-branches usually cannot be resolved; full mtDNA sequencing would.")
    else:
        L.append(f"Not called: {mt.get('error', 'unknown')}")
    L.append("\n## 6. Ancestry estimate (ezancestry)\n")
    if not ancestry:
        failed = [n for n in loaded.notes if n.startswith("ezancestry failed")]
        if failed:
            L.append(f"Requested but failed: {failed[0]}")
        else:
            L.append("Not run. Add `--ancestry` for a 1000 Genomes population classification "
                     "(needs `pip install git+https://github.com/arvkevi/ezancestry`, Python ≤3.11).")
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
    L.append("\n## 7. Multi-SNP readings\n")
    for name, val, note in combos:
        L.append(f"- **{name}:** {val}" + (f" — {note}" if note else ""))
    L.append("\n## 8. Marker panel\n")
    pop = extra.get("pop") or "EUR"
    L.append(f"Effect-allele counts are on the GRCh37 forward strand. `not on array` means the vendor's chip does not type that site. "
             f"`share` is the fraction of the {POP_LABEL.get(pop, pop)} 1000 Genomes reference population expected to carry the same genotype "
             f"(Hardy-Weinberg from the effect-allele frequency); the TSV has all five super-populations.\n")
    titles = dict(ancestry="Ancestry-informative markers", trait="Traits", health="Health-associated common variants",
                  carrier="Carrier / rare pathogenic sites typed on arrays", pgx="Pharmacogenomics")
    for cat in CATEGORY_ORDER:
        sub = panel[panel.category == cat]
        if sub.empty:
            continue
        L.append(f"\n### {titles.get(cat, cat)}\n")
        sub = sub.copy()
        col = f"gf_{pop}"
        sub["share"] = sub[col].map(lambda v: f"{v:.0%}" if isinstance(v, float) else "") if col in sub else ""
        L.append(_md_table(sub, ["rsid", "gene", "label", "genotype", "effect_count", "share", "interpretation"]))
    cv = extra.get("clinvar") or {}
    L.append("\n## 9. ClinVar scan of every typed site\n")
    L.append(extra.get("clinvar_text", "not run"))
    hits = cv.get("hits")
    if hits is not None and not hits.empty:
        L.append("\n" + _md_table(hits, ["rsid", "gene", "genotype", "zygosity", "significance", "stars", "confidence", "condition", "flags"]))
    pc = extra.get("pharmcat") or {}
    L.append("\n## 10. PharmCAT pharmacogenomics\n")
    L.append(extra.get("pharmcat_text", "not run"))
    calls = pc.get("calls")
    if calls is not None and not calls.empty:
        from .summary import PGX_DRUGS
        c = calls.copy()
        c["phenotype"] = c.phenotype.str.strip().str.rstrip("/").str.strip().replace("", "—")
        c["drugs"] = c.gene.map(PGX_DRUGS).fillna("")
        c["coverage"] = c.missing_positions.map({True: "partial", False: "complete"})
        L.append("\n" + _md_table(c, ["gene", "diplotype", "phenotype", "drugs", "typed", "coverage"]))
        L.append("\nRYR1/CACNA1S 'Uncertain Susceptibility' means a variant of uncertain significance for malignant hyperthermia was seen, not a positive result. "
                 "The CFTR row describes ivacaftor response, not cystic fibrosis status.")
    L.append("\n## 11. What this file cannot tell you\n")
    L.append("- Rare pathogenic variants (the chip targets common variation; BRCA, CFTR, etc. are covered at a handful of sites only).")
    L.append("- Copy-number genes such as CYP2D6, or HLA types beyond tag SNPs.")
    L.append("- Admixture percentages at the resolution of the vendor's own report (that needs their reference panel).")
    L.append("- Anything on sites the vendor did not type; the panel marks those `not on array`.")
    out.write_text("\n".join(L) + "\n")
    return out
