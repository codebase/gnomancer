"""Curated marker panel: genotype lookup with strand checking and interpretation.

The markers table lists, for each rsID, the two alleles on the GRCh37
forward strand (ref, effect) and optional text for 0/1/2 copies of the
effect allele. All major vendors report on the forward strand of build 37,
but a few probes are known to be flipped on some chips, so:
- if the observed alleles fit {ref, effect}: count directly
- if they fit the complement set and the SNP is not A/T or C/G: flip and note it
- otherwise: report a mismatch and do not interpret
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "I": "I", "D": "D"}
CATEGORY_ORDER = ["ancestry", "trait", "health", "carrier", "pgx"]
POPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
POP_LABEL = dict(AFR="African", AMR="Admixed American", EAS="East Asian", EUR="European", SAS="South Asian")
FREQ_FILE = Path(__file__).resolve().parent / "data" / "marker_frequencies.tsv"


def load_frequencies(path: Path = FREQ_FILE) -> dict:
    """rsid -> {pop: effect-allele frequency} from the bundled 1000 Genomes table."""
    out = {}
    if not path.exists():
        return out
    with open(path, newline="") as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        d = {p: float(r[p]) for p in POPS if r.get(p)}
        if d:
            out[r["rsid"]] = d
    return out


def genotype_frequency(effect_count: int, ploidy: int, p: float) -> float:
    """Share of a population expected to carry this genotype (Hardy-Weinberg)."""
    if ploidy == 1:
        return p if effect_count == 1 else 1 - p
    return {0: (1 - p) ** 2, 1: 2 * p * (1 - p), 2: p * p}[effect_count]


def load_markers(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    rows = [r for r in csv.DictReader(lines, delimiter="\t") if r.get("rsid")]
    if not rows:
        raise SystemExit(f"No markers loaded from {path}")
    return rows


def _lookup(df: pd.DataFrame, rsid: str):
    if rsid in df.index:
        r = df.loc[rsid]
        if isinstance(r, pd.DataFrame):
            r = r.iloc[0]
        return r
    return None


def evaluate(df: pd.DataFrame, markers: list[dict], frequencies: dict | None = None) -> pd.DataFrame:
    frequencies = load_frequencies() if frequencies is None else frequencies
    out = []
    for m in markers:
        row = _lookup(df, m["rsid"])
        rec = dict(rsid=m["rsid"], category=m["category"], gene=m["gene"], label=m["label"],
                   ref=m["ref"], effect=m["effect"], chrom="", pos="", genotype="", effect_count="",
                   interpretation="", note=m.get("note", ""), **{f"gf_{p}": None for p in POPS})
        if row is None:
            rec["genotype"] = "not on array"
            out.append(rec)
            continue
        rec["chrom"], rec["pos"] = row.chrom, int(row.pos)
        g = row.genotype
        if not isinstance(g, str) or g == "":
            rec["genotype"] = "no call"
            out.append(rec)
            continue
        rec["genotype"] = g
        alleles = set(g)
        expected = {m["ref"], m["effect"]}
        flipped = False
        if not alleles <= expected:
            comp_expected = {COMP.get(m["ref"], "?"), COMP.get(m["effect"], "?")}
            ambiguous = expected == comp_expected
            if alleles <= comp_expected and not ambiguous:
                g = "".join(COMP[a] for a in g)
                flipped = True
            else:
                rec["interpretation"] = f"allele mismatch: observed {g}, expected {m['ref']}/{m['effect']} (possible strand issue); not interpreted"
                out.append(rec)
                continue
        n = g.count(m["effect"])
        ploidy = len(g)
        rec["effect_count"] = f"{n}/{ploidy}"
        for pop, af in (frequencies.get(m["rsid"]) or {}).items():
            rec[f"gf_{pop}"] = round(genotype_frequency(n, ploidy, af), 4)
        key = {0: "interp0", 1: "interp1", 2: "interp2"}[n if ploidy == 2 else (2 if n == 1 else 0)]
        text = m.get(key) or f"{n} copy(ies) of effect allele {m['effect']}"
        if ploidy == 1:
            text += " (hemizygous)"
        if flipped:
            text += " [alleles were strand-flipped relative to the table; genotype shown as reported]"
        rec["interpretation"] = text
        out.append(rec)
    return pd.DataFrame(out)


def _gt(df, rsid):
    r = _lookup(df, rsid)
    return r.genotype if r is not None and isinstance(r.genotype, str) else None


def combos(df: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Multi-SNP readings: APOE and ABO."""
    res = []
    a, b = _gt(df, "rs429358"), _gt(df, "rs7412")
    if a and b and set(a) <= {"T", "C"} and set(b) <= {"C", "T"}:
        e4 = a.count("C")
        e2 = b.count("T")
        if e4 + e2 > 2:
            label = "unresolvable"
        else:
            alleles = ["e4"] * e4 + ["e2"] * e2
            alleles += ["e3"] * (2 - len(alleles))
            label = "/".join(sorted(alleles, key=lambda x: {"e2": 0, "e3": 1, "e4": 2}[x]))
        note = {"e3/e3": "most common genotype, reference risk for Alzheimer's and lipids",
                "e3/e4": "one e4 copy, ~2-3x population Alzheimer's risk; also higher LDL",
                "e4/e4": "two e4 copies, substantially higher Alzheimer's risk (roughly 8-12x); consider counseling",
                "e2/e3": "one e2 copy, modestly lower Alzheimer's risk; e2 relates to type III hyperlipidemia only when homozygous",
                "e2/e2": "two e2 copies; lower Alzheimer's risk, but predisposition to type III hyperlipidemia",
                "e2/e4": "one protective and one risk allele (could rarely be e1/e3); roughly e3/e4-like risk"}.get(label, "")
        res.append(("APOE", f"{label} (rs429358 {a}, rs7412 {b})", note))
    o, bb = _gt(df, "rs8176719"), _gt(df, "rs8176746")
    if o and bb:
        n_o = o.count("D")
        n_b = bb.count("T")
        if n_o == 2:
            grp = "O"
        elif n_b >= 1 and n_o == 1:
            grp = "B (B/O)"
        elif n_b == 2:
            grp = "B (B/B)"
        elif n_b == 1:
            grp = "AB"
        elif n_o == 1:
            grp = "A (A/O)"
        else:
            grp = "A (A/A)"
        res.append(("ABO blood group", grp + f" (rs8176719 {o}, rs8176746 {bb})",
                    "Rh factor is not typed on consumer arrays; A1 vs A2 subgroups not distinguished"))
    l1, l2 = _gt(df, "rs4988235"), _gt(df, "rs182549")
    if l1:
        persistent = "A" in l1 or (l2 and "T" in l2)
        res.append(("Lactase persistence", "persistent (can digest lactose as an adult)" if persistent else
                    "non-persistent (likely adult lactose intolerance)",
                    "Based on the European -13910 and -22018 variants; African and Middle Eastern persistence alleles "
                    "(rs41380347, rs145946881, rs41525747) are usually not on consumer arrays"))
    t = [_gt(df, r) for r in ("rs713598", "rs1726866", "rs10246939")]
    if all(t):
        pav = t[0].count("G") + t[1].count("G") + t[2].count("C")
        res.append(("TAS2R38 bitter taste (PTC/PROP)",
                    {6: "PAV/PAV: strong taster", 0: "AVI/AVI: non-taster"}.get(pav, "PAV/AVI or mixed: medium taster"),
                    f"haplotype score {pav}/6 across the three defining SNPs"))
    return res
