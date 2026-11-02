"""`gnomancer family KIT1 KIT2 [KIT3 ...]`: compare two or more raw files from relatives.

- Shared DNA segments per pair via apriha/lineage (HapMap genetic map), with a relationship estimate from total cM
- Which lines are shared: same mtDNA haplogroup = same maternal line; same Y = same paternal line (males)
- Side-by-side genotypes at the curated panel sites; for parent/parent/child trios, which parent each allele came from
- Mendelian-consistency check as a quality signal
Everything runs locally; nothing is uploaded.
"""
from __future__ import annotations

import argparse
import itertools
import sys
from html import escape as e
from pathlib import Path

import pandas as pd

from . import __version__
from .haplogroups import call_mt, call_y
from .load import load
from .panel import load_markers, evaluate
from .paths import DATA, HOME
from .qc import chromosome_table, summary

# Approximate total shared cM by relationship (Shared cM Project v4, rounded). Ranges overlap; this is a hint, not a verdict.
RELATIONSHIPS = [
    (2376, 3720, "parent / child", "the whole genome is shared along one chromosome of each pair"),
    (1613, 3488, "full siblings", "about a quarter of the genome is shared on both chromosomes"),
    (1201, 2282, "grandparent, half-sibling, aunt/uncle or niece/nephew", "one of these; DNA alone cannot tell them apart"),
    (396, 1397, "first cousin, great-grandparent, half-aunt/uncle or great-aunt/uncle", "one of these"),
    (102, 980, "first cousin once removed, half-first cousin or second cousin", "one of these"),
    (41, 592, "second cousin or similar", ""),
    (0, 234, "third cousin or more distant, or unrelated", ""),
]


def estimate_relationship(total_cm: float, two_chrom_cm: float) -> tuple[str, str]:
    if total_cm >= 2376 and two_chrom_cm < 200:
        return "parent / child", RELATIONSHIPS[0][3]
    if total_cm >= 1613 and two_chrom_cm >= 200:
        return "full siblings", RELATIONSHIPS[1][3]
    if total_cm < 41:
        return "unrelated (no significant shared segments)", ""
    for lo, hi, label, note in RELATIONSHIPS[2:]:
        if lo <= total_cm <= hi:
            return label, note
    return "distant relatives", ""


def mendelian_errors(child: pd.DataFrame, p1: pd.DataFrame, p2: pd.DataFrame) -> tuple[int, int]:
    """Sites where the child's genotype is impossible given both parents (ignores no-calls and haploid)."""
    common = child.index.intersection(p1.index).intersection(p2.index)
    c, a, b = child.loc[common, "genotype"], p1.loc[common, "genotype"], p2.loc[common, "genotype"]
    ok = c.notna() & a.notna() & b.notna() & (c.str.len() == 2) & (a.str.len() == 2) & (b.str.len() == 2)
    n = int(ok.sum())
    errs = 0
    for cg, ag, bg in zip(c[ok], a[ok], b[ok]):
        if not any(x in ag and y in bg or x in bg and y in ag for x in cg for y in cg if {x, y} == set(cg)):
            # child alleles must be one from each parent, in either order
            if not ((cg[0] in ag and cg[1] in bg) or (cg[1] in ag and cg[0] in bg)):
                errs += 1
    return errs, n


def shared_dna(paths: list[Path], names: list[str], outdir: Path) -> dict:
    from lineage import Lineage
    res_dir = HOME / "lineage"
    res_dir.mkdir(parents=True, exist_ok=True)
    L = Lineage(output_dir=str(outdir / "shared-dna"), resources_dir=str(res_dir), parallelize=False)
    inds = {n: L.create_individual(n, str(p)) for n, p in zip(names, paths)}
    pairs = {}
    for a, b in itertools.combinations(names, 2):
        # lineage's default of 1,100 SNPs per segment assumes a ~600k-site chip. Scale it to what the two kits share,
        # but never below 300: shorter matching runs happen by chance and get stitched into false segments.
        common = len(inds[a].snps.index.intersection(inds[b].snps.index))
        snp_threshold = max(300, min(1100, int(1100 * common / 600_000)))
        r = L.find_shared_dna([inds[a], inds[b]], cM_threshold=7.0, snp_threshold=snp_threshold, save_output=True)
        one, two = r.get("one_chrom_shared_dna"), r.get("two_chrom_shared_dna")
        one_cm = float(one["cMs"].sum()) if one is not None and len(one) else 0.0
        two_cm = float(two["cMs"].sum()) if two is not None and len(two) else 0.0
        segs = int(len(one)) if one is not None else 0
        label, note = estimate_relationship(one_cm, two_cm)
        pairs[(a, b)] = dict(one_cm=round(one_cm, 1), two_cm=round(two_cm, 1), segments=segs, relationship=label, note=note,
                             largest_cm=round(float(one["cMs"].max()), 1) if one is not None and len(one) else 0.0)
    return pairs


def run(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gnomancer family", description="Compare raw files from relatives: shared DNA, shared lineages, trait inheritance.")
    ap.add_argument("kits", nargs="+", help="two or more raw-data files")
    ap.add_argument("--names", default=None, help="comma-separated labels in the same order (default: file stems)")
    ap.add_argument("--out", default=None, help="output folder (default: family_analysis next to the first file)")
    ap.add_argument("--no-haplogroups", action="store_true", help="skip Y/mtDNA calls (faster)")
    a = ap.parse_args(argv)
    paths = [Path(k).expanduser() for k in a.kits]
    if len(paths) < 2:
        print("need at least two kits", file=sys.stderr)
        return 2
    names = a.names.split(",") if a.names else [p.stem.split("_")[0].split(".")[0] for p in paths]
    out = Path(a.out).expanduser() if a.out else paths[0].with_name("family_analysis")
    out.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] loading {len(paths)} kits ...", file=sys.stderr)
    loaded = {n: load(str(p)) for n, p in zip(names, paths)}
    dfs = {n: L.df for n, L in loaded.items()}
    sexes = {n: summary(dfs[n], chromosome_table(dfs[n]), loaded[n].sex)["sex_inferred"] for n in names}

    print("[2/4] shared DNA (lineage) ...", file=sys.stderr)
    pairs = shared_dna(paths, names, out)

    hg = {}
    if not a.no_haplogroups:
        print("[3/4] lineages ...", file=sys.stderr)
        for n in names:
            d = out / "haplogroups" / n
            d.mkdir(parents=True, exist_ok=True)
            y = call_y(dfs[n], d) if sexes[n] == "Male" else dict(error="not male")
            mt = call_mt(dfs[n], d, loaded[n].build_original)
            hg[n] = dict(y=y.get("haplogroup") if not y.get("error") else None, mt=mt.get("haplogroup") if not mt.get("error") else None)

    print("[4/4] panel and inheritance ...", file=sys.stderr)
    markers = load_markers(DATA / "markers.tsv")
    panels = {n: evaluate(dfs[n], markers).set_index("rsid") for n in names}
    # trios: child shares parent/child with two people who are unrelated to each other
    def rel(x, y):
        return (pairs.get((x, y)) or pairs.get((y, x)) or {}).get("relationship", "")
    trios = []
    for c in names:
        parents = [o for o in names if o != c and rel(c, o) == "parent / child"]
        for p1, p2 in itertools.combinations(parents, 2):
            if rel(p1, p2).startswith("unrelated"):
                trios.append((c, p1, p2))
    mend = {t: mendelian_errors(dfs[t[0]], dfs[t[1]], dfs[t[2]]) for t in trios}

    write_reports(out, names, sexes, pairs, hg, panels, markers, trios, mend)
    print(f"\nFamily report: {out / 'family-report.html'}", file=sys.stderr)
    return 0


def _phase(child_gt: str, p1_gt: str, p2_gt: str) -> str:
    """Explain which parent each child allele came from when it is determinable."""
    if any(not isinstance(g, str) or len(g) != 2 for g in (child_gt, p1_gt, p2_gt)):
        return ""
    opts = []
    for x, y in ((child_gt[0], child_gt[1]), (child_gt[1], child_gt[0])):
        if x in p1_gt and y in p2_gt:
            opts.append((x, y))
    opts = sorted(set(opts))
    if not opts:
        return "inconsistent with both parents"
    if len(opts) == 1:
        x, y = opts[0]
        return f"{x} from first parent, {y} from second"
    return "either way"


def write_reports(out, names, sexes, pairs, hg, panels, markers, trios, mend):
    md = [f"# Family comparison\n\nGenerated by gnomancer v{__version__}. Nothing was uploaded.\n", "## Who shares what\n",
          "| Pair | Shared (cM) | Fully identical (cM) | Segments | Largest (cM) | Likely relationship |", "|---|---|---|---|---|---|"]
    for (a, b), r in pairs.items():
        md.append(f"| {a} & {b} | {r['one_cm']:.0f} | {r['two_cm']:.0f} | {r['segments']} | {r['largest_cm']:.0f} | **{r['relationship']}**" + (f" ({r['note']})" if r['note'] else "") + " |")
    md.append("\nTotal shared centimorgans estimate relationships only roughly; overlapping ranges are listed as such. Plots of the shared segments per pair are in `shared-dna/`.\n")
    if hg:
        md += ["## Lineages\n", "| Person | Sex | Y haplogroup | mtDNA haplogroup |", "|---|---|---|---|"]
        md += [f"| {n} | {sexes[n]} | {hg[n]['y'] or '—'} | {hg[n]['mt'] or '—'} |" for n in names]
        lines = []
        for a, b in itertools.combinations(names, 2):
            if hg[a]["mt"] and hg[a]["mt"] == hg[b]["mt"]:
                lines.append(f"- {a} and {b} share the maternal line (mtDNA {hg[a]['mt']}).")
            if hg[a]["y"] and hg[a]["y"] == hg[b]["y"]:
                lines.append(f"- {a} and {b} share the paternal line (Y {hg[a]['y']}).")
        md += ["", *(lines or ["- No shared lineage labels among these kits (or too few kits were male for a Y comparison)."]), ""]
    if trios:
        md.append("## Parent-parent-child checks\n")
        for t in trios:
            errs, n = mend[t]
            md.append(f"- {t[0]} with parents {t[1]} and {t[2]}: {errs:,} Mendelian inconsistencies in {n:,} sites ({errs / max(n, 1):.2%}). Below 1% is normal array error; well above that means a mislabeled kit.")
        md.append("")
    md.append("## Panel sites side by side\n")
    md.append("| Site | " + " | ".join(names) + (" | Inheritance |" if trios else " |"))
    md.append("|---|" + "---|" * len(names) + ("---|" if trios else ""))
    rows_html = []
    for m in markers:
        gts = [panels[n].loc[m["rsid"], "genotype"] if m["rsid"] in panels[n].index else "" for n in names]
        if all(g in ("", "not on array", "no call") for g in gts):
            continue
        inh = ""
        if trios:
            c, p1, p2 = trios[0]
            inh = _phase(panels[c].loc[m["rsid"], "genotype"], panels[p1].loc[m["rsid"], "genotype"], panels[p2].loc[m["rsid"], "genotype"]) if m["rsid"] in panels[c].index else ""
            inh = inh.replace("first parent", p1).replace("second", p2)
        md.append(f"| {m['gene']} {m['label']} | " + " | ".join(gts) + (f" | {inh} |" if trios else " |"))
        rows_html.append((f"{m['gene']} · {m['label']}", gts, inh))
    (out / "family-report.md").write_text("\n".join(md) + "\n")

    css = ("body{font:15px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',Inter,Roboto,Helvetica,Arial,sans-serif;color:#111;background:#fbfbf9;margin:0}"
           "@media(prefers-color-scheme:dark){body{background:#141413;color:#f0efec}table,th,td{background:#1a1a19;color:#f0efec}th,td{border-color:#2e2e2c}table{border-color:#2e2e2c}.card{background:#1a1a19;border-color:#2e2e2c;color:#f0efec}.muted{color:#9a9890}img{border-color:#2e2e2c}}"
           ".wrap{max-width:1040px;margin:0 auto;padding:24px 16px 64px}h1{font-size:28px;margin:0 0 4px}h2{font-size:21px;margin:32px 0 8px}"
           "table{width:100%;border-collapse:collapse;font-size:14px;background:#fff;border:1px solid #e6e4df;border-radius:12px;overflow:hidden}th,td{text-align:left;padding:8px 10px;border-top:1px solid #e6e4df}th{border-top:0;font-weight:600}"
           ".gt{font-family:ui-monospace,Menlo,monospace}.card{background:#fff;border:1px solid #e6e4df;border-radius:12px;padding:16px 18px;margin:12px 0}img{max-width:100%;border-radius:12px;border:1px solid #e6e4df;margin:8px 0}.muted{color:#777;font-size:13px}")
    H = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Family comparison · gnomancer</title><style>{css}</style></head><body><div class="wrap">',
         f'<h1>Family comparison</h1><p class="muted">gnomancer v{__version__} · {len(names)} kits · nothing was uploaded</p>', "<h2>Who shares what</h2>",
         "<table><thead><tr><th>Pair</th><th>Shared (cM)</th><th>Fully identical (cM)</th><th>Segments</th><th>Largest (cM)</th><th>Likely relationship</th></tr></thead><tbody>"]
    for (a, b), r in pairs.items():
        H.append(f"<tr><td>{e(a)} &amp; {e(b)}</td><td>{r['one_cm']:.0f}</td><td>{r['two_cm']:.0f}</td><td>{r['segments']}</td><td>{r['largest_cm']:.0f}</td><td><b>{e(r['relationship'])}</b>" + (f"<div class='muted'>{e(r['note'])}</div>" if r['note'] else "") + "</td></tr>")
    H.append("</tbody></table><p class='muted'>Total shared centimorgans estimate relationships only roughly; overlapping ranges are listed as such.</p>")
    for (a, b), r in pairs.items():
        if r["segments"] == 0:
            continue
        pngs = sorted((out / "shared-dna").glob(f"shared_dna_{a}_{b}_*.png"), key=lambda p: p.stat().st_mtime)
        if pngs:
            H.append(f'<div class="card"><b>{e(a)} &amp; {e(b)}</b> · {r["one_cm"]:.0f} cM shared<br><img src="shared-dna/{e(pngs[-1].name)}" alt="shared DNA segments between {e(a)} and {e(b)}"></div>')
    if hg:
        H.append("<h2>Lineages</h2><table><thead><tr><th>Person</th><th>Sex</th><th>Y haplogroup</th><th>mtDNA haplogroup</th></tr></thead><tbody>")
        H += [f"<tr><td>{e(n)}</td><td>{e(sexes[n])}</td><td class='gt'>{e(hg[n]['y'] or '—')}</td><td class='gt'>{e(hg[n]['mt'] or '—')}</td></tr>" for n in names]
        H.append("</tbody></table><ul>")
        any_line = False
        for a, b in itertools.combinations(names, 2):
            if hg[a]["mt"] and hg[a]["mt"] == hg[b]["mt"]:
                H.append(f"<li>{e(a)} and {e(b)} share the maternal line (mtDNA {e(hg[a]['mt'])}).</li>"); any_line = True
            if hg[a]["y"] and hg[a]["y"] == hg[b]["y"]:
                H.append(f"<li>{e(a)} and {e(b)} share the paternal line (Y {e(hg[a]['y'])}).</li>"); any_line = True
        if not any_line:
            H.append("<li>No shared lineage labels among these kits.</li>")
        H.append("</ul>")
    if trios:
        H.append("<h2>Parent-parent-child checks</h2><ul>")
        for t in trios:
            errs, n = mend[t]
            H.append(f"<li>{e(t[0])} with parents {e(t[1])} and {e(t[2])}: {errs:,} Mendelian inconsistencies in {n:,} sites ({errs / max(n, 1):.2%}). Below 1% is normal array error.</li>")
        H.append("</ul>")
    H.append("<h2>Panel sites side by side</h2><table><thead><tr><th>Site</th>" + "".join(f"<th>{e(n)}</th>" for n in names) + ("<th>Inheritance</th>" if trios else "") + "</tr></thead><tbody>")
    for label, gts, inh in rows_html:
        H.append(f"<tr><td>{e(label)}</td>" + "".join(f"<td class='gt'>{e(g)}</td>" for g in gts) + (f"<td class='muted'>{e(inh)}</td>" if trios else "") + "</tr>")
    H.append("</tbody></table><p class='muted'>Genotypes on the GRCh37 forward strand. This is not medical advice; see each person's own report for interpretation.</p></div></body></html>")
    (out / "family-report.html").write_text("".join(H))
