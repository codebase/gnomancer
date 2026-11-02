"""One-page pharmacogenomics card for a pharmacist or prescriber.

Writes `pgx-card.html` (print-ready, one page) and, when reportlab is installed, `pgx-card.pdf`.
Content: the PharmCAT diplotypes and phenotypes with the drugs each gene matters for, plus the clotting and
abacavir single-site results, and a footer that says plainly this is array-based and provisional.
"""
from __future__ import annotations

import datetime as dt
from html import escape as e
from pathlib import Path

import pandas as pd

from . import __version__
from .summary import PGX_DRUGS

FOOTER = ("Derived from a consumer genotyping array, not a clinical test. PharmCAT assumed the reference allele at positions "
          "the chip did not type, so rare alleles can be missed and every call is provisional. Confirm with a clinical "
          "pharmacogenomic panel before changing any prescription. Guidelines: cpicpgx.org and pharmgkb.org.")


def _tidy_result(d: str, limit: int = 3) -> str:
    """'*1/*5 OR *1/*15 OR *1/*40 OR ...' -> first few alternatives plus a count."""
    alts = [a.strip() for a in str(d).split(" OR ")]
    if len(alts) <= limit:
        return " or ".join(alts)
    return " or ".join(alts[:limit]) + f" (+{len(alts) - limit} more possible)"


def _tidy_pheno(p: str) -> str:
    """'Decreased Function / Decreased Function / Poor Function' -> 'Decreased Function or Poor Function'."""
    seen, out = set(), []
    for x in [x.strip() for x in str(p).split("/")]:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return " or ".join(out) if out else "—"


def _rows(pc: dict, panel: pd.DataFrame) -> list[dict]:
    rows = []
    calls = pc.get("calls")
    if calls is not None and not calls.empty:
        for _, r in calls.iterrows():
            pheno = _tidy_pheno(r.phenotype)
            rows.append(dict(gene=r.gene, result=_tidy_result(r.diplotype), phenotype=pheno, drugs=PGX_DRUGS.get(r.gene, ""),
                             coverage="partial" if r.missing_positions else "complete"))
    have = {r["gene"] for r in rows}
    extra = {"HCP5": ("HLA-B*57:01 tag", "abacavir"), "F5": ("Factor V Leiden", PGX_DRUGS["F5"]), "F2": ("Prothrombin G20210A", PGX_DRUGS["F2"]),
             "VKORC1": ("warfarin sensitivity", "warfarin")}
    if panel is not None and not panel.empty:
        for _, r in panel[panel.category.isin(["pgx", "health"])].iterrows():
            if r.gene in extra and r.gene not in have and r.genotype not in ("not on array", "no call"):
                label, drugs = extra[r.gene]
                rows.append(dict(gene=r.gene, result=f"{label}: {r.genotype}", phenotype=str(r.interpretation), drugs=drugs, coverage="single site"))
                have.add(r.gene)
    return rows


def write_html(out: Path, rows: list[dict], source: str, input_name: str) -> Path:
    today = dt.date.today().isoformat()
    css = ("@page{size:A4;margin:14mm}body{font:11px/1.35 -apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;color:#111;margin:0}"
           ".wrap{max-width:720px;margin:0 auto;padding:16px}h1{font-size:18px;margin:0 0 2px}.sub{color:#555;font-size:11px;margin:0 0 10px}"
           ".id{display:flex;gap:18px;margin:8px 0 12px;font-size:11px}.id span{border-bottom:1px solid #999;min-width:180px;display:inline-block;padding:0 4px}"
           "table{width:100%;border-collapse:collapse;font-size:10.5px}th,td{text-align:left;padding:4px 6px;border-bottom:1px solid #ddd;vertical-align:top}th{background:#f3f1ee;font-weight:600}"
           "td.g{font-weight:600;white-space:nowrap}td.r{font-family:ui-monospace,Menlo,monospace;font-size:10px}.cov{color:#777;font-size:9.5px}"
           ".foot{margin-top:12px;font-size:9.5px;color:#444;border-top:1px solid #bbb;padding-top:6px}@media print{.wrap{padding:0}}")
    H = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Pharmacogenomics card</title><style>{css}</style></head><body><div class="wrap">',
         f'<h1>Pharmacogenomics summary</h1><p class="sub">Generated {today} by gnomancer v{__version__} from a {e(source)} raw-data file ({e(input_name)}). For a pharmacist or prescriber.</p>',
         '<div class="id"><div>Name: <span>&nbsp;</span></div><div>Date of birth: <span>&nbsp;</span></div></div>',
         '<table><thead><tr><th>Gene</th><th>Result</th><th>Phenotype</th><th>Drugs it matters for</th></tr></thead><tbody>']
    for r in rows:
        H.append(f'<tr><td class="g">{e(r["gene"])}</td><td class="r">{e(r["result"])}<div class="cov">{e(r["coverage"])}</div></td><td>{e(r["phenotype"])}</td><td>{e(r["drugs"])}</td></tr>')
    H.append("</tbody></table>")
    H.append(f'<p class="foot">{e(FOOTER)} "partial" means untyped positions were assumed reference.</p></div></body></html>')
    out.write_text("".join(H))
    return out


def write_pdf(out: Path, rows: list[dict], source: str, input_name: str) -> Path | None:
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
    except ImportError:
        return None
    styles = getSampleStyleSheet()
    small = ParagraphStyle("small", parent=styles["BodyText"], fontSize=8, leading=10)
    cell = ParagraphStyle("cell", parent=styles["BodyText"], fontSize=8.5, leading=10.5)
    doc = SimpleDocTemplate(str(out), pagesize=A4, leftMargin=14 * mm, rightMargin=14 * mm, topMargin=14 * mm, bottomMargin=12 * mm,
                            title="Pharmacogenomics summary", author="gnomancer")
    story = [Paragraph("Pharmacogenomics summary", styles["Title"]),
             Paragraph(f"Generated {dt.date.today().isoformat()} by gnomancer v{__version__} from a {source} raw-data file ({input_name}). For a pharmacist or prescriber.", small),
             Spacer(1, 4), Paragraph("Name: ____________________________&nbsp;&nbsp;&nbsp;&nbsp;Date of birth: ______________", styles["BodyText"]), Spacer(1, 8)]
    data = [["Gene", "Result", "Phenotype", "Drugs it matters for"]]
    for r in rows:
        data.append([Paragraph(f"<b>{r['gene']}</b>", cell), Paragraph(f"{r['result']}<br/><font size=7 color='#777777'>{r['coverage']}</font>", cell),
                     Paragraph(r["phenotype"], cell), Paragraph(r["drugs"], cell)])
    t = Table(data, colWidths=[22 * mm, 52 * mm, 42 * mm, 66 * mm], repeatRows=1)
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f3f1ee")), ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                           ("FONTSIZE", (0, 0), (-1, 0), 9), ("LINEBELOW", (0, 0), (-1, -1), 0.3, colors.HexColor("#dddddd")),
                           ("VALIGN", (0, 0), (-1, -1), "TOP"), ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
    story += [t, Spacer(1, 10), Paragraph(FOOTER + ' "partial" means untyped positions were assumed reference.', small)]
    doc.build(story)
    return out


def write(outdir: Path, pc: dict, panel: pd.DataFrame, source: str, input_name: str) -> dict:
    rows = _rows(pc, panel)
    if not rows:
        return dict(html=None, pdf=None, rows=0)
    html = write_html(outdir / "pgx-card.html", rows, source, input_name)
    pdf = write_pdf(outdir / "pgx-card.pdf", rows, source, input_name)
    return dict(html=str(html), pdf=str(pdf) if pdf else None, rows=len(rows))
