"""Self-contained HTML report: one file, inline CSS/JS, light and dark, printable, no external requests."""
from __future__ import annotations

import datetime as dt
import json
from html import escape as e
from pathlib import Path

import pandas as pd

from . import __version__
from .context import lookup
from .panel import CATEGORY_ORDER
from .report import DISCLAIMER

CAT_TITLES = dict(ancestry="Ancestry-informative markers", trait="Traits",
                  health="Health-associated common variants", carrier="Carrier and rare pathogenic sites",
                  pgx="Pharmacogenomics")
CAT_BLURB = dict(
    ancestry="Sites whose allele frequencies differ strongly between world regions. Together with the Y and mtDNA lineages they hint at ancestry; they are not an admixture estimate.",
    trait="Common variants with well-replicated effects on visible or sensory traits. Most are probabilistic, not deterministic.",
    health="Common variants that shift the odds of a condition. Odds ratios are from population studies; absolute risk depends on far more than one site.",
    carrier="A handful of well-known pathogenic variants that happen to be on consumer chips. A negative here covers only this one site, not the gene.",
    pgx="Variants that change how drugs are metabolized. Useful to mention to a prescriber; clinical decisions need a clinical-grade test.")
SENSITIVE = {"health", "carrier"}

CSS = """
:root{color-scheme:light;--bg:#fcfcfb;--surface:#ffffff;--line:#e6e4df;--ink:#0b0b0b;--ink-2:#52514e;--ink-3:#8a8781;
--accent:#7c3aed;--accent-soft:#efe9fd;--accent-ink:#5b21b6;--good:#0ca30c;--warn:#fab219;--serious:#ec835a;--critical:#d03b3b;
--radius:12px;--mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Inter,Roboto,Helvetica,Arial,sans-serif}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){color-scheme:dark;--bg:#141413;--surface:#1a1a19;--line:#2e2e2c;--ink:#f0efec;--ink-2:#c3c2b7;--ink-3:#8f8e86;--accent:#a78bfa;--accent-soft:#2e1a5e;--accent-ink:#c4b5fd}}
:root[data-theme=dark]{color-scheme:dark;--bg:#141413;--surface:#1a1a19;--line:#2e2e2c;--ink:#f0efec;--ink-2:#c3c2b7;--ink-3:#8f8e86;--accent:#a78bfa;--accent-soft:#2e1a5e;--accent-ink:#c4b5fd}
*{box-sizing:border-box}html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 var(--sans)}
.wrap{max-width:1040px;margin:0 auto;padding:0 16px 64px}
header{padding:40px 0 20px;display:flex;flex-wrap:wrap;gap:12px 24px;align-items:baseline;justify-content:space-between}
header h1{font-size:28px;margin:0;letter-spacing:-.01em}header h1 span{color:var(--accent)}
header .meta{color:var(--ink-2);font-size:14px}
.toolbar{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:4px 0 20px}
button,.btn{font:inherit;font-size:14px;padding:6px 12px;border-radius:8px;border:1px solid var(--line);background:var(--surface);color:var(--ink);cursor:pointer}
button:hover{border-color:var(--accent)}
nav.toc{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:14px;margin-bottom:24px}nav.toc a{color:var(--accent-ink);text-decoration:none}nav.toc a:hover{text-decoration:underline}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:0 0 28px}
.tile .value span.long{font-size:17px}.tile{background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:14px 16px;min-width:0}
.tile .label{font-size:12.5px;color:var(--ink-2);margin-bottom:4px}.tile .value{font-size:21px;font-weight:600;letter-spacing:-.01em;overflow-wrap:anywhere}.tile .value.long{font-size:17px}
.tile .sub{font-size:12.5px;color:var(--ink-3);margin-top:2px}
section{margin:36px 0}section h2{font-size:21px;margin:0 0 6px;letter-spacing:-.01em}section>p.lead{color:var(--ink-2);margin:0 0 14px;max-width:72ch}
.card{background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:18px 20px;margin:12px 0}
.card h3{margin:0 0 6px;font-size:17px}.card .hg{font-size:26px;font-weight:600;margin:2px 0 6px}.card .hg small{font-size:14px;font-weight:400;color:var(--ink-2);margin-left:8px}
.lineage{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}
.chain{margin:10px 0 0;padding:0;list-style:none}.chain li{padding:8px 0;border-top:1px solid var(--line);font-size:14px}.chain li b{display:block;color:var(--ink);margin-bottom:2px}
.chain li{color:var(--ink-2)}
.flag{display:flex;gap:8px;align-items:flex-start;padding:10px 12px;border-radius:8px;border:1px solid var(--line);margin:8px 0;font-size:14px}
.flag .ico{font-weight:700;flex:none;width:20px;height:20px;border-radius:50%;display:inline-grid;place-items:center;font-size:12px;color:#fff}
.flag.warn .ico{background:var(--warn);color:#000}.flag.good .ico{background:var(--good)}.flag.critical .ico{background:var(--critical)}
table{width:100%;border-collapse:collapse;font-size:14px;background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);overflow:hidden}
th,td{text-align:left;padding:8px 10px;border-top:1px solid var(--line);vertical-align:top}th{border-top:0;font-weight:600;color:var(--ink-2);font-size:13px;background:var(--surface)}
tr:first-child td{border-top:0}
.tablewrap{overflow-x:auto;border-radius:var(--radius)}
td.gt,td.num{font-family:var(--mono);font-variant-numeric:tabular-nums;white-space:nowrap}
.badge{display:inline-block;min-width:34px;text-align:center;padding:1px 8px;border-radius:999px;font-size:12.5px;font-family:var(--mono);border:1px solid var(--line);color:var(--ink-2)}
.badge.n1{border-color:var(--accent);color:var(--accent-ink)}.badge.n2{background:var(--accent);border-color:var(--accent);color:#fff}
.bar{fill:var(--accent)}
.badge.na{color:var(--ink-3);border-style:dashed}
.readings{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px}
.readings .card .value{font-size:18px;font-weight:600;margin:2px 0 4px}.readings .card .note{font-size:13.5px;color:var(--ink-2)}
details.sensitive{border:1px dashed var(--line);border-radius:var(--radius);padding:4px 16px;margin:12px 0}
details.sensitive>summary{cursor:pointer;padding:10px 0;font-weight:600}details.sensitive>summary small{display:block;font-weight:400;color:var(--ink-2);font-size:13.5px;margin-top:2px}
.filter{display:flex;gap:8px;align-items:center;margin:0 0 12px}.filter input{font:inherit;font-size:14px;padding:7px 10px;border:1px solid var(--line);border-radius:8px;background:var(--surface);color:var(--ink);width:min(100%,360px)}
.chart{background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:14px 12px 6px;margin:12px 0;position:relative}
.chart svg{display:block;width:100%;height:auto}
.chart .tt{position:absolute;pointer-events:none;background:var(--ink);color:var(--bg);font-size:12.5px;padding:6px 8px;border-radius:6px;white-space:nowrap;display:none;transform:translate(-50%,-110%)}
.bar{fill:var(--accent)}.bar.sex{fill:var(--accent-soft);stroke:var(--accent);stroke-width:1}
.axis{stroke:var(--line)}.lbl{fill:var(--ink-2);font-size:11px;font-family:var(--sans)}
.caption{font-size:13px;color:var(--ink-3);margin:6px 2px 0}
.limits li{margin:4px 0}
.glance{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}.glance .card ul{margin:6px 0 0;padding-left:18px}.glance .card li{margin:5px 0;font-size:14.5px}
.narrative .prose p{margin:0 0 10px;font-size:15.5px;line-height:1.6;max-width:72ch}
footer{margin-top:48px;padding-top:16px;border-top:1px solid var(--line);color:var(--ink-2);font-size:13.5px}
.disclaimer{background:var(--surface);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:8px;padding:12px 14px;font-size:14px;color:var(--ink-2);margin:0 0 20px}
.kv{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px;font-size:14px}.kv dt{color:var(--ink-2)}.kv dd{margin:0}
.hidden{display:none!important}
@media print{body{font-size:12px}.toolbar,nav.toc,.filter,.chart .tt{display:none}details.sensitive{border:0;padding:0}details.sensitive[open]>summary{display:none}section{break-inside:avoid-page}table{break-inside:auto}tr{break-inside:avoid}}
"""

JS = """
(function(){
  var root=document.documentElement;
  try{var t=localStorage.getItem('gnomancer-theme');if(t)root.setAttribute('data-theme',t);}catch(e){}
  var tb=document.getElementById('theme');if(tb)tb.onclick=function(){var cur=root.getAttribute('data-theme');var dark=cur?cur==='dark':matchMedia('(prefers-color-scheme: dark)').matches;var nxt=dark?'light':'dark';root.setAttribute('data-theme',nxt);try{localStorage.setItem('gnomancer-theme',nxt);}catch(e){}};
  var pb=document.getElementById('print');if(pb)pb.onclick=function(){window.print();};
  var rb=document.getElementById('reveal');if(rb)rb.onclick=function(){document.querySelectorAll('details.sensitive').forEach(function(d){d.open=true;});};
  var q=document.getElementById('q');if(q)q.oninput=function(){var v=q.value.trim().toLowerCase();document.querySelectorAll('table.panel tbody tr').forEach(function(tr){tr.classList.toggle('hidden',v&&tr.textContent.toLowerCase().indexOf(v)<0);});document.querySelectorAll('table.panel').forEach(function(t){var any=t.querySelectorAll('tbody tr:not(.hidden)').length>0;t.closest('.catblock').classList.toggle('hidden',!any);});};
  document.querySelectorAll('.chart').forEach(function(c){var tt=c.querySelector('.tt');c.querySelectorAll('[data-tt]').forEach(function(el){el.addEventListener('mousemove',function(ev){tt.textContent=el.getAttribute('data-tt');tt.style.display='block';var r=c.getBoundingClientRect();tt.style.left=(ev.clientX-r.left)+'px';tt.style.top=(ev.clientY-r.top-8)+'px';});el.addEventListener('mouseleave',function(){tt.style.display='none';});});});
})();
"""


def _split_detail(val: str):
    """'A (A/O) (rs8176719 DI, rs8176746 GG)' -> ('A (A/O)', 'rs8176719 DI, rs8176746 GG')"""
    i = val.find(" (rs")
    if i < 0:
        return val, ""
    return val[:i], val[i + 2:].rstrip(")")


def _tile(label, value, sub=""):
    return f'<div class="tile"><div class="label">{e(label)}</div><div class="value">{value}</div>' + (f'<div class="sub">{e(sub)}</div>' if sub else "") + "</div>"


def _het_chart(ct: pd.DataFrame) -> str:
    rows = ct.to_dict("records")
    n = len(rows)
    W, H, L, R, T, B = 1000, 260, 44, 10, 14, 34
    pw = W - L - R
    slot = pw / max(n, 1)
    bw = min(28, slot * 0.68)
    ymax = 0.5
    def y(v): return T + (H - T - B) * (1 - min(v, ymax) / ymax)
    parts = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Heterozygosity rate per chromosome">']
    for g in (0.1, 0.2, 0.3, 0.4, 0.5):
        yy = y(g)
        parts.append(f'<line class="axis" x1="{L}" x2="{W-R}" y1="{yy:.1f}" y2="{yy:.1f}"/><text class="lbl" x="{L-6}" y="{yy+4:.1f}" text-anchor="end">{g:.1f}</text>')
    base = y(0)
    parts.append(f'<line class="axis" x1="{L}" x2="{W-R}" y1="{base:.1f}" y2="{base:.1f}"/>')
    for i, r in enumerate(rows):
        v = r["het_rate"] if r["het_rate"] == r["het_rate"] else 0.0
        x = L + slot * i + (slot - bw) / 2
        top = y(v)
        h = max(base - top, 0)
        rr = min(4, h / 2)
        cls = "bar sex" if r["chrom"] in ("X", "Y", "MT") else "bar"
        tt = f"chr{r['chrom']}: heterozygosity {v:.3f} · {r['called']:,} called · {r['no_call']:,} no-call"
        path = (f"M{x:.1f},{base:.1f} V{top+rr:.1f} Q{x:.1f},{top:.1f} {x+rr:.1f},{top:.1f} H{x+bw-rr:.1f} "
                f"Q{x+bw:.1f},{top:.1f} {x+bw:.1f},{top+rr:.1f} V{base:.1f} Z")
        parts.append(f'<path class="{cls}" d="{path}" data-tt="{e(tt)}"/>')
        parts.append(f'<rect x="{L+slot*i:.1f}" y="{T}" width="{slot:.1f}" height="{H-T-B}" fill="transparent" data-tt="{e(tt)}"/>')
        parts.append(f'<text class="lbl" x="{x+bw/2:.1f}" y="{H-B+16}" text-anchor="middle">{e(r["chrom"])}</text>')
    parts.append("</svg>")
    return '<div class="chart">' + "".join(parts) + '<div class="tt"></div></div>' + \
        '<p class="caption">Share of called sites that are heterozygous, per chromosome. Autosomes typically sit between 0.25 and 0.35 on consumer arrays. ' \
        'X near zero with Y present indicates a male sample; MT is always homozygous.</p>'


def _panel_table(sub: pd.DataFrame) -> str:
    h = ['<div class="tablewrap"><table class="panel"><thead><tr><th>Marker</th><th>Gene</th><th>What it tags</th><th>Genotype</th><th>Effect allele</th><th>Reading</th></tr></thead><tbody>']
    for _, r in sub.iterrows():
        gt = r["genotype"]
        ec = str(r["effect_count"])
        if gt in ("not on array", "no call"):
            badge = f'<span class="badge na">{e(gt)}</span>'
            gtc = ""
        else:
            n = ec.split("/")[0] if ec else ""
            badge = f'<span class="badge n{n}">{e(ec)} {e(r["effect"])}</span>' if ec else ""
            gtc = e(gt)
        interp = e(str(r["interpretation"]))
        note = e(str(r.get("note", "") or ""))
        h.append(f'<tr><td class="num">{e(r["rsid"])}</td><td>{e(r["gene"])}</td><td>{e(r["label"])}</td><td class="gt">{gtc}</td><td>{badge}</td><td>{interp}' + (f'<div class="caption">{note}</div>' if note else "") + "</td></tr>")
    h.append("</tbody></table></div>")
    return "".join(h)


def write_html(out: Path, loaded, qc: dict, ct: pd.DataFrame, roh: pd.DataFrame, roh_text: str,
               panel: pd.DataFrame, combos: list, y: dict, mt: dict, ancestry: dict | None,
               input_name: str = "", extra: dict | None = None) -> Path:
    extra = extra or {}
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    combo_map = {k: (v, n) for k, v, n in combos}
    y_ok = bool(y.get("haplogroup")) and not y.get("error")
    mt_ok = bool(mt.get("haplogroup")) and not mt.get("error")
    H = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
         f'<title>gnomancer report</title><style>{CSS}</style></head><body><div class="wrap">']
    H.append(f'<header><div><h1><span>gnomancer</span> · raw DNA report</h1><div class="meta">{e(input_name)} · generated {now} · gnomancer v{__version__}</div></div>'
             f'<div class="toolbar"><button id="theme">Light / dark</button><button id="reveal">Reveal health sections</button><button id="print">Print / PDF</button></div></header>')
    H.append(f'<div class="disclaimer">{e(DISCLAIMER)}</div>')
    if extra.get("narrative"):
        paras = "".join(f"<p>{e(x.strip())}</p>" for x in extra["narrative"].split("\n\n") if x.strip())
        H.append(f'<section id="plain" class="narrative"><h2>In plain language</h2><div class="card prose">{paras}<div class="caption">Written by your AI assistant from the structured results below; it adds no new facts.</div></div></section>')
    g = extra.get("glance") or {}
    if g:
        H.append('<section id="glance"><h2>At a glance</h2><div class="glance">')
        for title, key in (("Who the data says you are", "who"), ("Worth knowing", "notable"), ("Reassuring", "reassuring"), ("Sensible next steps", "next")):
            if g.get(key):
                H.append(f'<div class="card"><h3>{title}</h3><ul>' + "".join(f"<li>{e(b)}</li>" for b in g[key]) + "</ul></div>")
        H.append("</div></section>")
    H.append('<nav class="toc"><a href="#glance">At a glance</a><a href="#overview">Overview</a><a href="#quality">Quality</a><a href="#roh">Homozygosity</a><a href="#lineage">Lineages</a>'
             '<a href="#readings">Readings</a><a href="#panel">Marker panel</a><a href="#clinvar">ClinVar</a><a href="#pharmcat">PharmCAT</a><a href="#limits">Limits</a></nav>')
    # Tiles
    vendor = e(loaded.source)
    tiles = [_tile("Vendor", f'<span class="{"long" if len(loaded.source) > 9 else ""}">{vendor}</span>', f"build GRCh{loaded.build}" + (" (remapped)" if loaded.remapped else "")),
             _tile("Sites", f"{loaded.count:,}", f"{qc['no_call_rate']:.2%} no-call"),
             _tile("Sex (X/Y)", e(qc["sex_inferred"])),
             _tile("Y haplogroup", e(y["haplogroup"]) if y_ok else "—", "paternal line" if y_ok else (y.get("error") or "")[:60]),
             _tile("mtDNA haplogroup", e(mt["haplogroup"]) if mt_ok else "—", "maternal line" if mt_ok else (mt.get("error") or "")[:60])]
    if "ABO blood group" in combo_map:
        tiles.append(_tile("Blood group", e(combo_map["ABO blood group"][0].split(" (")[0]), "Rh not typed"))
    H.append(f'<section id="overview"><h2>Overview</h2><div class="tiles">{"".join(tiles)}</div>')
    for f in qc["flags"]:
        H.append(f'<div class="flag warn"><span class="ico">!</span><span><b>Check:</b> {e(f)}</span></div>')
    if not qc["flags"]:
        H.append('<div class="flag good"><span class="ico">✓</span><span><b>Quality checks passed.</b> Call rate, heterozygosity and sex markers are all in the expected range.</span></div>')
    for n in loaded.notes:
        H.append(f'<div class="flag warn"><span class="ico">i</span><span>{e(n)}</span></div>')
    H.append("</section>")
    # Quality
    H.append('<section id="quality"><h2>Quality</h2><p class="lead">How complete and internally consistent the file is. These numbers say nothing about you; they say whether the rest of the report can be trusted.</p>')
    H.append(f'<dl class="kv"><dt>No-calls</dt><dd>{qc["no_call"]:,} of {qc["total"]:,} ({qc["no_call_rate"]:.2%})</dd>'
             f'<dt>Indel probes</dt><dd>{qc["indel_calls"]:,}</dd><dt>Autosomal heterozygosity</dt><dd>{qc["autosomal_het"]:.3f}</dd>'
             f'<dt>X heterozygosity</dt><dd>{qc["x_het_rate"]:.4f}</dd><dt>Y sites called</dt><dd>{qc["y_called"]:,} / {qc["y_total"]:,}</dd>'
             f'<dt>Sex inferred</dt><dd>{e(qc["sex_inferred"])} (parser: {e(qc["sex_library"] or "n/a")})</dd></dl>')
    H.append(_het_chart(ct))
    H.append("</section>")
    # ROH
    H.append(f'<section id="roh"><h2>Runs of homozygosity</h2><p class="lead">Long stretches where both copies of a chromosome are identical. They appear when parents share recent ancestors, and in populations with founder effects.</p><div class="card">{e(roh_text)}</div>')
    if not roh.empty:
        H.append('<div class="tablewrap"><table><thead><tr><th>Chromosome</th><th>Start</th><th>End</th><th>Length (Mb)</th><th>Sites</th></tr></thead><tbody>')
        for _, r in roh.head(20).iterrows():
            H.append(f'<tr><td>{e(str(r.chrom))}</td><td class="num">{int(r.start):,}</td><td class="num">{int(r.end):,}</td><td class="num">{r.length_mb:.2f}</td><td class="num">{int(r.snps):,}</td></tr>')
        H.append("</tbody></table></div>")
    H.append("</section>")
    # Lineages
    H.append('<section id="lineage"><h2>Lineages</h2><p class="lead">Two tiny slices of the genome are inherited along a single line: the Y chromosome from father to son, and mitochondrial DNA from mother to child. Each carries a lineage label (haplogroup) with a known geographic history. They describe one ancestral line each, not your overall ancestry.</p><div class="lineage">')
    H.append('<div class="card"><h3>Paternal line · Y chromosome</h3>')
    if y_ok:
        H.append(f'<div class="hg">{e(y["haplogroup"])}<small>{e(y["ycc"])}</small></div><div class="caption">yhaplo, {y["n_y_snps"]:,} Y sites. Arrays cannot resolve below this branch; Y sequencing would.</div><ul class="chain">')
        yl = y.get("yleaf") or {}
        labels = [y["haplogroup"], y["ycc"]]
        if yl.get("trees"):
            H.append('<div class="caption" style="margin-top:10px"><b>Second opinion, Yleaf on current trees:</b></div><dl class="kv" style="margin:4px 0 6px">')
            names = dict(yfull="YFull v14", ftdna="FamilyTreeDNA", isogg="ISOGG")
            for tree, t in yl["trees"].items():
                H.append(f'<dt>{e(names.get(tree, tree))}</dt><dd><span class="gt">{e(t["haplogroup"])}</span> <span class="caption">via {e(t["marker"])}, QC {e(t["qc"])}</span></dd>')
                labels.append(t["haplogroup"].split("*")[0].split("(")[0])
            H.append("</dl>")
        for r in lookup("Y", *labels):
            H.append(f'<li><b>{e(r["name"])}</b>{e(r["text"])}</li>')
        H.append("</ul>")
    else:
        H.append(f'<div class="caption">Not called: {e(y.get("error") or "unknown")}</div>')
    H.append('</div><div class="card"><h3>Maternal line · mitochondrial DNA</h3>')
    if mt_ok:
        H.append(f'<div class="hg">{e(mt["haplogroup"])}<small>quality {mt["quality"]:.2f}</small></div><div class="caption">Haplogrep3, {mt["n_mt_snps"]} typed sites. Differences from rCRS: {e(", ".join(mt["diffs"]) or "none")}.</div><ul class="chain">')
        for r in lookup("MT", mt["haplogroup"]):
            H.append(f'<li><b>{e(r["name"])}</b>{e(r["text"])}</li>')
        H.append("</ul>")
    else:
        H.append(f'<div class="caption">Not called: {e(mt.get("error") or "unknown")}</div>')
    H.append("</div></div>")
    if ancestry:
        from .summary import POP_NAMES
        pp, sp = ancestry.get("population_code"), ancestry.get("superpopulation_code")
        spp, ppp = ancestry.get("superpopulation_percent"), ancestry.get("population_percent")
        pct = (f" {spp:.0%} of the nearest reference samples are {e(POP_NAMES.get(sp, sp))}; {ppp:.0%} are specifically {e(POP_NAMES.get(pp, pp))}."
               if isinstance(spp, (int, float)) and isinstance(ppp, (int, float)) else "")
        H.append(f'<div class="card"><h3>Population estimate · ezancestry</h3><div class="hg">{e(POP_NAMES.get(sp, sp))}<small>nearest group: {e(POP_NAMES.get(pp, pp))}</small></div>'
                 f'<div class="caption">Nearest-neighbour classification against 1000 Genomes reference populations.{pct} This places you near a reference group; it is not an admixture breakdown.</div></div>')
    H.append("</section>")
    # Readings
    H.append('<section id="readings"><h2>Readings that need more than one site</h2><div class="readings">')
    for name, val, note in combos:
        if name == "APOE":
            continue
        val, detail = _split_detail(val)
        H.append(f'<div class="card"><h3>{e(name)}</h3><div class="value">{e(val)}</div><div class="note">{e(note)}</div>' + (f'<div class="caption">{e(detail)}</div>' if detail else "") + "</div>")
    H.append("</div>")
    if "APOE" in combo_map:
        v, n = combo_map["APOE"]
        v, detail = _split_detail(v)
        n = n + (f" · {detail}" if detail else "")
        H.append('<details class="sensitive"><summary>APOE (Alzheimer\'s and lipid genotype)<small>Collapsed on purpose. Some people prefer not to learn this result; open it only if you want to.</small></summary>'
                 f'<div class="readings"><div class="card"><h3>APOE</h3><div class="value">{e(v)}</div><div class="note">{e(n)}</div></div></div></details>')
    H.append("</section>")
    # Panel
    H.append('<section id="panel"><h2>Marker panel</h2><p class="lead">Each row is one well-studied site. The badge counts how many copies of the <em>effect allele</em> you carry (0, 1 or 2) on the GRCh37 forward strand. "Not on array" means this vendor\'s chip does not type the site.</p>'
             '<div class="filter"><input id="q" type="search" placeholder="Filter by gene, rsID or keyword…" aria-label="Filter markers"></div>')
    for cat in CATEGORY_ORDER:
        sub = panel[panel.category == cat]
        if sub.empty:
            continue
        block = f'<div class="catblock"><h3>{e(CAT_TITLES.get(cat, cat))}</h3><p class="lead">{e(CAT_BLURB.get(cat, ""))}</p>{_panel_table(sub)}</div>'
        if cat in SENSITIVE:
            H.append(f'<details class="sensitive"><summary>{e(CAT_TITLES[cat])}<small>Collapsed on purpose: contains disease-risk and carrier results. Open when you are ready, or use "Reveal health sections" above.</small></summary>{block}</details>')
        else:
            H.append(block)
    H.append("</section>")
    # ClinVar
    cv = extra.get("clinvar") or {}
    hits = cv.get("hits")
    H.append('<section id="clinvar"><h2>ClinVar scan of every typed site</h2><p class="lead">The curated panel above covers about a hundred sites. This pass checks every site on the chip against ClinVar, the public archive of clinically classified variants, and lists any pathogenic or likely-pathogenic allele you carry. Each hit is scored for plausibility, because arrays are unreliable at rare sites.</p>')
    body = f'<div class="card">{e(extra.get("clinvar_text", "not run"))}</div>'
    if hits is not None and not hits.empty:
        body += '<div class="tablewrap"><table><thead><tr><th>Marker</th><th>Gene</th><th>Genotype</th><th>Zygosity</th><th>ClinVar</th><th>Review</th><th>Confidence</th><th>Condition / flags</th></tr></thead><tbody>'
        for _, r in hits.iterrows():
            stars = "★" * int(r.stars) + "☆" * (4 - int(r.stars))
            body += (f'<tr><td class="num">{e(r.rsid)}</td><td>{e(r.gene)}</td><td class="gt">{e(r.genotype)}</td><td>{e(r.zygosity)}</td><td>{e(r.significance)}</td>'
                     f'<td title="{e(r.review)}">{stars}</td><td><span class="badge {"n2" if r.confidence == "high" else ("n1" if r.confidence == "medium" else "na")}">{e(r.confidence)}</span></td>'
                     f'<td>{e(r.condition)}' + (f'<div class="caption">{e(r["flags"])}</div>' if r["flags"] else "") + "</td></tr>")
        body += "</tbody></table></div>"
    H.append(f'<details class="sensitive"><summary>ClinVar findings<small>Collapsed on purpose: may contain carrier states and disease-risk variants. Open when ready.</small></summary>{body}</details></section>')
    # PharmCAT
    pc = extra.get("pharmcat") or {}
    calls = pc.get("calls")
    H.append('<section id="pharmcat"><h2>PharmCAT pharmacogenomics</h2><p class="lead">PharmCAT (PharmGKB) assigns star-allele diplotypes and CPIC phenotypes per gene from the sites the chip typed. It is the same engine used in clinical pipelines, but here it sees only a fraction of each gene\'s defining positions.</p>')
    H.append(f'<div class="card">{e(extra.get("pharmcat_text", "not run"))}</div>')
    if calls is not None and not calls.empty:
        from .summary import PGX_DRUGS
        H.append('<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Drugs it matters for</th><th>Positions typed</th><th>Coverage</th></tr></thead><tbody>')
        for _, r in calls.iterrows():
            pheno = str(r.phenotype).strip().rstrip("/").strip() or "—"
            drugs = PGX_DRUGS.get(r.gene, "")
            cov = ('<span class="badge na" title="untyped positions were assumed reference">partial</span>' if r.missing_positions
                   else '<span class="badge n2" title="every defining position was typed">complete</span>')
            H.append(f'<tr><td><b>{e(r.gene)}</b></td><td style="font-family:var(--mono);max-width:280px">{e(r.diplotype)}</td><td>{e(pheno)}</td><td><span class="caption">{e(drugs)}</span></td><td class="num">{e(r.typed)}</td><td>{cov}</td></tr>')
        H.append("</tbody></table></div>")
        H.append('<p class="caption">RYR1/CACNA1S "Uncertain Susceptibility" means a variant of uncertain significance for malignant hyperthermia was seen, not a positive result. The CFTR row describes ivacaftor response, not cystic fibrosis status.</p>')
    H.append("</section>")
    H.append('<section id="limits"><h2>What this file cannot tell you</h2><ul class="limits">'
             '<li>Rare pathogenic variants. The chip targets common variation; genes like BRCA1/2 and CFTR are covered at a handful of sites only.</li>'
             '<li>Copy-number genes such as CYP2D6, and full HLA types beyond tag SNPs.</li>'
             '<li>Admixture percentages at the resolution of the vendor\'s own report, which uses a proprietary reference panel.</li>'
             '<li>Anything at sites the vendor did not type.</li></ul></section>')
    H.append(f'<footer>Generated locally by gnomancer v{__version__}. No data left this computer. Source code: https://github.com/codebase/gnomancer</footer>')
    H.append(f"</div><script>{JS}</script></body></html>")
    out.write_text("".join(H))
    return out
