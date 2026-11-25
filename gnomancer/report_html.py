"""Self-contained HTML report: one file, inline CSS/JS, light and dark, printable, no external requests."""
from __future__ import annotations

import datetime as dt
from html import escape as e
from pathlib import Path

import pandas as pd

from . import __version__
from .context import lookup, page_url
from .panel import CATEGORY_ORDER, POPS, POP_LABEL, fmt_share
from .report import DISCLAIMER

CAT_TITLES = dict(ancestry="Ancestry-informative markers", trait="Traits",
                  health="Health-associated common variants", carrier="Carrier and rare pathogenic sites",
                  pgx="Pharmacogenomic sites (single markers)")
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
nav.toc{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:14px;margin-bottom:24px;position:sticky;top:0;z-index:5;background:var(--bg);padding:10px 0;border-bottom:1px solid var(--line)}nav.toc a{color:var(--accent-ink);text-decoration:none}nav.toc a:hover{text-decoration:underline}
section h2 a.top{font-size:12px;font-weight:400;color:var(--ink-3);text-decoration:none;margin-left:10px;vertical-align:middle}section h2 a.top:hover{color:var(--accent-ink)}
.rs{display:none}.scrollhint{display:none;font-size:12.5px;color:var(--ink-3);margin:0 0 6px}
.chart.mob{display:none}
.glossary dt{font-weight:600;margin-top:8px}.glossary dd{margin:0;color:var(--ink-2);font-size:14px}
.methods td:first-child{font-weight:600;white-space:nowrap}
@media (max-width:720px){table.panel{min-width:760px}table.panel th:first-child,table.panel td:first-child{display:none}.rs{display:block;font-family:var(--mono);font-size:11.5px;color:var(--ink-3)}.scrollhint{display:block}.chart.desk{display:none}.chart.mob{display:block}.card .hg small{display:block;margin:2px 0 0}nav.toc{position:static}}
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
.badge{display:inline-block;min-width:34px;text-align:center;padding:1px 8px;border-radius:999px;font-size:12.5px;font-family:var(--mono);border:1px solid var(--line);color:var(--ink-2);white-space:nowrap}
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
.axis{stroke:var(--line)}.lbl{fill:var(--ink-2);font-family:var(--sans)}
.caption{font-size:13px;color:var(--ink-3);margin:6px 2px 0}
.limits li{margin:4px 0}
.glance{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}.glance .wide{grid-column:1/-1}.glance .card ul{margin:6px 0 0;padding-left:18px}.glance .card li{margin:5px 0;font-size:14.5px}
.narrative .prose p{margin:0 0 10px;font-size:15.5px;line-height:1.6;max-width:72ch}
footer{margin-top:48px;padding-top:16px;border-top:1px solid var(--line);color:var(--ink-2);font-size:13.5px}
.disclaimer{background:var(--surface);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:8px;padding:12px 14px;font-size:14px;color:var(--ink-2);margin:0 0 20px}
.kv{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px;font-size:14px}.kv dt{color:var(--ink-2)}.kv dd{margin:0}
.hidden{display:none!important}
@media print{body{font-size:12px}.toolbar,nav.toc,.filter,.chart .tt,section h2 a.top,.scrollhint{display:none}details.sensitive{border:0;padding:0}details.sensitive[open]>summary{display:none}section{break-inside:avoid-page}table{break-inside:auto}tr{break-inside:avoid}}
"""

JS = """
(function(){
  var root=document.documentElement;
  try{var t=localStorage.getItem('gnomancer-theme');if(t)root.setAttribute('data-theme',t);}catch(e){}
  var tb=document.getElementById('theme');if(tb)tb.onclick=function(){var cur=root.getAttribute('data-theme');var dark=cur?cur==='dark':matchMedia('(prefers-color-scheme: dark)').matches;var nxt=dark?'light':'dark';root.setAttribute('data-theme',nxt);try{localStorage.setItem('gnomancer-theme',nxt);}catch(e){}};
  var pb=document.getElementById('print');if(pb)pb.onclick=function(){window.print();};
  var rb=document.getElementById('reveal');if(rb)rb.onclick=function(){document.querySelectorAll('details.sensitive').forEach(function(d){d.open=true;});};
  function openTarget(){var h=location.hash;if(!h||h.length<2)return;var el=document.getElementById(h.slice(1));if(!el)return;var d=el.closest('details');while(d){d.open=true;d=d.parentElement&&d.parentElement.closest('details');}el.querySelectorAll('details.sensitive').forEach(function(x){x.open=true;});}
  window.addEventListener('hashchange',openTarget);openTarget();
  document.querySelectorAll('section h2').forEach(function(h){var a=document.createElement('a');a.href='#top';a.className='top';a.textContent='\u2191 top';h.appendChild(a);});
  var wasOpen=[];window.addEventListener('beforeprint',function(){wasOpen=[];document.querySelectorAll('details.sensitive').forEach(function(d){wasOpen.push(d.open);d.open=true;});});
  window.addEventListener('afterprint',function(){document.querySelectorAll('details.sensitive').forEach(function(d,i){d.open=!!wasOpen[i];});});
  var q=document.getElementById('q');if(q)q.oninput=function(){var v=q.value.trim().toLowerCase();document.querySelectorAll('table.panel tbody tr').forEach(function(tr){tr.classList.toggle('hidden',v&&tr.textContent.toLowerCase().indexOf(v)<0);});document.querySelectorAll('table.panel').forEach(function(t){var any=t.querySelectorAll('tbody tr:not(.hidden)').length>0;t.closest('.catblock').classList.toggle('hidden',!any);});};
  document.querySelectorAll('.chart').forEach(function(c){var tt=c.querySelector('.tt');c.querySelectorAll('[data-tt]').forEach(function(el){el.addEventListener('mousemove',function(ev){tt.textContent=el.getAttribute('data-tt');tt.style.display='block';var r=c.getBoundingClientRect();tt.style.left=(ev.clientX-r.left)+'px';tt.style.top=(ev.clientY-r.top-8)+'px';});el.addEventListener('mouseleave',function(){tt.style.display='none';});});});
})();
"""


GLOSSARY = [
    ("Haplogroup", "A branch of the human family tree defined by shared mutations on the Y chromosome (paternal line) or mitochondrial DNA (maternal line). It describes one ancestral line, not overall ancestry."),
    ("Site / marker / SNP", "One position in the genome where people differ. A consumer chip types a fixed list of around 600,000 of them; this report calls them sites."),
    ("rsID", "The public catalogue number for a site, e.g. rs4988235. Useful for looking a site up in SNPedia or ClinVar."),
    ("Genotype", "The two letters you carry at a site, one from each parent (one letter on X or Y in males)."),
    ("Effect allele", "The version of a site that the published association is about. The badge counts how many copies you carry: 0, 1 or 2."),
    ("Heterozygous / homozygous / hemizygous", "Two different letters at a site / two of the same / only one copy present (X and Y in males, and mitochondrial DNA)."),
    ("Carrier", "One copy of a variant for a recessive condition. Carriers are usually unaffected; the result matters mainly for children, if the other parent carries one too."),
    ("Odds ratio (e.g. ~2x)", "How much a variant multiplies the odds of a condition in population studies. It is relative, not absolute: 2x of a 1% risk is 2%."),
    ("Run of homozygosity", "A long stretch where both chromosome copies match, which happens when parents share recent ancestors."),
    ("Diplotype / star allele", "A pharmacogenomics term for the pair of gene versions you carry, named *1, *2, *3 and so on; the pair determines a metabolizer phenotype."),
    ("Phenotype (metabolizer)", "What the diplotype predicts: normal, intermediate, poor or rapid metabolism of certain drugs."),
    ("ClinVar", "The public NIH archive of variants and their clinical classifications, with a review status from zero to four stars."),
    ("GRCh37 / build 37", "The version of the human reference genome that consumer vendors report positions against; files on build 38 are remapped."),
    ("1000 Genomes", "A public reference set of about 2,500 genomes from 26 populations; the source of the \"how common\" percentages."),
]


def _methods(loaded, y, mt, ancestry, cv: dict, pc: dict) -> list:
    """Rows for the methods table: (step, tool, reference)."""
    from importlib.metadata import version as _v
    def ver(pkg):
        try:
            return _v(pkg)
        except Exception:
            return "?"
    from .setup_cmd import HAPLOGREP_VERSION, PHARMCAT_VERSION
    seq = getattr(loaded, "kind", "array") != "array"
    if seq:
        from .seqvcf import HET_AB, HOM_AB, MIN_DP, MIN_GQ
        rows = [("Reading the VCF", f"gnomancer {__version__}", f"one pass, native GRCh{loaded.build}; alleles trimmed per alternate allele; calls need the caller's PASS, at least {MIN_DP} reads, "
                 f"genotype quality {MIN_GQ}+, and {HET_AB[0]:.0%}-{HET_AB[1]:.0%} of reads on a heterozygous allele ({HOM_AB:.0%}+ if homozygous)"),
                ("Marker panel and frequencies", f"gnomancer {__version__}", "curated sites genotyped by position (Ensembl and ClinVar coordinates); Hardy-Weinberg from 1000 Genomes phase 3 allele frequencies")]
    else:
        rows = [("Parsing and build detection", f"snps {ver('snps')}", "github.com/apriha/snps; remapped to GRCh37 when needed"),
                ("Marker panel and frequencies", f"gnomancer {__version__}", "111 curated sites; Hardy-Weinberg from 1000 Genomes phase 3 allele frequencies (Ensembl GRCh37)")]
    if seq and y.get("haplogroup") and not y.get("error"):
        rows.append(("Y haplogroup", "Yleaf", "YFull, FamilyTreeDNA and ISOGG trees"))
    elif y.get("haplogroup") and not y.get("error"):
        rows.append(("Y haplogroup", f"yhaplo {ver('yhaplo')}", "23andMe; ISOGG-based tree"))
        if (y.get("yleaf") or {}).get("trees"):
            rows.append(("Y second opinion", "Yleaf", "YFull v14, FamilyTreeDNA and ISOGG trees"))
    if mt.get("haplogroup") and not mt.get("error"):
        rows.append(("mtDNA haplogroup", f"Haplogrep3 {HAPLOGREP_VERSION}", f"PhyloTree 17.3 (rCRS); quality {mt.get('quality', 0):.2f}"))
    if ancestry:
        rows.append(("Population estimate", "ezancestry", "k-nearest neighbours on 1000 Genomes PCA"))
    if cv.get("available"):
        rows.append(("ClinVar scan", f"ClinVar GRCh{loaded.build if seq else 37} VCF, release {cv.get('release') or 'unknown'}",
                     "NCBI; exact position-and-allele match, indels included; actionable genes per ACMG SF v3.3" if seq else "NCBI; plausibility flags per Tandy-Connor et al. 2018"))
    if pc.get("available"):
        rows.append(("Pharmacogenomics", f"PharmCAT {PHARMCAT_VERSION}", "PharmGKB / CPIC; " + ("absent positions taken as reference, as PharmCAT's own --absent-to-ref does" if seq else "untyped positions assumed reference")))
    for x in (getattr(loaded, "methods", None) or []):
        rows.append(tuple(x))
    rows.append(("Report", f"gnomancer {__version__}", f"generated from {loaded.source} data, build GRCh{loaded.build}"))
    return rows


def _split_detail(val: str):
    """'A (A/O) (rs8176719 DI, rs8176746 GG)' -> ('A (A/O)', 'rs8176719 DI, rs8176746 GG')"""
    i = val.find(" (rs")
    if i < 0:
        return val, ""
    return val[:i], val[i + 2:].rstrip(")")


def _tile(label, value, sub=""):
    return f'<div class="tile"><div class="label">{e(label)}</div><div class="value">{value}</div>' + (f'<div class="sub">{e(sub)}</div>' if sub else "") + "</div>"


def _het_chart(ct: pd.DataFrame, seq: bool = False) -> str:
    """Two renderings of the same chart: a wide one for desktop and a narrower viewBox (so labels stay legible) for phones."""
    if seq:   # a VCF lists variants only, so the share is of variants, and sits near 0.6 instead of 0.3
        return (_het_svg(ct, "desk", 1000, 260, 11, 1.0, "variants") + _het_svg(ct, "mob", 520, 300, 13, 1.0, "variants")
                + '<p class="caption">Share of each chromosome\'s SNVs that are heterozygous. Autosomes typically sit between 0.55 and 0.67 in a whole genome; '
                  'a chromosome well below the rest carries a long homozygous stretch. X near zero indicates a male sample; MT is not diploid.</p>')
    return (_het_svg(ct, "desk", 1000, 260, 11) + _het_svg(ct, "mob", 520, 300, 13)
            + '<p class="caption">Share of called sites that are heterozygous, per chromosome. Autosomes typically sit between 0.25 and 0.35 on consumer arrays. '
              'X near zero with Y present indicates a male sample; MT is always homozygous.</p>')


def _het_svg(ct: pd.DataFrame, variant: str, W: int, H: int, font: int, ymax: float = 0.5, noun: str = "called") -> str:
    rows = [r for r in ct.to_dict("records") if r["chrom"] != "PAR"]  # the pseudo-autosomal rows stay in the TSV; they only crowd the axis
    n = len(rows)
    L, R, T, B = 44, 10, 14, 34
    pw = W - L - R
    slot = pw / max(n, 1)
    bw = min(28, slot * 0.68)
    def y(v): return T + (H - T - B) * (1 - min(v, ymax) / ymax)
    parts = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Heterozygosity rate per chromosome" style="font-size:{font}px">']
    for g in ((0.2, 0.4, 0.6, 0.8, 1.0) if ymax > 0.5 else (0.1, 0.2, 0.3, 0.4, 0.5)):
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
        tt = (f"chr{r['chrom']}: {v:.3f} of {r['called']:,} SNVs heterozygous · {r['no_call']:,} calls failed quality" if noun == "variants"
              else f"chr{r['chrom']}: heterozygosity {v:.3f} · {r['called']:,} called · {r['no_call']:,} no-call")
        path = (f"M{x:.1f},{base:.1f} V{top+rr:.1f} Q{x:.1f},{top:.1f} {x+rr:.1f},{top:.1f} H{x+bw-rr:.1f} "
                f"Q{x+bw:.1f},{top:.1f} {x+bw:.1f},{top+rr:.1f} V{base:.1f} Z")
        parts.append(f'<path class="{cls}" d="{path}" data-tt="{e(tt)}"/>')
        parts.append(f'<rect x="{L+slot*i:.1f}" y="{T}" width="{slot:.1f}" height="{H-T-B}" fill="transparent" data-tt="{e(tt)}"/>')
        parts.append(f'<text class="lbl" x="{x+bw/2:.1f}" y="{H-B+16}" text-anchor="middle">{e(r["chrom"])}</text>')
    parts.append("</svg>")
    return f'<div class="chart {variant}">' + "".join(parts) + '<div class="tt"></div></div>'


def _share_cell(r, pop: str) -> str:
    v = r.get(f"gf_{pop}")
    if not isinstance(v, float) or v != v:  # missing or NaN: no reference frequency for this site
        return '<td class="num" title="no reference frequency for this site"><span class="caption">—</span></td>'
    parts = [f"{POP_LABEL[p]}: {fmt_share(r.get(f'gf_{p}'))}" for p in POPS if isinstance(r.get(f"gf_{p}"), float) and r.get(f"gf_{p}") == r.get(f"gf_{p}")]
    return f'<td class="num" title="{e("; ".join(parts))}">{e(fmt_share(v))}</td>'


def _panel_table(sub: pd.DataFrame, pop: str = "EUR") -> str:
    absent = sub.genotype.isin(["not on array", "no call"])
    sub = pd.concat([sub[~absent], sub[absent]])  # typed sites first, untyped ones at the bottom
    h = [f'<p class="scrollhint">Swipe sideways for genotype, effect allele and reading.</p><div class="tablewrap"><table class="panel"><thead><tr><th>Marker</th><th>Gene</th><th>What it tags</th><th>Genotype</th><th>Effect allele</th><th title="Share of the {e(POP_LABEL.get(pop, pop))} 1000 Genomes reference population expected to carry the same genotype. Hover for all five super-populations.">How common ({e(pop)})</th><th>Reading</th></tr></thead><tbody>']
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
            if str(r.get("basis", "")).startswith("not in the VCF"):
                gtc += ' <span class="caption" title="No variant is recorded here, so the genotype is taken to be the reference on both copies. A gVCF would confirm it.">°</span>'
        interp = e(str(r["interpretation"]))
        if gt == "no call" and isinstance(r.get("basis"), str) and r.get("basis"):
            interp = e(r["basis"])
        note = e(str(r.get("note", "") or ""))
        h.append(f'<tr><td class="num">{e(r["rsid"])}</td><td>{e(r["gene"])}<div class="rs">{e(r["rsid"])}</div></td><td>{e(r["label"])}</td><td class="gt">{gtc}</td><td>{badge}</td>{_share_cell(r, pop)}<td>{interp}' + (f'<div class="caption">{note}</div>' if note else "") + "</td></tr>")
    h.append("</tbody></table></div>")
    return "".join(h)


SEQ_DISCLAIMER = (
    "Computed locally from a sequencing VCF. Sequencing reads far more than an array, but a VCF of small variants still cannot show "
    "large deletions and duplications, repeat expansions or genes with near-identical copies, and a variant no database has classified "
    "passes unremarked. Research-grade sequencing is not a clinical test. This is not medical advice: discuss anything concerning with a "
    "clinician or genetic counselor and confirm it with a clinical-grade test before acting on it."
)
SEQ_PANEL_NOTE = ('A small ° after a genotype means no variant is recorded at that position, so both copies are taken to match the reference; '
                  'a call that failed the quality gates is shown as "no call" with the reason.')


def seq_limits(loaded, qc: dict, extra: dict, pop: str) -> list:
    done = {x["id"] for x in extra.get("sections") or []}
    out = []
    if not qc.get("whole_genome"):
        out.append("Anything outside the regions this file covers. It is not a whole genome, so a missing variant means 'not looked at', never 'normal'.")
    elif qc.get("reads"):
        out.append("Whether positions beyond the named ones were read. The reads were checked at every marker-panel and PharmCAT position the VCF is silent about; "
                   "elsewhere (the ClinVar scan) absence of a variant is still taken to mean the reference.")
    elif not qc.get("gvcf"):
        out.append("Whether a position with no variant was actually read. The VCF lists only differences, so 'reference' is inferred from absence; a gVCF, or a coverage check on the BAM, would confirm it.")
    out.append("Variants nobody has classified. The ClinVar scan finds known pathogenic variants; a rare variant unique to your family is not in it, and predicting its effect needs a clinical lab.")
    if "sv" not in done:
        out.append("Large deletions, duplications and rearrangements (structural and copy-number variants). They are not in a small-variant VCF; pass the vendor's SV and CNV files with --sv and --cnv.")
    if "cyp2d6" not in done:
        out.append("CYP2D6, which is deleted, duplicated or fused with its pseudogene in many people, and other genes with near-identical copies (SMN1, PMS2, CYP21A2). These need the BAM or CRAM and a dedicated caller.")
    if "hla" not in done:
        out.append("HLA types. The marker panel's HLA rows are tag SNPs; real typing needs the BAM or CRAM.")
    if "repeats" not in done:
        out.append("Repeat expansions (Huntington disease, fragile X, myotonic dystrophy and others). They need the BAM or CRAM and a repeat caller, and short reads cannot size the largest ones.")
    if "smn" not in done:
        out.append("Spinal muscular atrophy carrier status. SMN1 has a near-identical neighbour, so its copy number needs the BAM or CRAM and a dedicated caller.")
    if done & {"cyp2d6", "hla", "repeats", "smn", "sv"}:
        out.append("Certainty about the read-level results above. Each comes from a research tool on 30x short reads and is right for most people, not all; "
                   "the clinical tests named in each section are what confirm them.")
    out.append("Low-level mosaic variants, and mitochondrial variants present in under about 3% of reads.")
    out.append('For a male sample, X-chromosome sites have one copy (hemizygous), so "carrier" and "homozygous" do not apply there in the usual sense.')
    out.append(f'"How common" uses the {POP_LABEL.get(pop, pop)} reference population unless --ancestry placed you in another 1000 Genomes group; a different reference changes those percentages, not your genotype.')
    return out


def _seq_clinvar_tables(hits: pd.DataFrame) -> str:
    """ClinVar hits from sequencing, grouped by meaning: actionable genes, both copies, single copies, set aside."""
    from .clinvar import seq_groups
    out = []
    for _, title, lead, df in seq_groups(hits):
        out.append(f'<h3 style="margin:18px 0 4px">{e(title)} <small style="font-weight:400;color:var(--ink-3)">({len(df)})</small></h3><p class="lead">{e(lead)}</p>'
                   '<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Change</th><th>Zygosity</th><th>ClinVar</th><th>Review</th><th title="Read depth, and the share of reads carrying the variant">Evidence</th><th>Condition / notes</th></tr></thead><tbody>')
        for _, r in df.iterrows():
            stars = "★" * int(r.stars) + "☆" * (4 - int(r.stars))
            share = f', {float(r.allele_share):.0%} of reads' if r.allele_share is not None and r.allele_share == r.allele_share else ""
            depth = f'{int(r.depth)} reads' if r.depth is not None and r.depth == r.depth else "depth not given"
            notes = "; ".join(x for x in (r.role if isinstance(r.role, str) else "", r["flags"] if isinstance(r["flags"], str) else "") if x)
            rs = f'<div class="rs">{e(r.rsid)}</div>' if str(r.rsid).startswith("rs") else ""
            cons = f'<div class="caption">{e(r.consequence)}</div>' if isinstance(r.consequence, str) and r.consequence else ""
            out.append(f'<tr><td><b>{e(r.gene)}</b>{rs}</td><td class="gt" title="{e(r.variant)}">{e(r.genotype)}{cons}</td><td>{e(r.zygosity)}</td><td>{e(r.significance)}</td>'
                       f'<td title="{e(r.review)}">{stars}</td><td><span class="caption">{e(depth + share)}</span></td>'
                       f'<td>{e(r.condition)}' + (f'<div class="caption">{e(notes)}</div>' if notes else "") + "</td></tr>")
        out.append("</tbody></table></div>")
    return "".join(out)


def _extra_section(x: dict) -> str:
    """A section supplied by an optional analysis: {id, title, lead, text, columns, rows, notes, sensitive}."""
    body = f'<div class="card">{e(x.get("text") or "")}</div>' if x.get("text") else ""
    if x.get("rows"):
        body += '<div class="tablewrap"><table><thead><tr>' + "".join(f"<th>{e(c)}</th>" for c in x["columns"]) + "</tr></thead><tbody>"
        for row in x["rows"]:
            body += "<tr>" + "".join(f"<td>{e(str(v))}</td>" for v in row) + "</tr>"
        body += "</tbody></table></div>"
    for n in x.get("notes") or []:
        body += f'<p class="caption">{e(n)}</p>'
    if x.get("sensitive"):
        body = f'<details class="sensitive"><summary>{e(x["title"])}<small>Collapsed on purpose: may contain disease-risk results. Open when ready.</small></summary>{body}</details>'
    return f'<section id="{e(x["id"])}"><h2>{e(x["title"])}</h2><p class="lead">{e(x.get("lead") or "")}</p>{body}</section>'


def write_html(
out: Path, loaded, qc: dict, ct: pd.DataFrame, roh: pd.DataFrame, roh_text: str,
               panel: pd.DataFrame, combos: list, y: dict, mt: dict, ancestry: dict | None,
               input_name: str = "", extra: dict | None = None) -> Path:
    extra = extra or {}
    seq = getattr(loaded, "kind", "array") != "array"       # sequencing VCF, not an array export
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    combo_map = {k: (v, n) for k, v, n in combos}
    y_ok = bool(y.get("haplogroup")) and not y.get("error")
    mt_ok = bool(mt.get("haplogroup")) and not mt.get("error")
    H = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
         f'<title>gnomancer report · {e(input_name or now)}</title><style>{CSS}</style></head><body><div class="wrap">']
    card = extra.get("card") or {}
    card_link = Path(card["pdf"] or card["html"]).name if card.get("rows") else ""
    H.append(f'<header id="top"><div><h1><span>gnomancer</span> · {"genome" if seq else "raw DNA"} report</h1><div class="meta">{e(input_name)} · generated {now} · gnomancer v{__version__}</div></div>'
             f'<div class="toolbar"><button id="theme">Light / dark</button><button id="reveal">Reveal health sections</button><button id="print">Print / PDF</button>'
             + (f'<a class="btn" href="{e(card_link)}">Pharmacist card</a>' if card_link else "") + '</div></header>')
    H.append(f'<div class="disclaimer">{e(SEQ_DISCLAIMER if seq else DISCLAIMER)}</div>')
    H.append('<nav class="toc">' + ('<a href="#plain">In plain language</a>' if extra.get("narrative") else "") + '<a href="#glance">At a glance</a><a href="#overview">Overview</a><a href="#lineage">Lineages</a>'
             '<a href="#readings">Readings</a><a href="#panel">Marker panel</a><a href="#clinvar">ClinVar</a><a href="#pharmcat">Drug-gene results</a>'
             + "".join(f'<a href="#{e(x["id"])}">{e(x.get("nav") or x["title"])}</a>' for x in extra.get("sections") or []) +
             '<a href="#quality">Quality</a><a href="#roh">Homozygosity</a><a href="#methods">Methods</a><a href="#glossary">Glossary</a><a href="#limits">Limits</a></nav>')
    if extra.get("narrative"):
        paras = "".join(f"<p>{e(x.strip())}</p>" for x in extra["narrative"].split("\n\n") if x.strip())
        H.append(f'<section id="plain" class="narrative"><h2>In plain language</h2><div class="card prose">{paras}<div class="caption">Narrative supplied with --summary or --ai, written from the structured results below. It should add no facts beyond them; if it does, trust the sections, not the prose.</div></div></section>')
    g = extra.get("glance") or {}
    if g:
        H.append('<section id="glance"><h2>At a glance</h2><div class="glance">')
        for title, key in (("Who the data says you are", "who"), ("Worth knowing", "notable"), ("Reassuring", "reassuring"), ("Sensible next steps", "next")):
            if g.get(key):
                H.append(f'<div class="card{" wide" if key == "next" else ""}"><h3>{title}</h3><ul>' + "".join(f"<li>{e(b)}</li>" for b in g[key]) + "</ul></div>")
        H.append("</div></section>")
    # Tiles
    vendor = e(loaded.source)
    tiles = [_tile("Source" if seq else "Vendor", f'<span class="{"long" if len(loaded.source) > 9 else ""}">{vendor}</span>', f"build GRCh{loaded.build}" + (" (remapped)" if loaded.remapped else "")),
             _tile("Variants", f"{loaded.count:,}", f"{qc['no_call_rate']:.1%} more failed quality gates") if seq else _tile("Sites", f"{loaded.count:,}", f"{qc['no_call_rate']:.2%} no-call"),
             _tile("Genetic sex", e(qc["sex_inferred"]), "from X heterozygosity" + ("" if seq else " and Y calls")),
             _tile("Y haplogroup", e(y["haplogroup"]) if y_ok else "—", "paternal line" if y_ok else (y.get("error") or "")[:60]),
             _tile("mtDNA haplogroup", e(mt["haplogroup"]) if mt_ok else "—", "maternal line" if mt_ok else (mt.get("error") or "")[:60])]
    if "ABO blood group" in combo_map:
        tiles.append(_tile("Blood group", e(combo_map["ABO blood group"][0].split(" (")[0]), "Rh not typed"))
    H.append(f'<section id="overview"><h2>Overview</h2><div class="tiles">{"".join(tiles)}</div>')
    for f in qc["flags"]:
        H.append(f'<div class="flag warn"><span class="ico">!</span><span><b>Check:</b> {e(f)}</span></div>')
    if not qc["flags"]:
        H.append('<div class="flag good"><span class="ico">✓</span><span><b>Quality checks passed.</b> ' + ("Variant counts, the transition/transversion ratio and heterozygosity are in the range expected for a whole genome." if seq else "Call rate, heterozygosity and sex markers are all in the expected range.") + '</span></div>')
    for n in loaded.notes:
        H.append(f'<div class="flag warn"><span class="ico">i</span><span>{e(n)}</span></div>')
    H.append("</section>")
    # Lineages
    H.append('<section id="lineage"><h2>Lineages</h2><p class="lead">Two tiny slices of the genome are inherited along a single line: the Y chromosome from father to son, and mitochondrial DNA from mother to child. Each carries a lineage label (haplogroup) with a known geographic history. They describe one ancestral line each, not your overall ancestry.</p><div class="lineage">')
    H.append('<div class="card"><h3>Paternal line · Y chromosome</h3>')
    if y_ok:
        if seq:
            H.append(f'<div class="hg">{e(y["haplogroup"])}<small>{e(y.get("ycc") or "")}</small></div><div class="caption">{e(y.get("basis") or "Yleaf, from the sequenced Y chromosome.")}</div><ul class="chain">')
        else:
            H.append(f'<div class="hg">{e(y["haplogroup"])}<small>{e(y["ycc"])}</small></div><div class="caption">yhaplo, {y["n_y_snps"]:,} Y sites. Arrays cannot resolve below this branch; Y sequencing would.</div><ul class="chain">')
        yl = y.get("yleaf") or {}
        labels = [y["haplogroup"], y.get("ycc") or ""]
        if yl.get("trees"):
            H.append('<div class="caption" style="margin-top:10px"><b>' + ("Yleaf on current trees:" if seq else "Second opinion, Yleaf on current trees:") + '</b></div><dl class="kv" style="margin:4px 0 6px">')
            names = dict(yfull="YFull v14", ftdna="FamilyTreeDNA", isogg="ISOGG")
            for tree, t in yl["trees"].items():
                H.append(f'<dt>{e(names.get(tree, tree))}</dt><dd><span class="gt">{e(t["haplogroup"])}</span> <span class="caption">via {e(t["marker"])}, QC {e(t["qc"])}</span></dd>')
                labels.append(t["haplogroup"].split("*")[0].split("(")[0])
            H.append('</dl><div class="caption">Each tree names branches its own way (ISOGG short name, FTDNA, YFull), so the labels differ even when they point at the same lineage. Read them together: ' + ("they describe one branch, not three." if seq else "at array resolution they describe one branch, not three.") + '</div>')
        for r in lookup("Y", *labels):
            H.append(f'<li><b>{e(r["name"])}</b>{e(r["text"])} <a href="{e(page_url("Y", r["name"]))}">More</a></li>')
        H.append("</ul>")
    else:
        H.append(f'<div class="caption">Not called: {e(y.get("error") or "unknown")}</div>')
    H.append('</div><div class="card"><h3>Maternal line · mitochondrial DNA</h3>')
    if mt_ok:
        if mt.get("sequenced"):
            het = mt.get("heteroplasmies") or []
            het_txt = (" Carried by only part of the reads (heteroplasmy), and left out of the haplogroup call: " + ", ".join(f'{h["pos"]}{h["ref"]}>{h["alt"]} at {h["share"]:.0%}' for h in het[:12]) + ("…" if len(het) > 12 else "") + ".") if het else ""
            H.append(f'<div class="hg">{e(mt["haplogroup"])}<small>quality {mt["quality"]:.2f}</small></div><div class="caption">Haplogrep3, from the whole mitochondrial sequence: {mt["n_mt_snps"]} differences from rCRS ({e(", ".join(mt["diffs"] + mt["indels"]) or "none")}).{e(het_txt)}</div><ul class="chain">')
        else:
            H.append(f'<div class="hg">{e(mt["haplogroup"])}<small>quality {mt["quality"]:.2f}</small></div><div class="caption">Haplogrep3, {mt["n_mt_snps"]} typed sites. Differences from rCRS: {e(", ".join(mt["diffs"]) or "none")}.</div><ul class="chain">')
        for r in lookup("MT", mt["haplogroup"]):
            H.append(f'<li><b>{e(r["name"])}</b>{e(r["text"])} <a href="{e(page_url("MT", r["name"]))}">More</a></li>')
        H.append("</ul>")
    else:
        H.append(f'<div class="caption">Not called: {e(mt.get("error") or "unknown")}</div>')
    H.append("</div></div>")
    if ancestry:
        from .summary import POP_NAMES
        pp, sp = ancestry.get("population_code"), ancestry.get("superpopulation_code")
        spp, ppp = ancestry.get("superpopulation_percent"), ancestry.get("population_percent")
        pct = (f" {spp:.0%} of the nearest reference samples are {e(POP_NAMES.get(sp, sp))}; {ppp:.0%} are specifically {e(POP_NAMES.get(pp, pp))}, the single most-represented group among those neighbours (related groups account for the rest)."
               if isinstance(spp, (int, float)) and isinstance(ppp, (int, float)) else "")
        H.append(f'<div class="card"><h3>Population estimate · ezancestry</h3><div class="hg">{e(POP_NAMES.get(sp, sp))}<small>nearest group: {e(POP_NAMES.get(pp, pp))}</small></div>'
                 f'<div class="caption">Nearest-neighbour classification against 1000 Genomes reference populations.{pct} This places you near a reference group; it is not an admixture breakdown.</div></div>')
    H.append("</section>")
    # Readings
    H.append('<section id="readings"><h2>Combined readings</h2><p class="lead">Results that depend on more than one site read together: blood group, lactase persistence, bitter taste and APOE.</p><div class="readings">')
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
        H.append('<details class="sensitive"><summary>APOE (Alzheimer\'s risk genotype)<small>Collapsed on purpose. Some people prefer not to learn this result; open it only if you want to.</small></summary>'
                 f'<div class="readings"><div class="card"><h3>APOE</h3><div class="value">{e(v)}</div><div class="note">{e(n)}</div></div></div></details>')
    H.append("</section>")
    # Panel
    pop = extra.get("pop") or "EUR"
    H.append(f'<section id="panel"><h2>Marker panel</h2><p class="lead">Each row is one well-studied site. The badge counts how many copies of the <em>effect allele</em> you carry (0, 1 or 2) on the forward strand. The share column says what fraction of the {e(POP_LABEL.get(pop, pop))} reference population (1000 Genomes) is expected to carry the same genotype; hover it for all five world regions. ' + (SEQ_PANEL_NOTE if seq else '"Not on array" means this vendor\'s chip does not type the site.') + '</p>'
             '<div class="filter"><input id="q" type="search" placeholder="Filter by gene, rsID or keyword…" aria-label="Filter markers"></div>')
    for cat in CATEGORY_ORDER:
        sub = panel[panel.category == cat]
        if sub.empty:
            continue
        n_absent = int(sub.genotype.isin(["not on array", "no call"]).sum())
        absent_note = f' <small style="font-weight:400;color:var(--ink-3)">({n_absent} site{"s" if n_absent != 1 else ""} {"without a confident call" if seq else "not on this chip"})</small>' if n_absent else ""
        block = f'<div class="catblock"><h3>{e(CAT_TITLES.get(cat, cat)).replace(" typed on arrays", "" if seq else " typed on arrays")}{absent_note}</h3><p class="lead">{e(CAT_BLURB.get(cat, ""))}</p>{_panel_table(sub, pop)}</div>'
        if cat in SENSITIVE:
            H.append(f'<details class="sensitive"><summary>{e(CAT_TITLES[cat]).replace(" typed on arrays", "" if seq else " typed on arrays")}<small>Collapsed on purpose: contains disease-risk and carrier results. Open when you are ready, or use "Reveal health sections" above.</small></summary>{block}</details>')
        else:
            H.append(block)
    H.append("</section>")
    # ClinVar
    cv = extra.get("clinvar") or {}
    hits = cv.get("hits")
    if seq:
        H.append('<section id="clinvar"><h2>ClinVar scan of every variant</h2><p class="lead">Every variant in the file is compared with ClinVar, the public archive of clinically classified variants, by exact position and alleles, insertions and deletions included. Anything classified pathogenic or likely pathogenic is listed, grouped by what it means. Each call shows its read depth and the share of reads carrying the variant, and a weak call or an unreviewed classification is set aside, not hidden.</p>')
        body = f'<div class="card">{e(extra.get("clinvar_text", "not run"))}</div>' + (_seq_clinvar_tables(hits) if hits is not None and not hits.empty else "")
    else:
        H.append('<section id="clinvar"><h2>ClinVar scan of every typed site</h2><p class="lead">The curated panel above covers about a hundred sites. This pass checks every site on the chip against ClinVar, the public archive of clinically classified variants, and lists any pathogenic or likely-pathogenic allele you carry. Each hit is scored for plausibility, because arrays are unreliable at rare sites.</p>')
        body = f'<div class="card">{e(extra.get("clinvar_text", "not run"))}</div>'
        if hits is not None and not hits.empty:
            body += '<div class="tablewrap"><table><thead><tr><th>Marker</th><th>Gene</th><th>Genotype</th><th>Zygosity</th><th>ClinVar</th><th>Review</th><th>Confidence</th><th>Condition / flags</th></tr></thead><tbody>'
            for _, r in hits.iterrows():
                stars = "★" * int(r.stars) + "☆" * (4 - int(r.stars))
                body += (f'<tr><td class="num">{e(r.rsid)}</td><td>{e(r.gene)}</td><td class="gt">{e(r.genotype)}</td><td>{e(r.zygosity)}</td><td>{e(r.significance)}</td>'
                         f'<td title="{e(r.review)}">{stars}<div class="caption">{e(r.review)}</div></td><td><span class="badge {"n2" if r.confidence == "high" else ("n1" if r.confidence == "medium" else "na")}">{e(r.confidence)}</span></td>'
                         f'<td>{e(r.condition)}' + (f'<div class="caption">{e(r["flags"])}</div>' if r["flags"] else "") + "</td></tr>")
            body += "</tbody></table></div>"
    H.append(f'<details class="sensitive"><summary>ClinVar findings<small>Collapsed on purpose: may contain carrier states and disease-risk variants. Open when ready.</small></summary>{body}</details></section>')
    # PharmCAT
    pc = extra.get("pharmcat") or {}
    calls = pc.get("calls")
    H.append('<section id="pharmcat"><h2>Drug-gene results (PharmCAT)</h2><p class="lead">PharmCAT (PharmGKB) assigns star-allele diplotypes and CPIC phenotypes per gene. It is the CPIC-aligned caller used in clinical research pipelines' + ("; from a whole genome it sees every defining position, insertions and deletions included." if seq else ", but from an array it sees only a fraction of each gene\'s defining positions.") + '</p>')
    H.append(f'<div class="card">{e(extra.get("pharmcat_text", "not run"))}</div>')
    if card_link:
        H.append(f'<p class="caption"><a class="btn" href="{e(card_link)}">Open the one-page card for your pharmacist</a> (same table, printable; "confirm clinically" footer included)</p>')
    if calls is not None and not calls.empty:
        from .summary import PGX_DRUGS, pgx_cells
        H.append('<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Result</th><th>Phenotype</th><th>Drugs it matters for</th><th>Positions ' + ("read" if seq else "typed") + '</th><th>Coverage</th></tr></thead><tbody>')
        for _, r in calls.iterrows():
            result, pheno = pgx_cells(r)
            drugs = PGX_DRUGS.get(r.gene, "")
            cov = ('<span class="badge na" title="untyped positions were assumed reference">partial</span>' if r.missing_positions
                   else '<span class="badge n2" title="every defining position was typed">complete</span>')
            H.append(f'<tr><td><b>{e(r.gene)}</b></td><td style="max-width:280px">{e(result)}</td><td>{e(pheno)}</td><td><span class="caption">{e(drugs)}</span></td><td class="num">{e(r.typed)}</td><td>{cov}</td></tr>')
        H.append("</tbody></table></div>")
        H.append('<p class="caption">RYR1/CACNA1S "Uncertain Susceptibility" means a variant of uncertain significance for malignant hyperthermia was seen, not a positive result. The CFTR row describes ivacaftor response, not cystic fibrosis status.</p>')
    H.append("</section>")
    for x in extra.get("sections") or []:
        H.append(_extra_section(x))
    # Quality
    H.append('<section id="quality"><h2>Quality</h2><p class="lead">How complete and internally consistent the file is. These numbers say nothing about you; they say whether the rest of the report can be trusted.</p>')
    if seq:
        f2 = lambda v, fmt: format(v, fmt) if isinstance(v, (int, float)) and v == v else "n/a"
        H.append(f'<dl class="kv"><dt>Layout</dt><dd>{"gVCF: reference blocks say where the sample was read and matched" if qc.get("gvcf") else "variant-only VCF: positions matching the reference are not listed"}'
                 + (f' · caller {e(qc["caller"])}' if qc.get("caller") else "") + '</dd>'
                 f'<dt>Passing variants</dt><dd>{qc["snvs"]:,} SNVs and {qc["indel_calls"]:,} insertions/deletions <span class="caption">(a whole genome has about 4-5 million in total)</span></dd>'
                 f'<dt title="caller filter, fewer than 8 reads, genotype quality under 20, or a read share that does not fit the genotype">Failed quality gates</dt><dd>{qc["no_call"]:,} ({qc["no_call_rate"]:.1%}) <span class="caption">(not used anywhere in this report)</span></dd>'
                 f'<dt title="A to G or C to T changes versus all others; about 2.0-2.1 across a whole genome. A low value means many false calls.">Transition / transversion</dt><dd>{f2(qc["ti_tv"], ".2f")}</dd>'
                 f'<dt title="Heterozygous SNVs for each homozygous one on the autosomes; typically 1.3-2.0 depending on ancestry.">Het / hom ratio</dt><dd>{f2(qc["het_hom_ratio"], ".2f")}</dd>'
                 f'<dt>Mean depth at variants</dt><dd>{f2(qc["mean_depth"], ".0f")}x</dd>'
                 f'<dt>X heterozygosity</dt><dd>{f2(qc["x_het_rate"], ".3f")} <span class="caption">(share of X variants that are heterozygous, outside the pseudo-autosomal ends)</span></dd>'
                 f'<dt>Genetic sex</dt><dd>{e(qc["sex_inferred"])}</dd></dl>')
        rd = qc.get("reads")
        if rd:
            H.append(f'<h3 style="margin:16px 0 6px">From the reads ({e(rd["file"])})</h3><dl class="kv">'
                     f'<dt>Mean depth</dt><dd>{rd["mean_depth"]:.0f}x <span class="caption">({rd["read_length"]}-base reads; mitochondrion about {rd["mt_depth"]:,.0f}x)</span></dd>'
                     f'<dt title="Depth of X and Y relative to the autosomes: about 1.0 and 0 for a female, 0.5 and 0.2-0.5 for a male">X and Y depth</dt><dd>{rd["x_ratio"]:.2f} and {rd["y_ratio"]:.2f} of autosomal depth: {e(rd["sex"])}</dd>'
                     f'<dt title="Named positions (marker panel and PharmCAT) where the VCF records no variant, checked in the reads">Silent positions checked</dt>'
                     f'<dd>{rd["sites_checked"]:,}: {rd["sites_confirmed"]:,} read as reference, {rd["sites_thin"]:,} too thinly covered, {rd["sites_disagree"]:,} where the reads show a variant the VCF lacks</dd></dl>')
    else:
        H.append(f'<dl class="kv"><dt>No-calls</dt><dd>{qc["no_call"]:,} of {qc["total"]:,} ({qc["no_call_rate"]:.2%})</dd>'
                 f'<dt title="insertion/deletion probes, reported as I/D rather than bases">Indel probes</dt><dd>{qc["indel_calls"]:,} <span class="caption">(insertion/deletion sites, reported as I/D)</span></dd><dt>Autosomal heterozygosity</dt><dd>{qc["autosomal_het"]:.3f}</dd>'
                 f'<dt>X heterozygosity</dt><dd>{qc["x_het_rate"]:.4f}</dd><dt>Y sites called</dt><dd>{qc["y_called"]:,} / {qc["y_total"]:,}</dd>'
                 f'<dt>Genetic sex</dt><dd>{e(qc["sex_inferred"])} (file header: {e(qc["sex_library"] or "n/a")})</dd></dl>')
    H.append(_het_chart(ct, seq))
    H.append("</section>")
    # ROH
    H.append(f'<section id="roh"><h2>Runs of homozygosity</h2><p class="lead">Long stretches where both copies of a chromosome are identical. They appear when parents share recent ancestors, and in populations with founder effects.</p><div class="card">{e(roh_text)}</div>')
    if not roh.empty:
        H.append('<div class="tablewrap"><table><thead><tr><th>Chromosome</th><th>Start</th><th>End</th><th>Length (Mb)</th><th>' + ("Homozygous SNVs" if seq else "Sites") + '</th></tr></thead><tbody>')
        for _, r in roh.head(20).iterrows():
            H.append(f'<tr><td>{e(str(r.chrom))}</td><td class="num">{int(r.start):,}</td><td class="num">{int(r.end):,}</td><td class="num">{r.length_mb:.2f}</td><td class="num">{int(r.snps):,}</td></tr>')
        H.append("</tbody></table></div>")
    H.append("</section>")
    # Methods
    H.append('<section id="methods"><h2>How this was computed</h2><p class="lead">Every step runs on this computer from the raw file. Tool versions and references, for reproducibility.</p><div class="tablewrap"><table class="methods"><thead><tr><th>Step</th><th>Tool</th><th>Reference</th></tr></thead><tbody>')
    for step, tool_, ref in _methods(loaded, y, mt, ancestry, cv, pc):
        H.append(f'<tr><td>{e(step)}</td><td>{e(tool_)}</td><td>{e(ref)}</td></tr>')
    H.append("</tbody></table></div></section>")
    # Glossary
    H.append('<section id="glossary"><h2>Glossary</h2><details class="sensitive"><summary>Terms used in this report<small>Open for one-line definitions.</small></summary><dl class="glossary">'
             + "".join(f"<dt>{e(t)}</dt><dd>{e(d)}</dd>" for t, d in GLOSSARY) + "</dl></details></section>")
    if seq:
        H.append('<section id="limits"><h2>What this file cannot tell you</h2><ul class="limits">' + "".join(f"<li>{e(x)}</li>" for x in seq_limits(loaded, qc, extra, pop)) + "</ul></section>")
    else:
        H.append('<section id="limits"><h2>What this file cannot tell you</h2><ul class="limits">'
                 '<li>Rare pathogenic variants. The chip targets common variation; genes like BRCA1/2 and CFTR are covered at a handful of sites only.</li>'
                 '<li>Copy-number genes such as CYP2D6, and full HLA types beyond tag SNPs.</li>'
                 '<li>Admixture percentages at the resolution of the vendor\'s own report, which uses a proprietary reference panel.</li>'
                 '<li>Anything at sites the vendor did not type. PharmCAT assumes the reference allele at untyped positions, so a "complete" diplotype can still miss a rare allele.</li>'
                 '<li>For a male sample, X-chromosome sites have one copy (hemizygous), so "carrier" and "homozygous" do not apply there in the usual sense.</li>'
                 f'<li>"How common" uses the {e(POP_LABEL.get(pop, pop))} reference population unless --ancestry placed you in another 1000 Genomes group; a different reference changes those percentages, not your genotype.</li></ul></section>')
    H.append(f'<footer>Generated locally by gnomancer v{__version__}. No data left this computer. To share with a clinician, print this page or hand over the pharmacist card ({e(card_link) if card_link else "pgx-card.pdf"}). '
             'Documentation: https://gnomancer.com · Source: https://github.com/codebase/gnomancer</footer>')
    H.append(f"</div><script>{JS}</script></body></html>")
    out.write_text("".join(H))
    return out
