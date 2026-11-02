"""Generate the haplogroup story pages under site/haplogroups/ from gnomancer/data/haplogroup_context.tsv.

Run from the repo root:  .venv/bin/python site/build_haplogroups.py   (then deploy the site as usual)
"""
from __future__ import annotations

import csv
import re
import sys
from html import escape as e
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from gnomancer.context import slug  # noqa: E402

TSV = ROOT / "gnomancer" / "data" / "haplogroup_context.tsv"
OUT = ROOT / "site" / "haplogroups"
INV = "https://discord.gg/aGjqzqErQH"
CODE = re.compile(r"^[A-Z][0-9a-z']*$")

CSS = open(ROOT / "site" / "install" / "index.html").read().split("<style>")[1].split("</style>")[0]
HEAD = ('<header><div class="wrap"><a class="brand" href="/"><img src="/assets/mark-64.png" alt="" width="30" height="30">gnomancer</a>'
        '<nav><a href="/haplogroups/">Haplogroups</a><a href="/get-your-data/">Get your data</a><a href="/install/">Install</a>'
        f'<a href="{INV}">Community</a><a href="https://github.com/codebase/gnomancer">GitHub</a></nav></div></header>')
FOOT = (f'<footer><div class="wrap">gnomancer · MIT licensed · <a href="https://github.com/codebase/gnomancer">source</a> · '
        f'<a href="{INV}">Discord</a> · no cookies, no analytics, no uploads.</div></footer>')
EXTRA_CSS = (".crumbs{display:flex;flex-wrap:wrap;gap:6px;align-items:center;font-size:14px;color:var(--ink-2);margin:0 0 18px}"
             ".crumbs a{border:1px solid var(--line);background:var(--surface);border-radius:999px;padding:3px 10px;color:var(--ink-2)}"
             ".crumbs a:hover{text-decoration:none;border-color:var(--accent)}.crumbs .here{border-color:var(--accent);color:var(--accent-ink);font-weight:600}"
             ".grid{display:grid;gap:12px;grid-template-columns:repeat(auto-fit,minmax(260px,1fr))}.tile{background:var(--surface);border:1px solid var(--line);border-radius:var(--r);padding:14px 16px}"
             ".tile h3{margin:0 0 4px;font-size:16px}.tile p{margin:0;color:var(--ink-2);font-size:14px}.tile a{color:inherit}.tile a:hover{text-decoration:none}"
             ".kicker{font-size:13px;letter-spacing:.06em;text-transform:uppercase;color:var(--ink-3);margin:48px 0 0}"
             ".filter input{font:inherit;font-size:15px;padding:9px 12px;border:1px solid var(--line);border-radius:10px;background:var(--surface);color:var(--ink);width:min(100%,420px);margin:0 0 16px}")
REFS = {"Y": [("ISOGG Y-DNA haplogroup tree", "https://isogg.org/tree/"), ("YFull YTree", "https://www.yfull.com/tree/"),
              ("FamilyTreeDNA Discover", "https://discover.familytreedna.com/")],
        "MT": [("PhyloTree", "https://www.phylotree.org/tree/"), ("Haplogrep", "https://haplogrep.i-med.ac.at/"), ("YFull MTree", "https://www.yfull.com/mtree/")]}


def load():
    with open(TSV, newline="") as f:
        rows = [r for r in csv.DictReader([l for l in f if l.strip() and not l.startswith("#")], delimiter="\t")]
    groups = {}
    for r in rows:
        g = groups.setdefault((r["type"], r["name"]), dict(type=r["type"], name=r["name"], text=r["text"], depth=int(r["depth"]), prefixes=[]))
        g["prefixes"].append(r["prefix"])
    for g in groups.values():
        codes = [p for p in g["prefixes"] if CODE.match(p)]
        g["code"] = max(codes, key=len) if codes else g["prefixes"][0]
    return list(groups.values())


def chain(g, all_groups):
    """Ancestors of g among the same type, by code prefix, broad to narrow."""
    anc = [h for h in all_groups if h["type"] == g["type"] and h is not g and g["code"].startswith(h["code"]) and CODE.match(h["code"])]
    return sorted(anc, key=lambda h: len(h["code"]))


def children(g, all_groups):
    return sorted([h for h in all_groups if h["type"] == g["type"] and h is not g and h["code"].startswith(g["code"]) and CODE.match(h["code"])
                   and len(h["code"]) > len(g["code"])], key=lambda h: (h["depth"], h["name"]))


def page(g, all_groups):
    kind = "Y chromosome (paternal line)" if g["type"] == "Y" else "mitochondrial DNA (maternal line)"
    url = f"/haplogroups/{'y' if g['type']=='Y' else 'mt'}/{slug(g['name'])}/"
    anc = chain(g, all_groups)
    kids = [h for h in children(g, all_groups) if len(chain(h, all_groups)) == len(anc) + 1]
    title = g["name"] if g["name"].lower().startswith("haplogroup") else f"Haplogroup {g['name']}"
    crumbs = "".join(f'<a href="/haplogroups/{"y" if h["type"]=="Y" else "mt"}/{slug(h["name"])}/">{e(h["name"].replace("Haplogroup ", ""))}</a><span>›</span>' for h in anc)
    crumbs += f'<span class="here crumbs-here"><a class="here" href="{url}">{e(g["name"].replace("Haplogroup ", ""))}</a></span>'
    desc = g["text"][:155].rsplit(" ", 1)[0] + ("…" if len(g["text"]) > 155 else "")
    refs = " · ".join(f'<a href="{u}">{e(n)}</a>' for n, u in REFS[g["type"]])
    kids_html = ("<h2>Branches</h2><div class='grid'>" + "".join(
        f'<div class="tile"><h3><a href="/haplogroups/{"y" if h["type"]=="Y" else "mt"}/{slug(h["name"])}/">{e(h["name"])}</a></h3><p>{e(h["text"][:140].rsplit(" ",1)[0])}…</p></div>' for h in kids) + "</div>") if kids else ""
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{e(title)} · {kind.split(' (')[0]} · gnomancer</title>
<meta name="description" content="{e(desc)}"><link rel="canonical" href="https://gnomancer.com{url}"><link rel="icon" href="/assets/mark-64.png" type="image/png">
<meta property="og:title" content="{e(title)}"><meta property="og:description" content="{e(desc)}"><meta property="og:image" content="https://gnomancer.com/assets/og.jpg"><meta property="og:image:width" content="1200"><meta property="og:image:height" content="628"><meta name="twitter:card" content="summary_large_image"><meta name="twitter:image" content="https://gnomancer.com/assets/og.jpg">
<style>{CSS}{EXTRA_CSS}</style></head><body>{HEAD}<main class="wrap">
<p class="kicker" style="margin-top:44px">{e(kind)}</p>
<h1 style="margin-top:6px">{e(title)}</h1>
<div class="crumbs">{crumbs}</div>
<p class="lead">{e(g["text"])}</p>
<div class="card"><h3 style="margin-top:0">How gnomancer finds it</h3>
<p>{"The Y haplogroup is called from the Y-chromosome sites on your chip by yhaplo (ISOGG tree) with a second opinion from Yleaf on the current YFull, FamilyTreeDNA and ISOGG trees. Arrays type a few thousand Y sites, so the call stops at roughly this depth; Y sequencing (e.g. Big Y-700) is needed for a terminal subclade." if g["type"]=="Y" else "The mitochondrial haplogroup is called by Haplogrep3 against PhyloTree 17 from the few hundred mtDNA sites on your chip, compared with the rCRS reference. Sub-branches usually need full mitochondrial sequencing."}</p>
<p>Run it on your own raw file: <a href="/install/">install gnomancer</a> (nothing is uploaded), or ask your AI assistant with the <a href="https://github.com/codebase/gnomancer#use-it-from-your-ai-assistant">MCP server</a>.</p></div>
{kids_html}
<h2>Compare notes</h2>
<p class="lead">People with this lineage talk in the <a href="{INV}">gnomancer Discord</a>, channel #haplogroups. Share your haplogroup, not your raw file.</p>
<p class="caption" style="color:var(--ink-3)">References: {refs}. Population statements are approximate summaries of published frequencies; corrections welcome in <a href="https://github.com/codebase/gnomancer/blob/main/gnomancer/data/haplogroup_context.tsv">the source table</a>.</p>
</main>{FOOT}</body></html>'''


def index(all_groups):
    def section(kind, label):
        gs = sorted([g for g in all_groups if g["type"] == kind], key=lambda g: (g["code"], g["name"]))
        tiles = "".join(f'<div class="tile" data-k="{e(g["name"].lower())} {e(g["code"].lower())}"><h3><a href="/haplogroups/{"y" if kind=="Y" else "mt"}/{slug(g["name"])}/">{e(g["name"])}</a></h3><p>{e(g["text"][:150].rsplit(" ",1)[0])}…</p></div>' for g in gs)
        return f'<p class="kicker">{label}</p><div class="grid">{tiles}</div>'
    js = ("document.getElementById('q').oninput=function(){var v=this.value.toLowerCase();document.querySelectorAll('.tile').forEach(function(t){t.style.display=(!v||t.dataset.k.indexOf(v)>-1)?'':'none'});};")
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Y-DNA and mtDNA haplogroups explained · gnomancer</title>
<meta name="description" content="Short, population-level guides to the major Y-chromosome and mitochondrial haplogroups: where each lineage comes from, where it is common today, and how a consumer DNA test finds yours.">
<link rel="canonical" href="https://gnomancer.com/haplogroups/"><link rel="icon" href="/assets/mark-64.png" type="image/png">
<style>{CSS}{EXTRA_CSS}</style></head><body>{HEAD}<main class="wrap">
<h1>Haplogroups, explained</h1>
<p class="lead">Two small parts of your genome are inherited along a single line: the Y chromosome from father to son, and mitochondrial DNA from mother to child. Each carries a lineage label with a known geographic history. These pages give the short version for every lineage gnomancer reports, written for people who just got their result.</p>
<div class="filter"><input id="q" type="search" placeholder="Find a haplogroup, e.g. R-M269, H1, E-M84…" aria-label="Find a haplogroup"></div>
{section("Y", "Y chromosome · paternal line")}
{section("MT", "Mitochondrial DNA · maternal line")}
<p class="caption" style="margin-top:36px;color:var(--ink-3)">Population statements are approximate summaries of published frequencies. Corrections welcome in <a href="https://github.com/codebase/gnomancer/blob/main/gnomancer/data/haplogroup_context.tsv">the source table</a>.</p>
</main>{FOOT}<script>{js}</script></body></html>'''


def main():
    groups = load()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(index(groups))
    urls = ["https://gnomancer.com/haplogroups/"]
    for g in groups:
        d = OUT / ("y" if g["type"] == "Y" else "mt") / slug(g["name"])
        d.mkdir(parents=True, exist_ok=True)
        (d / "index.html").write_text(page(g, groups))
        urls.append(f"https://gnomancer.com/haplogroups/{'y' if g['type']=='Y' else 'mt'}/{slug(g['name'])}/")
    # sitemap
    sm = ROOT / "site" / "sitemap.xml"
    base = ["https://gnomancer.com/", "https://gnomancer.com/get-your-data/", "https://gnomancer.com/install/"]
    sm.write_text('<?xml version="1.0" encoding="UTF-8"?>\n<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n' +
                  "".join(f"<url><loc>{u}</loc><changefreq>monthly</changefreq></url>\n" for u in base + urls) + "</urlset>\n")
    print(f"wrote {len(groups)} pages + index; sitemap has {len(base)+len(urls)} urls")


if __name__ == "__main__":
    main()
