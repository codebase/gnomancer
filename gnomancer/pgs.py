"""Polygenic scores, placed against a reference distribution that ships with gnomancer.

A polygenic score adds up many common variants, each with a small published weight. The sum has no meaning until
it is compared with the same sum in other people, so each bundled score carries the scores of the unrelated 1000
Genomes samples in five super-populations (built once by examples/build_pgs.py). A person's score is computed
the same way those were, and reported as a percentile within the population the sample sits closest to.

What this is and is not:
- Only small, published scores are bundled (tens to hundreds of variants). Genome-wide scores with millions of
  variants predict better and need a multi-gigabyte reference panel; they are out of scope.
- A percentile is a position in a distribution, not a probability of disease. Most people in the top tenth of a
  score never develop the condition, and a low score does not protect against a rare high-impact variant,
  family history or environment.
- Scores are developed mostly in European-ancestry studies and transfer imperfectly. The percentile is computed
  within the matching 1000 Genomes super-population, which corrects the shift in the distribution but not the
  loss of predictive power.
- GRCh38 whole genomes only: the score needs a genotype at every variant, including the reference ones.
"""
from __future__ import annotations

import bisect
import csv
import json

from .paths import DATA
from .seqsites import Site, genotype

PGS_DIR = DATA / "pgs"
POPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
POP_NAME = dict(AFR="African", AMR="Admixed American", EAS="East Asian", EUR="European", SAS="South Asian")
MAX_MISSING = 0.10


def load_index() -> dict:
    p = PGS_DIR / "index.json"
    return json.loads(p.read_text()) if p.exists() else {}


def load_variants(pid: str) -> list[dict]:
    with open(PGS_DIR / f"{pid}.tsv", newline="") as f:
        rows = list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))
    for r in rows:
        r["pos"], r["weight"] = int(r["pos"]), float(r["weight"])
        for p in POPS:
            r[f"af_{p}"] = float(r[f"af_{p}"])
    return rows


def positions() -> dict:
    """chrom -> sorted positions of every bundled score variant (GRCh38), for the main VCF pass to track."""
    out: dict = {}
    for pid in load_index():
        for r in load_variants(pid):
            out.setdefault(r["chrom"], set()).add(r["pos"])
    return {c: sorted(p) for c, p in out.items()}


def score_one(d, pid: str, meta: dict, whole_genome: bool, pop: str) -> dict:
    variants = load_variants(pid)
    total, missing, missing_weight, weight_sum = 0.0, 0, 0.0, 0.0
    for v in variants:
        site = Site(f"{pid}:{v['chrom']}:{v['pos']}", "pgs", v["chrom"], v["pos"], v["ref"], (v["alt"],))
        g, _ = genotype(d, site, whole_genome)
        weight_sum += abs(v["weight"])
        effect, allowed = v["effect"], {v["ref"], v["alt"]}
        if site.indel:                       # indel genotypes come back as I (the longer allele) and D
            other = v["alt"] if effect == v["ref"] else v["ref"]
            effect, allowed = ("I" if len(effect) > len(other) else "D"), {"I", "D"}
        if g is None or len(g) != 2 or not set(g) <= allowed:
            # No confident genotype: use the population's expected dosage, the usual treatment, and count it.
            missing += 1
            missing_weight += abs(v["weight"])
            total += v["weight"] * 2 * v[f"af_{pop}"]
            continue
        total += v["weight"] * g.count(effect)
    ref = meta["reference"][pop]
    res = dict(id=pid, label=meta["label"], citation=meta["citation"], note=meta.get("note", ""), variants=len(variants),
               variants_published=meta["variants_published"], missing=missing, score=round(total, 4), population=pop, usable=True, reason="")
    if weight_sum and missing_weight / weight_sum > MAX_MISSING:
        res.update(usable=False, reason=f"{missing} of {len(variants)} variants had no confident genotype, too many to place the score")
        return res
    scores = ref["scores"]
    res["percentile"] = round(100 * (bisect.bisect_left(scores, total) + bisect.bisect_right(scores, total)) / (2 * len(scores)), 1)
    res["z"] = round((total - ref["mean"]) / ref["sd"], 2) if ref["sd"] else 0.0
    res["reference_n"] = ref["n"]
    res["all_populations"] = {p: round(100 * bisect.bisect_left(meta["reference"][p]["scores"], total) / len(meta["reference"][p]["scores"]), 1) for p in POPS}
    return res


def ordinal(pct: float) -> str:
    """12.4 -> '12th', 22 -> '22nd', 99.6 -> '99th' (never '100th': nobody is above everyone)."""
    n = min(max(int(round(pct)), 1), 99)
    suffix = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def band(pct: float) -> str:
    if pct >= 95:
        return "in the top 5%"
    if pct >= 80:
        return "higher than about four in five people"
    if pct > 20:
        return "in the middle range"
    if pct > 5:
        return "lower than about four in five people"
    return "in the bottom 5%"


def score(d, whole_genome: bool, sex: str, pop: str | None) -> dict:
    """Score every bundled PGS that applies to this sample. `pop` is a 1000 Genomes super-population code, or None."""
    index = load_index()
    assumed = pop not in POPS
    pop = pop if pop in POPS else "EUR"
    if not index:
        return dict(available=False, error="no polygenic scores are bundled", results=[], section=None)
    if d.build != 38:
        return dict(available=False, error="polygenic scores are computed for GRCh38 genomes only", results=[], section=None)
    if not whole_genome:
        return dict(available=False, error="polygenic scores need a whole genome: every variant of the score must have a genotype", results=[], section=None)
    results = [score_one(d, pid, meta, whole_genome, pop) for pid, meta in index.items() if not meta.get("sex") or meta["sex"] == sex]
    return dict(available=True, error=None, results=results, population=pop, population_assumed=assumed, section=section(results, pop, assumed))


def section(results: list, pop: str, assumed: bool) -> dict:
    rows = []
    for r in results:
        if not r["usable"]:
            rows.append([r["label"], "not computed", r["reason"], f"{r['variants'] - r['missing']} of {r['variants_published']}", r["citation"]])
            continue
        rows.append([r["label"], f"{ordinal(r['percentile'])} percentile", band(r["percentile"]) + (f". {r['note']}" if r["note"] else ""),
                     f"{r['variants'] - r['missing']} of {r['variants_published']}", f"{r['citation']} ({r['id']})"])
    who = POP_NAME[pop]
    text = (f"Each score is compared with the unrelated {who} samples of the 1000 Genomes Project, "
            + ("because no ancestry estimate was available; if that is not your background the percentiles are off, so run with --ancestry. "
               if assumed else "the reference group the ancestry estimate placed this sample in. ")
            + "A percentile says where the score falls among those people, not how likely the condition is: polygenic scores capture a small part of risk, "
              "most people with a high score stay well, and a low score offers no protection against a rare high-impact variant (the ClinVar section), "
              "family history or lifestyle. The scores were developed mainly in people of European ancestry and predict less well elsewhere.")
    return dict(id="pgs", title="Polygenic scores", nav="Polygenic scores", sensitive=True,
                lead="Common variants that each nudge risk slightly, added up with published weights. Only small, published, openly licensed scores are used, "
                     "so each one is traceable to its paper.",
                text=text, columns=["Condition", "Where you fall", "Reading", "Variants used", "Source"], rows=rows,
                notes=["Method: effect-allele dosage times weight, summed; positions with no recorded variant count as reference; the reference distribution was "
                       "computed the same way on the 1000 Genomes 30x panel. Score files from the PGS Catalog (pgscatalog.org)."])


def bullets(res: dict) -> tuple[list, list]:
    if not res.get("available"):
        return [], []
    high = [r for r in res["results"] if r["usable"] and r["percentile"] >= 90]
    notable = [f"Polygenic score for {r['label'].lower()}: {ordinal(r['percentile'])} percentile among {POP_NAME[res['population']]} reference samples. "
               "A modest shift in odds, not a prediction; worth mentioning when screening is discussed." for r in high]
    rest = [r for r in res["results"] if r["usable"] and r["percentile"] < 90]
    good = ["Polygenic scores in the unremarkable range: " + ", ".join(f"{r['label'].lower()} ({ordinal(r['percentile'])} percentile)" for r in rest) + "."] if rest else []
    return notable, good
