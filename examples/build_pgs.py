"""Rebuild gnomancer/data/pgs/: a few published polygenic scores, with their distribution in 1000 Genomes.

A polygenic score is a weighted count of risk alleles across many sites. The number means nothing by itself:
it has to be placed in a reference distribution. The usual way is to download a 7-15 GB reference panel and
score it alongside the sample. gnomancer does that work once, here, and ships only the outcome, so scoring a
person needs no download:

  data/pgs/<ID>.tsv    the score's variants: chrom, pos (GRCh38), ref, alt, effect allele, weight, and the
                       effect-allele frequency in each 1000 Genomes super-population
  data/pgs/index.json  what each score is, where it was published, how many of its variants could be used, and
                       the sorted scores of the unrelated 1000 Genomes samples in each super-population

Only scores small enough to be transparent (tens to hundreds of variants, mostly genome-wide-significant hits)
are used. Genome-wide scores with millions of variants predict better but need the full panel and LD-aware
handling; they are out of scope by design.

How a reference sample is scored here is exactly how gnomancer scores a person: effect-allele dosage at each
variant, with the reference allele assumed where no variant is recorded. That keeps the two comparable.

Sources: PGS Catalog harmonized scoring files (pgscatalog.org, per-score citation in index.json) and the 1000
Genomes 30x phased panel (Byrska-Bishop et al., Cell 2022). Needs network and tabix (`gnomancer setup --wgs`).

Usage: python examples/build_pgs.py [PGS000004 ...]
"""
from __future__ import annotations

import gzip
import io
import json
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from gnomancer import wgstools  # noqa: E402

OUT = ROOT / "gnomancer" / "data" / "pgs"
KG = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
      "1kGP_high_coverage_Illumina.chr{}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
KG_SAMPLES = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
POPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
MIN_MATCHED = 0.85      # a score is bundled only if at least this share of its published variants is in the panel

# id -> (label shown in the report, sex the score applies to or None, one line on what a high score means)
SCORES = {
    "PGS000004": ("Breast cancer", "Female", "313 variants; in the original study the top 1% of women had roughly four times the risk of ER-positive disease of women in the middle fifth"),
    "PGS000662": ("Prostate cancer", "Male", "269 variants from a multi-ancestry study"),
    "PGS003438": ("Coronary artery disease", None, "241 variants"),
    "PGS020005": ("Atrial fibrillation", None, "166 variants; one of the scores the eMERGE network returns to patients"),
    "PGS005111": ("Type 2 diabetes", None, "150 variants"),
    "PGS000043": ("Venous thromboembolism", None, "297 variants; adds to, and does not replace, factor V Leiden and prothrombin status"),
}


def curl(url: str, binary: bool = False):
    p = subprocess.run(["curl", "-fsSL", "--retry", "3", url], capture_output=True)
    if p.returncode != 0:
        raise RuntimeError(f"download failed: {url}: {p.stderr.decode()[-200:]}")
    return p.stdout if binary else p.stdout.decode()


def scoring_file(pid: str) -> tuple[dict, list]:
    """Metadata and the usable variants of a score on GRCh38: [(chrom, pos, effect, other, weight)]."""
    meta = json.loads(curl(f"https://www.pgscatalog.org/rest/score/{pid}"))
    raw = curl(meta["ftp_harmonized_scoring_files"]["GRCh38"]["positions"], binary=True)
    lines = [l for l in gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode().splitlines() if l and not l.startswith("#")]
    head = lines[0].split("\t")
    col = {c: i for i, c in enumerate(head)}
    out, skipped = [], 0
    for line in lines[1:]:
        c = line.split("\t")
        get = lambda k: c[col[k]] if k in col and col[k] < len(c) else ""  # noqa: E731
        chrom, pos, effect = get("hm_chr"), get("hm_pos"), get("effect_allele").upper()
        other = (get("other_allele") or get("hm_inferOtherAllele")).upper()
        # Autosomal SNVs and simple indels with both alleles spelled out; X needs sex-aware dosage and is left out.
        if chrom not in map(str, range(1, 23)) or not pos.isdigit() or not effect or not other or set(effect + other) - set("ACGT") or effect == other:
            skipped += 1
            continue
        out.append((chrom, int(pos), effect, other, float(get("effect_weight"))))
    meta["_skipped_at_parse"] = skipped
    return meta, out


def samples() -> dict:
    """sample -> super-population, unrelated samples only (both parents unknown)."""
    out = {}
    for line in curl(KG_SAMPLES).splitlines()[1:]:
        f = line.split()
        if len(f) >= 7 and f[2] == "0" and f[3] == "0" and f[6] in POPS:
            out[f[1]] = f[6]
    return out


def panel_genotypes(variants: list, tabix: str) -> dict:
    """(chrom, pos) -> (ref, alt, [alt-allele count per sample]) and the sample order, from the 1000 Genomes panel."""
    by_chrom: dict = {}
    for chrom, pos, *_ in variants:
        by_chrom.setdefault(chrom, set()).add(pos)
    want = {(c, p): {e, o} for c, p, e, o, _ in variants}
    cache = Path(tempfile.gettempdir()) / "gnomancer-1kg-index"
    cache.mkdir(exist_ok=True)

    def fetch(chrom: str):
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as f:
            for p in sorted(by_chrom[chrom]):
                f.write(f"chr{chrom}\t{p - 1}\t{p}\n")
        r = subprocess.run([tabix, "-h", "-R", f.name, KG.format(chrom)], capture_output=True, text=True, cwd=str(cache))
        Path(f.name).unlink()
        if r.returncode != 0:
            raise RuntimeError(f"tabix failed on chr{chrom}: {r.stderr[-200:]}")
        got, names = {}, None
        for line in r.stdout.splitlines():
            if line.startswith("##"):
                continue
            c = line.split("\t")
            if line.startswith("#CHROM"):
                names = c[9:]
                continue
            k = (chrom, int(c[1]))
            if k not in want or {c[3], c[4]} != want[k]:
                continue                                    # another allele at the same position, or a different way of writing an indel
            got[k] = (c[3], c[4], [g[0:3].count("1") for g in c[9:]])
        print(f"  chr{chrom}: {len(got)}/{len(by_chrom[chrom])}", flush=True)
        return got, names

    from concurrent.futures import ThreadPoolExecutor
    out, names = {}, None
    with ThreadPoolExecutor(6) as pool:
        for got, n in pool.map(fetch, sorted(by_chrom, key=int)):
            out.update(got)
            names = names or n
    return out, names


def build(pid: str, pop_of: dict, tabix: str) -> dict:
    label, sex, note = SCORES.get(pid, (None, None, ""))
    meta, variants = scoring_file(pid)
    print(f"{pid} {meta['trait_reported']}: {len(variants)} usable of {meta['variants_number']} variants")
    geno, names = panel_genotypes(variants, tabix)
    keep = [i for i, n in enumerate(names) if n in pop_of]
    pops = [pop_of[names[i]] for i in keep]
    scores = [0.0] * len(keep)
    rows = []
    for chrom, pos, effect, other, weight in variants:
        g = geno.get((chrom, pos))
        if not g:
            continue
        ref, alt, counts = g
        dos = [counts[i] if effect == alt else 2 - counts[i] for i in keep]
        for j, d in enumerate(dos):
            scores[j] += weight * d
        freq = {p: sum(d for d, q in zip(dos, pops) if q == p) / (2 * pops.count(p)) for p in POPS}
        rows.append((chrom, pos, ref, alt, effect, weight, freq))
    rows.sort(key=lambda r: (int(r[0]), r[1]))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"{pid}.tsv", "w") as f:
        f.write(f"# {pid}: {meta['trait_reported']}. Variants found in the 1000 Genomes 30x panel, GRCh38. Generated by examples/build_pgs.py.\n")
        f.write("chrom\tpos\tref\talt\teffect\tweight\t" + "\t".join(f"af_{p}" for p in POPS) + "\n")
        for chrom, pos, ref, alt, effect, weight, freq in rows:
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{effect}\t{weight:.6g}\t" + "\t".join(f"{freq[p]:.4f}" for p in POPS) + "\n")
    pub = meta["publication"]
    dist = {p: sorted(round(s, 4) for s, q in zip(scores, pops) if q == p) for p in POPS}
    return dict(id=pid, label=label or meta["trait_reported"], trait=meta["trait_reported"], sex=sex, note=note, name=meta["name"],
                variants_published=meta["variants_number"], variants_used=len(rows),
                citation=f"{pub['firstauthor']} et al., {pub['journal']} {pub['date_publication'][:4]}", doi=pub.get("doi"), pmid=pub.get("PMID"),
                licence=meta.get("license"), development_ancestry=(meta.get("ancestry_distribution") or {}).get("gwas", {}).get("dist"),
                reference={p: dict(n=len(dist[p]), mean=round(statistics.mean(dist[p]), 4), sd=round(statistics.stdev(dist[p]), 4), scores=dist[p]) for p in POPS})


def main() -> None:
    tabix = wgstools.tool_bin("tabix")
    if not tabix:
        raise SystemExit("needs tabix; run `gnomancer setup --wgs`")
    ids = sys.argv[1:] or list(SCORES)
    pop_of = samples()
    print(f"{len(pop_of)} unrelated 1000 Genomes samples")
    index_path = OUT / "index.json"
    index = json.loads(index_path.read_text()) if index_path.exists() else {}
    for pid in ids:
        entry = build(pid, pop_of, tabix)
        if entry["variants_used"] / entry["variants_published"] < MIN_MATCHED:
            print(f"  {pid}: only {entry['variants_used']} of {entry['variants_published']} variants matched the panel; not bundled")
            (OUT / f"{pid}.tsv").unlink(missing_ok=True)
            index.pop(pid, None)
            continue
        index[pid] = entry
        OUT.mkdir(parents=True, exist_ok=True)
        index_path.write_text(json.dumps(dict(sorted(index.items())), indent=1) + "\n")
    print("done:", ", ".join(f"{k} {v['variants_used']}/{v['variants_published']}" for k, v in index.items()))


if __name__ == "__main__":
    main()
