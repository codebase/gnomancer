"""Quality summary and runs of homozygosity for a sequencing VCF.

Array QC asks 'how many probes failed'. For sequencing the useful questions are different: is this a whole
genome or a subset, do the variant counts and ratios look like a human genome, and how much did the quality
gates remove. Typical values for one 30x genome (GRCh38): 4-5 million variants, transition/transversion ratio
about 2.0-2.1, 1.5-2 heterozygous SNVs for each homozygous one.
"""
from __future__ import annotations

import pandas as pd

from .seqvcf import AUTOSOMES, CHROMS, SeqData

WHOLE_GENOME_MIN_SNVS = 2_000_000
BIN = 100_000
MIN_MB = 3.0


def is_whole_genome(d: SeqData) -> bool:
    """A whole genome has millions of SNVs; an exome has tens of thousands. Decides whether absent means reference."""
    return d.stats.get("snvs", 0) >= WHOLE_GENOME_MIN_SNVS


def chromosome_table(d: SeqData) -> pd.DataFrame:
    rows = []
    for ch in CHROMS:
        st = d.per_chrom.get(ch) or {}
        het, hom = st.get("snv_het", 0), st.get("snv_hom", 0)
        n = het + hom
        if not n and not st.get("failed"):
            continue
        rows.append(dict(chrom=ch, snps=n, called=n, no_call=st.get("failed", 0), het=het,
                         het_rate=(het / n if n else float("nan")), indels=st.get("indel_het", 0) + st.get("indel_hom", 0)))
    return pd.DataFrame(rows)


def summary(d: SeqData) -> dict:
    s = d.stats
    whole = is_whole_genome(d)
    total = s["passing_alleles"] + s["failed_alleles"]
    failed_rate = s["failed_alleles"] / total if total else 0.0
    x = s["x_het"] + s["x_hom"]
    flags = []
    if not whole:
        flags.append(f"Only {s['snvs']:,} SNVs: this looks like an exome, a gene panel or a subset, not a whole genome. "
                     "Sites missing from the file are treated as unknown, not as reference.")
    if whole and s["ti_tv"] is not None and not 1.85 <= s["ti_tv"] <= 2.3:
        flags.append(f"Transition/transversion ratio {s['ti_tv']:.2f} is outside the 1.9-2.2 expected for a whole genome; "
                     "a low value points to many false-positive calls.")
    r = s["het_hom_ratio"]
    if whole and r is not None and r > 3.0:
        flags.append(f"Heterozygous/homozygous ratio {r:.2f} is high; possible sample contamination or mixture.")
    if whole and r is not None and r < 1.0:
        flags.append(f"Heterozygous/homozygous ratio {r:.2f} is low: long homozygous stretches, see the runs-of-homozygosity section.")
    if failed_rate > 0.15:
        flags.append(f"{failed_rate:.0%} of variant calls failed the quality gates (caller filter, depth, genotype quality, allele balance).")
    md = s["mean_depth_at_variants"]
    if md is not None and md < 15:
        flags.append(f"Mean depth at variant sites is {md:.0f}x; below about 15x heterozygous variants are often missed.")
    if d.sex == "ambiguous":
        flags.append("Sex could not be inferred cleanly from X heterozygosity.")
    if d.sniff.mt_length and d.sniff.mt_length != 16569:
        flags.append(f"The mitochondrial contig is {d.sniff.mt_length} bases long, not the 16,569 of rCRS; mitochondrial positions cannot be compared and the maternal line is skipped.")
    return dict(kind="wgs", whole_genome=whole, gvcf=d.gvcf, caller=d.sniff.caller, total=total, no_call=s["failed_alleles"],
                no_call_rate=failed_rate, indel_calls=s["indels"], snvs=s["snvs"], ti_tv=s["ti_tv"], het_hom_ratio=r,
                mean_depth=md, autosomal_het=(r / (1 + r) if r is not None else float("nan")),
                x_het_rate=(s["x_het"] / x if x else float("nan")), y_called=s["y_calls"], y_total=s["y_calls"],
                sex_inferred=d.sex, sex_library="", flags=flags, records=s["records"], reference_blocks=s["reference_blocks"],
                symbolic_alleles=s["symbolic_alleles"], off_primary=s["off_primary_records"])


def find_roh(d: SeqData) -> pd.DataFrame:
    """Runs of homozygosity from heterozygous-SNV density.

    The genome is cut into 100 kb bins. Outside a run a bin holds dozens of heterozygous SNVs; inside one it holds
    almost none while homozygous SNVs continue. A bin is 'homozygous' when it has at least 5 homozygous SNVs and
    heterozygous ones make up under 10% (stray heterozygous calls come from mismapped reads). Bins with no data
    at all (centromeres, gaps) neither start nor break a run, and one stray bin inside a run is tolerated.
    """
    runs = []
    for ch in sorted(AUTOSOMES, key=int):
        het, hom = d.het_pos.get(ch), d.hom_pos.get(ch)
        if not hom:
            continue
        nb = max(max(hom), max(het) if het else 0) // BIN + 1
        h, o = [0] * nb, [0] * nb
        for p in het or ():
            h[p // BIN] += 1
        for p in hom:
            o[p // BIN] += 1
        state = []                       # 1 homozygous, 0 not, None no data
        for a, b in zip(h, o):
            state.append(None if (a == 0 and b < 5) else 1 if (b >= 5 and a / (a + b) < 0.10) else 0)
        i = 0
        while i < nb:
            if state[i] != 1:
                i += 1
                continue
            j = last = i
            strays = 0
            while j + 1 < nb:
                nxt = state[j + 1]
                if nxt == 1:
                    last = j = j + 1
                    strays = 0
                elif nxt is None:
                    j += 1
                elif strays == 0 and j + 2 < nb and state[j + 2] == 1:
                    strays = 1
                    j += 1
                else:
                    break
            good = sum(1 for s_ in state[i:last + 1] if s_ == 1)
            length = (last - i + 1) * BIN
            if length >= MIN_MB * 1e6 and good / (last - i + 1) >= 0.6:
                runs.append(dict(chrom=ch, start=i * BIN + 1, end=(last + 1) * BIN, length_mb=round(length / 1e6, 2), snps=sum(o[i:last + 1])))
            i = last + 1
    out = pd.DataFrame(runs, columns=["chrom", "start", "end", "length_mb", "snps"])
    return out.sort_values("length_mb", ascending=False) if not out.empty else out
