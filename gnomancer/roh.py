"""Runs of homozygosity (ROH) on autosomes.

Simple scan: a run is a stretch of called SNPs with at most MAX_HET
heterozygous calls (tolerates genotyping error). Reported if it spans at
least MIN_MB megabases and MIN_SNPS SNPs. Long runs (>5 Mb) arise from
recent shared ancestry between the parents; a sum of several hundred Mb is
the territory of first-cousin parents. A few runs of 3-5 Mb are normal in
many populations, especially ones with historical founder effects.
"""
from __future__ import annotations

import pandas as pd

MIN_MB = 3.0
MIN_SNPS = 150
MAX_HET = 2


def find_roh(df: pd.DataFrame) -> pd.DataFrame:
    runs = []
    for ch in [str(i) for i in range(1, 23)]:
        sub = df[(df.chrom == ch) & df.genotype.notna()].sort_values("pos")
        pos = sub.pos.to_numpy()
        gt = sub.genotype.to_numpy()
        het = [len(g) == 2 and g[0] != g[1] for g in gt]
        run_start = 0
        het_idx: list[int] = []

        def close(start, end):
            if end > start:
                length = int(pos[end] - pos[start])
                cnt = end - start + 1
                if length >= MIN_MB * 1e6 and cnt >= MIN_SNPS:
                    runs.append(dict(chrom=ch, start=int(pos[start]), end=int(pos[end]),
                                     length_mb=round(length / 1e6, 2), snps=cnt))

        for i, h in enumerate(het):
            if h:
                het_idx.append(i)
                if len(het_idx) > MAX_HET:
                    close(run_start, het_idx[0] - 1)
                    run_start = het_idx[0] + 1
                    het_idx = het_idx[1:]
        close(run_start, len(gt) - 1)
    out = pd.DataFrame(runs, columns=["chrom", "start", "end", "length_mb", "snps"])
    return out.sort_values("length_mb", ascending=False) if not out.empty else out


def interpret(roh: pd.DataFrame) -> str:
    if roh.empty:
        return "No runs of homozygosity ≥3 Mb: no sign of recent parental relatedness."
    total = roh.length_mb.sum()
    long_ = roh[roh.length_mb >= 5]
    if total > 200 or len(long_) >= 10:
        return (f"{len(roh)} runs totalling {total:.0f} Mb, {len(long_)} longer than 5 Mb: consistent with parents "
                f"who are related (roughly second-cousin level or closer). Worth discussing with a genetic counselor if unexpected.")
    if total > 50:
        return (f"{len(roh)} runs totalling {total:.0f} Mb: mild excess homozygosity, typical of populations with "
                f"founder effects or endogamy (e.g. island, Ashkenazi, Finnish, South Asian caste groups).")
    return f"{len(roh)} short run(s) totalling {total:.0f} Mb: within the normal range for an outbred sample."
