"""Per-chromosome quality metrics and sex inference."""
from __future__ import annotations

import pandas as pd

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "PAR", "MT"]


def _is_het(g) -> bool:
    return isinstance(g, str) and len(g) == 2 and g[0] != g[1]


def chromosome_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ch in CHR_ORDER:
        sub = df[df.chrom == ch]
        if sub.empty:
            continue
        called = sub.genotype.notna()
        n_called = int(called.sum())
        het = int(sub.genotype[called].map(_is_het).sum())
        rows.append(dict(chrom=ch, snps=len(sub), called=n_called, no_call=len(sub) - n_called,
                         het=het, het_rate=(het / n_called if n_called else float("nan")),
                         first_pos=int(sub.pos.min()), last_pos=int(sub.pos.max())))
    return pd.DataFrame(rows)


def summary(df: pd.DataFrame, chr_table: pd.DataFrame, vendor_sex: str) -> dict:
    total = len(df)
    nocall = int(df.genotype.isna().sum())
    indel = int(df.genotype.fillna("").str.contains("I|D").sum())
    auto = chr_table[chr_table.chrom.isin([str(i) for i in range(1, 23)])]
    auto_het = auto.het.sum() / max(auto.called.sum(), 1)
    x = chr_table[chr_table.chrom == "X"]
    y = chr_table[chr_table.chrom == "Y"]
    x_het = float(x.het_rate.iloc[0]) if not x.empty else float("nan")
    y_called = int(y.called.iloc[0]) if not y.empty else 0
    y_total = int(y.snps.iloc[0]) if not y.empty else 0
    # Independent sex inference: X heterozygosity and Y call rate.
    if x.empty and y.empty:
        inferred = "unknown"
    elif x_het < 0.02 and (y_total == 0 or y_called / y_total > 0.5):
        inferred = "Male"
    elif x_het > 0.1:
        inferred = "Female"
    else:
        inferred = "ambiguous"
    flags = []
    if nocall / total > 0.02:
        flags.append(f"High no-call rate ({nocall/total:.1%}); sample quality may be poor.")
    if auto_het < 0.20:
        flags.append(f"Low autosomal heterozygosity ({auto_het:.3f}); possible inbreeding, contamination-free but "
                     f"homozygous-rich sample, or array with many monomorphic probes.")
    if auto_het > 0.40:
        flags.append(f"High autosomal heterozygosity ({auto_het:.3f}); possible sample contamination or mixture.")
    if vendor_sex and inferred not in ("unknown", "ambiguous") and vendor_sex != inferred:
        flags.append(f"Sex inferred from X/Y ({inferred}) disagrees with library inference ({vendor_sex}).")
    return dict(total=total, no_call=nocall, no_call_rate=nocall / total, indel_calls=indel,
                autosomal_het=auto_het, x_het_rate=x_het, y_called=y_called, y_total=y_total,
                sex_inferred=inferred, sex_library=vendor_sex, flags=flags)
