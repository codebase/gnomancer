"""Load any consumer raw-data file and normalize it to GRCh37.

Relies on the `snps` library (apriha/snps) for format detection. It handles
AncestryDNA, 23andMe (v2-v5), MyHeritage, FamilyTreeDNA, LivingDNA, Genes for
Good, tellmeGen, Mapmygenome, Circle DNA, plain VCF, and generic CSV/TSV.

Normalization it gives us for free:
- chromosome labels 1-22, X, Y, PAR, MT (Ancestry's 23/24/25/26 are mapped)
- haploid single-letter calls on Y and MT
- no-calls as NaN regardless of vendor code ('0', '--', 'NN', '')
- CRLF line endings, CSV quoting, header comment lines
- build detection (36/37/38) and remapping to 37 when needed
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

TARGET_BUILD = 37


@dataclass
class Loaded:
    df: pd.DataFrame            # index rsid; columns chrom, pos, genotype (NaN = no-call)
    source: str
    build: int
    build_detected: bool
    build_original: int
    remapped: bool
    chip: str
    chip_version: str
    sex: str                    # 'Male', 'Female', or ''
    count: int
    notes: list[str] = field(default_factory=list)


def load(path: str, target_build: int = TARGET_BUILD) -> Loaded:
    from snps import SNPs

    s = SNPs(path)
    if s.count == 0:
        raise SystemExit(f"Could not parse any SNPs from {path}. Is it a raw-data export?")
    notes = []
    original_build = s.build
    remapped = False
    if s.build != target_build:
        notes.append(f"Input is build {s.build}; remapping to GRCh{target_build} (downloads chain file once).")
        chromosomes_remapped, chromosomes_not_remapped = s.remap(target_build)
        remapped = True
        if chromosomes_not_remapped:
            notes.append(f"Chromosomes not remapped: {chromosomes_not_remapped}")
        if original_build == 36:
            notes.append("Build 36 mitochondrial coordinates use the old NC_001807 reference, not rCRS; "
                         "the mtDNA haplogroup call will be skipped.")
    if not s.build_detected:
        notes.append("Build was not detected from the file and was assumed; verify with the vendor.")
    df = s.snps.copy()
    df["chrom"] = df["chrom"].astype(str)
    return Loaded(
        df=df, source=s.source or "unknown", build=s.build, build_detected=s.build_detected,
        build_original=original_build, remapped=remapped, chip=s.chip or "", chip_version=s.chip_version or "",
        sex=s.determine_sex() or "", count=int(s.count), notes=notes,
    )
