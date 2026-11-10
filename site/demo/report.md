> SYNTHETIC DEMO — No real person. Artificial genotypes and positions. Haplogroups, ClinVar and PharmCAT were skipped. Not clinical data.

# Raw DNA analysis report

Computed locally from a consumer genotyping array. The array reports a fixed list of common sites, not your whole genome; absence of a variant here never means a gene is normal, and single-marker readings are probabilistic. This is not medical advice. Discuss anything concerning with a clinician or genetic counselor and confirm with a clinical-grade test before acting on it.

## At a glance

**Who the data says you are**
- Blood group O (Rh not typed).
- Lactose: non-persistent.
- Bitter taste: PAV/AVI or mixed: medium taster.

**Worth knowing**
- G6PD A- (V68M): carrier / affected if male.
- NOD2 G908R (Crohn's): carrier.
- LPA High Lp(a) tag: elevated Lp(a) likely (~2x CAD odds); a direct Lp(a) blood test is cheap.
- GBA N370S (Gaucher / Parkinson's): carrier: ~5x Parkinson's odds; Gaucher carrier.
- One copy of a risk allele at: HLA-DQB1 (DQ8 tag (celiac)), TCF7L2 (Type 2 diabetes), CDKN2B-AS1 (9p21) (Coronary artery disease). Population-average territory, listed for completeness.
- CYP2C9 *2 (R144C): one *2: reduced warfarin/NSAID/phenytoin metabolism.
- CYP2C9 *3 (I359L): one *3: markedly reduced metabolism.
- VKORC1 1173C>T (proxy for -1639): sensitive.
- ABCG2 Q141K: homozygous reduced.
- HCP5 HLA-B*57:01 tag: HLA-B*57:01 carrier: abacavir hypersensitivity, avoid.

**Reassuring**
- CLU: two protective alleles.
- APOE e3/e3: most common genotype, reference risk for Alzheimer's and lipids.
- Negative at 28 checked health and carrier sites, including F5, F2, HFE, MTHFR, HLA-DQA1, HLA-DRB1, PTPN22, CFH.
- No runs of homozygosity ≥3 Mb.
- The file itself passed every quality check.

**Sensible next steps**
- Keep the pharmacogenomic table handy for any prescriber (the one-page card next to this report, pgx-card.pdf, is made for that); a clinical PGx panel confirms it and covers the sites the chip misses.

Sections: input · quality · homozygosity · paternal line · maternal line · ancestry · combined readings · marker panel · ClinVar · drug-gene results · limits

## Input

- Vendor/format detected: **AncestryDNA**
- Genome build: GRCh37
- Sites loaded: 20,342

## Quality

| Metric | Value |
|---|---|
| No-calls | 0 (0.00%) |
| Indel calls (I/D) | 6 |
| Autosomal heterozygosity | 0.326 |
| X heterozygosity | 0.2483 |
| Y sites called | 0 / 0 |
| Genetic sex | Female (file header: Female) |

Per-chromosome table: `chromosome-qc.tsv`.

## Runs of homozygosity

No runs of homozygosity ≥3 Mb: no sign of recent parental relatedness.

## Paternal line (Y chromosome)

Not called: skipped (sample is not male)

## Maternal line (mitochondrial DNA)

Not called: skipped (--no-mt)

## Ancestry estimate (ezancestry)

Not run. Add `--ancestry` for a 1000 Genomes population classification (`gnomancer setup` installs ezancestry into a side venv when Python 3.11 is available).

## Combined readings

- **APOE:** e3/e3 (rs429358 TT, rs7412 CC) — most common genotype, reference risk for Alzheimer's and lipids
- **ABO blood group:** O (rs8176719 DD, rs8176746 TG) — Rh factor is not typed on consumer arrays; A1 vs A2 subgroups not distinguished
- **Lactase persistence:** non-persistent (likely adult lactose intolerance) — Based on the European -13910 and -22018 variants; African and Middle Eastern persistence alleles (rs41380347, rs145946881, rs41525747) are usually not on consumer arrays
- **TAS2R38 bitter taste (PTC/PROP):** PAV/AVI or mixed: medium taster — haplotype score 5/6 across the three defining SNPs

## Marker panel

Effect-allele counts are on the GRCh37 forward strand. `not on array` means the vendor's chip does not type that site. `share` is the fraction of the European 1000 Genomes reference population expected to carry the same genotype (Hardy-Weinberg from the effect-allele frequency); the TSV has all five super-populations.


### Ancestry-informative markers

| rsid | gene | label | genotype | effect_count | share | interpretation |
|---|---|---|---|---|---|---|
| rs1426654 | SLC24A5 | Skin pigmentation A111T | AG | 1/2 | 1% | one light-skin allele |
| rs16891982 | SLC45A2 | Skin/hair pigmentation L374F | GC | 1/2 | 12% | one European light-pigmentation allele |
| rs3827760 | EDAR | V370A hair thickness / shovel incisors | GA | 1/2 | 2% | one copy of the East Asian/Native American allele |
| rs2814778 | ACKR1 (DARC) | Duffy null (FY*O) | CT | 1/2 | 1% | one Duffy-null allele (common in West/Central African ancestry) |
| rs4988235 | LCT | Lactase persistence -13910 | GG | 0/2 | 24% | non-persistent at this site |
| rs182549 | LCT | Lactase persistence -22018 | CC | 0/2 | 24% | non-persistent at this site |
| rs671 | ALDH2 | ALDH2*2 alcohol flush | GA | 1/2 | 0% | flush reaction, reduced acetaldehyde clearance (common in East Asia) |
| rs1229984 | ADH1B | ADH1B*2 (His48) | CC | 0/2 | 94% | *1/*1 typical of Europe and Africa |

### Traits

| rsid | gene | label | genotype | effect_count | share | interpretation |
|---|---|---|---|---|---|---|
| rs12913832 | HERC2 | Eye color (main determinant) | GA | 1/2 | 46% | brown or hazel likely, green possible |
| rs1800407 | OCA2 | Eye color modifier R419Q | TC | 1/2 | 14% | green/hazel modifier present |
| rs12896399 | SLC24A4 | Eye/hair color | GG | 0/2 | 32% | darker-associated |
| rs1042602 | TYR | S192Y pigmentation / freckling | CC | 0/2 | 39% | no lighter allele |
| rs1393350 | TYR | Eye/hair color | AG | 1/2 | 37% | intermediate |
| rs12203592 | IRF4 | Freckles, lighter hair, sun sensitivity | TT | 2/2 | 1% | freckles and sun sensitivity likely |
| rs1805007 | MC1R | R151C red hair | CT | 1/2 | 13% | red-hair carrier (R allele) |
| rs1805008 | MC1R | R160W red hair | CT | 1/2 | 12% | red-hair carrier (R allele) |
| rs1805009 | MC1R | D294H red hair | GG | 0/2 | 98% | no variant |
| rs1805006 | MC1R | D84E red hair | CC | 0/2 | 98% | no variant |
| rs11547464 | MC1R | R142H red hair | GG | 0/2 | 98% | no variant |
| rs1110400 | MC1R | I155T red hair | TT | 0/2 | 98% | no variant |
| rs885479 | MC1R | R163Q (weak variant) | GA | 1/2 | 13% | one weak r allele |
| rs11803731 | TCHH | Hair curl (European) | TT | 0/2 | 5% | curlier-associated |
| rs17822931 | ABCC11 | Earwax type / underarm odor | CC | 0/2 | 75% | wet earwax, typical body odor |
| rs4481887 | OR2M7 region | Asparagus metabolite smell | GG | 0/2 | 52% | less likely to smell it |
| rs10427255 | ZEB2 region | Photic sneeze reflex | TC | 1/2 | 50% | intermediate |
| rs1815739 | ACTN3 | R577X sprint/power muscle | CC | 0/2 | 32% | RR: two working copies, power-leaning |
| rs762551 | CYP1A2 | Caffeine metabolism *1F | CA | 1/2 | 44% | intermediate |
| rs2472297 | CYP1A1/2 region | Caffeine intake tendency | CC | 0/2 | 62% | lower-intake alleles |
| rs4410790 | AHR | Caffeine intake tendency | TT | 0/2 | 15% | lower-intake alleles |
| rs713598 | TAS2R38 | Bitter taste PTC (A49P) | GG | 2/2 |  | taster |
| rs1726866 | TAS2R38 | Bitter taste (V262A) | GA | 1/2 |  | taster |
| rs10246939 | TAS2R38 | Bitter taste (I296V) | CC | 2/2 |  | taster |
| rs72921001 | OR6A2 | Cilantro soapy taste | AC | 1/2 | 46% | one soapy-taste allele |
| rs9939609 | FTO | BMI/obesity risk | TT | 0/2 | 34% | no risk allele |
| rs17782313 | MC4R | BMI | TT | 0/2 | 58% | no risk allele |
| rs16969968 | CHRNA5 | Nicotine dependence D398N | GG | 0/2 | 40% | no risk allele |
| rs4680 | COMT | Val158Met | GA | 1/2 | 50% | Val/Met |
| rs6265 | BDNF | Val66Met | CT | 1/2 | 32% | Val/Met |
| rs1799971 | OPRM1 | A118G | GA | 1/2 | 27% | one G allele |
| rs601338 | FUT2 | Secretor status W154X | GG | 0/2 | 31% | secretor |
| rs2180439 | PAX1 region | Male pattern hair loss | TC | 1/2 | 49% | one risk allele |
| rs53576 | OXTR | Oxytocin receptor | GG | 0/2 | 42% | no variant allele |
| rs1800497 | ANKK1/DRD2 | TaqIA | GG | 0/2 | 66% | A2/A2 |
| rs8176719 | ABO | O allele (261delG) | DD | 2/2 | 37% | blood group O |
| rs8176746 | ABO | B allele | TG | 1/2 |  | one B allele |

### Health-associated common variants

| rsid | gene | label | genotype | effect_count | share | interpretation |
|---|---|---|---|---|---|---|
| rs1050828 | G6PD | G6PD A- (V68M) | TC | 1/2 |  | carrier / affected if male |
| rs429358 | APOE | e4-defining site (C112R) | TT | 0/2 | 71% | no e4 allele |
| rs7412 | APOE | e2-defining site (R158C) | CC | 0/2 | 88% | no e2 allele |
| rs6025 | F5 | Factor V Leiden R506Q | CC | 0/2 | 98% | negative |
| rs1799963 | F2 | Prothrombin G20210A | GG | 0/2 | 98% | negative |
| rs1800562 | HFE | Hemochromatosis C282Y | GG | 0/2 | 92% | negative |
| rs1799945 | HFE | Hemochromatosis H63D | CC | 0/2 | 69% | negative |
| rs1800730 | HFE | S65C | AA | 0/2 | 97% | negative |
| rs1801133 | MTHFR | C677T (A222V) | GG | 0/2 | 40% | normal activity |
| rs1801131 | MTHFR | A1298C | TT | 0/2 | 47% | normal |
| rs2187668 | HLA-DQA1 | DQ2.5 tag (celiac) | CC | 0/2 | 80% | no DQ2.5 |
| rs7454108 | HLA-DQB1 | DQ8 tag (celiac) | TC | 1/2 |  | DQ8 present |
| rs3135388 | HLA-DRB1 | DRB1*15:01 tag (multiple sclerosis) | GG | 0/2 | 77% | non-carrier |
| rs2476601 | PTPN22 | R620W autoimmunity | GG | 0/2 | 82% | negative |
| rs2066844 | NOD2 | R702W (Crohn's) | CC | 0/2 | 90% | negative |
| rs2066845 | NOD2 | G908R (Crohn's) | GC | 1/2 | 2% | carrier |
| rs7903146 | TCF7L2 | Type 2 diabetes | CT | 1/2 | 43% | one risk allele (~1.4x odds) |
| rs1061170 | CFH | Y402H macular degeneration | TT | 0/2 | 41% | no risk allele |
| rs10490924 | ARMS2 | A69S macular degeneration | GG | 0/2 | 65% | no risk allele |
| rs1333049 | CDKN2B-AS1 (9p21) | Coronary artery disease | GC | 1/2 | 50% | one copy (~1.25x) |
| rs10455872 | LPA | High Lp(a) tag | AG | 1/2 | 14% | elevated Lp(a) likely (~2x CAD odds); a direct Lp(a) blood test is cheap |
| rs3798220 | LPA | I4399M high Lp(a) | TT | 0/2 | 98% | no high-Lp(a) allele |
| rs11136000 | CLU | Alzheimer's (protective T) | TT | 2/2 | 15% | two protective alleles |
| rs75932628 | TREM2 | R47H Alzheimer's | CC | 0/2 | 99% | negative |
| rs76763715 | GBA | N370S (Gaucher / Parkinson's) | TC | 1/2 |  | carrier: ~5x Parkinson's odds; Gaucher carrier |
| rs34637584 | LRRK2 | G2019S Parkinson's | GG | 0/2 | 100% | negative |
| rs2802292 | FOXO3 | Longevity-associated | TT | 0/2 | 32% | no longevity allele |
| rs1800629 | TNF | -308 promoter | GG | 0/2 | 75% | no variant allele |

### Carrier / rare pathogenic sites typed on arrays

| rsid | gene | label | genotype | effect_count | share | interpretation |
|---|---|---|---|---|---|---|
| rs113993960 | CFTR | deltaF508 cystic fibrosis | II | 0/2 |  | negative |
| rs333 | CCR5 | delta32 | II | 0/2 |  | negative |
| rs80357906 | BRCA1 | 185delAG (Ashkenazi founder) | II | 0/2 |  | negative at this site |
| rs80359550 | BRCA2 | 6174delT (Ashkenazi founder) | II | 0/2 |  | negative at this site |
| rs11571833 | BRCA2 | K3326X | AA | 0/2 | 98% | negative |
| rs334 | HBB | Sickle cell HbS | TT | 0/2 | 100% | negative |
| rs33930165 | HBB | HbC | CC | 0/2 | 100% | negative |
| rs33950507 | HBB | HbE | CC | 0/2 | 100% | negative |

### Pharmacogenomic sites (single markers)

| rsid | gene | label | genotype | effect_count | share | interpretation |
|---|---|---|---|---|---|---|
| rs1799853 | CYP2C9 | *2 (R144C) | TC | 1/2 | 22% | one *2: reduced warfarin/NSAID/phenytoin metabolism |
| rs1057910 | CYP2C9 | *3 (I359L) | CA | 1/2 | 13% | one *3: markedly reduced metabolism |
| rs28371686 | CYP2C9 | *5 | GG | 2/2 | 0% | *5/*5 |
| rs9332131 | CYP2C9 | *6 | II | 0/2 |  | no *6 |
| rs28371685 | CYP2C9 | *11 | CC | 0/2 | 100% | no *11 |
| rs4244285 | CYP2C19 | *2 (loss of function) | GA | 1/2 | 25% | one *2: intermediate metabolizer (clopidogrel less effective) |
| rs4986893 | CYP2C19 | *3 (loss of function) | GA | 1/2 | 0% | one *3 |
| rs12248560 | CYP2C19 | *17 (increased function) | TC | 1/2 | 35% | one *17: rapid metabolizer |
| rs9923231 | VKORC1 | -1639G>A warfarin sensitivity | CC | 0/2 | 37% | normal warfarin dose requirement |
| rs9934438 | VKORC1 | 1173C>T (proxy for -1639) | AA | 2/2 | 15% | sensitive |
| rs2108622 | CYP4F2 | V433M | CT | 1/2 | 41% | slightly higher warfarin dose |
| rs3892097 | CYP2D6 | *4 (null) | CC | 0/2 |  | no *4 |
| rs1065852 | CYP2D6 | *10 (reduced) / also on *4 | GG | 0/2 |  | no *10 |
| rs28371725 | CYP2D6 | *41 (reduced) | CT | 1/2 |  | one *41 |
| rs776746 | CYP3A5 | *3 vs *1 (expressor) | TC | 1/2 | 11% | *1/*3 expressor: higher tacrolimus dose |
| rs3745274 | CYP2B6 | *6 (516G>T) | GG | 0/2 | 58% | *1/*1 |
| rs4149056 | SLCO1B1 | *5 (V174A) statin myopathy | TC | 1/2 | 27% | one C: ~4x simvastatin myopathy risk |
| rs2231142 | ABCG2 | Q141K | TT | 2/2 | 1% | homozygous reduced |
| rs3918290 | DPYD | *2A (IVS14+1G>A) | CC | 0/2 | 99% | normal |
| rs67376798 | DPYD | D949V | TT | 0/2 | 99% | normal |
| rs887829 | UGT1A1 | *28 proxy (Gilbert) | CT | 1/2 | 42% | one *28 |
| rs4148323 | UGT1A1 | *6 (G71R) | GA | 1/2 | 1% | one *6 |
| rs1801280 | NAT2 | *5 (I114T) | TT | 0/2 | 30% | no *5 |
| rs1799930 | NAT2 | *6 (R197Q) | GA | 1/2 | 41% | one slow allele |
| rs1799931 | NAT2 | *7 (G286E) | GA | 1/2 | 4% | one slow allele |
| rs1495741 | NAT2 | acetylator tag | AG | 1/2 | 37% | intermediate |
| rs1800460 | TPMT | *3B | CT | 1/2 | 5% | one *3B (thiopurine dose reduction) |
| rs1142345 | TPMT | *3C | CC | 2/2 | 0% | *3C/*3C |
| rs2395029 | HCP5 | HLA-B*57:01 tag | GT | 1/2 | 8% | HLA-B*57:01 carrier: abacavir hypersensitivity, avoid |
| rs1800566 | NQO1 | P187S (*2) | GA | 1/2 | 33% | one null allele |

## ClinVar scan of every typed site

skipped (--no-clinvar)

## Drug-gene results (PharmCAT)

skipped (--no-pharmcat)

One-page card for a pharmacist: `pgx-card.pdf`

## What this file cannot tell you

- Rare pathogenic variants (the chip targets common variation; BRCA, CFTR, etc. are covered at a handful of sites only).
- Copy-number genes such as CYP2D6, or HLA types beyond tag SNPs.
- Admixture percentages at the resolution of the vendor's own report (that needs their reference panel).
- Anything on sites the vendor did not type; the panel marks those `not on array`. PharmCAT assumes reference at untyped positions.
- For a male sample, X-chromosome sites are hemizygous (one copy), so carrier/homozygous language does not apply there.
- `share` uses the European reference population unless `--ancestry` placed the sample in another 1000 Genomes group.
