# Whole genomes

gnomancer reads a sequencing VCF (whole genome or exome, plain or gVCF, GRCh38 or GRCh37) natively. Nothing is
converted to an array-shaped file first, so the things sequencing adds are kept: insertions and deletions, exact
alleles, per-call evidence, and the difference between "matches the reference" and "was not looked at".

```bash
gnomancer genome.vcf.gz                                          # the VCF alone: about a minute
gnomancer genome.vcf.gz --sv genome.sv.vcf.gz --cnv genome.cnv.vcf.gz
gnomancer genome.vcf.gz --bam genome.bam                         # adds the read-level analyses below: about five minutes in all
```

The file type is detected from its contents; there is no separate command.

## What each input adds

| You have | You get |
|---|---|
| The VCF | Quality summary, runs of homozygosity, the marker panel, combined readings, a ClinVar scan of every variant (indels included) with the 84 ACMG actionable genes set apart, PharmCAT on every defining position, mtDNA haplogroup from the whole mitochondrial sequence, Y haplogroup, population estimate, polygenic scores |
| A gVCF instead of a VCF | The same, with every "reference" genotype observed in the file's reference blocks instead of inferred |
| `--sv`, `--cnv` | The vendor's structural-variant and copy-number calls checked against ClinGen's dosage-sensitive genes and recurrent regions, the ACMG genes, and a few drug and carrier genes |
| `--bam` (BAM or CRAM) | Depth and sex from the reads, every named position the VCF is silent about confirmed in the reads, Y haplogroup read marker by marker, low-level mitochondrial heteroplasmy, and the four callers below |

With `--bam`, four specialist callers run unless you limit them with `--callers` (a comma-separated list, or `none`):

| `--callers` name | Tool | Answers | Time on a 35x genome, Apple M-series |
|---|---|---|---|
| `repeats` | ExpansionHunter + stranger | 51 repeat-expansion loci (Huntington, fragile X, myotonic dystrophy, the ataxias, ...) against published normal and disease ranges | 2 min |
| `cyp2d6` | PyPGx | CYP2D6 star alleles including deletion, duplication and hybrids; PharmCAT then assigns the phenotype | under 1 min |
| `hla` | T1K | Two-field types for HLA-A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1, with the drug-reaction alleles, disease-linked alleles and celiac types named | 2 min |
| `smn` | SMNCopyNumberCaller | SMN1 and SMN2 copy number (spinal muscular atrophy carrier status) | under 1 min |

## Setup

The VCF path needs nothing beyond the normal first-run download plus ClinVar for the file's build (fetched
automatically, about 200 MB). The Y call from a VCF, and everything behind `--bam`, need:

```bash
gnomancer setup --wgs
```

This builds isolated conda environments under `~/.gnomancer` (samtools, bcftools and the four callers; an existing
conda, mamba or micromamba is used, otherwise micromamba is downloaded there), fetches the GRCh38 reference genome
(950 MB download, 3.1 GB unpacked), PyPGx's data bundle and the HLA index. About 2.5 GB to download and 6 GB on disk.
Nothing is installed system-wide. A caller that has no build for your platform is skipped and the report lists
its topic under "what this file cannot tell you" instead.

For a GRCh37 BAM or CRAM, point `GNOMANCER_REF_37` at an indexed FASTA. PharmCAT, polygenic scores and HLA typing
are GRCh38-only.

## Getting the files

Vendors differ in what they hand over. Sequencing.com delivers FASTQ, a GRCh38 VCF and separate structural-variant
and copy-number VCFs, and provides the BAM and a mitochondrial heteroplasmy VCF on request: ask for the BAM. Nebula /
DNA Complete and Dante Labs deliver a CRAM or BAM with the VCF. If you only have FASTQ, gnomancer does not align or
call variants; that takes hours and is a different job (DeepVariant, or nf-core/sarek, will produce the VCF and BAM).

## How the genotypes are decided

**Quality gates.** A call is used only if the caller passed it (FILTER is `PASS` or empty), at least 8 reads cover it,
genotype quality is 20 or more, and the share of reads carrying the variant fits the genotype: 20-80% for a
heterozygous call, 85% or more for a homozygous one. On the mitochondrion a low share is heteroplasmy, not an error.
Calls that fail are counted in the quality section and used nowhere.

**Alleles.** Multi-allelic records are split, and each alternate allele is trimmed of the bases it shares with the
reference. ClinVar's VCF is left-aligned (none of its 318,000 indels can be shifted left), and so is caller output,
so equal variants compare equal.

**Positions with no variant.** A plain VCF lists only differences. For a whole genome an absent position is taken as
homozygous reference; the marker panel marks those genotypes with a small ° and the TSV's `basis` column spells it out.
A gVCF or `--bam` replaces the inference with an observation. For anything smaller than a whole genome (under two
million SNVs) absence means unknown, and such positions are left as no-calls.

One exception guards against a false alarm. At a few sites the reference genome itself carries the allele of
interest (GRCh37 has factor V Leiden). Absence there would imply a positive finding with no evidence behind it, and a
coverage gap looks identical, so when that allele is rare the site is reported as a no-call unless a gVCF or the
reads settle it.

**Sex.** From the share of heterozygous calls on X outside the pseudo-autosomal regions; with `--bam`, also from X and
Y depth, and a disagreement between the two is flagged as a possible sample mix-up.

## The ClinVar scan

Every variant is matched to ClinVar by chromosome, position and both alleles. Matches classified pathogenic or likely
pathogenic are grouped by what they mean:

1. **In a medically actionable gene**: the ACMG secondary-findings list, v3.3 (84 genes), with the inheritance mode, so a
   single copy in a recessive gene reads as a carrier state.
2. **On both copies, or the only copy**: homozygous, hemizygous or mitochondrial.
3. **Single copies**: carrier states for recessive conditions.
4. **Common or low-penetrance alleles**: real calls that shift odds without predicting disease.
5. **Set aside**: a weak call or a classification nobody has reviewed. Listed, not hidden.

Two pathogenic variants in one gene are flagged, with their phase when the VCF provides one. A typical genome shows
one to three solid findings, nearly all of them carrier states.

Not covered: variants ClinVar has not classified, multi-nucleotide changes that callers split into separate SNVs, and
structural variants (use `--sv` / `--cnv`).

## Polygenic scores

Six small published scores (150 to 313 variants each) are bundled together with their distribution in the unrelated
1000 Genomes samples, so a score can be placed on a percentile with no download: breast cancer (for female samples),
prostate cancer (male), coronary artery disease, atrial fibrillation, type 2 diabetes and venous thromboembolism. Each
row names its paper and PGS Catalog ID. A score is bundled only if at least 85% of its published variants could be
matched in the reference panel; `examples/build_pgs.py` rebuilds them. The reference group is the
1000 Genomes super-population the ancestry estimate chooses. Genome-wide scores with millions of variants predict
better and need a multi-gigabyte reference panel; they are out of scope. A percentile is a position in a distribution,
not a probability, and the scores transfer imperfectly outside the ancestry they were developed in.

## Limits

- Research-grade 30x sequencing is not a clinical test. Each read-level section names the clinical test that confirms it.
- Short reads cannot size the largest repeat expansions, phase distant variants, or resolve every CYP2D6 hybrid.
- The reads are checked at the named positions (marker panel, PharmCAT, polygenic scores), not at every ClinVar position.
- Structural-variant calls from short reads contain many false positives, and overlap is measured against a gene's span,
  not specifically its coding sequence.

## Privacy

Nothing about the sample leaves the machine. Two details: reading a CRAM that was compressed against a reference
other than the bundled one makes samtools fetch the missing reference sequences by checksum from the European
Nucleotide Archive (the request contains the checksum of a public reference sequence, nothing from the sample); set
`REF_PATH` to a local cache to prevent it. And `gnomancer setup --wgs` downloads public tools and references only.

## Try it without your own data

```bash
python examples/make_wgs_demo.py demo.vcf.gz          # a synthetic whole genome with known variants planted
gnomancer demo.vcf.gz
```

The public Genome in a Bottle samples are real whole genomes with open consent; Google's DeepVariant case studies
publish a VCF, gVCF and BAM for HG003 that exercise every step here. On that sample (35x, 46 GB BAM, Apple M-series
laptop) the VCF alone takes about 40 seconds, the gVCF two and a half minutes, and the VCF with the BAM and all four
callers four and a half minutes. gnomancer recovers what is known about the sample: Y haplogroup J1-L816, mtDNA
K1a1b1a, two textbook Ashkenazi HLA haplotypes, and the factor XI founder variant among its carrier findings.
