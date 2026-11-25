# Tools used under the hood

| Tool | What it contributes | Source |
|---|---|---|
| `snps` | vendor detection, normalization, build remapping | github.com/apriha/snps |
| yhaplo | Y haplogroup (ISOGG 2016 tree) | github.com/23andMe/yhaplo |
| Haplogrep3 + PhyloTree 17 | mtDNA haplogroup | github.com/genepi/haplogrep3 |
| ClinVar | clinically classified variants at every typed site | ncbi.nlm.nih.gov/clinvar |
| PharmCAT | CPIC/PharmGKB star-allele calling | github.com/PharmGKB/PharmCAT |
| ezancestry | 1000 Genomes population classification | github.com/arvkevi/ezancestry |
| Yleaf | second-opinion Y haplogroup on current YFull / FTDNA / ISOGG trees | github.com/genid/Yleaf |
| rCRS (NC_012920.1) | mitochondrial reference | NCBI |
| samtools, bcftools | reads, pileups, indexing (whole-genome mode) | htslib.org |
| ExpansionHunter + stranger | repeat expansions against published ranges | github.com/Illumina/ExpansionHunter, github.com/Clinical-Genomics/stranger |
| PyPGx | CYP2D6 with deletions, duplications and hybrids | github.com/sbslee/pypgx |
| T1K | HLA typing against IPD-IMGT/HLA | github.com/mourisl/T1K |
| SMNCopyNumberCaller | SMN1 and SMN2 copy number | github.com/Illumina/SMNCopyNumberCaller |
| ClinGen dosage sensitivity | which genes matter for copy-number variants | dosage.clinicalgenome.org |
| PGS Catalog + 1000 Genomes 30x | small polygenic scores and their reference distribution | pgscatalog.org |

Evaluated and not used (yet): PLINK 2 (nothing an array file needs that is not already here), ADMIXTURE (needs a multi-GB reference panel),
pgsc_calc genome-wide polygenic scores (need a 7-15 GB reference panel; whole-genome mode ships small scores with a precomputed reference instead), SNPedia (non-commercial license, no bulk use).

Licenses and credits for everything above: [NOTICE.md](../NOTICE.md).

## Related tools

Local: [`apriha/lineage`](https://github.com/apriha/lineage) (shared segments between kits),
[PLINK 2](https://github.com/chrchang/plink-ng), ADMIXTURE + 1000 Genomes, [SNPedia](https://www.snpedia.com).
Web (uploads your data): GEDmatch, James Lick's mthap, Morley's Y-SNP predictor, Promethease, Genetic Genie.
