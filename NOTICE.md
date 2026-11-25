# Credits and third-party licenses

gnomancer is MIT licensed. It downloads or depends on the following, each under its own license. The tools are not bundled in the gnomancer package; `gnomancer setup` fetches them from their official releases. Small derived data tables (gene lists, score weights, site coordinates) ship inside the package and are marked as such below.

| Component | Role | License | Source |
|---|---|---|---|
| snps | parsing, normalization, build remapping | BSD-3-Clause | https://github.com/apriha/snps |
| lineage | shared-DNA segments in family mode | MIT | https://github.com/apriha/lineage |
| yhaplo | Y haplogroup | yhaplo license (non-commercial research use; see upstream) | https://github.com/23andMe/yhaplo |
| Yleaf | Y haplogroup on current trees | GPL-3.0 (run as a separate process in its own venv) | https://github.com/genid/Yleaf |
| Haplogrep3 | mtDNA haplogroup | MIT | https://github.com/genepi/haplogrep3 |
| PharmCAT | pharmacogenomic calling | MPL-2.0 | https://github.com/PharmGKB/PharmCAT |
| ezancestry | population estimate | MIT | https://github.com/arvkevi/ezancestry |
| ClinVar | variant classifications | public domain (NCBI) | https://www.ncbi.nlm.nih.gov/clinvar/ |
| 1000 Genomes allele frequencies (via Ensembl) | "how common" column | open data (Fort Lauderdale principles) | https://www.internationalgenome.org |
| rCRS NC_012920.1 | mitochondrial reference | public domain (NCBI) | https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1 |
| reportlab | pharmacist card PDF | BSD | https://www.reportlab.com/opensource/ |
| samtools, bcftools, htslib | reading BAM/CRAM, pileups (whole-genome mode) | MIT/Expat | https://www.htslib.org |
| ExpansionHunter | repeat-expansion sizing | see upstream (bioconda build of v5.0.0) | https://github.com/Illumina/ExpansionHunter |
| stranger | repeat catalog and normal / disease thresholds | MIT | https://github.com/Clinical-Genomics/stranger |
| PyPGx and pypgx-bundle | CYP2D6 from read depth | MIT | https://github.com/sbslee/pypgx |
| T1K | HLA typing | MIT | https://github.com/mourisl/T1K |
| IPD-IMGT/HLA | HLA allele database used by T1K | CC BY-ND 4.0 | https://www.ebi.ac.uk/ipd/imgt/hla/ |
| SMNCopyNumberCaller | SMN1/SMN2 copy number | PolyForm Strict 1.0.0 (non-commercial use; installed by the user's own setup, not redistributed) | https://github.com/Illumina/SMNCopyNumberCaller |
| GRCh38 analysis set | reference genome | public domain (UCSC / GRC) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/ |
| GENCODE | gene coordinates for the HLA index | open (EMBL-EBI terms) | https://www.gencodegenes.org |
| ClinGen dosage sensitivity | which genes and regions are copy-number sensitive (bundled as `data/dosage_genes.tsv`) | CC0 | https://dosage.clinicalgenome.org |
| ACMG SF v3.3 gene list | medically actionable genes (bundled as `data/acmg_sf.tsv`) | published list; Genet Med 2025;27:101454 | https://www.acmg.net |
| PGS Catalog scores | polygenic score variants and weights (bundled under `data/pgs/`, each with its citation) | per score; only scores under the Catalog's default terms are bundled | https://www.pgscatalog.org |
| 1000 Genomes 30x panel | reference distribution for polygenic scores | open data | https://www.internationalgenome.org |
| pandas, numpy | data handling | BSD-3-Clause | https://pandas.pydata.org, https://numpy.org |

This product includes software developed by 23andMe, Inc. (yhaplo).

The ClinVar false-positive figure quoted in the report is from Tandy-Connor et al., *Genetics in Medicine* 2018, doi:10.1038/gim.2018.38.
