# Vendors and quirks

Parsing is delegated to [`snps`](https://github.com/apriha/snps), which detects the vendor from the header and
normalizes everything. Quirks that are handled:

| Vendor | Quirk | Handling |
|---|---|---|
| AncestryDNA | CRLF line endings; chromosomes coded 23=X, 24=Y, 25=PAR, 26=MT; two allele columns; no-call `0`; indels `I`/`D`; diploid-looking `GG` on Y and MT | normalized to X/Y/PAR/MT, haploid single letters, NaN no-calls |
| AncestryDNA | V1 (~700k sites, 2012–2016) and V2 (~665k, 2016+) type different sites | report marks `not on array` per site |
| 23andMe | single genotype column; `--` no-call; haploid `A`; internal `i3000001` IDs; v3/v4/v5 chips differ (v5 dropped many PGx sites) | handled; missing sites reported |
| MyHeritage, FamilyTreeDNA | quoted CSV `RSID,CHROMOSOME,POSITION,RESULT` | handled |
| FamilyTreeDNA (old) | build 36 exports | remapped to 37; mtDNA call skipped (build 36 MT is not rCRS) |
| LivingDNA, Genes for Good, tellmeGen, Circle | tab/CSV variants | handled |
| Build 38 exports, VCF | | remapped to GRCh37 (chain file downloaded once) |
| Any | occasional probes reported on the reverse strand | panel auto-flips unambiguous SNPs and flags A/T and C/G ones |
| Any | MT numbering is rCRS (NC_012920) on builds 37/38 | diffed against the bundled rCRS |
